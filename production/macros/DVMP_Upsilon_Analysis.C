#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <numeric>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <memory>
#include <stdexcept>
#include <utility>

#include <TFile.h>
#include <TH1.h>
#include <TTree.h>
#include <TChain.h>
#include <TMath.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TROOT.h>

#include <podio/ObjectID.h>

#ifdef __CLING__
R__LOAD_LIBRARY(libpodio)
R__LOAD_LIBRARY(libpodioDict)
#endif

#include "ePICStyle.C"
#include "ElectronFinder.C" 

const double eMass = 0.000511;
const double pMass = 0.938272;
const double Upsilon1SMass = 9.46030;  
const double Upsilon2SMass = 10.02326; 
const double Upsilon3SMass = 10.35520; 

const double E_e_beam = 9.0;//10.0;
const double E_p_beam = 275.0;//100.0;

double GetMinUpsilonDiff(double mass) {
    double d1 = fabs(mass - Upsilon1SMass);
    double d2 = fabs(mass - Upsilon2SMass);
    double d3 = fabs(mass - Upsilon3SMass);
    return std::min({d1, d2, d3});
}

namespace UpsilonFEMC {
// ELKE_HALFPLANE_LEGACY_V1. Existing mode numbers remain unchanged.
enum Scenario { Full = 0, RemoveInner = 1, RemoveOuter = 2,
                KeepXPos = 3, KeepYNeg = 4 };

inline bool IsHalfPlane(int scenario) {
    return scenario == KeepXPos || scenario == KeepYNeg;
}

// Apply only to EcalEndcapPClusters, using signed position eta from the origin.
// Inner (near the beam pipe) means higher eta. The boundary belongs to the inner part.
inline bool KeepCluster(double eta, int scenario, double etaSplit) {
    if (scenario == Full) return true;
    return scenario == RemoveInner ? eta < etaSplit : eta >= etaSplit;
}

// Position-based dispatcher. ScopedEnergyMask validates the coordinates first.
// Use the stored EcalEndcapPClusters.position coordinates directly, with the
// split at x=0 or y=0. No rotation, recentering, momentum cut, or eta window is
// applied in the half-plane modes. "Positive x" avoids a viewing-dependent
// left/right name; negative y is the requested bottom-half convention.
inline bool KeepClusterAtPosition(double x, double y, double z,
                                  int scenario, double etaSplit) {
    switch (scenario) {
        case Full:      return true;
        case KeepXPos:  return x >= 0.0; // The x=0 boundary is retained.
        case KeepYNeg:  return y <  0.0; // The y=0 boundary is removed.
        case RemoveInner:
        case RemoveOuter:
            return KeepCluster(std::asinh(z / std::hypot(x, y)), scenario, etaSplit);
        default:
            throw std::invalid_argument("Unknown fEMC mode in position-based mask.");
    }
}

// Temporarily mask energies, preserving cluster indices and the existing association
// logic in ElectronFinder.C. Restore the exact input values after precomputation.
// Only the in-memory current entry is touched; the input ROOT files are read-only.
class ScopedEnergyMask {
public:
    ScopedEnergyMask(TTreeReaderArray<float>& energies,
                     const TTreeReaderArray<float>* posX,
                     const TTreeReaderArray<float>* posY,
                     const TTreeReaderArray<float>* posZ,
                     int scenario, double etaSplit) : energies_(energies) {
        if (scenario == Full) return;
        if (!posX || !posY || !posZ || posX->GetSetupStatus() < 0 ||
            posY->GetSetupStatus() < 0 || posZ->GetSetupStatus() < 0) {
            throw std::runtime_error("fEMC descoping requires readable float EcalEndcapPClusters.position.{x,y,z} branches.");
        }
        const std::size_t n = energies_.GetSize();
        if (posX->GetSize() != n || posY->GetSize() != n || posZ->GetSize() != n) {
            throw std::runtime_error("fEMC cluster energy and position arrays have different sizes.");
        }
        // Validate and save all entries before modifying any energy.
        for (std::size_t i = 0; i < n; ++i) {
            const double x = (*posX)[i], y = (*posY)[i], z = (*posZ)[i];
            const double r = std::hypot(x, y);
            if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z) ||
                r <= 0.0 || z <= 0.0) {
                throw std::runtime_error("Invalid fEMC cluster position at index " + std::to_string(i));
            }
            if (!KeepClusterAtPosition(x, y, z, scenario, etaSplit))
                saved_.emplace_back(i, energies_[i]);
        }
        for (const auto& entry : saved_) energies_[entry.first] = 0.0f;
    }

    ~ScopedEnergyMask() {
        for (const auto& entry : saved_) energies_[entry.first] = entry.second;
    }

    ScopedEnergyMask(const ScopedEnergyMask&) = delete;
    ScopedEnergyMask& operator=(const ScopedEnergyMask&) = delete;
    std::size_t RemovedCount() const { return saved_.size(); }
    // Expose the same geometric decision to the reconstructed-matching path.
    // A physical zero-energy cluster is not mistaken for a geometrically masked one.
    std::vector<unsigned char> ActiveMask() const {
        std::vector<unsigned char> active(energies_.GetSize(), 1);
        for (const auto& entry : saved_) active[entry.first] = 0;
        return active;
    }

private:
    TTreeReaderArray<float>& energies_;
    std::vector<std::pair<std::size_t, float>> saved_;
};
} // namespace UpsilonFEMC

// fEMC acceptance study on existing reconstructed events, not a geometry/reconstruction rerun.
//   full          : original calorimeter coverage (default)
//   remove_inner  : retain fEMC clusters with eta <  femc_eta_split
//   remove_outer  : retain fEMC clusters with eta >= femc_eta_split
//   keep_xpos     : retain fEMC clusters with position.x >= 0 (legacy ID only)
//   keep_yneg     : retain fEMC clusters with position.y <  0 (legacy ID only)
// The half-plane modes cover the full radial extent. They replace, rather than
// combine with, an inner/outer mask; femc_eta_split is not used by their decision.
// They are coordinate-based cluster masks, not an assertion of 50% of channels.
// The default 1.68 split follows the October 2026 email; it is configurable.
// The radial modes are complementary eta regions, not a calculation of removed physical area.
// Descoped modes require float EcalEndcapPClusters.position.{x,y,z} branches.
// Recomputing isolation and E/p can improve some candidates, so yields need not be monotonic.
// Electron-ID mode is an independent SIXTH argument, after the existing fEMC arguments.
//   legacy (default): preserve the original MC-assisted energy associations.
//   reco            : use reconstructed TrackClusterMatches and particle-to-track links.
// Both modes use the unchanged 0.8 <= E/p <= 1.2 and track-linked isolation cuts.
// The reco path filters stored matches; it does not rerun propagation or matching.
// Examples (use the updated ElectronFinder.C and your existing ePICStyle.C):
// root -l -b -q 'DVMP_Upsilon_Analysis.C("input.root","full.root","new")'
// root -l -b -q 'DVMP_Upsilon_Analysis.C("input.root","inner_removed.root","new","remove_inner")'
// root -l -b -q 'DVMP_Upsilon_Analysis.C("input.root","outer_removed.root","new","remove_outer")'
// root -l -b -q 'DVMP_Upsilon_Analysis.C("input.root","xpos.root","new","keep_xpos")'
// root -l -b -q 'DVMP_Upsilon_Analysis.C("input.root","yneg.root","new","keep_yneg")' 
// root -l -b -q 'DVMP_Upsilon_Analysis.C("input.root","reco_full.root","new","full",1.68,"reco")'
// root -l -b -q 'DVMP_Upsilon_Analysis.C("input.root","reco_inner.root","new","remove_inner",1.68,"reco")'
// This patch does not change the existing beam constants, weights, proton selection, or t estimators.
void DVMP_Upsilon_Analysis(TString infile, std::string outfilename, TString data_version = "new",
                          TString femc_scenario = "full", double femc_eta_split = 1.68,
                          TString ele_id_mode = "legacy")
{
    gROOT->SetBatch(kTRUE);

    ele_id_mode.ToLower();
    int b_ele_id_mode = UpsilonEleID::Legacy;
    if (ele_id_mode == "reco") b_ele_id_mode = UpsilonEleID::Reco;
    else if (ele_id_mode != "legacy")
        throw std::invalid_argument("Unknown ele_id_mode. Use legacy or reco.");

    femc_scenario.ToLower();
    int femc_mode = UpsilonFEMC::Full;
    if (femc_scenario == "remove_inner") femc_mode = UpsilonFEMC::RemoveInner;
    else if (femc_scenario == "remove_outer") femc_mode = UpsilonFEMC::RemoveOuter;
    else if (femc_scenario == "keep_xpos") femc_mode = UpsilonFEMC::KeepXPos;
    else if (femc_scenario == "keep_yneg") femc_mode = UpsilonFEMC::KeepYNeg;
    else if (femc_scenario != "full") {
        throw std::invalid_argument("Unknown fEMC scenario. Use full, remove_inner, remove_outer, keep_xpos, or keep_yneg.");
    }
    if (UpsilonFEMC::IsHalfPlane(femc_mode) && b_ele_id_mode != UpsilonEleID::Legacy) {
        throw std::invalid_argument("keep_xpos and keep_yneg are only implemented for legacy electron ID; use ele_id_mode=legacy.");
    }
    if (!std::isfinite(femc_eta_split) || femc_eta_split <= 0.0) {
        throw std::invalid_argument("femc_eta_split must be finite and positive.");
    }

    bool isNewFormat = false;
    data_version.ToLower(); 
    if (data_version == "new") isNewFormat = true;
    else if (data_version == "old") isNewFormat = false;
    else return; 

    double cross_section_fb = 0.0;
    TString infile_lower = infile;
    infile_lower.ToLower();
    if (infile_lower.Contains("1s")) cross_section_fb = 327.584; //, (104.793); 
    else if (infile_lower.Contains("2s")) cross_section_fb = 84.780;  
    else if (infile_lower.Contains("3s")) cross_section_fb = 68.929;  
    else cross_section_fb = 1.0;

    std::unique_ptr<TChain> mychain(new TChain("events"));
    mychain->Add(infile);
    const Long64_t total_entries = mychain->GetEntries();
    if (total_entries <= 0) {
        throw std::runtime_error("No input events found; no output file was created.");
    }
    // The default scenario keeps the original branch requirements.
    if (femc_mode != UpsilonFEMC::Full) {
        mychain->LoadTree(0);
        for (const char* branch : {"EcalEndcapPClusters.position.x",
                                   "EcalEndcapPClusters.position.y",
                                   "EcalEndcapPClusters.position.z"}) {
            if (!mychain->GetBranch(branch)) {
                throw std::runtime_error(std::string("Missing branch required for fEMC descoping: ") + branch);
            }
        }
    }
    // Validate the reconstruction-only input before creating the output file.
    UpsilonEleID::CollectionIDs firstFileCollectionIDs;
    if (b_ele_id_mode == UpsilonEleID::Reco) {
        mychain->LoadTree(0);
        UpsilonEleID::RequireRecoBranches(mychain->GetTree());
        firstFileCollectionIDs = UpsilonEleID::RequiredCollectionIDs(mychain->GetCurrentFile());
    }
    std::unique_ptr<TFile> ofile(TFile::Open(outfilename.c_str(), "RECREATE"));
    if (!ofile || ofile->IsZombie()) throw std::runtime_error("Cannot create output ROOT file: " + outfilename);
    ofile->cd();
    ResetElectronIDHistograms();
    TTreeReader tree_reader(mychain.get());

    std::cout << "fEMC scenario: " << femc_scenario << " (mode " << femc_mode << ")";
    if (femc_mode == UpsilonFEMC::KeepXPos)
        std::cout << "; keep stored cluster x >= 0; full radial extent; eta split unused";
    else if (femc_mode == UpsilonFEMC::KeepYNeg)
        std::cout << "; keep stored cluster y < 0; full radial extent; eta split unused";
    else
        std::cout << ", eta split = " << femc_eta_split;
    std::cout << std::endl;
    std::cout << "Electron ID: " << ele_id_mode << " (mode " << b_ele_id_mode << ")"
              << "; 0.8 <= E/p <= 1.2; track-linked isolation >= 0.9 in R=0.4" << std::endl;
    if (b_ele_id_mode == UpsilonEleID::Reco)
        std::cout << "Reco ID uses stored reconstructed matches only, with no MC-association fallback." << std::endl;

    double expected_lumi_invfb = 10.0;
    double total_gen_events = total_entries;
    double global_weight = (cross_section_fb * expected_lumi_invfb) / total_gen_events;
    
    TH1D* h_EventCount = new TH1D("h_EventCount", "Total Generated Events", 1, 0, 1);
    h_EventCount->SetBinContent(1, total_gen_events);

    const int nCuts = 10;
    TH1D* h_CutFlow = new TH1D("h_CutFlow", "Event Selection Cutflow", nCuts, 0.5, nCuts + 0.5);
    const char* cutNames[nCuts] = {
        "1. Total Processed", "2. Gen Upsilon Found", "3. Gen Leps in Acc.", 
        "4. Reco e+ Candidate", "5. Reco e- Decay Cand", "6. Reco Scat. e- Cand", 
        "7. Mass Window [8, 12]", "8. #Sigma(E-Pz) Exclusivity", "9. W Phase Space Cut", "10. Proton Tagged"
    };
    for (int i = 0; i < nCuts; i++) h_CutFlow->GetXaxis()->SetBinLabel(i+1, cutNames[i]);

    TTree *outTree = new TTree("UpsilonTree", "Flattened Kinematics for Post-Processing");

    int b_femc_nClustersRemoved = 0;
    // Keep configuration per event, including after hadd, without extra file-level keys.
    outTree->Branch("femc_mode", &femc_mode, "femc_mode/I");
    outTree->Branch("femc_eta_split", &femc_eta_split, "femc_eta_split/D");
    outTree->Branch("femc_nClustersRemoved", &b_femc_nClustersRemoved, "femc_nClustersRemoved/I");

    // Lightweight event-level provenance and reconstructed-match diagnostics.
    // Legacy mode uses -1 for diagnostics that do not apply. No new file-level keys.
    int b_ele_id_nMatched = -1, b_ele_id_nMissing = -1;
    int b_ele_id_nMasked = -1, b_ele_id_nAmbiguous = -1;
    outTree->Branch("ele_id_mode", &b_ele_id_mode, "ele_id_mode/I");
    outTree->Branch("ele_id_nMatched", &b_ele_id_nMatched, "ele_id_nMatched/I");
    outTree->Branch("ele_id_nMissing", &b_ele_id_nMissing, "ele_id_nMissing/I");
    outTree->Branch("ele_id_nMasked", &b_ele_id_nMasked, "ele_id_nMasked/I");
    outTree->Branch("ele_id_nAmbiguous", &b_ele_id_nAmbiguous, "ele_id_nAmbiguous/I");

    double b_weight, b_weight_proj; 
    double b_true_mass, b_true_y, b_true_pT, b_true_W, b_true_Q2, b_true_t;
    double b_true_ePlus_pT, b_true_ePlus_eta;
    double b_true_decay_eMinus_pT, b_true_decay_eMinus_eta;
    double b_true_scat_eMinus_pT, b_true_scat_eMinus_eta, b_true_scat_eMinus_E;
    int b_true_nTracks_Central, b_true_nTracks_Tagger;

    double b_reco_mass, b_reco_y, b_reco_pT, b_reco_W, b_reco_Q2;
    double b_reco_t_eX, b_reco_t_BABE, b_reco_t_eXBABE; 
    double b_reco_t_n, b_reco_t_x, b_reco_omega, b_reco_t_proj; 
    double b_reco_ePlus_pT, b_reco_ePlus_eta;
    double b_reco_decay_eMinus_pT, b_reco_decay_eMinus_eta;
    double b_reco_scat_eMinus_pT, b_reco_scat_eMinus_eta, b_reco_scat_eMinus_E;
    double b_reco_tagged_proton_p;
    double b_sigma_epz;
    int b_reco_nTracks_Central, b_reco_nTracks_Tagger;
    int b_reco_topology;

    outTree->Branch("weight", &b_weight, "weight/D");
    outTree->Branch("weight_proj", &b_weight_proj, "weight_proj/D"); 
    outTree->Branch("reco_topology", &b_reco_topology, "reco_topology/I");
    
    outTree->Branch("true_mass", &b_true_mass, "true_mass/D");
    outTree->Branch("true_y", &b_true_y, "true_y/D");
    outTree->Branch("true_pT", &b_true_pT, "true_pT/D");
    outTree->Branch("true_W", &b_true_W, "true_W/D");
    outTree->Branch("true_Q2", &b_true_Q2, "true_Q2/D");
    outTree->Branch("true_t", &b_true_t, "true_t/D");
    outTree->Branch("true_ePlus_pT", &b_true_ePlus_pT, "true_ePlus_pT/D");
    outTree->Branch("true_ePlus_eta", &b_true_ePlus_eta, "true_ePlus_eta/D");
    outTree->Branch("true_decay_eMinus_pT", &b_true_decay_eMinus_pT, "true_decay_eMinus_pT/D");
    outTree->Branch("true_decay_eMinus_eta", &b_true_decay_eMinus_eta, "true_decay_eMinus_eta/D");
    outTree->Branch("true_scat_eMinus_pT", &b_true_scat_eMinus_pT, "true_scat_eMinus_pT/D");
    outTree->Branch("true_scat_eMinus_eta", &b_true_scat_eMinus_eta, "true_scat_eMinus_eta/D");
    outTree->Branch("true_scat_eMinus_E", &b_true_scat_eMinus_E, "true_scat_eMinus_E/D");
    outTree->Branch("true_nTracks_Central", &b_true_nTracks_Central, "true_nTracks_Central/I");
    outTree->Branch("true_nTracks_Tagger", &b_true_nTracks_Tagger, "true_nTracks_Tagger/I");

    outTree->Branch("reco_mass", &b_reco_mass, "reco_mass/D");
    outTree->Branch("reco_y", &b_reco_y, "reco_y/D");
    outTree->Branch("reco_pT", &b_reco_pT, "reco_pT/D");
    outTree->Branch("reco_W", &b_reco_W, "reco_W/D");
    outTree->Branch("reco_Q2", &b_reco_Q2, "reco_Q2/D");
    outTree->Branch("reco_t_eX", &b_reco_t_eX, "reco_t_eX/D");
    outTree->Branch("reco_t_BABE", &b_reco_t_BABE, "reco_t_BABE/D");
    outTree->Branch("reco_t_eXBABE", &b_reco_t_eXBABE, "reco_t_eXBABE/D");
    outTree->Branch("reco_t_n", &b_reco_t_n, "reco_t_n/D");
    outTree->Branch("reco_t_x", &b_reco_t_x, "reco_t_x/D");
    outTree->Branch("reco_omega", &b_reco_omega, "reco_omega/D");
    outTree->Branch("reco_t_proj", &b_reco_t_proj, "reco_t_proj/D"); 
    outTree->Branch("reco_ePlus_pT", &b_reco_ePlus_pT, "reco_ePlus_pT/D");
    outTree->Branch("reco_ePlus_eta", &b_reco_ePlus_eta, "reco_ePlus_eta/D");
    outTree->Branch("reco_decay_eMinus_pT", &b_reco_decay_eMinus_pT, "reco_decay_eMinus_pT/D");
    outTree->Branch("reco_decay_eMinus_eta", &b_reco_decay_eMinus_eta, "reco_decay_eMinus_eta/D");
    outTree->Branch("reco_scat_eMinus_pT", &b_reco_scat_eMinus_pT, "reco_scat_eMinus_pT/D");
    outTree->Branch("reco_scat_eMinus_eta", &b_reco_scat_eMinus_eta, "reco_scat_eMinus_eta/D");
    outTree->Branch("reco_scat_eMinus_E", &b_reco_scat_eMinus_E, "reco_scat_eMinus_E/D");
    outTree->Branch("reco_tagged_proton_p", &b_reco_tagged_proton_p, "reco_tagged_proton_p/D");
    outTree->Branch("sigma_epz", &b_sigma_epz, "sigma_epz/D");
    outTree->Branch("reco_nTracks_Central", &b_reco_nTracks_Central, "reco_nTracks_Central/I");
    outTree->Branch("reco_nTracks_Tagger", &b_reco_nTracks_Tagger, "reco_nTracks_Tagger/I");

    // READERS SETUP
    TTreeReaderArray<int> partGenStat(tree_reader, "MCParticles.generatorStatus");
    TTreeReaderArray<double> partMomX(tree_reader, "MCParticles.momentum.x");
    TTreeReaderArray<double> partMomY(tree_reader, "MCParticles.momentum.y");
    TTreeReaderArray<double> partMomZ(tree_reader, "MCParticles.momentum.z");
    TTreeReaderArray<int> partPdg(tree_reader, "MCParticles.PDG");

    TTreeReaderArray<float> trackMomX(tree_reader, "ReconstructedChargedParticles.momentum.x");
    TTreeReaderArray<float> trackMomY(tree_reader, "ReconstructedChargedParticles.momentum.y");
    TTreeReaderArray<float> trackMomZ(tree_reader, "ReconstructedChargedParticles.momentum.z");
    TTreeReaderArray<float> trackCharge(tree_reader, "ReconstructedChargedParticles.charge");
    TTreeReaderArray<float> trackEng(tree_reader, "ReconstructedChargedParticles.energy");

    TTreeReaderArray<podio::ObjectID>* recoAssoc_new = nullptr;
    TTreeReaderArray<podio::ObjectID>* simuAssoc_new = nullptr;
    TTreeReaderArray<podio::ObjectID>* recoAssocEcalBarrel_new = nullptr;
    TTreeReaderArray<podio::ObjectID>* simuAssocEcalBarrel_new = nullptr;
    TTreeReaderArray<podio::ObjectID>* recoAssocEcalEndcapP_new = nullptr;
    TTreeReaderArray<podio::ObjectID>* simuAssocEcalEndcapP_new = nullptr;
    TTreeReaderArray<podio::ObjectID>* recoAssocEcalEndcapN_new = nullptr;
    TTreeReaderArray<podio::ObjectID>* simuAssocEcalEndcapN_new = nullptr;

    TTreeReaderArray<unsigned int>* recoAssoc_old = nullptr;
    TTreeReaderArray<unsigned int>* simuAssoc_old = nullptr;
    TTreeReaderArray<unsigned int>* recoAssocEcalBarrel_old = nullptr;
    TTreeReaderArray<unsigned int>* simuAssocEcalBarrel_old = nullptr;
    TTreeReaderArray<unsigned int>* recoAssocEcalEndcapP_old = nullptr;
    TTreeReaderArray<unsigned int>* simuAssocEcalEndcapP_old = nullptr;
    TTreeReaderArray<unsigned int>* recoAssocEcalEndcapN_old = nullptr;
    TTreeReaderArray<unsigned int>* simuAssocEcalEndcapN_old = nullptr;

    // Do not bind or read MC association branches in reco-ID mode.
    if (b_ele_id_mode == UpsilonEleID::Legacy) {
        if (isNewFormat) {
            recoAssoc_new = new TTreeReaderArray<podio::ObjectID>(tree_reader, "_ReconstructedChargedParticleAssociations_rec");
            simuAssoc_new = new TTreeReaderArray<podio::ObjectID>(tree_reader, "_ReconstructedChargedParticleAssociations_sim");
            recoAssocEcalBarrel_new = new TTreeReaderArray<podio::ObjectID>(tree_reader, "_EcalBarrelClusterAssociations_rec");
            simuAssocEcalBarrel_new = new TTreeReaderArray<podio::ObjectID>(tree_reader, "_EcalBarrelClusterAssociations_sim");
            recoAssocEcalEndcapP_new = new TTreeReaderArray<podio::ObjectID>(tree_reader, "_EcalEndcapPClusterAssociations_rec");
            simuAssocEcalEndcapP_new = new TTreeReaderArray<podio::ObjectID>(tree_reader, "_EcalEndcapPClusterAssociations_sim");
            recoAssocEcalEndcapN_new = new TTreeReaderArray<podio::ObjectID>(tree_reader, "_EcalEndcapNClusterAssociations_rec");
            simuAssocEcalEndcapN_new = new TTreeReaderArray<podio::ObjectID>(tree_reader, "_EcalEndcapNClusterAssociations_sim");
        } else {
            recoAssoc_old = new TTreeReaderArray<unsigned int>(tree_reader, "ReconstructedChargedParticleAssociations.recID");
            simuAssoc_old = new TTreeReaderArray<unsigned int>(tree_reader, "ReconstructedChargedParticleAssociations.simID");
            recoAssocEcalBarrel_old = new TTreeReaderArray<unsigned int>(tree_reader, "EcalBarrelClusterAssociations.recID");
            simuAssocEcalBarrel_old = new TTreeReaderArray<unsigned int>(tree_reader, "EcalBarrelClusterAssociations.simID");
            recoAssocEcalEndcapP_old = new TTreeReaderArray<unsigned int>(tree_reader, "EcalEndcapPClusterAssociations.recID");
            simuAssocEcalEndcapP_old = new TTreeReaderArray<unsigned int>(tree_reader, "EcalEndcapPClusterAssociations.simID");
            recoAssocEcalEndcapN_old = new TTreeReaderArray<unsigned int>(tree_reader, "EcalEndcapNClusterAssociations.recID");
            simuAssocEcalEndcapN_old = new TTreeReaderArray<unsigned int>(tree_reader, "EcalEndcapNClusterAssociations.simID");
        }
    }

    TTreeReaderArray<float> EcalBarrelEng(tree_reader, "EcalBarrelClusters.energy");
    TTreeReaderArray<float> EcalEndcapPEng(tree_reader, "EcalEndcapPClusters.energy");
    TTreeReaderArray<float> EcalEndcapNEng(tree_reader, "EcalEndcapNClusters.energy");
    std::unique_ptr<TTreeReaderArray<float>> EcalEndcapPPosX, EcalEndcapPPosY, EcalEndcapPPosZ;
    if (femc_mode != UpsilonFEMC::Full) {
        EcalEndcapPPosX.reset(new TTreeReaderArray<float>(tree_reader, "EcalEndcapPClusters.position.x"));
        EcalEndcapPPosY.reset(new TTreeReaderArray<float>(tree_reader, "EcalEndcapPClusters.position.y"));
        EcalEndcapPPosZ.reset(new TTreeReaderArray<float>(tree_reader, "EcalEndcapPClusters.position.z"));
    }
    
    TTreeReaderArray<float> RPEng(tree_reader, "ForwardRomanPotRecParticles.energy");
    TTreeReaderArray<float> RPMomX(tree_reader, "ForwardRomanPotRecParticles.momentum.x");
    TTreeReaderArray<float> RPMomY(tree_reader, "ForwardRomanPotRecParticles.momentum.y");
    TTreeReaderArray<float> RPMomZ(tree_reader, "ForwardRomanPotRecParticles.momentum.z");

    TTreeReaderArray<float> OffMEng(tree_reader, "ForwardOffMRecParticles.energy");
    TTreeReaderArray<float> OffMMomX(tree_reader, "ForwardOffMRecParticles.momentum.x");
    TTreeReaderArray<float> OffMMomY(tree_reader, "ForwardOffMRecParticles.momentum.y");
    TTreeReaderArray<float> OffMMomZ(tree_reader, "ForwardOffMRecParticles.momentum.z");

    TTreeReaderArray<float> reco_lowQ2_px    = {tree_reader, "TaggerTrackerReconstructedParticles.momentum.x"};
    TTreeReaderArray<float> reco_lowQ2_py    = {tree_reader, "TaggerTrackerReconstructedParticles.momentum.y"};
    TTreeReaderArray<float> reco_lowQ2_pz    = {tree_reader, "TaggerTrackerReconstructedParticles.momentum.z"}; 

    std::unique_ptr<UpsilonEleID::RecoMatchReader> recoIDReader;
    if (b_ele_id_mode == UpsilonEleID::Reco)
        recoIDReader.reset(new UpsilonEleID::RecoMatchReader(tree_reader, std::move(firstFileCollectionIDs)));
    int lastRecoInputTree = -1;
    Long64_t idMatchedTotal = 0, idMissingTotal = 0, idMaskedTotal = 0, idAmbiguousTotal = 0;

    int nEvents = 0;
    Long64_t femc_total_clusters = 0, femc_removed_clusters = 0, femc_affected_events = 0;
    TLorentzVector k_beam(0, 0, -E_e_beam, E_e_beam); 
    TLorentzVector P_beam(0, 0, E_p_beam, sqrt(E_p_beam*E_p_beam + pMass*pMass)); 

    while(tree_reader.Next()) 
    {
        nEvents++;
        b_femc_nClustersRemoved = 0;
        b_ele_id_nMatched = b_ele_id_nMissing = b_ele_id_nMasked = b_ele_id_nAmbiguous = -1;
        if (recoIDReader && mychain->GetTreeNumber() != lastRecoInputTree) {
            recoIDReader->SetInputTree(mychain->GetTree(), lastRecoInputTree < 0);
            lastRecoInputTree = mychain->GetTreeNumber();
        }
        if (nEvents % 1000 == 0) fprintf(stderr, "\rProcessing event %d...", nEvents);
        
        h_CutFlow->Fill(1);

        // 1. HARD RESET ALL VARIABLES TO -999.0
        b_weight = global_weight;
        b_weight_proj = 0.0;
        b_reco_topology = 0;
        b_true_mass = -999.; b_true_y = -999.; b_true_pT = -999.; 
        b_true_W = -999.; b_true_Q2 = -999.; b_true_t = -999.;
        b_true_ePlus_pT = -999.; b_true_ePlus_eta = -999.;
        b_true_decay_eMinus_pT = -999.; b_true_decay_eMinus_eta = -999.;
        b_true_scat_eMinus_pT = -999.; b_true_scat_eMinus_eta = -999.; b_true_scat_eMinus_E = -999.;
        b_true_nTracks_Central = 0; b_true_nTracks_Tagger = 0;

        b_reco_mass = -999.; b_reco_y = -999.; b_reco_pT = -999.; 
        b_reco_W = -999.; b_reco_Q2 = -999.; 
        b_reco_t_eX = -999.; b_reco_t_BABE = -999.; b_reco_t_eXBABE = -999.; 
        b_reco_t_n = -999.; b_reco_t_x = -999.; b_reco_omega = -999.; b_reco_t_proj = -999.;
        b_reco_ePlus_pT = -999.; b_reco_ePlus_eta = -999.;
        b_reco_decay_eMinus_pT = -999.; b_reco_decay_eMinus_eta = -999.;
        b_reco_scat_eMinus_pT = -999.; b_reco_scat_eMinus_eta = -999.; b_reco_scat_eMinus_E = -999.;
        b_reco_tagged_proton_p = -999.; b_sigma_epz = -999.;
        b_reco_nTracks_Central = 0; b_reco_nTracks_Tagger = 0;
        
        // -------------------------------------------------------------------
        // TRUTH ANALYSIS
        // -------------------------------------------------------------------
        int true_nTracks_Central = 0;
        int true_nTracks_Tagger = 0;
        std::vector<TLorentzVector> trueElectrons, truePositrons;

        for (unsigned int i = 0; i < partGenStat.GetSize(); i++) {
            if (partGenStat[i] != 1) continue; 
            int pdg = std::abs(partPdg[i]);
            bool isCharged = (pdg == 11 || pdg == 13 || pdg == 211 || pdg == 321 || 
                              pdg == 2212 || pdg == 3222 || pdg == 3112 || pdg == 3312 || pdg == 3334);
            if (isCharged) {
                TVector3 trueMom(partMomX[i], partMomY[i], partMomZ[i]);
                if (trueMom.Mag() > 0) {
                    double eta = trueMom.Eta();
                    if (std::abs(eta) < 3.5) true_nTracks_Central++;
                    else if (eta < -4.0) true_nTracks_Tagger++;
                }
            }
            if (partPdg[i] == 11) { 
                TLorentzVector v; v.SetXYZM(partMomX[i], partMomY[i], partMomZ[i], eMass);
                trueElectrons.push_back(v);
            } else if (partPdg[i] == -11) { 
                TLorentzVector v; v.SetXYZM(partMomX[i], partMomY[i], partMomZ[i], eMass);
                truePositrons.push_back(v);
            }
        }

        b_true_nTracks_Central = true_nTracks_Central;
        b_true_nTracks_Tagger = true_nTracks_Tagger;

        bool foundGenEvent = false;
        TLorentzVector v_Gen_Upsilon, v_Gen_ePlus, v_Gen_Decay_eMinus, v_Gen_Scat_eMinus;

        if (trueElectrons.size() >= 2 && truePositrons.size() >= 1) {
            double bestMetric = 9999.0;
            for (size_t i = 0; i < truePositrons.size(); i++) {
                for (size_t j = 0; j < trueElectrons.size(); j++) {
                    TLorentzVector ups = truePositrons[i] + trueElectrons[j];
                    int idx_other = -1;
                    for (size_t k = 0; k < trueElectrons.size(); k++) {
                        if (k != j) { idx_other = k; break; }
                    }
                    if (idx_other == -1) continue;
                    TLorentzVector scat = trueElectrons[idx_other];
                    if (scat.Pz() >= 0.0) continue; 

                    double diff = GetMinUpsilonDiff(ups.M());
                    if (diff < bestMetric) {
                        bestMetric = diff;
                        v_Gen_Upsilon = ups;
                        v_Gen_ePlus = truePositrons[i];
                        v_Gen_Decay_eMinus = trueElectrons[j];
                        v_Gen_Scat_eMinus = scat;
                        foundGenEvent = true;
                    }
                }
            }
        }

        if (foundGenEvent) {
            h_CutFlow->Fill(2);
            if (fabs(v_Gen_ePlus.Eta()) < 3.5 && fabs(v_Gen_Decay_eMinus.Eta()) < 3.5) h_CutFlow->Fill(3);

            TLorentzVector q = k_beam - v_Gen_Scat_eMinus;
            b_true_mass = v_Gen_Upsilon.M();
            b_true_y = v_Gen_Upsilon.Rapidity();
            b_true_pT = v_Gen_Upsilon.Pt();
            b_true_Q2 = -q.M2();
            b_true_W = (P_beam + q).M();
            b_true_t = fabs((q - v_Gen_Upsilon).M2()); 
            b_true_ePlus_pT = v_Gen_ePlus.Pt();
            b_true_ePlus_eta = v_Gen_ePlus.Eta();
            b_true_decay_eMinus_pT = v_Gen_Decay_eMinus.Pt();
            b_true_decay_eMinus_eta = v_Gen_Decay_eMinus.Eta();
            b_true_scat_eMinus_pT = v_Gen_Scat_eMinus.Pt();
            b_true_scat_eMinus_eta = v_Gen_Scat_eMinus.Eta();
            b_true_scat_eMinus_E = v_Gen_Scat_eMinus.E();
        }

        // -------------------------------------------------------------------
        // RECO ANALYSIS
        // -------------------------------------------------------------------
        b_reco_nTracks_Central = trackEng.GetSize();
        b_reco_nTracks_Tagger = reco_lowQ2_pz.GetSize();

        std::vector<double> trackCaloEnergies;
        std::vector<int> recoMatchStatus;
        {
            // Mask before BOTH E/p and isolation. No track or generated-event veto:
            // tracker counts, truth acceptance, and generated normalization stay fixed.
            UpsilonFEMC::ScopedEnergyMask femcMask(EcalEndcapPEng,
                EcalEndcapPPosX.get(), EcalEndcapPPosY.get(), EcalEndcapPPosZ.get(),
                femc_mode, femc_eta_split);
            b_femc_nClustersRemoved = static_cast<int>(femcMask.RemovedCount());
            if (recoIDReader) {
                auto result = recoIDReader->Compute(trackEng.GetSize(), EcalBarrelEng,
                    EcalEndcapPEng, EcalEndcapNEng, femcMask.ActiveMask());
                trackCaloEnergies = std::move(result.energy);
                recoMatchStatus = std::move(result.status);
                b_ele_id_nMatched = result.matched;
                b_ele_id_nMissing = result.missing;
                b_ele_id_nMasked = result.masked;
                b_ele_id_nAmbiguous = result.ambiguous;
                idMatchedTotal += result.matched;
                idMissingTotal += result.missing;
                idMaskedTotal += result.masked;
                idAmbiguousTotal += result.ambiguous;
            } else if (isNewFormat) {
                trackCaloEnergies = PrecomputeTrackEnergies_New(
                    trackEng.GetSize(), *recoAssoc_new, *simuAssoc_new, EcalBarrelEng, EcalEndcapPEng, EcalEndcapNEng,
                    *recoAssocEcalBarrel_new, *simuAssocEcalBarrel_new, *recoAssocEcalEndcapP_new, *simuAssocEcalEndcapP_new, *recoAssocEcalEndcapN_new, *simuAssocEcalEndcapN_new);
            } else {
                trackCaloEnergies = PrecomputeTrackEnergies_Old(
                    trackEng.GetSize(), *recoAssoc_old, *simuAssoc_old, EcalBarrelEng, EcalEndcapPEng, EcalEndcapNEng,
                    *recoAssocEcalBarrel_old, *simuAssocEcalBarrel_old, *recoAssocEcalEndcapP_old, *simuAssocEcalEndcapP_old, *recoAssocEcalEndcapN_old, *simuAssocEcalEndcapN_old);
            }
        } // Restore input cluster energies; trackCaloEnergies retains the selected scenario.
        femc_total_clusters += EcalEndcapPEng.GetSize();
        femc_removed_clusters += b_femc_nClustersRemoved;
        if (b_femc_nClustersRemoved > 0) ++femc_affected_events;

        std::vector<int> good_eMinus, good_ePlus, good_tagger_eMinus;
        for (unsigned int i = 0; i < trackEng.GetSize(); i++) {
            bool isElec = IsElectron_Fast(i, trackMomX, trackMomY, trackMomZ, trackCaloEnergies, false);
            // Missing/ambiguous matches may never enter the reco selection, even
            // if an ID threshold is changed later. Keep the legacy path unchanged.
            if (recoIDReader && recoMatchStatus[i] != UpsilonEleID::Matched) isElec = false;
            if (isElec) {
                if (trackCharge[i] > 0.5) good_ePlus.push_back(i);
                else if (trackCharge[i] < -0.5) good_eMinus.push_back(i);
            }
        }
        for (unsigned int i = 0; i < reco_lowQ2_px.GetSize(); i++) {
            TVector3 taggerMom(reco_lowQ2_px[i], reco_lowQ2_py[i], reco_lowQ2_pz[i]);
            if (reco_lowQ2_pz[i] < 0.0 && taggerMom.Eta() < -4.0) good_tagger_eMinus.push_back(i);
        }

        if (good_ePlus.size() >= 1) h_CutFlow->Fill(4);
        if (good_ePlus.size() >= 1 && good_eMinus.size() >= 1) h_CutFlow->Fill(5);
        bool hasScatElec = (good_eMinus.size() >= 2) || (good_eMinus.size() >= 1 && good_tagger_eMinus.size() >= 1);
        if (good_ePlus.size() >= 1 && hasScatElec) h_CutFlow->Fill(6);

        TLorentzVector v_ePlus, v_Decay_eMinus, v_Scattered_eMinus, v_UpsilonRec;
        double recoMass = -999.0;
        bool foundValidCombination = false;
        double bestMetric = 999.0;
        int reco_scenario = 0;

        if (!good_ePlus.empty() && !good_eMinus.empty()) {
            if (good_eMinus.size() >= 2) {
                for (int idx_plus : good_ePlus) {
                    TLorentzVector cand_ePlus;
                    cand_ePlus.SetXYZM(trackMomX[idx_plus], trackMomY[idx_plus], trackMomZ[idx_plus], eMass);
                    for (size_t i1 = 0; i1 < good_eMinus.size(); i1++) {
                        for (size_t i2 = 0; i2 < good_eMinus.size(); i2++) {
                            if (i1 == i2) continue; 
                            int idx_decay = good_eMinus[i1];
                            int idx_scat  = good_eMinus[i2];
                            if (trackMomZ[idx_scat] >= 0.0) continue; 

                            TLorentzVector cand_eMinus_decay, cand_eMinus_scat;
                            cand_eMinus_decay.SetXYZM(trackMomX[idx_decay], trackMomY[idx_decay], trackMomZ[idx_decay], eMass);
                            cand_eMinus_scat.SetXYZM(trackMomX[idx_scat], trackMomY[idx_scat], trackMomZ[idx_scat], eMass);

                            TLorentzVector cand_ups = cand_ePlus + cand_eMinus_decay;
                            double diff = GetMinUpsilonDiff(cand_ups.M());
                            
                            if (diff < bestMetric) {
                                bestMetric = diff; v_UpsilonRec = cand_ups; v_ePlus = cand_ePlus;
                                v_Decay_eMinus = cand_eMinus_decay; v_Scattered_eMinus = cand_eMinus_scat;
                                recoMass = cand_ups.M(); foundValidCombination = true; reco_scenario = 1;
                            }
                        }
                    }
                }
            }
            if (good_tagger_eMinus.size() >= 1) {
                for (int idx_plus : good_ePlus) {
                    for (int idx_decay : good_eMinus) {
                        TLorentzVector cand_ePlus, cand_eMinus_decay;
                        cand_ePlus.SetXYZM(trackMomX[idx_plus], trackMomY[idx_plus], trackMomZ[idx_plus], eMass);
                        cand_eMinus_decay.SetXYZM(trackMomX[idx_decay], trackMomY[idx_decay], trackMomZ[idx_decay], eMass);
                        
                        TLorentzVector cand_ups = cand_ePlus + cand_eMinus_decay;
                        double diffM = GetMinUpsilonDiff(cand_ups.M());
                        if (diffM < bestMetric) {
                            double ups_EPz = cand_ups.E() - cand_ups.Pz();
                            double expected_P = E_e_beam - (0.5 * ups_EPz);
                            double bestDiffP = 9999.0; int bestTaggerIdx = -1;
                            
                            for(int t_idx : good_tagger_eMinus) {
                                TVector3 taggerMom(reco_lowQ2_px[t_idx], reco_lowQ2_py[t_idx], reco_lowQ2_pz[t_idx]);
                                double diffP = std::abs(taggerMom.Mag() - expected_P);
                                if (diffP < bestDiffP) { bestDiffP = diffP; bestTaggerIdx = t_idx; }
                            }
                            
                            if (bestTaggerIdx != -1) {
                                bestMetric = diffM; v_UpsilonRec = cand_ups; v_ePlus = cand_ePlus;
                                v_Decay_eMinus = cand_eMinus_decay;
                                v_Scattered_eMinus.SetXYZM(reco_lowQ2_px[bestTaggerIdx], reco_lowQ2_py[bestTaggerIdx], reco_lowQ2_pz[bestTaggerIdx], eMass);
                                recoMass = cand_ups.M(); foundValidCombination = true; reco_scenario = 2;
                            }
                        }
                    }
                }
            }
        } // End Combinatorics

        // ===================================================================
        // 2. ASSIGN DATA (SAVE ALL CANDIDATES, REGARDLESS OF CUTS!)
        // ===================================================================
        if (foundValidCombination) {
            double SigmaEPz = (v_Scattered_eMinus.E() - v_Scattered_eMinus.Pz()) + 
                              (v_Decay_eMinus.E() - v_Decay_eMinus.Pz()) + 
                              (v_ePlus.E() - v_ePlus.Pz());
            
            TLorentzVector q_e = k_beam - v_Scattered_eMinus;
            double Q2_e = -q_e.M2();
            double y_e = (P_beam.Dot(q_e)) / (P_beam.Dot(k_beam));
            double W_e = (P_beam + q_e).M();
            double t_eX = (q_e - v_UpsilonRec).M2(); 

            double delta_h = (v_UpsilonRec.E() - v_UpsilonRec.Pz());
            double pt2_h = v_UpsilonRec.Pt() * v_UpsilonRec.Pt();
            double alpha_e = tan(v_Scattered_eMinus.Theta() / 2.0);
            double alpha_h = delta_h / sqrt(pt2_h);
            double y_DA = alpha_h / (alpha_e + alpha_h);
            double Q2_DA = (4.0 * E_e_beam * E_e_beam) / (alpha_e * (alpha_e + alpha_h));
            double x_DA = Q2_DA / (4.0 * E_e_beam * E_p_beam * y_DA);
            double W_DA = sqrt(pMass*pMass + Q2_DA * (1.0/x_DA - 1.0));

            double W_Best, Q2_Best;
            if (y_e > 0.045) { W_Best = W_e; Q2_Best = Q2_e; } 
            else { W_Best = W_DA; Q2_Best = Q2_DA; }

            // Proton Tagging
            TLorentzVector v_Scattered_p;
            bool foundProton = false;
            for (unsigned int i = 0; i < RPEng.GetSize(); i++) {
                if (RPEng[i] > E_p_beam * 0.5) { 
                    double pz = RPMomZ[i] < 0. ? -RPMomZ[i] : RPMomZ[i]; 
                    v_Scattered_p.SetXYZM(RPMomX[i], RPMomY[i], pz, pMass);
                    foundProton = true; break;
                }
            }
            if (!foundProton) {
                for (unsigned int i = 0; i < OffMEng.GetSize(); i++) {
                    if (OffMEng[i] > E_p_beam * 0.5) {
                        double pz = OffMMomZ[i] < 0. ? -OffMMomZ[i] : OffMMomZ[i];
                        v_Scattered_p.SetXYZM(OffMMomX[i], OffMMomY[i], pz, pMass);
                        foundProton = true; break;
                    }
                }
            }

            // Assign Branches
            b_reco_mass = recoMass;
            b_reco_y = v_UpsilonRec.Rapidity();
            b_reco_pT = v_UpsilonRec.Pt();
            b_reco_W = W_Best;
            b_reco_Q2 = Q2_Best;
            b_sigma_epz = SigmaEPz;
            b_reco_t_eX = fabs(t_eX);
            b_reco_ePlus_pT = v_ePlus.Pt();
            b_reco_ePlus_eta = v_ePlus.Eta();
            b_reco_decay_eMinus_pT = v_Decay_eMinus.Pt();
            b_reco_decay_eMinus_eta = v_Decay_eMinus.Eta();
            b_reco_scat_eMinus_pT = v_Scattered_eMinus.Pt();
            b_reco_scat_eMinus_eta = v_Scattered_eMinus.Eta();
            b_reco_scat_eMinus_E = v_Scattered_eMinus.E();

            if (foundProton) {
                b_reco_tagged_proton_p = v_Scattered_p.P();
                b_reco_t_BABE = fabs(-(P_beam - v_Scattered_p).M2());
                TLorentzVector P_Miss = k_beam + P_beam - v_Scattered_eMinus - v_UpsilonRec;
                double p_mag = P_Miss.P();
                TLorentzVector v_Scat_p_corr;
                v_Scat_p_corr.SetXYZM(
                    p_mag * sin(v_Scattered_p.Theta()) * cos(v_Scattered_p.Phi()),
                    p_mag * sin(v_Scattered_p.Theta()) * sin(v_Scattered_p.Phi()),
                    p_mag * cos(v_Scattered_p.Theta()), pMass
                );
                b_reco_t_eXBABE = fabs(-(P_beam - v_Scat_p_corr).M2());
            }
            
            if (reco_scenario == 1) b_reco_topology = foundProton ? 1 : 3;
            else if (reco_scenario == 2) b_reco_topology = foundProton ? 2 : 4;

            // Projective t
            TVector3 p3_k = k_beam.Vect();
            TVector3 p3_kp = v_Scattered_eMinus.Vect();
            TVector3 n_hat = (p3_k.Cross(p3_kp)).Unit(); 
            TVector3 p3_VM = v_UpsilonRec.Vect();
            b_reco_t_n = std::pow(p3_VM.Dot(n_hat), 2); 
            
            double abs_t = foundProton ? b_reco_t_eXBABE : b_reco_t_eX;
            b_reco_t_x = (abs_t > b_reco_t_n) ? (abs_t - b_reco_t_n) : 0.0;
            b_reco_omega = std::atan2(std::sqrt(b_reco_t_x), std::sqrt(b_reco_t_n));
            
            double omega_max = TMath::Pi() / 12.0;
            if (b_reco_omega < omega_max) {
                b_reco_t_proj = abs_t;
                b_weight_proj = global_weight * ((TMath::Pi() / 2.0) / omega_max); 
            }

            // ===================================================================
            // 3. BOOKKEEPING (STRICT CASCADE CUTS)
            // ===================================================================
            if (recoMass >= 8.0 && recoMass <= 12.0) {
                h_CutFlow->Fill(7);
                const double sigmaCenter = 2.0 * E_e_beam;
                const double sigmaHalfWidth = 0.2 * E_e_beam;
                
                if (std::abs(SigmaEPz - sigmaCenter) <= sigmaHalfWidth) {
                    h_CutFlow->Fill(8);
                    if (W_Best >= 0.0) { //20.0!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        h_CutFlow->Fill(9);
                        if (foundProton) h_CutFlow->Fill(10);
                    }
                }
            }
        } // End valid candidate checks

        // THE GOLDEN RULE: We output EVERY event to the TTree exactly once.
        outTree->Fill();
    } 

    if (nEvents != total_entries) {
        throw std::runtime_error("Event reading stopped before all input events were processed. Check branch/schema errors; output is incomplete.");
    }
    std::cout << "\nAnalysis Finished." << std::endl;
    std::cout << "fEMC clusters masked: " << femc_removed_clusters << " / " << femc_total_clusters
              << "; affected events: " << femc_affected_events << " / " << nEvents << std::endl;

    if (recoIDReader)
        std::cout << "Reco-ID match counts: usable=" << idMatchedTotal << ", missing=" << idMissingTotal
                  << ", masked-only=" << idMaskedTotal << ", ambiguous=" << idAmbiguousTotal << std::endl;

    ofile->cd();
    // Scenario details are in UpsilonTree; aggregate diagnostics are printed above.
    outTree->Write();
    h_EventCount->Write();
    h_CutFlow->Write(); 
    if (gDirectory->Get("h_EleID_EpRatio")) gDirectory->Get("h_EleID_EpRatio")->Write();
    if (gDirectory->Get("h_EleID_IsoFraction")) gDirectory->Get("h_EleID_IsoFraction")->Write();
    if (gDirectory->Get("h_EleID_CutFlow")) gDirectory->Get("h_EleID_CutFlow")->Write();
    ofile->Close();
    std::cout << "Flat Tree written to: " << outfilename << std::endl;
}
