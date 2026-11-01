#ifndef UPSILON_ELECTRON_FINDER_C
#define UPSILON_ELECTRON_FINDER_C
// ele_id_mode="legacy": original truth-assisted energy association, unchanged.
// ele_id_mode="reco": reconstructed TrackClusterMatches, no MC association fallback.
// The fEMC geometry mask remains in DVMP_Upsilon_Analysis.C and is applied first.

// This block has no ROOT or MC dependency. It is shared by the production reader
// below and the standalone regression tests. Inputs are reconstructed relations.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace UpsilonEleID {
enum Mode { Legacy = 0, Reco = 1 };
enum MatchStatus {
    Matched = 0, NoTrack = 1, NoMatch = 2, MaskedOnly = 3,
    MultipleTracks = 4, SharedTrack = 5, MultipleClusters = 6, SharedCluster = 7
};
struct ObjectKey {
    std::uint32_t collection = 0;
    int index = -1;
    bool operator<(const ObjectKey& other) const {
        return std::tie(collection, index) < std::tie(other.collection, other.index);
    }
    bool operator==(const ObjectKey& other) const {
        return collection == other.collection && index == other.index;
    }
};
struct RecoCluster {
    ObjectKey key;
    double energy = 0.0;
    bool active = true; // Already determined by the existing fEMC mask in DVMP.
};
struct RecoLink { ObjectKey track, cluster; };
struct RecoEnergyResult {
    std::vector<double> energy;       // Indexed by ReconstructedChargedParticles row.
    std::vector<int> status;          // MatchStatus, not an electron/PID decision.
    int matched = 0, missing = 0, masked = 0, ambiguous = 0;
};

inline RecoEnergyResult BuildRecoEnergyMap(
    const std::vector<std::vector<ObjectKey>>& particleTracks,
    const std::set<ObjectKey>& availableTracks,
    const std::vector<RecoCluster>& clusters,
    const std::vector<RecoLink>& links)
{
    // Validate complete ObjectIDs. Do not equate a particle row with a track index.
    std::map<ObjectKey, RecoCluster> clusterByID;
    for (const auto& c : clusters) {
        if (c.key.index < 0 || !std::isfinite(c.energy) || c.energy < 0.0 ||
            !clusterByID.emplace(c.key, c).second)
            throw std::runtime_error("Reco electron ID: invalid/duplicate cluster or energy.");
    }
    for (const auto& t : availableTracks)
        if (t.index < 0) throw std::runtime_error("Reco electron ID: invalid track identity.");

    std::vector<std::set<ObjectKey>> distinctTracks(particleTracks.size());
    std::map<ObjectKey, std::set<std::size_t>> particleOwners;
    for (std::size_t i = 0; i < particleTracks.size(); ++i) {
        for (const auto& t : particleTracks[i]) {
            if (!availableTracks.count(t))
                throw std::runtime_error("Reco electron ID: unresolved particle-to-track relation.");
            distinctTracks[i].insert(t);
            particleOwners[t].insert(i);
        }
    }

    std::map<ObjectKey, std::set<ObjectKey>> activeClusters, clusterOwners;
    std::set<ObjectKey> maskedTracks;
    for (const auto& link : links) {
        if (!availableTracks.count(link.track) || !clusterByID.count(link.cluster))
            throw std::runtime_error("Reco electron ID: unresolved track-cluster relation; no truth/index-only fallback.");
        if (!clusterByID.at(link.cluster).active) {
            maskedTracks.insert(link.track);
            continue;
        }
        // Repeated copies of the identical link do not duplicate cluster energy.
        activeClusters[link.track].insert(link.cluster);
        clusterOwners[link.cluster].insert(link.track);
    }

    RecoEnergyResult result;
    result.energy.assign(particleTracks.size(), 0.0);
    result.status.assign(particleTracks.size(), NoTrack);
    for (std::size_t i = 0; i < particleTracks.size(); ++i) {
        const auto& tracks = distinctTracks[i];
        int status = NoTrack;
        if (tracks.size() > 1) status = MultipleTracks;
        else if (tracks.size() == 1) {
            const auto& track = *tracks.begin();
            const auto& candidates = activeClusters[track];
            if (particleOwners.at(track).size() > 1) status = SharedTrack;
            else if (candidates.empty()) status = maskedTracks.count(track) ? MaskedOnly : NoMatch;
            else if (candidates.size() > 1) status = MultipleClusters;
            else {
                const auto& cluster = *candidates.begin();
                if (clusterOwners.at(cluster).size() > 1) status = SharedCluster;
                else {
                    status = Matched;
                    result.energy[i] = clusterByID.at(cluster).energy;
                }
            }
        }
        result.status[i] = status;
        if (status == Matched) ++result.matched;
        else if (status == NoTrack || status == NoMatch) ++result.missing;
        else if (status == MaskedOnly) ++result.masked;
        else ++result.ambiguous;
    }
    return result;
}
} // namespace UpsilonEleID

#ifndef UPSILON_ELEID_CORE_ONLY
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <numeric>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <vector>
#include <TFile.h>
#include <TH1D.h>
#include <TVector3.h>
#include <TTreeReaderArray.h>
#include <podio/ObjectID.h>

static TH1D* h_EleID_EpRatio = nullptr;
static TH1D* h_EleID_IsoFraction = nullptr;
static TH1D* h_EleID_CutFlow = nullptr;
static bool histograms_initialized = false;

// --- ID Configuration Parameters ---
const double mEoP_min = 0.8; 
const double mEoP_max = 1.2;
const double mIsoR = 0.4;  
const double mIsoE = 0.9;  

// Histogram objects are owned by the current output TFile. Reset pointers when
// starting a new analysis call; do not delete objects owned by another file.
void ResetElectronIDHistograms() {
    h_EleID_EpRatio = nullptr;
    h_EleID_IsoFraction = nullptr;
    h_EleID_CutFlow = nullptr;
    histograms_initialized = false;
}

void InitElectronIDHistograms() {
    if (histograms_initialized) return;
    if (gFile) {
        gFile->cd();
        h_EleID_EpRatio = new TH1D("h_EleID_EpRatio", "Electron ID: E/p Ratio;E_{cal}/p_{track};Entries", 100, 0.0, 2.0);
        h_EleID_IsoFraction = new TH1D("h_EleID_IsoFraction", "Electron ID: Isolation Fraction;E_{cand} / E_{cone};Entries", 100, 0.0, 1.1);
        h_EleID_CutFlow = new TH1D("h_EleID_CutFlow", "Electron ID Cut Flow", 5, 0, 5);
        h_EleID_CutFlow->GetXaxis()->SetBinLabel(1, "All Candidates");
        h_EleID_CutFlow->GetXaxis()->SetBinLabel(2, "Pass Pz < 0 (If Req)");
        h_EleID_CutFlow->GetXaxis()->SetBinLabel(3, "Pass 0.8 < E/p < 1.2");
        h_EleID_CutFlow->GetXaxis()->SetBinLabel(4, "Pass Iso > 0.9");
        h_EleID_CutFlow->GetXaxis()->SetBinLabel(5, "ACCEPTED");
        histograms_initialized = true;
    }
}

// ===================================================================================
// FAST PRE-COMPUTATION (Hash Maps)
// Builds a vector where index = TrackID, and value = Total ECal Energy
// ===================================================================================

// For OLD Format (unsigned int)
std::vector<double> PrecomputeTrackEnergies_Old(
    int numTracks,
    TTreeReaderArray<unsigned int>& recoAssoc, TTreeReaderArray<unsigned int>& simuAssoc,
    TTreeReaderArray<float>& EcalBarrelEng, TTreeReaderArray<float>& EcalEndcapPEng, TTreeReaderArray<float>& EcalEndcapNEng, 
    TTreeReaderArray<unsigned int>& recoAssocEcalBarrel, TTreeReaderArray<unsigned int>& simuAssocEcalBarrel, 
    TTreeReaderArray<unsigned int>& recoAssocEcalEndcapP, TTreeReaderArray<unsigned int>& simuAssocEcalEndcapP, 
    TTreeReaderArray<unsigned int>& recoAssocEcalEndcapN, TTreeReaderArray<unsigned int>& simuAssocEcalEndcapN) 
{
    std::unordered_map<int, double> simIdToEnergy;

    // 1. Build Dictionary: simID -> Total ECal Energy (Loops only once!)
    for (unsigned int j = 0; j < simuAssocEcalBarrel.GetSize(); j++) {
        if (recoAssocEcalBarrel[j] < EcalBarrelEng.GetSize()) simIdToEnergy[simuAssocEcalBarrel[j]] += EcalBarrelEng[recoAssocEcalBarrel[j]];
    }
    for (unsigned int j = 0; j < simuAssocEcalEndcapP.GetSize(); j++) {
        if (recoAssocEcalEndcapP[j] < EcalEndcapPEng.GetSize()) simIdToEnergy[simuAssocEcalEndcapP[j]] += EcalEndcapPEng[recoAssocEcalEndcapP[j]];
    }
    for (unsigned int j = 0; j < simuAssocEcalEndcapN.GetSize(); j++) {
        if (recoAssocEcalEndcapN[j] < EcalEndcapNEng.GetSize()) simIdToEnergy[simuAssocEcalEndcapN[j]] += EcalEndcapNEng[recoAssocEcalEndcapN[j]];
    }

    // 2. Map Track Index -> simID -> Energy
    std::vector<double> trackEnergies(numTracks, 0.0);
    for (unsigned int j = 0; j < recoAssoc.GetSize(); j++) {
        int trackIdx = recoAssoc[j];
        if (trackIdx >= 0 && trackIdx < numTracks) trackEnergies[trackIdx] = simIdToEnergy[simuAssoc[j]];
    }

    return trackEnergies;
}

// For NEW Format (podio::ObjectID)
std::vector<double> PrecomputeTrackEnergies_New(
    int numTracks,
    TTreeReaderArray<podio::ObjectID>& recoAssoc, TTreeReaderArray<podio::ObjectID>& simuAssoc,
    TTreeReaderArray<float>& EcalBarrelEng, TTreeReaderArray<float>& EcalEndcapPEng, TTreeReaderArray<float>& EcalEndcapNEng, 
    TTreeReaderArray<podio::ObjectID>& recoAssocEcalBarrel, TTreeReaderArray<podio::ObjectID>& simuAssocEcalBarrel, 
    TTreeReaderArray<podio::ObjectID>& recoAssocEcalEndcapP, TTreeReaderArray<podio::ObjectID>& simuAssocEcalEndcapP, 
    TTreeReaderArray<podio::ObjectID>& recoAssocEcalEndcapN, TTreeReaderArray<podio::ObjectID>& simuAssocEcalEndcapN) 
{
    std::unordered_map<int, double> simIdToEnergy;

    // 1. Build Dictionary: simID -> Total ECal Energy (Loops only once!)
    for (unsigned int j = 0; j < simuAssocEcalBarrel.GetSize(); j++) {
        if (recoAssocEcalBarrel[j].index < EcalBarrelEng.GetSize()) simIdToEnergy[simuAssocEcalBarrel[j].index] += EcalBarrelEng[recoAssocEcalBarrel[j].index];
    }
    for (unsigned int j = 0; j < simuAssocEcalEndcapP.GetSize(); j++) {
        if (recoAssocEcalEndcapP[j].index < EcalEndcapPEng.GetSize()) simIdToEnergy[simuAssocEcalEndcapP[j].index] += EcalEndcapPEng[recoAssocEcalEndcapP[j].index];
    }
    for (unsigned int j = 0; j < simuAssocEcalEndcapN.GetSize(); j++) {
        if (recoAssocEcalEndcapN[j].index < EcalEndcapNEng.GetSize()) simIdToEnergy[simuAssocEcalEndcapN[j].index] += EcalEndcapNEng[recoAssocEcalEndcapN[j].index];
    }

    // 2. Map Track Index -> simID -> Energy
    std::vector<double> trackEnergies(numTracks, 0.0);
    for (unsigned int j = 0; j < recoAssoc.GetSize(); j++) {
        int trackIdx = recoAssoc[j].index;
        if (trackIdx >= 0 && trackIdx < numTracks) trackEnergies[trackIdx] = simIdToEnergy[simuAssoc[j].index];
    }

    return trackEnergies;
}

// ===================================================================================
// MAIN ELECTRON ID (O(1) lookups using precomputed energies!)
// ===================================================================================
bool IsElectron_Fast(int candIndex, 
            TTreeReaderArray<float>& trackMomX, TTreeReaderArray<float>& trackMomY, TTreeReaderArray<float>& trackMomZ,
            const std::vector<double>& trackCaloEnergies, // Directly pass the fast map
            bool requireScattered = false)
{
    if (!histograms_initialized) InitElectronIDHistograms();

    TVector3 candMom(trackMomX[candIndex], trackMomY[candIndex], trackMomZ[candIndex]);
    
    // O(1) Fetch!
    double candCaloE = trackCaloEnergies[candIndex];
                                            
    double p = candMom.Mag();
    double EpRatio = (p > 0.0) ? (candCaloE / p) : 0.0;

    if (h_EleID_CutFlow) h_EleID_CutFlow->Fill("All Candidates", 1);
    if (h_EleID_EpRatio) h_EleID_EpRatio->Fill(EpRatio);

    bool passPz = requireScattered ? (candMom.Z() < 0.0) : true;
    if (passPz && h_EleID_CutFlow) h_EleID_CutFlow->Fill("Pass Pz < 0 (If Req)", 1);

    bool passEp = (EpRatio >= mEoP_min && EpRatio <= mEoP_max);
    if (passPz && passEp && h_EleID_CutFlow) h_EleID_CutFlow->Fill("Pass 0.8 < E/p < 1.2", 1);

    double coneEnergy = candCaloE; 
    for (int i = 0; i < trackMomX.GetSize(); i++) {
        if (i == candIndex) continue; 
        TVector3 otherMom(trackMomX[i], trackMomY[i], trackMomZ[i]);
        
        if (candMom.DeltaR(otherMom) < mIsoR) {
            // O(1) Fetch! No array scanning!
            coneEnergy += trackCaloEnergies[i]; 
        }
    }

    double isoFraction = (coneEnergy > 0.0) ? (candCaloE / coneEnergy) : 1.0;
    if (h_EleID_IsoFraction) h_EleID_IsoFraction->Fill(isoFraction);

    bool passIso = (isoFraction >= mIsoE);
    if (passPz && passEp && passIso && h_EleID_CutFlow) h_EleID_CutFlow->Fill("Pass Iso > 0.9", 1);

    bool accepted = (passPz && passEp && passIso);
    if (accepted && h_EleID_CutFlow) h_EleID_CutFlow->Fill("ACCEPTED", 1);

    return accepted;
}

// ROOT adapter for the reconstructed match collections in branches(1).txt.
// All particle kinematics remain those used by the existing DVMP macro.
// MCParticles are still used by DVMP for truth plots, never by this ID path.
#include <array>
#include <limits>
#include <memory>
#include <TChain.h>
#include <TDirectory.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TSystem.h>
#include <type_traits>
#include <utility>
#include <podio/CollectionIDTable.h>
// This header supplies the real on-disk metadata element type in recent PODIO.
// Do not copy its layout into the analysis or read its string member as an array.
#if __has_include(<podio/utilities/RootHelpers.h>)
#include <podio/utilities/RootHelpers.h>
#endif
// An incomplete declaration allows the old idTable path to keep compiling with
// older headers. The template below only instantiates the vector reader when
// the installed PODIO actually defines the named metadata record.
namespace podio { namespace root_utils { struct CollectionWriteInfo; } }

namespace UpsilonEleID {
inline ObjectKey Key(const podio::ObjectID& id) {
    return {static_cast<std::uint32_t>(id.collectionID), id.index};
}

inline void RequireRecoBranches(TTree* tree) {
    if (!tree) throw std::runtime_error("Reco electron ID: no current events tree.");
    for (const char* name : {
        "ReconstructedChargedParticles", "_ReconstructedChargedParticles_tracks", "CentralCKFTracks",
        "EcalBarrelClusters", "EcalEndcapNClusters", "EcalEndcapPClusters",
        "_EcalBarrelTrackClusterMatches_track", "_EcalBarrelTrackClusterMatches_cluster",
        "_EcalEndcapNTrackClusterMatches_track", "_EcalEndcapNTrackClusterMatches_cluster",
        "_EcalEndcapPTrackClusterMatches_track", "_EcalEndcapPTrackClusterMatches_cluster"}) {
        if (!tree->GetBranch(name))
            throw std::runtime_error(std::string("Reco electron ID: missing required branch ") + name +
                ". Use compatible reconstruction output; there is no automatic legacy fallback.");
    }
}

using CollectionIDs = std::map<std::string, std::uint32_t>;

namespace MetadataDetail {
template<class T, class = void> struct HasNamedCollectionInfo : std::false_type {};
template<class T> struct HasNamedCollectionInfo<T, std::void_t<
    decltype(std::declval<T>().name), decltype(std::declval<T>().collectionID)>>
    : std::true_type {};

// Read the complete vector branch. TTreeReaderArray<std::string> on the split
// ".name" member fails with TStreamerSTLstring on affected ROOT versions even
// when the input schema is valid. PODIO itself reads complete metadata records.
template<class Info, class AcceptTable>
Long64_t ReadNamedCollectionInfo(TTree* metadata, const AcceptTable& acceptTable) {
    if constexpr (HasNamedCollectionInfo<Info>::value) {
        // CollectionWriteInfo and vector<CollectionWriteInfo> are in PODIO's
        // ROOT-I/O dictionary, not just the core libpodioDict library. This
        // load is only requested for files using this metadata representation.
        if (!gSystem || gSystem->Load("libpodioRootIODict") < 0)
            throw std::runtime_error(
                "Reco electron ID: cannot load libpodioRootIODict for the complete "
                "CollectionTypeInfo vector. Use matching PODIO headers/libraries "
                "from the production environment; no collection-ID fallback is used.");

        TTreeReader reader(metadata);
        TTreeReaderValue<std::vector<Info>> records(reader, "events___CollectionTypeInfo");
        Long64_t count = 0;
        while (reader.Next()) {
            if (records.GetSetupStatus() < 0)
                throw std::runtime_error(
                    "Reco electron ID: cannot set up the complete CollectionTypeInfo "
                    "vector reader (ROOT setup status " +
                    std::to_string(records.GetSetupStatus()) + ").");
            const auto* values = records.Get();
            if (!values)
                throw std::runtime_error(
                    "Reco electron ID: failed to read CollectionTypeInfo metadata entry " +
                    std::to_string(count) + ".");
            std::vector<std::string> names;
            std::vector<std::uint32_t> ids;
            names.reserve(values->size());
            ids.reserve(values->size());
            for (const auto& info : *values) {
                names.push_back(info.name);
                ids.push_back(info.collectionID);
            }
            acceptTable(names, ids);
            ++count;
        }
        if (count != metadata->GetEntries())
            throw std::runtime_error(
                "Reco electron ID: incomplete complete-vector metadata read (" +
                std::to_string(count) + "/" + std::to_string(metadata->GetEntries()) +
                " entries; ROOT setup status " + std::to_string(records.GetSetupStatus()) +
                "). Check the installed PODIO ROOT-I/O dictionary and input file.");
        return count;
    } else {
        throw std::runtime_error(
            "Reco electron ID: this file uses named CollectionTypeInfo records, but "
            "the installed PODIO headers do not define them. Use PODIO headers and "
            "libraries compatible with this production. The older idTable path "
            "and legacy electron ID remain available for compatible input files.");
    }
}
} // namespace MetadataDetail
inline CollectionIDs RequiredCollectionIDs(TFile* input) {
    if (!input) throw std::runtime_error("Reco electron ID: no input ROOT file.");
    TDirectory::TContext directoryContext; // Do not change the active output directory.
    auto* metadata = dynamic_cast<TTree*>(input->Get("podio_metadata"));
    if (!metadata || metadata->GetEntries() < 1)
        throw std::runtime_error("Reco electron ID requires PODIO collection-ID metadata. Do not guess collection IDs.");

    CollectionIDs expected;
    const auto acceptTable = [&](const std::vector<std::string>& names,
                                 const std::vector<std::uint32_t>& ids) {
        if (names.size() != ids.size())
            throw std::runtime_error("Reco electron ID: malformed collection-ID metadata.");
        CollectionIDs all;
        std::map<std::uint32_t, std::string> reverse;
        for (std::size_t i = 0; i < names.size(); ++i) {
            if (!all.emplace(names[i], ids[i]).second || !reverse.emplace(ids[i], names[i]).second)
                throw std::runtime_error("Reco electron ID: duplicate name/ID in metadata.");
        }
        CollectionIDs selected;
        for (const char* name : {"CentralCKFTracks", "EcalBarrelClusters", "EcalEndcapNClusters", "EcalEndcapPClusters"}) {
            const auto it = all.find(name);
            if (it == all.end()) throw std::runtime_error(std::string("Reco electron ID: collection absent from metadata: ") + name);
            selected.emplace(*it);
        }
        if (expected.empty()) expected = std::move(selected);
        else if (expected != selected)
            throw std::runtime_error("Reco electron ID: inconsistent collection IDs within one ROOT file (unsafe merged metadata).");
    };

    // Keep the older CollectionIDTable path. For named CollectionWriteInfo,
    // read whole records through the installed PODIO type and ROOT-I/O dictionary.
    Long64_t tablesRead = 0;
    if (metadata->GetBranch("events___idTable")) {
        TTreeReader reader(metadata);
        TTreeReaderValue<podio::CollectionIDTable> table(reader, "events___idTable");
        while (reader.Next()) {
            if (table.GetSetupStatus() < 0)
                throw std::runtime_error("Reco electron ID: cannot read events___idTable; check PODIO dictionaries.");
            acceptTable(table->names(), table->ids());
            ++tablesRead;
        }
    } else if (metadata->GetBranch("events___CollectionTypeInfo")) {
        tablesRead = MetadataDetail::ReadNamedCollectionInfo<
            podio::root_utils::CollectionWriteInfo>(metadata, acceptTable);
    } else {
        throw std::runtime_error("Reco electron ID: unsupported PODIO metadata layout; refusing index-only matching.");
    }
    if (tablesRead != metadata->GetEntries() || expected.empty())
        throw std::runtime_error("Reco electron ID: incomplete collection-ID metadata read.");
    return expected;
}

class RecoMatchReader {
public:
    RecoMatchReader(TTreeReader& reader, CollectionIDs firstFileIDs)
        : ids_(std::move(firstFileIDs)),
          trackBegin_(reader, "ReconstructedChargedParticles.tracks_begin"),
          trackEnd_(reader, "ReconstructedChargedParticles.tracks_end"),
          particleTracks_(reader, "_ReconstructedChargedParticles_tracks"),
          centralTrackPx_(reader, "CentralCKFTracks.momentum.x"),
          barrelTrack_(reader, "_EcalBarrelTrackClusterMatches_track"),
          barrelCluster_(reader, "_EcalBarrelTrackClusterMatches_cluster"),
          backwardTrack_(reader, "_EcalEndcapNTrackClusterMatches_track"),
          backwardCluster_(reader, "_EcalEndcapNTrackClusterMatches_cluster"),
          forwardTrack_(reader, "_EcalEndcapPTrackClusterMatches_track"),
          forwardCluster_(reader, "_EcalEndcapPTrackClusterMatches_cluster") {}

    // Call when the TChain moves to a new input file. No cached numeric ID is
    // allowed to leak across files, and truth-seeded track collections are not substituted.
    void SetInputTree(TTree* tree, bool firstFile = false) {
        RequireRecoBranches(tree);
        if (!firstFile) ids_ = RequiredCollectionIDs(tree->GetCurrentFile());
    }

    RecoEnergyResult Compute(std::size_t nParticles,
                             TTreeReaderArray<float>& barrelEnergy,
                             TTreeReaderArray<float>& forwardEnergy,
                             TTreeReaderArray<float>& backwardEnergy,
                             const std::vector<unsigned char>& forwardActive) {
        Check(trackBegin_); Check(trackEnd_); Check(particleTracks_); Check(centralTrackPx_);
        Check(barrelTrack_); Check(barrelCluster_); Check(backwardTrack_); Check(backwardCluster_);
        Check(forwardTrack_); Check(forwardCluster_);
        Check(barrelEnergy); Check(forwardEnergy); Check(backwardEnergy);
        if (trackBegin_.GetSize() != nParticles || trackEnd_.GetSize() != nParticles ||
            forwardActive.size() != forwardEnergy.GetSize())
            throw std::runtime_error("Reco electron ID: inconsistent particle relation or fEMC-mask array size.");

        const auto trackID = ids_.at("CentralCKFTracks");
        const auto nTracks = centralTrackPx_.GetSize();
        if (nTracks > static_cast<std::size_t>(std::numeric_limits<int>::max()))
            throw std::runtime_error("Reco electron ID: too many tracks for PODIO ObjectID.");
        std::set<ObjectKey> tracks;
        for (std::size_t i = 0; i < nTracks; ++i) tracks.insert({trackID, static_cast<int>(i)});
        std::vector<std::vector<ObjectKey>> particleRefs(nParticles);
        for (std::size_t i = 0; i < nParticles; ++i) {
            const auto begin = trackBegin_[i], end = trackEnd_[i];
            if (begin > end || end > particleTracks_.GetSize())
                throw std::runtime_error("Reco electron ID: invalid tracks_begin/tracks_end.");
            for (std::size_t j = begin; j < end; ++j) {
                const auto key = Key(particleTracks_[j]);
                if (key.collection != trackID || !tracks.count(key))
                    throw std::runtime_error("Reco electron ID: particle track is not a valid CentralCKFTracks object.");
                particleRefs[i].push_back(key);
            }
        }

        std::vector<RecoCluster> clusters;
        std::vector<RecoLink> links;
        const auto add = [&](const char* collectionName, TTreeReaderArray<float>& energy,
                             TTreeReaderArray<podio::ObjectID>& tr,
                             TTreeReaderArray<podio::ObjectID>& cl,
                             const std::vector<unsigned char>* active) {
            if (tr.GetSize() != cl.GetSize())
                throw std::runtime_error("Reco electron ID: unequal track/cluster relation arrays.");
            const auto collectionID = ids_.at(collectionName);
            if (energy.GetSize() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
                throw std::runtime_error("Reco electron ID: too many clusters for PODIO ObjectID.");
            for (std::size_t i = 0; i < energy.GetSize(); ++i)
                clusters.push_back({{collectionID, static_cast<int>(i)}, energy[i], !active || (*active)[i] != 0});
            for (std::size_t i = 0; i < tr.GetSize(); ++i) {
                const auto track = Key(tr[i]), cluster = Key(cl[i]);
                if (track.collection != trackID || !tracks.count(track) ||
                    cluster.collection != collectionID || cluster.index < 0 ||
                    static_cast<std::size_t>(cluster.index) >= energy.GetSize())
                    throw std::runtime_error(std::string("Reco electron ID: invalid ObjectID in ") + collectionName + " matches.");
                links.push_back({track, cluster});
            }
        };
        add("EcalBarrelClusters", barrelEnergy, barrelTrack_, barrelCluster_, nullptr);
        add("EcalEndcapNClusters", backwardEnergy, backwardTrack_, backwardCluster_, nullptr);
        add("EcalEndcapPClusters", forwardEnergy, forwardTrack_, forwardCluster_, &forwardActive);
        return BuildRecoEnergyMap(particleRefs, tracks, clusters, links);
    }
private:
    template<class T> static void Check(const TTreeReaderArray<T>& a) {
        if (a.GetSetupStatus() < 0)
            throw std::runtime_error(std::string("Reco electron ID: unreadable input ") + a.GetBranchName());
    }
    CollectionIDs ids_;
    TTreeReaderArray<unsigned int> trackBegin_, trackEnd_;
    TTreeReaderArray<podio::ObjectID> particleTracks_;
    TTreeReaderArray<float> centralTrackPx_;
    TTreeReaderArray<podio::ObjectID> barrelTrack_, barrelCluster_, backwardTrack_, backwardCluster_, forwardTrack_, forwardCluster_;
};
} // namespace UpsilonEleID

#endif // UPSILON_ELEID_CORE_ONLY
#endif // UPSILON_ELECTRON_FINDER_C
