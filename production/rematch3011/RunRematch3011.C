// Batch entry point for the downstream EICrecon #3011 correction.
// jobs.py snapshots the beam-specific producer next to this file.
#define UPSILON_REMATCH_3011 1
#include "DVMP_Upsilon_Analysis.C"

void RunRematch3011(TString input, std::string output, TString data_version,
                    TString femc_scenario, double femc_eta_split) {
    std::cout << "Running reco-rematch3011. The archived producer prints its historical "
              << "stored-match description below. This run uses corrected links for ID "
              << "and checks the original algorithm against the stored links." << std::endl;
    DVMP_Upsilon_Analysis(input, output, data_version, femc_scenario, femc_eta_split, "reco");
    std::cout << "Completed reco-rematch3011 using corrected ECAL links for electron ID." << std::endl;
}
