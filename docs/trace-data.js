/* One source map for both reading directions. Source links identify the checked
   implementation. Run files are relative to a local replay or state workspace. */
(() => {
  const referenceRevision = "301c1450a3d2cff7af8f7e155b7a06f0790143ae";
  const git = "https://github.com/gkrintir/ePIC-Upsilon/";
  const src = (label, path) => ({label, url:git+"blob/"+referenceRevision+"/"+path});
  const guide = src;
  const page = (label, url) => ({label,url});
  const macro = "production/macros/DVMP_Upsilon_Analysis.C";
  const optimize = "analysis/scripts/reoptimization/";
  const plot = "plotting/upsilon_esr/";
  const sources = {
    upstream: page("Upstream detector reconstruction · EICrecon", "https://github.com/eic/EICrecon"),
    cards: page("July simulation cards · pinned source", "https://github.com/eic/Upsilon_Sim_Cards/tree/f62965a1612fc734dbab18e141f4809740d6b234"),
    production: guide("Run a new tuple production", "docs/production.md"),
    gap: guide("Which CERN production records are missing?", "docs/missing-production-sources.md"),
    entry: src("Tuple macro · entry point and ID mode",macro+"#L166-L188"),
    id: src("ElectronFinder · electron identification","production/macros/ElectronFinder.C#L243-L285"),
    match: src("Track–cluster associations","production/macros/ElectronFinder.C#L450-L529"),
    decay: src("Select the scattered electron and decay pair",macro+"#L596-L678"),
    kine: src("Calculate W and Q²",macro+"#L683-L705"),
    tag: src("Choose the proton candidate and fill reconstructed branches",macro+"#L707-L758"),
    schema: src("Define UpsilonTree and its branches",macro+"#L267-L345"),
    fill: src("Write the event row",macro+"#L790-L800"),
    merge: src("Merge explicit per-job files","production/jobs.py#L336-L399"),
    input: src("Merged input paths and SHA-256 values","inputs/paper_inputs.csv"),
    fetch: guide("Download or verify the merged inputs","inputs/README.md"),
    normalize: src("Reweight all generated rows","analysis/scripts/validation/prepare_latest_sample.py#L60-L135"),
    replayInputs: src("Verify inputs and prepare normalized copies","paper/replay_paper.py#L31-L65"),
    selection: src("Apply the reconstructed proton requirement","analysis/workflow/proton_selection.py#L41-L63"),
    calibrated: src("Build per-state mass spectra and calibrations","analysis/scripts/validation/plot_latest_sample_performance.py#L337-L366"),
    jointFit: src("Fit the combined spectrum",plot+"reconstruction.py#L63-L125"),
    intervals: src("Profile yields and core-width intervals",plot+"mass_fit_intervals.py#L147-L217"),
    nb07: src("Notebook 07 · combined fit and efficiency maps","plotting/07_mass_and_efficiency_w_xv.ipynb"),
    maps: src("Map producer · truth, pretag and tagged masks","extensions/efficiency/diagnostics/final_maps.py#L23-L44"),
    mapRun: src("Run the map producer and set the original 1S input","paper/maps.py#L16-L46"),
    mapCheck: src("Collect both map tables and check their cells","paper/maps.py#L57-L75"),
    load: src("Read the tuples and define the response sample",optimize+"data.py#L420-L448"),
    bins: src("Choose fit bins on training events",optimize+"response.py#L422-L455"),
    slope: src("Fit the slope and run its checks",optimize+"inference.py#L571-L651"),
    spectral: src("Infer the spectrum and check coverage",optimize+"inference.py#L921-L960"),
    profile: src("Paper plot settings and display-bin rule","configs/plots_9x275.paper.json"),
    rebin: src("Combine saved bins with their joint toys",plot+"rebinning.py#L56-L98"),
    wrender: src("Render the W spectrum",plot+"w_cross_section.py#L148-L207"),
    drender: src("Render the pT² and rapidity spectra",plot+"differential.py#L236-L294"),
    nbw: src("Notebook 02 · W figure","plotting/02_w_cross_section.ipynb"),
    nbd: src("Notebook 03 · differential figure","plotting/03_differential.ipynb"),
    rendercli: src("Run the same plotting functions from the command line","plotting/render_all.py#L35-L62"),
    paperCheck: src("Compare displayed bins, endpoints and slope","paper/check_spectra.py#L30-L66"),
    paperRun: guide("Reproduce the three paper figures","paper/README.md"),
    stats: guide("Validation by observable","docs/validation.md"),
    beamGuide: guide("Run the two-energy comparison","docs/figure-explorer.md#two-energy-kinematics"),
    beamMasks: src("Build truth and reconstructed selections","analysis/scripts/reporting/render_note_kinematics.py#L45-L82"),
    beamPlot: src("Build the two-energy histogram pages","analysis/scripts/reporting/render_note_kinematics.py#L464-L518"),
    beamTask: src("Inputs and output directory of the comparison task","analysis/config/pipeline_note_latest.toml#L1556-L1601"),
    beamExport: src("Read the original figure PDFs, numerical bins and audit","tools/export_explorer_beams.py"),
    flags: guide("Selections and normalization conventions","docs/selections-and-normalization.md"),
    viewMass: page("Open the joint mass fit · notebook and note page","./figures.html#mass"),
    viewEff: page("Open the 1S maps · notebook and note page","./figures.html#efficiency"),
    viewW: page("Open the W figure · notebook and note page","./figures.html#w"),
    viewDiff: page("Open pT² and rapidity · notebook and note page","./figures.html#differential"),
    viewBeams: page("Compare the two energies · original AN PDFs","./beams.html")
  };

  const states=["1s","2s","3s"];
  const run = state => "results/generated/runs/paper-"+state+"/artifacts/";
  const perf = state => run(state)+"performance/9x275/";
  const audits = states.map(s=>run(s)+"normalization/9x275_input_audit.json");
  const reopt = "results/generated/runs/paper-reopt/artifacts/reoptimization/";
  const mapbase = "results/generated/runs/paper-1s-tagged-maps/artifacts/final_efficiency_2d/1s_newid/";
  const wPdf = "spectra/w/upsilon_es_w_cross_section.pdf";
  const diffPdf = "spectra/differential/upsilon_es_differential.pdf";
  const beamPdfs = ["@compare/projected_@state_upsilon_kinematics.pdf","@compare/projected_@state_decay_leptons.pdf",
    "@compare/projected_@state_scattered_electron_tracking.pdf","@compare/t_estimators_projective_transverse.pdf"];
  const step = (id,title,input,output,detail,links,files=[]) => ({id,title,input,output,detail,links,files});
  const reco = (title,output,detail,links) => step("reconstruction",title,"Reconstructed ePIC events tree",output,detail,[...new Set([...links,"schema","production","upstream","gap"])]);
  const tuples = (detail,files=audits) => step("tuple","Start from the merged tuples","Per-job UpsilonTree and h_EventCount",
    "One merged ROOT file per state and energy",detail,["input","fetch","normalize","merge","cards"],files);
  const nominalTuple = "The paper replay verifies the three 9×275 input files against their SHA-256 values. It writes separate normalized state files, using every generated row for the rate. The recoil-proton requirement is applied later.";
  const sharedReopt = "Normalized state tuples with generator truth, reconstructed kinematics and the recoil-proton tag";
  const pub = (title,input,output,detail,links,files) => step("publication",title,input,output,detail,links,files);

  const routes = {
    mass: {
      label:"Three-state mass fit", section:"AN §2.5 · Fig. 7 (continued), Table 8", paper:"Paper Fig. 18(a)",figure:"viewMass",preview:"plots/mass.webp",
      description:"How a reconstructed decay pair becomes a fitted yield and core width for each Υ state.",
      steps:[
        reco("Build each dielectron candidate","reco_mass, reco_topology, reco_tagged_proton_p",
          "ElectronFinder matches tracks to calorimeter clusters. The tuple macro assigns the scattered electron and an opposite-sign decay pair, then stores the pair’s mass and the proton-candidate branches.",["match","id","decay","tag"]),
        tuples(nominalTuple),
        step("analysis","Calibrate the signal shape for each state","9x275_normalized.root for 1S, 2S and 3S","mass_bins.csv + mass_fit.json for each state",
          "Apply the tagged selection and fill a mass histogram for each state. The per-state fits supply the tails and relative widths for the joint fit.",["replayInputs","selection","calibrated"],states.flatMap(s=>[perf(s)+"mass_bins.csv",perf(s)+"mass_fit.json"])),
        step("plot","Fit the three peaks in notebook 07","Three saved mass histograms and their calibrations","Combined mass panel + mass_fit_intervals.json",
          "Notebook 07 combines the calibrated templates, fits the three yields with a common mass shift and width scale, and profiles their statistical intervals. The final fit happens here.",["nb07","jointFit","intervals"],["@fitDir/mass_fit_intervals.json","@figure"]),
        pub("Locate the fit in both documents","Reconstruction figure · panel (a)","AN Fig. 7 (continued) and Table 8 → paper Fig. 18(a)",
          "The AN also shows single-state mass calibrations and broad-window counts. Table 8 is the simultaneous tagged fit. The manuscript includes upsilon_es_reconstruction.pdf under fig:upsilon_es_reconstruction.",["viewMass","paperRun"],["@figure","paper_replay.json"])
      ]},
    efficiency: {
      label:"Tagged 1S efficiency maps",section:"AN §2.7 · Fig. 13",paper:"Paper Fig. 18(b,c)",figure:"viewEff",preview:"plots/mass.webp",
      description:"Follow an eligible simulated 1S decay into the denominator, then through reconstruction and the proton tag into the numerator.",
      steps:[
        reco("Record truth and reconstructed candidates","true_pT, true_y, true_Q2, true_W; reco_mass and proton tag",
          "Generator branches locate each decay in a truth bin. The reconstructed pair and proton branches determine whether that same row survives the cuts. Unreconstructed decays remain in the tuple.",["entry","tag","fill"]),
        tuples("The map producer reads the original merged 9×275 1S tuple. Each efficiency cell counts simulated events. The mass panel drawn alongside the maps uses normalized weights from the three state inputs.",[audits[0]]),
        step("analysis","Count eligible, pretag and tagged events","Original merged 1S tuple and the final-selection definition","final_efficiency_2d.csv and selection_audit.json",
          "The efficiency extension accumulates integer counts in truth coordinates, including pT and rapidity. A cell with zero eligible entries has undefined efficiency.",["maps","mapRun","selection"],[mapbase+"tables/final_efficiency_2d.csv",mapbase+"selection_audit.json"]),
        step("plot","Build the two displayed maps","Saved pT–rapidity counts and the original 1S tuple","reconstruction/efficiency_cells.csv + panels (b,c)",
          "Notebook 07 reads the saved pT–rapidity map. It reopens the original tuple to calculate W and xᵥ jointly. The replay collects both grids and checks the numerator, pretag count and denominator in every cell.",["nb07","mapCheck"],["reconstruction/efficiency_cells.csv","@figure"]),
        pub("Locate the efficiency panels","Reconstruction figure · panels (b,c)","AN Fig. 13 → paper Fig. 18(b,c)",
          "Panel (b) uses true pT and laboratory rapidity; panel (c) uses true W and xᵥ. The AN distinguishes the phase-space average from the conditional proton-tag efficiency. Both maps share paper label fig:upsilon_es_reconstruction.",["viewEff","paperRun"],["@figure","reconstruction/efficiency_cells.csv"])
      ]},
    w: {
      label:"Cross section versus W",section:"AN §2.9 · Fig. 14; Appendix B",paper:"Paper Fig. 19",figure:"viewW",preview:"plots/w.webp",
      description:"Follow the scattered-electron measurement into reco_W, the response matrix, and the flux-normalized spectrum.",
      steps:[
        reco("Reconstruct the photon–proton energy","true_W, reco_W, true_Q2, recoil-proton tag",
          "The recovered macro uses the electron method above yₑ = 0.045 and the double-angle estimate otherwise. It keeps the generated and reconstructed W for the response calculation.",["kine","tag"]),
        tuples(nominalTuple),
        step("analysis","Infer the W spectrum",sharedReopt,"Per-state w/result.json, response matrices and validation_status.csv",
          "Training events define bins and the response. Held-out events and toy studies check inference. The visible dielectron result is converted using the branching fraction and the photon flux in each state’s W bin.",["load","bins","spectral","flags"],[...states.map(s=>reopt+s.toUpperCase()+"/w/result.json"),reopt+"validation_status.csv"]),
        step("plot","Apply the paper’s display bins","Saved W inference, covariance, toys and plot profile","upsilon_es_w_cross_section.pdf and displayed_points.csv",
          "The replay calls the same renderer as notebook 02. For 1S it keeps the outer fit bins and combines the interior bins with the joint toys; the plotting step does not refit the response.",["wrender","nbw","profile","rebin"],[wPdf,"spectra/w/upsilon_es_w_cross_section_points.csv","spectra/displayed_points.csv"]),
        pub("Check what the W points mean","State-resolved W figure and its validation record","AN Fig. 14 and Table 14 → paper Fig. 19",
          "The projected W spectra reproduce the reference points and bands. All three states failed coverage for alternative injected spectra. The supplied manuscript uses fig:upsilon_es_w.",["viewW","stats","paperCheck"],[wPdf,reopt+"validation_status.csv","spectra/comparison.json"])
      ]},
    differential: {
      label:"pT² slope and rapidity spectra",section:"AN §2.9 · Fig. 15; Appendices B–C",paper:"Paper Fig. 20(a,b)",figure:"viewDiff",preview:"plots/differential.webp",
      description:"Distinguish the slope of the dielectron pT² spectrum from a fit targeting true momentum transfer.",
      steps:[
        reco("Reconstruct the Υ four-momentum","true_pT, reco_pT, true_y, reco_y, recoil-proton tag",
          "The decay-pair four-vector gives the meson pT and rapidity. The proton candidate selects tagged events; it is not the momentum used for the dielectron pT² spectrum.",["decay","tag"]),
        tuples(nominalTuple),
        step("analysis","Fit the slope and the rapidity bins",sharedReopt,"pt2/result.json, y/result.json and validation_status.csv",
          "The paper’s slope targets true pT² using reconstructed pT². The separate abs_t_via_pt2 case targets true |t|. Rapidity uses spectrum inference and the paper divides the ep spectrum by the state’s integrated photon flux.",["load","bins","slope","spectral"],[reopt+"1S/pt2/result.json",reopt+"1S/abs_t_via_pt2/result.json",reopt+"1S/y/result.json",reopt+"validation_status.csv"]),
        step("plot","Draw the available fits","Saved pT² and rapidity cases + 1S slope profile","upsilon_es_differential.pdf, slopes and per-bin points",
          "Notebook 03 and the command-line renderer use the saved fits. The rapidity panel merges the interior 1S bins for display. Excited states without a feasible pT² fit have no slope point.",["drender","nbd","profile","rendercli"],[diffPdf,"spectra/differential/upsilon_es_slopes.csv","spectra/differential/upsilon_es_differential_points.csv"]),
        pub("Check the scope of the slope and rapidity panels","Differential figure, profile limits and statuses","AN Fig. 15 → paper Fig. 20(a,b)",
          "The 1S pT² slope passed the conditional, signal-only checks. The rapidity spectra failed coverage; the true-|t| fit failed tail-shape closure. The supplied manuscript uses fig:upsilon_es_t.",["viewDiff","stats","paperCheck"],[diffPdf,reopt+"validation_status.csv","spectra/comparison.json"])
      ]},
    beams: {
      label:"Two-energy kinematics",section:"AN §§2.5–2.6 · Figs. 10–12",paper:"AN comparison only",figure:"viewBeams",
      description:"Compare 9×130 and 9×275 at 10 fb⁻¹ per energy, with no recoil-proton requirement.",
      steps:[
        reco("Record the decay pair and scattered electron","Reconstructed masses, momenta, track counts and t estimators",
          "The recovered source shows the branch definitions and electron-identification logic. Its proton-beam energy is fixed at 275 GeV. The original 130 GeV producer revision is still needed to reproduce that upstream job.",["entry","decay","kine","schema"]),
        tuples("Use all six merged input files. Normalize each state and energy independently, then pass the low- and high-energy tuples for the same state to an inclusive workspace.",[]),
        step("analysis","Select both energies at equal luminosity","Normalized tuples for the same Υ state at both energies","kinematic_bins.csv and kinematic_selection_audit.json",
          "The note renderer applies the inclusive reconstructed selection and beam-specific W windows. Its audit records the category, luminosity and masks. In Fig. 12 the projective estimator uses its angular weight.",["beamGuide","beamMasks","beamTask"],["@compare/kinematic_selection_audit.json","@compare/kinematic_bins.csv"]),
        step("plot","Overlay the two beam energies","The two selected samples and their histogram definitions","Original per-state PDFs in artifacts/note/compare_9x130_9x275/",
          "The same renderer writes the bin table and the figure PDFs. Figures 10 and 11 use visible yields per bin. Figure 12 compares projective and dielectron-pT² yield densities for 1S.",["beamPlot","beamExport"],beamPdfs),
        pub("Find the beam comparison in the AN","Per-state two-energy plots","AN Fig. 10 (meson), Fig. 11 (leptons and tracking), Fig. 12 (1S estimators)",
          "These plots answer how the two beams compare at the same luminosity. They are separate from the tagged 9×275 paper figures.",["viewBeams","beamGuide"],beamPdfs)
      ]}
  };
  const boundary = "The recovered producer fixes the proton beam at 275 GeV. The exact per-energy CERN producer, job environment and merge manifest are still missing. Reproduction has been verified starting from the checksum-listed merged tuples.";
  window.UPSiLON_TRACE = {referenceRevision, boundary, sources, routes};
})();
