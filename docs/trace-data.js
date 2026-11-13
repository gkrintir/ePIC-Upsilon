/* Stages are listed in production order. Source links identify the checked
   implementation. Local paper file templates retain the historical replay layout. */
(() => {
  const referenceRevision = "8e5d380f3d12351e9061bc41f660528c6ae163f4";
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
    gap: guide("Historical stored-matching producer comparisons", "docs/missing-production-sources.md"),
    eventCheck: guide("Compare two original inputs with the merged events", "docs/missing-production-sources.md#two-reconstructed-inputs"),
    cernCampaign: guide("Trace the 2S 9×130 Condor run into the merged tuple", "production/cern-2026-10/README.md"),
    producer130: src("9×130 working copy · beam setting", "production/macros/9x130/DVMP_Upsilon_Analysis.C#L35-L42"),
    entry: src("Tuple macro · entry point and ID mode",macro+"#L166-L188"),
    id: src("ElectronFinder · electron identification","production/macros/ElectronFinder.C#L243-L285"),
    match: src("Corrected downstream track–cluster matching","production/rematch3011/ElectronFinder.patch"),
    rematch: src("Enable the rematched reader", "production/rematch3011/RunRematch3011.C"),
    rematchGuide: guide("Reprocess the six samples with corrected matching", "docs/rematch-3011.md"),
    decay: src("Select the scattered electron and decay pair",macro+"#L593-L678"),
    kine: src("Calculate W and Q²",macro+"#L683-L705"),
    tag: src("Choose the proton candidate and fill reconstructed branches",macro+"#L707-L758"),
    schema: src("Define UpsilonTree and its branches",macro+"#L267-L345"),
    fill: src("Write the event row",macro+"#L793-L800"),
    merge: src("Merge validated per-job files","production/jobs.py#L691-L802"),
    input: page("Corrected input paths, SHA-256 values and generated counts","./data/paper-provenance.json"),
    fetch: guide("Prepare and audit rematched merged inputs","docs/rematch-3011.md"),
    normalize: src("Reweight all generated rows","analysis/scripts/validation/prepare_latest_sample.py#L60-L135"),
    replayInputs: src("Historical replay · verify stored-matching inputs","paper/replay_paper.py#L31-L65"),
    selection: src("Apply the reconstructed proton requirement","analysis/workflow/proton_selection.py#L41-L63"),
    calibrated: src("Build per-state mass spectra and calibrations","analysis/scripts/validation/plot_latest_sample_performance.py#L337-L366"),
    jointFit: src("Fit the combined spectrum",plot+"reconstruction.py#L63-L125"),
    combineBins: src("Read the three spectra and combine adjacent mass bins",plot+"reconstruction.py#L148-L213"),
    intervals: src("Profile yields and core-width intervals",plot+"mass_fit_intervals.py#L147-L217"),
    nb07: src("Notebook 07 · combined fit and efficiency maps","plotting/07_mass_and_efficiency_w_xv.ipynb"),
    maps: src("Map producer · truth, pretag and tagged masks","extensions/efficiency/diagnostics/final_maps.py#L23-L44"),
    mapCounts: src("Count the same events before and after selection in each truth bin","extensions/efficiency/diagnostics/final_maps.py#L106-L154"),
    mapRun: src("Historical replay · map input configuration","paper/maps.py#L16-L46"),
    mapCheck: page("Corrected map totals, grid counts and source hashes","./data/paper-provenance.json"),
    load: src("Read the tuples and define the response sample",optimize+"data.py#L408-L448"),
    bins: src("Choose fit bins on training events",optimize+"response.py#L422-L455"),
    response: src("Count selected events by generated and reconstructed bin",optimize+"response.py#L239-L300"),
    inferenceChoice: src("Pass the response and held-out sample to the fit","analysis/workflow/run_reoptimization.py#L292-L310"),
    foldedShape: src("Fold a trial exponential through the detector response",optimize+"inference.py#L119-L143"),
    fluxScaling: src("Convert fitted yields to cross sections",optimize+"data.py#L174-L224"),
    slope: src("Fit the slope and run its checks",optimize+"inference.py#L571-L651"),
    spectral: src("Infer the spectrum and check coverage",optimize+"inference.py#L921-L960"),
    profile: src("Paper plot settings and display-bin rule","configs/plots_9x275.paper.json"),
    rebin: src("Combine saved bins with their joint toys",plot+"rebinning.py#L56-L98"),
    wrender: src("Render the W spectrum",plot+"w_cross_section.py#L148-L207"),
    drender: src("Render the pT² and rapidity spectra",plot+"differential.py#L236-L294"),
    nbw: src("Notebook 02 · W figure","plotting/02_w_cross_section.ipynb"),
    nbd: src("Notebook 03 · differential figure","plotting/03_differential.ipynb"),
    rendercli: src("Run the same plotting functions from the command line","plotting/render_all.py#L35-L62"),
    paperCheck: src("Historical replay · compare stored-matching spectral values","paper/check_spectra.py#L30-L66"),
    paperRun: guide("Replay the historical stored-matching paper figures","paper/README.md"),
    stats: guide("Validation by observable","docs/validation.md"),
    beamGuide: guide("Run the two-energy comparison","docs/figure-explorer.md#two-energy-kinematics"),
    beamMasks: src("Build truth and reconstructed selections","analysis/scripts/reporting/render_note_kinematics.py#L51-L82"),
    beamBranches: src("Store the reconstructed meson and lepton kinematics",macro+"#L727-L741"),
    beamYields: src("Sum selected event weights in each bin","analysis/scripts/reporting/render_note_kinematics.py#L234-L260"),
    beamWindows: src("Reconstructed W range for each beam","analysis/config/pipeline_note_latest.toml#L40-L57"),
    beamPlot: src("Build the two-energy histogram pages","analysis/scripts/reporting/render_note_kinematics.py#L504-L553"),
    beamTask: src("Inputs and output directory of the comparison task","analysis/config/pipeline_note_latest.toml#L1556-L1601"),
    beamExport: src("Read the original figure PDFs, numerical bins and audit","tools/export_explorer_beams.py"),
    results: page("Corrected fit values, map cells and spectral points", "./data/paper.json"),
    provenance: page("Full-sample rematch provenance and validation", "./data/paper-provenance.json"),
    beamProvenance: page("Corrected two-energy figure sources", "./data/beams-provenance.json"),
    counters: page("AN Fig. 6 · corrected producer counters", "./plots/producer-counters.pdf"),
    counterTable: page("AN Fig. 6 · raw counts and defined adjacent losses", "./data/producer-counters.csv"),
    counterCode: src("Framework renderer · the same producer counter definitions", "analysis/scripts/reporting/render_note_kinematics.py#L395-L461"),
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
  const reco = (title,output,detail,links) => step("reconstruction",title,"Reconstructed ePIC events tree",output,detail,[...new Set([...links,"schema","production","rematchGuide","upstream","gap"])]);
  const tuples = (title,output,detail,files=audits) => step("tuple",title,"Per-job UpsilonTree and h_EventCount",
    output,detail,["input","fetch","normalize","merge","cards"],files);
  const mergedTuples = (branches, detail, links) => reco("Build and merge the event tuples",
    "Three corrected merged ROOT files with " + branches, detail + " The [rematching wrapper][rematch] enables the [corrected track–cluster assignment][match] before producing these tuples. The [production merge][merge] combines the per-job tuples separately for 1S, 2S and 3S.", [...links, "match", "rematch", "merge"]);
  const nominalTuples = next => step("tuple","Normalize the merged state samples",
    "Three merged 9×275 ROOT files, one per Υ state", "One 9x275_normalized.root per state, weighted at 10 fb⁻¹",
    "The current projection starts with the full rematched samples in the [verified input manifest][input]. The [normalization step][normalize] weights each generated event by the visible dielectron cross section times 10 fb⁻¹, divided by the actual total generated for that state. Each state then has the expected rate for that luminosity, regardless of how many jobs produced its tuples.\n\n" + next,
    ["input","fetch","normalize","provenance"], audits);
  const sharedReopt = "Normalized state tuples with generator truth, reconstructed kinematics and the recoil-proton tag";
  const pub = (title,input,output,detail,links,files) => step("publication",title,input,output,detail,links,files);

  const routes = {
    mass: {
      label:"Three-state mass fit", section:"AN §2.5 · Fig. 7 (continued), Table 8", paper:"Paper Fig. 18(a)",figure:"viewMass",preview:"plots/mass.webp?v=7a2b6cc265ba",
      question:"How many events can we assign to each Υ peak in the same dielectron spectrum?",
      description:"Simulating each state separately lets us calibrate its reconstructed line shape. Fitting those shapes to the summed mass spectrum gives a projected yield for each state in the proton-tagged 9×275 GeV sample at 2.5 fb⁻¹.",
      steps:[
        mergedTuples("reco_mass, reco_topology and reco_tagged_proton_p",
          "The final state contains both a scattered electron and an electron from the Υ decay. After [electron identification][id], the tuple macro [assigns the scattered electron and the decay pair][decay]. It saves the pair’s invariant mass and [records whether a recoil proton was reconstructed][tag], so the analysis can build the mass spectrum for the tagged sample.",["match","id","decay","tag"]),
        nominalTuples("The [mass analysis][mass/analysis] uses these weights to project each state’s selected spectrum to 2.5 fb⁻¹ for the combined fit."),
        step("analysis","Calibrate each Υ peak","9x275_normalized.root for 1S, 2S and 3S","mass_bins.csv + mass_fit.json for each state",
          "In a sample containing only one generated Υ state, the spread of reconstructed masses measures that state’s line shape. After the [proton tag][selection], the [per-state analysis][calibrated] fits the unweighted selected masses and separately fills a histogram with the expected yield at 2.5 fb⁻¹. For the [combined fit][mass/plot], mass_fit.json supplies the peak widths and tails, while mass_bins.csv supplies the expected events in each bin.",["provenance","selection","calibrated"],states.flatMap(s=>[perf(s)+"mass_bins.csv",perf(s)+"mass_fit.json"])),
        step("plot","Combined fit and figure","Three saved mass histograms and their calibrations","Mass panel, mass_fit_table.csv and mass_fit_intervals.json",
          "An event in the tail of one peak can fall under its neighbour. Notebook 07 [adds the three expected spectra in 50 MeV bins][combineBins] and uses the [single-state calibrations][mass/analysis] to describe the overlap. The [joint fit][jointFit] fixes the mass separations, tail shapes and relative widths. It lets the three yields vary independently, with one common mass shift and one common width scale.\n\nTo get the uncertainty on the 1S yield, the [interval calculation][intervals] holds that yield at a trial value and refits the other four parameters. The points where the Poisson deviance rises by one give the lower and upper errors beside 91.8 events in [Table 8][mass/publication].",["nb07","combineBins","jointFit","intervals"],["@fitDir/mass_fit_table.csv","@fitDir/mass_fit_intervals.json","@figure"]),
        pub("Analysis note and paper","Combined mass panel and fit table","AN Fig. 7 (continued), Table 8 and paper Fig. 18(a)",
          "The final page of AN Fig. 7 and the [Fig. 18(a) layout][viewMass] show the corrected fitted spectrum. From that same fit, Table 8 gives a Υ(1S) yield of 91.8 (+10.2/−9.5) expected events and a Gaussian core width of 54.4 (+6.9/−6.1) MeV. These conditional statistical intervals assume signal only and fixed MC tail shapes. All three widths depend on the same fitted scale, so their fractional uncertainties are fully correlated. The [validation record][provenance] retains the calibration-shape mismatch flags and the untested interval coverage.",["viewMass","results","provenance"],["@figure","paper_replay.json"])
      ]},
    efficiency: {
      label:"Tagged 1S efficiency maps",section:"AN §2.7 · Fig. 13",paper:"Paper Fig. 18(b,c)",figure:"viewEff",preview:"plots/mass.webp?v=7a2b6cc265ba",
      question:"How does the probability of reconstructing a proton-tagged Υ(1S) vary across phase space?",
      description:"In each bin of generated kinematics, count the eligible decays and how many survive each selection. The fraction that reaches the tagged sample gives the efficiency for that bin.",
      steps:[
        reco("Record truth and reconstructed candidates","true_pT, true_y, true_Q2, true_W; reco_mass and proton tag",
          "The denominator must include eligible decays even when reconstruction loses an electron. The tuple macro [writes an event row][fill] with the generated kinematics in either case. When it finds a candidate, it also [records the reconstructed mass and proton information][tag]. The [full corrected 1S input][input] contains 1,137,600 generated rows, of which 1,126,936 pass the truth eligibility cuts. Each selected decay contributes to the same rapidity–pT bin as its generated event. The [rematching wrapper][rematch] applies the corrected ECAL associations before the producer fills its reconstructed branches.",["entry","tag","fill","match","rematch"]),
        tuples("Merge the 1S tuples","Full rematched 9×275 1S ROOT file","The [merge][merge] collects the simulated 1S events across jobs. The [efficiency calculation][efficiency/analysis] needs integer counts in that sample, so it reads the rematched merged tuple. The [mass fit][mass/analysis] uses luminosity weights because its histogram predicts how many events the experiment would see.",[audits[0]]),
        step("analysis","Count eligible, pretag and tagged events","Full rematched 1S tuple and the final-selection definition","final_efficiency_2d.csv and selection_audit.json",
          "The [selection code][maps] first identifies eligible decays using generated mass, Q² and W. It then selects reconstructed candidates using the dielectron mass, W and the event’s E − pz balance, and finally requires a reconstructed recoil proton. For each truth bin the [accumulator][mapCounts] saves all three counts, so one can distinguish losses in reconstructing and selecting the leptons from losses at the proton tag.",["maps","mapCounts","provenance","selection"],[mapbase+"tables/final_efficiency_2d.csv",mapbase+"selection_audit.json"]),
        step("plot","Build the two displayed maps","Saved pT–rapidity counts and the rematched 1S tuple","efficiency_map_cells.csv + panels (b,c)",
          "[Notebook 07][nb07] can draw the pT–rapidity map directly from the saved bin counts. To find the efficiency at a particular pair of W and xᵥ values, it needs both quantities for the same decay. It therefore reopens the rematched 1S tuple and counts events on that grid with the same selections. The [current provenance record][mapCheck] records all 240 pT–rapidity and 288 xᵥ–W cells, including the three selected decays outside the displayed pT–rapidity range.",["nb07","mapCheck"],["reconstruction/efficiency_cells.csv","@figure"]),
        pub("Analysis note and paper","Reconstruction figure · panels (b,c)","AN Fig. 13 and paper Fig. 18(b,c)",
          "Across the eligible sample, 301,323 of 1,126,936 decays pass the lepton selection and 198,265 also have a recoil-proton tag. Thus 17.59% of eligible decays pass the full selection, and 65.80% of candidates already passing the lepton selection retain a tag. AN Fig. 13 and [paper Fig. 18(b,c)][viewEff] show the first fraction bin by bin, as a function of pT and rapidity or W and xᵥ.",["viewEff","results","provenance"],["@figure","reconstruction/efficiency_cells.csv"])
      ]},
    w: {
      label:"Cross section versus W",section:"AN §2.9 · Fig. 14; Appendix B",paper:"Paper Fig. 19",figure:"viewW",preview:"plots/w.webp?v=105f0a3231a1",
      question:"How does the Υ production cross section vary with photon–proton energy?",
      description:"Figure 19 gives a cross section averaged over the photon flux in each W interval. The [mass fit][mass/plot] gives one yield per state integrated over W. For the spectrum, the framework returns to the same simulated events and [infers a separate yield in each W bin][w/analysis].",
      steps:[
        mergedTuples("true_W, reco_W, true_Q2 and the recoil-proton tag",
          "An event generated in one W interval can be reconstructed in another. The tuple stores both values for the same event. To [estimate reconstructed W][kine], the macro uses the scattered electron when the electron-method inelasticity yₑ exceeds 0.045 and the double-angle estimate otherwise.",["kine","tag"]),
        nominalTuples("The [W analysis][w/response] scales expected counts to the paper exposure of 2.5 fb⁻¹. It retains generated and reconstructed W on the same event to estimate migration between bins."),
        step("response","Select tagged events and build the W response",sharedReopt,"W bin edges, detector response and held-out reconstructed counts for each state",
          "The [selection][load] identifies eligible generated decays and the subset with a reconstructed dielectron candidate and recoil-proton tag. The [bin optimization][bins] uses training events only. Within each generated W bin, the [response matrix][response] records the fraction reconstructed in each measured W bin. Its column sum gives the selection efficiency, and the off-diagonal entries describe migration. A separate sample supplies the expected reconstructed counts for the [spectrum fit][w/analysis].",["load","selection","bins","response","inferenceChoice"],states.map(s=>reopt+s.toUpperCase()+"/w/response_matrices.npz")),
        step("analysis","Infer the W spectrum","W response and held-out reconstructed counts for each state","Per-state w/result.json and validation_status.csv",
          "The [spectrum fit][spectral] varies the generated yields until applying the [W response][w/response] predicts the reconstructed counts in the held-out sample. To [obtain a γp cross section][fluxScaling], divide each fitted yield by the luminosity, the dielectron branching fraction and the photon flux integrated over that W bin.",["spectral","fluxScaling","flags"],[...states.map(s=>reopt+s.toUpperCase()+"/w/result.json"),reopt+"validation_status.csv"]),
        step("plot","Apply the paper’s display bins","Saved W inference, covariance, toys and plot profile","upsilon_es_w_cross_section.pdf and displayed_points.csv",
          "For 1S, the [display rule][profile] applied to the rematched fit gives 20–50, 50–65 and 65–100 GeV. To obtain a merged-bin cross section, the [plotting code][rebin] adds the fitted yields and divides by the combined flux with the same luminosity and branching fraction. It sums the bins in each simulated repeat of the fit before taking quantiles, preserving correlations in the displayed band.",["wrender","nbw","profile","rebin"],[wPdf,"spectra/w/upsilon_es_w_cross_section_points.csv","spectra/displayed_points.csv"]),
        pub("Analysis note and paper","State-resolved W figure and its validation record","AN Fig. 14, Table 14 and paper Fig. 19",
          "In the [corrected W spectrum][viewW], the middle 1S point is 117.2 pb for 50–65 GeV. Its 90.9–142.8 pb band shows the spread of the fitted estimates in repeated signal-only samples. The [numerical export][results] includes all nine displayed W points. The [validation record][provenance] reports failed injected-spectrum coverage and shape/combined closure for all three states. The adaptive-binning spectra are available as diagnostics. The separate fixed-bin 2S and 3S unfolding is withheld because a truth-populated bin has zero reconstruction efficiency.",["viewW","stats","results","provenance"],[wPdf,reopt+"validation_status.csv","spectra/comparison.json"])
      ]},
    differential: {
      label:"pT² slope and rapidity spectra",section:"AN §2.9 · Fig. 15; Appendices B–C",paper:"Paper Fig. 20(a,b)",figure:"viewDiff",preview:"plots/differential.webp?v=f3f1c863f1ba",
      question:"How steeply does the Υ spectrum fall with pT², and how does production vary with rapidity?",
      description:"Both panels use the reconstructed decay pair. Panel (a) extracts a slope from the pT² distribution. Panel (b) estimates the yield separately in each rapidity bin.",
      steps:[
        mergedTuples("true_pT, reco_pT, true_y, reco_y and the recoil-proton tag",
          "The macro [assigns the two decay tracks][decay] and sums their four-momenta. That sum gives the reconstructed Υ pT and laboratory rapidity. The macro also [records the recoil-proton candidate][tag] for the later tagged selection. Keeping the generated meson kinematics on the same row lets the analysis measure how often the reconstructed pair lands in a different bin.",["decay","tag"]),
        nominalTuples("The [differential analysis][differential/response] scales expected counts to 2.5 fb⁻¹. Each input pairs the generated Υ pT and rapidity with the reconstructed values, so the two detector responses can be built from the same events."),
        step("response","Select tagged events and build the two responses",sharedReopt,"pT² and rapidity bin edges, detector responses and held-out reconstructed counts",
          "The [selection][load] identifies eligible generated decays and the subset that passes the dielectron and recoil-proton requirements. The analysis uses the squared generated and reconstructed pT values to build the pT² response. On training events the [optimizer][bins] chooses bins separately for pT² and rapidity. For each observable, the [response matrix][response] records the fraction of decays in each generated bin that pass selection and where they appear in reconstructed bins. Separate events supply the expected reconstructed distributions for the [two fits][differential/analysis].",["load","selection","bins","response","inferenceChoice"],[reopt+"1S/pt2/response_matrices.npz",...states.map(s=>reopt+s.toUpperCase()+"/y/response_matrices.npz")]),
        step("analysis","Fit the pT² slope and the rapidity yields","The two detector responses and the held-out pT² and rapidity distributions","pt2/result.json, y/result.json and validation_status.csv",
          "The [analysis driver][inferenceChoice] fits the observables independently. For each trial slope b, it [folds an exponential in true pT² through the pT² response][foldedShape] to predict the reconstructed shape. The [slope fit][slope] allows the signal normalization to vary and includes a fixed MC estimate of signal migrating in from outside the fitted true-pT² range. The 1S and 2S cases have feasible binning. The 2S interval fails injected-slope coverage, and the 3S case has no feasible binning.\n\nFor rapidity, the [spectrum fit][spectral] varies the yield in each generated bin. Dividing by luminosity, the dielectron branching fraction and the bin width gives dσ(ep)/dy. The [paper’s rapidity distribution][fluxScaling] also divides by the photon flux integrated over that state’s W and Q² window.",["inferenceChoice","foldedShape","slope","spectral","fluxScaling"],[reopt+"1S/pt2/result.json",reopt+"1S/abs_t_via_pt2/result.json",...states.map(s=>reopt+s.toUpperCase()+"/y/result.json"),reopt+"validation_status.csv"]),
        step("plot","Draw the pT² and rapidity panels","Saved pT² and rapidity cases + 1S slope profile","upsilon_es_differential.pdf, slopes and per-bin points",
          "In panel (a), the [renderer][drender] reads the saved binwise pT² estimates for the spectrum and the separate exponential fit for the slope inset. The 1S slope passes conditional validation. The 2S slope is shown as a diagnostic with failed coverage, and the 3S slope is unavailable. In rapidity, the [display merge][rebin] sums the interior 1S yields and divides by the merged bin width. The 2S and 3S points keep their fitted bins.",["drender","nbd","profile","rebin","rendercli"],[diffPdf,"spectra/differential/upsilon_es_slopes.csv","spectra/differential/upsilon_es_differential_points.csv"]),
        pub("Analysis note and paper","Differential figure, profile limits and statuses","AN Fig. 15 and paper Fig. 20(a,b)",
          "The [corrected 1S slope][viewDiff] is 2.46 (+0.55/−0.53) GeV⁻² over 0–1 GeV² in pT². The slope fit profiles the normalization with a fixed detector response. Its interval passed the coverage check in signal-only pseudoexperiments at the tested slopes. The 2S slope is a diagnostic with failed injected-slope coverage. A separate 1S fit uses reconstructed pT² to infer true |t| and fails coverage, tail-shape and combined-closure checks. The [rapidity checks][stats] found inadequate interval coverage for alternative generated spectra in all three states.",["viewDiff","stats","results","provenance"],[diffPdf,reopt+"validation_status.csv","spectra/comparison.json"])
      ]},
    beams: {
      label:"Two-energy kinematics",section:"AN §§2.5–2.6 · Figs. 10–12",paper:"AN comparison only",figure:"viewBeams",
      question:"At the same luminosity, how do the reconstructed Υ and lepton spectra change with proton beam energy?",
      description:"Figures 10–12 compare the same Υ state at 9×130 and 9×275 GeV. Both samples represent 10 fb⁻¹, so the bin heights retain the different production rates and the fraction of decays selected at each energy.",
      steps:[
        reco("Record the decay pair and scattered electron","Meson rapidity, decay-lepton momenta and scattered-electron direction in each tuple",
          "A change in the meson boost can move its decay leptons into a different part of the detector. The recovered tuple macro [assigns the decay pair and scattered electron][decay], then [saves the pair rapidity and individual lepton directions][beamBranches]. Those branches supply the meson and lepton panels in Figs. 10 and 11. The [9×130 source][producer130] uses the lower proton energy for the event kinematics. The current six-sample [input manifest][input] records each corrected file’s checksum and generated-row count. The [rematching wrapper][rematch] enables the corrected track–cluster assignment at both energies.",["decay","beamBranches","match","rematch","input"]),
        tuples("Merge and normalize the tuples","One normalized ROOT file for each state and energy","The samples contain different numbers of simulated events. After the [per-job files are merged][merge], the [normalization step][normalize] assigns each generated row its state's visible cross section at that beam energy times 10 fb⁻¹, divided by the number generated. Summing selected rows then gives the expected yield for the same exposure at both energies.",[]),
        step("analysis","Select both energies at equal luminosity","Normalized tuples for the same Υ state at both energies","kinematic_bins.csv and kinematic_selection_audit.json",
          "The [comparison task][beamTask] pairs the two energies for one Υ state. After generated-mass and Q² cuts, [load_sample][beamMasks] keeps dielectron masses of 8–12 GeV and requires the three leptons’ Σ(E−p_z) to lie within 1.8 GeV of 18 GeV. The [reconstructed W ranges][beamWindows] are 20–70 GeV at 9×130 and 20–100 GeV at 9×275. The [inclusive selection][selection] accepts events with or without a reconstructed recoil proton. This comparison applies the reconstructed W windows and leaves generated W unrestricted.",["beamMasks","beamWindows","selection","beamTask"],["@compare/kinematic_selection_audit.json","@compare/kinematic_bins.csv"]),
        step("plot","Overlay the two beam energies","The two selected samples and their histogram definitions","Original per-state PDFs in artifacts/note/compare_9x130_9x275/",
          "For Figs. 10 and 11, [draw_pages][beamYields] sums the selected event weights in the same bins at both energies. For Fig. 12 it divides by the bin width and uses the angular weight only for the projective |t| estimate. The [figure exporter][beamExport] reads the plotted values from the renderer’s bin table and copies its PDFs. The [beam provenance record][beamProvenance] records the current source tables and PDF hashes.",["beamYields","beamPlot","beamExport"],beamPdfs),
        pub("Analysis note","Per-state two-energy plots","AN Fig. 10 (meson), Fig. 11 (leptons and tracking), Fig. 12 (1S estimators)",
          "Use [the beam comparison][viewBeams] to read the Υ rapidity distribution in Fig. 10 alongside the lepton directions in Fig. 11. A difference between the curves can reflect both the available phase space and which decays pass selection. The [joint mass fit][mass/plot] gives the expected yields and widths of the three peaks with a reconstructed proton at 9×275 GeV and 2.5 fb⁻¹.\n\n[Figure 6][counters] and its [counter table][counterTable] show the raw 1S producer counts at both energies. Truth stages 1–3 and reconstructed stages 4–10 form separate sequences. The [renderer][counterCode] breaks the curve between stages 3 and 4 and leaves that adjacent-loss entry undefined. The 54.3% loss from stage 5 to 6 at 9×275 measures the additional scattered-electron candidate requirement, which can be satisfied by a central candidate or a backward tagger candidate.",["viewBeams","beamGuide","beamProvenance","counters","counterTable","counterCode"],beamPdfs)
      ]}
  };
  const boundary = "The current projection uses all six full-sample rematched ROOT files with verified SHA-256 values and normalization from their actual generated-row counts. The 9×275 1S file has the same 1,137,600 truth records as its earlier stored-matching version. The linked historical producer comparisons document the original tuple campaign.";
  // The text of a reference reads as part of the explanation. Its key identifies
  // either a pinned source or another stage in this guide.
  function referenceParts(text) {
    const parts = [];
    let start = 0;
    for (const match of text.matchAll(/\[([^\]\n]+)\]\[([\w/-]+)\]/g)) {
      parts.push(text.slice(start, match.index), {label: match[1], target: match[2]});
      start = match.index + match[0].length;
    }
    parts.push(text.slice(start));
    return parts;
  }
  window.UPSiLON_TRACE = {referenceRevision, boundary, sources, routes, referenceParts};
})();
