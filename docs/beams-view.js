(() => {
  "use strict";
  const catalog = window.BEAM_FIGURES;
  let data = catalog;
  const byId = id => document.getElementById(id);
  const beams = ["130", "275"];
  const families = {
    quarkonium: ["Fig. 10 · Quarkonium", 10],
    decay: ["Fig. 11 · Decay leptons", 11],
    scattered: ["Fig. 11 · Scattered electron & tracks", 11],
    estimators: ["Fig. 12 · Momentum-transfer estimators", 12]
  };
  const meta = {
    reco_mass: ["Dielectron invariant mass", "m(e⁺e⁻)", "GeV"],
    reco_pT: ["Υ transverse momentum", "Υ transverse momentum", "GeV"],
    reco_y: ["Υ laboratory rapidity", "Υ laboratory rapidity", ""],
    reco_W: ["Photon–proton energy", "Wγp", "GeV"],
    reco_ePlus_pT: ["Decay positron transverse momentum", "Decay e⁺ transverse momentum", "GeV"],
    reco_ePlus_eta: ["Decay positron pseudorapidity", "Decay e⁺ pseudorapidity η", ""],
    reco_decay_eMinus_pT: ["Decay electron transverse momentum", "Decay e⁻ transverse momentum", "GeV"],
    reco_decay_eMinus_eta: ["Decay electron pseudorapidity", "Decay e⁻ pseudorapidity η", ""],
    reco_scat_eMinus_E: ["Scattered electron energy", "Scattered e⁻ energy", "GeV"],
    reco_scat_eMinus_eta: ["Scattered electron pseudorapidity", "Scattered e⁻ pseudorapidity η", ""],
    reco_nTracks_Central: ["Central track multiplicity", "Central reconstructed tracks", "tracks"],
    reco_nTracks_Tagger: ["Tagger track multiplicity", "Tagger reconstructed tracks", "tracks"],
    reco_t_proj: ["Projective momentum transfer", "Projective |t|", "GeV²"],
    reco_pT2: ["Squared Υ transverse momentum", "Υ transverse momentum squared", "GeV²"]
  };
  const state = byId("state"), group = byId("group"), observable = byId("observable");
  const scale = byId("scale"), slider = byId("bin"), svg = byId("chart");
  let series, index = 0, geometry, urlTimer, pendingBin;
  const currentFigure = () => data[state.value][group.value];
  const hasBins = () => Boolean(series?.["130"]?.length);
  const isDensity = () => group.value === "estimators";
  const isTracks = () => meta[observable.value][2] === "tracks";
  const yUnit = () => isDensity() ? "events/GeV²" : "events/bin";
  const figureNumber = () => families[group.value][1];

  function notePageNumber() {
    const n = Number(state.value[0]) - 1;
    return {quarkonium: 34 + n, decay: 37 + 2 * n, scattered: 38 + 2 * n, estimators: 44}[group.value];
  }

  function fillSelect(select, keys, labels, selected) {
    select.replaceChildren(...keys.map(key => new Option(labels[key][0], key)));
    select.value = keys.includes(selected) ? selected : keys[0];
  }

  function number(value, precision = 5) {
    if (value === null) return "Unavailable";
    if (value === 0) return "0";
    if (Math.abs(value) < 1e-4) return value.toExponential(2);
    return Number(value.toPrecision(precision)).toLocaleString("en-US", {maximumSignificantDigits: precision});
  }

  function updateHash() {
    clearTimeout(urlTimer);
    const hash = new URLSearchParams({state: state.value, group: group.value, observable: observable.value, scale: scale.value, ...(hasBins() ? {bin: index + 1} : pendingBin === undefined ? {} : {bin: pendingBin + 1})}).toString();
    for(const link of document.querySelectorAll('a[href*="#beams/plot/review"]')) link.href="./index.html#beams/plot/review?"+hash;
    if (location.hash.slice(1) !== hash) {
      try { history.replaceState(null, "", "#" + hash); }
      catch { if (location.hash.slice(1) !== hash) location.replace("#" + hash); }
    }
  }

  function restore() {
    const args = new URLSearchParams(location.hash.slice(1));
    state.value = Object.hasOwn(data, args.get("state")) ? args.get("state") : "1S";
    fillSelect(group, Object.keys(data[state.value]), families, args.get("group"));
    fillSelect(observable, Object.keys(currentFigure().observables), meta, args.get("observable"));
    scale.value = args.get("scale") === "log" ? "log" : "linear";
    const binNumber = Number(args.get("bin"));
    showObservable(Number.isInteger(binNumber) && binNumber > 0 ? binNumber - 1 : undefined);
  }

  function showObservable(binIndex) {
    const key = observable.value, [label, axis, unit] = meta[key];
    series = currentFigure().observables[key];
    const available = hasBins();
    byId("workspace").hidden = !available;
    byId("source-view").hidden = available;
    byId("import-box").hidden = available && data === catalog;
    if (available) {
      pendingBin = undefined;
      slider.max = series["130"].length;
      const rows = series["275"];
      const peak = rows.reduce((best, row, i) => row[2] > rows[best][2] ? i : best, 0);
      index = Math.min(binIndex ?? peak, rows.length - 1);
    } else {
      pendingBin = binIndex;
      index = 0;
      const preview = byId("pdf-preview");
      if (preview.getAttribute("src") !== currentFigure().pdf) preview.src = currentFigure().pdf;
      const figLabel = 'AN Fig. ' + figureNumber() + '  /  Υ(' + state.value + ')';
      byId("pdf-id").textContent = figLabel;
      byId("pdf-title").textContent = label;
      const panelIndex = Object.keys(currentFigure().observables).indexOf(key);
      const panelNames = isDensity() ? ["left", "right"] : ["upper left", "upper right", "lower left", "lower right"];
      byId("pdf-caption").textContent = 'See the ' + panelNames[panelIndex] + ' panel of the original figure. ' + (isDensity() ? "Fig. 12 displays reconstructed yield densities in events/GeV²." : "The vertical axis gives expected visible events per bin.");
      byId("pdf-open").href = currentFigure().pdf;
      preview.title = figLabel + ' · ' + families[group.value][0];
    }
    byId("chart-id").textContent = 'AN Fig. ' + figureNumber() + '  /  Υ(' + state.value + ')';
    byId("chart-title").textContent = label;
    byId("x-title").textContent = axis + (unit && !isTracks() ? ' [' + unit + ']' : '');
    const yieldLabel = !isDensity() ? "Expected visible events per bin" :
      key === "reco_t_proj" ? "Angular-weighted yield density [events/GeV²]" : "Visible yield density [events/GeV²]";
    byId("y-title").textContent = yieldLabel;
    byId("unit-caption").textContent = yieldLabel;
    byId("source").href = currentFigure().pdf;
    byId("source").textContent = 'Original Fig. ' + figureNumber() + ' PDF ↗';
    byId("note-jump").textContent = 'Open analysis note · PDF page ' + notePageNumber() + ' ↗';
    byId("selection").textContent = 'Each beam sample uses 8–12 GeV reconstructed dielectron mass and Σ(E − pz) = 18 ± 1.8 GeV. The reconstructed Wγp window is 20–70 GeV at 9 × 130 GeV and 20–100 GeV at 9 × 275 GeV. The generated Q² limit for Υ(' + state.value + ') is ' + (state.value === "1S" ? "1" : "0.01") + ' GeV².';
    let context = "Each bin sums luminosity-scaled event weights for the selected decay pair and reconstructed scattered electron. Both beam samples include events with or without a recoil-proton tag.";
    if (key === "reco_t_proj") context = "The projective estimator retains events with ω < π/12 and multiplies their weights by six. The source flags unavailable count bounds in some bins. The density uses the reconstructed estimator before unfolding.";
    if (key === "reco_pT2") context = "The squared transverse momentum comes from the reconstructed decay pair. This density uses ordinary visible dielectron weights, before unfolding or a branching-fraction correction.";
    byId("context").textContent = context;
    window.notePage.refresh();
    if (available) { drawChart(); selectBin(index); }
    updateHash();
  }

  function svgEl(tag, attributes, content) {
    const node = document.createElementNS("http://www.w3.org/2000/svg", tag);
    for (const [key, value] of Object.entries(attributes)) node.setAttribute(key, value);
    if (content !== undefined) node.textContent = content;
    return node;
  }

  function stepSize(span, n) {
    const size = span / n, factor = 10 ** Math.floor(Math.log10(size));
    const ratio = size / factor;
    return (ratio <= 1 ? 1 : ratio <= 2 ? 2 : ratio <= 2.5 ? 2.5 : ratio <= 5 ? 5 : 10) * factor;
  }

  function drawChart() {
    if (!hasBins() || byId("workspace").hidden) return;
    const width = Math.max(svg.clientWidth, 200), height = svg.clientHeight || 300;
    const left = 43, right = width - 12, top = 13, bottom = height - 26;
    const rows = series["130"], minX = rows[0][0], maxX = rows[rows.length - 1][1];
    const all = [...rows, ...series["275"]].flatMap(r => r.slice(2)).filter(n => n > 0);
    const log = scale.value === "log";
    const logFloor = Math.floor(Math.log10(Math.max(1e-300, all.length ? Math.min(...all) : 1)));
    const minY = log ? logFloor : 0;
    const maxData = Math.max(...all, 1e-12);
    const linearStep = stepSize(maxData, 5);
    const maxY = log ? Math.max(logFloor + 1, Math.ceil(Math.log10(maxData))) : Math.ceil(maxData / linearStep) * linearStep;
    const x = val => left + (val - minX) / (maxX - minX) * (right - left);
    const y = val => bottom - ((log ? Math.log10(val) : val) - minY) / (maxY - minY) * (bottom - top);
    geometry = {width, left, right, top, bottom, minX, maxX, x, y, log};
    svg.setAttribute("viewBox", '0 0 ' + width + ' ' + height);
    svg.replaceChildren(svgEl("desc", {}, "Step histograms for both beam energies. Use the bin inspector to read the displayed values and bounds."));
    const axis = svgEl("g", {"aria-hidden":"true"});
    const logStep = Math.max(1, Math.ceil((maxY - minY) / 5));
    const yTicks = log ? Array.from({length: Math.floor((maxY - minY) / logStep) + 1}, (_, i) => 10 ** (maxY - i * logStep)) :
      Array.from({length: Math.round(maxY / linearStep) + 1}, (_, i) => i * linearStep);
    for (const tick of yTicks) {
      const yp = y(tick);
      axis.append(svgEl("line", {class:"grid", x1:left, x2:right, y1:yp, y2:yp}), svgEl("text", {x:left - 7, y:yp, "text-anchor":"end", dy:".32em"}, tick > 0 && tick < .01 ? tick.toExponential(0) : number(tick, 3)));
    }
    const stride = stepSize(maxX - minX, width < 520 ? 4 : 7);
    const xStep = isTracks() ? Math.ceil(stride) : stride;
    for (let tick = Math.ceil(minX / xStep) * xStep; tick <= maxX + xStep * 1e-6; tick += xStep) {
      const xp = x(tick);
      axis.append(svgEl("line", {class:"tick", x1:xp,x2:xp,y1:bottom,y2:bottom + 4}),svgEl("text", {x:xp, y:bottom + 18,"text-anchor":"middle"}, number(Math.abs(tick) < xStep * 1e-7 ? 0 : tick, 4)));
    }
    axis.append(svgEl("line", {class:"tick", x1:left,x2:right,y1:bottom,y2:bottom}));
    svg.append(svgEl("rect",{id:"active-strip",y:top,height:bottom - top}), axis);
    for (const beam of beams) {
      let path = "", connected = false;
      for (const [lo, hi, val] of series[beam]) {
        if (log && val <= 0) { connected = false; continue; }
        path += (connected ? "L" : "M") + x(lo) + " " + y(val) + "H" + x(hi) + " ";
        connected = true;
      }
      svg.append(svgEl("path",{class:"trace beam" + beam, d:path}));
    }
    svg.append(svgEl("g", {id:"whiskers"}));
    drawSelection();
  }

  function drawSelection() {
    if (!geometry) return;
    const {x, y, log, bottom} = geometry;
    const [lo, hi] = series["130"][index];
    const strip = byId("active-strip"), whiskers = byId("whiskers");
    strip.setAttribute("x",x(lo));
    strip.setAttribute("width",Math.max(x(hi) - x(lo), 1));
    whiskers.replaceChildren();
    for (const [k, beam] of beams.entries()) {
      const [, , val, lower, upper] = series[beam][index];
      const xp = x(lo + (hi - lo) * (k === 0 ? .32 : .68));
      const g = svgEl("g",{class:"whisker beam" + beam});
      if (lower !== null && upper !== null && (!log || upper > 0)) {
        const upperY = y(upper), lowerY = log && lower === 0 ? bottom : y(lower);
        const path = 'M' + xp + ' ' + upperY + 'V' + lowerY + ' M' + (xp - 3) + ' ' + upperY + 'h6' + (log && lower === 0 ? '' : ' M' + (xp - 3) + ' ' + lowerY + 'h6');
        g.append(svgEl("path",{d:path}));
      }
      if (!log || val > 0) g.append(svgEl("circle",{cx:xp,cy:y(val),r:3.3}));
      whiskers.append(g);
    }
  }

  function selectBin(idx) {
    index = Math.max(0, Math.min(series["130"].length - 1, idx));
    const [lo, hi] = series["130"][index], unit = meta[observable.value][2];
    const range = isTracks() ? number((lo + hi) / 2) + ' tracks' : number(lo) + ' – ' + number(hi) + (unit ? ' ' + unit : '');
    byId("bin-range").textContent = range;
    byId("bin-num").textContent = 'Bin ' + (index + 1) + ' of ' + series["130"].length;
    byId("previous").disabled = index === 0;
    byId("next").disabled = index === series["130"].length - 1;
    slider.value = index + 1;
    let summary = range;
    const output = byId("bin-values");
    output.replaceChildren();
    for (const beam of beams) {
      const [, , central, lower, upper] = series[beam][index];
      const row = document.createElement("div");
      row.className = "value-row";
      const title = document.createElement("div"), val = document.createElement("strong"), bounds = document.createElement("span");
      title.className = "beam-label"; title.innerHTML = '<i class="line e' + beam + '"></i>9 × ' + beam + ' GeV';
      val.textContent = number(central); val.title = String(central);
      const boundsText = lower === null || upper === null ? "Bounds unavailable" : "68.3% bounds  " + number(lower) + " – " + number(upper);
      bounds.textContent = boundsText;
      if (lower !== null && upper !== null) bounds.title = String(lower) + " – " + String(upper);
      summary += ', 9 by ' + beam + ': ' + number(central) + ' ' + yUnit() + ', ' + boundsText;
      row.append(title, val, bounds); output.append(row);
    }
    slider.setAttribute("aria-valuetext", summary);
    drawSelection();
    clearTimeout(urlTimer);
    urlTimer = setTimeout(updateHash, 100);
  }

  function pointerBin(event) {
    if (!geometry) return;
    const rect = svg.getBoundingClientRect(), {width,left,right,top,bottom,minX,maxX} = geometry;
    const px = (event.clientX - rect.left) * width / rect.width, py = event.clientY - rect.top;
    if (px < left || px > right || py < top || py > bottom) return;
    const point = minX + (px - left) / (right - left) * (maxX - minX);
    const found = series["130"].findIndex(row => point < row[1]);
    const newIndex = found < 0 ? series["130"].length - 1 : found;
    if (newIndex !== index) selectBin(newIndex);
  }

  function exportCsv() {
    const columns = ["state","figure","observable","electron_GeV","proton_GeV","luminosity_invfb","bin","bin_low","bin_high","x_unit","display_value","display_lower","display_upper","y_unit","bounds_status"];
    const rows = [columns];
    for (const beam of beams) {
      series[beam].forEach(([lo,hi,val,lower,upper], i) => {
        rows.push([state.value, figureNumber(), observable.value, 9, beam, 10, i + 1, lo, hi, meta[observable.value][2] || "dimensionless", val, lower, upper, yUnit(), lower === null || upper === null ? "unavailable" : "available"]);
      });
    }
    const blob = new Blob([rows.map(row => row.join(",")).join("\r\n") + "\r\n"], {type:"text/csv;charset=utf-8"});
    const url = URL.createObjectURL(blob), link = document.createElement("a");
    link.href = url;
    link.download = ['AN-fig' + figureNumber(), state.value, observable.value, 'both-beams.csv'].join("-");
    document.body.append(link); link.click(); link.remove();
    setTimeout(() => URL.revokeObjectURL(url), 1000);
  }

  function validatedData(candidate) {
    const invalid = () => { throw new Error("This file does not match the beam-bin export. Choose beams.json from the replay."); };
    const exactKeys = (obj, expected) => obj !== null && typeof obj === "object" && !Array.isArray(obj) &&
      Object.keys(obj).length === expected.length && expected.every(key => Object.hasOwn(obj, key));
    const result = {};
    if (!exactKeys(candidate, Object.keys(catalog))) invalid();
    for (const [stateKey, groups] of Object.entries(catalog)) {
      const incomingGroups = candidate[stateKey];
      if (!exactKeys(incomingGroups, Object.keys(groups))) invalid();
      result[stateKey] = {};
      for (const [groupKey, trusted] of Object.entries(groups)) {
        const incoming = incomingGroups[groupKey];
        if (!incoming || !exactKeys(incoming.observables, Object.keys(trusted.observables))) invalid();
        if (incoming.sha256 !== undefined && incoming.sha256 !== trusted.sha256) invalid();
        const validated = {};
        for (const key of Object.keys(trusted.observables)) {
          const values = incoming.observables[key];
          if (!exactKeys(values, beams)) invalid();
          for (const beam of beams) {
            const bins = values[beam];
            if (!Array.isArray(bins) || !bins.length || bins.length > 2000 || bins.length !== values["130"].length) invalid();
            bins.forEach((row, i) => {
              if (!Array.isArray(row) || row.length !== 5) invalid();
              const [lo, hi, val, low, high] = row;
              if (![lo, hi, val].every(Number.isFinite) || lo >= hi || val < 0 ||
                  (i && lo !== bins[i - 1][1]) || lo !== values["130"][i][0] || hi !== values["130"][i][1]) invalid();
              if (low === null || high === null) {
                if (key !== "reco_t_proj" || low !== null || high !== null) invalid();
              } else if (![low, high].every(Number.isFinite) || low < 0 || high < low) invalid();
            });
          }
          validated[key] = values;
        }
        result[stateKey][groupKey] = {...trusted, observables: validated};
      }
    }
    return result;
  }

  let importRequest = 0;
  byId("bins-file").addEventListener("change", async event => {
    const file = event.target.files[0], request = ++importRequest, status = byId("import-status");
    if (!file) return;
    try {
      if (file.size > 10 * 1024 * 1024) throw new Error("Choose a beam-bin JSON file smaller than 10 MB.");
      const contents = await file.text();
      if (request !== importRequest) return;
      const checked = validatedData(JSON.parse(contents));
      data = checked;
      status.textContent = "Loaded all three states and both beam energies.";
      showObservable(pendingBin);
    } catch (err) {
      if (request === importRequest) status.textContent = err instanceof SyntaxError ? "The file contains invalid JSON. Choose beams.json from the replay." : err.message;
    }
  });

  state.addEventListener("change", () => {
    fillSelect(group, Object.keys(data[state.value]), families, group.value);
    fillSelect(observable, Object.keys(currentFigure().observables), meta, observable.value);
    showObservable();
  });
  group.addEventListener("change", () => {
    fillSelect(observable, Object.keys(currentFigure().observables), meta, observable.value);
    showObservable();
  });
  observable.addEventListener("change", () => showObservable());
  scale.addEventListener("change", () => {drawChart(); updateHash();});
  slider.addEventListener("input", () => selectBin(Number(slider.value) - 1));
  byId("previous").addEventListener("click", () => selectBin(index - 1));
  byId("next").addEventListener("click", () => selectBin(index + 1));
  byId("note-shortcut").addEventListener("click", () => {
    const jump = byId("note-jump");
    if (jump.getAttribute("href")) jump.click();
    else { byId("note-file").scrollIntoView({block:"center"}); byId("note-file").focus({preventScroll:true}); }
  });
  byId("csv").addEventListener("click", exportCsv);
  svg.addEventListener("pointermove", event => {if (event.pointerType !== "touch") pointerBin(event);});
  svg.addEventListener("pointerdown", pointerBin);
  window.addEventListener("hashchange", restore);
  window.notePage.attach({pick:byId("note-file"),link:byId("note-jump"),page:notePageNumber});
  restore();
  if (window.ResizeObserver) new ResizeObserver(drawChart).observe(svg);
  else window.addEventListener("resize", drawChart);
})();
