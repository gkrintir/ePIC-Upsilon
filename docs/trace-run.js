/* Read existing run products on the user's device. The picker does not upload files. */
(() => {
  "use strict";
  const files = new Map();
  const urls = [];
  let context = null;
  let view = null;
  let generation = 0;
  let drawing = 0;
  let panel, summary, instructions, result, picker;

  const node = (tag, text) => {
    const item = document.createElement(tag);
    if (text !== undefined) item.textContent = text;
    return item;
  };
  function safePath(path) {
    return typeof path === "string" && path.length < 1024 && !path.includes("\\") &&
      !path.includes(":") && !path.startsWith("/") &&
      path.split("/").every(part => part && part !== "." && part !== "..");
  }
  function disposeLinks() {
    urls.splice(0).forEach(url => URL.revokeObjectURL(url));
  }
  async function readJSON(path) {
    const file = files.get(path);
    if (!file || file.size > 5e6) throw new Error("Missing or unreadable " + path);
    return JSON.parse(await file.text());
  }
  async function readDirectory() {
    const revision = ++generation;
    context = null;
    disposeLinks();
    files.clear();
    for (const file of picker.files) {
      const path = file.webkitRelativePath.split("/").slice(1).join("/");
      if (safePath(path)) files.set(path, file);
    }
    if (!files.size) { draw(); return; }
    try {
      let next;
      if (files.has("paper_replay.json")) {
        const record = await readJSON("paper_replay.json");
        if (!safePath(record.figure) || !record.figure.endsWith(".pdf") || !record.cells_checked ||
            typeof record.spectra?.all_pass !== "boolean" || !Number.isInteger(record.cells_checked.pt_y) || record.cells_checked.pt_y < 1 ||
            !Number.isInteger(record.cells_checked.xv_w) || record.cells_checked.xv_w < 1 ||
            !Number.isInteger(record.mass_values_and_errors_checked) || record.mass_values_and_errors_checked < 1) {
          throw new Error("Invalid paper_replay.json. Choose a completed paper replay.");
        }
        next = {kind:"paper",record,fitDir:record.figure.split("/").slice(0,-1).join("/")};
      } else {
        const directory = ["artifacts/note/compare_9x130_9x275/", ""].find(prefix => files.has(prefix+"kinematic_selection_audit.json"));
        if (directory === undefined) throw new Error("Choose a directory containing paper_replay.json, or a two-energy state run with its kinematic selection audit.");
        const audit = await readJSON(directory+"kinematic_selection_audit.json");
        if (audit.proton_selection?.mode !== "inclusive" || audit.stat_mode !== "projected" ||
            audit.projection_luminosity_invfb !== 10 || !["1S","2S","3S"].includes(audit.state) ||
            !["9x130 GeV","9x275 GeV"].every(beam => audit.samples?.some(sample => sample.label === beam))) {
          throw new Error("This audit does not describe the inclusive two-energy comparison at 10 fb⁻¹.");
        }
        next = {kind:"beams",state:audit.state,directory};
      }
      if (revision === generation) context = next;
    } catch (err) {
      if (revision === generation) context = {kind:"error",message: err instanceof SyntaxError ? "The audit is not valid JSON." : err.message};
    }
    if (revision === generation) draw();
  }
  function resolvePath(template) {
    if (template === "@figure") return context.record?.figure;
    return template.replace("@fitDir/", (context.fitDir || "")+"/")
      .replace("@compare/", context.directory || "").replaceAll("@state", (context.state || "").toLowerCase());
  }
  function appendFileList(paths) {
    const list=node("ul"); list.className="trace-file-list";
    for (const template of paths) {
      const path=resolvePath(template);
      const row=node("li"), file=files.get(path);
      const name=node("code",path);
      if (file && file.size>0) {
        const ext=path.split(".").pop();
        if (["json","csv","md","pdf"].includes(ext)) {
          const link=node("a"); link.append(name);
          const type=ext==="pdf" ? "application/pdf" : "text/plain;charset=utf-8";
          const url=URL.createObjectURL(file.slice(0,file.size,type)); urls.push(url);
          link.href=url;link.target="_blank";link.rel="noopener";
          row.append(link,document.createTextNode(" · open from this run"));
        } else row.append(name,document.createTextNode(" · present"));
      } else row.append(name,document.createTextNode(" · missing"));
      list.append(row);
    }
    result.append(list);
  }
  function openInspector(key, message) {
    if (location.origin==="null") {
      message.textContent="Open figures.html and choose this replay folder to inspect the values.";
      return;
    }
    let child;
    const target = new URL("figures.html?from-run=1#"+key, location.href);
    const cleanup=() => {window.removeEventListener("message",ready);clearTimeout(timeout);};
    async function ready(event) {
      if (event.source!==child || event.origin!==location.origin || event.data?.type!=="epic-run-view-ready") return;
      cleanup();
      try {
        const payload=await window.paperDirectory.read(selected);
        child.postMessage({type:"epic-run-result",payload},location.origin);
        message.textContent="Opened this run in the figure inspector.";
      } catch(error) {
        message.textContent=error.message;
        if(child && !child.closed) child.close();
      }
    }
    const selected=new Map(files);
    const timeout=setTimeout(cleanup,45000);
    window.addEventListener("message",ready);
    child=window.open(target.href,"_blank");
    if(!child) {cleanup();message.textContent="The browser blocked the new tab. Open the figure and choose this folder there.";}
  }
  function draw() {
    const renderNumber=++drawing;
    if (!view || !panel) return;
    disposeLinks();
    result.replaceChildren();
    const [key, id]=view;
    const route=window.UPSiLON_TRACE.routes[key], step=route.steps.find(s=>s.id===id);
    const isBeam=key==="beams";
    instructions.textContent=isBeam ?
      "Choose a two-energy state run or its compare_9x130_9x275 directory. Only local files are opened." :
      "Choose the output directory containing paper_replay.json. Only local files are opened.";
    summary.textContent=context && context.kind!=="error" ? "Inspect files from your run · folder selected" : "Inspect files from your run";
    if (!context) return;
    if (context.kind==="error") {result.append(node("p",context.message));return;}
    if ((context.kind==="beams")!==isBeam) {
      result.append(node("p",isBeam ? "This folder is a paper replay. Select an inclusive two-energy run for this result." :
        "This folder is a two-energy comparison. Select a paper replay for this result."));return;
    }
    if (id==="reconstruction") {
      result.append(node("p","A downstream run cannot identify the historical producer and merge jobs. See the production records listed above."));return;
    }
    if (context.kind==="paper" && (id==="plot" || id==="publication")) {
      const r=context.record;
      let message;
      if (key==="mass") message="The replay records "+r.mass_values_and_errors_checked+" mass-fit values and errors checked.";
      else if (key==="efficiency") message="The replay records "+r.cells_checked.pt_y+" pT–rapidity cells and "+r.cells_checked.xv_w+" xᵥ–W cells checked.";
      else message=r.spectra.all_pass===true ? "The replay reports agreement with the reference values and validation statuses. Check the per-observable status below." : "The replay does not report a successful spectral comparison.";
      result.append(node("p",message));
    }
    if (context.kind==="beams" && id==="publication") result.append(node("p","Local plots are for Υ("+context.state+"). The note uses this state in Figs. 10–11; Fig. 12 uses 1S."));
    if (context.kind==="paper" && ["plot","publication"].includes(id)) {
      const button=node("button","Inspect the numbers from this folder");
      const message=node("p");
      button.type="button";
      button.addEventListener("click",()=>openInspector(key,message));
      result.append(button,message);
    }
    if (window.inspectTraceEvidence) {
      const notes=node("div");result.append(notes);
      window.inspectTraceEvidence({route:key,stage:id,kind:context.kind,files:new Map(files),
                                  record:context.record,referenceRevision:window.UPSiLON_TRACE.referenceRevision})
        .then(lines=>{if(renderNumber===drawing) notes.replaceChildren(...lines.map(line=>node("p",line)));})
        .catch(()=>{if(renderNumber===drawing) notes.textContent="The local run record could not be compared with the reference.";});
    }
    if (step.files?.length) appendFileList(step.files.filter(path=>context.kind!=="beams" || context.state==="1S" || !path.endsWith("t_estimators_projective_transverse.pdf")));
    else result.append(node("p","This step has no file in the selected downstream directory."));
  }
  window.showTraceRun = (key,id) => {
    const host=document.getElementById("trace-run");
    if (!host) return;
    if (!panel || !host.contains(panel)) {
      host.replaceChildren();
      panel=node("details");panel.className="run-details";
      summary=node("summary","Inspect files from your run");instructions=node("p");
      picker=node("input"); picker.type="file";picker.multiple=true;picker.setAttribute("webkitdirectory","");
      picker.id="trace-folder"; const label=node("label","Local run folder");label.htmlFor=picker.id;
      result=node("div");result.setAttribute("aria-live","polite");
      picker.addEventListener("change",readDirectory);
      panel.append(summary,instructions,label,picker,result);host.append(panel);
    }
    view=[key,id]; draw();
  };
})();
