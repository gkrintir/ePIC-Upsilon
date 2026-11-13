(() => {
  "use strict";
  const local = window.PAPER, images = {...local.images};
  let data = local, view = "mass", state = "1S", quantity = "pt2", map = "pt_y", bin = 0, xi = 0, yi = 0, csvOut;
  let fullFigure = false;
  const $ = id => document.getElementById(id);
  const states = ["1S", "2S", "3S"];
  const maps = ["pt_y","xv_w"];
  const cellFields = ["map","x_low","x_high","y_low","y_high","denominator","selected_before_proton","selected"];
  const pointFields = ["quantity","state","bin_low","bin_high","central","lower","upper","status"];
  const papers = {
    mass: {tag:"Signal extraction · Fig. 7 (continued)", title:"Three-state dielectron mass fit", lead:"The three fitted peaks determine the expected Υ(1S), Υ(2S) and Υ(3S) yields. The fit shares a mass shift and width scale across the states.", note:"Figure 7", page:27, about:"Figure 7 shows the three fitted Υ peaks. Table 8 gives the expected yield and Gaussian core width of each peak.", image:"mass", caption:"Paper panel (a) · expected events per 0.05 GeV vs. reconstructed dielectron mass [GeV].", inspector:"Fitted yields and core widths"},
    efficiency: {tag:"Detector corrections · Fig. 13", title:"Proton-tagged Υ(1S) efficiency", lead:"For each truth bin, the efficiency is the fraction of eligible simulated decays that pass reconstruction, selection and the recoil-proton tag.", note:"Figure 13", page:46, about:"Figure 13 shows the fraction of eligible Υ(1S) decays that pass the full selection, including the proton tag.",image:"mass",caption:"Paper panels (b,c) · Υ(1S) reconstruction and selection efficiency [%].",inspector:"Efficiency by truth bin"},
    w: {tag:"Energy dependence · Fig. 14", title:"Cross section versus photon–proton energy", lead:"The points give the adaptive-bin γp cross-section diagnostics for each Υ state. The horizontal bars span each displayed W bin.", note:"Figure 14",page:49,about:"Figure 14 shows the rematched adaptive-bin W projections. Their interval-coverage and closure checks are recorded with the numerical results.",image:"w",caption:"Bin-averaged σγp [pb] vs. Wγp [GeV].",inspector:"Energy-spectrum bins"},
    differential: {tag:"Differential projections · Fig. 15",title:"Dielectron pT² and laboratory rapidity",lead:"Panel (a) shows the Υ(1S) pT² projection and the available Υ(2S) diagnostic. Panel (b) uses the integrated photon flux of each state to normalize its rapidity spectrum.",note:"Figure 15",page:50,about:"Figure 15 shows the pT² spectra and slope profiles in panel (a), with the three rapidity spectra in panel (b). Each state retains its recorded validation status.",image:"differential",caption:"Panel (a) · dσγp/dpT² [pb/GeV²]. Panel (b) · flux-normalized ep yield density [pb per unit rapidity].",inspector:"Differential-spectrum bins"}
  };
  const make = (tag, content, cls) => {
    const node = document.createElement(tag);
    if (content !== undefined) node.textContent = content;
    if (cls) node.className = cls;
    return node;
  };
  const fmt = (n, digits=4) => n.toLocaleString("en-US",{maximumSignificantDigits:digits});
  const count = n => n.toLocaleString("en-US");
  const validN = Number.isFinite;
  function setRoute() {
    const params = new URLSearchParams();
    if (view === "efficiency") {
      params.set("map", map); if (data.efficiency.cells.length) {params.set("x",xi);params.set("y",yi);}
    }
    if (view === "w" || view === "differential") {
      params.set("state",state); if (view==="differential") params.set("observable",quantity); if(data.spectra.length)params.set("bin",bin+1);
    }
    if(fullFigure)params.set("panel","all");
    const suffix=params.toString();
    $("trace-link").href="./index.html#"+view+"/plot/review"+(suffix ? "?"+suffix : "");
    const hash="#"+view+(suffix?"?"+suffix:"");
    if (location.hash!==hash) {
      try {history.replaceState(null,"",hash);} catch {location.replace(hash);}
    }
  }
  function positiveIndex(x) { const n=Number(x); return Number.isSafeInteger(n)&&n>=0?n:0; }
  function restoreRoute() {
    const [name, qs=""] = location.hash.slice(1).split("?");
    if (name === "figure-content") {$("figure-content").focus();return;}
    view = Object.hasOwn(papers,name)?name:"mass";
    const args=new URLSearchParams(qs);
    fullFigure=args.get("panel")==="all";
    state=states.includes(args.get("state"))?args.get("state"):"1S";
    quantity=args.get("observable")==="y"?"y":"pt2";
    map=args.get("map")==="xv_w"?"xv_w":"pt_y";
    xi=positiveIndex(args.get("x")); yi=positiveIndex(args.get("y")); bin=Math.max(0,positiveIndex(args.get("bin"))-1);
    render();
  }
  function render() {
    const p=papers[view];
    document.title=p.title+" · ePIC Υ Atlas";
    for(const tab of document.querySelectorAll("[data-tab]")) {
      if(tab.dataset.tab===view) tab.setAttribute("aria-current","page"); else tab.removeAttribute("aria-current");
    }
    $("kicker").textContent=p.tag;$("figure-title").textContent=p.title;$("lead").textContent=p.lead;
    $("figure-caption").textContent=p.caption;
    $("figure-image").src=images[p.image];$("figure-image").alt=p.title+". "+p.caption;
    $("enlarge").href=images[p.image];$("enlarge").setAttribute("aria-label","Enlarge "+p.title);
    $("open-image").href=images[p.image];
    $("source-pdf").hidden=!local.pdfs?.[p.image];
    if(local.pdfs?.[p.image]) $("source-pdf").href=local.pdfs[p.image];
    $("note-description").textContent=p.about;
    $("note-link").textContent=p.note+" · Plot (PDF page "+p.page+")";
    if (window.ATLAS_HAS_LOCAL_NOTE===true) $("local-label").textContent="Use the bundled note below or choose a different copy.";
    $("inspector-title").textContent=p.inspector;
    $("import-box").hidden=false;
    const matchesReference = ["mass", "efficiency", "spectra"].every(key => JSON.stringify(data[key]) === JSON.stringify(local[key]));
    $("reference-artwork").hidden=!matchesReference;
    document.querySelector(".figure-grid").style.gridTemplateColumns=matchesReference ? "" : "1fr";
    $("reference-note").hidden=matchesReference;
    $("reference-note").textContent="The numerical inspector uses your selected local results. Reload this page to restore the bundled rematched figures. The note links below refer to the November 2026 reference.";
    window.notePage.refresh();
    renderResults();
  }
  function appendTable(columns,rows) {
    const wrap=make("div",undefined,"result-table-wrap"),t=make("table",undefined,"result-table"),thead=make("thead"),tr=make("tr");
    for(const col of columns) {const th=make("th",col);th.scope="col";tr.append(th);}
    thead.append(tr);t.append(thead);
    const tbody=make("tbody");
    for(const cells of rows) {const row=make("tr"); for(const text of cells)row.append(make("td",text));tbody.append(row);}
    t.append(tbody);wrap.append(t);$("results").append(wrap);
    return tbody;
  }
  function reportCsv(name,fields,rows) {
    csvOut={name,fields,rows};
    $("csv").hidden=!rows.length;
  }
  function note(text,cls="notes") {$("results").append(make("p",text,cls));}
  function fitText(values) {
    const [v,minus,plus]=values;
    return minus===null||plus===null ? v.toFixed(1) : v.toFixed(1)+"  (−"+minus.toFixed(1)+" / +"+plus.toFixed(1)+")";
  }
  function massView() {
    if (!data.mass) {
      note("Open your run’s paper.json below to see the fitted yields and uncertainties for all three states.");
      reportCsv("empty",[],[]);return;
    }
    appendTable(["State","Fitted expected events","Gaussian core width [MeV]"],states.map(s=>["Υ("+s+")",fitText(data.mass[s].yield),fitText(data.mass[s].width)]));
    const hasErrors=data.mass["1S"].yield[1]!==null;
    note(hasErrors?"Parentheses give the negative and positive 68.27% conditional profile-likelihood errors. The fit fixes the calibrated tails and relative widths, so the width errors share the same fitted scale. These signal-only intervals require separate checks of shape mismatch, backgrounds and low-count coverage.":"The table gives the rounded central values. For the conditional fit errors, open Table 8 or your run’s paper.json.");
    reportCsv(hasErrors?"mass_fit":"mass_fit_summary",["state","yield_events","yield_err_minus","yield_err_plus","core_width_MeV","width_err_minus_MeV","width_err_plus_MeV"],states.map(s=>[s,...data.mass[s].yield,...data.mass[s].width]));
  }
  function addSelect(box,id,label,choices,value,onChange) {
    const lab=make("label",label),selector=make("select");selector.id=id;lab.htmlFor=id;
    for(const [k,name] of choices)selector.add(new Option(name,k));selector.value=value;
    selector.addEventListener("change",()=>onChange(selector.value));lab.append(selector);box.append(lab);
    return selector;
  }
  function redrawn(selector) {renderResults();$(selector)?.focus();}
  function simpleStat(box,value,label) { const item=make("div");item.append(make("b",value),make("small",label));box.append(item); }
  function svgNode(tag,attrs,text) {
    const node=document.createElementNS("http://www.w3.org/2000/svg",tag);
    for(const [key,val] of Object.entries(attrs))node.setAttribute(key,val);
    if(text!==undefined)node.textContent=text;return node;
  }
  function color(t) {
    return "rgb("+[232-232*t,244-154*t,240-144*t].map(Math.round).join(",")+")";
  }
  function efficiencyView() {
    const e=data.efficiency, stats=make("div",undefined,"stats");
    if (!Number.isFinite(e.eligible)) {
      note("Open your run’s paper.json below to inspect the counts before and after proton tagging in each cell.");
      reportCsv("empty",[],[]);return;
    }
    simpleStat(stats,count(e.eligible),"Eligible simulated decays");
    simpleStat(stats,(100*e.pretag/e.eligible).toFixed(2)+"%","Selected before proton tag");
    simpleStat(stats,(100*e.tagged/e.eligible).toFixed(2)+"%","Selected with proton tag");$("results").append(stats);
    const stageBox = make("div", undefined, "filters");
    const rate = make("output");
    rate.id = "eff-rate";
    rate.setAttribute("aria-live", "polite");
    const stages = {
      tag_of_all: [e.tagged, e.eligible, "Tag and all preceding selections / eligible decays"],
      tag_of_pretag: [e.tagged, e.pretag, "Proton tag / events selected before the tag"],
      pretag_of_all: [e.pretag, e.eligible, "Selection before proton tagging / eligible decays"]
    };
    function showRate(key) {
      const [selected, total] = stages[key];
      rate.textContent = count(selected) + " / " + count(total) + " = " + (100*selected/total).toFixed(2) + "%";
    }
    addSelect(stageBox, "eff-stage", "Compare selection stages", Object.entries(stages).map(([key, item]) => [key, item[2]]), "tag_of_all", showRate);
    stageBox.append(rate);
    $("results").append(stageBox);
    showRate("tag_of_all");
    if(!e.cells.length) {
      note("Open your run’s paper.json below to inspect the numerator and denominator of each cell.");
      reportCsv("efficiency_summary",["eligible","selected_before_proton","selected"],[[e.eligible,e.pretag,e.tagged]]);return;
    }
    const filters=make("div",undefined,"filters");$("results").append(filters);
    addSelect(filters,"map","Truth coordinates",[["pt_y","Υ transverse momentum × rapidity"],["xv_w","xᵥ × photon–proton energy"]],map,value=>{map=value;xi=0;yi=0;redrawn("map");});
    const rows=e.cells.filter(row=>row.map===map);
    reportCsv("efficiency_"+map,cellFields,rows.map(row=>cellFields.map(key=>row[key])));
    const xs=[...new Set(rows.map(r=>r.x_low))].sort((a,b)=>a-b),ys=[...new Set(rows.map(r=>r.y_low))].sort((a,b)=>a-b);
    const grid=new Map(rows.map(row=>[[row.x_low,row.y_low].join("|"),row]));
    const get=(x,y)=>grid.get([xs[x],ys[y]].join("|"));
    xi=Math.min(xi,xs.length-1);yi=Math.min(yi,ys.length-1);
    const max=Math.max(...rows.map(r=>r.denominator?100*r.selected/r.denominator:0));
    const pt=map==="pt_y",xTitle=pt?"True Υ pT [GeV]":"True xᵥ (log scale)",yTitle=pt?"True lab rapidity":"True Wγp [GeV]";
    const layout=make("div",undefined,"heat-layout"),figure=make("div"),side=make("div",undefined,"heat-side");layout.append(figure,side);$("results").append(layout);
    const legend=make("div",undefined,"map-legend");legend.append("Tagged efficiency   0%",make("i"),fmt(max,3)+"%",make("em"),"No eligible events");figure.append(legend);
    const svg=svgNode("svg",{id:"heatmap",viewBox:"0 0 560 370",tabindex:"0",role:"group","aria-label":"Efficiency map. Arrow keys select a cell. Horizontal and vertical bin selectors follow."});
    const left=62,right=545,top=12,bottom=322;
    const xx=val=>left+(pt?(val-xs[0])/(get(xs.length-1,0).x_high-xs[0]):Math.log(val/xs[0])/Math.log(get(xs.length-1,0).x_high/xs[0]))*(right-left);
    const yy=val=>bottom-(val-ys[0])/(get(0,ys.length-1).y_high-ys[0])*(bottom-top);
    const rects=[];
    for(let x=0;x<xs.length;x++){
      rects[x]=[];
      for(let y=0;y<ys.length;y++){
        const c=get(x,y),rect=svgNode("rect",{x:xx(c.x_low),y:yy(c.y_high),width:xx(c.x_high)-xx(c.x_low),height:yy(c.y_low)-yy(c.y_high),fill:c.denominator?color(max?100*c.selected/c.denominator/max:0):"#d5dadd",stroke:"#fff","stroke-width":".6","data-x":x,"data-y":y});svg.append(rect);rects[x].push(rect);
      }
    }
    function ticks(vals,last,position,horizontal) {
      const edges=[...vals,last];
      for(const i of [...new Set([0,Math.round(vals.length/4),Math.round(vals.length/2),Math.round(3*vals.length/4),vals.length])]){
        const val=edges[i], pos=position(val);
        svg.append(svgNode("text",horizontal?{x:pos,y:bottom+19,"text-anchor":"middle"}:{x:left-7,y:pos+4,"text-anchor":"end"},fmt(val,pt||!horizontal?3:2)));
      }
    }
    ticks(xs,get(xs.length-1,0).x_high,xx,true);ticks(ys,get(0,ys.length-1).y_high,yy,false);
    svg.append(svgNode("text",{x:(left+right)/2,y:bottom+42,"text-anchor":"middle"},xTitle),
      svgNode("text",{transform:"translate(15 "+((top+bottom)/2)+") rotate(-90)","text-anchor":"middle"},yTitle));
    const focus=svgNode("rect",{fill:"none",stroke:"#102b3a","stroke-width":2.5,"pointer-events":"none"});svg.append(focus);
    figure.append(svg,make("p","Select a cell, use the bin menus, or focus the map and press the arrow keys.","keyhint"));
    side.append(make("p","Selected cell","eyebrow"));const score=make("h3"),readout=make("p"),bounds=make("p"),pretag=make("p");side.append(score,readout,pretag,bounds);
    const choices=make("div",undefined,"filters");side.append(choices);
    const range=(a,b)=>fmt(a)+" – "+fmt(b);
    const xSel=addSelect(choices,"map-x",xTitle,xs.map((lo,i)=>[i,range(lo,get(i,0).x_high)]),xi,value=>select(Number(value),yi));
    const ySel=addSelect(choices,"map-y",yTitle,ys.map((lo,i)=>[i,range(lo,get(0,i).y_high)]),yi,value=>select(xi,Number(value)));
    function select(x,y) {
      xi=Math.max(0,Math.min(xs.length-1,x));yi=Math.max(0,Math.min(ys.length-1,y));
      const c=get(xi,yi),r=rects[xi][yi];
      for(const k of ["x","y","width","height"])focus.setAttribute(k,r.getAttribute(k));
      score.textContent=c.denominator?fmt(100*c.selected/c.denominator)+"% tagged efficiency":"No eligible events";
      readout.textContent=count(c.selected)+" tagged / "+count(c.denominator)+" eligible simulated decays.";
      pretag.textContent=count(c.selected_before_proton)+" selected before the tag.";
      bounds.textContent=xTitle+"  "+range(c.x_low,c.x_high)+"  ·  "+yTitle+"  "+range(c.y_low,c.y_high);
      xSel.value=xi;ySel.value=yi;svg.setAttribute("aria-label",score.textContent+". "+bounds.textContent+". Arrow keys select a cell.");setRoute();
    }
    svg.addEventListener("click",ev=>{const r=ev.target;if(r.hasAttribute("data-x"))select(Number(r.dataset.x),Number(r.dataset.y));});
    svg.addEventListener("keydown",ev=>{
      const delta={ArrowLeft:[-1,0],ArrowRight:[1,0],ArrowUp:[0,1],ArrowDown:[0,-1]}[ev.key];
      if(delta){ev.preventDefault();select(xi+delta[0],yi+delta[1]);}
    });
    select(xi,yi);
    note(pt?"Both axes use true Υ kinematics. The integrated efficiency includes eligible events outside the displayed pT and rapidity range.":"xᵥ = (Q² + M₁²) / (Wγp² + Q² − mₚ²). The horizontal axis is logarithmic. This map includes the full eligible truth sample.");
  }
  function spectraView() {
    const q = view === "w" ? "W" : quantity;
    const filters = make("div", undefined, "filters");
    $("results").append(filters);
    if (view === "differential") {
      addSelect(filters, "quantity", "Observable", [["pt2", "Dielectron pT²"], ["y", "Laboratory rapidity"]],
        q, value => {quantity=value; bin=0; redrawn("quantity");});
    }
    addSelect(filters, "state", "Υ state", states.map(s => [s, "Υ("+s+")"]),
      state, value => {state=value; bin=0; redrawn("state");});
    const rows = data.spectra.filter(row => row.quantity === q && row.state === state);
    const availability = data.availability?.find(row => row.quantity === q && row.state === state);
    const slope = data.spectra.find(row => row.quantity === "slope" && row.state === state);
    const slopeStatus = data.slopes?.[state]?.status || slope?.status;
    const isW = q === "W", isPt = q === "pt2";
    if (isPt && !rows.length && (availability?.available === false || state !== "1S")) {
      note(availability?.reason || "The expected Υ("+state+") yield has no feasible pT² slope binning in this run. Its rapidity projection is available in the other observable view.", "status-note");
      reportCsv("empty", [], []);
      return;
    }
    if (isPt && slope) {
      const stat = make("div", undefined, "stats");
      const errors = " (−"+(slope.central-slope.lower).toFixed(2)+" / +"+(slope.upper-slope.central).toFixed(2)+") GeV⁻²";
      simpleStat(stat, slope.central.toFixed(2)+errors,
        "Υ("+state+") exponential slope b · "+fmt(slope.bin_low)+" < pT² < "+fmt(slope.bin_high)+" GeV²");
      $("results").append(stat);
    }
    const rangeUnit = isW ? "GeV" : isPt ? "GeV²" : "";
    const xLabel = isW ? "Wγp [GeV]" : isPt ? "Dielectron pT² [GeV²]" : "Lab rapidity";
    const yLabel = isW ? "σγp [pb]" : isPt ? "dσγp/dpT² [pb/GeV²]" : "(1/Φₙ) dσep/dy [pb]";
    if (rows.length) {
      bin = Math.min(bin, rows.length-1);
      const picker = addSelect(filters, "spectrum-bin", "Inspect bin",
        rows.map((r,i) => [i, fmt(r.bin_low)+" – "+fmt(r.bin_high)+(rangeUnit ? " "+rangeUnit : "")]),
        bin, value => select(Number(value)));
      const caption = make("p", undefined, "keyhint");
      caption.setAttribute("aria-live", "polite");
      $("results").append(caption);
      const tbody = appendTable([xLabel, yLabel, "Lower band", "Upper band", "Status"],
        rows.map(r => [fmt(r.bin_low)+" – "+fmt(r.bin_high), fmt(r.central), fmt(r.lower), fmt(r.upper),
          r.status === "PASS_CONDITIONAL" ? "Conditional pass" : "Diagnostic"]));
      Array.from(tbody.rows).forEach((row,i) => {
        const entry = row.cells[0], button = make("button", entry.textContent, "inline-link");
        button.type = "button";
        button.setAttribute("aria-label", "Inspect bin "+(i+1));
        button.onclick = () => {select(i); picker.value=i;};
        entry.replaceChildren(button);
      });
      function select(i) {
        bin = i;
        const row = rows[i];
        caption.textContent = "Selected bin "+(i+1)+" of "+rows.length+"  ·  "+yLabel+" = "+fmt(row.central,6)+
          "  ·  sampling band "+fmt(row.lower,6)+" – "+fmt(row.upper,6)+".";
        Array.from(tbody.rows).forEach((tr,j) => tr.classList.toggle("chosen-bin", j===i));
        setRoute();
      }
      select(bin);
      const csvRows = isPt ? [...rows, ...data.spectra.filter(r => r.quantity === "slope" && r.state === state)] : rows;
      reportCsv("spectrum_"+q+"_"+state, pointFields, csvRows.map(row => pointFields.map(key => row[key])));
    } else {
      note("Open your run’s paper.json below to inspect individual bins and download their full-precision values.");
      reportCsv("empty", [], []);
    }
    if (availability?.reason) {
      const interpretation = isPt
        ? "Slope errors are profile limits. The point bars are conditional toy sampling bands."
        : "Vertical bars describe conditional toy sampling bands in the fixed-response, signal-only model.";
      note(availability.reason+" "+interpretation, "status-note");
    } else if (isPt) {
      note(slopeStatus === "PASS_CONDITIONAL"
        ? "The Υ("+state+") pT² slope passed conditional validation in the fixed-response, signal-only tests. Its errors are profile limits. The point bars are toy sampling bands."
        : "The Υ("+state+") pT² result has diagnostic status. Its slope interval failed the injected-slope coverage check. The point bars are conditional toy sampling bands.", "status-note");
    } else {
      note("The "+(isW ? "adaptive-bin W" : "rapidity")+" projections have diagnostic status following the injected-spectrum coverage checks. Vertical bars describe conditional toy sampling bands in the fixed-response, signal-only model.", "status-note");
    }
    if (q === "y") note("The flux Φₙ is integrated across the full W range for each state. Positive laboratory rapidity follows the proton beam.");
  }
  function syncPanel() {
    const key = {mass:"mass",efficiency:map,w:state,differential:quantity}[view];
    const image = papers[view].image;
    window.showFigurePanel({view, image, src:images[image], selection:fullFigure ? "all" : key, onChange: selected => {
      fullFigure=selected==="all";
      if(!fullFigure) {
        if(view==="w") state=selected;
        if(view==="efficiency") {map=selected;xi=0;yi=0;}
        if(view==="differential") quantity=selected;
        bin=0;
      }
      renderResults();
      $("plot-panel").focus();
    }});
  }
  function renderResults() {
    $("results").replaceChildren();
    syncPanel();
    if(view==="mass")massView();else if(view==="efficiency")efficiencyView();else spectraView();
    setRoute();
  }
  function exportCsv() {
    if (!csvOut || !csvOut.rows.length)return;
    const lines=[csvOut.fields,...csvOut.rows].map(row=>row.join(","));
    const url=URL.createObjectURL(new Blob([lines.join("\r\n")+"\r\n"],{type:"text/csv;charset=utf-8"}));
    const a=make("a");a.href=url;a.download=csvOut.name+".csv";document.body.append(a);a.click();a.remove();setTimeout(()=>URL.revokeObjectURL(url),1000);
  }

  function validated(raw) {
    const fail=()=>{throw new Error("These files do not contain a complete set of paper results. Check the chosen run folder or exported paper.json.");};
    const record=x=>x!==null&&typeof x==="object"&&!Array.isArray(x);
    const n=x=>validN(x)&&Math.abs(x)<=1e15;
    const numArr=a=>Array.isArray(a)&&a.length===3&&a.every(x=>n(x)&&x>=0);
    const integer=x=>Number.isSafeInteger(x)&&x>=0;
    if(!record(raw)||!record(raw.mass)||!record(raw.efficiency)||!Array.isArray(raw.spectra)||!numArr(raw.slope))fail();
    const mass={};
    for(const s of states){const m=raw.mass[s];if(!record(m)||!numArr(m.yield)||!numArr(m.width))fail();mass[s]={yield:m.yield,width:m.width};}
    const e=raw.efficiency;
    if(![e.eligible,e.pretag,e.tagged].every(integer)||e.eligible<=0||e.pretag>e.eligible||e.tagged>e.pretag||!Array.isArray(e.cells)||!e.cells.length||e.cells.length>5000)fail();
    for(const c of e.cells){
      if(!record(c)||!maps.includes(c.map)||![c.x_low,c.x_high,c.y_low,c.y_high].every(n)||c.x_low>=c.x_high||c.y_low>=c.y_high||
          (c.map==="xv_w"&&c.x_low<=0)||![c.denominator,c.selected_before_proton,c.selected].every(integer)||
          c.selected_before_proton>c.denominator||c.selected>c.selected_before_proton)fail();
    }
    for(const m of maps){
      const rows=e.cells.filter(c=>c.map===m),xs=new Map(),ys=new Map(),pairs=new Set();
      for(const r of rows){
        const k=r.x_low+"|"+r.y_low;
        if(pairs.has(k)||(xs.has(r.x_low)&&xs.get(r.x_low)!==r.x_high)||(ys.has(r.y_low)&&ys.get(r.y_low)!==r.y_high))fail();
        pairs.add(k);xs.set(r.x_low,r.x_high);ys.set(r.y_low,r.y_high);
      }
      if(!rows.length||pairs.size!==xs.size*ys.size)fail();
      for(const d of [xs,ys]){
        const sorted=[...d.keys()].sort((a,b)=>a-b);if(sorted.some((v,i)=>i&&v!==d.get(sorted[i-1])))fail();
      }
    }
    if(!raw.spectra.length||raw.spectra.length>1000)fail();
    const seen=new Set();
    for(const r of raw.spectra){
      if(!record(r)||!states.includes(r.state)||!["pt2","W","y","slope"].includes(r.quantity)||
        ![r.bin_low,r.bin_high,r.central,r.lower,r.upper].every(n)||r.bin_low>=r.bin_high||r.lower<0||r.upper<r.lower||r.central<0||
        !["NOT_VALIDATED","PASS_CONDITIONAL"].includes(r.status))fail();
      const key=[r.quantity,r.state,r.bin_low,r.bin_high].join("|");if(seen.has(key))fail();seen.add(key);
    }
    for(const s of states) {
      for(const q of ["W","y"]) {
        const rows=raw.spectra.filter(r=>r.state===s&&r.quantity===q);
        if(!rows.length||rows.some((r,i)=>i>0&&r.bin_low!==rows[i-1].bin_high))fail();
      }
      const points=raw.spectra.filter(r=>r.state===s&&r.quantity==="pt2");
      const slopes=raw.spectra.filter(r=>r.state===s&&r.quantity==="slope");
      if(s==="1S"||points.length||slopes.length) {
        if(!points.length||slopes.length!==1||points.some((r,i)=>i>0&&r.bin_low!==points[i-1].bin_high))fail();
      }
    }
    const result={mass,efficiency:{eligible:e.eligible,pretag:e.pretag,tagged:e.tagged,cells:e.cells},spectra:raw.spectra,slope:raw.slope};
    if(raw.schema_version!==undefined) {
      if(raw.schema_version!==2||!Array.isArray(raw.availability)||raw.availability.length!==9||!record(raw.slopes))fail();
      const statuses=["PASS_CONDITIONAL","NOT_VALIDATED","NO_FEASIBLE_SLOPE_BINNING"];
      const validStatus=item=>record(item)&&typeof item.available==="boolean"&&statuses.includes(item.status)&&
        typeof item.reason==="string"&&Array.isArray(item.failed_checks)&&item.failed_checks.every(check=>typeof check==="string");
      const availableKeys=new Set();
      for(const item of raw.availability) {
        if(!validStatus(item)||!states.includes(item.state)||!["W","y","pt2"].includes(item.quantity))fail();
        const key=item.state+"|"+item.quantity;
        if(availableKeys.has(key))fail();
        availableKeys.add(key);
        const rows=raw.spectra.filter(r=>r.state===item.state&&r.quantity===item.quantity);
        if(item.available!==Boolean(rows.length)||rows.some(row=>row.status!==item.status))fail();
        if(!item.available&&item.status!=="NO_FEASIBLE_SLOPE_BINNING")fail();
      }
      const slopeFields=["value","error_low","error_high","lower","upper"];
      for(const state of states) {
        const item=raw.slopes[state];
        if(!validStatus(item))fail();
        const row=raw.spectra.find(point=>point.state===state&&point.quantity==="slope");
        if(item.available!==Boolean(row))fail();
        if(item.available) {
          if(!slopeFields.every(key=>n(item[key])&&item[key]>=0)||item.lower>item.upper||row.status!==item.status)fail();
        } else if(item.status!=="NO_FEASIBLE_SLOPE_BINNING"||slopeFields.some(key=>item[key]!==null))fail();
      }
      result.schema_version=2;
      result.availability=raw.availability;
      result.slopes=raw.slopes;
    }
    return result;
  }

  let req=0;
  $("data-file").addEventListener("change",async ev=>{
    const file=ev.target.files[0],seq=++req;if(!file)return;
    try{
      if(file.size>4*1024*1024)throw new Error("Choose a paper.json file smaller than 4 MB.");
      const text=await file.text();if(seq!==req)return;
      const result=validated(JSON.parse(text));data=result;
      $("load-status").textContent="Loaded fit intervals, efficiency maps and spectral bins.";render();
    }catch(err){if(seq===req)$("load-status").textContent=err instanceof SyntaxError?"The file contains invalid JSON. Choose the exported paper.json.":err.message;}
  });
  $("replay-folder").addEventListener("change", async () => {
    const seq=++req;
    try {
      const files=window.paperDirectory.indexDirectory($("replay-folder").files);
      if (!files.size) return;
      const imported=validated(await window.paperDirectory.read(files));
      if (seq!==req) return;
      data=imported;
      $("load-status").textContent="Loaded the fit, efficiency cells and spectral bins from your replay folder.";
      render();
    } catch (error) {
      if(seq===req) $("load-status").textContent=error.message;
    }
  });
  // The workflow can open a result with the folder its reader already selected.
  if (window.opener && location.origin!=="null" && new URLSearchParams(location.search).has("from-run")) {
    const accept = event => {
      if(event.source!==window.opener || event.origin!==location.origin || event.data?.type!=="epic-run-result") return;
      window.removeEventListener("message",accept);
      try {
        data=validated(event.data.payload);
        $("load-status").textContent="Inspecting the replay folder selected in the analysis workflow.";
        render();
      } catch(error) { $("load-status").textContent=error.message; }
    };
    window.addEventListener("message",accept);
    window.opener.postMessage({type:"epic-run-view-ready"},location.origin);
  }
  $("csv").addEventListener("click",exportCsv);
  window.notePage.attach({
    pick: $("local-note"),
    link: $("note-link"),
    page: () => papers[view].page,
    extraLinks: [{link: $("table-jump"), page: () => view === "mass" ? 24 : null}]
  });
  window.addEventListener("hashchange",restoreRoute);
  restoreRoute();
})();
