/* Adapt the files written by paper/replay_paper.py to the figure inspector. */
(() => {
  "use strict";
  const columns = {
    cells: ["map","x_low","x_high","y_low","y_high","denominator","selected_before_proton","selected"],
    points: ["quantity","state","bin_low","bin_high","central","lower","upper","status"]
  };
  function indexDirectory(selection) {
    const files=new Map();
    for (const file of selection) {
      const path=file.webkitRelativePath.split("/").slice(1).join("/");
      if (path && !path.split("/").some(part=>part===".." || part===".")) files.set(path,file);
    }
    return files;
  }
  async function textFile(files,path,limit=5000000) {
    const file=files.get(path);
    if (!file) throw new Error("Missing "+path+". Choose the replay output directory.");
    if (!file.size || file.size>limit) throw new Error("Unexpected file size for "+path);
    return file.text();
  }
  async function csv(files,path,fields,numeric) {
    // These two replay tables contain only controlled identifiers and numeric columns.
    const lines=(await textFile(files,path)).trim().split(/\r?\n/);
    if (lines.shift()!==fields.join(",")) throw new Error("Unexpected columns in "+path);
    return lines.map(line=>{
      const values=line.split(",");
      if (values.length!==fields.length) throw new Error("Malformed row in "+path);
      return Object.fromEntries(fields.map((name,i)=>{
        const value=values[i];
        if (numeric.includes(name) && (!value.trim() || !Number.isFinite(Number(value)))) throw new Error("Invalid "+name+" in "+path);
        return [name,numeric.includes(name)?Number(value):value];
      }));
    });
  }
  async function read(files) {
    const receipt=JSON.parse(await textFile(files,"paper_replay.json"));
    if (receipt.spectra?.all_pass !== true) throw new Error("The replay has no successful spectral comparison.");
    const figure=receipt.figure;
    if (typeof figure!=="string" || !/^reconstruction\/mass_and_efficiency_w_xv\/render-[\w.-]+\/upsilon_mass_and_efficiency_w_xv\.pdf$/.test(figure) || !files.has(figure)) {
      throw new Error("The figure recorded by paper_replay.json is missing or has an unexpected path.");
    }
    const fitPath=figure.slice(0,figure.lastIndexOf("/")+1)+"mass_fit_intervals.json";
    const fit=JSON.parse(await textFile(files,fitPath));
    const mass={};
    for (const state of ["1S","2S","3S"]) {
      const stateFit=fit.states?.[state];
      if (!stateFit) throw new Error("Missing fitted mass state "+state);
      mass[state]={};
      for (const [label,key] of [["yield","fitted_yield"],["width","fitted_sigma_core_MeV"]]) {
        mass[state][label]=["value","error_low","error_high"].map(k=>stateFit[key]?.[k]);
      }
    }
    const [cells,spectra]=await Promise.all([
      csv(files,"reconstruction/efficiency_cells.csv",columns.cells,columns.cells.slice(1)),
      csv(files,"spectra/displayed_points.csv",columns.points,columns.points.slice(2,7))
    ]);
    const slope=spectra.find(p=>p.quantity==="slope" && p.state==="1S");
    const totals=receipt.xv_w_totals;
    if (!slope || !totals || cells.filter(r=>r.map==="pt_y").length!==receipt.cells_checked?.pt_y ||
        cells.filter(r=>r.map==="xv_w").length!==receipt.cells_checked?.xv_w) throw new Error("The output tables disagree with the run record.");
    for (const key of ["denominator","selected_before_proton","selected"]) {
      if (cells.reduce((sum,c)=>sum+(c.map==="xv_w"?c[key]:0),0)!==totals[key]) throw new Error("The efficiency totals disagree with the run record.");
    }
    return {mass, efficiency:{eligible:totals.denominator,pretag:totals.selected_before_proton,tagged:totals.selected,cells},
      spectra,slope:[slope.central,slope.central-slope.lower,slope.upper-slope.central]};
  }
  window.paperDirectory={indexDirectory, read};
})();
