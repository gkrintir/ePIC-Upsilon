#!/usr/bin/env node
// Check the published numerical exports and the browser's supported import formats.
import assert from "node:assert/strict";
import {createHash} from "node:crypto";
import {readFileSync} from "node:fs";
import {resolve} from "node:path";
import {fileURLToPath} from "node:url";
import {runInNewContext} from "node:vm";

const root = fileURLToPath(new URL("../", import.meta.url));
const docs = resolve(root, "docs");
const readJson = path => JSON.parse(readFileSync(resolve(docs, path), "utf8"));
const paper = readJson("data/paper.json");
const beams = readJson("data/beams.json");

function checkAsset(record) {
  assert(!record.path.includes("..") && !record.path.startsWith("/"));
  const bytes = readFileSync(resolve(docs, record.path));
  assert.equal(createHash("sha256").update(bytes).digest("hex"), record.sha256, record.path);
  if (record.bytes !== undefined) assert.equal(bytes.length, record.bytes, record.path);
}

function checkDataScript(path, name, expected) {
  const context = {window: {}};
  runInNewContext(readFileSync(resolve(docs, path), "utf8"), context, {timeout: 1000});
  assert.deepEqual(JSON.parse(JSON.stringify(context.window[name])), expected, path);
}

function checkBrowserImports() {
  const source = readFileSync(resolve(docs, "app.js"), "utf8");
  const start = source.indexOf("  function validated(raw) {");
  const end = source.indexOf("\n  let req=0;", start);
  assert(start >= 0 && end > start, "Locate the browser result validator");
  const validate = runInNewContext("(" + source.slice(start, end).trim() + ")", {
    states: ["1S", "2S", "3S"], maps: ["pt_y", "xv_w"], validN: Number.isFinite
  }, {timeout: 1000});
  const current = validate(paper);
  assert.equal(current.slopes["2S"].status, "NOT_VALIDATED");
  assert.equal(current.slopes["3S"].available, false);
  assert(current.spectra.some(row => row.state === "2S" && row.quantity === "pt2"));

  // Historical exports contain only the original five top-level numerical fields.
  const legacy = Object.fromEntries(["mass", "efficiency", "spectra", "slope"].map(key => [key, paper[key]]));
  legacy.spectra = legacy.spectra.filter(row => row.state === "1S" || !["pt2", "slope"].includes(row.quantity));
  assert.equal(validate(legacy).spectra.length, legacy.spectra.length);

  const mutations = [
    value => {value.availability.find(row => row.state === "2S" && row.quantity === "pt2").available = false;},
    value => {value.spectra = value.spectra.filter(row => row.state !== "2S" || row.quantity !== "slope");},
    value => {value.slopes["3S"].value = 0;},
    value => {value.slopes["2S"].status = "PASS_CONDITIONAL";},
  ];
  for (const mutate of mutations) {
    const invalid = structuredClone(paper);
    mutate(invalid);
    assert.throws(() => validate(invalid), /complete set of paper results/);
  }
}

function main() {
  checkDataScript("data/paper.js", "PAPER", paper);
  checkDataScript("data/beams.js", "BEAM_FIGURES", beams);
  const provenance = readJson("data/paper-provenance.json");
  for (const record of provenance.exported_files) checkAsset(record);
  const beamProvenance = readJson("data/beams-provenance.json");
  for (const record of [...beamProvenance.figures, beamProvenance.data,
    beamProvenance.producer_counters.figure, beamProvenance.producer_counters.data]) checkAsset(record);
  assert.equal(provenance.verified_inputs.length, 6);
  assert.equal(paper.efficiency.cells.length, 528);
  assert.equal(paper.efficiency.eligible, 1126936);
  assert.equal(paper.efficiency.pretag, 301323);
  assert.equal(paper.efficiency.tagged, 198265);
  for (const key of ["denominator", "selected_before_proton", "selected"]) {
    const total = paper.efficiency.cells.filter(cell => cell.map === "xv_w")
      .reduce((sum, cell) => sum + cell[key], 0);
    const field = {denominator: "eligible", selected_before_proton: "pretag", selected: "tagged"}[key];
    assert.equal(total, paper.efficiency[field], key);
  }
  let beamRows = 0;
  for (const state of Object.values(beams)) {
    for (const figure of Object.values(state)) {
      checkAsset({path: figure.pdf, sha256: figure.sha256});
      for (const series of Object.values(figure.observables)) {
        assert(series["130"].length > 0);
        assert.deepEqual(series["130"].map(row => row.slice(0, 2)), series["275"].map(row => row.slice(0, 2)));
        for (const rows of Object.values(series)) {
          beamRows += rows.length;
          for (const [low, high, count, lower, upper] of rows) {
            assert([low, high, count].every(Number.isFinite) && low < high && count >= 0);
            assert(lower === null && upper === null || Number.isFinite(lower) && Number.isFinite(upper) && lower <= upper);
          }
        }
      }
    }
  }
  checkBrowserImports();
  console.log(`PASS Atlas exports, 528 efficiency cells, ${paper.spectra.length} spectrum rows, ${beamRows} beam rows and current/historical imports`);
}

try {main();}
catch (error) {console.error("FAIL " + error.message); process.exitCode = 1;}
