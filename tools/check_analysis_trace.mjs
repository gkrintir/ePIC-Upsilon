#!/usr/bin/env node
// From the repository root: node tools/check_analysis_trace.mjs [--replay PATH]
import assert from "node:assert/strict";
import {createHash} from "node:crypto";
import {readFileSync, statSync} from "node:fs";
import {dirname, isAbsolute, relative, resolve, sep} from "node:path";
import {fileURLToPath} from "node:url";
import {runInNewContext} from "node:vm";

const repo = fileURLToPath(new URL("../", import.meta.url));
const docs = resolve(repo, "docs");
const routeIds = ["mass", "efficiency", "w", "differential", "beams"];
const stageIds = ["reconstruction", "tuple", "analysis", "plot", "publication"];
const spectrumStageIds = ["reconstruction", "tuple", "response", "analysis", "plot", "publication"];
const publicRepo = "/gkrintir/ePIC-Upsilon/";

function nonempty(value, where) {
  assert(typeof value === "string" && value.trim(), where + " requires text");
}

function within(base, path) {
  const resolved = resolve(base, path);
  const rel = relative(base, resolved);
  assert(rel !== ".." && !rel.startsWith(".." + sep) && !isAbsolute(rel), "Path outside root: " + path);
  return resolved;
}

function fileAt(path) {
  assert(statSync(path).isFile(), "Expected a file: " + path);
}

function checkGuideAssets() {
  const html = readFileSync(resolve(docs, "index.html"), "utf8");
  const references = html.matchAll(/(?:src|href)=["'](\.\/[^"'?#]+\.(?:js|css))((?:\?[^"']*)?)["']/g);
  const checked = new Set();
  for (const [, path, query] of references) {
    const content = readFileSync(within(docs, path), "utf8").replaceAll("\r\n", "\n");
    const version = createHash("sha256").update(content).digest("hex").slice(0, 12);
    assert.equal(query, "?v=" + version, "In docs/index.html use " + path + "?v=" + version);
    checked.add(path);
  }
  assert(checked.has("./trace.js") && checked.has("./trace-data.js"), "Missing guide scripts in docs/index.html");
}

function checkSource({label, url}, name, pin) {
  nonempty(label, name + " link label");
  nonempty(url, name + " link URL");
  const link = new URL(url, "https://guide.invalid/index.html");
  assert(link.protocol === "https:", "Expected HTTPS or a site-relative link: " + url);
  if (link.hostname === "guide.invalid") {
    fileAt(within(docs, decodeURIComponent(link.pathname.slice(1))));
    return;
  }
  if (link.hostname !== "github.com" || !link.pathname.toLowerCase().startsWith(publicRepo.toLowerCase())) return;
  const [view, revision, ...parts] = link.pathname.slice(publicRepo.length).split("/");
  assert(["tree", "blob"].includes(view) && revision === pin, "Unpinned repository source: " + url);
  const path = within(repo, decodeURIComponent(parts.join("/")));
  if (view === "tree") {
    assert(statSync(path).isDirectory(), "Expected a repository directory: " + path);
    return;
  }
  fileAt(path);
  if (!link.hash.startsWith("#L")) return;
  const range = /^#L([1-9]\d*)(?:-L([1-9]\d*))?$/.exec(link.hash);
  assert(range, "Invalid line anchor: " + url);
  const source = readFileSync(path, "utf8");
  const count = source.split("\n").length - Number(source.endsWith("\n"));
  const start = Number(range[1]), end = Number(range[2] ?? range[1]);
  assert(start <= end && end <= count, "Line anchor outside repository file: " + url);
}

function replayPath(template, figure) {
  if (template === "@figure") return figure;
  if (template.startsWith("@fitDir/")) return dirname(figure) + "/" + template.slice("@fitDir/".length);
  assert(!template.includes("@"), "Unknown paper replay placeholder: " + template);
  return template;
}

function checkExplanation(text, data, where) {
  nonempty(text, where);
  const parts = data.referenceParts(text);
  for (const part of parts) {
    if (typeof part === "string") {
      assert(!part.includes("]["), where + " contains a malformed reference");
      continue;
    }
    if (Object.hasOwn(data.sources, part.target)) continue;
    const [route, stage, extra] = part.target.split("/");
    assert(!extra && data.routes[route]?.steps.some(step => step.id === stage),
      where + " refers to an unknown source or stage: " + part.target);
  }
  return parts.some(part => typeof part !== "string");
}

function main(args) {
  assert(args.length === 0 || (args.length === 2 && args[0] === "--replay" && args[1]),
    "Usage: node tools/check_analysis_trace.mjs [--replay PATH]");
  const replay = args.length ? resolve(args[1]) : null;
  const figure = replay ? JSON.parse(readFileSync(resolve(replay, "paper_replay.json"), "utf8")).figure : null;
  if (replay) {
    nonempty(figure, "paper_replay.json figure");
    assert(!isAbsolute(figure), "paper_replay.json figure must be relative");
    fileAt(within(replay, figure));
  }
  checkGuideAssets();
  const filename = resolve(docs, "trace-data.js");
  const sandbox = {window: {}};
  runInNewContext(readFileSync(filename, "utf8"), sandbox, {filename, timeout: 1000});
  const data = sandbox.window.UPSiLON_TRACE;
  const {referenceRevision, sources, routes} = data;
  assert.match(referenceRevision, /^[a-f0-9]{40}$/, "Pin a full Git commit for trace sources");
  assert.deepEqual(Object.keys(routes).sort(), [...routeIds].sort(), "Unexpected trace routes");
  for (const [key, source] of Object.entries(sources)) checkSource(source, key, referenceRevision);
  const replayFiles = new Set();
  for (const route of routeIds) {
    const {label, question, description, steps} = routes[route];
    nonempty(label, route + " label");
    nonempty(question, route + " question");
    checkExplanation(description, data, route + " description");
    const expected = ["w", "differential"].includes(route) ? spectrumStageIds : stageIds;
    assert.deepEqual(Array.from(steps, step => step.id), expected, route + " stages must follow calculation order");
    for (const step of steps) {
      const where = route + "/" + step.id;
      for (const key of ["title", "input", "output"]) nonempty(step[key], where + " " + key);
      assert(checkExplanation(step.detail, data, where), where + " needs a source or stage in its explanation");
      assert(Array.isArray(step.links) && step.links.length, where + " requires source links");
      for (const key of step.links) assert(Object.hasOwn(sources, key), where + " has an unknown source: " + key);
      assert(Array.isArray(step.files ?? []), where + " files must be an array");
      for (const template of step.files ?? []) {
        nonempty(template, where + " file");
        if (replay && route !== "beams") {
          const path = replayPath(template, figure);
          assert(!isAbsolute(path), where + " requires a relative replay path");
          replayFiles.add(within(replay, path));
        }
      }
    }
  }
  for (const path of replayFiles) fileAt(path);
  const stageCount = Object.values(routes).reduce((sum, route) => sum + route.steps.length, 0);
  console.log("PASS · " + routeIds.length + " results, " + stageCount + " stages, " + Object.keys(sources).length + " source links" +
    (replay ? ", " + replayFiles.size + " paper replay files" : ""));
}

try { main(process.argv.slice(2)); }
catch (error) { console.error("FAIL · " + error.message); process.exitCode = 1; }
