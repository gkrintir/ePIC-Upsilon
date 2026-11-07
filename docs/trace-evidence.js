/* Inspect small, existing run records in memory. Only the two published references are fetched. */
(() => {
  "use strict";
  const states = ["1S", "2S", "3S"];
  const runRoot = "results/generated/runs/";
  const manifest = "inputs/paper_inputs.csv";
  const checksums = "paper/reference/checksums.json";
  const publicReferences = new Map();
  const isSHA = value => typeof value === "string" && /^[0-9a-f]{64}$/i.test(value);

  // The status table contains quoted JSON with embedded commas and escaped quotes.
  function csv(source) {
    const input = source.replace(/^\uFEFF/, "");
    const field = /("(?:[^"]|"")*"|[^",\r\n]*)(,|\r\n|\n|\r|$)/gy;
    const rows = [];
    let row = [];
    while (true) {
      const match = field.exec(input);
      if (!match) throw new Error("Invalid CSV quoting");
      const [, raw, separator] = match;
      row.push(raw.startsWith('"') ? raw.slice(1, -1).replace(/""/g, '"') : raw);
      if (separator !== ",") {
        if (row.length !== 1 || row[0] !== "") rows.push(row);
        row = [];
      }
      if (!separator) break;
    }
    const header = rows.shift();
    if (!header || header.some(name => !name) || new Set(header).size !== header.length ||
        rows.some(values => values.length !== header.length)) throw new Error("Invalid CSV columns");
    return rows.map(values => Object.fromEntries(header.map((key, i) => [key, values[i]])));
  }

  async function localText(files, path) {
    const file = files.get(path);
    if (!file) throw new Error("Missing " + path);
    if (!file.size || file.size > 1024 * 1024) throw new Error("Empty or oversized record at " + path);
    return file.text();
  }
  async function reference(revision, path) {
    if (!/^[0-9a-f]{40}$/i.test(revision) || ![manifest, checksums].includes(path)) {
      throw new Error("Invalid public reference");
    }
    const url = "https://raw.githubusercontent.com/gkrintir/ePIC-Upsilon/" + revision + "/" + path;
    if (!publicReferences.has(url)) {
      const read = async () => {
        const controller = new AbortController();
        const timeout = setTimeout(() => controller.abort(), 6000);
        try {
          const response = await fetch(url, {credentials: "omit", referrerPolicy: "no-referrer", signal: controller.signal});
          if (!response.ok) throw new Error("Unavailable public reference");
          return await response.text();
        } finally { clearTimeout(timeout); }
      };
      publicReferences.set(url, read().catch(error => { publicReferences.delete(url); throw error; }));
    }
    return publicReferences.get(url);
  }

  async function tupleEvidence(files, revision) {
    const audits = await Promise.all(states.map(async state => {
      const path = runRoot + "paper-" + state.toLowerCase() + "/artifacts/normalization/9x275_input_audit.json";
      try {
        const audit = JSON.parse(await localText(files, path));
        if (!isSHA(audit?.sha256)) throw new Error();
        return {state, sha: audit.sha256.toLowerCase()};
      } catch {
        return {state, error: state + " · Missing or unreadable recorded source SHA-256 in " + path};
      }
    }));
    let rows;
    if (audits.some(audit => audit.sha)) {
      try { rows = csv(await reference(revision, manifest)); } catch { /* Offline or invalid CSV. */ }
    }
    return audits.map(({state, sha, error}) => {
      if (error) return error;
      const prefix = "Υ(" + state + "), 9×275 · Recorded source SHA-256 " + sha;
      const hits = rows?.filter(row => row.relative_path?.startsWith("Upsilon" + state + "/9x275/"));
      if (hits?.length !== 1 || !isSHA(hits[0].sha256)) {
        return prefix + ". Comparison unavailable because the pinned manifest is unavailable or invalid.";
      }
      return prefix + (sha === hits[0].sha256.toLowerCase() ? " matches" : " DIFFERS FROM") +
        " the pinned paper input manifest. The browser compared the audit entry without reading the ROOT file.";
    });
  }

  async function statusEvidence(route, files) {
    const path = runRoot + "paper-reopt/artifacts/reoptimization/validation_status.csv";
    if (!files.has(path)) return ["Missing " + path];
    try {
      const rows = csv(await localText(files, path));
      const messages = [];
      for (const observable of route === "w" ? ["w"] : ["pt2", "y"]) {
        for (const state of states) {
          const hits = rows.filter(row => row.state === state && row.observable === observable);
          const status = hits[0]?.overall_status ?? hits[0]?.status;
          const inference = hits[0]?.inference_status ?? hits[0]?.stage_status;
          const valid = value => typeof value === "string" && /^[A-Z][A-Z0-9_]{0,79}$/.test(value);
          const label = "Υ(" + state + "), " + ({pt2: "pT²", y: "rapidity", w: "W"})[observable];
          messages.push(hits.length === 1 && valid(status) && valid(inference)
            ? label + " · Recorded overall status " + status + ", inference " + inference + "."
            : label + " · Missing, duplicated or invalid validation status in " + path);
        }
      }
      return messages;
    } catch { return ["Cannot read the validation table at " + path + ". Check its CSV format and size."]; }
  }

  async function spectralDigest(files, record, revision) {
    const path = "spectra/comparison.json";
    if (!files.has(path)) return "Spectrum digest comparison unavailable. Missing " + path;
    let result;
    try { result = JSON.parse(await localText(files, path)); }
    catch { return "Cannot read the spectrum comparison record at " + path; }
    if (!isSHA(result?.sha256)) return "No valid recorded spectrum SHA-256 in " + path;
    const sha = result.sha256.toLowerCase();
    if (record?.spectra?.sha256 !== undefined &&
        (!isSHA(record.spectra.sha256) || record.spectra.sha256.toLowerCase() !== sha)) {
      return "The SHA-256 fields in paper_replay.json and spectra/comparison.json conflict.";
    }
    let expected;
    try { expected = JSON.parse(await reference(revision, checksums)).spectra; } catch { /* Offline or invalid JSON. */ }
    const prefix = "Recorded spectrum fingerprint " + sha;
    if (!isSHA(expected)) return prefix + ". Comparison unavailable because the pinned checksum record is unavailable or invalid.";
    return prefix + (sha === expected.toLowerCase() ? " matches" : " DIFFERS FROM") +
      " the pinned paper checksum. This comparison uses the digest recorded by the replay.";
  }

  window.inspectTraceEvidence = async ({route, stage, kind, files, record, referenceRevision}) => {
    if (kind !== "paper") return [];
    if (stage === "tuple" && ["mass", "efficiency", "w", "differential"].includes(route)) {
      return tupleEvidence(files, referenceRevision);
    }
    if (["w", "differential"].includes(route) && ["analysis", "plot", "publication"].includes(stage)) {
      const [statuses, digest] = await Promise.all([
        statusEvidence(route, files), spectralDigest(files, record, referenceRevision)
      ]);
      return [...statuses, digest];
    }
    return [];
  };
})();
