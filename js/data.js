"use strict";

/*
 * Data handling: default profiles (fetch), custom vehicles/payloads (session),
 * and file import/export (JSON, CSV, TSV – delimiter auto-detected).
 * No values are hardcoded here – the vehicle schema is the list below and all
 * defaults come from data/*.json.
 */
const Data = (() => {

  // ---------- Default profiles ----------
  async function loadDefaults() {
    const read = async (url) => {
      try {
        const r = await fetch(url, { cache: "no-store" });
        if (!r.ok) throw new Error("HTTP " + r.status);
        return await r.json();
      } catch (e) {
        console.warn("Could not load " + url, e);
        return null;
      }
    };
    const model = await read("data/model.json");
    const variables = await read("data/variables.json");
    const vehicles = await read("data/vehicles.json");
    const payloads = await read("data/payloads.json");
    const missing = [["model.json", model], ["variables.json", variables], ["vehicles.json", vehicles], ["payloads.json", payloads]].filter(([, v]) => !v).map(([n]) => n);
    return { ok: !missing.length, missing, model, variables, vehicles: vehicles || [], payloads: payloads || [] };
  }

  // ---------- CSV/TSV parsing (delimiter auto-detected) ----------
  function detectDelimiter(text) {
    const t = text.trim();
    if (t.indexOf("\t") >= 0) return "\t";
    const sem = (t.split(";").length - 1);
    const com = (t.split(",").length - 1);
    return sem > com ? ";" : ",";
  }

  function parseDelimited(text, requiredHeaders) {
    const delim = detectDelimiter(text);
    const lines = text.replace(/\r/g, "").split("\n").filter(l => l.trim() !== "");
    if (lines.length < 2) throw new Error("No data (needs a header row + at least one row).");
    const parseRow = (line) => {
      const out = [];
      let cur = "", inQ = false;
      for (let i = 0; i < line.length; i++) {
        const ch = line[i];
        if (ch === '"') { inQ = !inQ; continue; }
        if (ch === delim && !inQ) { out.push(cur.trim()); cur = ""; continue; }
        cur += ch;
      }
      out.push(cur.trim());
      return out;
    };
    const header = parseRow(lines[0]).map(h => h.trim());
    const lower = header.map(h => h.toLowerCase());
    (requiredHeaders || []).forEach(h => {
      if (!lower.includes(h.toLowerCase())) throw new Error("Missing column: " + h);
    });
    const rows = [];
    for (let i = 1; i < lines.length; i++) {
      const vals = parseRow(lines[i]);
      const obj = {};
      header.forEach((h, idx) => { obj[h] = vals[idx] != null ? vals[idx] : ""; });
      rows.push(obj);
    }
    return rows;
  }

  const num = (v) => {
    if (v == null || v === "") return null;
    const n = Number(String(v).replace(",", "."));
    return isNaN(n) ? null : n;
  };
  const bool = (v) => v === true || /^(true|1|yes|ja)$/i.test(String(v).trim());

  // ---------- Schemas ----------
  // Numeric vehicle fields (flat so that CSV/TSV round-trips cleanly).
  const vehicleNumeric = [
    "hullLength_m", "hullDiameter_m", "xcg_m", "volumeCoeff", "wettedCoeff", "pitchInertia_kgm2", "extraPitchDamping_Nms",
    "frontArea_m2", "frontAR", "frontX_m", "rearArea_m2", "rearAR", "rearX_m",
    "maxPower_kW", "jetEfficiency", "nozzleDiameter_m", "bsfc_kgpkWh", "emptyMass_kg", "fuelMass_kg"
  ];
  const vehicleHeaders = ["id", "name", "type", ...vehicleNumeric, "cruciform"];
  const vehicleRequired = ["hullLength_m", "hullDiameter_m", "frontArea_m2", "rearArea_m2"];
  const payloadHeaders = ["id", "name", "type", "weightKg"];

  function rowToVehicle(row) {
    const v = {
      id: String(row.id || "").trim() || "custom-" + Math.random().toString(36).slice(2, 8),
      name: String(row.name || "").trim() || "Unnamed vehicle",
      type: String(row.type || "").trim() || "custom",
      cruciform: row.cruciform == null || row.cruciform === "" ? true : bool(row.cruciform)
    };
    vehicleNumeric.forEach(k => { v[k] = num(row[k]); });
    if (v.extraPitchDamping_Nms == null) v.extraPitchDamping_Nms = 0;
    return v;
  }

  function rowToPayload(row) {
    return {
      id: String(row.id || "").trim() || "pl-" + Math.random().toString(36).slice(2, 8),
      name: String(row.name || "").trim() || "Unnamed payload",
      type: String(row.type || "").trim() || "other",
      weightKg: num(row.weightKg) || 0
    };
  }

  function parseVehiclesText(text, format) {
    const rows = format === "json" ? JSON.parse(text) : parseDelimited(text, vehicleRequired);
    return (Array.isArray(rows) ? rows : rows.vehicles || []).map(rowToVehicle);
  }

  function parsePayloadsText(text, format) {
    const rows = format === "json" ? JSON.parse(text) : parseDelimited(text, ["weightKg"]);
    return (Array.isArray(rows) ? rows : rows.payloads || []).map(rowToPayload);
  }

  // ---------- Serialisation ----------
  function toDelimited(rows, headers, delim) {
    const esc = (v) => {
      if (v == null) return "";
      const s = String(v);
      return (s.indexOf(delim) >= 0 || s.indexOf('"') >= 0 || s.indexOf("\n") >= 0) ? '"' + s.replace(/"/g, '""') + '"' : s;
    };
    const lines = [headers.join(delim)];
    rows.forEach(r => lines.push(headers.map(h => esc(r[h])).join(delim)));
    return lines.join("\n") + "\n";
  }

  function download(filename, content, mime) {
    const blob = new Blob([content], { type: mime || "application/octet-stream" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url; a.download = filename;
    document.body.appendChild(a); a.click();
    setTimeout(() => { URL.revokeObjectURL(url); a.remove(); }, 100);
  }

  const serializeVehicles = (list, format, delim) => format === "json" ? JSON.stringify(list, null, 2) : toDelimited(list, vehicleHeaders, delim);
  const serializePayloads = (list, format, delim) => format === "json" ? JSON.stringify(list, null, 2) : toDelimited(list, payloadHeaders, delim);

  return {
    loadDefaults, detectDelimiter, parseDelimited,
    parseVehiclesText, parsePayloadsText,
    serializeVehicles, serializePayloads,
    download, vehicleHeaders, payloadHeaders
  };
})();
