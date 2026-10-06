"use strict";

/*
 * Application logic – binds the UI, the physics engine (Model) and data (Data).
 * Nothing is hardcoded here: input fields, ranges and defaults are generated
 * from data/variables.json (fields + spans), data/model.json (defaults) and the
 * profiles in data/.
 */
(() => {

  const $ = (id) => document.getElementById(id);
  const fmt = (v, d = 1) => {
    if (v == null || Number.isNaN(v)) return "–";
    if (!Number.isFinite(v)) return "∞";
    return Number(v).toLocaleString("en-GB", { minimumFractionDigits: 0, maximumFractionDigits: d });
  };
  const DEG = 180 / Math.PI;
  const clamp = (v, lo, hi) => Math.max(lo, Math.min(hi, v));
  const show = (el) => el.classList.remove("hidden");
  const hide = (el) => el.classList.add("hidden");
  const numOr = (v, d) => { const n = parseFloat(v); return Number.isNaN(n) ? d : n; };
  const esc = (s) => String(s).replace(/[&<>"']/g, c => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]));

  // ---------- State ----------
  const state = {
    model: null,
    defaultVehicles: [], customVehicles: [],
    defaultPayloads: [], customPayloads: [],
    activeVehicleId: null,
    design: null,          // editable copy of the selected profile
    env: null,             // environment + operation inputs
    payloadId: null, payloadKg: 0,
    last: null,
    scenarios: [], history: [],
    map: null, mapConfig: null, mapCenter: null, mapTouched: false, mapLayer: null, marker: null,
    pending: false
  };

  const allVehicles = () => state.defaultVehicles.concat(state.customVehicles);
  const allPayloads = () => state.defaultPayloads.concat(state.customPayloads);
  const selectedProfile = () => allVehicles().find(v => v.id === state.activeVehicleId) || allVehicles()[0] || null;

  // ---------- Init ----------
  async function init() {
    bindStaticEvents();
    const res = await Data.loadDefaults();
    if (!res.ok) {
      show($("banner"));
      $("bannerText").textContent = "Could not load: " + res.missing.join(", ") + ". If you edited one of these files, check it is valid JSON (a missing comma or quote is the usual cause). If the page was opened directly from disk (file://), start the local server (start-server.bat) instead.";
      return;
    }
    state.model = res.model;
    state.model.sections = res.sections;
    state.vars = res.variables;
    try {
      state.model.presets = Model.presetsFromVariables(state.vars);
    } catch (e) {
      show($("banner")); $("bannerText").textContent = e.message; return;
    }
    const pv = state.vars.payload || { min: 0, max: 50, step: 0.1 };
    Object.assign($("payloadWeight"), { min: pv.min, max: pv.max, step: pv.step });
    state.defaultVehicles = res.vehicles;
    state.defaultPayloads = res.payloads;
    state.env = { ...res.model.environmentDefaults, turbulence: { ...res.model.turbulence } };
    $("vehicleColumns").textContent = Data.vehicleHeaders.join(", ");
    restoreHistory();
    buildVehicleSelect();
    buildPayloadSelect();
    loadProfile(selectedProfile());
    renderEnvBoxes();
    initMap();
    renderScenarios();
    schedule();
  }

  // ---------- Field generation from data/variables.json ----------
  // target: "design" | "env" | "turb"
  function fieldMax(def) {
    return typeof def.max === "string" ? state.design[def.max] : def.max;
  }
  function renderFields(group, containerId, target) {
    const box = $(containerId);
    box.innerHTML = "";
    (state.vars.groups[group] || []).forEach(def => {
      const obj = target === "design" ? state.design : target === "turb" ? state.env.turbulence : state.env;
      // Keep the computed value equal to what the slider shows when a span is narrowed.
      if (def.kind === "slider" && obj[def.key] != null) obj[def.key] = clamp(obj[def.key], def.min, fieldMax(def));
      const id = "f_" + def.key;
      const lab = document.createElement("label");
      lab.className = "field";
      lab.htmlFor = id;
      const title = document.createElement("span");
      title.textContent = def.label + (def.unit && !def.percent ? " [" + def.unit + "]" : "");
      lab.appendChild(title);
      const inp = document.createElement("input");
      inp.id = id; inp.className = "input";
      inp.dataset.key = def.key; inp.dataset.target = target;
      if (def.kind === "slider") {
        inp.type = "range"; inp.min = def.min; inp.max = fieldMax(def); inp.step = def.step;
        inp.value = obj[def.key];
        const out = document.createElement("span");
        out.className = "out"; out.id = id + "_out";
        lab.appendChild(inp); lab.appendChild(out);
      } else {
        inp.type = "number"; inp.step = def.step;
        inp.value = obj[def.key] == null ? "" : obj[def.key];
        if (def.nullable) inp.placeholder = "auto";
        lab.appendChild(inp);
      }
      inp.addEventListener("input", () => onField(def, inp, target));
      box.appendChild(lab);
      updateOut(def, target);
    });
  }

  function onField(def, inp, target) {
    const obj = target === "design" ? state.design : target === "turb" ? state.env.turbulence : state.env;
    if (inp.type === "number" && inp.value.trim() === "" && def.nullable) obj[def.key] = null;
    else obj[def.key] = numOr(inp.value, obj[def.key]);
    updateOut(def, target);
    if (def.key === "hullLength_m") syncPositionLimits();
    schedule();
  }

  function updateOut(def, target) {
    const out = $("f_" + def.key + "_out");
    if (!out) return;
    const obj = target === "design" ? state.design : target === "turb" ? state.env.turbulence : state.env;
    const v = obj[def.key];
    if (def.percent) out.textContent = fmt(v * 100, 0) + " %";
    else if (def.key === "windFromDeg") out.textContent = compass(v) + " (" + fmt(v, 0) + "°)";
    else out.textContent = fmt(v, def.digits) + (def.unit ? " " + def.unit : "");
    if (def.key === "cruiseSpeed_mps") out.textContent += " (" + fmt(v * 1.943844, 1) + " kn)";
  }

  // Position sliders (CG, fins) are bounded by the hull length.
  function syncPositionLimits() {
    ["hull", "front", "rear"].forEach(g => (state.vars.groups[g] || []).forEach(def => {
      if (typeof def.max !== "string") return;
      const el = $("f_" + def.key);
      if (!el) return;
      const mx = fieldMax(def);
      el.max = mx;
      if (state.design[def.key] > mx) { state.design[def.key] = mx; el.value = mx; }
      updateOut(def, "design");
    }));
  }

  function renderDesignBoxes() {
    renderFields("hull", "box_hull", "design");
    renderFields("front", "box_front", "design");
    renderSectionSelect("frontSection", "box_front");
    renderFields("rear", "box_rear", "design");
    renderSectionSelect("rearSection", "box_rear");
    renderFields("propulsion", "box_propulsion", "design");
    $("cruciform").checked = !!state.design.cruciform;
  }
  // 2D section (airfoil) picker – options come from data/sections.json.
  function renderSectionSelect(key, containerId) {
    const sc = state.model.sections;
    const list = (sc && sc.sections) || [];
    if (!list.length) return;
    const lab = document.createElement("label");
    lab.className = "field"; lab.htmlFor = "f_" + key;
    const title = document.createElement("span");
    title.textContent = "2D section (airfoil)";
    const sel = document.createElement("select");
    sel.id = "f_" + key; sel.className = "input";
    list.forEach(x => {
      const o = document.createElement("option");
      o.value = x.id; o.textContent = x.name + (x.confidence === "rough" ? " – rough data" : "");
      sel.appendChild(o);
    });
    if (!list.some(x => x.id === state.design[key])) state.design[key] = sc.default || list[0].id;
    sel.value = state.design[key];
    const info = document.createElement("span");
    info.className = "hint"; info.id = "f_" + key + "_info";
    const describe = () => { const x = list.find(y => y.id === sel.value); info.textContent = x ? (x.note ? x.note + " " : "") + "Source: " + x.source : ""; };
    describe();
    sel.addEventListener("change", () => { state.design[key] = sel.value; describe(); schedule(); });
    lab.appendChild(title); lab.appendChild(sel); lab.appendChild(info);
    $(containerId).appendChild(lab);
  }

  function renderEnvBoxes() {
    renderFields("environment", "box_environment", "env");
    renderFields("turbulence", "box_turbulence", "turb");
    $("useVehicleResponse").checked = state.env.useVehicleResponse !== false;
    $("includeQPenalty").checked = !!state.env.includeQPenalty;
    $("turbPlaceholder").classList.toggle("hidden", !state.model.turbulence.placeholder);
  }

  const COMPASS = ["N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE", "S", "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW"];
  const compass = (deg) => COMPASS[Math.round((((deg % 360) + 360) % 360) / 22.5) % 16];

  // ---------- Profiles & payloads ----------
  function buildVehicleSelect() {
    const sel = $("vehicleSelect");
    sel.innerHTML = "";
    allVehicles().forEach(v => {
      const o = document.createElement("option");
      o.value = v.id; o.textContent = v.name + " (" + (v.type || "custom") + ")";
      sel.appendChild(o);
    });
    if (!allVehicles().some(v => v.id === state.activeVehicleId)) state.activeVehicleId = allVehicles()[0] ? allVehicles()[0].id : null;
    sel.value = state.activeVehicleId;
  }

  function loadProfile(p) {
    if (!p) return;
    state.activeVehicleId = p.id;
    state.design = JSON.parse(JSON.stringify(p));
    renderDesignBoxes();
    renderCustomLists();
    fillCalibTheory();
  }

  function buildPayloadSelect() {
    const sel = $("payloadSelect");
    const prev = state.payloadId;
    sel.innerHTML = "";
    allPayloads().forEach(l => {
      const o = document.createElement("option");
      o.value = l.id; o.textContent = l.name + " (" + fmt(l.weightKg, 1) + " kg)";
      sel.appendChild(o);
    });
    const custom = document.createElement("option");
    custom.value = "__custom"; custom.textContent = "Custom mass (slider)";
    sel.appendChild(custom);
    if (allPayloads().some(l => l.id === prev)) sel.value = prev;
    else if (allPayloads().length) { sel.value = allPayloads()[0].id; state.payloadKg = allPayloads()[0].weightKg; }
    else sel.value = "__custom";
    state.payloadId = sel.value;
    $("payloadWeight").value = state.payloadKg;
    $("payloadWeightLabel").textContent = fmt(state.payloadKg, 1) + " kg";
  }

  const currentPayloadName = () => { const l = allPayloads().find(x => x.id === state.payloadId); return l ? l.name : "Custom"; };

  // ---------- Calibration ----------
  function calibInput() {
    const use = (id) => $(id + "_use").checked ? numOr($(id).value, null) : null;
    const c = {
      hullCDwet: use("cal_hull"),
      clAlphaFront: use("cal_claf"),
      clAlphaRear: use("cal_clar"),
      etaJet: use("cal_eta"),
      extraDragArea_m2: use("cal_extra")
    };
    $("cal_runInfo").textContent = "";
    if ($("cal_run_use").checked && c.hullCDwet == null) {
      const P = numOr($("cal_runP").value, null), V = numOr($("cal_runV").value, null);
      if (P > 0 && V > 0) {
        const derived = deriveHullCD(P * 1000, V, c);
        if (derived != null && derived > 0) {
          c.hullCDwet = derived;
          c.hullSrc = "derived from run";
          $("cal_runInfo").textContent = "Derived hull C_D,wet = " + fmt(derived, 5) + " (calm water, run at " + fmt(V, 1) + " m/s, " + fmt(P, 1) + " kW).";
        } else {
          $("cal_runInfo").textContent = "The run implies less drag than the fins alone – check the measurement.";
        }
      }
    }
    return c;
  }

  // Calm-water run: thrust from measured engine power, minus non-hull drag → hull C_D,wet.
  function deriveHullCD(P_W, V, c) {
    const calmEnv = { ...state.env, cruiseSpeed_mps: V, payload_kg: state.payloadKg, turbulence: { ...state.env.turbulence, sigma_u: 0, sigma_v: 0, sigma_w: 0 } };
    const r = Model.analyze(state.design, calmEnv, state.model, { ...c, hullCDwet: 0 });
    if (!r) return null;
    const eta = c.etaJet != null ? c.etaJet : state.design.jetEfficiency;
    const T = Model._internal.thrustForJetPower(P_W * eta, V, r.ctx.rho, r.ctx.An).T;
    const other = r.chosen.drag.total;
    return (T - other) / (r.chosen.qbar * r.ctx.hull.wetted);
  }

  function fillCalibTheory() {
    if (!state.design) return;
    // Theory = Helmbold with the selected section's cℓα (ignoring any active override).
    const r = state.last, I = Model._internal;
    const theory = (n) => r ? I.helmbold(state.design[n + "AR"], r.ctx[n].section.a0) : I.helmbold(state.design[n + "AR"]);
    $("cc_claf").textContent = fmt(theory("front"), 2);
    $("cc_clar").textContent = fmt(theory("rear"), 2);
    $("cc_eta").textContent = fmt(state.design.jetEfficiency, 2);
    $("cc_extra").textContent = fmt(state.model.hydro.extraDragArea_m2 || 0, 4);
  }

  // ---------- Compute (coalesced; a timer rather than rAF so it also runs in background tabs) ----------
  function schedule() {
    if (state.pending) return;
    state.pending = true;
    setTimeout(() => { state.pending = false; compute(); }, 16);
  }

  function compute() {
    if (!state.model || !state.design) return null;
    const env = { ...state.env, payload_kg: state.payloadKg };
    const r = Model.analyze(state.design, env, state.model, calibInput());
    state.last = r;
    fillCalibTheory();
    renderResults(r);
    return r;
  }

  // ---------- Results ----------
  function renderResults(r) {
    if (!r) {
      ["rangeBig", "endureBig", "smBig"].forEach(id => $(id).textContent = "–");
      return;
    }
    const c = r.chosen, ctx = r.ctx;
    $("rangeBig").textContent = fmt(c.rangeKm, 1);
    $("rangeSub").textContent = "at " + fmt(c.V, 2) + " m/s" + (c.feasible ? "" : " – NOT reachable with available power");
    const hrs = c.enduranceH, h = Math.floor(hrs), m = Math.round((hrs - h) * 60);
    $("endureBig").textContent = Number.isFinite(hrs) ? h + " h " + String(m).padStart(2, "0") + " min" : "–";
    $("smBig").textContent = fmt(r.staticMargin_pct, 1);
    $("smSub").textContent = r.staticMargin_m > 0 ? "of hull length · stable" : "of hull length · UNSTABLE";

    chip("cf_stab", r.staticMargin_m > 0 && c.vert.stable && (!ctx.cruciform || c.lat.stable) ? "✓ stable" : "✕ unstable", r.staticMargin_m > 0 && c.vert.stable && ctx.cruciform && c.lat.stable ? "good" : "bad");
    chip("cf_feas", (c.feasible ? "✓ " : "✕ ") + fmt(c.throttle * 100, 0) + " % power", c.feasible ? "good" : "bad");
    chip("cf_cla", "CLα " + fmt(ctx.front.clAlpha, 2) + " / " + fmt(ctx.rear.clAlpha, 2), ctx.front.clAlphaSrc === "empirical" || ctx.rear.clAlphaSrc === "empirical" ? "emp" : "teo");
    chip("cf_hull", "C_D,wet " + fmt(c.CDwet, 5), ctx.hullCDwetOverride != null ? "emp" : "teo");
    chip("cf_eta", "η_jet " + fmt(ctx.etaJet, 2), $("cal_eta_use").checked ? "emp" : "teo");
    chip("cf_mode", c.useDynV ? "gust: with vehicle response" : "gust: fixed attitude", "teo");

    // Trajectory stability
    const modeTxt = (sys) => !sys.stable ? "unstable" : (sys.zeta >= 1 ? "overdamped" : "ω_n " + fmt(sys.wn, 2) + " rad/s, ζ " + fmt(sys.zeta, 2));
    $("tsNP").textContent = fmt(r.xnp, 3) + " m from nose";
    $("tsSM").textContent = fmt(r.staticMargin_m, 3) + " m (" + fmt(r.staticMargin_pct, 1) + " % L)";
    $("tsPitch").textContent = modeTxt(c.vert) + (c.vert.stable ? " · eig " + c.vert.eig.map(e => fmt(e.re, 1) + (Math.abs(e.im) > 1e-9 ? (e.im > 0 ? "+" : "−") + fmt(Math.abs(e.im), 1) + "i" : "")).join(", ") : "");
    $("tsYaw").textContent = ctx.cruciform ? modeTxt(c.lat) : "no vertical surfaces";
    $("tsAlpha").textContent = fmt(c.sv.sigFinAlpha[0] * DEG, 2) + "° / " + fmt(c.sv.sigFinAlpha[1] * DEG, 2) + "°";
    $("tsN").textContent = "1 ± " + fmt(c.sv.sigL / ctx.W, 2) + " g";
    $("tsQ").textContent = fmt(c.sv.sigQ * DEG, 2) + " °/s";
    $("tsGamma").textContent = fmt(r.pathAngle.vert_rad * DEG, 2) + "° / " + fmt(r.pathAngle.lat_rad * DEG, 2) + "°";
    $("tsD1").textContent = fmt(r.per1km.depth_m, 1) + " m / " + fmt(r.per1km.lateral_m, 1) + " m";
    $("tsDR").textContent = r.atRange ? fmt(r.atRange.depth_m, 0) + " m / " + fmt(r.atRange.lateral_m, 0) + " m (at " + fmt(r.atRange.km, 0) + " km)" : "–";

    // Propulsion & drag
    const D = c.drag;
    $("pdD").textContent = fmt(D.total, 1) + " N (calm water " + fmt(D.calm, 1) + " N)";
    $("pdP").textContent = fmt(c.P_engine / 1000, 2) + " kW of " + fmt(state.design.maxPower_kW, 1) + " kW (" + fmt(c.throttle * 100, 0) + " %)";
    $("pdJet").textContent = fmt(c.jet.Vj, 2) + " m/s · η_F " + fmt(c.jet.etaF, 3);
    $("pdEta").textContent = fmt(c.etaTot, 3) + " (η_jet × η_F)";
    $("pdFuel").textContent = fmt(c.fuelFlow_kgph, 2) + " kg/h";
    $("pdC").textContent = fmt(c.cEq_perHour, 3) + " 1/h";
    $("pdLD").textContent = fmt(c.LD, 3);
    $("pdVmax").textContent = r.vMax != null ? fmt(r.vMax, 2) + " m/s" : "none in scan";
    $("pdBest").textContent = r.best ? fmt(r.best.V, 2) + " m/s → " + fmt(r.best.rangeKm, 1) + " km" : "–";
    const parts = [["Hull", D.hull], ["Fin profile", D.finProfile], ["Trim induced", D.trim], ["Gust vert.", D.gustV], ["Gust lat.", D.gustL], ["Extra", D.extra]];
    $("dragBreakdown").innerHTML = parts.map(([n, v]) => "<div class='stat'><span>" + n + "</span><b>" + fmt(v, 1) + " N</b><em>" + fmt(v / D.total * 100, 0) + " %</em></div>").join("");

    // Warnings
    const wb = $("warningsBox");
    if (r.warnings.length) {
      show(wb); wb.innerHTML = "";
      r.warnings.forEach(w => { const d = document.createElement("div"); d.className = "warn-item"; d.textContent = "⚠ " + w; wb.appendChild(d); });
    } else hide(wb);

    // Steps
    const sbox = $("stepsBox"); sbox.textContent = "";
    r.steps.forEach(s => {
      const row = document.createElement("div"); row.className = "step-row";
      const t = document.createElement("div"); t.className = "step-t"; t.textContent = s.t;
      const d = document.createElement("div"); d.className = "step-d"; d.textContent = s.d;
      row.appendChild(t); row.appendChild(d); sbox.appendChild(row);
    });

    drawSideView(r);
    drawCharts(r);
    drawMap(Model.envelope(r, state.env));
  }

  function chip(id, text, kind) {
    const el = $(id);
    el.textContent = text;
    el.className = "chip " + ({ good: "chip-good", bad: "chip-bad", emp: "chip-emp", teo: "chip-teo" }[kind] || "");
  }

  // ---------- Side view (inline SVG) ----------
  function drawSideView(r) {
    const ctx = r.ctx, d = state.design;
    const W = 640, Hmax = 220, pad = 28;
    const bMax = Math.max(ctx.front.span, ctx.rear.span);
    const s = Math.min((W - 2 * pad) / ctx.hull.L, (Hmax - 2 * pad - 24) / (ctx.hull.D + bMax));
    const H = Math.max(120, (ctx.hull.D + bMax) * s + 2 * pad + 24);
    const x0 = (W - ctx.hull.L * s) / 2, yc = pad + (bMax / 2 + ctx.hull.D / 2) * s;
    const X = (x) => x0 + x * s;
    const hr = Math.max(3, ctx.hull.D * s / 2);
    let svg = `<svg viewBox="0 0 ${W} ${H}" role="img" aria-label="Side view of the AUV design"><title>Side view</title>`;
    // Hull (ogive-ish nose, tapered tail)
    const xn = X(0), xt = X(ctx.hull.L), nose = Math.min(hr * 2.2, (xt - xn) * 0.2), tail = Math.min(hr * 3, (xt - xn) * 0.25);
    svg += `<path class="sv-hull" d="M ${xn} ${yc} Q ${xn} ${yc - hr} ${xn + nose} ${yc - hr} L ${xt - tail} ${yc - hr} L ${xt} ${yc - hr * 0.45} L ${xt} ${yc + hr * 0.45} L ${xt - tail} ${yc + hr} L ${xn + nose} ${yc + hr} Q ${xn} ${yc + hr} ${xn} ${yc} Z"/>`;
    // Nozzle
    const nz = Math.max(2, d.nozzleDiameter_m * s / 2);
    svg += `<rect class="sv-nozzle" x="${xt}" y="${yc - nz}" width="${Math.max(4, 0.04 * s)}" height="${2 * nz}"/>`;
    // Fins: leading edge at x_cl − c/4, semispan b/2 above and below the hull
    const fin = (p, cls, label) => {
      const le = X(p.xcl - p.chord / 4), te = X(p.xcl + 0.75 * p.chord), semi = p.span / 2 * s, sweep = 0.3 * p.chord * s;
      let out = "";
      [-1, 1].forEach(dir => {
        const yr = yc + dir * hr, yt = yr + dir * semi;
        out += `<path class="${cls}" d="M ${le} ${yr} L ${le + sweep} ${yt} L ${te} ${yt} L ${te} ${yr} Z"/>`;
      });
      out += `<text class="sv-label" x="${X(p.xcl)}" y="${yc - hr - semi - 6}" text-anchor="middle">${label} ${fmt(p.xcl, 2)} m</text>`;
      return out;
    };
    svg += fin(ctx.front, "sv-fin sv-fin-front", "front");
    svg += fin(ctx.rear, "sv-fin sv-fin-rear", "rear");
    // CG
    const cgx = X(ctx.xcg);
    svg += `<circle class="sv-cg" cx="${cgx}" cy="${yc}" r="6"/><path class="sv-cg-q" d="M ${cgx} ${yc} L ${cgx + 6} ${yc} A 6 6 0 0 1 ${cgx} ${yc + 6} Z M ${cgx} ${yc} L ${cgx - 6} ${yc} A 6 6 0 0 1 ${cgx} ${yc - 6} Z"/>`;
    // Neutral point
    if (Number.isFinite(r.xnp)) {
      const npx = clamp(X(r.xnp), 4, W - 4), yb = yc + hr + 4;
      svg += `<path class="${r.staticMargin_m > 0 ? "sv-np" : "sv-np bad"}" d="M ${npx} ${yb} L ${npx - 6} ${yb + 10} L ${npx + 6} ${yb + 10} Z"/>`;
      svg += `<text class="sv-label" x="${npx}" y="${yb + 22}" text-anchor="middle">NP ${fmt(r.xnp, 2)} m</text>`;
    }
    svg += `<text class="sv-label" x="${cgx}" y="${yc - hr - 4 < pad ? yc + hr + 34 : yc + hr + 34}" text-anchor="middle">CG ${fmt(ctx.xcg, 2)} m</text>`;
    // Ruler
    const yR = H - 8, step = ctx.hull.L > 4 ? 1 : ctx.hull.L > 1.5 ? 0.5 : 0.25;
    svg += `<line class="axis" x1="${X(0)}" y1="${yR - 6}" x2="${X(ctx.hull.L)}" y2="${yR - 6}"/>`;
    for (let x = 0; x <= ctx.hull.L + 1e-9; x += step) svg += `<line class="axis" x1="${X(x)}" y1="${yR - 9}" x2="${X(x)}" y2="${yR - 3}"/><text class="tick" x="${X(x)}" y="${yR + 4}" text-anchor="middle">${fmt(x, 2)}</text>`;
    svg += "</svg>";
    $("sideView").innerHTML = svg;
  }

  // ---------- Charts (inline SVG; crosshair + tooltip) ----------
  // series: [{ name, cls, data:[{x,y}] }] – y may be non-finite (gap).
  function makeChart(container, opts) {
    const W = 600, H = 250, padL = 54, padR = 14, padT = 26, padB = 34;
    const vw = W - padL - padR, vh = H - padT - padB;
    const tx = (v) => opts.logX ? Math.log10(v) : v, ty = (v) => opts.logY ? Math.log10(v) : v;
    const pts = opts.series.flatMap(s => s.data).filter(p => Number.isFinite(p.y) && Number.isFinite(p.x) && (!opts.logY || p.y > 0) && (!opts.logX || p.x > 0));
    if (!pts.length) { container.innerHTML = "<p class='hint'>" + esc(opts.empty || "No data to plot.") + "</p>"; return; }
    let xMin = Math.min(...pts.map(p => tx(p.x))), xMax = Math.max(...pts.map(p => tx(p.x)));
    let yMin = Math.min(...pts.map(p => ty(p.y))), yMax = Math.max(...pts.map(p => ty(p.y)));
    (opts.hlines || []).forEach(h => { yMin = Math.min(yMin, h); yMax = Math.max(yMax, h); });
    if (!opts.logY && opts.zeroBase !== false) yMin = Math.min(0, yMin);
    if (yMax - yMin < 1e-12) { yMax = yMin + 1; }
    if (xMax - xMin < 1e-12) { xMax = xMin + 1; }
    const X = (x) => padL + (tx(x) - xMin) / (xMax - xMin) * vw;
    const Y = (y) => padT + (yMax - ty(y)) / (yMax - yMin) * vh;
    const tickTxt = (v, log) => log ? "10^" + Math.round(v) : fmt(v, Math.abs(v) < 10 ? 2 : 0);
    if (opts.logY) { yMin = Math.floor(yMin); yMax = Math.ceil(yMax); }

    // Round tick steps (1/2/5·10^n); integer decades on log axes.
    const ticks = (lo, hi, n, log) => {
      if (log) { const out = []; for (let v = Math.ceil(lo); v <= Math.floor(hi); v++) out.push(v); return out.length > 1 ? out : [lo, hi]; }
      const raw = (hi - lo) / n, mag = Math.pow(10, Math.floor(Math.log10(raw)));
      const step = [1, 2, 5, 10].map(f => f * mag).find(s => s >= raw);
      const out = []; for (let v = Math.ceil(lo / step) * step; v <= hi + 1e-9 * step; v += step) out.push(Math.abs(v) < 1e-12 ? 0 : v);
      return out;
    };
    const Yt = (t) => padT + (yMax - t) / (yMax - yMin) * vh, Xt = (t) => padL + (t - xMin) / (xMax - xMin) * vw;

    let svg = `<svg viewBox="0 0 ${W} ${H}" role="img" aria-label="${esc(opts.title)}"><title>${esc(opts.title)}</title>`;
    ticks(yMin, yMax, 4, opts.logY).forEach(t => {
      const yy = Yt(t);
      svg += `<line class="grid" x1="${padL}" y1="${yy}" x2="${padL + vw}" y2="${yy}"/><text class="tick" x="${padL - 6}" y="${yy + 3}" text-anchor="end">${tickTxt(t, opts.logY)}</text>`;
    });
    ticks(xMin, xMax, 5, opts.logX).forEach(t => {
      svg += `<text class="tick" x="${Xt(t)}" y="${padT + vh + 14}" text-anchor="middle">${tickTxt(t, opts.logX)}</text>`;
    });
    svg += `<line class="axis" x1="${padL}" y1="${padT + vh}" x2="${padL + vw}" y2="${padT + vh}"/>`;
    (opts.hlines || []).forEach(h => { svg += `<line class="chart-ref" x1="${padL}" y1="${Y(opts.logY ? Math.pow(10, h) : h)}" x2="${padL + vw}" y2="${Y(opts.logY ? Math.pow(10, h) : h)}"/>`; });
    (opts.markers || []).forEach(m => {
      if (!Number.isFinite(m.x)) return;
      const xx = X(m.x);
      if (xx < padL || xx > padL + vw) return;
      svg += `<line class="chart-hl" x1="${xx}" y1="${padT}" x2="${xx}" y2="${padT + vh}"/><text class="chart-label" x="${Math.min(xx + 4, padL + vw - 2)}" y="${padT + 10 + (m.row || 0) * 12}" text-anchor="${xx > padL + vw * 0.7 ? "end" : "start"}">${esc(m.label)}</text>`;
    });
    opts.series.forEach(s => {
      let d = "", pen = false;
      s.data.forEach(p => {
        const ok = Number.isFinite(p.y) && (!opts.logY || p.y > 0) && (!opts.logX || p.x > 0);
        if (!ok) { pen = false; return; }
        d += (pen ? " L" : " M") + X(p.x).toFixed(1) + " " + Y(p.y).toFixed(1); pen = true;
      });
      svg += `<path d="${d}" class="chart-line ${s.cls}"/>`;
    });
    // Legend (≥ 2 series) – line keys, text in text colour
    if (opts.series.length > 1) {
      let lx = padL + 4;
      opts.series.forEach(s => {
        svg += `<line class="chart-line ${s.cls}" x1="${lx}" y1="${12}" x2="${lx + 16}" y2="${12}"/><text class="legend" x="${lx + 20}" y="${16}">${esc(s.name)}</text>`;
        lx += 30 + s.name.length * 6.2;
      });
    } else {
      svg += `<text class="legend" x="${padL + 4}" y="${16}">${esc(opts.series[0].name)}</text>`;
    }
    svg += `<text class="chart-axislabel" x="${padL + vw / 2}" y="${H - 4}" text-anchor="middle">${esc(opts.xLabel)}</text>`;
    svg += `<text class="chart-axislabel" x="12" y="${padT + vh / 2}" text-anchor="middle" transform="rotate(-90 12 ${padT + vh / 2})">${esc(opts.yLabel)}</text>`;
    svg += `<line class="crosshair hidden" x1="0" y1="${padT}" x2="0" y2="${padT + vh}"/>`;
    svg += `<rect class="hit" x="${padL}" y="${padT}" width="${vw}" height="${vh}" fill="transparent"/>`;
    svg += "</svg>";
    container.innerHTML = svg;

    // Hover: crosshair snaps to the nearest x of the first series; tooltip lists every series.
    const el = container.querySelector("svg"), hit = el.querySelector(".hit"), ch = el.querySelector(".crosshair"), tip = $("chartTip");
    const xs = opts.series[0].data.map(p => p.x);
    const move = (ev) => {
      const box = el.getBoundingClientRect();
      const px = (ev.clientX - box.left) / box.width * W;
      let best = 0, bd = Infinity;
      xs.forEach((x, i) => { const dd = Math.abs(X(x) - px); if (dd < bd) { bd = dd; best = i; } });
      const xv = xs[best];
      ch.setAttribute("x1", X(xv)); ch.setAttribute("x2", X(xv)); ch.classList.remove("hidden");
      tip.textContent = "";
      const head = document.createElement("div"); head.className = "tip-head"; head.textContent = opts.xLabel + ": " + fmt(xv, opts.logX ? 4 : 2); tip.appendChild(head);
      opts.series.forEach(s => {
        const p = s.data.reduce((a, b) => Math.abs(b.x - xv) < Math.abs(a.x - xv) ? b : a, s.data[0]);
        const row = document.createElement("div"); row.className = "tip-row";
        const key = document.createElement("i"); key.className = "tip-key " + s.cls;
        const val = document.createElement("b"); val.textContent = opts.fmtY ? opts.fmtY(p.y) : fmt(p.y, 2);
        const nm = document.createElement("span"); nm.textContent = " " + s.name;
        row.appendChild(key); row.appendChild(val); row.appendChild(nm); tip.appendChild(row);
      });
      show(tip);
      tip.style.left = Math.min(window.innerWidth - tip.offsetWidth - 8, ev.clientX + 14) + "px";
      tip.style.top = (ev.clientY + 14) + "px";
    };
    hit.addEventListener("pointermove", move);
    hit.addEventListener("pointerleave", () => { ch.classList.add("hidden"); hide(tip); });
  }

  function drawCharts(r) {
    const c = r.chosen;
    makeChart($("chartSpeed"), {
      title: "One-way range vs cruise speed",
      series: [{ name: "One-way range (km)", cls: "s1", data: r.scan.map(p => ({ x: p.V, y: p.feasible ? p.rangeKm : NaN })) }],
      xLabel: "Cruise speed (m/s)", yLabel: "Range (km)",
      markers: [
        { x: c.V, label: "chosen " + fmt(c.V, 2) + " m/s" },
        r.best ? { x: r.best.V, label: "best " + fmt(r.best.V, 2) + " m/s", row: 1 } : null,
        r.vMax != null ? { x: r.vMax, label: "max " + fmt(r.vMax, 2) + " m/s", row: 2 } : null
      ].filter(Boolean),
      fmtY: (y) => fmt(y, 1) + " km",
      empty: "No speed in the scan is reachable with the available power."
    });
    makeChart($("chartDisp"), {
      title: "Track wander vs distance (no steering)",
      series: [
        { name: "Depth", cls: "s1", data: r.disp.map(p => ({ x: p.km, y: p.depth_m })) },
        { name: "Lateral", cls: "s2", data: r.disp.map(p => ({ x: p.km, y: p.lateral_m })) }
      ],
      xLabel: "Distance (km)", yLabel: "RMS wander (m)",
      fmtY: (y) => fmt(y, 1) + " m",
      empty: "Unbounded – the vehicle is unstable without control."
    });
    makeChart($("chartSM"), {
      title: "Static margin vs rear-pair position",
      series: [{ name: "Static margin (% of L)", cls: "s1", data: r.smSweep.map(p => ({ x: p.x, y: p.sm })) }],
      xLabel: "Rear pair x (m from nose)", yLabel: "Static margin (% L)",
      markers: [{ x: state.design.rearX_m, label: "current " + fmt(state.design.rearX_m, 2) + " m" }],
      hlines: [0], zeroBase: false,
      fmtY: (y) => fmt(y, 1) + " %"
    });
    makeChart($("chartSpec"), {
      title: "Incidence spectrum: raw gust vs felt at rear pair",
      series: [
        { name: "Raw (w/V)²", cls: "s1", data: r.spectrum.map(p => ({ x: p.Om, y: p.raw })) },
        { name: "Felt at rear pair", cls: "s2", data: r.spectrum.map(p => ({ x: p.Om, y: p.felt })) }
      ],
      xLabel: "Spatial frequency Ω (rad/m)", yLabel: "PSD (rad²·m)", logX: true, logY: true,
      fmtY: (y) => y.toExponential(2)
    });
  }

  // ---------- Map (Leaflet) ----------
  function initMap() {
    if (typeof L === "undefined") { $("mapDiv").textContent = "The Leaflet map library could not be loaded (needs internet)."; return; }
    fetch("data/map.json", { cache: "no-store" }).then(r => r.json()).then(cfg => {
      state.mapConfig = cfg;
      state.mapCenter = { lat: cfg.center.lat, lng: cfg.center.lng };
      const map = L.map("mapDiv", { center: [cfg.center.lat, cfg.center.lng], zoom: cfg.center.zoom });
      state.map = map;
      map.on("dragstart zoomstart", () => { state.mapTouched = true; });
      const base = {};
      Object.values(cfg.layers).forEach((ly, i) => {
        const t = L.tileLayer(ly.url, { attribution: ly.attribution, maxZoom: 19 });
        if (i === 0) t.addTo(map);
        base[ly.name] = t;
      });
      L.control.layers(base, null, { collapsed: false }).addTo(map);
      state.marker = L.marker([cfg.center.lat, cfg.center.lng], { draggable: true, zIndexOffset: 1000 })
        .addTo(map).bindTooltip("Start: " + cfg.center.name, { permanent: true, direction: "top" });
      state.marker.on("dragstart", () => { state.mapTouched = true; });
      state.marker.on("dragend", (e) => { const ll = e.target.getLatLng(); state.mapCenter = { lat: ll.lat, lng: ll.lng }; schedule(); });
      schedule();
    }).catch(() => { $("mapDiv").textContent = "data/map.json could not be loaded – run via the local server."; });
  }

  function drawMap(env) {
    if (!state.map || !state.mapCenter) return;
    const ctr = state.mapCenter, mPerDegLat = 111320;
    const point = (deg, km) => {
      const b = deg * Math.PI / 180, d = km * 1000;
      const lat = ctr.lat + d * Math.cos(b) / mPerDegLat;
      return [lat, ctr.lng + d * Math.sin(b) / (mPerDegLat * Math.cos(lat * Math.PI / 180))];
    };
    const poly = env.oneWayKm.map((km, i) => point(env.headings[i], km));
    if (state.mapLayer) state.map.removeLayer(state.mapLayer);
    state.mapLayer = null;
    if (poly.length >= 3 && env.max > 0) {
      state.mapLayer = L.polygon(poly.concat([poly[0]]), { color: "#d95926", weight: 2, fillColor: "#d95926", fillOpacity: 0.08 })
        .addTo(state.map).bindTooltip("One-way reach: " + fmt(env.min, 1) + "–" + fmt(env.max, 1) + " km", { sticky: true });
      if (!state.mapTouched) state.map.fitBounds(state.mapLayer.getBounds().pad(0.15));
    }
  }

  // ---------- Scenarios / history ----------
  function snapshot(r) {
    return {
      ts: new Date().toISOString(),
      name: state.design.name, payloadKg: state.payloadKg, payloadName: currentPayloadName(),
      V: r.chosen.V, rho: state.env.density_kgpm3,
      rangeKm: r.chosen.rangeKm, enduranceH: r.chosen.enduranceH,
      sm: r.staticMargin_pct, wander1km: r.per1km.depth_m, wanderLat1km: r.per1km.lateral_m,
      L: state.design.hullLength_m, D: state.design.hullDiameter_m,
      Sf: state.design.frontArea_m2, Sr: state.design.rearArea_m2, P: state.design.maxPower_kW
    };
  }
  function addScenario() {
    const r = compute(); if (!r) return;
    state.scenarios.push(snapshot(r));
    renderScenarios();
  }
  function renderScenarios() {
    const box = $("scenarioList"); box.innerHTML = "";
    state.scenarios.forEach(s => {
      const card = document.createElement("div"); card.className = "scenario-card";
      card.innerHTML =
        "<div class='sc-title'>" + esc(s.name) + "</div>" +
        "<div class='sc-sub'>L " + fmt(s.L, 2) + " m · D " + fmt(s.D, 3) + " m · S " + fmt(s.Sf, 3) + "/" + fmt(s.Sr, 3) + " m² · " + fmt(s.P, 1) + " kW · " + fmt(s.V, 2) + " m/s · ρ " + fmt(s.rho, 1) + " · " + esc(s.payloadName) + " " + fmt(s.payloadKg, 1) + " kg</div>" +
        "<div class='sc-nums'><div>Range<br><b>" + fmt(s.rangeKm, 1) + " km</b></div><div>Static margin<br><b>" + fmt(s.sm, 1) + " %</b></div><div>Wander/1 km<br><b>" + fmt(s.wander1km, 1) + " / " + fmt(s.wanderLat1km, 1) + " m</b></div></div>";
      box.appendChild(card);
    });
    $("scenarioCount").textContent = state.scenarios.length ? state.scenarios.length + " scenario(s)" : "None saved";
  }
  function addHistoryEntry(snap) {
    state.history.unshift(snap);
    if (state.history.length > 20) state.history.pop();
    try { sessionStorage.setItem("jetify_history", JSON.stringify(state.history)); } catch (e) { /* storage unavailable */ }
    renderHistory();
  }
  function restoreHistory() {
    try { const raw = sessionStorage.getItem("jetify_history"); if (raw) state.history = JSON.parse(raw) || []; } catch (e) { state.history = []; }
    renderHistory();
  }
  function renderHistory() {
    const box = $("historyList"); box.innerHTML = "";
    state.history.forEach(h => {
      const row = document.createElement("div"); row.className = "hist-row";
      row.innerHTML = "<span>" + esc(h.name) + " · " + fmt(h.V, 2) + " m/s · SM " + fmt(h.sm, 1) + " %</span><span>" + new Date(h.ts).toLocaleTimeString("en-GB") + "</span><b>" + fmt(h.rangeKm, 1) + " km</b>";
      box.appendChild(row);
    });
    $("historyCount").textContent = state.history.length ? state.history.length + " computation(s) this session" : "No history yet";
  }

  function summaryText(r) {
    const d = state.design, c = r.chosen, e = state.env, t = e.turbulence;
    return [
      "JETIFY – SUMMARY",
      "────────────────────────────",
      "Design: " + d.name,
      "Hull: L " + fmt(d.hullLength_m, 2) + " m, D " + fmt(d.hullDiameter_m, 3) + " m, CG " + fmt(d.xcg_m, 2) + " m from nose",
      "Front pair: S " + fmt(d.frontArea_m2, 3) + " m², AR " + fmt(d.frontAR, 2) + ", x " + fmt(d.frontX_m, 2) + " m",
      "Rear pair:  S " + fmt(d.rearArea_m2, 3) + " m², AR " + fmt(d.rearAR, 2) + ", x " + fmt(d.rearX_m, 2) + " m" + (d.cruciform ? " (cruciform)" : ""),
      "Waterjet: " + fmt(d.maxPower_kW, 1) + " kW max, η_jet " + fmt(r.ctx.etaJet, 2) + ", nozzle " + fmt(d.nozzleDiameter_m * 1000, 0) + " mm, BSFC " + fmt(d.bsfc_kgpkWh, 2) + " kg/kWh",
      "Mass: m0 " + fmt(r.ctx.m0, 1) + " kg (fuel " + fmt(d.fuelMass_kg, 1) + " kg, payload " + fmt(state.payloadKg, 1) + " kg), reserve " + fmt(e.reserveFraction * 100, 0) + " %",
      "Water density " + fmt(e.density_kgpm3, 1) + " kg/m³ · cruise " + fmt(c.V, 2) + " m/s",
      "Turbulence: σ_u/v/w " + t.sigma_u + "/" + t.sigma_v + "/" + t.sigma_w + " m/s, L_u/v/w " + t.L_u + "/" + t.L_v + "/" + t.L_w + " m" + (state.model.turbulence.placeholder ? "  [ATMOSPHERIC PLACEHOLDER]" : ""),
      "",
      "One-way range: " + fmt(c.rangeKm, 1) + " km · endurance " + fmt(c.enduranceH, 2) + " h",
      "Best range: " + (r.best ? fmt(r.best.rangeKm, 1) + " km at " + fmt(r.best.V, 2) + " m/s" : "–") + " · max speed " + (r.vMax != null ? fmt(r.vMax, 2) + " m/s" : "–"),
      "Drag " + fmt(c.drag.total, 1) + " N · engine " + fmt(c.P_engine / 1000, 2) + " kW (" + fmt(c.throttle * 100, 0) + " %) · fuel " + fmt(c.fuelFlow_kgph, 2) + " kg/h",
      "Static margin " + fmt(r.staticMargin_m, 3) + " m (" + fmt(r.staticMargin_pct, 1) + " % L), NP " + fmt(r.xnp, 3) + " m",
      "RMS gust incidence front/rear " + fmt(c.sv.sigFinAlpha[0] * DEG, 2) + "°/" + fmt(c.sv.sigFinAlpha[1] * DEG, 2) + "° · load factor 1 ± " + fmt(c.sv.sigL / r.ctx.W, 2),
      "Track wander per 1 km (depth/lateral): " + fmt(r.per1km.depth_m, 1) + " / " + fmt(r.per1km.lateral_m, 1) + " m",
      "Start point: " + (state.mapCenter ? fmt(state.mapCenter.lat, 4) + ", " + fmt(state.mapCenter.lng, 4) : "default"),
      "Warnings: " + (r.warnings.length ? r.warnings.join(" | ") : "none")
    ].join("\n");
  }

  // ---------- Custom lists ----------
  function renderCustomLists() {
    const vl = $("customVehicleList"); vl.innerHTML = "";
    state.customVehicles.forEach(v => {
      const row = document.createElement("div"); row.className = "cust-row";
      row.innerHTML = "<span><b>" + esc(v.name) + "</b> · L " + fmt(v.hullLength_m, 2) + " m · " + fmt(v.maxPower_kW, 1) + " kW</span><button class='btn-link' data-remove-vehicle='" + esc(v.id) + "'>Remove</button>";
      vl.appendChild(row);
    });
  }

  // ---------- Import / export ----------
  function handleImport(text) {
    const content = text.trim().replace(/^﻿/, "");
    const format = (content.startsWith("[") || content.startsWith("{")) ? "json" : "csv";
    let added = 0;
    if (content.indexOf("hullLength_m") >= 0) {
      const list = Data.parseVehiclesText(content, format);
      list.forEach(v => { const i = state.customVehicles.findIndex(x => x.id === v.id); if (i >= 0) state.customVehicles[i] = v; else state.customVehicles.push(v); added++; });
      buildVehicleSelect();
    }
    if (content.indexOf("weightKg") >= 0) {
      let list = [];
      if (format === "json") { const j = JSON.parse(content); list = Data.parsePayloadsText(JSON.stringify(Array.isArray(j) ? j.filter(x => x.weightKg != null) : (j.payloads || [])), "json"); }
      else list = Data.parsePayloadsText(content, format);
      list.forEach(l => { const i = state.customPayloads.findIndex(x => x.id === l.id); if (i >= 0) state.customPayloads[i] = l; else state.customPayloads.push(l); added++; });
      buildPayloadSelect();
    }
    if (!added) throw new Error("Could not tell whether the file contains vehicles or payloads.");
    renderCustomLists(); schedule(); flash("Imported " + added + " item(s).");
  }
  function exportEverything() {
    const f = $("dlg_format").value || "json", delim = f === "tsv" ? "\t" : ";", ext = f === "tsv" ? "tsv" : "csv";
    if (f === "json") Data.download("jetify-data.json", JSON.stringify({ vehicles: allVehicles(), payloads: allPayloads() }, null, 2), "application/json");
    else {
      Data.download("jetify-vehicles." + ext, Data.serializeVehicles(allVehicles(), f, delim), "text/csv;charset=utf-8");
      Data.download("jetify-payloads." + ext, Data.serializePayloads(allPayloads(), f, delim), "text/csv;charset=utf-8");
    }
  }
  function downloadTemplates() {
    const f = $("dlg_format").value || "json", delim = f === "tsv" ? "\t" : ";", ext = f === "tsv" ? "tsv" : "csv";
    if (f === "json") {
      const blank = {}; Data.vehicleHeaders.forEach(h => blank[h] = h === "cruciform" ? true : (["id", "name", "type"].includes(h) ? "" : 0));
      Data.download("template-vehicles.json", JSON.stringify([blank], null, 2), "application/json");
      Data.download("template-payloads.json", JSON.stringify([{ id: "", name: "", type: "", weightKg: 0 }], null, 2), "application/json");
    } else {
      Data.download("template-vehicles." + ext, Data.serializeVehicles([], f, delim), "text/csv;charset=utf-8");
      Data.download("template-payloads." + ext, Data.serializePayloads([], f, delim), "text/csv;charset=utf-8");
    }
    flash("Templates downloaded.");
  }

  function flash(msg) {
    const el = $("flash");
    el.textContent = msg; show(el);
    clearTimeout(el._t); el._t = setTimeout(() => hide(el), 2500);
  }

  // ---------- Events ----------
  function bindStaticEvents() {
    const on = (id, ev, fn) => { const el = $(id); if (el) el.addEventListener(ev, fn); };

    on("vehicleSelect", "change", () => { loadProfile(allVehicles().find(v => v.id === $("vehicleSelect").value)); schedule(); });
    on("resetDesignBtn", "click", () => { loadProfile(selectedProfile()); schedule(); });
    on("saveDesignBtn", "click", () => {
      const name = $("saveName").value.trim() || (state.design.name + " (edited)");
      const v = { ...JSON.parse(JSON.stringify(state.design)), id: "custom-" + Date.now(), name, type: "custom" };
      state.customVehicles.push(v);
      state.activeVehicleId = v.id;
      buildVehicleSelect(); loadProfile(v); schedule();
      flash("Saved: " + name);
    });
    on("cruciform", "change", () => { state.design.cruciform = $("cruciform").checked; schedule(); });
    on("useVehicleResponse", "change", () => { state.env.useVehicleResponse = $("useVehicleResponse").checked; schedule(); });
    on("includeQPenalty", "change", () => { state.env.includeQPenalty = $("includeQPenalty").checked; schedule(); });
    on("resetTurbBtn", "click", () => { state.env.turbulence = { ...state.model.turbulence }; renderFields("turbulence", "box_turbulence", "turb"); schedule(); });

    on("payloadSelect", "change", () => {
      state.payloadId = $("payloadSelect").value;
      const l = allPayloads().find(x => x.id === state.payloadId);
      if (l) state.payloadKg = l.weightKg;
      $("payloadWeight").value = state.payloadKg;
      $("payloadWeightLabel").textContent = fmt(state.payloadKg, 1) + " kg";
      schedule();
    });
    on("payloadWeight", "input", () => {
      state.payloadKg = parseFloat($("payloadWeight").value);
      state.payloadId = "__custom"; $("payloadSelect").value = "__custom";
      $("payloadWeightLabel").textContent = fmt(state.payloadKg, 1) + " kg";
      schedule();
    });

    ["cal_hull", "cal_claf", "cal_clar", "cal_eta", "cal_extra", "cal_runP", "cal_runV"].forEach(id => on(id, "input", schedule));
    ["cal_hull_use", "cal_claf_use", "cal_clar_use", "cal_eta_use", "cal_extra_use", "cal_run_use"].forEach(id => on(id, "change", schedule));
    on("cal_resetBtn", "click", () => {
      ["cal_hull", "cal_claf", "cal_clar", "cal_eta", "cal_extra", "cal_runP", "cal_runV"].forEach(id => { $(id).value = ""; });
      ["cal_hull_use", "cal_claf_use", "cal_clar_use", "cal_eta_use", "cal_extra_use", "cal_run_use"].forEach(id => { $(id).checked = false; });
      schedule();
    });

    on("computeBtn", "click", () => { const r = compute(); if (r) addHistoryEntry(snapshot(r)); });
    on("addScenarioBtn", "click", addScenario);
    on("clearScenariosBtn", "click", () => { state.scenarios = []; renderScenarios(); });
    on("clearHistoryBtn", "click", () => { state.history = []; try { sessionStorage.removeItem("jetify_history"); } catch (e) { /* ignore */ } renderHistory(); });
    on("copySummaryBtn", "click", async () => {
      const r = compute(); if (!r) return;
      try { await navigator.clipboard.writeText(summaryText(r)); flash("Summary copied."); } catch (e) { flash("Could not copy."); }
    });
    on("downloadResultsBtn", "click", () => { const r = compute(); if (r) Data.download("jetify-results.txt", summaryText(r), "text/plain;charset=utf-8"); });

    on("resetStartBtn", "click", () => {
      if (!state.map || !state.mapConfig) return;
      state.mapTouched = false;
      state.mapCenter = { lat: state.mapConfig.center.lat, lng: state.mapConfig.center.lng };
      state.marker.setLatLng([state.mapCenter.lat, state.mapCenter.lng]);
      schedule();
    });

    on("btnDownloadTemplates", "click", downloadTemplates);
    on("btnExportAll", "click", exportEverything);
    on("fileInput", "change", (ev) => {
      const files = Array.from(ev.target.files); if (!files.length) return;
      let pending = files.length; const errors = [];
      files.forEach(f => {
        const reader = new FileReader();
        reader.onload = () => {
          try { handleImport(String(reader.result)); } catch (e) { errors.push(f.name + ": " + e.message); }
          if (--pending === 0) { if (errors.length) flash("Error: " + errors.join(" | ")); ev.target.value = ""; }
        };
        reader.readAsText(f);
      });
    });
    on("retryBtn", "click", () => location.reload());

    document.addEventListener("click", (e) => {
      const rv = e.target.closest("[data-remove-vehicle]");
      if (!rv) return;
      const id = rv.getAttribute("data-remove-vehicle");
      state.customVehicles = state.customVehicles.filter(x => x.id !== id);
      if (state.activeVehicleId === id) state.activeVehicleId = null;
      buildVehicleSelect(); loadProfile(selectedProfile()); schedule();
    });
  }

  document.addEventListener("DOMContentLoaded", init);
})();
