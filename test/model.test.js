"use strict";
// Run: node --test test/
const test = require("node:test");
const assert = require("node:assert");
const fs = require("node:fs");
const path = require("node:path");
const Model = require("../js/model.js");
const I = Model._internal;

const model = JSON.parse(fs.readFileSync(path.join(__dirname, "../data/model.json"), "utf8"));
const variables = JSON.parse(fs.readFileSync(path.join(__dirname, "../data/variables.json"), "utf8"));
model.presets = Model.presetsFromVariables(variables);
const vehicles = JSON.parse(fs.readFileSync(path.join(__dirname, "../data/vehicles.json"), "utf8"));
const envOf = (over) => ({ ...model.environmentDefaults, payload_kg: 0, turbulence: { ...model.turbulence }, ...over });

const integrate = (f, n = 4000) => { const g = I.logGrid(1e-7, 1e5, n); let s = 0; for (let k = 0; k < g.n; k++) s += f(g.om[k]) * g.w[k]; return s; };

test("von Kármán spectra integrate to σ²", () => {
  const sT = integrate(Om => I.phiTransverse(Om, 3, 50));
  const sL = integrate(Om => I.phiLongitudinal(Om, 5, 200));
  assert.ok(Math.abs(sT / 9 - 1) < 0.01, "transverse " + sT);
  assert.ok(Math.abs(sL / 25 - 1) < 0.01, "longitudinal " + sL);
});

test("Lamb added-mass coefficients match the prolate-spheroid table (λ = 2, 4)", () => {
  const a = I.lambCoefficients(2), b = I.lambCoefficients(4);
  assert.ok(Math.abs(a.k1 - 0.210) < 0.002 && Math.abs(a.k2 - 0.702) < 0.003 && Math.abs(a.kp - 0.240) < 0.003, JSON.stringify(a));
  assert.ok(Math.abs(b.k1 - 0.082) < 0.002 && Math.abs(b.k2 - 0.860) < 0.003 && Math.abs(b.kp - 0.608) < 0.003, JSON.stringify(b));
});

test("waterjet: thrust ↔ jet power are inverse", () => {
  const rho = 1005, An = Math.PI / 4 * 0.06 ** 2, V = 10, T = 400;
  const j = I.jetForThrust(T, V, rho, An);
  const back = I.thrustForJetPower(j.Pjet, V, rho, An);
  assert.ok(Math.abs(back.T / T - 1) < 1e-6);
  assert.ok(Math.abs(j.etaF - 2 * V / (V + j.Vj)) < 1e-9);
});

test("closed-form gust drag penalty equals the MATLAB script's 120 s Monte-Carlo mean", () => {
  // Reference script values (air, interceptor) – used only to validate the statistics.
  const rho = 1.225, g = 9.81, S = 0.20, AR = 2.5, e = 0.8, CD0 = 0.12, m0 = 50, mf = 35;
  // dt finer than the script's 0.02 s: Euler at dt/τ ≈ 0.07 inflates the variance by ~4 %.
  const V0 = 180, sw = 3, Lw = 50, CLa = 4.5, dt = 0.002;
  const W = m0 * g, CLt = W / (0.5 * rho * V0 * V0 * S);
  // Seeded Gaussian RNG
  let s = 12345;   // mulberry32
  const rnd = () => { s = (s + 0x6D2B79F5) | 0; let t = Math.imul(s ^ (s >>> 15), 1 | s); t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t; return ((t ^ (t >>> 14)) >>> 0) / 4294967296; };
  const randn = () => Math.sqrt(-2 * Math.log(rnd() + 1e-12)) * Math.cos(2 * Math.PI * rnd());
  const tau = Lw / V0, N = Math.round(4000 / dt);   // longer than 120 s to beat the noise
  let w = 0, sumCL = 0, sumCD = 0;
  for (let k = 0; k < N; k++) {
    w = w - dt / tau * w + Math.sqrt(2 * sw * sw * dt / tau) * randn();
    const CL = CLt + CLa * Math.atan2(w, V0);
    sumCL += CL; sumCD += CD0 + CL * CL / (Math.PI * e * AR);
  }
  const LDmc = (sumCL / N) / (sumCD / N);
  // Frequency-domain / closed form: E[CD] = CD0 + (CLt² + CLα²σ_α²)/(π e AR), σ_α = σ_w/V0
  const sa = sw / V0;
  const LDcf = CLt / (CD0 + (CLt * CLt + CLa * CLa * sa * sa) / (Math.PI * e * AR));
  assert.ok(Math.abs(LDcf / LDmc - 1) < 0.01, `L/D closed form ${LDcf} vs MC ${LDmc}`);
  // Breguet equivalence: script R = V/c·L/D·ln(m0/mf) with c = 2/h
  const R = V0 / (2 / 3600) * LDcf * Math.log(m0 / mf) / 1000;
  assert.ok(R > 0);
});

test("baseline vehicle produces a finite, stable result", () => {
  const r = Model.analyze(vehicles[0], envOf({}), model, {});
  assert.ok(r && isFinite(r.chosen.rangeKm) && r.chosen.rangeKm > 0);
  assert.ok(r.staticMargin_m > 0, "baseline should be statically stable");
  assert.ok(r.chosen.vert.stable);
  assert.ok(isFinite(r.per1km.depth_m));
});

test("moving the rear pair aft increases the static margin", () => {
  const d = vehicles[0];
  const a = Model.analyze({ ...d, rearXle_m: 1.4 }, envOf({}), model, {});
  const b = Model.analyze({ ...d, rearXle_m: 1.7 }, envOf({}), model, {});
  assert.ok(b.staticMargin_m > a.staticMargin_m);
});

test("zero turbulence → no gust drag, and penalty grows with σ_w", () => {
  const d = vehicles[0];
  const calm = Model.analyze(d, envOf({ turbulence: { ...model.turbulence, sigma_u: 0, sigma_v: 0, sigma_w: 0 } }), model, {});
  assert.ok(calm.chosen.drag.gustV < 1e-9 && calm.chosen.drag.gustL < 1e-9);
  const base = Model.analyze(d, envOf({}), model, {});
  assert.ok(base.chosen.rangeKm < calm.chosen.rangeKm);
});

test("variables.json: every field is well formed and every slider span contains the defaults", () => {
  const all = Object.values(variables.groups).flat();
  const keys = new Set(all.map(d => d.key));
  assert.strictEqual(keys.size, all.length, "duplicate key");
  const envDefaults = { ...model.environmentDefaults };
  all.forEach(d => {
    assert.ok(d.key && d.label && ["slider", "number"].includes(d.kind), "bad entry " + JSON.stringify(d));
    assert.ok(d.step > 0, d.key + ": step must be > 0");
    if (d.kind !== "slider") return;
    const max = typeof d.max === "string" ? null : d.max;
    if (typeof d.max === "string") assert.ok(keys.has(d.max), d.key + ": max refers to unknown key " + d.max);
    else assert.ok(d.min < max, d.key + ": min must be < max");
    const values = [envDefaults[d.key], ...vehicles.map(v => v[d.key])].filter(v => v != null);
    values.forEach(v => {
      const hi = max != null ? max : Infinity;
      assert.ok(v >= d.min && v <= hi, `${d.key}: default ${v} is outside the span ${d.min}–${d.max}`);
    });
  });
});

test("generalised Helmbold: a0 = 2π gives the classic form, large AR tends to a0", () => {
  const AR = 3;
  assert.ok(Math.abs(I.helmbold(AR) - 2 * Math.PI * AR / (2 + Math.sqrt(AR * AR + 4))) < 1e-12);
  assert.ok(Math.abs(I.helmbold(AR, 2 * Math.PI) - I.helmbold(AR)) < 1e-12);
  assert.ok(Math.abs(I.helmbold(1e5, 5.7) / 5.7 - 1) < 1e-3);
  assert.ok(I.helmbold(AR, 5.7) < I.helmbold(AR));
});

test("section data: log-Re interpolation and clamping", () => {
  const sec = { points: [{ Re: 1e5, clAlpha_perRad: 5, cdMin: 0.01, clMax: 1 }, { Re: 1e7, clAlpha_perRad: 7, cdMin: null, clMax: 2 }] };
  const mid = I.sectionAt(sec, 1e6);
  assert.ok(Math.abs(mid.clAlpha - 6) < 1e-12 && Math.abs(mid.clMax - 1.5) < 1e-12 && mid.cdMin === null);
  assert.strictEqual(I.sectionAt(sec, 1e3).clAlpha, 5);
  assert.strictEqual(I.sectionAt(sec, 1e9).clAlpha, 7);
});

test("sections.json is valid and the ideal section reproduces the old model", () => {
  const sections = JSON.parse(fs.readFileSync(path.join(__dirname, "../data/sections.json"), "utf8"));
  assert.ok(sections.sections.some(s => s.id === sections.default));
  sections.sections.forEach(s => {
    assert.ok(s.id && s.name && s.points.length, s.id);
    s.points.forEach(p => assert.ok(p.Re > 0 && p.clAlpha_perRad > 0, s.id));
  });
  const v = vehicles[0], env = envOf({});
  const old = Model.analyze(v, env, model);
  const withIdeal = Model.analyze({ ...v, frontSection: "ideal", rearSection: "ideal" }, env, { ...model, sections });
  assert.ok(Math.abs(withIdeal.chosen.rangeKm / old.chosen.rangeKm - 1) < 1e-12);
  const n12 = Model.analyze({ ...v, frontSection: "naca0012", rearSection: "naca0012" }, env, { ...model, sections });
  assert.ok(n12.ctx.rear.clAlpha < old.ctx.rear.clAlpha, "real section has a lower slope than 2π");
  assert.notStrictEqual(n12.ctx.rear.cd0, env.finCD0);
});

test("planform: rectangle and straight-TE delta geometry", () => {
  const r = I.planform(0.2, 4, 1, 0, 1);
  const c = Math.sqrt(0.2 / 4);
  assert.ok(Math.abs(r.mac - c) < 1e-12 && Math.abs(r.xAc25 - (1 + c / 4)) < 1e-12);
  // Delta, λ = 0, tanΛ_LE = 4/AR (straight trailing edge): MAC = ⅔c_r, ¼-MAC at ½c_r, centroid at ⅔c_r.
  const AR = 2.3, d = I.planform(0.3, AR, 0, Math.atan(4 / AR) * 180 / Math.PI, 0);
  assert.ok(Math.abs(d.mac - 2 / 3 * d.cr) < 1e-12);
  assert.ok(Math.abs(d.xAc25 - d.cr / 2) < 1e-9, "a.c. " + d.xAc25 / d.cr);
  assert.ok(Math.abs(d.xCentroid - 2 / 3 * d.cr) < 1e-9);
  assert.ok(Math.abs(d.tanHc - (4 / AR - 2 / AR)) < 1e-12);
});

test("legacy centre-of-lift x converts without changing the result", () => {
  const v = vehicles[0], env = envOf({});
  const pf = (n) => I.planform(v[n + "Area_m2"], v[n + "AR"], v[n + "Taper"], v[n + "SweepLE_deg"], v[n + "Xle_m"]);
  const legacy = { ...v, frontX_m: pf("front").xAc25, rearX_m: pf("rear").xAc25 };
  delete legacy.frontXle_m; delete legacy.rearXle_m;
  const a = Model.analyze(v, env, model), b = Model.analyze(legacy, env, model);
  assert.ok(Math.abs(a.chosen.rangeKm / b.chosen.rangeKm - 1) < 1e-9);
  assert.ok(Math.abs(a.staticMargin_m - b.staticMargin_m) < 1e-9);
});

test("sweep lowers the DATCOM lift slope", () => {
  assert.ok(I.helmbold(6, 2 * Math.PI, Math.tan(30 / 57.2958)) < I.helmbold(6, 2 * Math.PI, 0));
});

test("delta (Polhamus): small-α slope = Kp, vortex lift adds slope, CDi = CL·tanα", () => {
  const pair = { kind: "delta", clAlpha: 2.4, Kv: Math.PI, stallAlpha: 25 / 57.2958, pf: I.planform(0.3, 2.3, 0, 60, 0) };
  const s0 = I.liftState(pair, 1e-6);
  assert.ok(Math.abs(s0.slope / 2.4 - 1) < 1e-3, "slope at 0: " + s0.slope);
  const s = I.liftState(pair, 0.4);
  assert.ok(Math.abs(I.deltaCL(pair, s.alpha) - 0.4) < 1e-9);
  assert.ok(s.slope > 2.4, "vortex lift should raise the slope");
  assert.ok(Math.abs(s.CDi - 0.4 * Math.tan(s.alpha)) < 1e-12);
  assert.ok(s.xAc > pair.pf.xAc25 && s.xAc < pair.pf.xCentroid, "a.c. moves aft towards the centroid");
});

test("downwash: plausible gradient for a conventional wing + tail, and it lowers the static margin", () => {
  const ctx = { downwash: true, tailHeight: 0, front: { AR: 7, clAlpha: 4.8, pf: I.planform(0.5, 7, 0.5, 0, 0) } };
  const b = ctx.front.pf.b;
  const E = I.downwashGradient(ctx, { xAc: 0, slope: 4.8 }, { xAc: 0.5 * b }).E;
  assert.ok(E > 0.25 && E < 0.6, "dε/dα " + E);
  const v = vehicles[0], env = envOf({});
  const on = Model.analyze({ ...v, downwash: true }, env, model), off = Model.analyze({ ...v, downwash: false }, env, model);
  assert.ok(on.chosen.trim.dw.E > 0);
  assert.ok(on.staticMargin_m < off.staticMargin_m);
});

test("delta wing + aft tail produces a finite result", () => {
  const v = { ...vehicles[0], hullLength_m: 1.2, hullDiameter_m: 0.09, xcg_m: 0.62, emptyMass_kg: 1.8, fuelMass_kg: 0.2,
    frontArea_m2: 0.3, frontAR: 2.3, frontTaper: 0, frontSweepLE_deg: 60, frontXle_m: 0.15, frontPlanform: "delta",
    rearArea_m2: 0.06, rearAR: 4, rearTaper: 0.6, rearSweepLE_deg: 20, rearXle_m: 1.05, rearHeight_m: 0.05, cruciform: false };
  const r = Model.analyze(v, envOf({ cruiseSpeed_mps: 20 }), model);
  assert.ok(r && Number.isFinite(r.staticMargin_m) && Number.isFinite(r.chosen.drag.total));
  assert.strictEqual(r.ctx.front.kind, "delta");
});
