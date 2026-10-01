"use strict";
// Run: node --test test/
const test = require("node:test");
const assert = require("node:assert");
const fs = require("node:fs");
const path = require("node:path");
const Model = require("../js/model.js");
const I = Model._internal;

const model = JSON.parse(fs.readFileSync(path.join(__dirname, "../data/model.json"), "utf8"));
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
  const a = Model.analyze({ ...d, rearX_m: 1.4 }, envOf({}), model, {});
  const b = Model.analyze({ ...d, rearX_m: 1.7 }, envOf({}), model, {});
  assert.ok(b.staticMargin_m > a.staticMargin_m);
});

test("zero turbulence → no gust drag, and penalty grows with σ_w", () => {
  const d = vehicles[0];
  const calm = Model.analyze(d, envOf({ turbulence: { ...model.turbulence, sigma_u: 0, sigma_v: 0, sigma_w: 0 } }), model, {});
  assert.ok(calm.chosen.drag.gustV < 1e-9 && calm.chosen.drag.gustL < 1e-9);
  const base = Model.analyze(d, envOf({}), model, {});
  assert.ok(base.chosen.rangeKm < calm.chosen.rangeKm);
});
