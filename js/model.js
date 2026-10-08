"use strict";

/*
 * Jetify physics engine – one-way range and trajectory stability for a
 * waterjet-propelled, jet-fuel AUV. Quasi-steady; no time simulation.
 *
 * Every coefficient is read from data/model.json, the vehicle design and the
 * environment inputs (and can be overridden by the calibration panel).
 * See docs/MODEL.md for the full derivation.
 *
 *   Static balance (as in the MATLAB reference script, buoyancy deliberately
 *   excluded):  L_front + L_rear = W = m0·g, split by moment balance about CG.
 *
 *   Drag:  hull  = q·S_wet·Cf(Re)·FF(λ)            (ITTC-57 + Hoerner form factor)
 *          fins  = q·S_i·[CD0 + E[CL_i²]/(π·e·AR_i)]
 *          E[CL_i²] = CL_i,trim² + CLα_i²·σ²_α,i     (von Kármán gust penalty)
 *
 *   σ²_α,i = ∫ |G_α,i(Ω)|² Φ_vK(Ω) dΩ  – integrated in the frequency domain.
 *   G includes fin-chord averaging, the fore/aft gust delay Δx and (optionally)
 *   the rigid-body pitch/heave response of the vehicle.
 *
 *   Waterjet:  T = ρ·A_n·Vj·(Vj − V),  P_jet = ½·ρ·A_n·Vj·(Vj² − V²),
 *              P_engine = P_jet / η_jet,  η_Froude = 2V/(V + Vj)
 *   Range (Breguet, fuel burn):  R = η_jet·η_F/(BSFC·g) · (L/D) · ln(m0/m_end)
 *   (equivalent to the script's  R = V/c · L/D · ln(m0/mf)  with c = BSFC·g·V/η)
 */
const Model = (() => {

  const round = (v, d) => { const m = Math.pow(10, d); return Math.round(v * m) / m; };
  const DEG = 180 / Math.PI;

  // ---------- Small complex helpers ----------
  const cx = (re, im) => ({ re, im: im || 0 });
  const cadd = (a, b) => cx(a.re + b.re, a.im + b.im);
  const cmul = (a, b) => cx(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re);
  const cscale = (a, s) => cx(a.re * s, a.im * s);
  const cdiv = (a, b) => { const d = b.re * b.re + b.im * b.im; return cx((a.re * b.re + a.im * b.im) / d, (a.im * b.re - a.re * b.im) / d); };
  const cabs2 = (a) => a.re * a.re + a.im * a.im;
  const cexpi = (phi) => cx(Math.cos(phi), Math.sin(phi));
  const sinc = (x) => Math.abs(x) < 1e-6 ? 1 - x * x / 6 : Math.sin(x) / x;

  // ---------- Von Kármán spectra (one-sided, spatial frequency Ω [rad/m]) ----------
  // ∫₀^∞ Φ dΩ = σ²  for both forms (MIL-F-8785C / MIL-HDBK-1797).
  const VK = 1.339;
  function phiLongitudinal(Om, sigma, L) {
    const x = VK * L * Om;
    return sigma * sigma * (2 * L / Math.PI) / Math.pow(1 + x * x, 5 / 6);
  }
  function phiTransverse(Om, sigma, L) {
    const x = VK * L * Om;
    return sigma * sigma * (L / Math.PI) * (1 + (8 / 3) * x * x) / Math.pow(1 + x * x, 11 / 6);
  }

  // Log-spaced grid with trapezoid weights in ln Ω  (∫f dΩ ≈ Σ f(Ω_k)·w_k).
  function logGrid(omMin, omMax, n) {
    const a = Math.log(omMin), b = Math.log(omMax), h = (b - a) / (n - 1);
    const om = new Float64Array(n), w = new Float64Array(n);
    for (let k = 0; k < n; k++) {
      om[k] = Math.exp(a + k * h);
      w[k] = om[k] * h * ((k === 0 || k === n - 1) ? 0.5 : 1);
    }
    return { om, w, n };
  }

  // ---------- Hull: geometry, added mass (Lamb's k-factors), drag ----------
  // Prolate spheroid of the same fineness ratio λ = L/D (Lamb, Hydrodynamics §373).
  function lambCoefficients(lam) {
    if (lam <= 1.0001) return { k1: 0.5, k2: 0.5, kp: 0 };
    const e = Math.sqrt(1 - 1 / (lam * lam));
    const lnT = Math.log((1 + e) / (1 - e));
    const e3 = e * e * e;
    const a0 = 2 * (1 - e * e) / e3 * (0.5 * lnT - e);
    const b0 = 1 / (e * e) - (1 - e * e) / (2 * e3) * lnT;
    const k1 = a0 / (2 - a0);
    const k2 = b0 / (2 - b0);
    const kp = Math.pow(e, 4) * (b0 - a0) / ((2 - e * e) * (2 * e * e - (2 - e * e) * (b0 - a0)));
    return { k1, k2, kp };
  }

  function hullGeometry(d) {
    const L = d.hullLength_m, D = d.hullDiameter_m;
    const lam = L / D;
    const volume = d.volumeCoeff * Math.PI / 4 * D * D * L;
    const wetted = d.wettedCoeff * Math.PI * D * L;
    const lamb = lambCoefficients(lam);
    const formFactor = 1 + 1.5 / Math.pow(lam, 1.5) + 7 / Math.pow(lam, 3);   // Hoerner, bodies of revolution
    return { L, D, lam, volume, wetted, lamb, formFactor };
  }

  // ITTC-1957 friction line (valid for turbulent flow; Re clamped to 1e5 at very low speed).
  function ittcCf(Re) {
    const l = Math.log10(Math.max(Re, 1e5)) - 2;
    return 0.075 / (l * l);
  }

  // ---------- Atmosphere (ISA, 0–20 km) ----------
  // Troposphere T = 288.15 − 0.0065·h, isothermal 216.65 K above 11 km; Sutherland viscosity.
  function isa(h) {
    const R = 287.05, g0 = 9.80665;
    const hh = Math.max(0, Math.min(h || 0, 20000));
    let T, p;
    if (hh <= 11000) { T = 288.15 - 0.0065 * hh; p = 101325 * Math.pow(T / 288.15, g0 / (0.0065 * R)); }
    else { T = 216.65; p = 22632.06 * Math.exp(-g0 * (hh - 11000) / (R * T)); }
    const rho = p / (R * T);
    const mu = 1.458e-6 * Math.pow(T, 1.5) / (T + 110.4);
    return { h: hh, T, p, rho, mu, nu: mu / rho, a: Math.sqrt(1.4 * R * T), sigma: rho / 1.225 };
  }

  // useISA (default): density, viscosity and speed of sound from the altitude; otherwise the
  // manual values (water: ρ ≈ 1005, ν ≈ 1.3e-6, a ≈ 1480).
  function resolveAtmosphere(env) {
    if (env.useISA === false) return { ...env, speedOfSound_mps: env.speedOfSound_mps || 1480, atm: null };
    const atm = isa(env.altitude_m);
    return { ...env, density_kgpm3: atm.rho, kinematicViscosity_m2ps: atm.nu, speedOfSound_mps: atm.a, atm };
  }

  // ---------- Control-surface pairs ----------
  // DATCOM lift-curve slope [1/rad] with Prandtl–Glauert compressibility, β = √(1 − M²):
  //   CLα = 2πAR / (2 + √(AR²β²/κ²·(1 + tan²Λ½c/β²) + 4)),  κ = a0/2π (a0 = incompressible section slope).
  // M = 0, Λ = 0, a0 = 2π gives the classic Helmbold 2πAR/(2+√(AR²+4)). M is capped at MACH_CAP.
  const MACH_CAP = 0.95;
  const helmbold = (AR, a0, tanHalfChord, M) => {
    const k = (a0 != null ? a0 : 2 * Math.PI) / (2 * Math.PI);
    const t = tanHalfChord || 0;
    const m = Math.min(Math.max(M || 0, 0), MACH_CAP), b2 = 1 - m * m;
    return 2 * Math.PI * AR / (2 + Math.sqrt(AR * AR * b2 / (k * k) * (1 + t * t / b2) + 4));
  };

  // Drag divergence (Korn equation, Mason's form) and wave drag (Lock's 4th-power law):
  //   M_dd = κ_A/cosΛ − (t/c)/cos²Λ − CL/(10·cos³Λ),  M_crit = M_dd − (0.1/80)^(1/3),
  //   CD_wave = 20·(M − M_crit)⁴ for M > M_crit.   κ_A ≈ 0.87 conventional, 0.95 supercritical.
  function waveDrag(pair, CL, M) {
    const c = 1 / Math.sqrt(1 + pair.pf.tanQc * pair.pf.tanQc);
    const Mdd = pair.kappaA / c - pair.tc / (c * c) - Math.abs(CL) / (10 * c * c * c);
    const Mcrit = Mdd - Math.cbrt(0.1 / 80);
    return { Mdd, Mcrit, CD: M > Mcrit ? 20 * Math.pow(M - Mcrit, 4) : 0 };
  }

  // Section thickness and Korn factor (ideal 2π section: thickness from model.json).
  function compressData(sec, model) {
    const cm = (model && model.compressibility) || {};
    const tc = sec && sec.thickness_pct > 0 ? sec.thickness_pct / 100 : (cm.idealSectionThickness != null ? cm.idealSectionThickness : 0.1);
    const kappaA = sec && sec.kornKappa != null ? sec.kornKappa : (cm.defaultKornKappa != null ? cm.defaultKornKappa : 0.87);
    return { tc, kappaA };
  }

  // 2D section data at chord Reynolds number Re: linear in log10(Re), clamped to the data range.
  function sectionAt(sec, Re) {
    const pts = sec.points.slice().sort((a, b) => a.Re - b.Re);
    const pick = (p) => ({ clAlpha: p.clAlpha_perRad, cdMin: p.cdMin, clMax: p.clMax });
    if (pts.length === 1 || !(Re > pts[0].Re)) return pick(pts[0]);
    if (Re >= pts[pts.length - 1].Re) return pick(pts[pts.length - 1]);
    const i = pts.findIndex(p => p.Re > Re), a = pts[i - 1], b = pts[i];
    const t = (Math.log10(Re) - Math.log10(a.Re)) / (Math.log10(b.Re) - Math.log10(a.Re));
    const lerp = (x, y) => x == null || y == null ? null : x + t * (y - x);
    return { clAlpha: lerp(a.clAlpha_perRad, b.clAlpha_perRad), cdMin: lerp(a.cdMin, b.cdMin), clMax: lerp(a.clMax, b.clMax) };
  }

  // 3D CLmax ≈ 0.9·cℓmax·cosΛ¼c (Raymer); conservative for low-AR fins, which stall later.
  const CLMAX_3D = 0.9;

  // Trapezoidal planform of one pair (both halves together): area S, aspect ratio AR,
  // taper λ = c_tip/c_root, leading-edge sweep Λ_LE, root leading edge (apex) at x_le.
  function planform(S, AR, taper, sweepLEdeg, xle) {
    const lam = Math.max(0, Math.min(1, taper));
    const b = Math.sqrt(S * AR);
    const cr = 2 * S / (b * (1 + lam)), ct = lam * cr;
    const mac = 2 / 3 * cr * (1 + lam + lam * lam) / (1 + lam);
    const yMac = b / 6 * (1 + 2 * lam) / (1 + lam);
    const tanLE = Math.tan(sweepLEdeg / DEG);
    const tanAt = (n) => tanLE - 4 * n / AR * (1 - lam) / (1 + lam);   // sweep of the n-chord line
    const xAc25 = xle + yMac * tanLE + 0.25 * mac;                        // ¼-MAC (attached-flow a.c.)
    const xCentroid = xle + b / 2 * tanLE * (1 + 2 * lam) / (3 * (1 + lam)) + (cr * cr + cr * ct + ct * ct) / (3 * (cr + ct));
    return { b, cr, ct, mac, yMac, taper: lam, tanLE, sweepLEdeg, tanQc: tanAt(0.25), tanHc: tanAt(0.5), xle, xAc25, xCentroid };
  }

  const TAIL_TYPES = ["traditional", "vtail", "tailless", "cruciform"];

  // Fill planform defaults (model.json designDefaults) and convert legacy frontX_m / rearX_m
  // (old profiles gave the centre of lift of a rectangular pair) to the root leading edge.
  function normalizeDesign(design, model) {
    const defs = (model && model.designDefaults) || {};
    const d = { ...design };
    Object.keys(defs).forEach(k => { if (d[k] == null) d[k] = defs[k]; });
    // Older profiles: cruciform true → "cruciform"; false → "traditional" with no vertical fin.
    if (!TAIL_TYPES.includes(d.tailType)) d.tailType = d.cruciform === true ? "cruciform" : "traditional";
    if (d.finXle_m == null) d.finXle_m = d.rearXle_m;
    delete d.cruciform;
    ["front", "rear"].forEach(n => {
      if (d[n + "Xle_m"] == null && d[n + "X_m"] != null) {
        const pf = planform(d[n + "Area_m2"], d[n + "AR"], d[n + "Taper"] != null ? d[n + "Taper"] : 1, d[n + "SweepLE_deg"] || 0, 0);
        d[n + "Xle_m"] = d[n + "X_m"] - pf.xAc25;
      }
      delete d[n + "X_m"];
    });
    return d;
  }

  // Pair = planform + section + lift model. n = "front" | "rear".
  function finPair(n, design, clAlphaOverride, sec, env, model) {
    const area = design[n + "Area_m2"], AR = design[n + "AR"];
    const pf = planform(area, AR, design[n + "Taper"] != null ? design[n + "Taper"] : 1, design[n + "SweepLE_deg"] || 0, design[n + "Xle_m"]);
    const nu = env.kinematicViscosity_m2ps, V = env.cruiseSpeed_mps;
    const Re = V > 0 && nu > 0 ? V * pf.mac / nu : null;
    const s2 = sec ? sectionAt(sec, Re) : { clAlpha: 2 * Math.PI, cdMin: null, clMax: null };
    const M = env.speedOfSound_mps > 0 ? V / env.speedOfSound_mps : 0;
    const clAlpha = clAlphaOverride != null ? clAlphaOverride : helmbold(AR, s2.clAlpha, pf.tanHc, M);
    const clAlpha0 = clAlphaOverride != null ? clAlphaOverride : helmbold(AR, s2.clAlpha, pf.tanHc, 0);
    const cd0 = s2.cdMin != null ? s2.cdMin : env.finCD0;
    const e = design[n + "OswaldE"] != null && design[n + "OswaldE"] > 0 ? design[n + "OswaldE"] : env.oswaldE;
    const kind = design[n + "Planform"] === "delta" ? "delta" : "conventional";
    const dm = (model && model.delta) || {};
    const pair = {
      n, area, AR, pf, span: pf.b, chord: pf.mac, xcl: pf.xAc25, clAlpha, clAlpha0, M, cd0, e, kind, ...compressData(sec, model),
      Kv: dm.Kv != null ? dm.Kv : Math.PI, stallAlpha: (dm.stallAlphaDeg != null ? dm.stallAlphaDeg : 30) / DEG,
      section: { id: sec ? sec.id : "ideal", name: sec ? sec.name : "Ideal thin airfoil (2π)", confidence: sec ? sec.confidence : "typical", Re, a0: s2.clAlpha, cdMin: s2.cdMin, clMax2D: s2.clMax },
      clAlphaSrc: clAlphaOverride != null ? "empirical" : "DATCOM(AR, Λ½c, a0, M)"
    };
    const cosQc = 1 / Math.sqrt(1 + pf.tanQc * pf.tanQc);
    pair.clMax = kind === "delta" ? deltaCL(pair, pair.stallAlpha)
      : (s2.clMax != null ? CLMAX_3D * s2.clMax * cosQc : env.finCLmax);
    return pair;
  }

  // Single vertical fin (traditional / tailless). Geometry from its mirror-image pair
  // (area 2S, AR 2·AR_geo → same MAC, a.c. and sweep lines); lift slope from the effective
  // AR = factor·AR_geo, AR_geo = h²/S (endplate effect of fuselage / tailplane, Raymer ≈ 1.55).
  function verticalFin(design, sec, env, model) {
    const S = design.finArea_m2, ARg = design.finAR;
    if (!(S > 0) || !(ARg > 0)) return null;
    const vt = (model && model.verticalTail) || {};
    const ARe = (vt.effectiveARFactor != null ? vt.effectiveARFactor : 1.55) * ARg;
    const pf = planform(2 * S, 2 * ARg, design.finTaper != null ? design.finTaper : 1, design.finSweepLE_deg || 0, design.finXle_m);
    const nu = env.kinematicViscosity_m2ps, V = env.cruiseSpeed_mps;
    const Re = V > 0 && nu > 0 ? V * pf.mac / nu : null;
    const s2 = sec ? sectionAt(sec, Re) : { clAlpha: 2 * Math.PI, cdMin: null, clMax: null };
    const cosQc = 1 / Math.sqrt(1 + pf.tanQc * pf.tanQc);
    const M = env.speedOfSound_mps > 0 ? V / env.speedOfSound_mps : 0;
    return {
      n: "fin", area: S, AR: ARe, ARgeo: ARg, pf, span: pf.b / 2, chord: pf.mac, xcl: pf.xAc25, M, ...compressData(sec, model),
      clAlpha: helmbold(ARe, s2.clAlpha, pf.tanHc, M), clAlpha0: helmbold(ARe, s2.clAlpha, pf.tanHc, 0), cd0: s2.cdMin != null ? s2.cdMin : env.finCD0, e: env.oswaldE, kind: "conventional",
      clMax: s2.clMax != null ? CLMAX_3D * s2.clMax * cosQc : env.finCLmax,
      section: { id: sec ? sec.id : "ideal", name: sec ? sec.name : "Ideal thin airfoil (2π)", confidence: sec ? sec.confidence : "typical", Re, a0: s2.clAlpha, cdMin: s2.cdMin, clMax2D: s2.clMax },
      clAlphaSrc: "DATCOM(AR_eff, Λ½c, a0)"
    };
  }

  // Polhamus leading-edge-suction analogy (sharp LE, delta):
  //   CL = Kp·sinα·cos²α + Kv·cosα·sin²α     (Kp = attached-flow slope, Kv = vortex-lift factor)
  function deltaParts(pair, a) {
    const sn = Math.sin(a), cs = Math.cos(a);
    return { p: pair.clAlpha * sn * cs * cs, v: pair.Kv * cs * sn * Math.abs(sn) };
  }
  function deltaCL(pair, a) { const t = deltaParts(pair, a); return t.p + t.v; }

  // Lift state of a pair at trim lift coefficient CL:
  //   alpha, slope dCL/dα, lift centre xLift, perturbation a.c. xAc,
  //   CDi (lift-dependent drag) and d2 = d²CDi/dα² (gust penalty ½·d2·σα²).
  function liftState(pair, CL) {
    const pf = pair.pf;
    if (pair.kind !== "delta") {
      const a = pair.clAlpha, k = 1 / (Math.PI * pair.e * pair.AR);
      return { alpha: CL / a, slope: a, xLift: pf.xAc25, xAc: pf.xAc25, CDi: k * CL * CL, d2: 2 * k * a * a, stalled: Math.abs(CL) > pair.clMax };
    }
    // Delta: invert CL(α) by bisection (monotonic below ~45°); no leading-edge suction → CDi = CL·tanα.
    const aMax = 45 / DEG;
    let lo = -aMax, hi = aMax;
    for (let i = 0; i < 80; i++) { const m = (lo + hi) / 2; if (deltaCL(pair, m) < CL) lo = m; else hi = m; }
    const alpha = (lo + hi) / 2, h = 1e-4;
    const g = (a) => deltaCL(pair, a) * Math.tan(a);
    const slope = (deltaCL(pair, alpha + h) - deltaCL(pair, alpha - h)) / (2 * h);
    const d2 = (g(alpha + h) - 2 * g(alpha) + g(alpha - h)) / (h * h);
    const t = deltaParts(pair, alpha), tp = deltaParts(pair, alpha + h), tm = deltaParts(pair, alpha - h);
    const dp = (tp.p - tm.p) / (2 * h), dv = (tp.v - tm.v) / (2 * h);
    // Potential lift at ¼-MAC, vortex lift at the planform centroid (≈ ⅔ c_root for a delta).
    const xLift = Math.abs(t.p + t.v) > 1e-12 ? (t.p * pf.xAc25 + t.v * pf.xCentroid) / (t.p + t.v) : pf.xAc25;
    const xAc = Math.abs(dp + dv) > 1e-12 ? (dp * pf.xAc25 + dv * pf.xCentroid) / (dp + dv) : pf.xAc25;
    return { alpha, slope, xLift, xAc, CDi: CL * Math.tan(alpha), d2, stalled: Math.abs(alpha) > pair.stallAlpha };
  }

  // Downwash gradient dε/dα at the rear pair from the front pair (DATCOM / Raymer, low speed):
  //   dε/dα = 4.44·[K_A·K_λ·K_H·√cosΛ¼c]^1.19
  //   K_A = 1/AR − 1/(1 + AR^1.7),  K_λ = (10 − 3λ)/7,  K_H = (1 − |h_H/b|)/(2·l_H/b)^(1/3)
  // scaled by (actual slope / attached slope) so vortex lift on a delta strengthens the wake.
  function downwashGradient(ctx, fs, rs) {
    const lH = rs.xAc - fs.xAc;
    if (!ctx.downwash || !(lH > 0)) return { E: 0, lH };
    const f = ctx.front, pf = f.pf, AR = f.AR;
    const KA = 1 / AR - 1 / (1 + Math.pow(AR, 1.7));
    const Kl = (10 - 3 * pf.taper) / 7;
    const KH = (1 - Math.abs(ctx.tailHeight) / pf.b) / Math.cbrt(2 * lH / pf.b);
    if (!(KH > 0) || !(KA > 0)) return { E: 0, lH };
    const cosQc = 1 / Math.sqrt(1 + pf.tanQc * pf.tanQc);
    const E0 = 4.44 * Math.pow(KA * Kl * KH * Math.sqrt(cosQc), 1.19) * (fs.slope / f.clAlpha) * (f.clAlpha / f.clAlpha0);
    return { E: Math.min(E0, 0.95), E0, lH, capped: E0 > 0.95 };
  }

  // ---------- Waterjet ----------
  // Required thrust T at speed V → jet velocity, jet power, Froude efficiency.
  function jetForThrust(T, V, rho, An) {
    const Vj = (V + Math.sqrt(V * V + 4 * T / (rho * An))) / 2;
    const Pjet = 0.5 * rho * An * Vj * (Vj * Vj - V * V);
    return { Vj, Pjet, etaF: Pjet > 0 ? T * V / Pjet : 1 };
  }
  // Available thrust for a given jet power (monotone in Vj → bisection).
  function thrustForJetPower(Pjet, V, rho, An) {
    if (Pjet <= 0) return { T: 0, Vj: V };
    const P = (Vj) => 0.5 * rho * An * Vj * (Vj * Vj - V * V);
    let lo = V, hi = Math.max(2 * V, 1);
    while (P(hi) < Pjet) hi *= 2;
    for (let i = 0; i < 80; i++) { const m = (lo + hi) / 2; if (P(m) < Pjet) lo = m; else hi = m; }
    const Vj = (lo + hi) / 2;
    return { T: rho * An * Vj * (Vj - V), Vj };
  }

  // ---------- Rigid-body gust response in one plane (pitch/heave or yaw/sway) ----------
  // States α, q (incidence & body rate); γ̇ = ΔL/(m'V) is the path-angle rate.
  // Per unit reference gust w(Ω) measured at the nose, each lifting element i
  // sees  w_i = sinc(Ω·c_i/2)·e^{−iΩx_i}·w  (chord averaging + transport delay);
  // the hull sees the hull-length-averaged gust (Munk moment).
  //
  // Underwater the accelerating water also pushes the hull (Froude–Krylov +
  // added mass, Morison): F = (1 + k2)·ρ·Vol·∂w_hull/∂t, acting at mid-hull.
  // Negligible in air, first-order in water – it makes the body drift with
  // the water parcel. Dynamic term only; the static L = W balance is untouched.
  //
  // Downwash: a fin with dw = {src, E, tau} feels  −E·e^{−iωτ}·(incidence of fin src),
  // τ = l_H/V (wake transport lag). Eigenvalues use the quasi-steady limit (τ → 0, no
  // downwash-lag damping – conservative); the gust transfer functions keep the lag.
  function planeSystem(p) {
    // p: { V, qbar, fins:[{K, l, chord, x, dw?}], Mah, Mqx, mPrime, Iprime, hullL, fkMass, fkArm }
    const V = p.V;
    // Incidence of fin i = aA_i·α + aQ_i·q + (gust part), with complex aA, aQ at temporal frequency w.
    function coeffs(w) {
      let La = 0, Lai = 0, Lq = 0, Lqi = 0, Ma = p.Mah, Mai = 0, Mq = p.Mqx, Mqi = 0;
      const rows = new Array(p.fins.length);
      for (let i = 0; i < p.fins.length; i++) {
        const f = p.fins[i];
        let cr = 0, ci = 0, ls = 0;
        if (f.dw) { cr = f.dw.E * Math.cos(w * f.dw.tau); ci = -f.dw.E * Math.sin(w * f.dw.tau); ls = p.fins[f.dw.src].l; }
        const aAr = 1 - cr, aAi = -ci, aQr = (f.l - cr * ls) / V, aQi = -ci * ls / V;
        rows[i] = { c: cx(cr, ci), aA: cx(aAr, aAi), aQ: cx(aQr, aQi) };
        La += f.K * aAr; Lai += f.K * aAi; Lq += f.K * aQr; Lqi += f.K * aQi;
        Ma -= f.K * f.l * aAr; Mai -= f.K * f.l * aAi; Mq -= f.K * f.l * aQr; Mqi -= f.K * f.l * aQi;
      }
      return { La: cx(La, Lai), Lq: cx(Lq, Lqi), Ma: cx(Ma, Mai), Mq: cx(Mq, Mqi), rows };
    }
    const s0 = coeffs(0);
    const hasDw = p.fins.some(f => f.dw);   // without downwash the coefficients are frequency-independent
    const La = s0.La.re, Lq = s0.Lq.re, Ma = s0.Ma.re, Mq = s0.Mq.re;
    const mV = p.mPrime * V;
    // A = [[-La/mV, 1 - Lq/mV], [Ma/I, Mq/I]]
    const a11 = -La / mV, a12 = 1 - Lq / mV, a21 = Ma / p.Iprime, a22 = Mq / p.Iprime;
    const tr = a11 + a22, det = a11 * a22 - a12 * a21;
    const disc = tr * tr / 4 - det;
    let eig;
    if (disc >= 0) eig = [cx(tr / 2 + Math.sqrt(disc)), cx(tr / 2 - Math.sqrt(disc))];
    else eig = [cx(tr / 2, Math.sqrt(-disc)), cx(tr / 2, -Math.sqrt(-disc))];
    const stable = tr < 0 && det > 0 && La > 0;
    const wn = det > 0 ? Math.sqrt(det) : null;
    const zeta = wn ? -tr / (2 * wn) : null;

    // Transfer functions at spatial frequency Ω (rigid = vehicle attitude frozen).
    function transfer(Om, rigid) {
      const w = Om * V;                     // temporal frequency
      const sc = hasDw ? coeffs(w) : s0;
      const fw = p.fins.map(f => cscale(cexpi(-Om * f.x), sinc(Om * f.chord / 2) / V)); // w_i/V per unit w
      const gw = p.fins.map((f, i) => f.dw ? cadd(fw[i], cscale(cmul(sc.rows[i].c, fw[f.dw.src]), -1)) : fw[i]);
      const hw = cscale(cexpi(-Om * p.hullL / 2), sinc(Om * p.hullL / 2) / V);
      const Ffk = cmul(cx(0, w), cscale(hw, p.fkMass * V));          // (1+k2)ρVol · iω · w_hull
      let Lw = Ffk, Mw = cadd(cscale(hw, p.Mah), cscale(Ffk, -p.fkArm));
      p.fins.forEach((f, i) => { Lw = cadd(Lw, cscale(gw[i], f.K)); Mw = cadd(Mw, cscale(gw[i], -f.K * f.l)); });
      let alpha = cx(0), q = cx(0);
      if (!rigid) {
        // (iω + La/mV)α + (Lq/mV − 1)q = −Lw/mV ;  −(Ma/I)α + (iω − Mq/I)q = Mw/I
        const A11 = cadd(cx(0, w), cscale(sc.La, 1 / mV)), A12 = cadd(cscale(sc.Lq, 1 / mV), cx(-1));
        const A21 = cscale(sc.Ma, -1 / p.Iprime), A22 = cadd(cx(0, w), cscale(sc.Mq, -1 / p.Iprime));
        const b1 = cscale(Lw, -1 / mV), b2 = cscale(Mw, 1 / p.Iprime);
        const D = cadd(cmul(A11, A22), cscale(cmul(A12, A21), -1));
        alpha = cdiv(cadd(cmul(b1, A22), cscale(cmul(A12, b2), -1)), D);
        q = cdiv(cadd(cmul(A11, b2), cscale(cmul(A21, b1), -1)), D);
      }
      // Effective incidence at each fin, lift perturbation, path-rate.
      const finAlpha = p.fins.map((f, i) => f.dw
        ? cadd(cadd(cmul(sc.rows[i].aA, alpha), cmul(sc.rows[i].aQ, q)), gw[i])
        : cadd(cadd(alpha, cscale(q, f.l / V)), gw[i]));
      let dL = cx(0);
      p.fins.forEach((f, i) => { dL = cadd(dL, cscale(finAlpha[i], f.K)); });
      const gammaDot = cscale(cadd(dL, rigid ? cx(0) : Ffk), 1 / mV);
      return { finAlpha, dL, q, alpha, gammaDot };
    }
    return { La, Lq, Ma, Mq, eig, stable, wn, zeta, transfer };
  }

  // Integrate response statistics of one plane over the gust spectrum.
  function planeStats(sys, phi, grid, rigid) {
    const nF = sys.transfer(1, true).finAlpha.length;
    const varFin = new Array(nF).fill(0);
    let varL = 0, varQ = 0, varAlpha = 0;
    for (let k = 0; k < grid.n; k++) {
      const Om = grid.om[k], ph = phi(Om) * grid.w[k];
      const t = sys.transfer(Om, rigid);
      for (let i = 0; i < nF; i++) varFin[i] += cabs2(t.finAlpha[i]) * ph;
      varL += cabs2(t.dL) * ph;
      varQ += cabs2(t.q) * ph;
      varAlpha += cabs2(t.alpha) * ph;
    }
    return { sigFinAlpha: varFin.map(Math.sqrt), sigL: Math.sqrt(varL), sigQ: Math.sqrt(varQ), sigAlpha: Math.sqrt(varAlpha) };
  }

  // ---------- Lateral-directional dynamics (sideslip, roll, yaw, bank) ----------
  // Small helpers: characteristic polynomial (Faddeev–LeVerrier) and its roots (Durand–Kerner).
  function charPoly(A) {
    const n = A.length, c = new Array(n + 1).fill(0);
    c[n] = 1;
    let M = A.map(r => r.map(() => 0));
    const mul = (X, Y) => X.map((r, i) => Y[0].map((_, j) => r.reduce((a, _x, k) => a + X[i][k] * Y[k][j], 0)));
    for (let k = 1; k <= n; k++) {
      M = mul(A, M); for (let i = 0; i < n; i++) M[i][i] += c[n - k + 1];
      const AM = mul(A, M);
      c[n - k] = -AM.reduce((a, r, i) => a + r[i], 0) / k;
    }
    return c;                       // c[i] = coefficient of λ^i
  }
  function polyRoots(c) {
    const n = c.length - 1;
    const scale = Math.max(1, ...c.slice(0, n).map(Math.abs));
    let z = Array.from({ length: n }, (_, i) => cscale(cexpi(0.4 + 2 * Math.PI * i / n), scale));
    const ev = (x) => { let r = cx(1); for (let i = n - 1; i >= 0; i--) r = cadd(cmul(r, x), cx(c[i])); return r; };
    for (let it = 0; it < 2000; it++) {
      let moved = 0;
      z = z.map((zi, i) => {
        let den = cx(1);
        z.forEach((zj, j) => { if (j !== i) den = cmul(den, cadd(zi, cscale(zj, -1))); });
        const d = cdiv(ev(zi), den);
        moved = Math.max(moved, Math.sqrt(cabs2(d)) / (1 + Math.sqrt(cabs2(zi))));
        return cadd(zi, cscale(d, -1));
      });
      if (moved < 1e-13) break;
    }
    return z.map(e => Math.abs(e.im) < 1e-9 * (1 + Math.abs(e.re)) ? cx(e.re) : e);
  }
  // Complex Gaussian elimination: solves M·x = b (M n×n of {re,im}).
  function csolve(M, b) {
    const n = b.length, A = M.map((r, i) => [...r, b[i]]);
    for (let k = 0; k < n; k++) {
      let piv = k;
      for (let i = k + 1; i < n; i++) if (cabs2(A[i][k]) > cabs2(A[piv][k])) piv = i;
      [A[k], A[piv]] = [A[piv], A[k]];
      for (let i = k + 1; i < n; i++) {
        const f = cdiv(A[i][k], A[k][k]);
        for (let j = k; j <= n; j++) A[i][j] = cadd(A[i][j], cscale(cmul(f, A[k][j]), -1));
      }
    }
    const x = new Array(n);
    for (let i = n - 1; i >= 0; i--) {
      let sum = A[i][n];
      for (let j = i + 1; j < n; j++) sum = cadd(sum, cscale(cmul(A[i][j], x[j]), -1));
      x[i] = cdiv(sum, A[i][i]);
    }
    return x;
  }

  // Strip-theory roll/yaw derivatives of a horizontal (or panel) pair, dimensional:
  //   C_lβ = −a·Γ·(1+2λ)/(6(1+λ)) − CL·tanΛ¼c·(1+2λ)/(3(1+λ))     (dihedral + sweep)
  //   L_p  = −q·a·c_r·b³·(1+3λ)/(48V)                              (C_lp = −a(1+3λ)/(12(1+λ)))
  //   C_lr = CL/4,  C_np = −CL/8,  C_nr = −CD0/4   (per p·b/2V, r·b/2V; Nelson)
  function pairRollDerivs(pair, slope, CL, dihedralRad, qbar, V) {
    const f = pair.pf, b = f.b, lam = f.taper, qSb = qbar * pair.area * b, k = b / (2 * V);
    const Clb = -slope * dihedralRad * (1 + 2 * lam) / (6 * (1 + lam)) - CL * f.tanQc * (1 + 2 * lam) / (3 * (1 + lam));
    return {
      Lb: qSb * Clb,
      Lp: -qbar * slope * f.cr * b * b * b * (1 + 3 * lam) / (48 * V),
      Lr: qSb * CL / 4 * k, Np: -qSb * CL / 8 * k, Nr: -qSb * pair.cd0 / 4 * k
    };
  }

  // States [β, p, r, φ] (sideslip, roll rate, yaw rate, bank; stability axes, I_xz neglected):
  //   m'V(β̇ + r) = Y + W·φ,   I_x·ṗ = L − K_φ·φ,   I_z'·ṙ = N,   φ̇ = p
  // Side-force surfaces i (fin, V-tail yaw share, cruciform pairs) at arm l_i (aft) and height
  // z_i above the roll axis feel β_i = β − r·l_i/V + p·z_i/V (+ gust) and give Y_i = −K_i·β_i,
  // L_i = Y_i·z_i, N_i = −Y_i·l_i. Horizontal pairs add the strip-theory derivatives above.
  // Inputs: lateral gust v (delayed/chord-averaged per surface, Froude–Krylov on the hull) and
  // the rolling gust p_g (antisymmetric w, acts on the wing like a roll rate).
  function lateralSystem(p) {
    const V = p.V, mV = p.mPrime * V, Ix = p.Ix, Iz = p.Iz;
    const H = p.horiz;   // summed horizontal-pair derivatives
    let Yb = 0, Yp = 0, Yr = 0, Lb = H.Lb, Lp = H.Lp, Lr = H.Lr, Nb = -p.Mah, Np = H.Np, Nr = H.Nr;
    p.fins.forEach(f => {
      Yb -= f.K; Yp -= f.K * f.z / V; Yr += f.K * f.l / V;
      Lb -= f.K * f.z; Lp -= f.K * f.z * f.z / V; Lr += f.K * f.l * f.z / V;
      Nb += f.K * f.l; Np += f.K * f.l * f.z / V; Nr -= f.K * f.l * f.l / V;
    });
    const Afull = [
      [Yb / mV, Yp / mV, Yr / mV - 1, p.W / mV],
      [Lb / Ix, Lp / Ix, Lr / Ix, -p.Kphi / Ix],
      [Nb / Iz, Np / Iz, Nr / Iz, 0],
      [0, 1, 0, 0]
    ];
    // Roll held level (p = φ = 0, e.g. AUV roll control): only sideslip and yaw remain.
    const idx = p.rollLocked ? [0, 2] : [0, 1, 2, 3];
    const A = idx.map(i => idx.map(j => Afull[i][j]));
    const eig = polyRoots(charPoly(A)).sort((a, b) => a.re - b.re);
    const stable = eig.every(e => e.re < -1e-9);
    // Mode identification: complex pair = Dutch roll (the faster one if two pairs);
    // real roots: most negative = roll subsidence, smallest |λ| = spiral.
    const cplx = eig.filter(e => e.im > 0), real = eig.filter(e => e.im === 0);
    const wnOf = (e) => Math.sqrt(e.re * e.re + e.im * e.im);
    const dr = cplx.length ? cplx.reduce((a, e) => wnOf(e) > wnOf(a) ? e : a) : null;
    const dutch = dr ? { wn: wnOf(dr), zeta: -dr.re / wnOf(dr), eig: dr } : null;
    const others = eig.filter(e => !(dr && Math.abs(e.re - dr.re) < 1e-12 && Math.abs(Math.abs(e.im) - dr.im) < 1e-12));
    let roll = null, spiral = null, coupled = null;
    if (p.rollLocked) { /* no roll or spiral mode */ }
    else if (others.every(e => e.im === 0) && others.length) {
      const srt = others.slice().sort((a, b) => Math.abs(a.re) - Math.abs(b.re));
      spiral = { lambda: srt[0].re };
      if (srt.length > 1) roll = { lambda: srt[srt.length - 1].re, tau: -1 / srt[srt.length - 1].re };
    } else if (others.length) {
      const e = others.find(x => x.im > 0);
      coupled = { wn: wnOf(e), zeta: -e.re / wnOf(e) };   // roll–spiral ("lateral phugoid")
    }
    if (spiral) spiral.T = spiral.lambda !== 0 ? Math.LN2 / Math.abs(spiral.lambda) : Infinity;   // time to double (unstable) or halve

    function transfer(Om, rigid, input) {
      const w = Om * V;
      const isV = input !== "p";
      const fw = p.fins.map(f => isV ? cscale(cexpi(-Om * f.x), sinc(Om * f.chord / 2) / V) : cx(0));
      const hw = isV ? cscale(cexpi(-Om * p.hullL / 2), sinc(Om * p.hullL / 2) / V) : cx(0);
      const ww = isV ? cscale(cexpi(-Om * p.wingX), sinc(Om * p.wingChord / 2) / V) : cx(0);
      const Ffk = cmul(cx(0, w), cscale(hw, -p.fkMass * V));
      let Yg = Ffk, Lg = isV ? cscale(ww, H.Lb) : cx(-H.Lp), Ng = isV ? cscale(hw, -p.Mah) : cx(-H.Np);
      Ng = cadd(Ng, cscale(Ffk, -p.fkArm));
      p.fins.forEach((f, i) => {
        Yg = cadd(Yg, cscale(fw[i], -f.K));
        Lg = cadd(Lg, cscale(fw[i], -f.K * f.z));
        Ng = cadd(Ng, cscale(fw[i], f.K * f.l));
      });
      const x = [cx(0), cx(0), cx(0), cx(0)];
      if (!rigid) {
        const rhs = [cscale(Yg, 1 / mV), cscale(Lg, 1 / Ix), cscale(Ng, 1 / Iz), cx(0)];
        const M = A.map((r, i) => r.map((a, j) => cx(-a, i === j ? w : 0)));   // iωI − A
        const xs = csolve(M, idx.map(i => rhs[i]));
        idx.forEach((i, k) => { x[i] = xs[k]; });
      }
      const [beta, pr, r, phi] = x;
      const finAlpha = p.fins.map((f, i) => cadd(cadd(cadd(beta, cscale(r, -f.l / V)), cscale(pr, f.z / V)), fw[i]));
      let dY = cx(0);
      p.fins.forEach((f, i) => { dY = cadd(dY, cscale(finAlpha[i], -f.K)); });
      const gammaDot = cscale(cadd(cadd(dY, rigid ? cx(0) : Ffk), cscale(phi, p.W)), 1 / mV);   // course rate
      return { finAlpha, dL: dY, q: r, alpha: beta, gammaDot, phi, p: pr };
    }
    const La = -Yb;   // side-force slope (kept for compatibility with the pitch-plane object)
    return {
      La, eig, stable, dutch, roll, spiral, coupled, transfer,
      wn: dutch ? dutch.wn : null, zeta: dutch ? dutch.zeta : null,
      deriv: { Yb, Yp, Yr, Lb, Lp, Lr, Nb, Np, Nr },
      // Views with one gust input each (same interface as planeSystem for the integrators).
      byInput: (inp) => ({ stable, transfer: (Om, rigid) => transfer(Om, rigid, inp) })
    };
  }

  // RMS statistics of the lateral system over both gust inputs (v and rolling gust p_g).
  function lateralStats(sys, phiV, phiP, grid, rigid) {
    const nF = sys.transfer(1, true, "v").finAlpha.length;
    const vF = new Array(nF).fill(0);
    let vY = 0, vR = 0, vB = 0, vPhi = 0, vP = 0;
    [["v", phiV], ["p", phiP]].forEach(([inp, phi]) => {
      for (let k = 0; k < grid.n; k++) {
        const Om = grid.om[k], ph = phi(Om) * grid.w[k];
        const t = sys.transfer(Om, rigid, inp);
        for (let i = 0; i < nF; i++) vF[i] += cabs2(t.finAlpha[i]) * ph;
        vY += cabs2(t.dL) * ph; vR += cabs2(t.q) * ph; vB += cabs2(t.alpha) * ph;
        vPhi += cabs2(t.phi) * ph; vP += cabs2(t.p) * ph;
      }
    });
    return { sigFinAlpha: vF.map(Math.sqrt), sigL: Math.sqrt(vY), sigQ: Math.sqrt(vR), sigAlpha: Math.sqrt(vB), sigPhi: Math.sqrt(vPhi), sigP: Math.sqrt(vP) };
  }

  // Rolling-gust spectrum (MIL-F-8785C), spatial: Φ_p(Ω) = σ_w²/L_w · 0.8(πL_w/4b)^{1/3} / (1 + (4bΩ/π)²).
  function phiRolling(Om, sigma, L, b) {
    const x = 4 * b * Om / Math.PI;
    return sigma * sigma / L * 0.8 * Math.cbrt(Math.PI * L / (4 * b)) / (1 + x * x);
  }

  // RMS track wander about the mean course after distance X, no steering.
  // y(T) − y(0) = V ∫₀ᵀ γ dt, γ = γ̇/(iω) (stationary path angle; linear theory
  // returns the path to its original direction once a gust has passed):
  //   Var = V² ∫ |G_γ̇|² Φ · 4·sin²(ωT/2)/ω⁴ dΩ
  // (Forcing γ(0) = 0 instead would turn the launch-instant path angle into a
  //  permanent course offset – reported separately as σ_γ.)
  function dispersion(sys, phi, V, X, cfg) {
    if (!sys.stable) return Infinity;
    const T = X / V;
    const g = logGrid(Math.min(cfg.omegaMin_radpm, 0.01 / X), cfg.omegaMax_radpm, cfg.points);
    let v = 0;
    for (let k = 0; k < g.n; k++) {
      const Om = g.om[k], w = Om * V, s = Math.sin(w * T / 2);
      v += cabs2(sys.transfer(Om, false).gammaDot) * phi(Om) * 4 * s * s / (w * w * w * w) * g.w[k];
    }
    return V * Math.sqrt(v);
  }

  // Same, measured from the INITIAL course (γ(0) = 0):
  //   y(T) = V ∫₀ᵀ ∫₀ᵗ γ̇ dt' dt  →  Var = V² ∫ |G_γ̇|² Φ · |(e^{iωT} − 1 − iωT)/(iω)|²/ω² dΩ
  // Finite even when a steady input gives a steady turn rate (heading random walk), as the
  // rolling gust does; then the "about the mean course" form above has no finite limit.
  function dispersionFromStart(sys, phi, V, X, cfg) {
    if (!sys.stable) return Infinity;
    const T = X / V;
    const g = logGrid(Math.min(cfg.omegaMin_radpm, 0.01 / X), cfg.omegaMax_radpm, cfg.points);
    let v = 0;
    for (let k = 0; k < g.n; k++) {
      const Om = g.om[k], w = Om * V, a = w * T;
      // |(e^{ia} − 1 − ia)/(ia)|²·T²/ω²  (series for small a: T⁴/4)
      let kern;
      if (a < 1e-3) kern = T * T * T * T / 4 * (1 - a * a / 9);
      else { const re = Math.cos(a) - 1, im = Math.sin(a) - a; kern = (re * re + im * im) / (a * a) * T * T / (w * w); }
      v += cabs2(sys.transfer(Om, false).gammaDot) * phi(Om) * kern * g.w[k];
    }
    return V * Math.sqrt(v);
  }

  // RMS path-angle wander σ_γ (rad).
  function pathAngleRms(sys, phi, V, cfg) {
    if (!sys.stable) return Infinity;
    const g = logGrid(cfg.omegaMin_radpm, cfg.omegaMax_radpm, cfg.points);
    let v = 0;
    for (let k = 0; k < g.n; k++) {
      const Om = g.om[k], w = Om * V;
      v += cabs2(sys.transfer(Om, false).gammaDot) * phi(Om) / (w * w) * g.w[k];
    }
    return Math.sqrt(v);
  }

  // ---------- Resolve design + environment into one computational context ----------
  function resolve(design, env, model, calib) {
    const ph = model.physics, hy = model.hydro;
    const c = calib || {};
    const g = ph.gravity_mps2;
    const rho = env.density_kgpm3;
    design = normalizeDesign(design, model);
    const hull = hullGeometry(design);
    // Section data evaluated at the fin chord Re for the cruise speed (fixed over the speed scan).
    const secs = (model.sections && model.sections.sections) || [];
    const secOf = (id) => secs.find(x => x.id === (id || (model.sections && model.sections.default) || "ideal")) || null;
    const front = finPair("front", design, c.clAlphaFront, secOf(design.frontSection), env, model);
    const tailType = design.tailType;
    const rear = tailType === "tailless" ? null : finPair("rear", design, c.clAlphaRear, secOf(design.rearSection), env, model);
    const fin = tailType === "traditional" || tailType === "tailless" ? verticalFin(design, secOf(design.finSection), env, model) : null;
    // Share of each surface's force in pitch / yaw. V-tail (Purser–Campbell): cos Γ / sin Γ per
    // panel, so the effective slope in each plane scales with cos²Γ / sin²Γ.
    const G = (design.rearDihedral_deg || 0) / DEG;
    const proj = {
      front: { pitch: 1, yaw: tailType === "cruciform" ? 1 : 0 },
      rear: tailType === "vtail" ? { pitch: Math.cos(G), yaw: Math.sin(G) } : { pitch: 1, yaw: tailType === "cruciform" ? 1 : 0 },
      fin: { pitch: 0, yaw: 1 }
    };
    const m0 = design.emptyMass_kg + (env.payload_kg || 0) + design.fuelMass_kg;
    const mEnd = m0 - design.fuelMass_kg * (1 - env.reserveFraction);
    const W = m0 * g;
    const xcg = design.xcg_m;
    // Pitch inertia: override or uniform solid cylinder; added mass/inertia from Lamb.
    const I0 = design.pitchInertia_kgm2 != null && design.pitchInertia_kgm2 > 0
      ? design.pitchInertia_kgm2 : m0 * (hull.L * hull.L / 12 + hull.D * hull.D / 16);
    // Roll / yaw inertia: override or estimate. Wing part from the non-dimensional radius of
    // gyration R̄x = 2k_x/b (model.json lateral.rollGyrationFactor, ≈ 0.25 for light aircraft).
    const Rx = (model.lateral && model.lateral.rollGyrationFactor != null) ? model.lateral.rollGyrationFactor : 0.25;
    const IxWing = m0 * Math.pow(Rx * front.pf.b / 2, 2);
    const Ix = design.rollInertia_kgm2 > 0 ? design.rollInertia_kgm2 : m0 * hull.D * hull.D / 8 + IxWing;
    const Iz0 = design.yawInertia_kgm2 > 0 ? design.yawInertia_kgm2 : I0 + IxWing;
    const a = hull.L / 2, b = hull.D / 2;
    const Iadd = hull.lamb.kp * rho * hull.volume * (a * a + b * b) / 5;
    const mAdd = hull.lamb.k2 * rho * hull.volume;
    const An = Math.PI / 4 * design.nozzleDiameter_m * design.nozzleDiameter_m;
    const etaJet = c.etaJet != null ? c.etaJet : design.jetEfficiency;
    return {
      g, rho, hull, front, rear, fin, tailType, proj, m0, mEnd, W, xcg, I0, Iadd, mAdd, An, etaJet,
      a: env.speedOfSound_mps || 340.29, atm: env.atm || null,
      // Available engine power ∝ σ^n with altitude (ISA only; model.json engine.powerLapseExponent).
      powerFactor: env.atm ? Math.pow(env.atm.sigma, (model.engine && model.engine.powerLapseExponent != null) ? model.engine.powerLapseExponent : 0.7) : 1,
      Ix, Iz: Iz0 + Iadd, Kphi: design.extraRollStiffness_Nmprad || 0, rollFree: design.rollFree !== false, wingDihedral: (design.frontDihedral_deg || 0) / DEG,
      vtailDihedral: tailType === "vtail" ? (design.rearDihedral_deg || 0) / DEG : 0,
      nu: env.kinematicViscosity_m2ps,
      extraDragArea: (c.extraDragArea_m2 != null ? c.extraDragArea_m2 : (hy.extraDragArea_m2 || 0)),
      hullCDwetOverride: c.hullCDwet != null ? c.hullCDwet : null,
      cruciform: tailType === "cruciform",
      downwash: design.downwash !== false && !!rear, tailHeight: design.rearHeight_m || 0,
      extraDamping: design.extraPitchDamping_Nms || 0
    };
  }

  // Lift split between the two pairs (moment balance about CG, L_f + L_r = W), lift acting at xf / xr.
  function liftSplit(ctx, xf, xr) {
    const dx = xr - xf;
    if (Math.abs(dx) < 1e-6) return { Lf: ctx.W / 2, Lr: ctx.W / 2, degenerate: true };
    const Lf = ctx.W * (xr - ctx.xcg) / dx;
    return { Lf, Lr: ctx.W - Lf, degenerate: false };
  }

  // Trim at dynamic pressure qbar. A delta's lift centre moves with α (vortex lift), so the
  // split and the lift states are iterated to a fixed point.
  // Rear CL is the panel CL (V-tail panels carry L_rear/cos Γ). Tailless: the wing carries W
  // and elevon trim is assumed (its moment balances; its trim drag is not modelled).
  function trim(ctx, qbar) {
    if (!ctx.rear) {
      const split = { Lf: ctx.W, Lr: 0, degenerate: false, elevonTrim: true };
      const CLf = ctx.W / (qbar * ctx.front.area);
      const fs = liftState(ctx.front, CLf);
      return { split, CLf, CLr: 0, fs, rs: null, dw: { E: 0, lH: 0 } };
    }
    const pr = ctx.proj.rear.pitch;
    let xf = ctx.front.pf.xAc25, xr = ctx.rear.pf.xAc25, split, CLf, CLr, fs, rs;
    for (let it = 0; it < 8; it++) {
      split = liftSplit(ctx, xf, xr);
      CLf = split.Lf / (qbar * ctx.front.area);
      CLr = split.Lr / (qbar * ctx.rear.area * pr);
      fs = liftState(ctx.front, CLf); rs = liftState(ctx.rear, CLr);
      if (Math.abs(fs.xLift - xf) < 1e-7 && Math.abs(rs.xLift - xr) < 1e-7) break;
      xf = fs.xLift; xr = rs.xLift;
    }
    return { split, CLf, CLr, fs, rs, dw: downwashGradient(ctx, fs, rs) };
  }

  // Pitch/heave (vertical) and yaw/sway (lateral) systems at speed V.
  // vTags / lTags describe the fins of each system: { name, pair, proj, st } (st = lift state about trim).
  function planeSystems(ctx, V, qbar, tr) {
    const Mah = (ctx.hull.lamb.k2 - ctx.hull.lamb.k1) * ctx.rho * ctx.hull.volume * V * V;
    const base = {
      V, qbar, Mah, Mqx: -Math.abs(ctx.extraDamping), mPrime: ctx.m0 + ctx.mAdd, Iprime: ctx.I0 + ctx.Iadd, hullL: ctx.hull.L,
      fkMass: (1 + ctx.hull.lamb.k2) * ctx.rho * ctx.hull.volume, fkArm: ctx.hull.L / 2 - ctx.xcg
    };
    const fin = (t, slope, x, dw) => ({ K: qbar * t.pair.area * slope * t.proj * t.proj, l: x - ctx.xcg, chord: t.pair.chord, x, dw });
    const vTags = [{ name: "front", pair: ctx.front, proj: 1, st: tr.fs }];
    if (ctx.rear && ctx.proj.rear.pitch > 0) vTags.push({ name: "rear", pair: ctx.rear, proj: ctx.proj.rear.pitch, st: tr.rs });
    const dw = tr.dw.E > 0 ? { src: 0, E: tr.dw.E, tau: tr.dw.lH / V } : null;
    const vert = planeSystem({ ...base, fins: vTags.map((t, i) => fin(t, t.st.slope, t.st.xAc, i === 1 ? dw : null)) });
    // Lateral surfaces fly at zero side-force trim: attached-flow slope at ¼-MAC, no sidewash.
    const lTags = [];
    if (ctx.proj.front.yaw > 0) lTags.push({ name: "front", pair: ctx.front, proj: ctx.proj.front.yaw });
    if (ctx.rear && ctx.proj.rear.yaw > 0) lTags.push({ name: ctx.tailType === "vtail" ? "V-tail" : "rear", pair: ctx.rear, proj: ctx.proj.rear.yaw });
    if (ctx.fin) lTags.push({ name: "fin", pair: ctx.fin, proj: 1 });
    lTags.forEach(t => { t.st = liftState(t.pair, 0); });
    // Height of each side-force surface's a.c. above the roll axis (hull centreline).
    const zOf = (t) => t.pair === ctx.fin ? ctx.hull.D / 2 + ctx.fin.pf.yMac
      : (ctx.tailType === "vtail" && t.pair === ctx.rear ? (ctx.hull.D / 2 + ctx.rear.pf.yMac) * Math.sin(ctx.vtailDihedral) : 0);
    // Horizontal pairs: wing (with dihedral), rear tail / V-tail panels (roll damping), and in
    // cruciform the vertical pairs as well (roll damping only).
    const H = { Lb: 0, Lp: 0, Lr: 0, Np: 0, Nr: 0 };
    const addH = (d) => Object.keys(H).forEach(k => { H[k] += d[k]; });
    addH(pairRollDerivs(ctx.front, tr.fs.slope, tr.CLf, ctx.wingDihedral, qbar, V));
    if (ctx.rear) addH(pairRollDerivs(ctx.rear, tr.rs.slope, ctx.tailType === "vtail" ? 0 : tr.CLr, 0, qbar, V));
    if (ctx.cruciform) {
      addH(pairRollDerivs(ctx.front, ctx.front.clAlpha, 0, 0, qbar, V));
      addH(pairRollDerivs(ctx.rear, ctx.rear.clAlpha, 0, 0, qbar, V));
    }
    const lat = lateralSystem({
      V, W: ctx.W, mPrime: ctx.m0 + ctx.mAdd, Ix: ctx.Ix, Iz: ctx.Iz, Kphi: ctx.Kphi, rollLocked: !ctx.rollFree,
      Mah, fkMass: base.fkMass, fkArm: base.fkArm, hullL: ctx.hull.L,
      wingX: tr.fs.xAc, wingChord: ctx.front.chord, horiz: H,
      fins: lTags.map(t => ({ ...fin(t, t.pair.clAlpha, t.pair.pf.xAc25), z: zOf(t) }))
    });
    return { vert, lat, vTags, lTags };
  }

  // Turbulent skin friction compressibility factor (Raymer): Cf ∝ (1 + 0.144·M²)^−0.65.
  const compressibleCf = (M) => Math.pow(1 + 0.144 * M * M, -0.65);

  // Wave drag of all lifting surfaces (CD·S summed), per surface Korn/Lock data.
  function waveTerms(ctx, tr, M) {
    const list = [{ name: "front", pair: ctx.front, CL: tr.CLf }];
    if (ctx.rear) list.push({ name: ctx.tailType === "vtail" ? "V-tail" : "rear", pair: ctx.rear, CL: tr.CLr });
    if (ctx.cruciform) { list.push({ name: "front (vert.)", pair: ctx.front, CL: 0 }); if (ctx.rear) list.push({ name: "rear (vert.)", pair: ctx.rear, CL: 0 }); }
    if (ctx.fin) list.push({ name: "fin", pair: ctx.fin, CL: 0 });
    const surfaces = list.map(t => ({ ...t, ...waveDrag(t.pair, t.CL, M) }));
    return { surfaces, CDA: surfaces.reduce((a, t) => a + t.CD * t.pair.area, 0) };
  }

  // Calm-air drag of steady level flight (no turbulence) – used for the altitude envelope.
  function calmDrag(ctx, V, tr) {
    const qbar = 0.5 * ctx.rho * V * V, M = V / ctx.a;
    const Cf = ittcCf(V * ctx.hull.L / ctx.nu) * compressibleCf(M);
    const CDwet = ctx.hullCDwetOverride != null ? ctx.hullCDwetOverride : Cf * ctx.hull.formFactor;
    const prof = (ctx.front.cd0 * ctx.front.area + (ctx.rear ? ctx.rear.cd0 * ctx.rear.area : 0)) * (ctx.cruciform ? 2 : 1) + (ctx.fin ? ctx.fin.cd0 * ctx.fin.area : 0);
    const induced = ctx.front.area * tr.fs.CDi + (ctx.rear ? ctx.rear.area * tr.rs.CDi : 0);
    return qbar * (ctx.hull.wetted * CDwet + prof + induced + ctx.extraDragArea + waveTerms(ctx, tr, M).CDA);
  }

  // Flight envelope (ISA, calm air, take-off mass): for each speed, the band of altitudes where
  // level flight is possible – L = W without exceeding CLmax, engine power required ≤ available
  // (P_max·σ^n) and M ≤ MACH_CAP. Also the stall-limited and power-limited ceilings separately.
  function altitudeEnvelope(design, env, model, calib, opts) {
    if (env.useISA === false) return null;
    const am = model.atmosphere || {};
    const hTop = am.envelopeMaxAltitude_m || 20000, nH = am.envelopeAltitudeSteps || 40, nV = am.envelopeSpeedSteps || 40;
    const pr = model.presets, v0 = pr.speedScanMin_mps, v1 = pr.speedScanMax_mps;
    const check = (V, h) => {
      const e = resolveAtmosphere({ ...env, altitude_m: h, cruiseSpeed_mps: V });
      const c = resolve(design, e, model, calib);
      const tr = trim(c, 0.5 * c.rho * V * V);
      const stall = !!(tr.fs.stalled || (tr.rs && tr.rs.stalled));
      const P = jetForThrust(calmDrag(c, V, tr), V, c.rho, c.An).Pjet / c.etaJet;
      const power = P <= opts.maxPower_W * c.powerFactor && V / c.a <= MACH_CAP;
      return { stall, power, ok: !stall && power };
    };
    // Boundary between a passing altitude `pass` and a failing altitude `fail` (bisection, dh/64).
    const edge = (V, pass, fail, test) => {
      for (let k = 0; k < 6; k++) { const m = (pass + fail) / 2; if (test(check(V, m))) pass = m; else fail = m; }
      return pass;
    };
    const dh = hTop / nH, points = [];
    for (let i = 0; i <= nV; i++) {
      const V = v0 + (v1 - v0) * i / nV;
      let hMin = null, hMax = null, hStall = null, hPower = null;
      const res = [];
      for (let j = 0; j <= nH; j++) res.push(check(V, dh * j));
      res.forEach((c, j) => {
        if (!c.stall) hStall = dh * j;
        if (c.power) hPower = dh * j;
        if (c.ok) { if (hMin == null) hMin = dh * j; hMax = dh * j; }
      });
      // Refine the coarse grid where a boundary lies inside the range.
      if (hMax != null && hMax < hTop) hMax = edge(V, hMax, hMax + dh, c => c.ok);
      if (hMin != null && hMin > 0) hMin = edge(V, hMin, hMin - dh, c => c.ok);
      if (hStall != null && hStall < hTop) hStall = edge(V, hStall, hStall + dh, c => !c.stall);
      if (hPower != null && hPower < hTop) hPower = edge(V, hPower, hPower + dh, c => c.power);
      points.push({ V, hMin, hMax, hStall, hPower });
    }
    return { points, hTop };
  }

  // Everything at one speed. opts.withDynamics: use vehicle response in drag/load penalty.
  function pointAt(V, ctx, env, model, opts) {
    const sp = model.spectral;
    const grid = logGrid(sp.omegaMin_radpm, sp.omegaMax_radpm, sp.points);
    const qbar = 0.5 * ctx.rho * V * V;
    const tr = trim(ctx, qbar);
    const split = tr.split, CLf = tr.CLf, CLr = tr.CLr;
    const { vert, lat, vTags, lTags } = planeSystems(ctx, V, qbar, tr);
    const hasYaw = lTags.length > 0;
    const tu = env.turbulence;
    const phiW = (Om) => phiTransverse(Om, tu.sigma_w, tu.L_w);
    const phiV = (Om) => phiTransverse(Om, tu.sigma_v, tu.L_v);
    const phiU = (Om) => phiLongitudinal(Om, tu.sigma_u, tu.L_u);

    const useDynV = opts.withDynamics && vert.stable;
    const useDynL = opts.withDynamics && lat.stable;
    const sv = planeStats(vert, phiW, grid, !useDynV);
    const phiP = (Om) => phiRolling(Om, tu.sigma_w, tu.L_w, ctx.front.pf.b);
    const sl = lateralStats(lat, phiV, phiP, grid, !useDynL);
    // Per-surface panel incidence σ (plane incidence × projection).
    const surfV = vTags.map((t, i) => ({ ...t, CL: i === 0 ? CLf : CLr, sig: sv.sigFinAlpha[i] * t.proj }));
    const surfL = lTags.map((t, i) => ({ ...t, sig: sl.sigFinAlpha[i] * t.proj }));

    // Dynamic-pressure penalty from u/v/w (optional; the reference script ignores it).
    let qPar = qbar;
    if (env.includeQPenalty) {
      let vu = 0, vv = 0, vw = 0;
      for (let k = 0; k < grid.n; k++) {
        const Om = grid.om[k], hf = sinc(Om * ctx.hull.L / 2) ** 2 * grid.w[k];
        vu += phiU(Om) * hf; vv += phiV(Om) * hf; vw += phiW(Om) * hf;
      }
      qPar = 0.5 * ctx.rho * (V * V + vu + vv + vw);
    }

    // Drag build-up
    const Re = V * ctx.hull.L / ctx.nu;
    const Mach = V / ctx.a;
    const Cf = ittcCf(Re) * compressibleCf(Mach);
    const CDwet = ctx.hullCDwetOverride != null ? ctx.hullCDwetOverride : Cf * ctx.hull.formFactor;
    const D_hull = qPar * ctx.hull.wetted * CDwet;
    const pairsCD0A = (ctx.front.cd0 * ctx.front.area + (ctx.rear ? ctx.rear.cd0 * ctx.rear.area : 0)) * (ctx.cruciform ? 2 : 1);
    const finCD0A = pairsCD0A + (ctx.fin ? ctx.fin.cd0 * ctx.fin.area : 0);
    const D_finProfile = qPar * finCD0A;
    // Lift-dependent drag: E[CDi(α0 + α_g)] ≈ CDi(α0) + ½·CDi''(α0)·σα²  (exact for the parabolic polar).
    const D_trim = qbar * (ctx.front.area * tr.fs.CDi + (ctx.rear ? ctx.rear.area * tr.rs.CDi : 0));
    const gust = (t) => qbar * t.pair.area * 0.5 * t.st.d2 * t.sig * t.sig;
    const D_gustV = surfV.reduce((a, t) => a + gust(t), 0);
    const D_gustL = surfL.reduce((a, t) => a + gust(t), 0);
    const D_extra = qPar * ctx.extraDragArea;
    const wave = waveTerms(ctx, tr, Mach);
    const D_wave = qbar * wave.CDA;
    const D = D_hull + D_finProfile + D_trim + D_gustV + D_gustL + D_extra + D_wave;
    const D_calm = qbar * ctx.hull.wetted * CDwet + qbar * finCD0A + D_trim + qbar * ctx.extraDragArea + D_wave;

    // Propulsion
    const jet = jetForThrust(D, V, ctx.rho, ctx.An);
    const P_engine = jet.Pjet / ctx.etaJet;
    const etaTot = ctx.etaJet * jet.etaF;
    const Pmax = opts.maxPower_W * ctx.powerFactor;
    const avail = thrustForJetPower(Pmax * ctx.etaJet, V, ctx.rho, ctx.An);
    const throttle = Pmax > 0 ? P_engine / Pmax : Infinity;

    // Breguet range (L = W at m0, as in the script)
    const LD = ctx.W / D;
    const bsfc = opts.bsfc_kgpJ;
    const lnM = Math.log(ctx.m0 / ctx.mEnd);
    const rangeM = etaTot / (bsfc * ctx.g) * LD * lnM;
    const cEq_perHour = bsfc * ctx.g * V / etaTot * 3600;
    const fuelFlow_kgph = bsfc * P_engine * 3600;

    return {
      V, qbar, qPar, Re, Cf, CDwet, split, CLf, CLr, trim: tr, hasYaw, surfV, surfL,
      vert, lat, sv, sl, useDynV, useDynL,
      drag: { hull: D_hull, finProfile: D_finProfile, trim: D_trim, gustV: D_gustV, gustL: D_gustL, extra: D_extra, wave: D_wave, total: D, calm: D_calm },
      Mach, wave, Pavail: Pmax,
      jet, P_engine, etaTot, avail, throttle, feasible: throttle <= 1,
      LD, lnM, rangeKm: rangeM / 1000, enduranceH: rangeM / V / 3600, cEq_perHour, fuelFlow_kgph,
      phiW, phiV, phiP
    };
  }

  // ---------- Main analysis ----------
  function analyze(design, env, model, calib) {
    if (!design || !(design.hullLength_m > 0) || !(design.hullDiameter_m > 0)) return null;
    if (!(design.frontArea_m2 > 0) || !(design.rearArea_m2 > 0) || !(design.nozzleDiameter_m > 0)) return null;
    design = normalizeDesign(design, model);
    env = resolveAtmosphere(env);
    const ctx = resolve(design, env, model, calib);
    const opts = {
      withDynamics: env.useVehicleResponse !== false,
      maxPower_W: design.maxPower_kW * 1000,
      bsfc_kgpJ: (calib && calib.bsfc != null ? calib.bsfc : design.bsfc_kgpkWh) / 3.6e6
    };
    const chosen = pointAt(env.cruiseSpeed_mps, ctx, env, model, opts);

    // Speed sweep: best range & max speed.
    const pr = model.presets;
    const scan = [];
    let best = null, vMax = null;
    for (let V = pr.speedScanMin_mps; V <= pr.speedScanMax_mps + 1e-9; V += pr.speedScanStep_mps) {
      // Section Re and Mach depend on speed → context per speed.
      const p = pointAt(V, resolve(design, { ...env, cruiseSpeed_mps: V }, model, calib), env, model, opts);
      scan.push({ V: round(V, 3), rangeKm: p.rangeKm, feasible: p.feasible, P_kW: p.P_engine / 1000, D: p.drag.total });
      if (p.feasible) {
        vMax = V;
        if (!best || p.rangeKm > best.rangeKm) best = { V, rangeKm: p.rangeKm };
      }
    }

    // Static stability (speed-independent ratios; evaluated at the chosen point).
    const v = chosen.vert;
    const sumK = v.La;
    const xnp = ctx.xcg + (sumK > 0 ? -v.Ma / sumK : -Infinity);
    const staticMargin_m = xnp - ctx.xcg;
    const staticMargin_pct = staticMargin_m / ctx.hull.L * 100;
    const qSb = chosen.qbar * ctx.front.area * ctx.front.pf.b;
    const lateralStatic = { Cnb: chosen.lat.deriv.Nb / qSb, Clb: chosen.lat.deriv.Lb / qSb };

    // Trajectory dispersion with no steering.
    const sp = model.spectral;
    const distances = [];
    const nD = sp.dispersionPoints;
    const xEnd = Math.max(0.1, Math.min(chosen.rangeKm, sp.dispersionMaxKm)) * 1000;
    for (let i = 1; i <= nD; i++) distances.push(xEnd * i / nD);
    const tu = env.turbulence;
    // Lateral wander: lateral gust v and rolling gust p_g are independent → variances add.
    const latWander = (X) => Math.hypot(dispersion(chosen.lat.byInput("v"), chosen.phiV, chosen.V, X, sp), dispersionFromStart(chosen.lat.byInput("p"), chosen.phiP, chosen.V, X, sp));
    const disp = distances.map(X => ({
      km: X / 1000,
      depth_m: dispersion(chosen.vert, chosen.phiW, chosen.V, X, sp),
      lateral_m: latWander(X)
    }));
    const per1km = {
      depth_m: dispersion(chosen.vert, chosen.phiW, chosen.V, 1000, sp),
      lateral_m: latWander(1000)
    };
    const atRange = disp.length ? disp[disp.length - 1] : null;
    const pathAngle = {
      vert_rad: pathAngleRms(chosen.vert, chosen.phiW, chosen.V, sp),
      lat_rad: pathAngleRms(chosen.lat.byInput("v"), chosen.phiV, chosen.V, sp)   // rolling-gust heading drift is a random walk → in the track wander
    };

    // Static-margin design sweep vs the rear pair's root-LE position (tailless: the wing's).
    const smSweep = [];
    const smSweepOf = ctx.rear ? "rear" : "front";
    const nS = 40;
    for (let i = 0; i <= nS; i++) {
      const x = ctx.hull.L * i / nS;
      const c2 = resolve({ ...design, [smSweepOf + "Xle_m"]: x }, env, model, calib);
      const tr2 = trim(c2, chosen.qbar);
      const v2 = planeSystems(c2, chosen.V, chosen.qbar, tr2).vert;
      smSweep.push({ x, sm: v2.La > 0 ? -v2.Ma / v2.La / c2.hull.L * 100 : NaN });
    }

    // Spectrum chart data: raw w-spectrum vs incidence felt at the rear pair.
    const grid = logGrid(sp.omegaMin_radpm * 10, sp.omegaMax_radpm / 10, 120);
    const spectrum = [];
    for (let k = 0; k < grid.n; k++) {
      const Om = grid.om[k];
      const raw = chosen.phiW(Om) / (chosen.V * chosen.V);
      const t = chosen.vert.transfer(Om, !chosen.useDynV);
      spectrum.push({ Om, raw, felt: cabs2(t.finAlpha[t.finAlpha.length - 1]) * chosen.phiW(Om) });
    }

    // ---------- Warnings ----------
    const warnings = [];
    const lim = model.limits;
    if (model.turbulence && model.turbulence.placeholder)
      warnings.push("Turbulence presets are still the ATMOSPHERIC values from the reference script (σ_w = " + tu.sigma_w + " m/s, L_w = " + tu.L_w + " m). Replace them with measured underwater values before trusting the turbulence penalty.");
    const sigAdeg = Math.max(...chosen.sv.sigFinAlpha) * DEG;
    if (sigAdeg > lim.maxGustAlphaDeg)
      warnings.push(`RMS gust incidence at the fins ≈ ${round(sigAdeg, 1)}° exceeds ${lim.maxGustAlphaDeg}° – small-angle (linear) assumption is stretched.`);
    if (chosen.sl && chosen.lat.stable && ctx.rollFree && chosen.sl.sigPhi * DEG > 30)
      warnings.push(`RMS bank angle ≈ ${round(chosen.sl.sigPhi * DEG, 0)}° – far beyond small-angle theory. Uncontrolled, the vehicle would roll off into a spiral; the lateral wander figures are only an indication. More dihedral effect, a rolling-gust-tolerant layout or roll control is needed.`);
    [["depth", per1km.depth_m], ["lateral", per1km.lateral_m]].forEach(([n, y]) => {
      if (Number.isFinite(y) && y > 300)
        warnings.push(`RMS ${n} track wander after 1 km ≈ ${round(y, 0)} m exceeds 30 % of the distance – beyond small-angle theory: without control the vehicle does not hold its course. Use the mode data (stability, damping, time constants) rather than the wander figure.`);
    });
    if (chosen.Mach > MACH_CAP)
      warnings.push(`Mach ${round(chosen.Mach, 2)} > ${MACH_CAP}: transonic/supersonic – Prandtl–Glauert and the Korn/Lock wave-drag model are not valid here.`);
    else if (chosen.wave.CDA > 0)
      warnings.push(`Wave drag active at M ${round(chosen.Mach, 3)}: ${chosen.wave.surfaces.filter(t => t.CD > 0).map(t => t.name + " (M_crit " + round(t.Mcrit, 3) + ")").join(", ")} – ${round(chosen.drag.wave, 1)} N.`);
    if (!chosen.feasible)
      warnings.push(`Cruise speed ${round(chosen.V, 1)} m/s needs ${round(chosen.P_engine / 1000, 1)} kW engine power – above the ${round(design.maxPower_kW, 1)} kW available.`);
    if (staticMargin_m <= 0)
      warnings.push(`Statically UNSTABLE in pitch: neutral point (${round(xnp, 3)} m) is ahead of CG (${round(ctx.xcg, 3)} m). Move fins aft, enlarge the rear pair or move CG forward.`);
    else if (!chosen.vert.stable)
      warnings.push("Pitch/heave mode is dynamically unstable – trajectory dispersion is unbounded without control.");
    if (!chosen.hasYaw)
      warnings.push("No vertical surfaces (fin area 0): the hull alone is unstable in yaw (Munk moment) – lateral dispersion is unbounded.");
    else if (!chosen.lat.stable) {
      const L = chosen.lat;
      if (L.dutch && L.dutch.zeta <= 0) warnings.push(`Dutch roll is unstable (ζ ${round(L.dutch.zeta, 3)}) – more fin area/arm or less dihedral effect.`);
      if (L.spiral && L.spiral.lambda >= 0) warnings.push(`Spiral mode diverges (time to double ${round(L.spiral.T, 1)} s) – without roll control the path slowly banks away; lateral dispersion is unbounded. More dihedral effect (C_lβ more negative) or less fin stabilises it.`);
      if (L.roll && L.roll.lambda >= 0) warnings.push("Roll mode is unstable.");
      if (!L.dutch && !L.spiral && !L.roll) warnings.push("Lateral-directional motion is unstable – lateral dispersion is unbounded without control.");
    }
    if (chosen.CLf < 0 || chosen.CLr < 0)
      warnings.push(`CG (${round(ctx.xcg, 2)} m) lies outside the fin pairs – one pair carries negative lift (CL_front ${round(chosen.CLf, 3)}, CL_rear ${round(chosen.CLr, 3)}).`);
    const T = chosen.trim;
    chosen.surfV.map(t => [t.name, t.CL, t.sig, t.pair, t.st]).forEach(([n, CL, s, pair, st]) => {
      if (st.stalled)
        warnings.push(`The ${n} pair is beyond its stall limit already at trim (α ≈ ${round(st.alpha * DEG, 1)}°, CL ${round(CL, 2)}, CLmax ${round(pair.clMax, 2)}).`);
      else if (Math.abs(CL) + 2 * st.slope * s > pair.clMax)
        warnings.push(`The ${n} pair reaches CL ≈ ${round(Math.abs(CL) + 2 * st.slope * s, 2)} in 2σ gusts (CLmax ${round(pair.clMax, 2)}, ${pair.kind === "delta" ? "delta, stall α " + round(pair.stallAlpha * DEG, 0) + "°" : pair.section.name}) – risk of stall/ventilation.`);
      if (pair.kind === "delta" && pair.pf.sweepLEdeg < 45)
        warnings.push(`The ${n} pair uses the delta (vortex-lift) model with only ${round(pair.pf.sweepLEdeg, 0)}° LE sweep – the Polhamus analogy assumes a sharp, highly swept leading edge (≳ 45–50°).`);
    });
    if (ctx.downwash && ctx.rear && ctx.rear.span > ctx.front.span)
      warnings.push(`Downwash: the rear pair (b ${round(ctx.rear.span, 2)} m) is wider than the front pair (b ${round(ctx.front.span, 2)} m). The DATCOM formula assumes the tail lies inside the wing's vortex span; the outer tail sees upwash, so dε/dα (${round(T.dw.E, 2)}) is overestimated and the static margin underestimated. Consider switching downwash off for small front fins.`);
    if (T.dw.capped)
      warnings.push(`Downwash gradient dε/dα ≈ ${round(T.dw.E0, 2)} from the DATCOM formula was capped at 0.95 – the rear pair sits very close behind the front pair (l_H ${round(T.dw.lH, 2)} m), outside the formula's range.`);
    [["front", ctx.front], ["rear", ctx.rear], ["fin", ctx.fin]].filter(([, p]) => p).forEach(([n, pair]) => {
      const q = pair.section, sec = ((model.sections && model.sections.sections) || []).find(x => x.id === q.id);
      if (q.confidence === "rough")
        warnings.push(`The ${n} section (${q.name}) uses ROUGH placeholder coefficients – verify against the cited source in data/sections.json before using the result as justification.`);
      if (sec && sec.points.length > 1 && q.Re != null) {
        const lo = Math.min(...sec.points.map(p => p.Re)), hi = Math.max(...sec.points.map(p => p.Re));
        if (q.Re < lo || q.Re > hi)
          warnings.push(`The ${n} fin chord Re ≈ ${q.Re.toExponential(1)} is outside the ${q.name} data (${lo.toExponential(0)}–${hi.toExponential(0)}); the nearest data point is used. Check the kinematic viscosity (water ≈ 1.3e-6, air ≈ 1.5e-5 m²/s).`);
      }
    });
    [["Front", ctx.front.xcl], ["Rear", ctx.rear && ctx.rear.xcl], ["Fin", ctx.fin && ctx.fin.xcl], ["CG", ctx.xcg]].filter(([, x]) => x != null).forEach(([n, x]) => {
      if (x < 0 || x > ctx.hull.L) warnings.push(`${n} position ${round(x, 2)} m is outside the hull (0–${round(ctx.hull.L, 2)} m).`);
    });
    if (ctx.mEnd >= ctx.m0 || design.fuelMass_kg <= 0) warnings.push("No usable fuel – check fuel mass and reserve.");

    // ---------- Pedagogical steps ----------
    const D = chosen.drag;
    const steps = [
      { t: "Atmosphere", d: (ctx.atm ? `ISA h = ${round(ctx.atm.h, 0)} m: T ${round(ctx.atm.T - 273.15, 1)} °C, ρ ${round(ctx.rho, 4)} kg/m³, ν ${ctx.nu.toExponential(3)} m²/s, a ${round(ctx.a, 1)} m/s, σ ${round(ctx.atm.sigma, 3)} → available power ×${round(ctx.powerFactor, 3)}` : `manual: ρ ${round(ctx.rho, 2)} kg/m³, ν ${ctx.nu.toExponential(3)} m²/s, a ${round(ctx.a, 0)} m/s`) + ` · M = ${round(chosen.Mach, 3)}, β = ${round(Math.sqrt(Math.max(0, 1 - Math.min(chosen.Mach, MACH_CAP) ** 2)), 3)}` },
      { t: "Compressibility", d: chosen.wave.surfaces.map(t => `${t.name}: t/c ${round(t.pair.tc, 3)}, κ_A ${t.pair.kappaA}, M_dd ${round(t.Mdd, 3)}, M_crit ${round(t.Mcrit, 3)}, CD_wave ${t.CD.toExponential(2)}`).join(" · ") + ` · CLα front ${round(ctx.front.clAlpha0, 2)} → ${round(ctx.front.clAlpha, 2)}/rad (PG) · hull Cf ×${round(compressibleCf(chosen.Mach), 3)}` },
      { t: "Mass", d: `m0 = ${round(ctx.m0, 2)} kg (empty ${round(design.emptyMass_kg, 2)} + payload ${round(env.payload_kg || 0, 2)} + fuel ${round(design.fuelMass_kg, 2)}); m_end = ${round(ctx.mEnd, 2)} kg (reserve ${round(env.reserveFraction * 100, 0)} % of fuel)` },
      { t: "Hull", d: `L/D = ${round(ctx.hull.lam, 2)}, S_wet = ${round(ctx.hull.wetted, 3)} m², Vol = ${round(ctx.hull.volume * 1000, 1)} L, Re = ${chosen.Re.toExponential(2)}, Cf = ${round(chosen.Cf, 5)}, FF = ${round(ctx.hull.formFactor, 3)}` },
      { t: "Tail", d: { traditional: "traditional – horizontal tail (rear pair) + vertical fin", vtail: `V-tail – rear pair at dihedral Γ ${round(Math.acos(ctx.proj.rear.pitch) * DEG, 1)}° (pitch cos²Γ = ${round(ctx.proj.rear.pitch ** 2, 3)}, yaw sin²Γ = ${round(ctx.proj.rear.yaw ** 2, 3)})`, tailless: "tailless – wing only (elevon trim) + optional vertical fin", cruciform: "cruciform – both pairs also act in yaw" }[ctx.tailType] + (ctx.fin ? `; fin S ${round(ctx.fin.area, 4)} m², h ${round(ctx.fin.span, 3)} m, AR_geo ${round(ctx.fin.ARgeo, 2)} → AR_eff ${round(ctx.fin.AR, 2)}, CLα ${round(ctx.fin.clAlpha, 2)}/rad, a.c. x ${round(ctx.fin.xcl, 3)} m` : "") },
      { t: "Planform", d: ["front", "rear"].filter(n => ctx[n]).map(n => { const p = ctx[n], f = p.pf; return `${n} (${p.kind}): b = ${round(f.b, 3)} m, c_root ${round(f.cr, 3)} / c_tip ${round(f.ct, 3)} m, MAC ${round(f.mac, 3)} m, Λ_LE ${round(f.sweepLEdeg, 1)}°, Λ½c ${round(Math.atan(f.tanHc) * DEG, 1)}°, apex x ${round(f.xle, 3)} m → ¼-MAC x ${round(f.xAc25, 3)} m`; }).join(" · ") },
      { t: "Lift", d: ["front", "rear"].filter(n => ctx[n]).map(n => { const p = ctx[n], q = p.section, st = chosen.trim[n === "front" ? "fs" : "rs"]; return `${n}: ${q.name}, Re_MAC = ${q.Re != null ? q.Re.toExponential(2) : "–"}, cℓα = ${round(q.a0, 2)} → CLα = ${round(p.clAlpha, 2)}/rad` + (p.kind === "delta" ? ` (Kp), Kv = ${round(p.Kv, 2)}; trim α ${round(st.alpha * DEG, 2)}°, dCL/dα ${round(st.slope, 2)}/rad, a.c. x ${round(st.xAc, 3)} m` : `, e = ${round(p.e, 2)}`) + `, CD0 = ${round(p.cd0, 4)}, CLmax = ${round(p.clMax, 2)}`; }).join(" · ") },
      { t: "Downwash", d: !ctx.rear ? "no rear pair" : ctx.downwash ? `dε/dα at rear = ${round(chosen.trim.dw.E, 3)} (l_H ${round(chosen.trim.dw.lH, 3)} m, tail height ${round(ctx.tailHeight, 3)} m) → rear pair effectiveness ${round((1 - chosen.trim.dw.E) * 100, 0)} %` : "off" },
      { t: "Static balance (L = W)", d: ctx.rear
        ? `q = ${round(chosen.qbar, 0)} Pa · L_front = ${round(chosen.split.Lf, 1)} N (CL ${round(chosen.CLf, 3)}), L_rear = ${round(chosen.split.Lr, 1)} N (${ctx.tailType === "vtail" ? "panel CL " + round(chosen.CLr, 3) + ", Γ " + round(Math.acos(ctx.proj.rear.pitch) * DEG, 1) + "°" : "CL " + round(chosen.CLr, 3)})`
        : `q = ${round(chosen.qbar, 0)} Pa · tailless: wing carries W = ${round(chosen.split.Lf, 1)} N (CL ${round(chosen.CLf, 3)}); pitch trim by elevons assumed (trim drag not modelled)` },
      { t: "Gust incidence (freq. domain)", d: `σ_α ${chosen.surfV.map(t => t.name + " " + round(t.sig * DEG, 2) + "°").join(", ")}${chosen.surfL.length ? " · lateral " + chosen.surfL.map(t => t.name + " " + round(t.sig * DEG, 2) + "°").join(", ") : ""} (${chosen.useDynV ? "with vehicle response" : "fixed attitude, as script"})` },
      { t: "Drag build-up", d: `hull ${round(D.hull, 1)} + fin profile ${round(D.finProfile, 1)} + trim induced ${round(D.trim, 1)} + gust (vert) ${round(D.gustV, 1)} + gust (lat) ${round(D.gustL, 1)} + extra ${round(D.extra, 1)} = ${round(D.total, 1)} N` },
      { t: "Waterjet", d: `T = D → Vj = ${round(chosen.jet.Vj, 2)} m/s, η_F = ${round(chosen.jet.etaF, 3)}, P_jet = ${round(chosen.jet.Pjet / 1000, 2)} kW, P_engine = P_jet/η_jet = ${round(chosen.P_engine / 1000, 2)} kW (${round(chosen.throttle * 100, 0)} % of max)` },
      { t: "Breguet", d: `R = η_jet·η_F/(BSFC·g)·(L/D)·ln(m0/m_end) = ${round(chosen.etaTot, 3)}/(BSFC·g)·${round(chosen.LD, 3)}·${round(chosen.lnM, 4)} = ${round(chosen.rangeKm, 2)} km  (equiv. c = ${round(chosen.cEq_perHour, 3)} 1/h)` },
      { t: "Lateral", d: (ctx.rollFree ? "" : "roll held level (p = φ = 0) · ") + `C_nβ = ${round(lateralStatic.Cnb, 4)} (weathercock, > 0 stable), C_lβ = ${round(lateralStatic.Clb, 4)} (dihedral effect, < 0 stable), wing Γ ${round(ctx.wingDihedral * DEG, 1)}° · I_x ${round(ctx.Ix, 4)}, I_z ${round(ctx.Iz, 4)} kg·m²` + (ctx.Kphi ? `, roll stiffness ${round(ctx.Kphi, 2)} N·m/rad` : "") + ` · ` + [chosen.lat.dutch ? `Dutch roll ω_n ${round(chosen.lat.dutch.wn, 2)} rad/s ζ ${round(chosen.lat.dutch.zeta, 3)}` : "no Dutch-roll oscillation", chosen.lat.roll ? `roll τ ${round(chosen.lat.roll.tau, 3)} s` : null, chosen.lat.spiral ? `spiral λ ${round(chosen.lat.spiral.lambda, 4)} 1/s` : null, chosen.lat.coupled ? `roll–spiral coupled ω_n ${round(chosen.lat.coupled.wn, 2)}` : null].filter(Boolean).join(", ") },
      { t: "Stability", d: `x_NP = ${round(xnp, 3)} m, static margin = ${round(staticMargin_m, 3)} m (${round(staticMargin_pct, 1)} % L), Munk k2−k1 = ${round(ctx.hull.lamb.k2 - ctx.hull.lamb.k1, 3)}` }
    ];

    return {
      design, env, ctx, chosen, scan, best, vMax,
      xnp, staticMargin_m, staticMargin_pct, lateralStatic,
      disp, per1km, atRange, pathAngle, smSweep, smSweepOf, spectrum,
      altEnvelope: altitudeEnvelope(design, env, model, calib, opts),
      warnings, steps
    };
  }

  /*
   * One-way reach per course (180 headings) for the map. The "wind" inputs are
   * kept from the original tool (atmospheric placeholder; to be replaced by a
   * current model later). Time on task is fixed by the Breguet result:
   * t = R / V; ground distance = ground speed · t (crab angle included).
   */
  function envelope(r, params) {
    const COURSES = 180, STEP = 360 / COURSES;
    const V = r.chosen.V, W = params.windSpeed_mps || 0;
    const fromDeg = (params.windFromDeg != null ? params.windFromDeg : 0) % 360;
    const blowRad = ((fromDeg + 180) % 360) * Math.PI / 180;
    const t = r.chosen.rangeKm * 1000 / V;
    const headings = [], oneWayKm = [];
    for (let i = 0; i < COURSES; i++) {
      const th = i * STEP * Math.PI / 180;
      const along = W * Math.cos(th - blowRad), cross = W * Math.sin(th - blowRad);
      const a2 = V * V - cross * cross;
      const gs = a2 > 0 ? Math.sqrt(a2) + along : 0;
      headings.push(i * STEP); oneWayKm.push(gs > 0 ? gs * t / 1000 : 0);
    }
    const pos = oneWayKm.filter(x => x > 0);
    return { headings, oneWayKm, min: pos.length ? Math.min(...pos) : 0, max: pos.length ? Math.max(...pos) : 0 };
  }

  // Speed-scan range = the cruise-speed slider span in data/variables.json.
  function presetsFromVariables(vars) {
    const all = Object.values(vars.groups).flat();
    const v = all.find(d => d.key === "cruiseSpeed_mps");
    if (!v) throw new Error("variables.json: cruiseSpeed_mps is missing");
    return { speedScanMin_mps: v.min, speedScanMax_mps: v.max, speedScanStep_mps: (vars.speedScan && vars.speedScan.step_mps) || v.step };
  }

  return {
    analyze, envelope, round, presetsFromVariables,
    // exposed for tests / docs
    normalizeDesign, isa,
    _internal: { waveDrag, altitudeEnvelope, resolveAtmosphere, phiLongitudinal, phiTransverse, phiRolling, charPoly, polyRoots, csolve, pairRollDerivs, lateralSystem, logGrid, lambCoefficients, helmbold, sectionAt, finPair, planform, liftState, deltaCL, downwashGradient, trim, jetForThrust, thrustForJetPower, ittcCf, hullGeometry, planeSystem, planeStats }
  };
})();

if (typeof module !== "undefined") module.exports = Model;
