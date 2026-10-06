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

  // ---------- Control-surface pairs ----------
  // Generalised Helmbold / DATCOM lift-curve slope [1/rad] (incompressible):
  //   CLα = 2πAR / (2 + √(AR²/κ²·(1 + tan²Λ½c) + 4)),  κ = a0/2π.
  // With a0 = 2π and Λ = 0 this is the classic 2πAR/(2+√(AR²+4)).
  const helmbold = (AR, a0, tanHalfChord) => {
    const k = (a0 != null ? a0 : 2 * Math.PI) / (2 * Math.PI);
    const t = tanHalfChord || 0;
    return 2 * Math.PI * AR / (2 + Math.sqrt(AR * AR / (k * k) * (1 + t * t) + 4));
  };

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

  // Fill planform defaults (model.json designDefaults) and convert legacy frontX_m / rearX_m
  // (old profiles gave the centre of lift of a rectangular pair) to the root leading edge.
  function normalizeDesign(design, model) {
    const defs = (model && model.designDefaults) || {};
    const d = { ...design };
    Object.keys(defs).forEach(k => { if (d[k] == null) d[k] = defs[k]; });
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
    const clAlpha = clAlphaOverride != null ? clAlphaOverride : helmbold(AR, s2.clAlpha, pf.tanHc);
    const cd0 = s2.cdMin != null ? s2.cdMin : env.finCD0;
    const e = design[n + "OswaldE"] != null && design[n + "OswaldE"] > 0 ? design[n + "OswaldE"] : env.oswaldE;
    const kind = design[n + "Planform"] === "delta" ? "delta" : "conventional";
    const dm = (model && model.delta) || {};
    const pair = {
      n, area, AR, pf, span: pf.b, chord: pf.mac, xcl: pf.xAc25, clAlpha, cd0, e, kind,
      Kv: dm.Kv != null ? dm.Kv : Math.PI, stallAlpha: (dm.stallAlphaDeg != null ? dm.stallAlphaDeg : 30) / DEG,
      section: { id: sec ? sec.id : "ideal", name: sec ? sec.name : "Ideal thin airfoil (2π)", confidence: sec ? sec.confidence : "typical", Re, a0: s2.clAlpha, cdMin: s2.cdMin, clMax2D: s2.clMax },
      clAlphaSrc: clAlphaOverride != null ? "empirical" : "DATCOM(AR, Λ½c, a0)"
    };
    const cosQc = 1 / Math.sqrt(1 + pf.tanQc * pf.tanQc);
    pair.clMax = kind === "delta" ? deltaCL(pair, pair.stallAlpha)
      : (s2.clMax != null ? CLMAX_3D * s2.clMax * cosQc : env.finCLmax);
    return pair;
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
    const E0 = 4.44 * Math.pow(KA * Kl * KH * Math.sqrt(cosQc), 1.19) * (fs.slope / f.clAlpha);
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
    const rear = finPair("rear", design, c.clAlphaRear, secOf(design.rearSection), env, model);
    const m0 = design.emptyMass_kg + (env.payload_kg || 0) + design.fuelMass_kg;
    const mEnd = m0 - design.fuelMass_kg * (1 - env.reserveFraction);
    const W = m0 * g;
    const xcg = design.xcg_m;
    // Pitch inertia: override or uniform solid cylinder; added mass/inertia from Lamb.
    const I0 = design.pitchInertia_kgm2 != null && design.pitchInertia_kgm2 > 0
      ? design.pitchInertia_kgm2 : m0 * (hull.L * hull.L / 12 + hull.D * hull.D / 16);
    const a = hull.L / 2, b = hull.D / 2;
    const Iadd = hull.lamb.kp * rho * hull.volume * (a * a + b * b) / 5;
    const mAdd = hull.lamb.k2 * rho * hull.volume;
    const An = Math.PI / 4 * design.nozzleDiameter_m * design.nozzleDiameter_m;
    const etaJet = c.etaJet != null ? c.etaJet : design.jetEfficiency;
    return {
      g, rho, hull, front, rear, m0, mEnd, W, xcg, I0, Iadd, mAdd, An, etaJet,
      nu: env.kinematicViscosity_m2ps,
      extraDragArea: (c.extraDragArea_m2 != null ? c.extraDragArea_m2 : (hy.extraDragArea_m2 || 0)),
      hullCDwetOverride: c.hullCDwet != null ? c.hullCDwet : null,
      cruciform: !!design.cruciform,
      downwash: design.downwash !== false, tailHeight: design.rearHeight_m || 0,
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
  function trim(ctx, qbar) {
    let xf = ctx.front.pf.xAc25, xr = ctx.rear.pf.xAc25, split, CLf, CLr, fs, rs;
    for (let it = 0; it < 8; it++) {
      split = liftSplit(ctx, xf, xr);
      CLf = split.Lf / (qbar * ctx.front.area);
      CLr = split.Lr / (qbar * ctx.rear.area);
      fs = liftState(ctx.front, CLf); rs = liftState(ctx.rear, CLr);
      if (Math.abs(fs.xLift - xf) < 1e-7 && Math.abs(rs.xLift - xr) < 1e-7) break;
      xf = fs.xLift; xr = rs.xLift;
    }
    return { split, CLf, CLr, fs, rs, dw: downwashGradient(ctx, fs, rs) };
  }

  // Pitch/heave (vertical) and yaw/sway (lateral) systems at speed V.
  function planeSystems(ctx, V, qbar, tr) {
    const Mah = (ctx.hull.lamb.k2 - ctx.hull.lamb.k1) * ctx.rho * ctx.hull.volume * V * V;
    const base = {
      V, qbar, Mah, Mqx: -Math.abs(ctx.extraDamping), mPrime: ctx.m0 + ctx.mAdd, Iprime: ctx.I0 + ctx.Iadd, hullL: ctx.hull.L,
      fkMass: (1 + ctx.hull.lamb.k2) * ctx.rho * ctx.hull.volume, fkArm: ctx.hull.L / 2 - ctx.xcg
    };
    const fin = (pair, slope, x, dw) => ({ K: qbar * pair.area * slope, l: x - ctx.xcg, chord: pair.chord, x, dw });
    const dw = tr.dw.E > 0 ? { src: 0, E: tr.dw.E, tau: tr.dw.lH / V } : null;
    const vert = planeSystem({ ...base, fins: [fin(ctx.front, tr.fs.slope, tr.fs.xAc), fin(ctx.rear, tr.rs.slope, tr.rs.xAc, dw)] });
    // Lateral surfaces fly at zero side-force trim: attached-flow slope at ¼-MAC, no sidewash.
    const lat = planeSystem({ ...base, fins: ctx.cruciform ? [fin(ctx.front, ctx.front.clAlpha, ctx.front.pf.xAc25), fin(ctx.rear, ctx.rear.clAlpha, ctx.rear.pf.xAc25)] : [] });
    return { vert, lat };
  }

  // Everything at one speed. opts.withDynamics: use vehicle response in drag/load penalty.
  function pointAt(V, ctx, env, model, opts) {
    const sp = model.spectral;
    const grid = logGrid(sp.omegaMin_radpm, sp.omegaMax_radpm, sp.points);
    const qbar = 0.5 * ctx.rho * V * V;
    const tr = trim(ctx, qbar);
    const split = tr.split, CLf = tr.CLf, CLr = tr.CLr;
    const { vert, lat } = planeSystems(ctx, V, qbar, tr);
    const tu = env.turbulence;
    const phiW = (Om) => phiTransverse(Om, tu.sigma_w, tu.L_w);
    const phiV = (Om) => phiTransverse(Om, tu.sigma_v, tu.L_v);
    const phiU = (Om) => phiLongitudinal(Om, tu.sigma_u, tu.L_u);

    const useDynV = opts.withDynamics && vert.stable;
    const useDynL = opts.withDynamics && lat.stable;
    const sv = planeStats(vert, phiW, grid, !useDynV);
    const sl = ctx.cruciform ? planeStats(lat, phiV, grid, !useDynL) : null;

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
    const Cf = ittcCf(Re);
    const CDwet = ctx.hullCDwetOverride != null ? ctx.hullCDwetOverride : Cf * ctx.hull.formFactor;
    const D_hull = qPar * ctx.hull.wetted * CDwet;
    const planes = ctx.cruciform ? 2 : 1;
    const finCD0A = ctx.front.cd0 * ctx.front.area + ctx.rear.cd0 * ctx.rear.area;
    const D_finProfile = qPar * finCD0A * planes;
    // Lift-dependent drag: E[CDi(α0 + α_g)] ≈ CDi(α0) + ½·CDi''(α0)·σα²  (exact for the parabolic polar).
    const D_trim = qbar * (ctx.front.area * tr.fs.CDi + ctx.rear.area * tr.rs.CDi);
    const gust = (pair, st, sig) => qbar * pair.area * 0.5 * st.d2 * sig * sig;
    const D_gustV = gust(ctx.front, tr.fs, sv.sigFinAlpha[0]) + gust(ctx.rear, tr.rs, sv.sigFinAlpha[1]);
    const D_gustL = sl ? gust(ctx.front, liftState(ctx.front, 0), sl.sigFinAlpha[0]) + gust(ctx.rear, liftState(ctx.rear, 0), sl.sigFinAlpha[1]) : 0;
    const D_extra = qPar * ctx.extraDragArea;
    const D = D_hull + D_finProfile + D_trim + D_gustV + D_gustL + D_extra;
    const D_calm = qbar * ctx.hull.wetted * CDwet + qbar * finCD0A * planes + D_trim + qbar * ctx.extraDragArea;

    // Propulsion
    const jet = jetForThrust(D, V, ctx.rho, ctx.An);
    const P_engine = jet.Pjet / ctx.etaJet;
    const etaTot = ctx.etaJet * jet.etaF;
    const Pmax = opts.maxPower_W;
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
      V, qbar, qPar, Re, Cf, CDwet, split, CLf, CLr, trim: tr,
      vert, lat, sv, sl, useDynV, useDynL,
      drag: { hull: D_hull, finProfile: D_finProfile, trim: D_trim, gustV: D_gustV, gustL: D_gustL, extra: D_extra, total: D, calm: D_calm },
      jet, P_engine, etaTot, avail, throttle, feasible: throttle <= 1,
      LD, lnM, rangeKm: rangeM / 1000, enduranceH: rangeM / V / 3600, cEq_perHour, fuelFlow_kgph,
      phiW, phiV
    };
  }

  // ---------- Main analysis ----------
  function analyze(design, env, model, calib) {
    if (!design || !(design.hullLength_m > 0) || !(design.hullDiameter_m > 0)) return null;
    if (!(design.frontArea_m2 > 0) || !(design.rearArea_m2 > 0) || !(design.nozzleDiameter_m > 0)) return null;
    design = normalizeDesign(design, model);
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
      const p = pointAt(V, ctx, env, model, opts);
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

    // Trajectory dispersion with no steering.
    const sp = model.spectral;
    const distances = [];
    const nD = sp.dispersionPoints;
    const xEnd = Math.max(0.1, Math.min(chosen.rangeKm, sp.dispersionMaxKm)) * 1000;
    for (let i = 1; i <= nD; i++) distances.push(xEnd * i / nD);
    const tu = env.turbulence;
    const disp = distances.map(X => ({
      km: X / 1000,
      depth_m: dispersion(chosen.vert, chosen.phiW, chosen.V, X, sp),
      lateral_m: ctx.cruciform ? dispersion(chosen.lat, chosen.phiV, chosen.V, X, sp) : Infinity
    }));
    const per1km = {
      depth_m: dispersion(chosen.vert, chosen.phiW, chosen.V, 1000, sp),
      lateral_m: ctx.cruciform ? dispersion(chosen.lat, chosen.phiV, chosen.V, 1000, sp) : Infinity
    };
    const atRange = disp.length ? disp[disp.length - 1] : null;
    const pathAngle = {
      vert_rad: pathAngleRms(chosen.vert, chosen.phiW, chosen.V, sp),
      lat_rad: ctx.cruciform ? pathAngleRms(chosen.lat, chosen.phiV, chosen.V, sp) : Infinity
    };

    // Static-margin design sweep vs rear-pair root-LE position (same speed, trim and downwash model).
    const smSweep = [];
    const nS = 40;
    for (let i = 0; i <= nS; i++) {
      const xr = ctx.hull.L * i / nS;
      const c2 = resolve({ ...design, rearXle_m: xr }, env, model, calib);
      const tr2 = trim(c2, chosen.qbar);
      const v2 = planeSystems(c2, chosen.V, chosen.qbar, tr2).vert;
      smSweep.push({ x: xr, sm: v2.La > 0 ? -v2.Ma / v2.La / c2.hull.L * 100 : NaN });
    }

    // Spectrum chart data: raw w-spectrum vs incidence felt at the rear pair.
    const grid = logGrid(sp.omegaMin_radpm * 10, sp.omegaMax_radpm / 10, 120);
    const spectrum = [];
    for (let k = 0; k < grid.n; k++) {
      const Om = grid.om[k];
      const raw = chosen.phiW(Om) / (chosen.V * chosen.V);
      const t = chosen.vert.transfer(Om, !chosen.useDynV);
      spectrum.push({ Om, raw, felt: cabs2(t.finAlpha[1]) * chosen.phiW(Om) });
    }

    // ---------- Warnings ----------
    const warnings = [];
    const lim = model.limits;
    if (model.turbulence && model.turbulence.placeholder)
      warnings.push("Turbulence presets are still the ATMOSPHERIC values from the reference script (σ_w = " + tu.sigma_w + " m/s, L_w = " + tu.L_w + " m). Replace them with measured underwater values before trusting the turbulence penalty.");
    const sigAdeg = Math.max(...chosen.sv.sigFinAlpha) * DEG;
    if (sigAdeg > lim.maxGustAlphaDeg)
      warnings.push(`RMS gust incidence at the fins ≈ ${round(sigAdeg, 1)}° exceeds ${lim.maxGustAlphaDeg}° – small-angle (linear) assumption is stretched.`);
    if (!chosen.feasible)
      warnings.push(`Cruise speed ${round(chosen.V, 1)} m/s needs ${round(chosen.P_engine / 1000, 1)} kW engine power – above the ${round(design.maxPower_kW, 1)} kW available.`);
    if (staticMargin_m <= 0)
      warnings.push(`Statically UNSTABLE in pitch: neutral point (${round(xnp, 3)} m) is ahead of CG (${round(ctx.xcg, 3)} m). Move fins aft, enlarge the rear pair or move CG forward.`);
    else if (!chosen.vert.stable)
      warnings.push("Pitch/heave mode is dynamically unstable – trajectory dispersion is unbounded without control.");
    if (!ctx.cruciform)
      warnings.push("No vertical surfaces (cruciform off): the hull alone is unstable in yaw (Munk moment) – lateral dispersion is unbounded.");
    else if (!chosen.lat.stable)
      warnings.push("Yaw/sway mode is unstable – lateral dispersion is unbounded without control.");
    if (chosen.CLf < 0 || chosen.CLr < 0)
      warnings.push(`CG (${round(ctx.xcg, 2)} m) lies outside the fin pairs – one pair carries negative lift (CL_front ${round(chosen.CLf, 3)}, CL_rear ${round(chosen.CLr, 3)}).`);
    const T = chosen.trim;
    [["front", chosen.CLf, chosen.sv.sigFinAlpha[0], ctx.front, T.fs], ["rear", chosen.CLr, chosen.sv.sigFinAlpha[1], ctx.rear, T.rs]].forEach(([n, CL, s, pair, st]) => {
      if (st.stalled)
        warnings.push(`The ${n} pair is beyond its stall limit already at trim (α ≈ ${round(st.alpha * DEG, 1)}°, CL ${round(CL, 2)}, CLmax ${round(pair.clMax, 2)}).`);
      else if (Math.abs(CL) + 2 * st.slope * s > pair.clMax)
        warnings.push(`The ${n} pair reaches CL ≈ ${round(Math.abs(CL) + 2 * st.slope * s, 2)} in 2σ gusts (CLmax ${round(pair.clMax, 2)}, ${pair.kind === "delta" ? "delta, stall α " + round(pair.stallAlpha * DEG, 0) + "°" : pair.section.name}) – risk of stall/ventilation.`);
      if (pair.kind === "delta" && pair.pf.sweepLEdeg < 45)
        warnings.push(`The ${n} pair uses the delta (vortex-lift) model with only ${round(pair.pf.sweepLEdeg, 0)}° LE sweep – the Polhamus analogy assumes a sharp, highly swept leading edge (≳ 45–50°).`);
    });
    if (T.dw.capped)
      warnings.push(`Downwash gradient dε/dα ≈ ${round(T.dw.E0, 2)} from the DATCOM formula was capped at 0.95 – the rear pair sits very close behind the front pair (l_H ${round(T.dw.lH, 2)} m), outside the formula's range.`);
    [["front", ctx.front], ["rear", ctx.rear]].forEach(([n, pair]) => {
      const q = pair.section, sec = ((model.sections && model.sections.sections) || []).find(x => x.id === q.id);
      if (q.confidence === "rough")
        warnings.push(`The ${n} section (${q.name}) uses ROUGH placeholder coefficients – verify against the cited source in data/sections.json before using the result as justification.`);
      if (sec && sec.points.length > 1 && q.Re != null) {
        const lo = Math.min(...sec.points.map(p => p.Re)), hi = Math.max(...sec.points.map(p => p.Re));
        if (q.Re < lo || q.Re > hi)
          warnings.push(`The ${n} fin chord Re ≈ ${q.Re.toExponential(1)} is outside the ${q.name} data (${lo.toExponential(0)}–${hi.toExponential(0)}); the nearest data point is used. Check the kinematic viscosity (water ≈ 1.3e-6, air ≈ 1.5e-5 m²/s).`);
      }
    });
    [["Front", ctx.front.xcl], ["Rear", ctx.rear.xcl], ["CG", ctx.xcg]].forEach(([n, x]) => {
      if (x < 0 || x > ctx.hull.L) warnings.push(`${n} position ${round(x, 2)} m is outside the hull (0–${round(ctx.hull.L, 2)} m).`);
    });
    if (ctx.mEnd >= ctx.m0 || design.fuelMass_kg <= 0) warnings.push("No usable fuel – check fuel mass and reserve.");

    // ---------- Pedagogical steps ----------
    const D = chosen.drag;
    const steps = [
      { t: "Mass", d: `m0 = ${round(ctx.m0, 2)} kg (empty ${round(design.emptyMass_kg, 2)} + payload ${round(env.payload_kg || 0, 2)} + fuel ${round(design.fuelMass_kg, 2)}); m_end = ${round(ctx.mEnd, 2)} kg (reserve ${round(env.reserveFraction * 100, 0)} % of fuel)` },
      { t: "Hull", d: `L/D = ${round(ctx.hull.lam, 2)}, S_wet = ${round(ctx.hull.wetted, 3)} m², Vol = ${round(ctx.hull.volume * 1000, 1)} L, Re = ${chosen.Re.toExponential(2)}, Cf = ${round(chosen.Cf, 5)}, FF = ${round(ctx.hull.formFactor, 3)}` },
      { t: "Planform", d: ["front", "rear"].map(n => { const p = ctx[n], f = p.pf; return `${n} (${p.kind}): b = ${round(f.b, 3)} m, c_root ${round(f.cr, 3)} / c_tip ${round(f.ct, 3)} m, MAC ${round(f.mac, 3)} m, Λ_LE ${round(f.sweepLEdeg, 1)}°, Λ½c ${round(Math.atan(f.tanHc) * DEG, 1)}°, apex x ${round(f.xle, 3)} m → ¼-MAC x ${round(f.xAc25, 3)} m`; }).join(" · ") },
      { t: "Lift", d: ["front", "rear"].map(n => { const p = ctx[n], q = p.section, st = chosen.trim[n === "front" ? "fs" : "rs"]; return `${n}: ${q.name}, Re_MAC = ${q.Re != null ? q.Re.toExponential(2) : "–"}, cℓα = ${round(q.a0, 2)} → CLα = ${round(p.clAlpha, 2)}/rad` + (p.kind === "delta" ? ` (Kp), Kv = ${round(p.Kv, 2)}; trim α ${round(st.alpha * DEG, 2)}°, dCL/dα ${round(st.slope, 2)}/rad, a.c. x ${round(st.xAc, 3)} m` : `, e = ${round(p.e, 2)}`) + `, CD0 = ${round(p.cd0, 4)}, CLmax = ${round(p.clMax, 2)}`; }).join(" · ") },
      { t: "Downwash", d: ctx.downwash ? `dε/dα at rear = ${round(chosen.trim.dw.E, 3)} (l_H ${round(chosen.trim.dw.lH, 3)} m, tail height ${round(ctx.tailHeight, 3)} m) → rear pair effectiveness ${round((1 - chosen.trim.dw.E) * 100, 0)} %` : "off" },
      { t: "Static balance (L = W)", d: `q = ${round(chosen.qbar, 0)} Pa · L_front = ${round(chosen.split.Lf, 1)} N (CL ${round(chosen.CLf, 3)}), L_rear = ${round(chosen.split.Lr, 1)} N (CL ${round(chosen.CLr, 3)})` },
      { t: "Gust incidence (freq. domain)", d: `σ_α front ${round(chosen.sv.sigFinAlpha[0] * DEG, 2)}°, rear ${round(chosen.sv.sigFinAlpha[1] * DEG, 2)}° (${chosen.useDynV ? "with vehicle response" : "fixed attitude, as script"})` },
      { t: "Drag build-up", d: `hull ${round(D.hull, 1)} + fin profile ${round(D.finProfile, 1)} + trim induced ${round(D.trim, 1)} + gust (vert) ${round(D.gustV, 1)} + gust (lat) ${round(D.gustL, 1)} + extra ${round(D.extra, 1)} = ${round(D.total, 1)} N` },
      { t: "Waterjet", d: `T = D → Vj = ${round(chosen.jet.Vj, 2)} m/s, η_F = ${round(chosen.jet.etaF, 3)}, P_jet = ${round(chosen.jet.Pjet / 1000, 2)} kW, P_engine = P_jet/η_jet = ${round(chosen.P_engine / 1000, 2)} kW (${round(chosen.throttle * 100, 0)} % of max)` },
      { t: "Breguet", d: `R = η_jet·η_F/(BSFC·g)·(L/D)·ln(m0/m_end) = ${round(chosen.etaTot, 3)}/(BSFC·g)·${round(chosen.LD, 3)}·${round(chosen.lnM, 4)} = ${round(chosen.rangeKm, 2)} km  (equiv. c = ${round(chosen.cEq_perHour, 3)} 1/h)` },
      { t: "Stability", d: `x_NP = ${round(xnp, 3)} m, static margin = ${round(staticMargin_m, 3)} m (${round(staticMargin_pct, 1)} % L), Munk k2−k1 = ${round(ctx.hull.lamb.k2 - ctx.hull.lamb.k1, 3)}` }
    ];

    return {
      design, env, ctx, chosen, scan, best, vMax,
      xnp, staticMargin_m, staticMargin_pct,
      disp, per1km, atRange, pathAngle, smSweep, spectrum,
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
    normalizeDesign,
    _internal: { phiLongitudinal, phiTransverse, logGrid, lambCoefficients, helmbold, sectionAt, finPair, planform, liftState, deltaCL, downwashGradient, trim, jetForThrust, thrustForJetPower, ittcCf, hullGeometry, planeSystem, planeStats }
  };
})();

if (typeof module !== "undefined") module.exports = Model;
