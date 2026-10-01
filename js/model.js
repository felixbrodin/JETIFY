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
  // Helmbold low-AR lift-curve slope [1/rad].
  const helmbold = (AR) => 2 * Math.PI * AR / (2 + Math.sqrt(AR * AR + 4));

  function finPair(area, AR, xcl, clAlphaOverride) {
    const span = Math.sqrt(area * AR);
    const chord = Math.sqrt(area / AR);
    const clAlpha = clAlphaOverride != null ? clAlphaOverride : helmbold(AR);
    return { area, AR, xcl, span, chord, clAlpha, clAlphaSrc: clAlphaOverride != null ? "empirical" : "Helmbold(AR)" };
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
  function planeSystem(p) {
    // p: { V, qbar, fins:[{K, l, chord, x}], Mah, Mqx, mPrime, Iprime, hullL, fkMass, fkArm }
    const V = p.V;
    let La = 0, Lq = 0, Ma = p.Mah, Mq = p.Mqx;
    p.fins.forEach(f => { La += f.K; Lq += f.K * f.l / V; Ma -= f.K * f.l; Mq -= f.K * f.l * f.l / V; });
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
      const fw = p.fins.map(f => cscale(cexpi(-Om * f.x), sinc(Om * f.chord / 2) / V)); // w_i/V per unit w
      const hw = cscale(cexpi(-Om * p.hullL / 2), sinc(Om * p.hullL / 2) / V);
      const Ffk = cmul(cx(0, w), cscale(hw, p.fkMass * V));          // (1+k2)ρVol · iω · w_hull
      let Lw = Ffk, Mw = cadd(cscale(hw, p.Mah), cscale(Ffk, -p.fkArm));
      p.fins.forEach((f, i) => { Lw = cadd(Lw, cscale(fw[i], f.K)); Mw = cadd(Mw, cscale(fw[i], -f.K * f.l)); });
      let alpha = cx(0), q = cx(0);
      if (!rigid) {
        // (iω − a11)α − a12 q = −Lw/mV ;  −a21 α + (iω − a22) q = Mw/I
        const A11 = cx(-a11, w), A12 = cx(-a12), A21 = cx(-a21), A22 = cx(-a22, w);
        const b1 = cscale(Lw, -1 / mV), b2 = cscale(Mw, 1 / p.Iprime);
        const D = cadd(cmul(A11, A22), cscale(cmul(A12, A21), -1));
        alpha = cdiv(cadd(cmul(b1, A22), cscale(cmul(A12, b2), -1)), D);
        q = cdiv(cadd(cmul(A11, b2), cscale(cmul(A21, b1), -1)), D);
      }
      // Effective incidence at each fin, lift perturbation, path-rate.
      const finAlpha = p.fins.map((f, i) => cadd(cadd(alpha, cscale(q, f.l / V)), fw[i]));
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
    const hull = hullGeometry(design);
    const front = finPair(design.frontArea_m2, design.frontAR, design.frontX_m, c.clAlphaFront);
    const rear = finPair(design.rearArea_m2, design.rearAR, design.rearX_m, c.clAlphaRear);
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
      e: env.oswaldE, cd0Fin: env.finCD0, clMax: env.finCLmax, nu: env.kinematicViscosity_m2ps,
      extraDragArea: (c.extraDragArea_m2 != null ? c.extraDragArea_m2 : (hy.extraDragArea_m2 || 0)),
      hullCDwetOverride: c.hullCDwet != null ? c.hullCDwet : null,
      cruciform: !!design.cruciform,
      extraDamping: design.extraPitchDamping_Nms || 0
    };
  }

  // Lift split between the two pairs (moment balance about CG, L_f + L_r = W).
  function liftSplit(ctx) {
    const xf = ctx.front.xcl, xr = ctx.rear.xcl, dx = xr - xf;
    if (Math.abs(dx) < 1e-6) return { Lf: ctx.W / 2, Lr: ctx.W / 2, degenerate: true };
    const Lf = ctx.W * (xr - ctx.xcg) / dx;
    return { Lf, Lr: ctx.W - Lf, degenerate: false };
  }

  // Everything at one speed. opts.withDynamics: use vehicle response in drag/load penalty.
  function pointAt(V, ctx, env, model, opts) {
    const sp = model.spectral;
    const grid = logGrid(sp.omegaMin_radpm, sp.omegaMax_radpm, sp.points);
    const qbar = 0.5 * ctx.rho * V * V;
    const split = liftSplit(ctx);
    const CLf = split.Lf / (qbar * ctx.front.area);
    const CLr = split.Lr / (qbar * ctx.rear.area);

    // Plane systems (vertical uses w-gusts; lateral uses v-gusts if cruciform).
    const Mah = (ctx.hull.lamb.k2 - ctx.hull.lamb.k1) * ctx.rho * ctx.hull.volume * V * V;
    const mk = (finsOn) => planeSystem({
      V, qbar,
      fins: finsOn ? [
        { K: qbar * ctx.front.area * ctx.front.clAlpha, l: ctx.front.xcl - ctx.xcg, chord: ctx.front.chord, x: ctx.front.xcl },
        { K: qbar * ctx.rear.area * ctx.rear.clAlpha, l: ctx.rear.xcl - ctx.xcg, chord: ctx.rear.chord, x: ctx.rear.xcl }
      ] : [],
      Mah, Mqx: -Math.abs(ctx.extraDamping), mPrime: ctx.m0 + ctx.mAdd, Iprime: ctx.I0 + ctx.Iadd, hullL: ctx.hull.L,
      fkMass: (1 + ctx.hull.lamb.k2) * ctx.rho * ctx.hull.volume, fkArm: ctx.hull.L / 2 - ctx.xcg
    });
    const vert = mk(true);
    const lat = mk(ctx.cruciform);
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
    const D_finProfile = qPar * ctx.cd0Fin * (ctx.front.area + ctx.rear.area) * planes;
    const ind = (pair, CL) => qbar * pair.area * CL * CL / (Math.PI * ctx.e * pair.AR);
    const D_trim = ind(ctx.front, CLf) + ind(ctx.rear, CLr);
    const gust = (pair, sig) => qbar * pair.area * pair.clAlpha * pair.clAlpha * sig * sig / (Math.PI * ctx.e * pair.AR);
    const D_gustV = gust(ctx.front, sv.sigFinAlpha[0]) + gust(ctx.rear, sv.sigFinAlpha[1]);
    const D_gustL = sl ? gust(ctx.front, sl.sigFinAlpha[0]) + gust(ctx.rear, sl.sigFinAlpha[1]) : 0;
    const D_extra = qPar * ctx.extraDragArea;
    const D = D_hull + D_finProfile + D_trim + D_gustV + D_gustL + D_extra;
    const D_calm = qbar * ctx.hull.wetted * CDwet + qbar * ctx.cd0Fin * (ctx.front.area + ctx.rear.area) * planes + D_trim + qbar * ctx.extraDragArea;

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
      V, qbar, qPar, Re, Cf, CDwet, split, CLf, CLr,
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

    // Static-margin design sweep vs rear-fin position.
    const smSweep = [];
    const nS = 40;
    for (let i = 0; i <= nS; i++) {
      const xr = ctx.hull.L * i / nS;
      const d2 = { ...design, rearX_m: xr };
      const c2 = resolve(d2, env, model, calib);
      const qb = 0.5 * c2.rho;   // V² cancels in the ratio
      const Kf = qb * c2.front.area * c2.front.clAlpha, Kr = qb * c2.rear.area * c2.rear.clAlpha;
      const Mah = (c2.hull.lamb.k2 - c2.hull.lamb.k1) * c2.rho * c2.hull.volume;
      const num = Kf * (c2.front.xcl - c2.xcg) + Kr * (c2.rear.xcl - c2.xcg) - Mah;
      smSweep.push({ x: xr, sm: num / (Kf + Kr) / c2.hull.L * 100 });
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
    [["front", chosen.CLf, chosen.sv.sigFinAlpha[0], ctx.front], ["rear", chosen.CLr, chosen.sv.sigFinAlpha[1], ctx.rear]].forEach(([n, CL, s, pair]) => {
      if (Math.abs(CL) + 2 * pair.clAlpha * s > ctx.clMax)
        warnings.push(`The ${n} pair reaches CL ≈ ${round(Math.abs(CL) + 2 * pair.clAlpha * s, 2)} in 2σ gusts (CLmax ${ctx.clMax}) – risk of stall/ventilation.`);
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
      { t: "Fins", d: `front: b = ${round(ctx.front.span, 3)} m, c = ${round(ctx.front.chord, 3)} m, CLα = ${round(ctx.front.clAlpha, 2)}/rad · rear: b = ${round(ctx.rear.span, 3)} m, c = ${round(ctx.rear.chord, 3)} m, CLα = ${round(ctx.rear.clAlpha, 2)}/rad` },
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

  return {
    analyze, envelope, round,
    // exposed for tests / docs
    _internal: { phiLongitudinal, phiTransverse, logGrid, lambCoefficients, helmbold, jetForThrust, thrustForJetPower, ittcCf, hullGeometry, planeSystem, planeStats }
  };
})();

if (typeof module !== "undefined") module.exports = Model;
