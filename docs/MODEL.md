# Jetify – physics model

Quasi-steady mission-planning model for a waterjet-propelled, jet-fuel AUV.
No time simulation: turbulence statistics are obtained by integrating the
von Kármán spectrum in the frequency domain. Implementation: `js/model.js`.
All coefficients and defaults: `data/model.json`; slider spans: `data/variables.json`.

Origin: the MATLAB reference script "VON KARMAN TURBULENCE MODEL" (first-order
shaping filter, 120 s Monte-Carlo, Breguet range). Jetify keeps its static
equations and replaces the time simulation with spectral integration.

---

## 1. Static balance (unchanged from the script)

`L_front + L_rear = W = m0·g` – buoyancy is deliberately excluded (handled in a
separate calculation). The lift split follows from moment balance about CG:

```
L_front = W·(x_rear − x_CG)/(x_rear − x_front)
L_rear  = W − L_front
CL_i    = L_i / (q·S_i),   q = ½ρV0²
```

x_front / x_rear are the lift centres: ¼-MAC for conventional pairs; for a delta
the lift centre moves aft with α (vortex lift), so split and α are iterated to a
fixed point. With one pair at CG this reduces to the script's `CLtrim = W/(½ρV0²S)`.
Trim is evaluated at take-off mass m0 (as in the script).

## 2. Geometry

| Item | Relation |
|---|---|
| Planform (trapezoid, S = both halves) | `b = √(S·AR)`, `c_r = 2S/(b(1+λ))`, `c_t = λc_r`, `MAC = ⅔c_r(1+λ+λ²)/(1+λ)`, `ȳ_MAC = b/6·(1+2λ)/(1+λ)` |
| Sweep of the n-chord line | `tanΛ_n = tanΛ_LE − 4n(1−λ)/(AR(1+λ))` |
| Position input | root leading edge (apex) x_le; `x_ac = x_le + ȳ_MAC·tanΛ_LE + ¼MAC`; planform centroid `x_le + b/2·tanΛ_LE(1+2λ)/(3(1+λ)) + (c_r²+c_rc_t+c_t²)/(3(c_r+c_t))` |
| Lift slope | DATCOM (incompressible) `CLα = 2πAR/(2+√(AR²/κ²·(1+tan²Λ½c)+4))`, `κ = cℓα/2π` from the selected 2D section (script used fixed 4.5 – overridable in calibration). The sweep term also accounts for the section working normal to the sweep line |
| Oswald e | per pair; blank = environment value |
| Hull volume | `Vol = C_vol · πD²L/4` |
| Hull wetted area | `S_wet = C_wet · πDL` |
| Added mass / inertia | Lamb's k₁, k₂, k′ for a prolate spheroid of the same L/D |

### 2D sections (`data/sections.json`)

Each fin pair has its own section. A section stores data points
(cℓα, cd_min, cℓmax) at chord Reynolds numbers; the model evaluates them at
`Re = V_cruise·MAC/ν` (cruise speed, fixed over the speed scan), linear in
log10 Re and clamped to the data range. The pair then uses

```
CLα   = DATCOM(AR, Λ½c, cℓα)
CD0   = cd_min                    (environment "Fin CD0" if null – the ideal 2π section)
CLmax = 0.9·cℓmax·cosΛ¼c          (Raymer; conservative for low-AR fins)
```

Not modelled: compressibility (Prandtl–Glauert, matters above M ≈ 0.3 in
air), cd rise with cℓ away from the drag bucket, cavitation.

### Delta planform (sharp leading edge) – Polhamus suction analogy

Selected per pair ("Planform model: Delta"). NASA TN D-3767 (1966):

```
CL(α) = Kp·sinα·cos²α + Kv·cosα·sin²α      Kp = DATCOM CLα, Kv = model.json delta.Kv (default π)
CDi   = CL·tanα                            (no leading-edge suction)
CLmax = CL(α_stall)                        α_stall = model.json delta.stallAlphaDeg (default 25°)
```

Trim α is found by inverting CL(α). For stability and gusts the pair is
linearised about trim: slope `dCL/dα`, and the perturbation a.c. is the
slope-weighted mix of ¼-MAC (potential lift) and the planform centroid
(vortex lift, ⅔c_r for a delta). Gust drag uses the curvature of CDi(α):
`ΔCD = ½·CDi''(α0)·σα²` (identical to the parabolic-polar term for a
conventional pair). Polhamus is valid for sharp, highly swept leading edges
(Λ ≳ 45–50°) and overpredicts lift after vortex breakdown; a warning is shown
below 45° sweep.

### Downwash on the rear pair

DATCOM / Raymer low-speed downwash gradient from the front pair:

```
dε/dα = 4.44·[K_A·K_λ·K_H·√cosΛ¼c]^1.19 · (dCL/dα / CLα)_front
K_A = 1/AR − 1/(1 + AR^1.7),  K_λ = (10 − 3λ)/7,  K_H = (1 − |h_H/b|)/(2l_H/b)^(1/3)
```

l_H = distance between the a.c.s, h_H = rear-pair height above the front-pair
plane (input), b = front span. The last factor scales the wake for delta
vortex lift. Capped at 0.95 (warning) when the tail sits too close behind the
wing. Applied in the vertical plane only (no sidewash in yaw); can be switched
off per design.
The bundled coefficients are approximate literature-typical values; sections
marked `rough` raise a warning. Verify against the cited source before using a
result as justification.

## 3. Drag build-up

```
D = D_hull + D_fin,profile + D_trim + D_gust,vert + D_gust,lat + D_extra
D_hull        = q_p·S_wet·Cf(Re)·FF            Cf = 0.075/(log10 Re − 2)²   (ITTC-57)
                                               FF = 1 + 1.5/λ^1.5 + 7/λ³    (Hoerner)
D_fin,profile = q_p·ΣCD0_i·S_i·(2 if cruciform)
D_trim        = Σ q·S_i·CDi_i(α_i)          conventional: CL²/(π·e_i·AR_i); delta: CL·tanα
D_gust        = Σ q·S_i·½·CDi_i''·σ²_α,i     conventional: CLα²σ²/(π·e·AR) ← the script's turbulent drag growth
```

`E[CL²] = CL_trim² + CLα²σ_α²` is exactly the mean of the script's
`CL = CLtrim + CLα·α_g` squared, so the expected drag is obtained without a
time series. `test/model.test.js` verifies against a Monte-Carlo run of the
script (agreement < 1 %).

`q_p = q` unless "dynamic-pressure penalty" is enabled, then
`q_p = ½ρ(V0² + σ²_u,eff + σ²_v,eff + σ²_w,eff)` (hull-length filtered).
The script ignores this effect; it is off by default.

## 4. Turbulence – frequency domain

One-sided von Kármán spectra in spatial frequency Ω [rad/m] (MIL-HDBK-1797),
each integrating to σ²:

```
Φ_u(Ω) = σ_u²·(2L_u/π) / (1 + (1.339·L_u·Ω)²)^(5/6)
Φ_w(Ω) = σ_w²·(L_w/π)·(1 + 8/3·(1.339·L_w·Ω)²) / (1 + (1.339·L_w·Ω)²)^(11/6)   (same form for v)
```

Frozen turbulence: temporal frequency ω = Ω·V0. Each lifting element sees the
gust at its own position, averaged over its chord:

```
w_i(Ω) = sinc(Ω·c_i/2) · e^(−iΩ·x_i) · w(Ω)        (hull: sinc(ΩL/2)·e^(−iΩL/2))
σ²_α,i = ∫ |G_α,i(Ω)|² Φ_w(Ω) dΩ
```

*Fixed-attitude mode* (`G_α,i = w_i/V`) reproduces the script's assumption.
*Vehicle-response mode* (default) includes the rigid-body response below.
Integration: 600-point log grid, Ω = 10⁻⁵ … 10³ rad/m, trapezoid in ln Ω.

## 5. Trajectory stability (no steering)

Linear pitch/heave model (yaw/sway identical with v-gusts when cruciform),
states α and q, flight-path rate γ̇:

```
K_i = q·S_i·(dCL/dα)_i,   l_i = x_ac,i − x_CG (positive aft)
α_i = α + q·l_i/V + w_i/V                        (front)
α_r = α + q·l_r/V + w_r/V − E·e^{−iωτ}·α_f        (rear, downwash E = dε/dα, τ = l_H/V)
ΔL  = Σ K_i α_i + F_FK
ΔM  = −Σ K_i l_i α_i + M_Munk(α + w_h/V) + M_q,extra·q − F_FK·(L/2 − x_CG)
(m + k₂ρVol)·V·γ̇ = ΔL,   α̇ = q − γ̇,   (I_yy + I_added)·q̇ = ΔM
M_Munk = (k₂ − k₁)·ρ·Vol·V²                (destabilising hull moment)
F_FK   = (1 + k₂)·ρ·Vol·∂w_hull/∂t         (Froude–Krylov + added mass; negligible in air, first-order in water)
```

Outputs:

* **Neutral point / static margin** – `x_NP = x_CG − M_α/L_α` with the rear pair's K reduced by (1 − dε/dα).
  Eigenvalues use the quasi-steady limit τ → 0 (downwash-lag damping omitted – conservative);
  the gust transfer functions keep the lag.
* **Short-period eigenvalues** – stability, ω_n, ζ.
* **RMS gust incidence, load factor, pitch rate** – spectral integrals of the transfer functions.
* **Track wander** – RMS deviation from the mean straight course after distance X (T = X/V):
  `Var = V² ∫ |G_γ̇|² Φ · 4 sin²(ωT/2)/ω⁴ dΩ` (grows ∝ √X at long range).
* **Path-angle wander** σ_γ – note the classic "rotary gust" effect: in long-wave gusts
  the nose meets the gust before the tail, so the vehicle pitches with the gust gradient
  and its path angle can exceed σ_w/V.

Pitch inertia defaults to a uniform solid cylinder `m0(L²/12 + D²/16)` unless given.

## 6. Waterjet and range

```
T = D,  Vj = (V + √(V² + 4T/(ρA_n)))/2
P_jet    = ½ρA_n·Vj·(Vj² − V²),   η_F = 2V/(V + Vj)
P_engine = P_jet / η_jet
R = η_jet·η_F/(BSFC·g) · (W/D) · ln(m0/m_end),   m_end = m0 − (1 − reserve)·m_fuel
```

Equivalent to the script's `R = V/c · L/D · ln(m0/mf)` with `c = BSFC·g·V/(η_jet·η_F)`;
the app displays this equivalent c (1/h) for comparison with the script's 2.0 1/h.
Max speed = highest scanned speed where `P_engine ≤ P_max`.

## 7. Known simplifications (flag list – open architecture)

| # | Simplification | Where to change |
|---|---|---|
| 1 | **Turbulence presets are atmospheric placeholders** (σ_u,v,w = 5/5/3 m/s, L = 200/200/50 m) | `data/model.json → turbulence` (set `placeholder: false` once replaced) |
| 2 | Linear (small-angle) gust response; warning above `limits.maxGustAlphaDeg` | `model.js planeSystem` |
| 3 | CG fixed during fuel burn; trim/L-D evaluated at m0 | `model.js analyze` |
| 4 | Hull lift neglected (Munk moment kept); fin–hull interference ignored | `model.js pointAt` |
| 5 | Same fin pairs used in yaw when cruciform | `design.cruciform` |
| 6 | "Wind" on the map is the original tool's atmospheric placeholder (future: currents) | `model.js envelope` |
| 7 | Pitch damping only from fins + optional extra term | `extraPitchDamping_Nms` |
| 8 | ITTC-57 + Hoerner form factor for hull drag | calibration panel / `model.js` |
| 9 | No compressibility; trapezoidal planforms only (no cranked/double delta); no wing–body carry-over | `model.js planform` |
| 10 | Delta: Polhamus Kv and stall α are single global values | `data/model.json → delta` |
| 11 | Downwash: DATCOM empirical formula, vertical plane only, no downwash lag in the eigenvalues | `model.js downwashGradient` |
| 12 | No separate vertical tail – "cruciform" reuses the same pairs in yaw | `design.cruciform` |
