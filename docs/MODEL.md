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

With one pair at CG this reduces to the script's `CLtrim = W/(½ρV0²S)`.
Trim is evaluated at take-off mass m0 (as in the script).

## 2. Geometry

| Item | Relation |
|---|---|
| Pair span / chord | `b = √(S·AR)`, `c = √(S/AR)` (S = both fins of the pair) |
| Lift slope | Helmbold `CLα = 2πAR/(2+√(AR²+4))` (script used fixed 4.5 – overridable in calibration) |
| Hull volume | `Vol = C_vol · πD²L/4` |
| Hull wetted area | `S_wet = C_wet · πDL` |
| Added mass / inertia | Lamb's k₁, k₂, k′ for a prolate spheroid of the same L/D |

## 3. Drag build-up

```
D = D_hull + D_fin,profile + D_trim + D_gust,vert + D_gust,lat + D_extra
D_hull        = q_p·S_wet·Cf(Re)·FF            Cf = 0.075/(log10 Re − 2)²   (ITTC-57)
                                               FF = 1 + 1.5/λ^1.5 + 7/λ³    (Hoerner)
D_fin,profile = q_p·CD0·ΣS_i·(2 if cruciform)
D_trim        = Σ q·S_i·CL_i²/(π·e·AR_i)
D_gust        = Σ q·S_i·CLα_i²·σ²_α,i/(π·e·AR_i)        ← the script's turbulent drag growth
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
K_i = q·S_i·CLα_i,   l_i = x_i − x_CG (positive aft)
ΔL  = Σ K_i (α + q·l_i/V + w_i/V) + F_FK
ΔM  = −Σ K_i l_i (α + q·l_i/V + w_i/V) + M_Munk(α + w_h/V) + M_q,extra·q − F_FK·(L/2 − x_CG)
(m + k₂ρVol)·V·γ̇ = ΔL,   α̇ = q − γ̇,   (I_yy + I_added)·q̇ = ΔM
M_Munk = (k₂ − k₁)·ρ·Vol·V²                (destabilising hull moment)
F_FK   = (1 + k₂)·ρ·Vol·∂w_hull/∂t         (Froude–Krylov + added mass; negligible in air, first-order in water)
```

Outputs:

* **Neutral point / static margin** – `x_NP = x_CG + (Σ K_i l_i − M_Munk)/Σ K_i` (speed-independent).
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
