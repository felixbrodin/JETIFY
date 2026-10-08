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

### Tail configuration (`tailType`)

| Option | Pitch plane | Yaw plane |
|---|---|---|
| `traditional` | wing + horizontal tail (rear pair) | vertical fin |
| `vtail` | wing + rear pair, K·cos²Γ | rear pair, K·sin²Γ |
| `tailless` | wing only – elevon trim assumed | vertical fin (optional) |
| `cruciform` (AUV) | front + rear pair | front + rear pair (second identical pair) |

* **Vertical fin** (single surface, height h, area S, AR_geo = h²/S): geometry from its
  mirror-image pair; lift slope from `AR_eff = factor·AR_geo`
  (`model.json → verticalTail.effectiveARFactor`, default 1.55 – endplate effect, Raymer).
  Fin area 0 = no fin (yaw then unstable, warning).
* **V-tail** (Purser–Campbell): each panel's normal force has a vertical share cos Γ and a
  side share sin Γ, and a vertical/lateral gust gives panel incidence ·cos Γ / ·sin Γ, so
  `K_pitch = q·S·CLα·cos²Γ`, `K_yaw = q·S·CLα·sin²Γ`. Panel CL at trim = `L_rear/(q·S·cos Γ)`
  (used for induced drag and the stall check). S, AR and span are measured along the panels.
* **Tailless**: the wing carries W; the pitching moment is assumed trimmed by elevons
  (reflex) and **elevon trim drag is not modelled**. Static margin = neutral point (wing a.c.
  with the hull Munk moment) − CG.
* Old profiles: `cruciform: true` → `cruciform`; `cruciform: false` → `traditional` with fin area 0.

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
off per design. The formula assumes the tail span lies inside the wing's vortex
span – a warning is shown when the rear pair is wider than the front pair (typical
AUV fins: there the outer tail sees upwash and dε/dα is overestimated; the bundled
AUV profiles therefore have downwash off).
The bundled coefficients are approximate literature-typical values; sections
marked `rough` raise a warning. Verify against the cited source before using a
result as justification.

## 3. Drag build-up

```
D = D_hull + D_fin,profile + D_trim + D_gust,vert + D_gust,lat + D_wave + D_extra
D_hull        = q_p·S_wet·Cf(Re)·FF·(1 + 0.144M²)^−0.65
                                               Cf = 0.075/(log10 Re − 2)²   (ITTC-57)
                                               FF = 1 + 1.5/λ^1.5 + 7/λ³    (Hoerner)
                                               last factor: compressible turbulent friction (Raymer)
D_wave        = Σ q·S_i·CD_wave,i            Korn/Lock, see "Compressibility" below
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

### Atmosphere and compressibility

* **ISA** (default, `useISA`): from the altitude input, troposphere `T = 288.15 − 0.0065h`,
  isothermal 216.65 K above 11 km (to 20 km); ρ from the ideal-gas law, ν from Sutherland,
  `a = √(1.4·R·T)`. ISA off = manual ρ, ν and speed of sound (water: a ≈ 1480 m/s → M ≈ 0).
* **Prandtl–Glauert** in the DATCOM slope, β = √(1 − M²) (M capped at 0.95):
  `CLα = 2πAR / (2 + √(AR²β²/κ²·(1 + tan²Λ½c/β²) + 4))`, κ = a0/2π. Applies to every pair and
  the fin (and so to the delta's Kp); the downwash gradient is scaled by CLα,M/CLα,0 (DATCOM).
* **Drag divergence** per surface (Korn equation, Mason's form) and **wave drag** (Lock):
  `M_dd = κ_A/cosΛ¼c − (t/c)/cos²Λ¼c − CL/(10cos³Λ¼c)`, `M_crit = M_dd − (0.1/80)^{1/3}`,
  `CD_wave = 20(M − M_crit)⁴` for M > M_crit. t/c and κ_A from `data/sections.json`
  (κ_A 0.87 conventional, 0.95 supercritical; ideal 2π section: t/c from model.json).
* Section Re and Mach are evaluated at each speed of the speed scan (re-resolved per speed).
* **Engine power lapse**: `P_avail = P_max·σ^n`, σ = ρ/ρ0, n = `model.json → engine.powerLapseExponent`
  (0.7, turbine-like). Manual (ISA off): no lapse.

### Flight envelope (altitude vs speed chart)

Calm air, take-off mass, ISA. For each speed (`atmosphere.envelopeSpeedSteps`) the altitudes
0–20 km (`envelopeAltitudeSteps`, edges refined by bisection) are checked for steady level
flight: trim L = W without stall (CLmax, delta stall α), engine power for the calm-air drag
≤ P_max·σ^n, and M ≤ 0.95. Plotted: ceiling and floor of the flyable band, plus the
stall-limited and power/Mach-limited ceilings separately.

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

### Lateral-directional: sideslip, roll, yaw, bank

States `[β, p, r, φ]` (stability axes, I_xz neglected):

```
m'V(β̇ + r) = Y + W·φ          I_x·ṗ = L − K_φ·φ          I_z'·ṙ = N          φ̇ = p
```

* **Side-force surfaces** (vertical fin, V-tail yaw share, cruciform pairs) with slope
  `K_i = q·S_i·CLα_i·proj²`, arm l_i (aft of CG) and height z_i of the a.c. above the roll
  axis (fin: D/2 + ȳ_MAC; V-tail: (D/2 + ȳ_MAC)·sinΓ; cruciform: 0) feel
  `β_i = β − r·l_i/V + p·z_i/V (+ gust)` and give `Y_i = −K_i·β_i`, `L_i = Y_i·z_i`, `N_i = −Y_i·l_i`.
* **Horizontal pairs** (strip theory / Nelson, per pair, slope a at trim):
  `C_lβ = −a·Γ·(1+2λ)/(6(1+λ)) − CL·tanΛ¼c·(1+2λ)/(3(1+λ))` (dihedral + sweep),
  `C_lp = −a(1+3λ)/(12(1+λ))`, `C_lr = CL/4`, `C_np = −CL/8`, `C_nr = −CD0/4`.
  Dihedral Γ is an input for the wing (front pair); tail/V-tail panels add roll damping.
* **Hull**: Munk yaw moment `N_β −= (k₂−k₁)ρVol·V²`, Froude–Krylov side force on the hull.
* **Inertia** (blank = estimate): `I_x = m0(D²/8 + (R̄x·b/2)²)`, `I_z = I_y + m0(R̄x·b/2)²`,
  R̄x = `model.json → lateral.rollGyrationFactor` (0.25). Optional roll stiffness K_φ
  (e.g. hydrostatic W·BG for an AUV – buoyancy itself stays outside the model).
* **Roll held level** (`rollFree: false`, used by the bundled AUV profiles): only β and r
  remain – the previous yaw/sway model.
* **Modes**: eigenvalues of the 4×4 matrix (characteristic polynomial + Durand–Kerner).
  Complex pair = Dutch roll; real roots: most negative = roll subsidence (τ = −1/λ),
  smallest |λ| = spiral (T½ or T₂ = ln2/|λ|). Static: `C_nβ` (> 0) and `C_lβ` (< 0),
  referenced to the wing (q·S·b).
* **Gust inputs** (independent, variances add): lateral gust v (delay + chord averaging per
  surface, as in pitch) and the MIL-F-8785C rolling gust
  `Φ_p(Ω) = σ_w²/L_w · 0.8(πL_w/4b)^{1/3} / (1 + (4bΩ/π)²)` acting on the horizontal pairs.
* **Lateral track wander**: course rate `χ̇ = (Y + W·φ)/(m'V)`. A steady rolling moment gives
  a steady turn (heading random walk), so the rolling-gust part is measured from the
  initial course: `Var = V²∫|G_χ̇|²Φ·|(e^{iωT}−1−iωT)/(iω)|²/ω² dΩ`. A warning flags wander
  above 30 % of the distance (beyond small-angle theory – the uncontrolled vehicle does not
  hold course) and RMS bank above 30°.

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
| 5 | Cruciform reuses the same pairs in yaw | `design.tailType` |
| 5b | Lateral: strip-theory derivatives; no I_xz, no wing vertical position (high/low wing) term, no sidewash, rolling gust on horizontal pairs only; linear (small bank) | `model.js lateralSystem / pairRollDerivs` |
| 6 | "Wind" on the map is the original tool's atmospheric placeholder (future: currents) | `model.js envelope` |
| 7 | Pitch damping only from fins + optional extra term | `extraPitchDamping_Nms` |
| 8 | ITTC-57 + Hoerner form factor for hull drag | calibration panel / `model.js` |
| 9 | Compressibility: Prandtl–Glauert + Korn/Lock only (no shocks, no transonic CLmax/buffet, no hull/body wave drag, invalid above M ≈ 0.95); trapezoidal planforms only (no cranked/double delta); no wing–body carry-over | `model.js planform` |
| 10 | Delta: Polhamus Kv and stall α are single global values | `data/model.json → delta` |
| 11 | Downwash: DATCOM empirical formula, vertical plane only, no downwash lag in the eigenvalues | `model.js downwashGradient` |
| 12 | Tailless: elevon trim drag and Cm0 (reflex) not modelled; V-tail: no panel interference | `model.js trim` |
