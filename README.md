# Jetify

A design tool for the overall dimensions of a **waterjet-propelled, jet-fuel AUV**.
It predicts **one-way range** and **trajectory stability** (how steadily the vehicle
holds its course in turbulence) for mission planning. There is no time simulation.
The turbulence penalty comes from integrating the **von Kármán spectrum in the
frequency domain**.

Jetify grew out of the *drönarräckvidd* range tool and reuses its app structure,
data import/export, map, scenarios and calibration. The rotor/battery physics was
replaced with the model in [`docs/MODEL.md`](docs/MODEL.md).

## What you can design

| Box | Variables (sliders) | Presets (editable) |
|---|---|---|
| Hull | length, diameter, centre of gravity | volume and wetted-area coefficients, pitch inertia, extra pitch damping |
| Control surfaces | front and rear pair: area S, aspect ratio AR, centre-of-lift position x | cruciform on/off |
| Waterjet, fuel & mass | max engine power | η_jet (engine → jet-stream power), nozzle diameter, BSFC, empty mass, fuel mass, payload |
| Environment | water density ρ, cruise speed V0, fuel reserve, wind (placeholder) | viscosity, Oswald e, fin CD0, CLmax |
| Turbulence | — | σ_u, σ_v, σ_w, L_u, L_v, L_w (**atmospheric placeholders**) |

## Outputs

- One-way range (Breguet with fuel burn), endurance, best-range speed, max speed
- Drag build-up (hull, fin profile, trim-induced, gust-induced) and required engine power
- Static margin and neutral point (fins plus the hull's destabilising Munk moment)
- Pitch/yaw mode stability, RMS gust incidence, load factor and pitch rate
- Track wander (depth and lateral) vs distance with no steering
- Side view, map reach, and charts: range vs speed, wander vs distance, static margin vs rear-fin position, gust spectrum

## Adjusting slider spans

Every adjustable variable (label, unit, slider min/max, step and decimals) is defined in
**`data/variables.json`**. Edit the file and reload the page; no code changes are needed.
The file starts with a short `_readme` explaining each field. Run the tests afterwards to
catch typos and spans that exclude a default value.

| File | Holds |
|---|---|
| `data/variables.json` | spans of all sliders and inputs (+ payload slider, speed-scan step) |
| `data/model.json` | physics coefficients, environment and turbulence default values |
| `data/vehicles.json` | vehicle profiles (default design values) |

## Run

The page loads its data with `fetch`, so it needs a local web server:

```
start-server.bat                 (Windows; serves on http://localhost:8000)
python -m http.server 8000       (alternative)
```

Tests (Node ≥ 18, no dependencies):

```
node --test test/model.test.js
```

The tests check the spectra (they integrate to σ²), the Lamb added-mass table and the
waterjet relations, and compare the closed-form gust drag penalty against a Monte-Carlo
run of the original MATLAB script.

## ⚠ Before trusting the turbulence numbers

The turbulence presets in `data/model.json` are still the **atmospheric** values
from the MATLAB reference script. Replace them with measured underwater values
for the operating area (e.g. the Baltic), then set `"placeholder": false`.

## Layout

```
index.html, styles.css
js/model.js     physics engine (no DOM; also loadable in Node)
js/data.js      profile loading, CSV/TSV/JSON import & export
js/app.js       UI – input fields are generated from data/variables.json
data/           variables.json (slider spans), model.json (coefficients), vehicles.*, payloads.*, map.json
docs/MODEL.md   derivation and list of simplifications
test/           node:test suite
```
