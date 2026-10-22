# orbit-state-vector-toolkit

Python tools for describing Earth orbits from their classical orbital elements, and for analysing an orbit from its state vector (position and velocity).

| Part | Tool | Input | Status |
|---|---|---|---|
| 1 | Orbit Describer (`describe.py`) | classical orbital elements | complete |
| 2 | Vector invariants (`vectors.py`) | state vectors (r, v) | complete, 23 tests |
| 3 | Elements ↔ state vector conversions | both | next |
| 4 | Kepler solver and elements explorer | | planned |

## Repository structure

```
orbit-state-vector-toolkit/
├── describe.py              Part 1: orbit describer
├── vectors.py               Part 2: vector invariants from (r, v)
├── pytest.ini               pytest configuration
├── data/
│   ├── orbits.csv           Part 1 input: 6 real orbit types + 2 invalid rows
│   └── states.csv           Part 2 input: 3 test states + 2 invalid rows
├── output/
│   ├── descriptions.csv     Part 1 results
│   ├── skipped.txt          Part 1 rejected rows
│   ├── invariants.csv       Part 2 results
│   └── states_skipped.txt   Part 2 rejected rows
├── tests/
│   └── test_vectors.py      23 tests for Part 2
```

## Setup and run

Requires Python 3 with NumPy (Part 2) and pytest (tests). Part 1 uses the standard library only.

```
git clone git@github.com:dkumar2016/orbit-state-vector-toolkit.git
cd orbit-state-vector-toolkit

python describe.py     # Part 1: describe the orbits in data/orbits.csv
python vectors.py      # Part 2: invariants of the states in data/states.csv
pytest -v              # run the tests
```

Both programs follow the same pipeline: **read** a CSV, **process** each valid row once, **output** to the console and to files. Bad rows are skipped with their line number and a reason; they never stop the program.

## Conventions

| | |
|---|---|
| Frame | Earth-centred inertial (ECI): I towards the vernal equinox (♈), K towards the north pole, J completing the right-handed set |
| Units | km, km/s, degrees (radians only inside calculations) |
| Constants | μ⊕ = 398,600.4418 km³/s², R⊕ = 6,378.137 km |
| Model | two-body (Keplerian): no J2, drag or other perturbations |

---

## Part 1: Orbit Describer

`describe.py` reads classical orbital elements, rejects physically impossible orbits with a reason, and describes each valid orbit in one plain-English sentence:

```
Molniya: highly elliptical HEO, altitude 528-39,836 km, prograde at 63.4 deg, period 11.97 h (718.0 min). Ground track reaches +/-63.4 deg latitude.
```

### Input: `data/orbits.csv`

| Column | Element | Units |
|---|---|---|
| `name` | orbit label | text |
| `a_km` | semi-major axis *a* | km |
| `e` | eccentricity *e* | – |
| `i_deg` | inclination *i* | deg |
| `raan_deg` | right ascension of the ascending node *Ω* | deg |
| `argp_deg` | argument of perigee *ω* | deg |
| `nu_deg` | true anomaly *ν* | deg |

Six real orbit types (ISS, GPS, GEO, Molniya, sun-synchronous, retrograde LEO) and two deliberately invalid rows. Values are representative of each orbit type, not current tracking data.

### Pipeline

```
read_orbits      CSV row  -> Elements      (validated on creation, or skipped with a reason)
summarize        Elements -> summary dict  (every label and derived value, computed once)
describe / write summary  -> sentence, descriptions.csv
```

The sentence and the CSV row are both built from the same summary, so they can't disagree.

### Validation

| Check | Why |
|---|---|
| 0 ≤ e < 1 | e ≥ 1 is a parabolic or hyperbolic trajectory, not a closed orbit |
| periapsis radius a(1 − e) > R⊕ | the orbit would pass through the Earth |
| 0 ≤ i ≤ 180° | inclination is only defined on this range |

Ω, ω and ν are wrapped into [0°, 360°). Non-numeric values and missing columns are also skipped with a reason.

### Classification rules

| Label | Rule |
|---|---|
| Shape | circular e < 0.001 · near-circular < 0.1 · elliptical < 0.5 · highly elliptical ≥ 0.5 |
| Regime | HEO if e ≥ 0.5; otherwise by mean altitude a − R⊕: LEO ≤ 2,000 km · GEO within ±200 km of 35,786 km · MEO below GEO · Beyond GEO |
| Inclination | equatorial < 5° · prograde < 85° · polar 85–95° · retrograde ≤ 175° · retrograde equatorial; LEO at 96–102° flagged as sun-synchronous band |
| Max latitude | i for prograde, 180° − i for retrograde |

HEO is checked first: Molniya's mean altitude falls in the MEO range, but its altitude swings from 528 to 39,836 km, so "MEO" would misdescribe it.

### Results

| Orbit | Shape | Regime | Inclination class | Perigee–apogee alt (km) | Period |
|---|---|---|---|---|---|
| ISS | circular | LEO | prograde | 396–403 | 92.6 min |
| GPS IIF | near-circular | MEO | prograde | 19,916–20,448 | 11.97 h |
| GOES-16 | circular | GEO | equatorial | 35,782–35,790 | 23.93 h |
| Molniya | highly elliptical | HEO | prograde | 528–39,836 | 11.97 h |
| Sentinel-2A | circular | LEO | retrograde (SSO band) | 785–787 | 100.6 min |
| Retrograde LEO | near-circular | LEO | retrograde | 431–569 | 94.6 min |

Skipped: line 8 (e = 1.2, not a closed orbit) and line 9 (periapsis 3,500 km, inside the Earth).

### Checks on the results

- **GEO period = one sidereal day:** GOES-16's 1436.1 min (23.93 h) is why it stays over one longitude.
- **Period depends only on a:** GPS and Molniya share a = 26,560 km and so the same 11.97 h period, despite very different shapes.
- **ISS:** period 92.6 min and altitude around 400 km match published values.

---

## Part 2: Vector invariants from a state vector

`vectors.py` takes a state vector — position **r** and velocity **v** in ECI at an epoch — and computes the vectors and scalars that define the orbit. These are the building blocks for converting a state vector into orbital elements (Part 3).

### What it computes

| Quantity | Formula | Gives |
|---|---|---|
| angular momentum **h** | **r** × **v** | the orbit plane; i = angle between **K** and **h** |
| node vector **n** | **K** × **h** | direction of the ascending node (Ω in Part 3) |
| eccentricity vector **e** | [(v² − μ/r) **r** − (**r·v**) **v**] / μ | points to periapsis; its length is e |
| specific energy ε | v²/2 − μ/r | a = −μ / (2ε) |
| **r · v** | dot product | climbing (> 0) or falling (< 0); the quadrant check for ν |

It also flags the two singular cases, where some elements are undefined:

| Flag | Condition | Element lost |
|---|---|---|
| `circular` | \|**e**\| < 10⁻¹⁰ | ω (no periapsis) |
| `equatorial` | \|**n**\| < 10⁻¹⁰ | Ω (no node) |

### Input: `data/states.csv`

```
name,epoch_utc,x_km,y_km,z_km,vx_km_s,vy_km_s,vz_km_s
```

Every state carries an **epoch**: a state vector is only valid at one instant, so a file without times would be incomplete.

| State | r (km) | v (km/s) | Purpose |
|---|---|---|---|
| ISS-like | [6778, 0, 0] | [0, 4.759, 6.013] | matches a hand calculation |
| Equatorial circular | [7000, 0, 0] | [0, √(μ/7000), 0] | both singular cases at once |
| Eccentric | [7000, 0, 0] | [0, 7.4511, 4.3019] | periapsis of an orbit with e = 0.3, i = 30°, built from vis-viva |
| Escape | [7000, 0, 0] | [0, 11.0, 0] | **invalid**: above escape speed, ε > 0 |
| Inside Earth | [6000, 0, 0] | [0, 8.0, 0] | **invalid**: \|r\| < R⊕ |

The test velocities are stored at full double precision. Rounding the circular case to six decimals gives e ≈ 8 × 10⁻⁸, which is a genuinely (slightly) elliptical orbit and is correctly *not* flagged circular.

### Pipeline

```
read_states      CSV row -> (line, name, epoch, r, v)   (bad numbers skipped)
invariants       r, v    -> dict                         (rejects non-orbits; all values computed once)
output           dict    -> console table, invariants.csv
```

Reading only converts numbers; `invariants` decides whether a state is a real orbit. Rejections from both stages are reported together, in file order.

### Results

```
Read 5 rows from states.csv: 3 processed, 2 skipped

Name                  |h| km2/s         e      a km   i deg  circular equatorial
ISS-like                51976.4  0.000063    6777.6   51.64        no         no
Equatorial circular     52822.4  0.000000    7000.0    0.00       yes        yes
Eccentric               60226.8  0.300000   10000.0   30.00        no         no

Skipped rows:
  line 5  Escape: not a closed orbit: energy 3.557 km2/s2 >= 0
  line 6  Inside Earth: position 6000.0 km is inside the Earth
```

| | ISS-like | Equatorial circular | Eccentric |
|---|---|---|---|
| **h** | [0, −40,756.1, 32,256.5] | [0, 0, 52,822.4] | [0, −30,113.4, 52,157.9] |
| **n** | [40,756.1, 0, 0] | [0, 0, 0] | [30,113.4, 0, 0] |
| **e** vector | [−6.3 × 10⁻⁵, 0, 0] | ~0 | [0.3, 0, 0] |
| ε (km²/s²) | −29.406 | −28.471 | −19.930 |

`output/invariants.csv` holds every value, with vectors split into x, y, z columns.

### What the results show

- **The hand calculation is reproduced exactly** (**h**, i = 51.64°, a = 6,777.6 km).
- **The ISS-like state is at apoapsis, not on a circle.** Its velocity, rounded to 3 decimals, is 0.24 m/s below circular speed. That makes e = 6.3 × 10⁻⁵, with **e** pointing *opposite* **r**: the satellite is at its highest point. The direction of a tiny eccentricity vector is set by rounding, which is why ω is unreliable for near-circular orbits and the argument of latitude (ω + ν) is used instead.
- **The eccentric state's e vector points along r**, confirming it sits at periapsis, as it was built.
- **The equatorial circular state loses both n and e**, so Ω and ω are both undefined; Part 3 must handle this.

### Tests

`pytest -v` runs **23 tests** in four groups:

| Group | Tests | What it checks |
|---|---|---|
| Known answers | 5 | ISS **h**, i and **n** against the hand calculation; eccentric e = 0.3, a = 10,000 km, i = 30°; **e** points to periapsis |
| Edge cases | 5 | circular and equatorial flags set (and *not* set for the ISS); escape and inside-Earth states rejected with the right message |
| Identities | 12 (4 × 3 states) | **h** ⊥ **r** and **v**; **e** lies in the orbit plane; \|h\| = √(μa(1 − e²)); Lagrange's identity h² + (r·v)² = r²v² |
| File reading | 1 | name, epoch and line number kept; a typo row skipped |

The identity tests hold for any orbit, so they check the maths itself, not just the cases written down. Writing them found one real bug: `is_zero` returned NumPy's `bool_`, which is not Python's `True` (`np.bool_(True) is True` is `False`). It now returns a plain `bool`.

---

## Limitations

- Two-body model only: no J2, drag, third bodies or solar pressure.
- Part 1's classification thresholds are conventions, not standards; values near a boundary can flip label.
- The sun-synchronous flag is an inclination band only. True sun-synchronism depends on a, e and i together (J2 nodal precession).
- The zero tolerance (10⁻¹⁰) separates rounding noise from real values. Deciding when a "nearly circular" orbit should be *treated* as circular is a mission choice, not set here.
- Epochs are carried through but not yet used; time scales and propagation come later.


## Coming next

- `coe_to_rv`: elements → state vector (perifocal frame and three rotations)
- `rv_to_coe`: state vector → elements, using **h**, **n** and **e**, with quadrant checks and singular-case handling; round-trip tests below 10⁻⁶
- Kepler's equation solver and an orbital-elements explorer notebook
