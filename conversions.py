"""Convert between classical orbital elements and a state vector (r, v) in ECI.

coe_to_rv: elements -> r, v        rv_to_coe: r, v -> elements

Pipeline (main):  read (orbits.csv -> elements)  ->  process (elements -> r, v)  ->  output
"""

import csv
from pathlib import Path

import numpy as np

from describe import read_orbits
from vectors import MU_KM3_S2, angle_deg, invariants

# --- Constants ---------------------------------------------------------------

I_HAT = np.array([1.0, 0.0, 0.0])   # ECI X axis (vernal equinox direction)

BASE_DIR = Path(__file__).resolve().parent
INPUT_PATH = BASE_DIR / "data" / "orbits.csv"
OUTPUT_DIR = BASE_DIR / "output"
STATES_PATH = OUTPUT_DIR / "states_from_elements.csv"
SKIPPED_PATH = OUTPUT_DIR / "elements_skipped.txt"

EPOCH_UTC = "2026-10-21T00:00:00Z"   # orbits.csv has no epoch, so one fixed epoch is used

R_FIELDS = ("x_km", "y_km", "z_km")                  # position columns
V_FIELDS = ("vx_km_s", "vy_km_s", "vz_km_s")         # velocity columns
STATE_FIELDS = ("name", "epoch_utc") + R_FIELDS + V_FIELDS   # same columns as data/states.csv



# --- Helpers -----------------------------------------------------------------

def Rz(theta):
    """Rotate a vector by theta (rad) about the Z axis."""
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c,  -s,  0.0],
                     [s,   c,  0.0],
                     [0.0, 0.0, 1.0]])


def Rx(theta):
    """Rotate a vector by theta (rad) about the X axis."""
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[1.0, 0.0, 0.0],
                     [0.0,  c,  -s],
                     [0.0,  s,   c]])


def rotation_matrix(raan_rad, i_rad, argp_rad):
    """Perifocal-to-ECI rotation Q = Rz(raan) @ Rx(i) @ Rz(argp).

    Read right to left: turn by argp, then tilt by i, then swing by raan.
    """
    return Rz(raan_rad) @ Rx(i_rad) @ Rz(argp_rad)


def angle_360(a, b, sign):
    """Angle from a to b, 0-360 deg. sign < 0 means the second half."""
    ang = angle_deg(a, b)
    return ang if sign >= 0 else 360.0 - ang


# --- 1. Read -----------------------------------------------------------------
# Elements are read by read_orbits() from describe.py (imported above), which
# validates each row (0 <= e < 1, periapsis above the Earth, 0 <= i <= 180)
# and returns (orbits, skipped).


# --- 2. Process --------------------------------------------------------------

def perifocal(a_km, e, nu_rad):
    """Position (km) and velocity (km/s) in the perifocal frame."""
    p = a_km * (1 - e**2)                   # semi-latus rectum
    r = p / (1 + e * np.cos(nu_rad))        # orbit equation
    k = np.sqrt(MU_KM3_S2 / p)              # speed factor (= mu / h)

    r_pqw = r * np.array([np.cos(nu_rad), np.sin(nu_rad), 0.0])
    v_pqw = k * np.array([-np.sin(nu_rad), e + np.cos(nu_rad), 0.0])

    return r_pqw, v_pqw


def coe_to_rv(a_km, e, i_deg, raan_deg, argp_deg, nu_deg):
    """Classical orbital elements -> position (km) and velocity (km/s) in ECI.

    Angles in degrees. Raises ValueError if a <= 0 or e is outside 0 <= e < 1.
    """
    if a_km <= 0:
        raise ValueError(f"semi-major axis must be positive, got {a_km}")
    if not (0 <= e < 1):
        raise ValueError(f"eccentricity must be 0 <= e < 1, got {e}")

    i, raan, argp, nu = np.radians([i_deg, raan_deg, argp_deg, nu_deg])

    r_pqw, v_pqw = perifocal(a_km, e, nu)
    Q = rotation_matrix(raan, i, argp)

    return Q @ r_pqw, Q @ v_pqw


def rv_to_coe(r, v):
    """State vector (km, km/s) in ECI -> classical orbital elements (angles in deg).

    Singular cases (undefined angles set to 0, replacement angle stored instead):
      circular, inclined     : argp = 0, nu = u (argument of latitude, n -> r)
      equatorial, elliptical : raan = 0, argp = longitude of periapsis (I -> e)
      circular, equatorial   : raan = argp = 0, nu = true longitude (I -> r)
    u_deg is returned for every inclined orbit (None when equatorial).
    Raises ValueError (from invariants) if the state is not a closed orbit.
    """
    inv = invariants(r, v)                  # validates; h, n, e_vec, a, e, i, flags
    n = inv["n_vec"]
    e_vec = inv["e_vec"]
    circular = inv["circular"]
    equatorial = inv["equatorial"]

    # Argument of latitude: angle n -> r, defined whenever the orbit has a node
    u = None if equatorial else angle_360(n, r, r[2])         # test r_z

    if not circular and not equatorial:                       # general orbit
        raan = angle_360(I_HAT, n, n[1])                      # Omega: test n_y
        argp = angle_360(n, e_vec, e_vec[2])                  # omega: test e_z
        nu = angle_360(e_vec, r, inv["r_dot_v"])              # nu:    test r.v
    elif circular and not equatorial:                         # no periapsis
        raan = angle_360(I_HAT, n, n[1])
        argp = 0.0
        nu = u
    elif equatorial and not circular:                         # no node
        raan = 0.0
        argp = angle_360(I_HAT, e_vec, e_vec[1])              # longitude of periapsis: test e_y
        nu = angle_360(e_vec, r, inv["r_dot_v"])
    else:                                                     # no node, no periapsis
        raan = 0.0
        argp = 0.0
        nu = angle_360(I_HAT, r, r[1])                        # true longitude: test r_y

    return {
        "a_km": inv["a_km"],
        "e": inv["e"],
        "i_deg": inv["i_deg"],
        "raan_deg": raan,
        "argp_deg": argp,
        "nu_deg": nu,
        "u_deg": u,
        "circular": circular,
        "equatorial": equatorial,
    }


# --- 3. Output ---------------------------------------------------------------

TABLE_HEADER = f"{'Name':<20}{'|r| km':>10}{'|v| km/s':>10}   {'r_eci (km)':<30}{'v_eci (km/s)'}"


def fmt_vec(vec, decimals):
    """'[ 1.0, -2.0, 3.0]' with a fixed number of decimals."""
    return "[" + ", ".join(f"{x:.{decimals}f}" for x in vec) + "]"


def format_row(name, r, v):
    """One console table row for one converted orbit."""
    return (f"{name:<20}{np.linalg.norm(r):>10.1f}{np.linalg.norm(v):>10.3f}   "
            f"{fmt_vec(r, 1):<30}{fmt_vec(v, 3)}")


def format_skipped(entry):
    """One line for a rejected row, used by both the console and the skipped file."""
    line_number, name, reason = entry
    return f"line {line_number}  {name}: {reason}"


def to_csv_row(name, r, v):
    """One CSV row in the states.csv format."""
    row = {"name": name, "epoch_utc": EPOCH_UTC}
    row.update(zip(R_FIELDS, r))
    row.update(zip(V_FIELDS, v))
    return row


def write_states(results, path):
    """Write one row per converted orbit, in the same format as data/states.csv."""
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=STATE_FIELDS)
        writer.writeheader()
        writer.writerows(to_csv_row(name, r, v) for name, r, v in results)


def write_skipped(skipped, path):
    """Write one line per rejected row, or say that none were rejected."""
    lines = [format_skipped(entry) for entry in skipped] or ["No rows skipped."]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


# --- Main --------------------------------------------------------------------

def main():
    # 1. Read
    orbits, skipped = read_orbits(INPUT_PATH)

    # 2. Process
    results = []
    for el in orbits:
        r, v = coe_to_rv(el.a_km, el.e, el.i_deg, el.raan_deg, el.argp_deg, el.nu_deg)
        results.append((el.name, r, v))

    # 3. Output
    print(f"Read {len(results) + len(skipped)} rows from {INPUT_PATH.name}: "
          f"{len(results)} converted, {len(skipped)} skipped\n")
    print(TABLE_HEADER)
    for name, r, v in results:
        print(format_row(name, r, v))
    if skipped:
        print("\nSkipped rows:")
        for entry in skipped:
            print("  " + format_skipped(entry))

    OUTPUT_DIR.mkdir(exist_ok=True)
    write_states(results, STATES_PATH)
    write_skipped(skipped, SKIPPED_PATH)
    print(f"\nWrote {len(results)} state vectors to {STATES_PATH}")
    print(f"Wrote {len(skipped)} skipped rows to {SKIPPED_PATH}")


if __name__ == "__main__":
    main()