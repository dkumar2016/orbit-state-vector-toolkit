"""Vector invariants of an orbit from a state vector (r, v) in ECI.

Pipeline:  read (CSV -> r, v)  ->  process (r, v -> invariants)  ->  output
"""

import csv
from pathlib import Path

import numpy as np

# --- Constants ---------------------------------------------------------------

MU_KM3_S2 = 398600.4418             # Earth's gravitational parameter
R_EARTH_KM = 6378.137               # Earth's equatorial radius
K_HAT = np.array([0.0, 0.0, 1.0])   # ECI Z axis (north pole)
ZERO_TOL = 1e-10                    # vector length treated as zero (circular / equatorial)

BASE_DIR = Path(__file__).resolve().parent
INPUT_PATH = BASE_DIR / "data" / "states.csv"
OUTPUT_DIR = BASE_DIR / "output"
INVARIANTS_PATH = OUTPUT_DIR / "invariants.csv"
SKIPPED_PATH = OUTPUT_DIR / "states_skipped.txt"

R_FIELDS = ("x_km", "y_km", "z_km")                  # position columns in states.csv
V_FIELDS = ("vx_km_s", "vy_km_s", "vz_km_s")         # velocity columns in states.csv
RESULT_FIELDS = ("name", "epoch_utc",
                 "h_x", "h_y", "h_z", "n_x", "n_y", "n_z", "e_x", "e_y", "e_z",
                 "h", "e", "energy", "a_km", "i_deg", "r_dot_v",
                 "circular", "equatorial")            # columns in invariants.csv


# --- Helpers -----------------------------------------------------------------

def angle_deg(a, b):
    """Angle between two vectors, 0-180 deg (clipped against rounding)."""
    c = np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b))
    return np.degrees(np.arccos(np.clip(c, -1.0, 1.0)))


def is_zero(x, tol=ZERO_TOL):
    """True if a vector's length is below the tolerance."""
    return bool(np.linalg.norm(x) < tol)


# --- 1. Read -----------------------------------------------------------------

def read_states(path):
    """Read state vectors from a CSV file.

    Returns (states, skipped):
      states  - list of (line_number, name, epoch, r, v) for every row with valid numbers
      skipped - list of (line_number, name, reason) for every row that could not be read
    """
    states = []
    skipped = []
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for line_number, row in enumerate(reader, start=2):   # line 1 is the header
            name = row["name"]
            epoch = row["epoch_utc"]
            try:
                r = np.array([float(row[field]) for field in R_FIELDS])
                v = np.array([float(row[field]) for field in V_FIELDS])
                states.append((line_number, name, epoch, r, v))
            except (ValueError, TypeError) as err:
                skipped.append((line_number, name, str(err)))
    return states, skipped


# --- 2. Process --------------------------------------------------------------

def angular_momentum(r, v):
    """Specific angular momentum h = r x v (km^2/s), normal to the orbit plane."""
    return np.cross(r, v)


def node_vector(h):
    """Node vector n = K x h, pointing to the ascending node."""
    return np.cross(K_HAT, h)


def eccentricity_vector(r, v):
    """Eccentricity vector: points to periapsis, length = e."""
    r_mag = np.linalg.norm(r)
    return ((np.dot(v, v) - MU_KM3_S2 / r_mag) * r - np.dot(r, v) * v) / MU_KM3_S2


def specific_energy(r, v):
    """Specific orbital energy eps = v^2/2 - mu/r (km^2/s^2)."""
    return np.dot(v, v) / 2 - MU_KM3_S2 / np.linalg.norm(r)


def semi_major_axis(r, v):
    """Semi-major axis a = -mu / (2 eps) (km)."""
    return -MU_KM3_S2 / (2 * specific_energy(r, v))


def inclination_deg(h):
    """Inclination: angle between K and h (deg)."""
    return angle_deg(K_HAT, h)


def invariants(r, v):
    """Validate one state vector and compute all its invariants.

    Returns a dict. Raises ValueError if the state is not a closed orbit
    around the Earth.
    """
    r_mag = np.linalg.norm(r)
    if r_mag <= R_EARTH_KM:
        raise ValueError(f"position {r_mag:.1f} km is inside the Earth")

    energy = specific_energy(r, v)
    if energy >= 0:
        raise ValueError(f"not a closed orbit: energy {energy:.3f} km2/s2 >= 0")

    h_vec = angular_momentum(r, v)
    n_vec = node_vector(h_vec)
    e_vec = eccentricity_vector(r, v)

    return {
        "h_vec": h_vec,
        "n_vec": n_vec,
        "e_vec": e_vec,
        "h": np.linalg.norm(h_vec),
        "e": np.linalg.norm(e_vec),
        "energy": energy,
        "a_km": semi_major_axis(r, v),
        "i_deg": inclination_deg(h_vec),
        "r_dot_v": np.dot(r, v),
        "circular": is_zero(e_vec),
        "equatorial": is_zero(n_vec),
    }


def process_states(states):
    """Compute invariants for every state that was read.

    Returns (results, rejected):
      results  - list of dicts: name, epoch_utc and all invariants
      rejected - list of (line_number, name, reason) for states that are not real orbits
    """
    results = []
    rejected = []
    for line_number, name, epoch, r, v in states:
        try:
            inv = invariants(r, v)
        except ValueError as err:
            rejected.append((line_number, name, str(err)))
            continue
        results.append({"name": name, "epoch_utc": epoch, **inv})
    return results, rejected


# --- 3. Output ---------------------------------------------------------------

TABLE_HEADER = (f"{'Name':<22}{'|h| km2/s':>12}{'e':>11}{'a km':>11}{'i deg':>9}"
                f"{'circular':>10}{'equatorial':>12}")


def yes_no(flag):
    """'yes' / 'no' for the console table."""
    return "yes" if flag else "no"


def format_row(res):
    """One console table row for one result."""
    return (f"{res['name']:<22}{res['h']:>12.1f}{res['e']:>11.6f}{res['a_km']:>11.1f}"
            f"{res['i_deg']:>9.2f}{yes_no(res['circular']):>10}{yes_no(res['equatorial']):>12}")


def format_skipped(entry):
    """One line for a rejected row, used by both the console and the skipped file."""
    line_number, name, reason = entry
    return f"line {line_number}  {name}: {reason}"


def to_csv_row(res):
    """Flatten one result for the CSV: each vector becomes x, y, z columns."""
    row = {"name": res["name"], "epoch_utc": res["epoch_utc"]}
    for vec, prefix in (("h_vec", "h"), ("n_vec", "n"), ("e_vec", "e")):
        for axis, value in zip("xyz", res[vec]):
            row[f"{prefix}_{axis}"] = value
    for key in ("h", "e", "energy", "a_km", "i_deg", "r_dot_v", "circular", "equatorial"):
        row[key] = res[key]
    return row


def write_invariants(results, path):
    """Write one CSV row per result."""
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=RESULT_FIELDS)
        writer.writeheader()
        writer.writerows(to_csv_row(res) for res in results)


def write_skipped(skipped, path):
    """Write one line per rejected row, or say that none were rejected."""
    lines = [format_skipped(entry) for entry in skipped] or ["No rows skipped."]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


# --- Main --------------------------------------------------------------------

def main():
    # 1. Read
    states, skipped = read_states(INPUT_PATH)

    # 2. Process
    results, rejected = process_states(states)
    skipped = sorted(skipped + rejected)   # both kinds of rejection, in file order

    # 3. Output
    print(f"Read {len(results) + len(skipped)} rows from {INPUT_PATH.name}: "
          f"{len(results)} processed, {len(skipped)} skipped\n")
    print(TABLE_HEADER)
    for res in results:
        print(format_row(res))
    if skipped:
        print("\nSkipped rows:")
        for entry in skipped:
            print("  " + format_skipped(entry))

    OUTPUT_DIR.mkdir(exist_ok=True)
    write_invariants(results, INVARIANTS_PATH)
    write_skipped(skipped, SKIPPED_PATH)
    print(f"\nWrote {len(results)} results to {INVARIANTS_PATH}")
    print(f"Wrote {len(skipped)} skipped rows to {SKIPPED_PATH}")


if __name__ == "__main__":
    main()