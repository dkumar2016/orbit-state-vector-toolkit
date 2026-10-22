import csv
from pathlib import Path
import numpy as np


# --- Constants ---------------------------------------------------------------

MU_KM3_S2 = 398600.4418   # Earth's gravitational parameter
R_EARTH_KM = 6378.137     # Earth's equatorial radius
K_HAT = np.array([0.0, 0.0, 1.0])   # ECI Z axis (north pole)
ZERO_TOL = 1e-10          # vector length treated as zero (circular / equatorial)

R_FIELDS = ("x_km", "y_km", "z_km")  # position
V_FIELDS = ("vx_km_s", "vy_km_s", "vz_km_s")  # velocity


BASE_DIR = Path(__file__).resolve().parent
INPUT_PATH = BASE_DIR / "data" / "states.csv"
OUTPUT_DIR = BASE_DIR / "output"
INVARIANTS_PATH = OUTPUT_DIR / "invariants.csv"
SKIPPED_PATH = OUTPUT_DIR / "states_skipped.txt"


# --- Helpers -----------------------------------------------------------------

def angle_deg(a, b):
    """Angle between two vectors, 0-180 deg (clipped against rounding)."""
    c = np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b))
    return np.degrees(np.arccos(np.clip(c, -1.0, 1.0)))


def is_zero(x, tol=ZERO_TOL):
    """True if a vector's length is below the tolerance."""
    return np.linalg.norm(x) < tol


# --- Vector invariants -------------------------------------------------------

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


#-------------Compute Invariants-----------------

def invariants(r, v):
    """Validate a state vector and compute its vector invariants.

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


#----------------------Read CSV File-------------------------

def read_states(path):
    """Read state vectors from a CSV file.

    Returns (states, skipped):
      states  - list of (r, v) tuples for every valid row
      skipped - list of (line_number, name, reason) for every rejected row
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

if __name__ == "__main__":
    states, skipped = read_states(INPUT_PATH)
    for line_number, name, epoch, r, v in states:
        try:
            inv = invariants(r, v)
            print(f"{name} @ {epoch}: e={inv['e']:.6f}  a={inv['a_km']:.1f} km  i={inv['i_deg']:.2f} deg")
        except ValueError as err:
            skipped.append((line_number, name, str(err)))
    for line_number, name, reason in skipped:
        print(f"Skipped line {line_number}  {name}: {reason}")