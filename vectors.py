import csv
from pathlib import Path
import numpy as np


# --- Constants ---------------------------------------------------------------

MU_KM3_S2 = 398600.4418   # Earth's gravitational parameter
R_EARTH_KM = 6378.137     # Earth's equatorial radius
K_HAT = np.array([0.0, 0.0, 1.0])   # ECI Z axis (north pole)
ZERO_TOL = 1e-10          # vector length treated as zero (circular / equatorial)

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


if __name__ == "__main__":
    r = np.array([6778.0, 0.0, 0.0])
    v = np.array([0.0, 4.759, 6.013])
    h = angular_momentum(r, v)
    print("h     =", h)
    print("n     =", node_vector(h))
    print("e vec =", eccentricity_vector(r, v))
    print("eps   =", specific_energy(r, v))
    print("a     =", semi_major_axis(r, v))
    print("i     =", inclination_deg(h))