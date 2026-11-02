"""Physical constants and repository paths, defined once for the whole toolkit."""

from pathlib import Path

import numpy as np

# --- Earth --------------------------------------------------------------------

MU_KM3_S2 = 398600.4418             # Earth's gravitational parameter (km^3/s^2)
R_EARTH_KM = 6378.137               # Earth's equatorial radius (km)

# --- ECI unit vectors -----------------------------------------------------------

I_HAT = np.array([1.0, 0.0, 0.0])   # ECI X axis (vernal equinox direction)
K_HAT = np.array([0.0, 0.0, 1.0])   # ECI Z axis (north pole)

# --- Numerical tolerances -------------------------------------------------------

ZERO_TOL = 1e-10                    # vector length treated as zero (circular / equatorial)

# --- Repository paths -------------------------------------------------------------

REPO_DIR = Path(__file__).resolve().parent.parent   # toolkit/ -> repository root
DATA_DIR = REPO_DIR / "data"
OUTPUT_DIR = REPO_DIR / "output"
FIGURES_DIR = REPO_DIR / "figures"
