"""Orbit describer: reads classical orbital elements and describes each orbit
in plain English (shape, regime, inclination class, altitudes, period)."""

import math
from dataclasses import dataclass
from pathlib import Path

MU_KM3_S2 = 398600.4418   # Earth's gravitational parameter
R_EARTH_KM = 6378.137     # Earth's equatorial radius

BASE_DIR = Path(__file__).resolve().parent
INPUT_PATH = BASE_DIR / "data" / "orbits.csv"
OUTPUT_DIR = BASE_DIR / "output"


def wrap_360(angle_deg):
    """Wrap any angle into the range [0, 360)."""
    return angle_deg % 360


def max_latitude_deg(i_deg):
    """Highest latitude the ground track reaches: i for prograde, 180 - i for retrograde."""
    return i_deg if i_deg <= 90 else 180 - i_deg


@dataclass
class Elements:
    """Classical orbital elements for one orbit. Angles in degrees, distances in km."""
    name: str
    a_km: float
    e: float
    i_deg: float
    raan_deg: float
    argp_deg: float
    nu_deg: float

    def __post_init__(self):
        if not (0 <= self.e < 1):
            raise ValueError(f"{self.name}: eccentricity must be 0 <= e < 1, got {self.e}")
        r_p = self.a_km * (1 - self.e)
        if r_p <= R_EARTH_KM:
            raise ValueError(f"{self.name}: periapsis {r_p:.1f} km is inside the Earth")
        if not (0 <= self.i_deg <= 180):
            raise ValueError(f"{self.name}: inclination must be 0-180 deg, got {self.i_deg}")
        self.raan_deg = wrap_360(self.raan_deg)
        self.argp_deg = wrap_360(self.argp_deg)
        self.nu_deg = wrap_360(self.nu_deg)

    @property
    def perigee_alt_km(self):
        return self.a_km * (1 - self.e) - R_EARTH_KM

    @property
    def apogee_alt_km(self):
        return self.a_km * (1 + self.e) - R_EARTH_KM

    @property
    def period_min(self):
        return 2 * math.pi * math.sqrt(self.a_km**3 / MU_KM3_S2) / 60


if __name__ == "__main__":
    iss = Elements("ISS", 6778, 0.0005, 51.64, 0, 0, 0)
    print(iss)
    print(f"Perigee {iss.perigee_alt_km:.1f} km, apogee {iss.apogee_alt_km:.1f} km, "
          f"period {iss.period_min:.2f} min, max latitude {max_latitude_deg(iss.i_deg)} deg")
    print(f"CSV found: {INPUT_PATH.is_file()}")