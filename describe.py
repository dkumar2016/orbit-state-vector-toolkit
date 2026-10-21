"""Orbit describer: reads classical orbital elements from a CSV file and
describes each orbit in plain English (shape, regime, inclination class,
altitudes, period).

Pipeline:  read (CSV -> Elements)  ->  process (Elements -> summary)  ->  output
"""

import csv
import math
from dataclasses import dataclass
from pathlib import Path

# --- Constants ---------------------------------------------------------------

MU_KM3_S2 = 398600.4418   # Earth's gravitational parameter
R_EARTH_KM = 6378.137     # Earth's equatorial radius

GEO_ALT_KM = 35786        # geostationary altitude
GEO_BAND_KM = 200         # +/- band counted as GEO
LEO_MAX_ALT_KM = 2000     # upper limit of LEO
HEO_MIN_E = 0.5           # eccentricity at which an orbit counts as highly elliptical
SSO_BAND_DEG = (96, 102)  # inclinations flagged as sun-synchronous for LEO

BASE_DIR = Path(__file__).resolve().parent
INPUT_PATH = BASE_DIR / "data" / "orbits.csv"
OUTPUT_DIR = BASE_DIR / "output"
DESCRIPTIONS_PATH = OUTPUT_DIR / "descriptions.csv"
SKIPPED_PATH = OUTPUT_DIR / "skipped.txt"

FLOAT_FIELDS = ("a_km", "e", "i_deg", "raan_deg", "argp_deg", "nu_deg")
DESCRIPTION_FIELDS = ("name", "shape", "regime", "inclination_class", "i_deg",
                      "perigee_alt_km", "apogee_alt_km", "period_min", "max_lat_deg")


# --- Data model --------------------------------------------------------------

def wrap_360(angle_deg):
    """Wrap any angle into the range [0, 360)."""
    return angle_deg % 360


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
            raise ValueError(f"eccentricity must be 0 <= e < 1, got {self.e}")
        r_p = self.a_km * (1 - self.e)
        if r_p <= R_EARTH_KM:
            raise ValueError(f"periapsis {r_p:.1f} km is inside the Earth")
        if not (0 <= self.i_deg <= 180):
            raise ValueError(f"inclination must be 0-180 deg, got {self.i_deg}")
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

    @property
    def max_lat_deg(self):
        """Highest latitude the ground track reaches."""
        return self.i_deg if self.i_deg <= 90 else 180 - self.i_deg


# --- 1. Read -----------------------------------------------------------------

def read_orbits(path):
    """Read orbits from a CSV file.

    Returns (orbits, skipped):
      orbits  - list of Elements for every valid row
      skipped - list of (line_number, name, reason) for every rejected row
    """
    orbits = []
    skipped = []
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for line_number, row in enumerate(reader, start=2):   # line 1 is the header
            name = row["name"]
            try:
                values = {field: float(row[field]) for field in FLOAT_FIELDS}
                orbits.append(Elements(name=name, **values))
            except (ValueError, TypeError) as err:
                skipped.append((line_number, name, str(err)))
    return orbits, skipped


# --- 2. Process --------------------------------------------------------------

def shape_label(e):
    """Describe the orbit's shape from its eccentricity."""
    if e < 0.001:
        return "circular"
    if e < 0.1:
        return "near-circular"
    if e < HEO_MIN_E:
        return "elliptical"
    return "highly elliptical"


def regime_label(a_km, e):
    """Orbit regime. Highly elliptical orbits are HEO, whatever their mean altitude."""
    if e >= HEO_MIN_E:
        return "HEO"
    mean_alt_km = a_km - R_EARTH_KM
    if mean_alt_km <= LEO_MAX_ALT_KM:
        return "LEO"
    if abs(mean_alt_km - GEO_ALT_KM) <= GEO_BAND_KM:
        return "GEO"
    if mean_alt_km < GEO_ALT_KM:
        return "MEO"
    return "Beyond GEO"


def inclination_label(i_deg, regime):
    """Inclination class, flagging the sun-synchronous band for LEO orbits."""
    if i_deg < 5:
        label = "equatorial"
    elif i_deg < 85:
        label = "prograde"
    elif i_deg <= 95:
        label = "polar"
    elif i_deg <= 175:
        label = "retrograde"
    else:
        label = "retrograde equatorial"
    low, high = SSO_BAND_DEG
    if regime == "LEO" and low <= i_deg <= high:
        label += " (sun-synchronous band)"
    return label


def summarize(el):
    """Every derived label and value for one orbit, as a dict (one CSV row)."""
    regime = regime_label(el.a_km, el.e)
    return {
        "name": el.name,
        "shape": shape_label(el.e),
        "regime": regime,
        "inclination_class": inclination_label(el.i_deg, regime),
        "i_deg": el.i_deg,
        "perigee_alt_km": round(el.perigee_alt_km, 1),
        "apogee_alt_km": round(el.apogee_alt_km, 1),
        "period_min": round(el.period_min, 2),
        "max_lat_deg": round(el.max_lat_deg, 2),
    }


# --- 3. Output ---------------------------------------------------------------

def format_period(minutes):
    """'92.6 min' for short periods, '11.97 h (718.0 min)' for long ones."""
    if minutes < 120:
        return f"{minutes:.1f} min"
    return f"{minutes / 60:.2f} h ({minutes:.1f} min)"


def describe(s):
    """One plain-English sentence from a summary dict."""
    return (
        f"{s['name']}: {s['shape']} {s['regime']}, "
        f"altitude {s['perigee_alt_km']:,.0f}-{s['apogee_alt_km']:,.0f} km, "
        f"{s['inclination_class']} at {s['i_deg']:g} deg, "
        f"period {format_period(s['period_min'])}. "
        f"Ground track reaches +/-{s['max_lat_deg']:g} deg latitude."
    )


def format_skipped(entry):
    """One line for a rejected row, used by both the console and skipped.txt."""
    line_number, name, reason = entry
    return f"line {line_number}  {name}: {reason}"


def write_descriptions(summaries, path):
    """Write one CSV row per summary."""
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=DESCRIPTION_FIELDS)
        writer.writeheader()
        writer.writerows(summaries)


def write_skipped(skipped, path):
    """Write one line per rejected row, or say that none were rejected."""
    lines = [format_skipped(entry) for entry in skipped] or ["No rows skipped."]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    # 1. Read
    orbits, skipped = read_orbits(INPUT_PATH)

    # 2. Process: summarize each orbit once; everything below uses these
    summaries = [summarize(el) for el in orbits]

    # 3. Output
    print(f"Read {len(orbits) + len(skipped)} rows from {INPUT_PATH.name}: "
          f"{len(orbits)} described, {len(skipped)} skipped\n")
    for s in summaries:
        print(describe(s))
    if skipped:
        print("\nSkipped rows:")
        for entry in skipped:
            print("  " + format_skipped(entry))

    OUTPUT_DIR.mkdir(exist_ok=True)
    write_descriptions(summaries, DESCRIPTIONS_PATH)
    write_skipped(skipped, SKIPPED_PATH)
    print(f"\nWrote {len(summaries)} descriptions to {DESCRIPTIONS_PATH}")
    print(f"Wrote {len(skipped)} skipped rows to {SKIPPED_PATH}")


if __name__ == "__main__":
    main()