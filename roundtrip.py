"""Round-trip every orbit in data/orbits.csv: elements -> (r, v) -> elements.

Prints the recovered elements and the worst error for each orbit, then the
maximum over all orbits, against the target of 1e-6.
"""

from conversions import INPUT_PATH, format_skipped, roundtrip
from describe import read_orbits

TARGET = 1e-6

TABLE_HEADER = (f"{'Name':<20}{'a km':>10}{'e':>10}{'i deg':>9}{'Ω deg':>9}"
                f"{'ω deg':>9}{'ν deg':>9}{'max error':>12}")


def format_row(name, out, error):
    """One table row: recovered elements and the worst error for one orbit."""
    return (f"{name:<20}{out['a_km']:>10.1f}{out['e']:>10.6f}{out['i_deg']:>9.3f}"
            f"{out['raan_deg']:>9.3f}{out['argp_deg']:>9.3f}{out['nu_deg']:>9.3f}"
            f"{error:>12.1e}")


def main():
    # 1. Read
    orbits, skipped = read_orbits(INPUT_PATH)

    # 2. Process
    results = []
    for el in orbits:
        out, error = roundtrip(el.a_km, el.e, el.i_deg, el.raan_deg, el.argp_deg, el.nu_deg)
        results.append((el.name, out, error))

    # 3. Output
    print(f"Round trip of {len(results)} orbits from {INPUT_PATH.name} "
          f"({len(skipped)} rows skipped)\n")
    print(TABLE_HEADER)
    for name, out, error in results:
        print(format_row(name, out, error))
    if skipped:
        print("\nSkipped rows:")
        for entry in skipped:
            print("  " + format_skipped(entry))

    worst = max(error for _, _, error in results)
    verdict = "PASS" if worst < TARGET else "FAIL"
    print(f"\nMaximum error over {len(results)} orbits: {worst:.1e}  "
          f"(target {TARGET:.0e})  {verdict}")


if __name__ == "__main__":
    main()