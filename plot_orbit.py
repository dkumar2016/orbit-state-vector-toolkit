"""3D plot of orbits computed with coe_to_rv: Molniya and the ISS in ECI."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from conversions import coe_to_rv
from vectors import R_EARTH_KM

BASE_DIR = Path(__file__).resolve().parent
FIGURE_PATH = BASE_DIR / "figures" / "orbit_3d.png"

# (name, a_km, e, i_deg, raan_deg, argp_deg, colour)
ORBITS = [
    ("Molniya", 26560.0, 0.74, 63.4, 0.0, 270.0, "tab:red"),
    ("ISS", 6778.0, 0.0005, 51.64, 0.0, 0.0, "tab:blue"),
]


def orbit_track(a_km, e, i_deg, raan_deg, argp_deg, n_points=361):
    """Positions (n_points x 3 array, km) for one full orbit, nu from 0 to 360 deg."""
    nus = np.linspace(0.0, 360.0, n_points)
    return np.array([coe_to_rv(a_km, e, i_deg, raan_deg, argp_deg, nu)[0] for nu in nus])


def draw_earth(ax):
    """Earth as a wireframe sphere, plus the equator."""
    u, w = np.meshgrid(np.linspace(0, 2 * np.pi, 36), np.linspace(0, np.pi, 18))
    x = R_EARTH_KM * np.cos(u) * np.sin(w)
    y = R_EARTH_KM * np.sin(u) * np.sin(w)
    z = R_EARTH_KM * np.cos(w)
    ax.plot_wireframe(x, y, z, color="tab:blue", linewidth=0.3, alpha=0.4)
    t = np.linspace(0, 2 * np.pi, 200)
    ax.plot(R_EARTH_KM * np.cos(t), R_EARTH_KM * np.sin(t), 0 * t,
            color="gray", linewidth=0.8, linestyle="--")


def draw_eci_axes(ax, length):
    """I, J, K axes from Earth's centre."""
    for vec, label in (((1, 0, 0), "I"), ((0, 1, 0), "J"), ((0, 0, 1), "K")):
        ax.quiver(0, 0, 0, *(length * np.array(vec)), color="gray", linewidth=1,
                  arrow_length_ratio=0.08)
        ax.text(*(1.08 * length * np.array(vec)), label, color="gray")


def main():
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(projection="3d")

    draw_earth(ax)
    draw_eci_axes(ax, length=30000)

    for name, a, e, i, raan, argp, colour in ORBITS:
        track = orbit_track(a, e, i, raan, argp)
        ax.plot(track[:, 0], track[:, 1], track[:, 2], color=colour, linewidth=1.5,
                label=f"{name} (a = {a:,.0f} km, e = {e}, i = {i} deg)")
        if e > 0.01:
            peri = coe_to_rv(a, e, i, raan, argp, 0.0)[0]
            apo = coe_to_rv(a, e, i, raan, argp, 180.0)[0]
            ax.scatter(*peri, color=colour, marker="o", s=30)
            ax.scatter(*apo, color=colour, marker="^", s=30)
            ax.text(*peri, "  periapsis", color=colour, fontsize=8)
            ax.text(*apo, "  apoapsis", color=colour, fontsize=8)

    limit = 40000
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.set_zlim(-limit, limit)
    ax.set_box_aspect((1, 1, 1))          # equal scale on all axes, so the shapes are true
    ax.set_xlabel("x (km)")
    ax.set_ylabel("y (km)")
    ax.set_zlabel("z (km)")
    ax.set_title("Orbits from coe_to_rv, in ECI")
    ax.view_init(elev=20, azim=-60)
    ax.legend(loc="upper left", fontsize=8)

    FIGURE_PATH.parent.mkdir(exist_ok=True)
    fig.savefig(FIGURE_PATH, dpi=150, bbox_inches="tight")
    print(f"Saved {FIGURE_PATH}")
    plt.show()


if __name__ == "__main__":
    main()