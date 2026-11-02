"""Anomalies, Kepler's equation and two-body propagation.

    nu (true)  <->  E (eccentric)  <->  M (mean)  =  n (t - t_p)

nu <-> E is a closed-form half-angle formula. E -> M is Kepler's equation,
M = E - e sin E. M -> E has no closed form and is solved by Newton-Raphson.

The anomaly functions work in radians (they are maths); time_since_periapsis
and propagate take degrees, like coe_to_rv.
"""

import numpy as np

from toolkit.constants import MU_KM3_S2
from toolkit.conversions import coe_to_rv

TWO_PI = 2.0 * np.pi


# --- Mean motion and period ---------------------------------------------------

def mean_motion(a_km):
    """Mean motion n = sqrt(mu / a^3), in rad/s."""
    return np.sqrt(MU_KM3_S2 / a_km**3)


def period_s(a_km):
    """Orbital period T = 2 pi / n, in seconds."""
    return TWO_PI / mean_motion(a_km)


# --- Kepler's equation ----------------------------------------------------------

def solve_kepler(M, e, tol=1e-12, max_iter=50):
    """Solve Kepler's equation M = E - e sin E for E (rad), by Newton-Raphson.

    M in radians (any value; wrapped to 0-2pi), 0 <= e < 1.
    Starts from E = M, or from E = pi when e > 0.8, where E = M can overshoot.
    Raises RuntimeError if it has not converged after max_iter steps.
    """
    if not (0 <= e < 1):
        raise ValueError(f"eccentricity must be 0 <= e < 1, got {e}")
    M = M % TWO_PI
    E = np.pi if e > 0.8 else M
    for _ in range(max_iter):
        step = (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
        E -= step
        if abs(step) < tol:
            return E
    raise RuntimeError(f"Kepler's equation did not converge in {max_iter} iterations "
                       f"(M = {M}, e = {e})")


# --- Converting between anomalies (radians) ----------------------------------------

def true_to_eccentric(nu, e):
    """True anomaly -> eccentric anomaly (rad), using half angles (no quadrant test)."""
    E = 2.0 * np.arctan2(np.sqrt(1 - e) * np.sin(nu / 2), np.sqrt(1 + e) * np.cos(nu / 2))
    return E % TWO_PI


def eccentric_to_true(E, e):
    """Eccentric anomaly -> true anomaly (rad), using half angles (no quadrant test)."""
    nu = 2.0 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2), np.sqrt(1 - e) * np.cos(E / 2))
    return nu % TWO_PI


def eccentric_to_mean(E, e):
    """Kepler's equation: M = E - e sin E (rad)."""
    return (E - e * np.sin(E)) % TWO_PI


def true_to_mean(nu, e):
    """True anomaly -> mean anomaly (rad). Direct: nu -> E -> M."""
    return eccentric_to_mean(true_to_eccentric(nu, e), e)


def mean_to_true(M, e, tol=1e-12):
    """Mean anomaly -> true anomaly (rad). Needs the solver: M -> E -> nu."""
    return eccentric_to_true(solve_kepler(M, e, tol), e)


# --- Time and propagation ----------------------------------------------------------

def time_since_periapsis(a_km, e, nu_deg):
    """Time (s) from periapsis to true anomaly nu (deg), 0 <= t < T."""
    return true_to_mean(np.radians(nu_deg), e) / mean_motion(a_km)


def propagate(a_km, e, i_deg, raan_deg, argp_deg, nu_deg, dt_s):
    """Two-body propagation: elements at t0 -> (r, v) in ECI at t0 + dt_s.

    Only the position along the orbit changes: M advances by n * dt, Kepler's
    equation gives the new nu, and coe_to_rv gives r and v. Angles in degrees.
    """
    M0 = true_to_mean(np.radians(nu_deg), e)
    M = M0 + mean_motion(a_km) * dt_s
    nu_new_deg = np.degrees(mean_to_true(M, e))
    return coe_to_rv(a_km, e, i_deg, raan_deg, argp_deg, nu_new_deg)


if __name__ == "__main__":
    a, e = 26560.0, 0.74                                     # Molniya
    T_min = period_s(a) / 60
    print(f"Molniya: T = {T_min:.2f} min")
    for nu in (90, 180, 270):
        print(f"  nu = {nu:3d} deg reached {time_since_periapsis(a, e, nu) / 60:6.1f} min "
              f"after periapsis")
    M = mean_motion(a) * 3600
    print(f"  1 h after periapsis: nu = {np.degrees(mean_to_true(M, e)):.2f} deg")
