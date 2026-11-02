"""Tests for kepler.py: the solver, anomaly conversions, timing and propagation."""

import numpy as np
import pytest

from toolkit.conversions import coe_to_rv, rv_to_coe
from toolkit.kepler import (
    eccentric_to_mean,
    eccentric_to_true,
    mean_motion,
    mean_to_true,
    period_s,
    propagate,
    solve_kepler,
    time_since_periapsis,
    true_to_eccentric,
    true_to_mean,
)

ECCENTRICITIES = [0.0, 0.01, 0.3, 0.74, 0.9, 0.95]
MEAN_ANOMALIES = np.radians([0.0, 1.0, 30.0, 90.0, 179.0, 180.0, 181.0, 270.0, 359.0])

MOLNIYA = (26560.0, 0.74, 63.4, 0.0, 270.0)          # a, e, i, raan, argp


# --- A. The solver -------------------------------------------------------------------

def test_circular_orbit_gives_E_equal_M():
    for M in MEAN_ANOMALIES:
        assert solve_kepler(M, 0.0) == pytest.approx(M % (2 * np.pi), abs=1e-12)


@pytest.mark.parametrize("e", ECCENTRICITIES)
def test_solution_satisfies_keplers_equation(e):
    for M in MEAN_ANOMALIES:
        E = solve_kepler(M, e)
        assert E - e * np.sin(E) == pytest.approx(M, abs=1e-12)


@pytest.mark.parametrize("e", ECCENTRICITIES)
def test_converges_in_fewer_than_10_iterations(e):
    for M in MEAN_ANOMALIES:
        solve_kepler(M, e, max_iter=10)          # raises RuntimeError if it needs more


def test_newton_example_from_theory():
    """Molniya, 1 h after periapsis (Block 1, Example 4)."""
    M = mean_motion(26560.0) * 3600.0
    assert solve_kepler(M, 0.74) == pytest.approx(1.2200243, abs=1e-7)


def test_rejects_hyperbolic_eccentricity():
    with pytest.raises(ValueError, match="eccentricity"):
        solve_kepler(1.0, 1.2)


# --- B. Anomaly conversions -----------------------------------------------------------

@pytest.mark.parametrize("e", ECCENTRICITIES)
def test_true_mean_round_trip(e):
    for nu in np.radians([0.0, 45.0, 90.0, 179.0, 181.0, 270.0, 350.0]):
        nu_back = mean_to_true(true_to_mean(nu, e), e)
        assert np.cos(nu_back - nu) == pytest.approx(1.0, abs=1e-12)


def test_anomalies_agree_at_apsides():
    for e in ECCENTRICITIES:
        for nu in (0.0, np.pi):
            E = true_to_eccentric(nu, e)
            assert E == pytest.approx(nu, abs=1e-12)
            assert eccentric_to_mean(E, e) == pytest.approx(nu, abs=1e-12)


def test_worked_example_day8_orbit():
    """a = 10,000 km, e = 0.3, nu = 45 deg (Block 1, Example 2)."""
    E = true_to_eccentric(np.radians(45.0), 0.3)
    assert np.degrees(E) == pytest.approx(33.813, abs=1e-3)
    assert np.degrees(eccentric_to_mean(E, 0.3)) == pytest.approx(24.248, abs=1e-3)
    assert eccentric_to_true(E, 0.3) == pytest.approx(np.radians(45.0))


# --- C. Timing ---------------------------------------------------------------------------

def test_molniya_period():
    assert period_s(26560.0) / 60 == pytest.approx(717.96, abs=0.01)


@pytest.mark.parametrize("nu_deg, minutes", [(90.0, 27.42), (180.0, 358.98), (270.0, 690.54)])
def test_molniya_timing(nu_deg, minutes):
    assert time_since_periapsis(26560.0, 0.74, nu_deg) / 60 == pytest.approx(minutes, abs=0.01)


@pytest.mark.parametrize("e", ECCENTRICITIES)
def test_apoapsis_at_half_period(e):
    assert time_since_periapsis(10000.0, e, 180.0) == pytest.approx(period_s(10000.0) / 2)


# --- D. Propagation -----------------------------------------------------------------------

def test_one_full_period_returns_to_start():
    r0, v0 = coe_to_rv(*MOLNIYA, 37.0)
    r1, v1 = propagate(*MOLNIYA, 37.0, period_s(MOLNIYA[0]))
    assert r1 == pytest.approx(r0, rel=1e-9)
    assert v1 == pytest.approx(v0, rel=1e-9)


def test_half_period_from_periapsis_reaches_apoapsis():
    a, e = MOLNIYA[0], MOLNIYA[1]
    r, _ = propagate(*MOLNIYA, 0.0, period_s(a) / 2)
    assert np.linalg.norm(r) == pytest.approx(a * (1 + e))


def test_propagation_keeps_the_orbit():
    """Two-body motion changes only nu: a, e, i, raan, argp must not change."""
    r, v = propagate(*MOLNIYA, 10.0, 5000.0)
    out = rv_to_coe(r, v)
    for key, value in zip(("a_km", "e", "i_deg", "raan_deg", "argp_deg"), MOLNIYA):
        assert out[key] == pytest.approx(value, abs=1e-6)


def test_one_hour_after_molniya_periapsis():
    """Block 1, Example 4: nu = 122.11 deg after 1 h."""
    r, v = propagate(*MOLNIYA, 0.0, 3600.0)
    assert rv_to_coe(r, v)["nu_deg"] == pytest.approx(122.109, abs=1e-3)
