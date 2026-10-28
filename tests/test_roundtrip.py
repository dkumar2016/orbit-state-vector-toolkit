"""Tests for rv_to_coe: sign tests, round trips, singular cases, invalid input."""

import numpy as np
import pytest

from conversions import angle_diff_deg, coe_to_rv, roundtrip, rv_to_coe

TOL = 1e-6        # plan target for the round-trip error
ANGLE_TOL = 1e-5  # degrees; arccos loses precision near 0 and 180 deg (about 1e-6 deg)


# --- A. The three sign tests ------------------------------------------------------
# (elements in, expected Omega, omega, nu). Each case flips at least one test.

SIGN_CASES = [
    pytest.param((10000.0, 0.3, 30.0, 40.0, 60.0, 45.0), (40.0, 60.0, 45.0), id="no-flip"),
    pytest.param((7000.0, 0.1, 141.0, 250.0, 30.0, 200.0), (250.0, 30.0, 200.0), id="flip-ny-and-rv"),
    pytest.param((26560.0, 0.74, 63.4, 0.0, 270.0, 0.0), (0.0, 270.0, 0.0), id="flip-ez"),
]


@pytest.mark.parametrize("coe, expected", SIGN_CASES)
def test_sign_tests_recover_angles(coe, expected):
    out = rv_to_coe(*coe_to_rv(*coe))
    for key, value in zip(("raan_deg", "argp_deg", "nu_deg"), expected):
        assert abs(angle_diff_deg(out[key], value)) < ANGLE_TOL


# --- B. Round trip, general orbits -------------------------------------------------
# Between them, Omega, omega and nu each appear above and below 180 deg.

GENERAL_ORBITS = [
    pytest.param((10000.0, 0.3, 30.0, 40.0, 60.0, 45.0), id="worked-example"),
    pytest.param((7000.0, 0.1, 141.0, 250.0, 30.0, 200.0), id="retrograde"),
    pytest.param((26560.0, 0.74, 63.4, 0.0, 270.0, 0.0), id="molniya"),
    pytest.param((8000.0, 0.2, 98.0, 300.0, 200.0, 100.0), id="raan-argp-past-180"),
    pytest.param((7200.0, 0.001, 90.0, 120.0, 10.0, 300.0), id="polar-near-circular"),
    pytest.param((42164.0, 0.0001, 0.05, 0.0, 0.0, 0.0), id="near-geo"),
]


@pytest.mark.parametrize("coe", GENERAL_ORBITS)
def test_general_round_trip(coe):
    out, error = roundtrip(*coe)
    assert error < TOL
    assert out["circular"] is False
    assert out["equatorial"] is False


@pytest.mark.parametrize("coe", GENERAL_ORBITS)
def test_argument_of_latitude_is_argp_plus_nu(coe):
    out = rv_to_coe(*coe_to_rv(*coe))
    assert abs(angle_diff_deg(out["u_deg"], coe[4] + coe[5])) < ANGLE_TOL


# --- C. Singular orbits ---------------------------------------------------------------
# (elements in, circular, equatorial, expected Omega, omega, nu)

SINGULAR_CASES = [
    pytest.param((7000.0, 0.0, 30.0, 40.0, 0.0, 75.0), True, False, (40.0, 0.0, 75.0), id="circular-inclined"),
    pytest.param((7000.0, 0.0, 30.0, 40.0, 0.0, 250.0), True, False, (40.0, 0.0, 250.0), id="circular-inclined-u-past-180"),
    pytest.param((9000.0, 0.2, 0.0, 70.0, 50.0, 30.0), False, True, (0.0, 120.0, 30.0), id="equatorial-lonperi-is-sum"),
    pytest.param((9000.0, 0.2, 0.0, 0.0, 250.0, 300.0), False, True, (0.0, 250.0, 300.0), id="equatorial-past-180"),
    pytest.param((7000.0, 0.0, 0.0, 0.0, 0.0, 123.0), True, True, (0.0, 0.0, 123.0), id="circular-equatorial"),
    pytest.param((7000.0, 0.0, 0.0, 0.0, 0.0, 300.0), True, True, (0.0, 0.0, 300.0), id="circular-equatorial-past-180"),
]


@pytest.mark.parametrize("coe, circular, equatorial, expected", SINGULAR_CASES)
def test_singular_flags_and_replacement_angles(coe, circular, equatorial, expected):
    out = rv_to_coe(*coe_to_rv(*coe))
    assert out["circular"] is circular
    assert out["equatorial"] is equatorial
    for key, value in zip(("raan_deg", "argp_deg", "nu_deg"), expected):
        assert abs(angle_diff_deg(out[key], value)) < ANGLE_TOL


@pytest.mark.parametrize("coe, circular, equatorial, expected", SINGULAR_CASES)
def test_singular_round_trip_reproduces_state(coe, circular, equatorial, expected):
    _, error = roundtrip(*coe)
    assert error < TOL


@pytest.mark.parametrize("coe, circular, equatorial, expected", SINGULAR_CASES)
def test_argument_of_latitude_only_when_inclined(coe, circular, equatorial, expected):
    out = rv_to_coe(*coe_to_rv(*coe))
    assert (out["u_deg"] is None) == equatorial


def test_no_nan_for_singular_orbits():
    for case in SINGULAR_CASES:
        out = rv_to_coe(*coe_to_rv(*case.values[0]))
        for key in ("raan_deg", "argp_deg", "nu_deg"):
            assert not np.isnan(out[key])


# --- D. Invalid input -------------------------------------------------------------------

def test_rejects_escape_state():
    with pytest.raises(ValueError, match="not a closed orbit"):
        rv_to_coe(np.array([7000.0, 0.0, 0.0]), np.array([0.0, 11.0, 0.0]))


def test_rejects_state_inside_earth():
    with pytest.raises(ValueError, match="inside the Earth"):
        rv_to_coe(np.array([6000.0, 0.0, 0.0]), np.array([0.0, 8.0, 0.0]))


# --- E. angle_diff_deg ------------------------------------------------------------------

@pytest.mark.parametrize("a, b, expected", [
    (359.99999915, 0.0, -8.5e-7),
    (10.0, 350.0, 20.0),
    (350.0, 10.0, -20.0),
    (45.0, 45.0, 0.0),
])
def test_angle_diff_wraps_round_360(a, b, expected):
    assert angle_diff_deg(a, b) == pytest.approx(expected, abs=1e-9)