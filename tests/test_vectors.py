"""Tests for vectors.py: known answers, edge cases, and identities that hold for any orbit."""

import numpy as np
import pytest

from vectors import (
    MU_KM3_S2,
    angular_momentum,
    eccentricity_vector,
    inclination_deg,
    invariants,
    node_vector,
    read_states,
)

# --- Test states (typed here so the tests never depend on data/states.csv) ---

ISS_R = np.array([6778.0, 0.0, 0.0])
ISS_V = np.array([0.0, 4.759, 6.013])

CIRC_R = np.array([7000.0, 0.0, 0.0])
CIRC_V = np.array([0.0, 7.546053290107541, 0.0])

ECC_R = np.array([7000.0, 0.0, 0.0])
ECC_V = np.array([0.0, 7.451130602178055, 4.301912258934557])

GOOD_STATES = [
    pytest.param(ISS_R, ISS_V, id="iss"),
    pytest.param(CIRC_R, CIRC_V, id="equatorial-circular"),
    pytest.param(ECC_R, ECC_V, id="eccentric"),
]


# --- A. Known answers ---------------------------------------------------------

def test_iss_angular_momentum():
    h = angular_momentum(ISS_R, ISS_V)
    assert h == pytest.approx([0.0, -40756.114, 32256.502])


def test_iss_inclination():
    assert inclination_deg(angular_momentum(ISS_R, ISS_V)) == pytest.approx(51.64, abs=0.01)


def test_iss_node_points_to_plus_i():
    n = node_vector(angular_momentum(ISS_R, ISS_V))
    assert n == pytest.approx([40756.114, 0.0, 0.0], abs=1e-6)


def test_eccentric_e_a_and_i():
    inv = invariants(ECC_R, ECC_V)
    assert inv["e"] == pytest.approx(0.3)
    assert inv["a_km"] == pytest.approx(10000.0)
    assert inv["i_deg"] == pytest.approx(30.0)


def test_eccentricity_vector_points_to_periapsis():
    # The eccentric state was built at periapsis, so e must point along r.
    e_vec = eccentricity_vector(ECC_R, ECC_V)
    assert e_vec / np.linalg.norm(e_vec) == pytest.approx(ECC_R / np.linalg.norm(ECC_R))


# --- B. Edge cases ------------------------------------------------------------

def test_circular_orbit_is_flagged():
    inv = invariants(CIRC_R, CIRC_V)
    assert inv["e"] == pytest.approx(0.0, abs=1e-10)
    assert inv["circular"] is True


def test_equatorial_orbit_is_flagged():
    inv = invariants(CIRC_R, CIRC_V)
    assert np.linalg.norm(inv["n_vec"]) == pytest.approx(0.0, abs=1e-10)
    assert inv["equatorial"] is True
    assert inv["i_deg"] == pytest.approx(0.0, abs=1e-10)


def test_iss_is_not_flagged():
    inv = invariants(ISS_R, ISS_V)
    assert inv["circular"] is False
    assert inv["equatorial"] is False


def test_rejects_escape_trajectory():
    with pytest.raises(ValueError, match="not a closed orbit"):
        invariants(np.array([7000.0, 0.0, 0.0]), np.array([0.0, 11.0, 0.0]))


def test_rejects_position_inside_earth():
    with pytest.raises(ValueError, match="inside the Earth"):
        invariants(np.array([6000.0, 0.0, 0.0]), np.array([0.0, 8.0, 0.0]))


# --- C. Identities (must hold for any orbit) -----------------------------------

@pytest.mark.parametrize("r, v", GOOD_STATES)
def test_h_is_perpendicular_to_r_and_v(r, v):
    h = angular_momentum(r, v)
    scale = np.linalg.norm(h)
    assert np.dot(h, r) / (scale * np.linalg.norm(r)) == pytest.approx(0.0, abs=1e-12)
    assert np.dot(h, v) / (scale * np.linalg.norm(v)) == pytest.approx(0.0, abs=1e-12)


@pytest.mark.parametrize("r, v", GOOD_STATES)
def test_e_vector_lies_in_orbit_plane(r, v):
    h = angular_momentum(r, v)
    e_vec = eccentricity_vector(r, v)
    assert np.dot(e_vec, h) / np.linalg.norm(h) == pytest.approx(0.0, abs=1e-12)


@pytest.mark.parametrize("r, v", GOOD_STATES)
def test_h_matches_elements(r, v):
    inv = invariants(r, v)
    expected = np.sqrt(MU_KM3_S2 * inv["a_km"] * (1 - inv["e"] ** 2))
    assert inv["h_km2_s"] == pytest.approx(expected, rel=1e-9)


@pytest.mark.parametrize("r, v", GOOD_STATES)
def test_lagrange_identity(r, v):
    h = angular_momentum(r, v)
    lhs = np.dot(h, h) + np.dot(r, v) ** 2
    rhs = np.dot(r, r) * np.dot(v, v)
    assert lhs == pytest.approx(rhs, rel=1e-12)


# --- D. Reading the file --------------------------------------------------------

def test_read_states_skips_bad_numbers(tmp_path):
    path = tmp_path / "states.csv"
    path.write_text(
        "name,epoch_utc,x_km,y_km,z_km,vx_km_s,vy_km_s,vz_km_s\n"
        "Good,2026-10-21T00:00:00Z,7000,0,0,0,7.5,0\n"
        "Typo,2026-10-21T00:00:00Z,70O0,0,0,0,7.5,0\n",
        encoding="utf-8",
    )
    states, skipped = read_states(path)
    assert len(states) == 1
    line_number, name, epoch, r, v = states[0]
    assert (line_number, name, epoch) == (2, "Good", "2026-10-21T00:00:00Z")
    assert r == pytest.approx([7000.0, 0.0, 0.0])
    assert skipped[0][:2] == (3, "Typo")