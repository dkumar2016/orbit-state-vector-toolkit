"""Tests for conversions.py: known answers, bad input, properties, and cross-checks with vectors.py."""

import numpy as np
import pytest

from conversions import (
    MU_KM3_S2,
    STATE_FIELDS,
    coe_to_rv,
    rotation_matrix,
    to_csv_row,
    write_states,
)
from vectors import (
    angular_momentum,
    eccentricity_vector,
    inclination_deg,
    node_vector,
    read_states,
    semi_major_axis,
)

# Element sets (a_km, e, i_deg, raan_deg, argp_deg, nu_deg) used by the property tests.
# Chosen to cover prograde and retrograde, low and high e, and nu on both halves of the orbit.
ORBITS = [
    pytest.param((10000.0, 0.3, 30.0, 40.0, 60.0, 45.0), id="worked-example"),
    pytest.param((26560.0, 0.74, 63.4, 0.0, 270.0, 0.0), id="molniya"),
    pytest.param((7000.0, 0.1, 141.0, 250.0, 30.0, 200.0), id="retrograde"),
    pytest.param((7200.0, 0.001, 90.0, 120.0, 10.0, 300.0), id="polar"),
]


# --- A. Known answers (the two hand calculations) ------------------------------

def test_circular_hand_check():
    r, v = coe_to_rv(6778.0, 0.0, 51.64, 0.0, 0.0, 0.0)
    assert r == pytest.approx([6778.0, 0.0, 0.0], abs=1e-9)
    assert v == pytest.approx([0.0, 4.75916, 6.01318], abs=1e-5)


def test_elliptical_worked_example():
    r, v = coe_to_rv(10000.0, 0.3, 30.0, 40.0, 60.0, 45.0)
    assert r == pytest.approx([-5525.2395, 3561.8490, 3625.8117], abs=1e-3)
    assert v == pytest.approx([-5.81348, -5.69228, -0.36008], abs=1e-4)


def test_rotation_matrix_worked_example():
    Q = rotation_matrix(*np.radians([40.0, 30.0, 60.0]))
    expected = [[-0.09907, -0.94175, 0.32139],
                [0.89593, -0.22496, -0.38302],
                [0.43301, 0.25000, 0.86603]]
    assert Q == pytest.approx(np.array(expected), abs=1e-5)


# --- B. Bad input ----------------------------------------------------------------

def test_rejects_hyperbolic_eccentricity():
    with pytest.raises(ValueError, match="eccentricity"):
        coe_to_rv(10000.0, 1.2, 30.0, 40.0, 60.0, 45.0)


def test_rejects_negative_semi_major_axis():
    with pytest.raises(ValueError, match="semi-major axis"):
        coe_to_rv(-7000.0, 0.1, 30.0, 40.0, 60.0, 45.0)


# --- C. Properties (must hold for any orbit) --------------------------------------

@pytest.mark.parametrize("coe", ORBITS)
def test_rotation_matrix_is_orthonormal(coe):
    _, _, i, raan, argp, _ = coe
    Q = rotation_matrix(*np.radians([raan, i, argp]))
    assert Q.T @ Q == pytest.approx(np.eye(3), abs=1e-12)
    assert np.linalg.det(Q) == pytest.approx(1.0)


@pytest.mark.parametrize("coe", ORBITS)
def test_radius_matches_orbit_equation(coe):
    a, e, _, _, _, nu = coe
    r, _ = coe_to_rv(*coe)
    p = a * (1 - e**2)
    assert np.linalg.norm(r) == pytest.approx(p / (1 + e * np.cos(np.radians(nu))))


@pytest.mark.parametrize("coe", ORBITS)
def test_speed_matches_vis_viva(coe):
    a = coe[0]
    r, v = coe_to_rv(*coe)
    expected = np.sqrt(MU_KM3_S2 * (2 / np.linalg.norm(r) - 1 / a))
    assert np.linalg.norm(v) == pytest.approx(expected)


# --- D. Cross-checks with vectors.py (Day 7) ---------------------------------------

@pytest.mark.parametrize("coe", ORBITS)
def test_inclination_recovered(coe):
    r, v = coe_to_rv(*coe)
    assert inclination_deg(angular_momentum(r, v)) == pytest.approx(coe[2])


@pytest.mark.parametrize("coe", ORBITS)
def test_semi_major_axis_recovered(coe):
    r, v = coe_to_rv(*coe)
    assert semi_major_axis(r, v) == pytest.approx(coe[0])


@pytest.mark.parametrize("coe", ORBITS)
def test_eccentricity_recovered(coe):
    r, v = coe_to_rv(*coe)
    assert np.linalg.norm(eccentricity_vector(r, v)) == pytest.approx(coe[1], abs=1e-9)


@pytest.mark.parametrize("coe", ORBITS)
def test_raan_recovered_from_node_vector(coe):
    r, v = coe_to_rv(*coe)
    n = node_vector(angular_momentum(r, v))
    raan = np.degrees(np.arctan2(n[1], n[0])) % 360
    assert raan == pytest.approx(coe[3])


# --- E. Output format -------------------------------------------------------------

def test_csv_row_has_states_csv_columns():
    r, v = coe_to_rv(10000.0, 0.3, 30.0, 40.0, 60.0, 45.0)
    assert tuple(to_csv_row("Test", r, v)) == STATE_FIELDS


def test_written_file_reads_back_with_vectors_reader(tmp_path):
    r, v = coe_to_rv(10000.0, 0.3, 30.0, 40.0, 60.0, 45.0)
    path = tmp_path / "states.csv"
    write_states([("Test", r, v)], path)
    states, skipped = read_states(path)
    assert skipped == []
    _, name, _, r_read, v_read = states[0]
    assert name == "Test"
    assert r_read == pytest.approx(r)
    assert v_read == pytest.approx(v)