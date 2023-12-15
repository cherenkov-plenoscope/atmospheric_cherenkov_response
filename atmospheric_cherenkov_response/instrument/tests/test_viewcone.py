import atmospheric_cherenkov_response as acr
import numpy as np

TAU = 2 * np.pi


def assert_cone_edge(az_rad, ha_rad, d, margin=1e-9):
    _d = acr.instrument.viewcone.make_direction_on_edge_of_viewcone_when_viewcone_pointing_to_z(
        viewcone_azimuth_rad=az_rad, viewcone_half_angle_rad=ha_rad
    )
    assert np.linalg.norm(_d - d) < 1e-9


def test_edge_of_viewcone():
    assert_cone_edge(az_rad=0.0, ha_rad=0.0, d=[0, 0, 1])
    assert_cone_edge(az_rad=0.0, ha_rad=1 / 4 * TAU, d=[1, 0, 0])
    assert_cone_edge(az_rad=1 / 4 * TAU, ha_rad=1 / 4 * TAU, d=[0, 1, 0])
    assert_cone_edge(az_rad=1 / 2 * TAU, ha_rad=1 / 4 * TAU, d=[-1, 0, 0])
    assert_cone_edge(az_rad=3 / 4 * TAU, ha_rad=1 / 4 * TAU, d=[0, -1, 0])
    assert_cone_edge(
        az_rad=0.0, ha_rad=1 / 8 * TAU, d=[np.sqrt(1 / 2), 0, np.sqrt(1 / 2)]
    )


def assert_pointing(az_rad, zd_rad, d, margin=1e-9):
    _d = acr.instrument.viewcone.make_direction_on_edge_of_viewcone_when_viewcone_pointing(
        pointing_azimuth_rad=az_rad,
        pointing_zenith_rad=zd_rad,
        viewcone_azimuth_rad=0.0,
        viewcone_half_angle_rad=0.0,
    )
    assert np.linalg.norm(_d - d) < margin


def test_azimuth_zenith():
    assert_pointing(az_rad=0.0, zd_rad=0.0, d=[0, 0, 1])
    assert_pointing(az_rad=0.0, zd_rad=1 / 4 * TAU, d=[1, 0, 0])
    assert_pointing(az_rad=1 / 4 * TAU, zd_rad=1 / 4 * TAU, d=[0, 1, 0])
    assert_pointing(az_rad=1 / 2 * TAU, zd_rad=1 / 4 * TAU, d=[-1, 0, 0])
