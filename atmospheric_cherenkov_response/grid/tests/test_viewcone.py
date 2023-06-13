import atmospheric_cherenkov_response as acr
import numpy as np


def assert_cone_edge(az_deg, ha_deg, d, margin=1e-9):
    _d = acr.grid.viewcone.make_direction_on_edge_of_viewcone_when_viewcone_pointing_to_z(
        viewcone_azimuth_deg=az_deg, viewcone_half_angle_deg=ha_deg
    )
    assert np.linalg.norm(_d - d) < 1e-9


def test_edge_of_viewcone():
    assert_cone_edge(az_deg=0.0, ha_deg=0.0, d=[0, 0, 1])
    assert_cone_edge(az_deg=0.0, ha_deg=90.0, d=[1, 0, 0])
    assert_cone_edge(az_deg=90.0, ha_deg=90.0, d=[0, 1, 0])
    assert_cone_edge(az_deg=180.0, ha_deg=90.0, d=[-1, 0, 0])
    assert_cone_edge(az_deg=270.0, ha_deg=90.0, d=[0, -1, 0])
    assert_cone_edge(
        az_deg=0.0, ha_deg=45.0, d=[np.sqrt(1 / 2), 0, np.sqrt(1 / 2)]
    )


def assert_pointing(az_deg, zd_deg, d, margin=1e-9):
    _d = acr.grid.viewcone.make_direction_on_edge_of_viewcone_when_viewcone_pointing(
        pointing_azimuth_deg=az_deg,
        pointing_zenith_deg=zd_deg,
        viewcone_azimuth_deg=0.0,
        viewcone_half_angle_deg=0.0,
    )
    assert np.linalg.norm(_d - d) < margin


def test_azimuth_zenith():
    assert_pointing(az_deg=0.0, zd_deg=0.0, d=[0, 0, 1])
    assert_pointing(az_deg=0.0, zd_deg=90.0, d=[1, 0, 0])
    assert_pointing(az_deg=90.0, zd_deg=90.0, d=[0, 1, 0])
    assert_pointing(az_deg=180.0, zd_deg=90.0, d=[-1, 0, 0])
