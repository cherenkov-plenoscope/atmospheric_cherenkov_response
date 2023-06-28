import atmospheric_cherenkov_response as acr
import homogeneous_transformation
import numpy as np


def assert_close(a, b, eps=1e-6):
    assert np.linalg.norm(a - b) <= eps


def transform_xyz_altitude_azimuth_mount(pointing):
    rot_civil = acr.pointing.make_civil_rotation_for_altitude_azimuth_mount(
        **pointing
    )

    origin = np.array([0, 0, 0])

    t_civil = {"rot": rot_civil, "pos": origin}

    t = homogeneous_transformation.compile(t_civil)

    x = np.array([1, 0, 0])
    y = np.array([0, 1, 0])
    z = np.array([0, 0, 1])

    x_t = homogeneous_transformation.transform_orientation(t, x)
    y_t = homogeneous_transformation.transform_orientation(t, y)
    z_t = homogeneous_transformation.transform_orientation(t, z)

    return x_t, y_t, z_t


def test_altitude_azimuth_mount():
    x, y, z = transform_xyz_altitude_azimuth_mount(
        acr.pointing.init(azimuth_deg=0.0, zenith_deg=0.0)
    )
    assert_close(x, [1, 0, 0])
    assert_close(x, [1, 0, 0])
    assert_close(x, [1, 0, 0])
