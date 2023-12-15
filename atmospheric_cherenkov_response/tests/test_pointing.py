import atmospheric_cherenkov_response as acr
import homogeneous_transformation
import numpy as np

ISQ = 1 / np.sqrt(2)
TAU = 2.0 * np.pi

CASES_AZIMUTH_ZERO = [
    {"az": 0, "zd": 0, "x": [1, 0, 0], "y": [0, 1, 0], "z": [0, 0, 1]},
    {
        "az": 0,
        "zd": (1 / 4) * TAU,
        "x": [0, 0, -1],
        "y": [0, 1, 0],
        "z": [1, 0, 0],
    },
    {
        "az": 0,
        "zd": (1 / 8) * TAU,
        "x": [ISQ, 0, -ISQ],
        "y": [0, 1, 0],
        "z": [ISQ, 0, ISQ],
    },
    {
        "az": 0,
        "zd": -(1 / 8) * TAU,
        "x": [ISQ, 0, ISQ],
        "y": [0, 1, 0],
        "z": [-ISQ, 0, ISQ],
    },
]


def assert_close(a, b, eps=1e-6):
    assert np.linalg.norm(a - b) <= eps


def transform_xyz(rot_civil):
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


def assert_case(mount, case):
    c = case
    rot_civil = acr.pointing.pointing_to_civil_rotation(
        pointing=acr.pointing.Pointing(
            azimuth_rad=c["az"], zenith_rad=c["zd"]
        ),
        mount=mount,
    )
    x, y, z = transform_xyz(rot_civil)
    assert_close(x, c["x"])
    assert_close(y, c["y"])
    assert_close(z, c["z"])


def test_altitude_azimuth_mount():
    for case in CASES_AZIMUTH_ZERO:
        assert_case(mount="altitude_azimuth", case=case)

    cases = [
        {
            "az": (1 / 4) * TAU,
            "zd": 0,
            "x": [0, 1, 0],
            "y": [-1, 0, 0],
            "z": [0, 0, 1],
        },
        {
            "az": (1 / 2) * TAU,
            "zd": 0,
            "x": [-1, 0, 0],
            "y": [0, -1, 0],
            "z": [0, 0, 1],
        },
        {
            "az": (3 / 4) * TAU,
            "zd": 0,
            "x": [0, -1, 0],
            "y": [1, 0, 0],
            "z": [0, 0, 1],
        },
    ]

    for case in cases:
        assert_case(mount="altitude_azimuth", case=case)


def test_cable_robot_mount():
    for case in CASES_AZIMUTH_ZERO:
        assert_case(mount="cable_robot", case=case)

    # rotations in azimuth do nothing when zenith is zero
    for az in np.linspace(-137, 137, 111):
        c = {"az": az, "zd": 0, "x": [1, 0, 0], "y": [0, 1, 0], "z": [0, 0, 1]}
        assert_case(mount="cable_robot", case=c)

    # rotate around x-axis into +y
    case = {
        "az": (1 / 4) * TAU,
        "zd": (1 / 8) * TAU,
        "x": [1, 0, 0],
        "y": [0, ISQ, -ISQ],
        "z": [0, ISQ, ISQ],
    }
    assert_case(mount="cable_robot", case=case)

    # rotate around x-axis into -y
    case = {
        "az": (1 / 4) * TAU,
        "zd": -(1 / 8) * TAU,
        "x": [1, 0, 0],
        "y": [0, ISQ, ISQ],
        "z": [0, -ISQ, ISQ],
    }
    assert_case(mount="cable_robot", case=case)
