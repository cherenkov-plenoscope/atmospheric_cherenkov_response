import numpy as np
import homogeneous_transformation
from .. import pointing

EXAMPLE = {
    "half_angle_deg": 3.25,
    "azimuth_deg": 0.0,
    "zenith_deg": 45.0,
}


def pointing_to_civil_rotation(azimuth_deg, zenith_deg):
    """
    X-axis is north where azimuth is 0deg.
    """
    return {
        "repr": "tait_bryan",
        "xyz_deg": np.array([0.0, -zenith_deg, -azimuth_deg]),
    }


def make_direction_on_edge_of_viewcone_when_viewcone_pointing_to_z(
    viewcone_azimuth_deg, viewcone_half_angle_deg
):
    sin = np.sin
    cos = np.cos
    theta = np.deg2rad(viewcone_half_angle_deg)
    phi = np.deg2rad(viewcone_azimuth_deg)

    return np.array([sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)])


def make_direction_on_edge_of_viewcone_when_viewcone_pointing(
    pointing_azimuth_deg,
    pointing_zenith_deg,
    viewcone_azimuth_deg,
    viewcone_half_angle_deg,
):
    dir_z = make_direction_on_edge_of_viewcone_when_viewcone_pointing_to_z(
        viewcone_azimuth_deg=viewcone_azimuth_deg,
        viewcone_half_angle_deg=viewcone_half_angle_deg,
    )

    rot_civil = pointing.make_civil_rotation_for_altitude_azimuth_mount(
        azimuth_deg=pointing_azimuth_deg, zenith_deg=pointing_zenith_deg,
    )
    t_civil = {
        "pos": [0, 0, 0],
        "rot": rot_civil,
    }
    t = homogeneous_transformation.compile(t_civil)
    return homogeneous_transformation.transform_orientation(t=t, d=dir_z)
