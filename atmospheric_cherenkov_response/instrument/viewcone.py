import numpy as np
import homogeneous_transformation
from .. import pointing

EXAMPLE = {
    "half_angle_rad": 3.25,
    "azimuth_rad": 0.0,
    "zenith_rad": 45.0,
}


def make_direction_on_edge_of_viewcone_when_viewcone_pointing_to_z(
    viewcone_azimuth_rad, viewcone_half_angle_rad
):
    sin = np.sin
    cos = np.cos
    theta = viewcone_half_angle_rad
    phi = viewcone_azimuth_rad

    return np.array([sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)])


def make_direction_on_edge_of_viewcone_when_viewcone_pointing(
    pointing_azimuth_rad,
    pointing_zenith_rad,
    viewcone_azimuth_rad,
    viewcone_half_angle_rad,
):
    dir_z = make_direction_on_edge_of_viewcone_when_viewcone_pointing_to_z(
        viewcone_azimuth_rad=viewcone_azimuth_rad,
        viewcone_half_angle_rad=viewcone_half_angle_rad,
    )

    rot_civil = pointing.pointing_to_civil_rotation(
        {
            "azimuth_rad": pointing_azimuth_rad,
            "zenith_rad": pointing_zenith_rad,
        },
        mount="cable_robot",
    )
    t_civil = {
        "pos": [0, 0, 0],
        "rot": rot_civil,
    }
    t = homogeneous_transformation.compile(t_civil)
    return homogeneous_transformation.transform_orientation(t=t, d=dir_z)
