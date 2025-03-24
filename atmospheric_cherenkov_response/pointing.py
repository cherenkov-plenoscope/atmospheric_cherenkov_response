"""
Defines how the pointing of a telescope's mount effects the
rotation of its principal aperture plane.
"""

import numpy as np
import spherical_coordinates


def Pointing(azimuth_rad, zenith_rad):
    """
    Defines a dict for  the pointing of an instrument.
    Ensures: -PI <= azimuth_rad < +PI

    Parameters
    ----------
    azimuth_rad : float
        Direction of optical axis.
    zenith_rad : float
        Direction of optical axis.

    Returns
    -------
    pointing : dict
    """
    return {
        "azimuth_rad": spherical_coordinates.azimuth_range(azimuth_rad),
        "zenith_rad": zenith_rad,
    }


def pointing_to_civil_rotation(pointing, mount):
    """
    Defines the goemetry of a telescope mount.

    Parameters
    ----------
    poining : dict {azimuth_rad, zenith_rad}
        Direction of the instrument's optical axis.
    mount : str
        Defines the geometry of the mount.
            - ["cable_robot", "hexapod", "steward"]
            - "altitude_azimuth"

    Returns
    -------
    civil rotation : dict
        Can be compiled with homogeneous_transformation.compile to obtain
        rotation matrices.
    """
    if mount in ["cable_robot", "hexapod", "steward"]:
        return make_civil_rotation_for_mounts_without_z_rotation(**pointing)
    elif mount == "altitude_azimuth":
        return make_civil_rotation_for_altitude_azimuth_mount(**pointing)
    else:
        raise KeyError("No such mount: '{:s}'".format(mount))


def make_civil_rotation_for_altitude_azimuth_mount(azimuth_rad, zenith_rad):
    """
    x-axis is magnetic north where azimuth is 0.
    z-axis goes up.
    """
    return {
        "repr": "tait_bryan",
        "xyz_deg": np.array(
            [0.0, np.rad2deg(-zenith_rad), np.rad2deg(-azimuth_rad)]
        ),
    }


def make_civil_rotation_for_mounts_without_z_rotation(
    azimuth_rad,
    zenith_rad,
    eps_rad=1e-7,
):
    """
    Does not rotate in z-axis
    """
    if np.abs(zenith_rad) > eps_rad:
        zenith_direction = np.array([0, 0, 1])
        pointing_direction = np.array(
            spherical_coordinates.az_zd_to_cx_cy_cz(
                azimuth_rad=azimuth_rad,
                zenith_rad=zenith_rad,
            )
        )
        rot_axis = np.cross(zenith_direction.T, pointing_direction)
        rot_axis /= np.linalg.norm(rot_axis)
        rot_angle = spherical_coordinates.angle_between_xyz(
            zenith_direction,
            pointing_direction,
        )
        rot = {
            "repr": "axis_angle",
            "axis": rot_axis,
            "angle_deg": np.rad2deg(rot_angle),
        }
    else:
        rot = {
            "repr": "tait_bryan",
            "xyz_deg": [0, 0, 0],
        }
    return rot
