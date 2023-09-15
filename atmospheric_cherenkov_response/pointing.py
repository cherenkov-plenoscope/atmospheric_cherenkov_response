import numpy as np
import homogeneous_transformation


def init(azimuth_deg, zenith_deg):
    assert not np.isnan(azimuth_deg)
    assert not np.isnan(zenith_deg)
    return {"azimuth_deg": azimuth_deg, "zenith_deg": zenith_deg}


def make_pointing_key(pointing):
    return "az{:06d}mdeg_zd{:06d}mdeg".format(
        int(pointing["azimuth_deg"] * 1e3),
        int(pointing["zenith_deg"] * 1e3),
    )


def make_civil_rotation_of_principal_aperture_plane(pointing, mount):
    """ """
    if mount == "cable_robot_mount":
        return make_civil_rotation_for_cable_robot_mount(**pointing)
    elif mount == "altitude_azimuth_mount":
        return make_civil_rotation_for_altitude_azimuth_mount(**pointing)
    else:
        raise KeyError("No such mount: '{:s}'".format(mount))


def make_civil_rotation_for_altitude_azimuth_mount(azimuth_deg, zenith_deg):
    """
    x-axis is magnetic north where azimuth is 0deg.
    z-axis goes up.

    """
    return {
        "repr": "tait_bryan",
        "xyz_deg": np.array([0.0, -zenith_deg, -azimuth_deg]),
    }


def make_civil_rotation_for_cable_robot_mount(
    azimuth_deg, zenith_deg, eps_deg=1e-5, eps_1=1e-5
):
    """
    does not rotate in z-axis
    """
    if np.abs(zenith_deg) > eps_deg:
        rot_civil_altaz = make_civil_rotation_for_altitude_azimuth_mount(
            azimuth_deg=azimuth_deg,
            zenith_deg=zenith_deg,
        )
        t_civil_altaz = {"rot": rot_civil_altaz, "pos": [0, 0, 0]}
        t_altaz = homogeneous_transformation.compile(t_civil_altaz)
        z = np.array([0, 0, 1])
        z_t = homogeneous_transformation.transform_orientation(t_altaz, z)
        assert np.abs(np.linalg.norm(z_t) - 1) < eps_1
        pointing_direction = z_t
        rot_axis = np.cross(z.T, pointing_direction)
        rot_axis /= np.linalg.norm(rot_axis)
        assert np.abs(np.linalg.norm(rot_axis) - 1) < eps_1
        angle_z_and_pointing_rad = np.arccos(np.dot(pointing_direction, z))
        angle_z_and_pointing_deg = np.rad2deg(angle_z_and_pointing_rad)
        assert np.abs(zenith_deg) - angle_z_and_pointing_deg < eps_1
        rot = {
            "repr": "axis_angle",
            "axis": rot_axis,
            "angle_deg": angle_z_and_pointing_deg,
        }
    else:
        rot = {
            "repr": "tait_bryan",
            "xyz_deg": [0, 0, 0],
        }

    return rot
