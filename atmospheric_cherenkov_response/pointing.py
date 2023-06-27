import numpy as np
import homogeneous_transformation


def init_pointing(azimuth_deg, zenith_deg):
    assert not np.isnan(azimuth_deg)
    assert not np.isnan(zenith_deg)
    return {"azimuth_deg": azimuth_deg, "zenith_deg": zenith_deg}


def make_pointing_key(pointing):
    return "az{:06d}mdeg_zd{:06d}mdeg".format(
        int(pointing["azimuth_deg"] * 1e3), int(pointing["zenith_deg"] * 1e3),
    )


def init_civil_rotation_of_principal_aperture_plane(pointing, mount):
    """
    cable-robot-mount
    -----------------
    no rotation in z-axis

    altitude-azimuth-mount
    ----------------------
    rotates in z-axis
    """
    if mount == "cable-robot-mount":
        return {"repr": None}
    elif mount == "altitude-azimuth-mount":
        return {"repr": None}
    else:
        raise KeyError("No such mount: '{:s}'".format(mount))
