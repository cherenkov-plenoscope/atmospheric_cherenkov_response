import numpy as np
import homogeneous_transformation


def init_pointing(azimuth_deg, zenith_deg):
    assert not np.isnan(azimuth_deg)
    assert not np.isnan(zenith_deg)
    return {"azimuth_deg": azimuth_deg, "zenith_deg": zenith_deg}
