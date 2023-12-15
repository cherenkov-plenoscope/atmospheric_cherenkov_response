import numpy as np
import corsika_primary
from . import pointing

"""
Defines the solid angle where an instrument can point to.
For now, this is a simple cone.
"""


def PointingRange_from_cone(
    azimuth_rad=0.0,
    zenith_rad=0.0,
    half_angel_rad=np.deg2rad(60),
):
    return {
        "type": "cone",
        "cone": {
            "azimuth_rad": azimuth_rad,
            "zenith_rad": zenith_rad,
            "half_angel_rad": half_angel_rad,
        },
    }


def draw_pointing(pointing_range, prng):
    """
    Parameters
    ----------
    pointing_range : dict
        Defiens the solid angle where an instrument can point to.
    prng : numpy.random.Generator
        Pseudo random number generator.

    Returns
    -------
    pointing : dict
        azimuth_rad and zenith_rad
    """
    cor = corsika_primary.random.distributions

    cone = pointing_range["cone"]
    az, zd = cor.draw_azimuth_zenith_in_viewcone(
        prng=prng,
        azimuth_rad=cone["azimuth_rad"],
        zenith_rad=cone["zenith_rad"],
        min_scatter_opening_angle_rad=0.0,
        max_scatter_opening_angle_rad=cone["half_angel_rad"],
        max_zenith_rad=np.pi,
        max_iterations=1000000,
    )
    return pointing.Pointing(azimuth_rad=az, zenith_rad=zd)
