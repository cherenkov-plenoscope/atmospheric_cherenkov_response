import numpy as np
import corsika_primary as cpw


def refractive_index(altitude_asl_m):
    """
    Parameters
    ----------
    altitude_asl_m : float
        Altitude above sea level.

    Returns
    -------
    index : float
        Refractive index of atmosphere.
    """
    return cpw.particles.identification.refractive_index_atmosphere(
        altitude_asl_m=altitude_asl_m
    )
