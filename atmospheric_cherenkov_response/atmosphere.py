import numpy as np


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
    REFRACTION_AIR_273KELVIN_1ATM = 1.00027357
    DENSITY_LENGTH_AIR_M = 8435.0
    # Rise in altitude by which atmospheric density is reduced by 1/e.

    # refractive index of atmosphere vs. altitude
    return 1.0 + (REFRACTION_AIR_273KELVIN_1ATM - 1.0) * np.exp(
        -altitude_asl_m / DENSITY_LENGTH_AIR_M
    )
