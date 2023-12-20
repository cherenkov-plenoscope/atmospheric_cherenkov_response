import copy
import numpy as np
from . import atmosphere


def _all():
    sites = {}
    sites["namibia"] = {
        "observation_level_asl_m": 2300,
        "earth_magnetic_field_x_muT": 12.5,
        "earth_magnetic_field_z_muT": -25.9,
        "corsika_atmosphere_id": 10,
        "geomagnetic_cutoff_rigidity_GV": 12.5,
        "coordinates_wgs1984": [-23.3425, 16.225556],
        "comment": "The Gamsberg-mesa in Khoma, Namibia, southern Africa.",
        "plotting": {"label": "Gamsberg", "marker": "+", "linestyle": "--"},
    }

    sites["chile"] = {
        "observation_level_asl_m": 5000,
        "earth_magnetic_field_x_muT": 20.815,
        "earth_magnetic_field_z_muT": -11.366,
        "corsika_atmosphere_id": 26,
        "geomagnetic_cutoff_rigidity_GV": 10.0,
        "coordinates_wgs1984": [-23.0193, -67.7532],
        "comment": "Llano de Chajnantor in Chile, southern America.",
        "plotting": {"label": "Chajnantor", "marker": "*", "linestyle": ":"},
    }

    sites["lapalma"] = {
        "observation_level_asl_m": 2200,
        "earth_magnetic_field_x_muT": 30.419,
        "earth_magnetic_field_z_muT": 23.856,
        "corsika_atmosphere_id": 8,
        "geomagnetic_cutoff_rigidity_GV": 10.0,
        "coordinates_wgs1984": [28.7615, -17.8906],
        "comment": "Roque de los Muchachos Observatory on La Palma, Spain",
        "plotting": {
            "label": "Roque",
            "marker": "^",
            "linestyle": "-",
        },
    }

    PRACTICALLY_ZERO = (
        1e-9  # Corsika's magnetic field must not be exactly zero.
    )

    sites["namibiaOff"] = copy.deepcopy(sites["namibia"])
    sites["namibiaOff"][
        "comment"
    ] += " Here earth's magnetic field is zero. Only for demonstration!"
    sites["namibiaOff"]["earth_magnetic_field_x_muT"] = PRACTICALLY_ZERO
    sites["namibiaOff"]["earth_magnetic_field_z_muT"] = PRACTICALLY_ZERO
    sites["namibiaOff"]["geomagnetic_cutoff_rigidity_GV"] = 0.0
    sites["namibiaOff"]["plotting"] = {
        "label": "Gamsberg-Off",
        "marker": ".",
        "linestyle": "-.",
    }

    # assert all have same keys
    # -------------------------
    for sk in sites:
        for key in sites[sk]:
            assert (
                key in sites["namibia"]
            ), "The key '{:s}' from '{:s}' not in 'namibia'".format(key, sk)
        for key in sites["namibia"]:
            assert key in sites[sk], "{:s} not in '{:s}'".format(key, sk)

    # add refractive index at observation level
    for sk in sites:
        sites[sk][
            "atmosphere_refractive_index_at_observation_level"
        ] = atmosphere.refractive_index(
            altitude_asl_m=sites[sk]["observation_level_asl_m"]
        )

    # put key into dict
    # -----------------
    for key in sites:
        sites[key]["key"] = key

    return sites


def keys():
    return list(_all().keys())


def init(key):
    return _all()[key]


def assert_valid(site):
    assert (
        site["observation_level_asl_m"] >= -500
    )  # already unreasonable, but hey!
    assert site["corsika_atmosphere_id"] >= 0
    assert not np.isnan(float(site["earth_magnetic_field_x_muT"]))
    assert not np.isnan(float(site["earth_magnetic_field_z_muT"]))
