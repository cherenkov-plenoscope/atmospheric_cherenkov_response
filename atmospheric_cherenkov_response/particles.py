import numpy as np
import corsika_primary
import binning_utils


def _all():
    par = {}

    ebin_GeV = binning_utils.power10.lower_bin_edge
    stop_GeV = ebin_GeV(decade=3, bin=2, num_bins_per_decade=5)

    stop_GeV = {
        "repr": "binning",
        "decade": 3,
        "bin": 2,
        "num_bins_per_decade": 5,
    }

    par["gamma"] = {
        "corsika_particle_id": 1,
        "electric_charge_qe": 0.0,
        "population": {
            "energy": {
                "start_GeV": {
                    "repr": "binning",
                    "decade": -1,
                    "bin": 2,
                    "num_bins_per_decade": 5,
                },
                "stop_GeV": stop_GeV,
                "power_law_slope": -1.5,
            },
            "direction": {"scatter_cone_half_angle_deg": 3.25},
        },
    }

    par["electron"] = {
        "corsika_particle_id": 3,
        "electric_charge_qe": -1.0,
        "population": {
            "energy": {
                "start_GeV": {
                    "repr": "binning",
                    "decade": -1,
                    "bin": 3,
                    "num_bins_per_decade": 5,
                },
                "stop_GeV": stop_GeV,
                "power_law_slope": -1.5,
            },
            "direction": {"scatter_cone_half_angle_deg": 6.5},
        },
    }

    MIN_PROTON_ENERGY_GEV = 5.0
    par["proton"] = {
        "corsika_particle_id": 14,
        "electric_charge_qe": +1.0,
        "population": {
            "energy": {
                "start_GeV": {
                    "repr": "explicit",
                    "value": MIN_PROTON_ENERGY_GEV,
                },
                "stop_GeV": stop_GeV,
                "power_law_slope": -1.5,
            },
            "direction": {"scatter_cone_half_angle_deg": 18.3},
        },
    }

    MIN_HELIUM_ENERGY_GEV = 10.0
    par["helium"] = {
        "corsika_particle_id": 402,
        "electric_charge_qe": +2.0,
        "population": {
            "energy": {
                "start_GeV": {
                    "repr": "explicit",
                    "value": MIN_HELIUM_ENERGY_GEV,
                },
                "stop_GeV": stop_GeV,
                "power_law_slope": -1.5,
            },
            "direction": {"scatter_cone_half_angle_deg": 18.3},
        },
    }

    # assert all have same keys
    # -------------------------
    for pk in par:
        for key in par[pk]:
            assert (
                key in par["gamma"]
            ), "The key '{:s}' from '{:s}' not in 'gamma'".format(key, pk)
        for key in par["gamma"]:
            assert key in par[pk], "{:s} not in '{:s}'".format(key, pk)

    # put key into dict
    # -----------------
    for key in par:
        par[key]["key"] = key

    return par


def keys():
    return list(_all().keys())


def init(key):
    return _all()[key]


def compile_energy(energy):
    _ecpy = energy.copy()
    _repr = _ecpy.pop("repr")
    if _repr == "explicit":
        return energy["value"]
    elif _repr == "binning":
        return binning_utils.power10.lower_bin_edge(**_ecpy)
    else:
        raise KeyError()


def assert_valid(particle):
    assert particle["corsika_particle_id"] >= 0
    assert (
        particle["population"]["direction"]["scatter_cone_half_angle_deg"]
        >= 0.0
    )
    assert compile_energy(particle["population"]["energy"]["start_GeV"]) > 0.0
    assert compile_energy(particle["population"]["energy"]["stop_GeV"]) > 0.0
