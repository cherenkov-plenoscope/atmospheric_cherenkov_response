import numpy as np
import corsika_primary
import binning_utils


def scatter_cone(key):
    out = {}
    out["gamma"] = {
        "half_angle_rad": np.deg2rad([3.25, 3.25]),
        "energy_GeV": [1e-1, 1e4],
    }
    out["electron"] = {
        "half_angle_rad": np.deg2rad([6.5, 6.5]),
        "energy_GeV": [1e-1, 1e4],
    }
    out["proton"] = {
        "half_angle_rad": np.deg2rad([15.0, 15.0]),
        "energy_GeV": [1e0, 1e4],
    }
    out["helium"] = {
        "half_angle_rad": np.deg2rad([15.0, 15.0]),
        "energy_GeV": [1e0, 1e4],
    }
    return out[key]


def energy_stop_GeV():
    stop_GeV = {
        "repr": "binning",
        "decade": 3,
        "bin": 2,
        "num_bins_per_decade": 5,
    }
    return stop_GeV


def _all():
    par = {}

    par["gamma"] = {
        "corsika": {"particle_id": 1, "min_energy_GeV": None},
        "electric_charge_qe": 0.0,
    }

    par["electron"] = {
        "corsika": {"particle_id": 3, "min_energy_GeV": None},
        "electric_charge_qe": -1.0,
    }

    par["proton"] = {
        "corsika": {"particle_id": 14, "min_energy_GeV": 5.0},
        "electric_charge_qe": +1.0,
    }

    par["helium"] = {
        "corsika": {"particle_id": 402, "min_energy_GeV": 10.0},
        "electric_charge_qe": +2.0,
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


"""
def get_scatter_cone_half_angle_rad(particle, energy_GeV):
    sc = particle["population"]["direction"]["scatter_cone"]
    return interpolate_scatter_cone_half_angle(
        energy_GeV=energy_GeV,
        scatter_cone_energy_GeV=sc["energy_GeV"],
        scatter_cone_half_angle_rad=sc["half_angle_rad"],
    )


def get_energy_start_GeV(particle):
    return compile_energy(energy=particle["population"]["energy"]["start_GeV"])


def get_energy_stop_GeV(particle):
    return compile_energy(energy=particle["population"]["energy"]["stop_GeV"])


def compile_energy(energy):
    _ecpy = energy.copy()
    _repr = _ecpy.pop("repr")
    if _repr == "explicit":
        return energy["value"]
    elif _repr == "binning":
        return binning_utils.power10.lower_bin_edge(**_ecpy)
    else:
        raise KeyError()
"""


def interpolate_scatter_solid_angle(
    energy_GeV,
    scatter_energy_GeV,
    scatter_solid_angle_sr,
):
    energies_GeV = np.asarray(scatter_energy_GeV)
    solid_angles_sr = np.asarray(scatter_solid_angle_sr)

    assert binning_utils.is_strictly_monotonic_increasing(energies_GeV)
    assert np.all(energies_GeV > 0.0)
    assert energy_GeV > 0.0
    assert np.all(solid_angles_sr >= 0.0)
    assert np.all(solid_angles_sr <= 4 * np.pi)
    assert min(energies_GeV) <= energy_GeV <= max(energies_GeV)

    solid_angle_sr = np.interp(
        x=energy_GeV, xp=energies_GeV, fp=solid_angles_sr
    )
    return solid_angle_sr


def assert_valid(particle):
    assert particle["corsika"]["particle_id"] >= 0
