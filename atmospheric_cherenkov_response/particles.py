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
            "direction": {
                "scatter_cone": {
                    "half_angle_rad": np.deg2rad([3.25, 3.25]),
                    "energy_GeV": [1e-1, 1e4],
                },
            },
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
            "direction": {
                "scatter_cone": {
                    "half_angle_rad": np.deg2rad([6.5, 6.5]),
                    "energy_GeV": [1e-1, 1e4],
                },
            },
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
            "direction": {
                "scatter_cone": {
                    "half_angle_rad": np.deg2rad([18.3, 18.3]),
                    "energy_GeV": [1e0, 1e4],
                },
            },
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
            "direction": {
                "scatter_cone": {
                    "half_angle_rad": np.deg2rad([18.3, 18.3]),
                    "energy_GeV": [1e0, 1e4],
                },
            },
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


def interpolate_scatter_cone_half_angle(
    energy_GeV,
    scatter_cone_energy_GeV,
    scatter_cone_half_angle_rad,
):
    energies_GeV = np.asarray(scatter_cone_energy_GeV)
    half_angles_rad = np.asarray(scatter_cone_half_angle_rad)

    assert binning_utils.is_strictly_monotonic_increasing(energies_GeV)
    assert np.all(energies_GeV > 0.0)
    assert energy_GeV > 0.0
    assert np.all(half_angles_rad >= 0.0)
    assert np.all(half_angles_rad <= np.pi)
    assert min(energies_GeV) <= energy_GeV <= max(energies_GeV)

    half_angle_rad = np.interp(
        x=energy_GeV, xp=energies_GeV, fp=half_angles_rad
    )
    return half_angle_rad


def assert_valid(particle):
    assert particle["corsika_particle_id"] >= 0
    e_start = get_energy_start_GeV(particle=particle)
    e_stop = get_energy_stop_GeV(particle=particle)
    assert e_start < e_stop

    sc_e_start = get_scatter_cone_half_angle_rad(
        particle=particle, energy_GeV=e_start
    )
    sc_e_stop = get_scatter_cone_half_angle_rad(
        particle=particle, energy_GeV=e_stop
    )

    assert sc_e_start >= 0.0
    assert sc_e_stop >= 0.0
