import numpy as np
import corsika_primary
import magnetic_deflection

from .. import sites
from .. import particles
from .. import pointing_range
from .. import pointing


def draw_event(
    run_id,
    site_particle_magnetic_deflection,
    instrument_pointing_range,
    instrument_field_of_view_half_angle_rad,
    num_events,
):
    """
    Draw the random distribution of particles to induce showers and emitt
    Cherenkov-light which is relevant for our instrument.

    Parameters
    ----------
    run_id : int
        The run-number/run-id of the corsika-run. Must be > 0.
    site_particle_magnetic_deflection : magnetic_deflection.allsky.Random()
        Describes from what direction the given particle must be thrown in
        order to see its Cherenkov-light. Must match 'particle' and 'site'.
    num_events : int
        The number of events in the run.

    Returns
    -------
    steering_dict : dict
        To be given to CORSIKA-primary.
        Describes explicitly how each particle shall be thrown in CORSIKA.
    """

    # assertion checks
    # ----------------
    assert run_id > 0
    assert num_events > 0
    assert instrument_field_of_view_half_angle_rad > 0.0
    prng = np.random.Generator(np.random.PCG64(run_id))

    site = site_particle_magnetic_deflection.config["site"]
    sites.assert_valid(site)
    particle = site_particle_magnetic_deflection.config["particle"]
    particles.assert_valid(particle)

    # energies
    # --------
    start_energy_GeV = particles.compile_energy(
        particle["population"]["energy"]["start_GeV"]
    )
    stop_energy_GeV = particles.compile_energy(
        particle["population"]["energy"]["stop_GeV"]
    )
    energies_GeV = corsika_primary.random.distributions.draw_power_law(
        prng=prng,
        lower_limit=start_energy_GeV,
        upper_limit=stop_energy_GeV,
        power_slope=particle["population"]["energy"]["power_law_slope"],
        num_samples=num_events,
    )

    # instrument pointings
    # --------------------
    instrument_pointings = []
    for i in range(num_events):
        instrument_pointing = pointing_range.draw_pointing(
            pointing_range=instrument_pointing_range,
            prng=prng,
        )
        instrument_pointings.append(instrument_pointing)

    # primary directions
    # ------------------
    rnd = magnetic_deflection.allsky.random.Random(
        allsky_deflection=site_particle_magnetic_deflection
    )
    _shower_spread_half_angle_rad = np.deg2rad(
        particle["population"]["direction"]["scatter_cone_half_angle_deg"]
    )
    primary_directions = []
    for i in range(num_events):
        res, dbg = rnd.draw_particle_direction(
            prng=prng,
            method="grid",
            azimuth_rad=instrument_pointings[i]["azimuth_rad"],
            zenith_rad=instrument_pointings[i]["zenith_rad"],
            half_angle_rad=instrument_field_of_view_half_angle_rad,
            energy_GeV=energies_GeV[i],
            shower_spread_half_angle_rad=_shower_spread_half_angle_rad,
            min_num_cherenkov_photons=1e3,
        )
        primary_directions.append(res)

    i8 = np.int64
    f8 = np.float64

    run = {
        "run_id": i8(run_id),
        "event_id_of_first_event": i8(1),
        "observation_level_asl_m": f8(site["observation_level_asl_m"]),
        "earth_magnetic_field_x_muT": f8(site["earth_magnetic_field_x_muT"]),
        "earth_magnetic_field_z_muT": f8(site["earth_magnetic_field_z_muT"]),
        "atmosphere_id": i8(site["corsika_atmosphere_id"]),
        "energy_range": {
            "start_GeV": f8(start_energy_GeV),
            "stop_GeV": f8(stop_energy_GeV),
        },
        "random_seed": corsika_primary.random.seed.make_simple_seed(run_id),
    }

    primaries = []
    for e in range(num_events):
        prm = {}
        prm["particle_id"] = f8(particle["corsika_particle_id"])
        prm["energy_GeV"] = f8(energies_GeV[e])
        prm["azimuth_rad"] = f8(primary_directions[e]["particle_azimuth_rad"])
        prm["zenith_rad"] = f8(primary_directions[e]["particle_zenith_rad"])
        prm["depth_g_per_cm2"] = f8(0.0)
        primaries.append(prm)

    corsika_primary_steering = {"run": run, "primaries": primaries}

    out = {}
    out["corsika_primary_steering"] = corsika_primary_steering
    out["primary_directions"] = []
    for x in primary_directions:
        y = {}
        for key in ["cutoff", "solid_angle_thrown_sr"]:
            y[key] = x[key]
        out["primary_directions"].append(y)
    return out