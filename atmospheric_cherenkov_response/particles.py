import numpy as np
import corsika_primary


def _assert_deflection(site_particle_deflection):
    spd = site_particle_deflection

    required_keys = [
        "particle_energy_GeV",
        "particle_azimuth_deg",
        "particle_zenith_deg",
        "cherenkov_x_m",
        "cherenkov_y_m",
    ]

    for required_key in required_keys:
        assert required_key in spd

    assert len(spd["particle_energy_GeV"]) >= 2
    for _key in spd:
        assert len(spd["particle_energy_GeV"]) == len(spd[_key])
    for energy in spd["particle_energy_GeV"]:
        assert energy > 0.0
    for zenith_deg in spd["particle_zenith_deg"]:
        assert zenith_deg >= 0.0
    assert np.all(np.diff(spd["particle_energy_GeV"]) >= 0)


def _assert_site(site):
    required_keys = [
        "observation_level_asl_m",
        "earth_magnetic_field_x_muT",
        "earth_magnetic_field_z_muT",
        "atmosphere_id",
    ]
    for required_key in required_keys:
        assert required_key in site


def _assert_particle(particle):
    required_keys = [
        "particle_id",
        "energy_bin_edges_GeV",
        "max_scatter_angle_deg",
        "energy_power_law_slope",
    ]
    for required_key in required_keys:
        assert required_key in particle

    assert np.all(np.diff(particle["energy_bin_edges_GeV"]) >= 0)
    assert len(particle["energy_bin_edges_GeV"]) == 2


def draw_corsika_primary_steering(
    run_id, site, particle, site_particle_deflection, num_events, prng,
):
    """
    Draw the random distribution of particles to induce showers and emitt
    Cherenkov-light which is relevant for our instrument.

    Parameters
    ----------
    run_id : int
        The run-number/run-id of the corsika-run. Must be > 0.
    site : dict
        Describes the site's magnetic field, altitude and atmosphere.
    particle : dict
        Describes the particle's type and distribution in energy and
        solid angle.
    site_particle_deflection : dict
        Describes from what direction the given particle must be thrown in
        order to see its Cherenkov-light. Must match 'particle' and 'site'.
    num_events : int
        The number of events in the run.
    prng : numpy.random.Generator
        Instance to draw random numbers from.

    Returns
    -------
    steering_dict : dict
        To be given to CORSIKA-primary.
        Describes explicitly how each particle shall be thrown in CORSIKA.
    """

    # Assertions shall prevent 'late' errors in large paralel productions on
    # cluster-computers which are more difficult to debug
    assert run_id > 0
    _assert_site(site)
    _assert_particle(particle)
    _assert_deflection(site_particle_deflection)
    assert num_events > 0

    max_scatter_rad = np.deg2rad(particle["max_scatter_angle_deg"])

    start_energy_GeV = np.max(
        [
            np.min(particle["energy_bin_edges_GeV"]),
            np.min(site_particle_deflection["particle_energy_GeV"]),
        ]
    )
    stop_energy_GeV = np.max(particle["energy_bin_edges_GeV"])

    energies_GeV = corsika_primary.random.distributions.draw_power_law(
        prng=prng,
        lower_limit=start_energy_GeV,
        upper_limit=stop_energy_GeV,
        power_slope=particle["energy_power_law_slope"],
        num_samples=num_events,
    )

    i8 = np.int64
    f8 = np.float64

    run = {
        "run_id": i8(run_id),
        "event_id_of_first_event": i8(1),
        "observation_level_asl_m": f8(site["observation_level_asl_m"]),
        "earth_magnetic_field_x_muT": f8(site["earth_magnetic_field_x_muT"]),
        "earth_magnetic_field_z_muT": f8(site["earth_magnetic_field_z_muT"]),
        "atmosphere_id": i8(site["atmosphere_id"]),
        "energy_range": {
            "start_GeV": f8(start_energy_GeV),
            "stop_GeV": f8(stop_energy_GeV),
        },
        "random_seed": corsika_primary.random.seed.make_simple_seed(run_id),
    }

    primaries = []
    for e in range(num_events):
        prm = {}
        prm["particle_id"] = f8(particle["particle_id"])
        prm["energy_GeV"] = f8(energies_GeV[e])
        prm["magnet_azimuth_rad"] = np.deg2rad(
            np.interp(
                x=prm["energy_GeV"],
                xp=site_particle_deflection["particle_energy_GeV"],
                fp=site_particle_deflection["particle_azimuth_deg"],
            )
        )
        prm["magnet_zenith_rad"] = np.deg2rad(
            np.interp(
                x=prm["energy_GeV"],
                xp=site_particle_deflection["particle_energy_GeV"],
                fp=site_particle_deflection["particle_zenith_deg"],
            )
        )
        prm["magnet_cherenkov_pool_x_m"] = np.interp(
            x=prm["energy_GeV"],
            xp=site_particle_deflection["particle_energy_GeV"],
            fp=site_particle_deflection["cherenkov_x_m"],
        )
        prm["magnet_cherenkov_pool_y_m"] = np.interp(
            x=prm["energy_GeV"],
            xp=site_particle_deflection["particle_energy_GeV"],
            fp=site_particle_deflection["cherenkov_y_m"],
        )
        (
            az,
            zd,
        ) = corsika_primary.random.distributions.draw_azimuth_zenith_in_viewcone(
            prng=prng,
            azimuth_rad=prm["magnet_azimuth_rad"],
            zenith_rad=prm["magnet_zenith_rad"],
            min_scatter_opening_angle_rad=0.0,
            max_scatter_opening_angle_rad=max_scatter_rad,
        )
        prm["max_scatter_rad"] = f8(max_scatter_rad)
        prm["zenith_rad"] = f8(zd)
        prm["azimuth_rad"] = f8(az)
        prm["depth_g_per_cm2"] = f8(0.0)
        primaries.append(prm)

    return {
        "run": run,
        "primaries": primaries,
    }
