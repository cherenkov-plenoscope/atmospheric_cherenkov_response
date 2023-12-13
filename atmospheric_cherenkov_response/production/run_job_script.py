import atmospheric_cherenkov_response as acr
import numpy as np
import corsika_primary as cpw
import tempfile
import os
import sys
import tarfile
import json_utils
import magnetic_deflection as mdfl
import sparse_numeric_table as spt
import solid_angle_utils

op = os.path


def init_pointing_cone():
    out = {
        "type": "cone",
        "azimuth_rad": 0.0,
        "zenith_rad": 0.0,
        "half_angel_rad": np.deg2rad(45),
    }


def draw_corsika_primary_steering(
    run_id,
    site,
    particle,
    site_particle_magnetic_deflection_allsky,
    pointing_cone_azimuth_rad,
    pointing_cone_zenith_rad,
    pointing_cone_half_angel_rad,
    instrument_field_of_view_half_angle_rad,
    num_events,
    prng,
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
    site_particle_deflection : magnetic_deflection.allsky.Random()
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
    assert num_events > 0
    _assert_site(site)
    _assert_particle(particle)
    _assert_deflection(site_particle_deflection)

    # energies
    # --------
    start_energy_GeV = compile_energy(
        particle["population"]["energy"]["start_GeV"]
    )
    stop_energy_GeV = compile_energy(
        particle["population"]["energy"]["stop_GeV"]
    )
    energies_GeV = corsika_primary.random.distributions.draw_power_law(
        prng=prng,
        lower_limit=start_energy_GeV,
        upper_limit=stop_energy_GeV,
        power_slope=particle["population"]["energy"]["power_law_slope"],
        num_samples=num_events,
    )

    # pointings
    # ---------
    assert pointing_cone_zenith_deg >= 0.0
    pointings_azimuth_rad = np.nan * np.ones(num_events)
    pointings_zenith_rad = np.nan * np.ones(num_events)
    for i in range(num_events):
        (
            pointings_azimuth_rad[i],
            pointings_zenith_rad[i],
        ) = corsika_primary.random.distributions.draw_azimuth_zenith_in_viewcone(
            prng=prng,
            azimuth_rad=pointing_cone_azimuth_rad,
            zenith_rad=pointing_cone_zenith_rad,
            min_scatter_opening_angle_rad=0.0,
            max_scatter_opening_angle_rad=pointing_cone_half_angel_rad,
        )

    # primary directions
    # ------------------
    rnd = mdfl.allsky.random.Random(
        allsky_deflection=site_particle_magnetic_deflection_allsky
    )
    for i in range(num_events):
        res, dbg = rnd.draw_particle_direction(
            prng=prng,
            method="cone",
            azimuth_deg=np.rad2deg(pointings_azimuth_rad[i]),
            zenith_deg=np.rad2deg(pointings_zenith_rad[i]),
            half_angle_deg=np.rad2deg(instrument_field_of_view_half_angle_rad),
            energy_GeV=energies_GeV[i],
            shower_spread_half_angle_deg=np.rad2deg(),
            min_num_cherenkov_photons=1e3,
        )

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

        """
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
        """

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

    corsika_primary_steering = {
        "run": run,
        "primaries": primaries,
    }

    thrown = []

    return corsika_primary_steering


def make_cherenkovsize_record(uid, cherenkov_bunches):
    cb = cherenkov_bunches
    ase = {spt.IDX: uid}
    ase["num_bunches"] = cb.shape[0]
    ase["num_photons"] = np.sum(cb[:, cpw.I.BUNCH.BUNCH_SIZE_1])
    return ase


def estimate_cherenkov_bunch_statistics(cherenkov_bunches):
    cb = cherenkov_bunches
    ase = {}
    assert cb.shape[0] > 0
    ase["maximum_asl_m"] = cpw.CM2M * np.median(
        cb[:, cpw.I.BUNCH.EMISSOION_ALTITUDE_ASL_CM]
    )
    ase["wavelength_median_nm"] = np.abs(
        np.median(cb[:, cpw.I.BUNCH.WAVELENGTH_NM])
    )
    ase["cx_median_rad"] = np.median(cb[:, cpw.I.BUNCH.CX_RAD])
    ase["cy_median_rad"] = np.median(cb[:, cpw.I.BUNCH.CY_RAD])
    ase["x_median_m"] = cpw.CM2M * np.median(cb[:, cpw.I.BUNCH.X_CM])
    ase["y_median_m"] = cpw.CM2M * np.median(cb[:, cpw.I.BUNCH.Y_CM])
    ase["bunch_size_median"] = np.median(cb[:, cpw.I.BUNCH.BUNCH_SIZE_1])
    return ase


def init_grid_geometry_from_job(job):
    assert job["plenoscope_pointing"]["zenith_deg"] == 0.0
    assert job["plenoscope_pointing"]["azimuth_deg"] == 0.0
    plenoscope_pointing_direction = np.array([0, 0, 1])  # For now this is fix.

    _scenery_path = op.join(job["plenoscope_scenery_path"], "scenery.json")
    _light_field_sensor_geometry = production.merlict.read_plenoscope_geometry(
        merlict_scenery_path=_scenery_path
    )
    plenoscope_diameter = (
        2.0
        * _light_field_sensor_geometry[
            "expected_imaging_system_aperture_radius"
        ]
    )
    plenoscope_radius = 0.5 * plenoscope_diameter
    plenoscope_field_of_view_radius_deg = (
        0.5 * _light_field_sensor_geometry["max_FoV_diameter_deg"]
    )

    grid_geometry = atmospheric_cherenkov_response.grid.init_geometry(
        instrument_aperture_outer_diameter=plenoscope_diameter,
        bin_width_overhead=job["grid"]["bin_width_overhead"],
        instrument_field_of_view_outer_radius_deg=(
            plenoscope_field_of_view_radius_deg
        ),
        instrument_pointing_direction=plenoscope_pointing_direction,
        field_of_view_overhead=job["grid"]["field_of_view_overhead"],
        num_bins_radius=job["grid"]["num_bins_radius"],
    )
    return grid_geometry


job = {}
job["tmp_dir"] = "~/Desktop/2023-acr-tmp"
job["run_id"] = 1337
job["num_events"] = 100
job["site_key"] = "lapalma"
job["particle_key"] = "electron"
job["magnetic_deflection"] = {}
job["magnetic_deflection"]["method"] = "cone"
job["magnetic_deflection"]["half_angle_deg"] = 6.5


os.makedirs(job["tmp_dir"], exist_ok=True)
prng = np.random.Generator(np.random.PCG64(job["run_id"]))

site = acr.sites.init(job["site_key"])
particle = acr.sites.init(job["particle_key"])

site_particle_magnetic_deflection_allsky = mdfl.allsky.AllSky(
    op.join(
        "~/Desktop/2023-12-08_magsky", job["site_key"], job["particle_key"]
    )
)

corsika_primary_steering = acr.particles.draw_corsika_primary_steering(
    run_id=job["production"]["run_id"],
    site=site,
    particle=particle,
    site_particle_deflection=site_particle_deflection,
    num_events=job["num_events"],
    prng=prng,
)

# loop over air-showers
# ---------------------
cherenkov_pools_path = op.join(tmp_dir, "cherenkov_pools.tar")
particle_pools_dat_path = op.join(tmp_dir, "particle_pools.dat")
particle_pools_tar_path = op.join(tmp_dir, "particle_pools.tar.gz")
grid_histogram_path = op.join(tmp_dir, "grid.tar")
grid_roi_histogram_path = op.join(tmp_dir, "grid_roi.tar")


def read_all_cherenkov_bunches(cherenkov_reader):
    return np.vstack([b for b in cherenkov_reader])


def _is_close(a, b, eps=1e-6):
    return np.abs(a - b) < eps


def make_and_assert_uid_and_uid_str(corsika_evth, run_id, event_id):
    assert _is_close(run_id, corsika_evth[cpw.I.EVTH.RUN_NUMBER])
    assert _is_close(event_id, corsika_evth[cpw.I.EVTH.EVENT_NUMBER])
    uid = unique.make_uid(run_id=run_id, event_id=event_id)
    uid_str = unique.make_uid_str(run_id=run_id, event_id=event_id)


def make_primary_record(uid, corsika_evth, corsika_primary_steering):
    EVTH = cpw.I.EVTH

    run_id, event_id = unique.split_uid(uid=uid)
    event_idx = event_id - 1
    prmst = corsika_primary_steering["primaries"][event_idx]

    # unique id of shower
    prim = {spt.IDX: uid}

    prim["particle_id"] = prmst["particle_id"]
    assert _is_close(prim["particle_id"], corsika_evth[EVTH.PARTICLE_ID])

    prim["energy_GeV"] = prmst["energy_GeV"]
    assert _is_close(prim["energy_GeV"], corsika_evth[EVTH.TOTAL_ENERGY_GEV])

    prim["azimuth_rad"] = prmst["azimuth_rad"]
    assert _is_close(prim["azimuth_rad"], corsika_evth[EVTH.AZIMUTH_RAD])

    prim["zenith_rad"] = prmst["zenith_rad"]
    assert _is_close(prim["zenith_rad"], corsika_evth[EVTH.ZENITH_RAD])

    prim["depth_g_per_cm2"] = prmst["depth_g_per_cm2"]
    assert _is_close(
        prim["depth_g_per_cm2"], corsika_evth[EVTH.STARTING_DEPTH_G_PER_CM2]
    )

    prim["max_scatter_rad"] = prmst["max_scatter_rad"]
    prim["solid_angle_thrown_sr"] = solid_angle_utils.cone.solid_angle(
        half_angle_rad=prim["max_scatter_rad"]
    )

    prim["momentum_x_GeV_per_c"] = corsika_evth[EVTH.PX_MOMENTUM_GEV_PER_C]
    prim["momentum_y_GeV_per_c"] = corsika_evth[EVTH.PY_MOMENTUM_GEV_PER_C]
    prim["momentum_z_GeV_per_c"] = (
        -1.0 * corsika_evth[EVTH.PZ_MOMENTUM_GEV_PER_C]
    )
    prim["first_interaction_height_asl_m"] = (
        -1.0 * cpw.CM2M * corsika_evth[EVTH.Z_FIRST_INTERACTION_CM]
    )
    prim["starting_height_asl_m"] = (
        cpw.CM2M * corsika_evth[EVTH.STARTING_HEIGHT_CM]
    )
    obs_lvl_intersection = acr.utils.ray_plane_x_y_intersection(
        support=[0, 0, prim["starting_height_asl_m"]],
        direction=[
            prim["momentum_x_GeV_per_c"],
            prim["momentum_y_GeV_per_c"],
            prim["momentum_z_GeV_per_c"],
        ],
        plane_z=job["site"]["observation_level_asl_m"],
    )
    prim["starting_x_m"] = -1.0 * obs_lvl_intersection[0]
    prim["starting_y_m"] = -1.0 * obs_lvl_intersection[1]
    prim["magnet_azimuth_rad"] = prmst["magnet_azimuth_rad"]
    prim["magnet_zenith_rad"] = prmst["magnet_zenith_rad"]
    prim["magnet_cherenkov_pool_x_m"] = prmst["magnet_cherenkov_pool_x_m"]
    prim["magnet_cherenkov_pool_y_m"] = prmst["magnet_cherenkov_pool_y_m"]
    return prim


def append_to_tabrec(tabrec, key, record):
    if key not in tabrec:
        tabrec[key] = []
    tabrec[key].append(record)


tabrec = {}

with cpw.cherenkov.CherenkovEventTapeWriter(
    path=cherenkov_pools_path
) as evttar, tarfile.open(grid_histogram_path, "w") as imgtar, tarfile.open(
    grid_roi_histogram_path, "w"
) as imgroitar:
    with cpw.CorsikaPrimary(
        corsika_path=config["executables"]["corsika_primary"]["path"],
        steering_dict=corsika_primary_steering,
        stdout_path=op.join(tmp_dir, "corsika.stdout"),
        stderr_path=op.join(tmp_dir, "corsika.stderr"),
        particle_output_path=particle_pools_dat_path,
    ) as corsika_run:
        evttar.write_runh(runh=corsika_run.runh)

        for event_idx, corsika_event in enumerate(corsika_run):
            corsika_evth, cherenkov_reader = corsika_event

            cherenkov_bunches = read_all_cherenkov_bunches(
                cherenkov_reader=cherenkov_reader
            )

            uid, uid_str = make_and_assert_uid_and_uid_str(
                corsika_evth=corsika_evth,
                event_id=event_idx + 1,
                run_id=corsika_primary_steering["run"]["run_id"],
            )

            tabrec = append_to_tabrec(
                tabrec=tabrec,
                key="primary",
                record=make_primary_record(
                    uid=uid,
                    corsika_evth=corsika_evth,
                    corsika_primary_steering=corsika_primary_steering,
                ),
            )

            tabrec = append_to_tabrec(
                tabrec=tabrec,
                key="cherenkovsize",
                record=make_cherenkovsize_record(
                    uid=uid,
                    cherenkov_bunches=cherenkov_bunches,
                ),
            )

            (
                grid_random_shift_x,
                grid_random_shift_y,
            ) = acr.grid.draw_random_shift_x_y(
                grid_geometry=grid_geometry, prng=prng
            )

            print("size", len(cherenkov_bunches))
