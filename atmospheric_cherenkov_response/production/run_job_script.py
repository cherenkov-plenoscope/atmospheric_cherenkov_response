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
import merlict_development_kit_python

op = os.path


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
    _light_field_sensor_geometry = (
        merlict_development_kit_python.plenoscope_propagator.read_plenoscope_geometry(
            merlict_scenery_path=_scenery_path
        )
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
