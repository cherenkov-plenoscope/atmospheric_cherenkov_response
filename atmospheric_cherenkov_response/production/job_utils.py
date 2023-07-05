from .. import sites
from .. import particles
import copy


def make_example_site():
    job = {}
    job["site"] = copy.deepcopy(sites.init["lapalma"])
    return job


def make_job(
    production_key,
    run_dir,
    run_id,
    site_key,
    particle_key,
    config,
    deflection_table,
    num_air_showers,
    corsika_primary_path,
    merlict_plenoscope_propagator_path,
    tmp_dir,
    keep_tmp_dir,
    date_dict,
):
    job = {
        "run_id": run_id,
        "production_key": production_key,
        "site_key": site_key,
        "particle_key": particle_key,
        "num_air_showers": num_air_showers,
        "plenoscope_pointing": config["plenoscope_pointing"],
        "particle": config["particles"][particle_key],
        "site": config["sites"][site_key],
        "site_particle_deflection": deflection_table[site_key][particle_key],
        "grid": config["grid"],
        "raw_sensor_response": config["raw_sensor_response"],
        "sum_trigger": config["sum_trigger"],
        "cherenkov_classification": config["cherenkov_classification"],
        "reconstruction": config["reconstruction"],
        "corsika_primary_path": str(corsika_primary_path),
        "merlict_plenoscope_propagator_path": str(
            merlict_plenoscope_propagator_path
        ),
        "plenoscope_scenery_path": os.path.join(
            run_dir, "light_field_geometry", "input", "scenery"
        ),
        "light_field_geometry_path": os.path.join(
            run_dir, "light_field_geometry"
        ),
        "trigger_geometry_path": os.path.join(run_dir, "trigger_geometry"),
        "merlict_plenoscope_propagator_config_path": os.path.join(
            run_dir, "input", "merlict_propagation_config.json"
        ),
        "log_dir": os.path.join(
            run_dir, production_key, site_key, particle_key, "log.map"
        ),
        "past_trigger_dir": os.path.join(
            run_dir, production_key, site_key, particle_key, "past_trigger.map"
        ),
        "particles_dir": os.path.join(
            run_dir, production_key, site_key, particle_key, "particles.map"
        ),
        "past_trigger_reconstructed_cherenkov_dir": os.path.join(
            run_dir,
            production_key,
            site_key,
            particle_key,
            "past_trigger_reconstructed_cherenkov_dir.map",
        ),
        "feature_dir": os.path.join(
            run_dir, production_key, site_key, particle_key, "features.map"
        ),
        "keep_tmp": keep_tmp_dir,
        "tmp_dir": tmp_dir,
        "date": date_dict,
        "artificial_core_limitation": config["artificial_core_limitation"][
            particle_key
        ],
    }
    return job
