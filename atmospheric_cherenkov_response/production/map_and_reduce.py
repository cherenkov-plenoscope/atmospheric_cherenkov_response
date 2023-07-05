from .. import sites
from .. import particles
import magnetic_deflection
import copy


def make_example_job():
    job = {}

    job["work_dir"] = "./2023-06-15_acr"
    job["site_key"] = "lapalma"
    job["particle_key"] = "electron"
    job["pointing_key"] = "B"

    # deflection will be searched for automatically in the work_dir

    job["production"] = {}
    job["production"]["key"] = "test"
    job["production"]["run_id"] = 13
    job["production"]["num_showers"] = 100

    job["grid"] = {
        "instrument_aperture_outer_diameter",
        "bin_width_overhead",
        "instrument_field_of_view_outer_radius_deg",
        "instrument_pointing_direction",
        "field_of_view_overhead",
        "num_bins_radius",
    }
    job["instrument"] = {}

    job["tmp"] = {}
    job["tmp"]["dir"] = "tmp-production-acr"
    job["tmp"]["keep"] = True

    return job
