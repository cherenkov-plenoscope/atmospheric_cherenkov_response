from .. import sites
from .. import particles
import magnetic_deflection
import copy


def make_example_job():
    job = {}
    job["site"] = copy.deepcopy(sites.init("lapalma"))
    job["particle"] = copy.deepcopy(particles.init("electron"))

    # deflection will be searched for automatically in the work_dir

    job["production"] = {}
    job["production"]["key"] = "test"
    job["production"]["run_id"] = 13
    job["production"]["num_showers"] = 100

    job["executables"] = {}
    job["executables"][
        "corsika_primary_path"
    ] = magnetic_deflection.examples.CORSIKA_PRIMARY_MOD_PATH

    job["grid"] = {}
    job["instrument"] = {}

    return job
