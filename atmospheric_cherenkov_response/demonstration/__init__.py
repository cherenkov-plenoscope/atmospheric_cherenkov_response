import numpy as np
import corsika_primary
import json_numpy
import os
import magnetic_deflection
import json_line_logger


GRID = {
    "num_bins_radius": 512,
    "threshold_num_photons": 50,
    "field_of_view_overhead": 1.1,
    "bin_width_overhead": 1.1,
    "output_after_num_events": 25,
}


def init(work_dir, config=None):
    if config == None:
        config = {
            "corsika": magnetic_deflection.examples.CORSIKA_PRIMARY_MOD_PATH,
            "particles": magnetic_deflection.examples.PARTICLES,
            "sites": magnetic_deflection.examples.SITES,
            "pointings": [
                {"azimuth_deg": 0.0, "zenith_deg": 0.0},
                {"azimuth_deg": 0.0, "zenith_deg": 22.5},
                {"azimuth_deg": 0.0, "zenith_deg": 45.0},
            ],
            "magnetic_deflection": {
                "num_energy_supports": 512,
                "max_energy_GeV": 64,
            },
            "runs": {
                "gamma": {"num": 64, "first_run_id": 1},
                "electron": {"num": 64, "first_run_id": 1},
                "proton": {"num": 64, "first_run_id": 1},
                "helium": {"num": 64, "first_run_id": 1},
            },
            "num_airshowers_per_run": 100,
            "grid": GRID,
            "classic": {},
        }

    with open(os.path.join(work_dir, "config.json"), "wt") as f:
        f.write(json_numpy.dumps(config, indent=4))


def run(work_dir, pool, logger=json_line_logger.LoggerStdout()):
    config = json_numpy.read(os.path.join(work_dir, "config.json"))

    # magnetic deflections
    # --------------------

    jobs = []
    for pointing in config["pointings"]:
        pointing_key = make_pointing_key(pointing)
        pointing_dir = os.path.join(work_dir, "map", pointing_key)

        if not os.path.exists(pointing_dir):
            logger.info("Adding magnetic-deflection for " + pointing_key)

            magnetic_deflection.init(
                work_dir=pointing_dir,
                particles=config["particles"],
                sites=config["sites"],
                pointing=pointing,
                max_energy=config["magnetic_deflection"]["max_energy_GeV"],
                num_energy_supports=config["magnetic_deflection"][
                    "num_energy_supports"
                ],
            )

            jobs += magnetic_deflection.make_jobs(work_dir=pointing_dir)

    if len(jobs):
        logger.info(
            "Running {:d} jobs for magnetic-deflection".format(len(jobs))
        )

    _ = pool.map(magnetic_deflection.map_and_reduce.run_job, jobs)

    for pointing in config["pointings"]:
        pointing_key = make_pointing_key(pointing)
        pointing_dir = os.path.join(work_dir, "map", pointing_key)

        try:
            _ = magnetic_deflection.read_deflection(
                work_dir=pointing_dir, style="dict",
            )
        except:
            logger.info("Reducing magnetic-deflection " + pointing_key)
            magnetic_deflection.reduce(work_dir=pointing_dir, logger=logger)


def make_pointing_key(pointing):
    return "az{:06d}mdeg_zd{:06d}mdeg".format(
        int(pointing["azimuth_deg"] * 1e3), int(pointing["zenith_deg"] * 1e3),
    )
