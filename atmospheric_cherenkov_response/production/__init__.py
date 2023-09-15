import numpy as np
import corsika_primary
import json_utils
import os
import magnetic_deflection
import json_line_logger
import rename_after_writing as rnw
from . import map_and_reduce
from .. import pointing
from .. import sites
from .. import particles
from .. import grid
from . import table


def init(work_dir, config=None):
    join = os.path.join
    if config == None:
        cfg_dir = join(work_dir, "config")
        os.makedirs(cfg_dir, exist_ok=True)

        cfg_exe_dir = join(cfg_dir, "executables")
        os.makedirs(cfg_exe_dir, exist_ok=True)

        with rnw.open(join(cfg_exe_dir, "corsika_primary.json"), "wt") as f:
            f.write(
                json_utils.dumps(
                    {
                        "path": magnetic_deflection.examples.CORSIKA_PRIMARY_MOD_PATH
                    },
                    indent=4,
                )
            )

        with rnw.open(join(cfg_dir, "sites.json"), "wt") as f:
            f.write(json_utils.dumps(sites._all(), indent=4))

        with rnw.open(join(cfg_dir, "particles.json"), "wt") as f:
            f.write(json_utils.dumps(particles._all(), indent=4))

        with rnw.open(join(cfg_dir, "grid.json"), "wt") as f:
            f.write(json_utils.dumps(grid.EXAMPLE, indent=4))

        cfg_prd_dir = join(cfg_dir, "production")
        os.makedirs(cfg_prd_dir, exist_ok=True)

        with rnw.open(
            join(cfg_prd_dir, "magnetic_deflection.json"), "wt"
        ) as f:
            f.write(
                json_utils.dumps(
                    {
                        "num_energy_supports": 512,
                        "max_energy_GeV": 64,
                    },
                    indent=4,
                )
            )

        with rnw.open(
            join(cfg_prd_dir, "instrument_response.json"), "wt"
        ) as f:
            f.write(
                json_utils.dumps(
                    {
                        "num_airshowers_per_run": 100,
                        "first_run_id": 1,
                        "export_details": {"grid": {"probability": 0.01}},
                    },
                    indent=4,
                )
            )

        with rnw.open(join(cfg_dir, "pointings.json"), "wt") as f:
            f.write(
                json_utils.dumps(
                    {
                        "A": pointing.init(azimuth_deg=0.0, zenith_deg=0.0),
                        "B": pointing.init(azimuth_deg=0.0, zenith_deg=22.5),
                        "C": pointing.init(azimuth_deg=0.0, zenith_deg=45.0),
                    },
                    indent=4,
                )
            )

        with rnw.open(join(cfg_dir, "toy.json"), "wt") as f:
            f.write(
                json_utils.dumps(
                    {
                        "sites": sites.keys(),
                        "particles": particles.keys(),
                        "pointings": ["A", "B", "C"],
                    },
                    indent=4,
                )
            )


def run_magnetic_deflection(
    work_dir, pool, logger=json_line_logger.LoggerStdout()
):
    join = os.path.join
    config = json_utils.tree.read(join(work_dir, "config"))

    mdf_dir = join(work_dir, "magnetic_deflection")

    jobs = []
    for ptg in config["pointings"]:
        ptg_dir = join(mdf_dir, ptg)

        if not os.path.exists(ptg_dir):
            logger.info("Adding magnetic-deflection for pointing ", ptg)

            magnetic_deflection.init(
                work_dir=ptg_dir,
                particles=config["particles"],
                sites=config["sites"],
                pointing=config["pointings"][ptg],
                max_energy=config["production"]["magnetic_deflection"][
                    "max_energy_GeV"
                ],
                num_energy_supports=config["production"][
                    "magnetic_deflection"
                ]["num_energy_supports"],
            )

            jobs += magnetic_deflection.make_jobs(work_dir=ptg_dir)

    if len(jobs):
        logger.info(
            "Running {:d} jobs for magnetic-deflection".format(len(jobs))
        )

    _ = pool.map(magnetic_deflection.map_and_reduce.run_job, jobs)

    for ptg in config["pointings"]:
        ptg_dir = os.path.join(mdf_dir, ptg)

        try:
            _ = magnetic_deflection.read_deflection(
                work_dir=ptg_dir,
                style="dict",
            )
        except:
            logger.info("Reducing magnetic-deflection for pointing " + ptg)
            magnetic_deflection.reduce(work_dir=ptg_dir, logger=logger)


def run(work_dir, pool, logger=json_line_logger.LoggerStdout()):
    run_magnetic_deflection(work_dir=work_dir, pool=pool, logger=logger)
