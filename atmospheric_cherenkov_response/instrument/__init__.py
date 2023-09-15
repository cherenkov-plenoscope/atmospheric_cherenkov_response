import numpy as np
from . import toy
from . import viewcone
from . import bounding_box


EXAMPLE = {
    "viewcone": viewcone.EXAMPLE,
    "bounding_box": bounding_box.EXAMPLE,
}


def intersection_of_ray_on_ground(support, direction):
    """
    Returns the intersection-vector of a ray and the ground.
    Ground is the x-y-plane where z = 0.

    Parameters
    ----------
    support : array 3
        Ray's support-vector
    direction : array 3
        Ray's direction-vector
    """
    direction = direction / np.linalg.norm(direction)
    lam = -support[2] / direction[2]
    return support + lam * direction


def estimate_projection_on_ground(instrument, num_rays=1000):
    """
    Returns the outer bounds ((xmin, xmax), (ymin, ymax)) of an instrument on
    ground. Ground is the x-y-plane where z = 0.
    The outer bounds mark the area where the ground collects photons which
    can potentially be seen by the instrument.

    Parameters
    ----------
    instrument : dict
        Contains the instrument's bounding box and viewcone.
    num_rays : int, default = 1000
        The number of rays emitted to project the instrument's viewcone.
    """
    xmin = -0.5 * instrument["bounding_box"]["width_x"]
    xmax = +0.5 * instrument["bounding_box"]["width_x"]
    ymin = -0.5 * instrument["bounding_box"]["width_y"]
    ymax = +0.5 * instrument["bounding_box"]["width_y"]

    for ix in [0, 1]:
        for iy in [0, 1]:
            corner = bounding_box.get_corner(
                bounding_box=instrument["bounding_box"], ix=ix, iy=iy, iz=1
            )

            for viewcone_azimuth_deg in np.linspace(0, 360, num_rays):
                viewing_direction = viewcone.make_direction_on_edge_of_viewcone_when_viewcone_pointing(
                    pointing_azimuth_deg=instrument["viewcone"]["azimuth_deg"],
                    pointing_zenith_deg=instrument["viewcone"]["zenith_deg"],
                    viewcone_azimuth_deg=viewcone_azimuth_deg,
                    viewcone_half_angle_deg=instrument["viewcone"][
                        "half_angle_deg"
                    ],
                )

                incoming_direction = -1.0 * viewing_direction

                isec = intersection_of_ray_on_ground(
                    support=corner,
                    direction=incoming_direction,
                )
                px = isec[0]
                py = isec[1]

                xmin = np.min([xmin, px])
                xmax = np.max([xmax, px])

                ymin = np.min([ymin, py])
                ymax = np.max([ymax, py])

    return ((xmin, xmax), (ymin, ymax))
