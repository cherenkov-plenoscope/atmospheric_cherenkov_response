import numpy as np
from . import bounding_box
from . import viewcone
from . import instrument


def project_bounding_box_onto_ground(
    bounding_box_width_x,
    bounding_box_width_y,
    bounding_box_height_z,
    viewcone_half_angle_deg,
    viewcone_azimuth_deg,
    viewcone_zenith_distance_deg,
):
    bb = np.array([bounding_box_x, bounding_box_y, bounding_box_z,])

    """
    project all eight corners of the instrument's bounding box
    on the ground (x-y-plane)
    """
