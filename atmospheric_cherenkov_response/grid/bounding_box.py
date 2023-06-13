import numpy as np


EXAMPLE = {
    "width_x": 1.0,
    "width_y": 1.0,
    "height_z": 1.0,
}


def init(width_x, width_y, height_z):
    assert width_x > 0
    assert width_y > 0
    assert height_z >= 0
    return {
        "width_x": width_x,
        "width_y": width_y,
        "height_z": height_z,
    }


def get_corner(bounding_box, ix=0, iy=0, iz=0):
    dx = bounding_box["width_x"] / 2.0
    dy = bounding_box["width_y"] / 2.0
    hz = bounding_box["height_z"]
    return np.array([dx if ix else -dx, dy if iy else -dy, hz if iz else 0.0])
