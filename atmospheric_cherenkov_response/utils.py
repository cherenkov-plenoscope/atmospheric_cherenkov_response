import numpy as np


def ray_plane_x_y_intersection(support, direction, plane_z):
    direction = np.array(direction)
    support = np.array(support)
    direction_norm = direction / np.linalg.norm(direction)
    ray_parameter = -(support[2] - plane_z) / direction_norm[2]
    intersection = support + ray_parameter * direction_norm
    assert np.abs(intersection[2] - plane_z) < 1e-3
    return intersection
