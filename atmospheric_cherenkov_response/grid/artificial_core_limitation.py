import numpy as np
import json_numpy


# artificial core limitation
# --------------------------

def where_grid_idxs_within_radius(grid_geometry, radius, center_x, center_y):
    """
    Returns the idxs in x, and y of the grid where the bin is within a circle
    of radius and center_x/y.
    Same format as np.where()
    """
    pgg = grid_geometry
    radius_bins = int(np.ceil(1.1 * radius / pgg["bin_width"]))

    x_bin_center = int(np.round(center_x / pgg["bin_width"]))
    x_bin_range = (
        x_bin_center
        + pgg["num_bins_radius"]
        + np.arange(-radius_bins, radius_bins + 1)
    )
    y_bin_center = int(np.round(center_y / pgg["bin_width"]))
    y_bin_range = (
        y_bin_center
        + pgg["num_bins_radius"]
        + np.arange(-radius_bins, radius_bins + 1)
    )

    num_bins_thrown = 0
    grid_idxs_x = []
    grid_idxs_y = []
    for x_bin in x_bin_range:
        for y_bin in y_bin_range:
            if x_bin < 0 or x_bin >= pgg["num_bins_diameter"]:
                continue
            if y_bin < 0 or y_bin >= pgg["num_bins_diameter"]:
                continue
            delta_x = pgg["xy_bin_centers"][x_bin] - center_x
            delta_y = pgg["xy_bin_centers"][y_bin] - center_y
            if delta_x ** 2 + delta_y ** 2 > radius ** 2:
                continue
            grid_idxs_x.append(x_bin)
            grid_idxs_y.append(y_bin)

    return np.array(grid_idxs_x), np.array(grid_idxs_y)


def intersection_of_bin_idxs(a_bin_idxs, b_bin_idxs):
    """
    Returns only bin-idxs which are both in a_bins and b_bins.
    """
    a = a_bin_idxs
    b = b_bin_idxs
    assert len(a[0]) == len(a[1])
    assert len(b[0]) == len(b[1])
    a_set = set()
    for ia in range(len(a[0])):
        a_set.add((a[0][ia], a[1][ia]))
    b_set = set()
    for ib in range(len(b[0])):
        b_set.add((b[0][ib], b[1][ib]))
    ab_set = a_set.intersection(b_set)
    ab = list(ab_set)
    x_idxs = np.array([idx[0] for idx in ab])
    y_idxs = np.array([idx[1] for idx in ab])
    return (x_idxs, y_idxs)


def apply_bin_limitation_and_warn(
    bin_idxs_above_threshold, bin_idxs_limitation
):
    bin_idxs_above_threshold_and_in_limits = intersection_of_bin_idxs(
        a_bin_idxs=bin_idxs_above_threshold, b_bin_idxs=bin_idxs_limitation
    )
    msg = {}
    msg["artificial_core_limitation"] = {
        "num_grid_bins_in_limits": len(bin_idxs_limitation[0]),
        "num_grid_bins_above_threshold": len(bin_idxs_above_threshold[0]),
        "num_grid_bins_above_threshold_and_in_limits": len(
            bin_idxs_above_threshold_and_in_limits[0]
        ),
    }
    print(json_numpy.dumps(msg))

    return bin_idxs_above_threshold_and_in_limits
