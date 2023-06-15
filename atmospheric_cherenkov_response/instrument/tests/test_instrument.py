import atmospheric_cherenkov_response as acr
import numpy as np


def test_projection_no_cone_straight_dwon():
    ins = {
        "bounding_box": {"width_x": 100, "width_y": 100, "height_z": 100},
        "viewcone": {"azimuth_deg": 0, "zenith_deg": 0, "half_angle_deg": 0},
    }
    xlim, ylim = acr.grid.instrument.estimate_projection_on_ground(ins)
    assert xlim == (-50, 50)
    assert ylim == (-50, 50)


def test_projection_no_cone_45deg_flat():
    ins = {
        "bounding_box": {"width_x": 100, "width_y": 100, "height_z": 0},
        "viewcone": {"azimuth_deg": 0, "zenith_deg": 45, "half_angle_deg": 0},
    }
    xlim, ylim = acr.grid.instrument.estimate_projection_on_ground(ins)
    assert xlim == (-50, 50)
    assert ylim == (-50, 50)


def test_projection_no_cone_45deg():
    S = 100
    ins = {
        "bounding_box": {"width_x": S, "width_y": S, "height_z": S},
        "viewcone": {"azimuth_deg": 0, "zenith_deg": 45, "half_angle_deg": 0},
    }
    xlim, ylim = acr.grid.instrument.estimate_projection_on_ground(ins)
    assert xlim == (-S / 2 - S, S / 2)
    assert ylim == (-S / 2, S / 2)


def test_projection_cone5deg_straight_dwon():
    S = 100
    ins = {
        "bounding_box": {"width_x": S, "width_y": S, "height_z": S},
        "viewcone": {"azimuth_deg": 0, "zenith_deg": 0, "half_angle_deg": 5},
    }
    xlim, ylim = acr.grid.instrument.estimate_projection_on_ground(ins)

    off = S * np.tan(np.deg2rad(5))

    np.testing.assert_almost_equal(xlim[0], -S / 2 - off, 3)
    np.testing.assert_almost_equal(xlim[1], S / 2 + off, 3)
    np.testing.assert_almost_equal(ylim[0], -S / 2 - off, 3)
    np.testing.assert_almost_equal(ylim[1], S / 2 + off, 3)
