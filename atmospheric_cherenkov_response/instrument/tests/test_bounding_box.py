import atmospheric_cherenkov_response as acr
import numpy as np


def assert_corner(bb, ix, iy, iz, c, margin=1e-9):
    _c = acr.instrument.bounding_box.get_corner(
        bounding_box=bb, ix=ix, iy=iy, iz=iz
    )
    assert np.linalg.norm(_c - c) < margin


def test_corners_shape():
    bb = acr.instrument.bounding_box.init(width_x=2, width_y=2, height_z=1)
    c = acr.instrument.bounding_box.get_corner(
        bounding_box=bb, ix=0, iy=0, iz=0
    )
    assert c.shape == (3,)


def test_corners():
    bb = acr.instrument.bounding_box.init(width_x=2, width_y=2, height_z=1)

    assert_corner(bb, 0, 0, 0, [-1, -1, 0])
    assert_corner(bb, 0, 1, 0, [-1, 1, 0])
    assert_corner(bb, 1, 0, 0, [1, -1, 0])
    assert_corner(bb, 1, 1, 0, [1, 1, 0])

    assert_corner(bb, 0, 0, 1, [-1, -1, 1])
    assert_corner(bb, 0, 1, 1, [-1, 1, 1])
    assert_corner(bb, 1, 0, 1, [1, -1, 1])
    assert_corner(bb, 1, 1, 1, [1, 1, 1])
