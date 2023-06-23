import atmospheric_cherenkov_response as acr
from corsika_primary.I import BUNCH
import numpy as np


d2r = np.deg2rad


def sind(a):
    return np.sin(d2r(a))


def tand(a):
    return np.tan(d2r(a))


def approx(a, b, eps=1e-6):
    return np.abs(a - b) < eps


def test_normalize_axis_1():
    v = np.array([[1.0, 1, 1], [2.0, 2, 2],])
    vn = acr.instrument.dummy.normalize_axis_1(v=v)
    assert vn.shape == v.shape

    for row in range(len(vn)):
        nn = np.linalg.norm(vn[row])
        assert nn == 1.0


def test_calculate_photon_x_y_intersections_on_screen():
    # on axis
    x, y = acr.instrument.dummy.calculate_photon_x_y_intersections_on_screen(
        focal_length=1,
        screen_distance=1,
        cx_Tpap=np.array([0]),
        cy_Tpap=np.array([0]),
        x_Tpap=np.array([0]),
        y_Tpap=np.array([0]),
    )
    assert (x, y) == (0.0, 0.0)

    # small angle
    for x_Tpap in [-1, 0, 1]:
        for y_Tpap in [-1, 0, 1]:
            (
                x,
                y,
            ) = acr.instrument.dummy.calculate_photon_x_y_intersections_on_screen(
                focal_length=1,
                screen_distance=1,
                cx_Tpap=np.array([d2r(1)]),
                cy_Tpap=np.array([0]),
                x_Tpap=np.array([x_Tpap]),
                y_Tpap=np.array([y_Tpap]),
            )
            assert approx(x, -tand(1))

    # set focus different from infinity
    for x_Tpap in [-1, 0, 1]:
        for y_Tpap in [-1, 0, 1]:
            (
                x,
                y,
            ) = acr.instrument.dummy.calculate_photon_x_y_intersections_on_screen(
                focal_length=1,
                screen_distance=1.01,
                cx_Tpap=np.array([0]),
                cy_Tpap=np.array([0]),
                x_Tpap=np.array([x_Tpap]),
                y_Tpap=np.array([y_Tpap]),
            )

            x_exp = -0.01 * x_Tpap
            assert approx(x, x_exp)


def test_calculate_relative_path_length_for_isochor_imagen():
    # on axis
    d = acr.instrument.dummy.calculate_relative_path_length_for_isochor_imagen(
        cx_Tpap=np.array([0]),
        cy_Tpap=np.array([0]),
        x_Tpap=np.array([0]),
        y_Tpap=np.array([0]),
    )
    assert d.shape == (1,)
    assert d[0] == 0.0

    # smalle angle off axis but impacts in center of mirror
    d = acr.instrument.dummy.calculate_relative_path_length_for_isochor_imagen(
        cx_Tpap=np.array([0]),
        cy_Tpap=np.array([d2r(3)]),
        x_Tpap=np.array([0]),
        y_Tpap=np.array([0]),
    )
    assert d.shape == (1,)
    assert d[0] == 0.0

    # away from center of mirror
    #              z
    #              |
    #              |                       / incoming trajectory
    #              |                      /
    #              |                     /
    #              |                  1 /
    #              |           hyp    |/
    # -------------0--__--------------x-----> y, principal-aperture-plane (pap)
    #      origin       --__         /
    #                       --__    / opp = d
    #                           --_x

    d = acr.instrument.dummy.calculate_relative_path_length_for_isochor_imagen(
        cx_Tpap=np.array([0]),
        cy_Tpap=np.array([d2r(9)]),
        x_Tpap=np.array([0]),
        y_Tpap=np.array([1]),
    )
    # hyp = 1
    # opp = d
    # opp/hyp = sin(a)
    assert d.shape == (1,)
    assert approx(d[0], sind(9), eps=1e-3)


def test_calculate_photon_x_y_t_on_screen():
    c = 3e8
    f = 106.5
    pseudo_infinity = 1e99

    x, y, t = acr.instrument.dummy.calculate_photon_x_y_t_on_screen(
        focus_depth=pseudo_infinity,
        focal_length=f,
        cx_Tpap=np.array([0.0, d2r(1.0), 0.0]),
        cy_Tpap=np.array([d2r(3), 0.0, d2r(3)]),
        x_Tpap=np.array([0.0, 0.0, 0.0]),
        y_Tpap=np.array([30.0, 0.0, -30]),
        t_Tpap=np.array([0.0, 0.0, 0.0]),
    )
    assert x.shape == (3,)
    assert y.shape == (3,)
    assert t.shape == (3,)

    assert approx(x[0], 0.0)
    assert approx(y[0], -f * tand(3.0), eps=1e-5)
    assert approx(t[0], (30 * sind(3)) / c, eps=1e-10)

    assert approx(x[1], -f * tand(1.0), eps=1e-5)
    assert approx(y[1], 0.0)
    assert approx(
        t[1], 0.0, eps=1e-10
    )  # because impact is in center of mirror

    # photon 3: same direction as photon 1 but opposite site of mirror
    assert approx(x[2], x[0])
    assert approx(y[2], y[0])
    assert approx(t[2], -t[0])


def test_make_pixel_x_y():
    field_of_view_radius = 1.0
    pixel_radius = 0.05

    x, y = acr.instrument.dummy.make_pixel_x_y(
        field_of_view_radius=field_of_view_radius, pixel_radius=pixel_radius
    )

    A_fov = np.pi * field_of_view_radius ** 2
    A_pix = np.pi * pixel_radius ** 2

    num = int(np.round(A_fov / A_pix))

    assert 0.9 < len(x) / num < 1.1

    assert x.shape == y.shape
    assert np.abs(np.mean(x)) < 1e-2
    assert np.abs(np.mean(y)) < 1e-2

    assert 0.5 < np.std(x) < 1.0
    assert 0.5 < np.std(y) < 1.0


def test_get_cherenkov_bunches_which_cause_response_in_dummy_instrument():
    prng = np.random.Generator(np.random.PCG64(13))

    scenarios = [
        {"refl": 1.0, "pde": 0.5, "atmo": 1.0},
        {"refl": 0.5, "pde": 0.5, "atmo": 1.0},
        {"refl": 0.5, "pde": 0.5, "atmo": 0.5},
        {"refl": 1.0, "pde": 1.0, "atmo": 1.0},
    ]

    for scn in scenarios:
        toy = {}
        toy["mirror"] = {}
        toy["mirror"]["diameter_m"] = 2.0
        toy["mirror"]["reflectivity"] = {}
        toy["mirror"]["reflectivity"]["wavelength_m"] = [200e-9, 1200e-9]
        toy["mirror"]["reflectivity"]["value"] = [scn["refl"], scn["refl"]]
        toy["camera"] = {}
        toy["camera"]["efficiency"] = {}
        toy["camera"]["efficiency"]["wavelength_m"] = [200e-9, 1200e-9]
        toy["camera"]["efficiency"]["value"] = [scn["pde"], scn["pde"]]

        cherenkov_bunches_Tpap = np.array(
            # x/cm, y/cm, cx/1, cy/1, time/ns, zem/cm, bsize/1, wvl/nm
            [[0.0, 0.0, 0.0, 0.0, 0.0, 1e5, scn["atmo"], 433],],
            dtype=np.float32,
        )

        expected_ratio = scn["refl"] * scn["pde"] * scn["atmo"]

        NN = 1000
        nn = 0
        for i in range(NN):
            cer = acr.instrument.dummy.get_cherenkov_bunches_which_cause_response(
                cherenkov_bunches_Tpap=cherenkov_bunches_Tpap,
                dummy_instrument=toy,
                prng=prng,
            )
            nn += len(cer)

        assert (0.8 * expected_ratio) < nn / NN <= (1.2 * expected_ratio)


def test_mirror():
    prng = np.random.Generator(np.random.PCG64(13))

    toy = {}
    toy["mirror"] = {}
    toy["mirror"]["diameter_m"] = 2.0
    toy["mirror"]["reflectivity"] = {}
    toy["mirror"]["reflectivity"]["wavelength_m"] = [200e-9, 1200e-9]
    toy["mirror"]["reflectivity"]["value"] = [1, 1]
    toy["camera"] = {}
    toy["camera"]["efficiency"] = {}
    toy["camera"]["efficiency"]["wavelength_m"] = [200e-9, 1200e-9]
    toy["camera"]["efficiency"]["value"] = [1, 1]

    Nr = 100
    cherenkov_bunches_Tpap = []
    cm_over_m = 1e2
    for xx in np.linspace(-1, 1, 2 * Nr):
        for yy in np.linspace(-1, 1, 2 * Nr):
            ccc = [
                xx * cm_over_m,
                yy * cm_over_m,
                0.0,
                0.0,
                0.0,
                1e5,
                1.0,
                433,
            ]
            cherenkov_bunches_Tpap.append(ccc)
    cherenkov_bunches_Tpap = np.array(cherenkov_bunches_Tpap)

    assert len(cherenkov_bunches_Tpap) == 4 * Nr * Nr

    cer = acr.instrument.dummy.get_cherenkov_bunches_which_cause_response(
        cherenkov_bunches_Tpap=cherenkov_bunches_Tpap,
        dummy_instrument=toy,
        prng=prng,
    )

    actual_density = len(cer)
    expected_density = np.pi * (Nr * Nr)

    assert 0.95 < actual_density / expected_density < 1.05


def test_bin_t_into_time_slice():
    tx = acr.instrument.dummy.bin_t_into_time_slice(
        t_Tscreen=[-0.5, 0.5, 1.9, 3.4], time_bin_edges=[0, 1, 2]
    )

    assert tx[0] == -1
    assert tx[1] == 0
    assert tx[2] == 1
    assert tx[3] == -1


def test_bin_x_y_into_pixel():
    portal = acr.instrument.dummy.init_portal()

    px = acr.instrument.dummy.bin_x_y_into_pixel(
        x_Tscreen=[0, 100],
        y_Tscreen=[0, 100],
        pixel_x_y_tree=portal["camera"]["pixel"]["x_y_tree"],
        pixel_r=portal["camera"]["pixel"]["r_m"],
    )

    assert px[0] != -1
    assert px[1] == -1


def test_draw_night_sky_background():
    prng = np.random.Generator(np.random.PCG64(42))

    num_time_slices = 100
    num_pixel = 200

    dt = 1e-9
    R = 60e9

    nsb = acr.instrument.dummy.draw_night_sky_background(
        prng=prng,
        num_time_slices=num_time_slices,
        num_pixel=num_pixel,
        time_slice_duration_s=dt,
        rate_per_s=R,
    )

    I = R * dt

    assert nsb.shape == (num_time_slices, num_pixel)
    assert 0.95 < np.mean(nsb) / I < 1.05
    assert 0.95 < np.std(nsb) / np.sqrt(I) < 1.05


def test_make_cherenkov_bunches():
    prng = np.random.Generator(np.random.PCG64(42))

    cer = acr.instrument.dummy.make_cherenkov_bunches(
        prng=prng,
        size=1337,
        emission_point_m=[0, 0, 1e4],
        wavelength_m=433e-9,
        bunch_size=1.0,
        mirror_radius_m=10,
    )

    assert len(cer) == 1337
    assert approx(np.mean(cer[:, BUNCH.X]), 0.0, eps=10)  # cm
    assert approx(np.mean(cer[:, BUNCH.Y]), 0.0, eps=10)  # cm

    assert d2r(0.02) < np.std(cer[:, BUNCH.CX]) < d2r(0.1)  # rad
    assert d2r(0.02) < np.std(cer[:, BUNCH.CY]) < d2r(0.1)  # rad

    assert np.std(cer[:, BUNCH.TIME]) < 1  # ns

    assert np.all(cer[:, BUNCH.ZEM] == 1e4 * 1e2)  # cm
    assert np.all(cer[:, BUNCH.BSIZE] == 1.0)
    assert np.all(cer[:, BUNCH.WVL] == 433)


def test_response():
    prng = np.random.Generator(np.random.PCG64(42))

    emission_depth_m = 10e3

    num_cer = 42 * 1337
    cherenkov_bunches_Tpap = acr.instrument.dummy.make_cherenkov_bunches(
        prng=prng,
        size=num_cer,
        emission_point_m=[30, 400, emission_depth_m],
        wavelength_m=433e-9,
        bunch_size=1.0,
        mirror_radius_m=40,
    )

    assert len(cherenkov_bunches_Tpap) == num_cer

    portal_in_focus = acr.instrument.dummy.init_portal(
        focus_depth_m=emission_depth_m
    )

    resi, truth = acr.instrument.dummy.estimate_response_to_cherenkov(
        prng=prng,
        dummy_instrument=portal_in_focus,
        cherenkov_bunches_Tpap=cherenkov_bunches_Tpap,
    )

    assert np.sum(resi) > 0.2 * num_cer
    assert np.median(resi) == 0
    assert np.max(resi) < num_cer

    time_resi = np.sum(resi, axis=1)
    assert time_resi.shape == (
        portal_in_focus["camera"]["time_slices"]["num"],
    )
    assert np.sum(time_resi > 0) < 3

    img_resi = np.sum(resi, axis=0)
    assert 0 < np.sum(img_resi > 0) < 5

    # acr.instrument.plot.plot_response(
    #    path="portal_in_focus.jpg",
    #    toy_instrument=portal_in_focus,
    #    response=resi,
    # )

    assert len(cherenkov_bunches_Tpap) == num_cer

    # out of focus
    # ------------
    portal_out_of_focus = acr.instrument.dummy.init_portal(focus_depth_m=3e3)

    reso, truth = acr.instrument.dummy.estimate_response_to_cherenkov(
        prng=prng,
        dummy_instrument=portal_out_of_focus,
        cherenkov_bunches_Tpap=cherenkov_bunches_Tpap,
    )

    time_reso = np.sum(reso, axis=1)
    assert time_reso.shape == (
        portal_out_of_focus["camera"]["time_slices"]["num"],
    )
    assert np.sum(time_reso > 0) < 3

    img_reso = np.sum(reso, axis=0)
    assert 30 < np.sum(img_reso > 0) < 60

    # acr.instrument.plot.plot_response(
    #    path="portal_out_of_focus.jpg",
    #    toy_instrument=portal_out_of_focus,
    #    response=reso,
    # )
