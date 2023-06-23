import numpy as np
import photon_spectra
import binning_utils
import solid_angle_utils
import plenopy
from optic_object_wavefronts.geometry.grid import hexagonal as hexgrid
from corsika_primary.I import BUNCH
from .. import night_sky_background
import scipy
from scipy import spatial


def assert_wavelength_values(wavelength, values):
    assert len(wavelength) == len(values)
    assert np.all(wavelength > 0)
    assert np.all(np.gradient(wavelength) > 0)
    assert np.all(values >= 0)


def init(
    mirror_diameter_m,
    mirror_focal_length_m,
    field_of_view_half_angle_deg,
    focus_depth_m,
    pixel_half_angle_deg,
    trigger_threshold_num_photons,
    trigger_integration_num_time_slices,
    time_slice_duration_s,
    num_time_slices,
    reflectivity_wavelength_m,
    reflectivity,
    photon_detection_efficiency_wavelength_m,
    photon_detection_efficiency,
    background_wavelength_m,
    background_differential_flux_per_m2_per_sr_per_s_per_m,
):
    assert mirror_diameter_m > 0
    assert mirror_focal_length_m > 0
    assert field_of_view_half_angle_deg > 0
    assert focus_depth_m > 0
    assert pixel_half_angle_deg > 0
    assert trigger_threshold_num_photons > 0
    assert trigger_integration_num_time_slices > 0
    assert time_slice_duration_s > 0
    assert num_time_slices > 0

    assert_wavelength_values(
        wavelength=reflectivity_wavelength_m, values=reflectivity
    )
    assert_wavelength_values(
        wavelength=photon_detection_efficiency_wavelength_m,
        values=photon_detection_efficiency,
    )

    f = mirror_focal_length_m
    pixel_x_Tscreen, pixel_y_Tscreen = make_pixel_x_y(
        field_of_view_radius=f * tandeg(field_of_view_half_angle_deg),
        pixel_radius=f * tandeg(pixel_half_angle_deg),
    )

    pixel = {
        "x_m": pixel_x_Tscreen,
        "y_m": pixel_y_Tscreen,
        "num": len(pixel_y_Tscreen),
        "r_m": f * tandeg(pixel_half_angle_deg),
        "r_deg": pixel_half_angle_deg,
    }
    pixel["x_y_tree"] = scipy.spatial.cKDTree(
        data=np.c_[pixel["x_m"], pixel["y_m"]]
    )

    dt = time_slice_duration_s
    num_t = num_time_slices
    time = {
        "num": num_t,
        "duration_s": dt,
        "edges_s": np.linspace(0.0, dt * num_t, num_t + 1) - 0.5 * dt * num_t,
    }

    dummy_instrument = {
        "mirror": {
            "diameter_m": mirror_diameter_m,
            "focal_length_m": mirror_focal_length_m,
            "reflectivity": {
                "wavelength_m": reflectivity_wavelength_m,
                "value": reflectivity,
            },
        },
        "camera": {
            "field_of_view_half_angle_deg": field_of_view_half_angle_deg,
            "pixel": pixel,
            "time_slices": time,
            "efficiency": {
                "wavelength_m": photon_detection_efficiency_wavelength_m,
                "value": photon_detection_efficiency,
            },
            "focus_depth_m": focus_depth_m,
        },
        "trigger": {
            "integration_num_time_slices": trigger_integration_num_time_slices,
        },
    }

    dummy_instrument["camera"]["pixel"][
        "average_background_rate_per_s"
    ] = estimate_rate_in_trigger_pixel_from_night_sky_background(
        dummy_instrument=dummy_instrument,
        background_wavelength_m=background_wavelength_m,
        background_differential_flux_per_m2_per_sr_per_s_per_m=background_differential_flux_per_m2_per_sr_per_s_per_m,
        num_wavelength_steps=1337,
    )

    return dummy_instrument


def init_portal(focus_depth_m=1e4):
    opt = photon_spectra.cta_mirrors.init("cta_mst_dielectric_after")
    pde = photon_spectra.hamamatsu_r11920_100_05.init()

    nsb = (
        night_sky_background.init_night_sky_background_differential_flux_per_m2_per_sr_per_s_per_m()
    )
    background_wavelength_m = nsb["wavelength_m"]
    background_differential_flux_per_m2_per_sr_per_s_per_m = nsb["value"]

    return init(
        mirror_diameter_m=71,
        mirror_focal_length_m=71 * 1.5,
        field_of_view_half_angle_deg=3.25,
        focus_depth_m=focus_depth_m,
        pixel_half_angle_deg=2.2 * 0.5 * 0.067,
        trigger_threshold_num_photons=50,
        trigger_integration_num_time_slices=10,
        time_slice_duration_s=0.5e-9,
        num_time_slices=200,
        reflectivity_wavelength_m=opt["wavelength"],
        reflectivity=opt["value"],
        photon_detection_efficiency_wavelength_m=pde["wavelength"],
        photon_detection_efficiency=pde["value"],
        background_wavelength_m=background_wavelength_m,
        background_differential_flux_per_m2_per_sr_per_s_per_m=background_differential_flux_per_m2_per_sr_per_s_per_m,
    )


def estimate_response(prng, dummy_instrument, cherenkov_bunches_Tpap):
    cer, cer_truth = estimate_response_to_cherenkov(
        prng=prng,
        dummy_instrument=dummy_instrument,
        cherenkov_bunches_Tpap=cherenkov_bunches_Tpap,
    )
    nsb = estimate_response_to_night_sky_background(
        prng=prng, dummy_instrument=dummy_instrument,
    )
    res = nsb + cer
    return res, {"cherenkov": cer_truth}


def estimate_response_to_night_sky_background(prng, dummy_instrument):
    dum = dummy_instrument
    return draw_night_sky_background(
        prng=prng,
        num_time_slices=dum["camera"]["time_slices"]["num"],
        num_pixel=dum["camera"]["pixel"]["num"],
        time_slice_duration_s=dum["camera"]["time_slices"]["duration_s"],
        rate_per_s=dum["camera"]["pixel"]["average_background_rate_per_s"],
    )


def draw_night_sky_background(
    prng, num_time_slices, num_pixel, time_slice_duration_s, rate_per_s
):
    mean = rate_per_s * time_slice_duration_s
    assert mean > 10
    std = np.sqrt(mean)
    size = num_time_slices * num_pixel
    nnn = np.round(prng.normal(loc=mean, scale=std, size=size))
    nnn = nnn.astype(np.int)
    nnn[nnn < 0] = 0
    nsb = nnn.reshape((num_time_slices, num_pixel))
    return nsb


def estimate_response_to_cherenkov(
    prng,
    dummy_instrument,
    cherenkov_bunches_Tpap,
    speed_of_light_m_per_s=299792458,
):
    dum = dummy_instrument

    cer_Tpap = get_cherenkov_bunches_which_cause_response(
        cherenkov_bunches_Tpap=cherenkov_bunches_Tpap,
        dummy_instrument=dum,
        prng=prng,
    )

    m_over_cm = 1e-2
    s_over_ns = 1e9
    x_Tscreen, y_Tscreen, t_Tscreen = calculate_photon_x_y_t_on_screen(
        focus_depth=dum["camera"]["focus_depth_m"],
        focal_length=dum["mirror"]["focal_length_m"],
        cx_Tpap=cer_Tpap[:, BUNCH.CX],
        cy_Tpap=cer_Tpap[:, BUNCH.CY],
        x_Tpap=cer_Tpap[:, BUNCH.X] * m_over_cm,
        y_Tpap=cer_Tpap[:, BUNCH.Y] * m_over_cm,
        t_Tpap=cer_Tpap[:, BUNCH.TIME] * s_over_ns,
        speed_of_light=speed_of_light_m_per_s,
    )

    true_arrival_time_Tscreen = np.median(t_Tscreen)
    t_rel = t_Tscreen - true_arrival_time_Tscreen

    pixel_id, tslice_id = bin_x_y_t_into_pixel_and_time_slice(
        x_Tscreen=x_Tscreen,
        y_Tscreen=y_Tscreen,
        t_Tscreen=t_rel,
        pixel_x_y_tree=dum["camera"]["pixel"]["x_y_tree"],
        pixel_r=dum["camera"]["pixel"]["r_m"],
        time_bin_edges=dum["camera"]["time_slices"]["edges_s"],
    )
    image_sequence = np.zeros(
        shape=(
            dum["camera"]["time_slices"]["num"],
            dum["camera"]["pixel"]["num"],
        )
    )
    mask_in_pixel = pixel_id != -1
    mask_in_time_slice = np.logical_and(
        tslice_id < dum["camera"]["time_slices"]["num"], tslice_id > -1
    )
    mask = np.logical_and(mask_in_pixel, mask_in_time_slice)

    pixel_id = pixel_id[mask]
    tslice_id = tslice_id[mask]

    for iii in range(len(pixel_id)):
        image_sequence[tslice_id[iii], pixel_id[iii]] += 1

    truth = {"arrival_time_s": true_arrival_time_Tscreen}
    return image_sequence, truth


def tandeg(a):
    return np.tan(np.deg2rad(a))


def estimate_rate_in_trigger_pixel_from_night_sky_background(
    dummy_instrument,
    background_wavelength_m,
    background_differential_flux_per_m2_per_sr_per_s_per_m,
    num_wavelength_steps=1337,
):
    dum = dummy_instrument

    common_wavelength_bin = binning_utils.Binning(
        bin_edges=np.linspace(240e-9, 700e-9, num_wavelength_steps + 1)
    )

    q_opt = photon_spectra.utils.make_values_for_common_wavelength(
        wavelength=dum["mirror"]["reflectivity"]["wavelength_m"],
        value=dum["mirror"]["reflectivity"]["value"],
        common_wavelength=common_wavelength_bin["centers"],
    )

    q_pde = photon_spectra.utils.make_values_for_common_wavelength(
        wavelength=dum["camera"]["efficiency"]["wavelength_m"],
        value=dum["camera"]["efficiency"]["value"],
        common_wavelength=common_wavelength_bin["centers"],
    )

    q_nsb_per_m2_per_sr_per_s_per_m = photon_spectra.utils.make_values_for_common_wavelength(
        wavelength=background_wavelength_m,
        value=background_differential_flux_per_m2_per_sr_per_s_per_m,
        common_wavelength=common_wavelength_bin["centers"],
    )

    # integrate
    detected_nsb_flux_per_m2_per_sr_per_s = 0.0
    for i, delta_wavelength_m in enumerate(common_wavelength_bin["widths"]):
        detected_nsb_flux_per_m2_per_sr_per_s += (
            q_opt[i]
            * q_pde[i]
            * q_nsb_per_m2_per_sr_per_s_per_m[i]
            * delta_wavelength_m
        )

    A_m2 = np.pi * (dum["mirror"]["diameter_m"] / 2) ** 2
    O_sr = solid_angle_utils.cone.solid_angle(
        half_angle_rad=np.deg2rad(dum["camera"]["pixel"]["r_deg"])
    )

    detected_nsb_rate_per_s = (
        detected_nsb_flux_per_m2_per_sr_per_s * A_m2 * O_sr
    )

    return detected_nsb_rate_per_s


def evaluate_by_wavelength(
    photon_wavelength_m, function_wavelength_m, function_value
):
    return np.interp(
        x=photon_wavelength_m, xp=function_wavelength_m, fp=function_value
    )


def get_cherenkov_bunches_which_cause_response(
    cherenkov_bunches_Tpap, dummy_instrument, prng,
):
    dum = dummy_instrument
    cer = cherenkov_bunches_Tpap

    # remove bunches outside of mirror
    # --------------------------------
    m_over_cm = 1e-2
    cer_x_m = m_over_cm * cer[:, BUNCH.X]
    cer_y_m = m_over_cm * cer[:, BUNCH.Y]

    bunch_r2_m2 = cer_x_m ** 2 + cer_y_m ** 2
    mirror_r2_m2 = (0.5 * dum["mirror"]["diameter_m"]) ** 2
    mask_on_mirror = bunch_r2_m2 <= mirror_r2_m2
    cer = cer[mask_on_mirror]

    # remove bunches (atmosphere, mirror, photo-sensor)
    # -------------------------------------------------
    p_passing_atmosphere = cer[:, BUNCH.BSIZE]
    assert np.all(p_passing_atmosphere <= 1.0)
    assert np.all(p_passing_atmosphere >= 0.0)

    p_reflected_by_mirror = evaluate_by_wavelength(
        photon_wavelength_m=1e-9 * cer[:, BUNCH.WVL],
        function_wavelength_m=dum["mirror"]["reflectivity"]["wavelength_m"],
        function_value=dum["mirror"]["reflectivity"]["value"],
    )
    assert np.all(p_reflected_by_mirror <= 1.0)
    assert np.all(p_reflected_by_mirror >= 0.0)

    p_detected_by_photo_sensor = evaluate_by_wavelength(
        photon_wavelength_m=1e-9 * cer[:, BUNCH.WVL],
        function_wavelength_m=dum["camera"]["efficiency"]["wavelength_m"],
        function_value=dum["camera"]["efficiency"]["value"],
    )
    assert np.all(p_detected_by_photo_sensor <= 1.0)
    assert np.all(p_detected_by_photo_sensor >= 0.0)

    p_total = (
        p_passing_atmosphere
        * p_reflected_by_mirror
        * p_detected_by_photo_sensor
    )

    prob = prng.uniform(low=0, high=1, size=len(cer))

    mask = p_total >= prob
    return cer[mask]


def make_pixel_x_y(
    field_of_view_radius, pixel_radius,
):
    spacing = 2 * pixel_radius
    num_pixel_on_diagonal = int(
        np.ceil(2 * field_of_view_radius / pixel_radius)
    )
    grid = hexgrid.init_from_spacing(
        spacing=spacing, fN=num_pixel_on_diagonal,
    )
    out = []
    for key in grid:
        ppos = grid[key]
        if ppos[0] ** 2 + ppos[1] ** 2 <= field_of_view_radius ** 2:
            out.append(ppos)
    pos_xyz = np.array(out)
    return pos_xyz[:, 0], pos_xyz[:, 1]


def normalize_axis_1(v):
    no = np.linalg.norm(v, axis=1)
    v[:, 0] /= no
    v[:, 1] /= no
    v[:, 2] /= no
    return v


def calculate_photon_x_y_intersections_on_screen(
    focal_length, screen_distance, cx_Tpap, cy_Tpap, x_Tpap, y_Tpap
):
    """
    Based on the thin-lens.

    Parameters
    ----------
    focal_length : float
        The focal-length of the imaging optics.
    screen_distance : float
        The distance where the screen is placed.
        If screen_distance == focal_length, the screen is focused to infinity.
    cx_Tpap, cy_Tpap, x_Tpap, y_Tpap : arrays of floats
        The photons directions (cx, cy) and impacts (x, y) relative to
        the aperture's principal plane (PAP).
    """
    assert focal_length > 0.0
    assert screen_distance > 0.0
    num = len(cx_Tpap)
    assert num == len(cy_Tpap)
    assert num == len(x_Tpap)
    assert num == len(y_Tpap)

    f = focal_length
    d = screen_distance
    tan = np.tan

    intersections_on_focal_plane = np.array(
        [-f * tan(cx_Tpap), -f * tan(cy_Tpap), f * np.ones(num)]
    ).T

    supports_on_aperture = np.array([x_Tpap, y_Tpap, np.zeros(num)]).T

    directions_leaving_aperture = (
        intersections_on_focal_plane - supports_on_aperture
    )
    directions_leaving_aperture = normalize_axis_1(directions_leaving_aperture)
    scale_factors = screen_distance / directions_leaving_aperture[:, 2]

    intersections_on_screen = (
        supports_on_aperture
        + (scale_factors * directions_leaving_aperture.T).T
    )

    return intersections_on_screen[:, 0], intersections_on_screen[:, 1]


def calculate_relative_path_length_for_isochor_imagen(
    cx_Tpap, cy_Tpap, x_Tpap, y_Tpap
):
    """
    The ray r(d) = (x,y,0)^T + d*(cx,cy,cz)^T has a point closest to the
    aperture's principal plane's origin (0,0,0)^T.
    The 'd' to reach this point on the ray r(d) is the path-length we are
    looking for.
    """
    d = cx_Tpap * x_Tpap + cy_Tpap * y_Tpap
    return d


def calculate_photon_x_y_t_on_screen(
    focus_depth,
    focal_length,
    cx_Tpap,
    cy_Tpap,
    x_Tpap,
    y_Tpap,
    t_Tpap,
    speed_of_light=299792458,
):
    screen_distance = plenopy.thin_lens.object_distance_2_image_distance(
        object_distance=focus_depth, focal_length=focal_length,
    )
    x_Tscreen, y_Tscreen = calculate_photon_x_y_intersections_on_screen(
        focal_length=focal_length,
        screen_distance=screen_distance,
        cx_Tpap=cx_Tpap,
        cy_Tpap=cy_Tpap,
        x_Tpap=x_Tpap,
        y_Tpap=y_Tpap,
    )
    dd = calculate_relative_path_length_for_isochor_imagen(
        cx_Tpap=cx_Tpap, cy_Tpap=cy_Tpap, x_Tpap=x_Tpap, y_Tpap=y_Tpap,
    )

    dt = dd / speed_of_light
    t_Tscreen = t_Tpap + dt

    return x_Tscreen, y_Tscreen, t_Tscreen


def bin_x_y_t_into_pixel_and_time_slice(
    x_Tscreen, y_Tscreen, t_Tscreen, pixel_x_y_tree, pixel_r, time_bin_edges,
):
    px = bin_x_y_into_pixel(
        x_Tscreen=x_Tscreen,
        y_Tscreen=y_Tscreen,
        pixel_x_y_tree=pixel_x_y_tree,
        pixel_r=pixel_r,
    )

    tx = bin_t_into_time_slice(
        t_Tscreen=t_Tscreen, time_bin_edges=time_bin_edges
    )

    return px, tx


def bin_x_y_into_pixel(x_Tscreen, y_Tscreen, pixel_x_y_tree, pixel_r):
    """
    returns -1 when missed
    """
    dist, pixel_i = pixel_x_y_tree.query(np.c_[x_Tscreen, y_Tscreen])
    pixel_i[dist >= pixel_r] = -1
    return pixel_i


def bin_t_into_time_slice(t_Tscreen, time_bin_edges):
    """
    returns -1 when missed
    """
    tslice_i = np.digitize(x=t_Tscreen, bins=time_bin_edges, right=False) - 1
    num_bins = len(time_bin_edges) - 1

    mask_miss = tslice_i >= num_bins
    tslice_i[mask_miss] = -1

    return tslice_i


def make_cherenkov_bunches(
    prng,
    size=1000,
    emission_point_m=[0, 0, 1e4],
    wavelength_m=433e-9,
    bunch_size=1.0,
    mirror_radius_m=10,
):
    cer = []

    cm_over_m = 1e2
    ns_over_s = 1e-9
    nm_over_m = 1e9
    speed_of_light = 299792458

    mirror_impacts_x = prng.uniform(
        low=-mirror_radius_m, high=mirror_radius_m, size=size
    )
    mirror_impacts_y = prng.uniform(
        low=-mirror_radius_m, high=mirror_radius_m, size=size
    )

    mirror_impacts = np.array(
        [mirror_impacts_x, mirror_impacts_y, np.zeros(size)]
    ).T

    incidents = -mirror_impacts + emission_point_m

    distances = np.linalg.norm(incidents, axis=1)

    cxcycz = np.array(incidents)
    cxcycz[:, 0] /= distances
    cxcycz[:, 1] /= distances
    cxcycz[:, 2] /= distances

    cer = np.zeros(shape=(size, 8), dtype=np.float32)
    times = distances / speed_of_light

    # x/cm, y/cm, cx/1, cy/1, time/ns, zem/cm, bsize/1, wvl/nm
    # --------------------------------------------------------
    cer[:, BUNCH.X] = mirror_impacts_x * cm_over_m
    cer[:, BUNCH.Y] = mirror_impacts_y * cm_over_m
    cer[:, BUNCH.CX] = cxcycz[:, 0]
    cer[:, BUNCH.CY] = cxcycz[:, 1]
    cer[:, BUNCH.TIME] = times * ns_over_s
    cer[:, BUNCH.ZEM] = np.ones(size) * emission_point_m[2] * cm_over_m
    cer[:, BUNCH.BSIZE] = np.ones(size) * bunch_size
    cer[:, BUNCH.WVL] = np.ones(size) * wavelength_m * nm_over_m
    return cer
