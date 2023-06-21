import numpy as np
import photon_spectra
import binning_utils
import solid_angle_utils
from optic_object_wavefronts.geometry.grid import hexagonal as hexgrid
from corsika_primary.I import BUNCH


def assert_wavelength_values(wavelength, values):
    assert len(wavelength) == len(values)
    assert np.all(wavelength > 0)
    assert np.all(np.gradient(wavelength) > 0)
    assert np.all(values >= 0)


def init(
    mirror_diameter_m,
    field_of_view_half_angle_deg,
    focus_depth_m,
    trigger_pixel_half_angle_deg,
    trigger_threshold_num_photons,
    trigger_integration_time_s,
    trigger_time_sample_duration_s,
    optical_efficiency_wavelength_m,
    optical_efficiency,
    photon_detection_efficiency_wavelength_m,
    photon_detection_efficiency,
):
    assert mirror_diameter_m > 0
    assert field_of_view_half_angle_deg > 0
    assert focus_depth_m > 0
    assert trigger_pixel_half_angle_deg > 0
    assert trigger_threshold_num_photons > 0
    assert trigger_integration_time_s > 0
    assert trigger_time_sample_duration_s > 0
    assert trigger_time_sample_duration_s <= trigger_integration_time_s

    assert_wavelength_values(
        wavelength=optical_efficiency_wavelength_m, values=optical_efficiency
    )
    assert_wavelength_values(
        wavelength=photon_detection_efficiency_wavelength_m,
        values=photon_detection_efficiency,
    )

    return {
        "mirror_diameter_m": mirror_diameter_m,
        "field_of_view_half_angle_deg": field_of_view_half_angle_deg,
        "focus_depth_m": focus_depth_m,
        "trigger_pixel_half_angle_deg": trigger_pixel_half_angle_deg,
        "trigger_threshold_num_photons": trigger_threshold_num_photons,
        "trigger_integration_time_s": trigger_integration_time_s,
        "trigger_time_sample_duration_s": trigger_time_sample_duration_s,
        "optical_efficiency": {
            "wavelength_m": optical_efficiency_wavelength_m,
            "value": optical_efficiency,
        },
        "photon_detection_efficiency": {
            "wavelength_m": photon_detection_efficiency_wavelength_m,
            "value": photon_detection_efficiency,
        },
    }


def init_portal():
    opt = photon_spectra.cta_mirrors.init("cta_mst_dielectric_after")
    pde = photon_spectra.hamamatsu_r11920_100_05.init()

    return init(
        mirror_diameter_m=71,
        field_of_view_half_angle_deg=3.25,
        focus_depth_m=10e3,
        trigger_pixel_half_angle_deg=0.033,
        trigger_threshold_num_photons=50,
        trigger_integration_time_s=5e-9,
        trigger_time_sample_duration_s=0.5e-9,
        optical_efficiency_wavelength_m=opt["wavelength"],
        optical_efficiency=opt["value"],
        photon_detection_efficiency_wavelength_m=pde["wavelength"],
        photon_detection_efficiency=pde["value"],
    )


def estimate_rate_in_trigger_pixel_from_night_sky_background(
    dummy_instrument,
    night_sky_background_differential_flux_per_m2_per_sr_per_s_per_m,
    num_wavelength_steps=1337,
):
    di = dummy_instrument
    nsb_per_m2_per_sr_per_s_per_m = (
        night_sky_background_differential_flux_per_m2_per_sr_per_s_per_m
    )

    common_wavelength_bin = binning_utils.Binning(
        bin_edges=np.linspace(240e-9, 700e-9, num_wavelength_steps + 1)
    )

    q_opt = photon_spectra.utils.make_values_for_common_wavelength(
        wavelength=di["optical_efficiency"]["wavelength_m"],
        value=di["optical_efficiency"]["value"],
        common_wavelength=common_wavelength_bin["centers"],
    )

    q_pde = photon_spectra.utils.make_values_for_common_wavelength(
        wavelength=di["photon_detection_efficiency"]["wavelength_m"],
        value=di["photon_detection_efficiency"]["value"],
        common_wavelength=common_wavelength_bin["centers"],
    )

    q_nsb_per_m2_per_sr_per_s_per_m = photon_spectra.utils.make_values_for_common_wavelength(
        wavelength=nsb_per_m2_per_sr_per_s_per_m["wavelength_m"],
        value=nsb_per_m2_per_sr_per_s_per_m["value"],
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

    A_m2 = np.pi * (di["mirror_diameter_m"] / 2) ** 2
    O_sr = solid_angle_utils.cone.solid_angle(
        half_angle_rad=np.deg2rad(di["trigger_pixel_half_angle_deg"])
    )

    detected_nsb_rate_per_s = (
        detected_nsb_flux_per_m2_per_sr_per_s * A_m2 * O_sr
    )

    return detected_nsb_rate_per_s


def estimate_trigger_response(
    dummy_instrument,
    night_sky_background_differential_flux_per_m2_per_sr_per_s_per_m,
    cherenkov_bunches,
    prng,
):
    # remove bunches outside of mirror
    # --------------------------------
    bunch_r2 = (
        cherenkov_bunches[:, BUNCH.X] ** 2 + cherenkov_bunches[:, BUNCH.Y] ** 2
    )
    mirror_r2 = (0.5 * dummy_instrument["mirror_diameter_m"]) ** 2
    mask_on_mirror = r2 <= mirror_r2

    cherenkov_bunches = cherenkov_bunches[mask_on_mirror]

    # remove bunches absorbed in atmosphere
    # -------------------------------------
    prob = prng.uniform(low=0, high=1, size=len(cherenkov_bunches))
    mask_passing_atmosphere = cherenkov_bunches[:, BUNCH.BSIZE] >= prob
    cherenkov_bunches = cherenkov_bunches[prob]

    # assign photon to pixel
    # ----------------------


def make_trigger_pixel_centers(
    field_of_view_half_angle_deg, trigger_pixel_half_angle_deg,
):
    spacing_deg = trigger_pixel_half_angle_deg  # overlap is intended
    num_pixel_on_diagonal = int(
        np.ceil(
            2 * field_of_view_half_angle_deg / trigger_pixel_half_angle_deg
        )
    )
    grid = hexgrid.init_from_spacing(
        spacing=spacing_deg, fN=num_pixel_on_diagonal,
    )
    out = []
    for key in grid:
        ppos = grid[key]
        if ppos[0] ** 2 + ppos[1] ** 2 <= field_of_view_half_angle_deg ** 2:
            out.append(ppos)
    pos_xyz = np.array(out)
    return pos_xyz[:, 0:2]


def normalize_axis_1(v):
    no = np.linalg.norm(v, axis=1)
    v[:, 0] /= no
    v[:, 1] /= no
    v[:, 2] /= no
    return v


def compute_photon_x_y_intersections_on_screen(
    focal_length, screen_distance, cx, cy, x, y
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
    cx, cy, x, y : arrays of floats
        The photons directions (cx, cy) and impacts (x, y) relative to
        the aperture's principal plane.
    """
    assert focal_length > 0.0
    assert screen_distance > 0.0
    num = len(cx)
    assert num == len(cy)
    assert num == len(x)
    assert num == len(y)

    f = focal_length
    d = screen_distance
    tan = np.tan

    intersections_on_focal_plane = np.array(
        [-f * tan(cx), -f * tan(cy), f * np.ones(num)]
    )

    supports_on_aperture = np.array([x, y, np.zeros(num)])

    directions_leaving_aperture = (
        intersections_on_focal_plane - supports_on_aperture
    )
    directions_leaving_aperture = normalize_axis_1(directions_leaving_aperture)

    scale_factor = screen_distance / directions_leaving_aperture[:, 2]

    intersections_on_screen = (
        supports_on_aperture
        + (scale_factors * directions_leaving_aperture.T).T
    )

    return intersections_on_screen
