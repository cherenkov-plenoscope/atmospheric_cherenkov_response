import numpy as np
import photon_spectra
import binning_utils
import solid_angle_utils


def assert_wavelength_values(wavelength, values):
    assert len(wavelength) == len(values)
    assert np.all(wavelength > 0)
    assert np.all(np.gradient(wavelength) > 0)
    assert np.all(values >= 0)


def init(
    mirror_diameter_m,
    field_of_view_half_angle_deg,
    trigger_pixel_half_angle_deg,
    trigger_threshold_num_photons,
    trigger_integration_time_s,
    optical_efficiency_wavelength_m,
    optical_efficiency,
    photon_detection_efficiency_wavelength_m,
    photon_detection_efficiency,
    night_sky_background_wavelength_m,
    night_sky_background_differential_flux_per_m2_per_sr_per_s_per_m,
):
    assert mirror_diameter_m > 0
    assert field_of_view_half_angle_deg > 0
    assert trigger_pixel_half_angle_deg > 0
    assert trigger_threshold_num_photons > 0
    assert trigger_integration_time_s > 0

    assert_wavelength_values(
        wavelength=optical_efficiency_wavelength_m, values=optical_efficiency
    )
    assert_wavelength_values(
        wavelength=photon_detection_efficiency_wavelength_m,
        values=photon_detection_efficiency,
    )
    assert_wavelength_values(
        wavelength=night_sky_background_wavelength_m,
        values=night_sky_background_differential_flux_per_m2_per_sr_per_s_per_m,
    )

    return {
        "mirror_diameter_m": mirror_diameter_m,
        "field_of_view_half_angle_deg": field_of_view_half_angle_deg,
        "trigger_pixel_half_angle_deg": trigger_pixel_half_angle_deg,
        "trigger_threshold_num_photons": trigger_threshold_num_photons,
        "trigger_integration_time_s": trigger_integration_time_s,
        "optical_efficiency": {
            "wavelength_m": optical_efficiency_wavelength_m,
            "value": optical_efficiency,
        },
        "photon_detection_efficiency": {
            "wavelength_m": photon_detection_efficiency_wavelength_m,
            "value": photon_detection_efficiency,
        },
        "night_sky_background_differential_flux_per_m2_per_sr_per_s_per_m": {
            "wavelength_m": night_sky_background_wavelength_m,
            "value": night_sky_background_differential_flux_per_m2_per_sr_per_s_per_m,
        },
    }


def init_portal():
    opt = photon_spectra.cta_mirrors.init("cta_mst_dielectric_after")
    pde = photon_spectra.hamamatsu_r11920_100_05.init()
    nsb = photon_spectra.nsb_la_palma_2013_benn.init()

    return init(
        mirror_diameter_m=71,
        field_of_view_half_angle_deg=3.25,
        trigger_pixel_half_angle_deg=0.033,
        trigger_threshold_num_photons=50,
        trigger_integration_time_s=5e-9,
        optical_efficiency_wavelength_m=opt["wavelength"],
        optical_efficiency=opt["value"],
        photon_detection_efficiency_wavelength_m=pde["wavelength"],
        photon_detection_efficiency=pde["value"],
        night_sky_background_wavelength_m=nsb["wavelength"],
        night_sky_background_differential_flux_per_m2_per_sr_per_s_per_m=nsb[
            "value"
        ],
    )


def estimate_rate_in_trigger_pixel_from_night_sky_background(
    dummy_instrument, num_wavelength_steps=1337
):
    di = dummy_instrument

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

    q_nsb = photon_spectra.utils.make_values_for_common_wavelength(
        wavelength=di[
            "night_sky_background_differential_flux_per_m2_per_sr_per_s_per_m"
        ]["wavelength_m"],
        value=di[
            "night_sky_background_differential_flux_per_m2_per_sr_per_s_per_m"
        ]["value"],
        common_wavelength=common_wavelength_bin["centers"],
    )

    # integrate
    detected_nsb_flux_per_m2_per_sr_per_s = 0.0
    for i, delta_wavelength_m in enumerate(common_wavelength_bin["widths"]):
        detected_nsb_flux_per_m2_per_sr_per_s += (
            q_opt[i] * q_pde[i] * q_nsb[i] * delta_wavelength_m
        )

    A_m2 = np.pi * (di["mirror_diameter_m"] / 2) ** 2
    O_sr = solid_angle_utils.cone.solid_angle(
        half_angle_rad=np.deg2rad(di["trigger_pixel_half_angle_deg"])
    )

    detected_nsb_rate_per_s = (
        detected_nsb_flux_per_m2_per_sr_per_s * A_m2 * O_sr
    )

    return detected_nsb_rate_per_s


def estimate_trigger_response(instrument, cherenkov_bunches):
    pass
