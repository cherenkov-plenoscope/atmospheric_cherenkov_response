import numpy as np
import photon_spectra
import binning_utils


def init_night_sky_background_differential_flux_per_m2_per_sr_per_s_per_m():
    nsb = photon_spectra.nsb_la_palma_2013_benn.init()
    return {
        "wavelength_m": nsb["wavelength"],
        "value": nsb["value"],
    }
