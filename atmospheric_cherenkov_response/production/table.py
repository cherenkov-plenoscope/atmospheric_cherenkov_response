import numpy as np

STRUCTURE = {}

STRUCTURE["primary"] = {
    "particle_id": {"dtype": "<i8", "comment": "CORSIKA particle-id"},
    "energy_GeV": {"dtype": "<f8", "comment": ""},
    "azimuth_rad": {
        "dtype": "<f8",
        "comment": "Direction of the primary particle w.r.t. magnetic north.",
    },
    "zenith_rad": {
        "dtype": "<f8",
        "comment": "Direction of the primary particle.",
    },
    "max_scatter_rad": {"dtype": "<f8", "comment": ""},
    "magnet_azimuth_rad": {
        "dtype": "<f8",
        "comment": "The azimuth direction that the primary particle needs "
        "to have in order to induce an air-shower that emits its "
        "Cherenkov-light head on the pointing of the instrument.",
    },
    "magnet_zenith_rad": {
        "dtype": "<f8",
        "comment": "The zenith direction that the primary particle needs "
        "to have in order to induce an air-shower that emits its "
        "Cherenkov-light head on the pointing of the instrument.",
    },
    "magnet_cherenkov_pool_x_m": {
        "dtype": "<f8",
        "comment": "This offset must be added to the core-position, where "
        "the trajectory of the primary particle intersects the "
        "observation-level, in order for the instrument to stand in "
        "the typical center of the Cherenkov-pool.",
    },
    "magnet_cherenkov_pool_y_m": {"dtype": "<f8", "comment": ""},
    "solid_angle_thrown_sr": {"dtype": "<f8", "comment": ""},
    "depth_g_per_cm2": {"dtype": "<f8", "comment": ""},
    "momentum_x_GeV_per_c": {"dtype": "<f8", "comment": ""},
    "momentum_y_GeV_per_c": {"dtype": "<f8", "comment": ""},
    "momentum_z_GeV_per_c": {"dtype": "<f8", "comment": ""},
    "first_interaction_height_asl_m": {"dtype": "<f8", "comment": ""},
    "starting_height_asl_m": {"dtype": "<f8", "comment": ""},
    "starting_x_m": {"dtype": "<f8", "comment": ""},
    "starting_y_m": {"dtype": "<f8", "comment": ""},
}

STRUCTURE["cherenkovsize"] = {
    "num_bunches": {"dtype": "<i8", "comment": ""},
    "num_photons": {"dtype": "<f8", "comment": ""},
}

STRUCTURE["particlepool"] = {
    "num_water_cherenkov": {
        "dtype": "<i8",
        "comment": "The number of particles which reach the observation-level "
        "and will emitt Cherenkov-light in water",
    },
    "num_air_cherenkov": {
        "dtype": "<i8",
        "comment": "Same as 'num_water_cherenkov' but for the air at the "
        "instruments altitude.",
    },
    "num_unknown": {
        "dtype": "<i8",
        "comment": "Particles which are not (yet) in our"
        "corsika-particle-zoo.",
    },
}

STRUCTURE["grid"] = {
    "num_bins_thrown": {
        "dtype": "<i8",
        "comment": "The number of all grid-bins which can collect "
        "Cherenkov-photons.",
    },
    "bin_width_m": {"dtype": "<f8", "comment": ""},
    "field_of_view_radius_deg": {"dtype": "<f8", "comment": ""},
    "field_of_view_azimuth_deg": {
        "dtype": "<f8",
        "comment": "Pointing azimuth w.r.t. magnetic north.",
    },
    "field_of_view_zenith_deg": {
        "dtype": "<f8",
        "comment": "Pointing zenith-distance.",
    },
    "pointing_direction_x": {"dtype": "<f8", "comment": ""},
    "pointing_direction_y": {"dtype": "<f8", "comment": ""},
    "pointing_direction_z": {"dtype": "<f8", "comment": ""},
    "random_shift_x_m": {"dtype": "<f8", "comment": ""},
    "random_shift_y_m": {"dtype": "<f8", "comment": ""},
    "magnet_shift_x_m": {"dtype": "<f8", "comment": ""},
    "magnet_shift_y_m": {"dtype": "<f8", "comment": ""},
    "total_shift_x_m": {
        "dtype": "<f8",
        "comment": "Sum of random and magnetic shift.",
    },
    "total_shift_y_m": {
        "dtype": "<f8",
        "comment": "Sum of random and magnetic shift.",
    },
    "num_bins_above_threshold": {"dtype": "<i8", "comment": ""},
    "overflow_x": {"dtype": "<i8", "comment": ""},
    "overflow_y": {"dtype": "<i8", "comment": ""},
    "underflow_x": {"dtype": "<i8", "comment": ""},
    "underflow_y": {"dtype": "<i8", "comment": ""},
    "area_thrown_m2": {"dtype": "<f8", "comment": ""},
    "artificial_core_limitation": {"dtype": "<i8", "comment": "Flag"},
    "artificial_core_limitation_radius_m": {"dtype": "<f8", "comment": ""},
}

STRUCTURE["cherenkovpool"] = {
    "maximum_asl_m": {"dtype": "<f8", "comment": ""},
    "wavelength_median_nm": {"dtype": "<f8", "comment": ""},
    "cx_median_rad": {"dtype": "<f8", "comment": ""},
    "cy_median_rad": {"dtype": "<f8", "comment": ""},
    "x_median_m": {"dtype": "<f8", "comment": ""},
    "y_median_m": {"dtype": "<f8", "comment": ""},
    "bunch_size_median": {"dtype": "<f8", "comment": ""},
}

STRUCTURE["cherenkovsizepart"] = STRUCTURE["cherenkovsize"].copy()
STRUCTURE["cherenkovpoolpart"] = STRUCTURE["cherenkovpool"].copy()

STRUCTURE["core"] = {
    "bin_idx_x": {"dtype": "<i8", "comment": ""},
    "bin_idx_y": {"dtype": "<i8", "comment": ""},
    "core_x_m": {"dtype": "<f8", "comment": ""},
    "core_y_m": {"dtype": "<f8", "comment": ""},
}

STRUCTURE["particlepoolonaperture"] = {
    "num_air_cherenkov_on_aperture": {
        "dtype": "<i8",
        "comment": "Same as 'num_air_cherenkov' but also run through the "
        "instrument's aperture for Cherenkov-light.",
    },
}

STRUCTURE["instrument"] = {
    "start_time_of_exposure_s": {
        "dtype": "<f8",
        "comment": "The start-time of the instrument's exposure-window"
        "relative to the clock in CORSIKA.",
    },
}
