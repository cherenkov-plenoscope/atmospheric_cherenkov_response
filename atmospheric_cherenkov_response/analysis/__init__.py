import numpy as np


def effective_quantity_for_grid(
    mask_detected,
    quantity_scatter,
    num_grid_cells_above_lose_threshold,
    total_num_grid_cells,
):
    """
    Returns the effective quantity and its absolute uncertainty.

    Parameters
    ----------
    mask_detected : array(num. thrown airshower)
        A flag/weight for each airshower marking its detection.
    quantity_scatter : array(num. thrown airshower)
        The scatter-quantity for each airshower. This is area/m^2 for point
        like sources, or acceptance/m^2 sr for diffuse sources.
    num_grid_cells_above_lose_threshold : array(num. thrown airshower)
        Num. of grid cells passing the lose threshold of the grid for each
        airshower.
    total_num_grid_cells : array(num. thrown airshower)
        The total number of grid-cells thrown for each airshower.

    Returns
    -------
    (effective_quantity, effective_quantity abs. uncertainty) : (float, float)

    Formula
    -------
    Q_effective &=& /frac{ Q_detected }{ C_thrown }

    Q_detected &=& /sum_m^M { f_{detected,m} N_{S,m} Q_{scatter,m} }

    C_thrown &=& /sum_m^M { f_{thrown,m} N_G }

    /frac{
        /Delta_Q_effective
    }{
        Q_effective
    } &=& /frac{
        /sqrt{ /sum_m^M { f_{detected,m} N_{S,m}^2 } }
    }{
        /sum_m^M { f_{detected,m} N_{S,m} }
    }

    Variables
    ---------
    N_G             Num. bins in grid.

    N_{S,m}         Num. bins in grid with cherenkov-photons
                    above losethreshold for m-th air-shower.

    Q_{scatter,m}   Scatter-quantity of m-th air-shower.
                    Scatter-area times scatter-solid-angle

    f_{thrown,m}    Flag marking that m-th air-shower is a valid thrown.

    f_{detected,m}  Flag marking that m-th air-shower is a valid detection.

    """

    # sanity checks
    # -------------
    assert np.all(mask_detected >= 0)
    assert np.all(quantity_scatter >= 0.0)
    assert np.all(num_grid_cells_above_lose_threshold >= 0)
    assert np.all(total_num_grid_cells > 0)

    _num_events = len(mask_detected)
    assert _num_events == len(mask_detected)
    assert _num_events == len(quantity_scatter)
    assert _num_events == len(num_grid_cells_above_lose_threshold)
    assert _num_events == len(total_num_grid_cells)

    f_detected = np.asarray(mask_detected, dtype=bool)
    Q_scatter = np.asarray(quantity_scatter, dtype=float)
    N_S = np.asarray(num_grid_cells_above_lose_threshold, dtype=np.int64)
    N_G = np.asarray(total_num_grid_cells, dtype=np.int64)

    # mean
    # ----
    Q_detected = np.sum(f_detected * N_S * Q_scatter)
    C_thrown = np.sum(N_G)

    if C_thrown == 0:
        Q_effective = 0.0
    else:
        Q_effective = Q_detected / C_thrown

    # uncertainty
    # -----------
    # according to Werner EffAreaComment.pdf 2020-03-21 17:35

    Sum_f_x_N_S_x_N_S = np.sum(f_detected * N_S**2)
    assert Sum_f_x_N_S_x_N_S >= 0.0

    Sum_f_x_N_S = np.sum(f_detected * N_S)
    assert Sum_f_x_N_S >= 0.0

    if Sum_f_x_N_S == 0.0:
        Q_effective_relunc = 0.0
    else:
        Q_effective_relunc = np.sqrt(Sum_f_x_N_S_x_N_S) / Sum_f_x_N_S

    Q_effective_absunc = Q_effective * Q_effective_relunc
    return Q_effective, Q_effective_absunc


def effective_quantity_for_grid_vs_energy(
    energy_bin_edges_GeV,
    energy_GeV,
    mask_detected,
    quantity_scatter,
    num_grid_cells_above_lose_threshold,
    total_num_grid_cells,
):
    """
    Returns the effective quantity and its absolute uncertainty vs. the energy.

    Parameters
    ----------
    energy_bin_edges_GeV : array like

    energy_GeV : array like (num. thrown airshower)
        The energy of each airshower.
    mask_detected : array(num. thrown airshower)
        A flag/weight for each airshower marking its detection.
    quantity_scatter : array(num. thrown airshower)
        The scatter-quantity for each airshower. This is area/m^2 for point
        like sources, or acceptance/m^2 sr for diffuse sources.
    num_grid_cells_above_lose_threshold : array(num. thrown airshower)
        Num. of grid cells passing the lose threshold of the grid for each
        airshower.
    total_num_grid_cells : array(num. thrown airshower)
        The total number of grid-cells thrown for each airshower.

    Returns
    -------
    (effective_quantity, effective_quantity abs. uncertainty) : (array, array)
        Number of energy bins is len(energy_bin_edges_GeV) - 1.

    Formula
    -------
    Q_effective &=& /frac{ Q_detected }{ C_thrown }

    Q_detected &=& /sum_m^M { f_{detected,m} N_{S,m} Q_{scatter,m} }

    C_thrown &=& /sum_m^M { f_{thrown,m} N_G }

    /frac{
        /Delta_Q_effective
    }{
        Q_effective
    } &=& /frac{
        /sqrt{ /sum_m^M { f_{detected,m} N_{S,m}^2 } }
    }{
        /sum_m^M { f_{detected,m} N_{S,m} }
    }

    Variables
    ---------
    N_G             Num. bins in grid.

    N_{S,m}         Num. bins in grid with cherenkov-photons
                    above losethreshold for m-th air-shower.

    Q_{scatter,m}   Scatter-quantity of m-th air-shower.
                    Scatter-area times scatter-solid-angle

    f_{thrown,m}    Flag marking that m-th air-shower is a valid thrown.

    f_{detected,m}  Flag marking that m-th air-shower is a valid detection.

    """
    assert np.all(energy_bin_edges_GeV > 0)
    assert np.all(np.gradient(energy_bin_edges_GeV) > 0)
    assert np.all(energy_GeV > 0.0)

    num_energy_bins = len(energy_bin_edges_GeV) - 1
    assert num_energy_bins >= 1

    Q_effective_vs_energy = np.zeros(num_energy_bins)
    Q_effective_absunc_vs_energy = np.zeros(num_energy_bins)

    for enbin in range(num_energy_bins):
        energy_start_GeV = energy_bin_edges_GeV[enbin]
        energy_stop_GeV = energy_bin_edges_GeV[ebin + 1]
        energy_mask = np.logical_and(
            energy_start_GeV <= energy_GeV, energy_stop_GeV > energy_GeV
        )

        (Q_effective_vs_energy[enbin], Q_effective_absunc_vs_energy[enbin]) = (
            effective_quantity_for_grid(
                mask_detected=mask_detected[energy_mask],
                quantity_scatter=quantity_scatter[energy_mask],
                num_grid_cells_above_lose_threshold=num_grid_cells_above_lose_threshold[
                    energy_mask
                ],
                total_num_grid_cells=total_num_grid_cells[energy_mask],
            )
        )

    return Q_effective_vs_energy, Q_effective_absunc_vs_energy

    """
    assert np.all(mask_detected >= 0)

    assert np.all(quantity_scatter >= 0.0)

    assert np.all(num_grid_cells_above_lose_threshold >= 0)

    assert np.all(total_num_grid_cells > 0)

    num_grid_cells_above_lose_threshold = np.asarray(
        num_grid_cells_above_lose_threshold, dtype=np.int64
    )
    total_num_grid_cells = np.asarray(total_num_grid_cells, dtype=np.int64)

    quantity_detected = np.histogram(
        energy_GeV,
        bins=energy_bin_edges_GeV,
        weights=(
            mask_detected
            * num_grid_cells_above_lose_threshold
            * quantity_scatter
        ),
    )[0]
    assert np.all(quantity_detected >= 0.0)

    count_thrown = np.histogram(
        energy_GeV,
        weights=total_num_grid_cells,
        bins=energy_bin_edges_GeV,
    )[0]
    assert np.all(count_thrown >= 0.0)

    effective_quantity = divide_silent(
        numerator=quantity_detected, denominator=count_thrown, default=0.0
    )

    # uncertainty
    # according to Werner EffAreaComment.pdf 2020-03-21 17:35

    A_square = np.histogram(
        energy_GeV,
        bins=energy_bin_edges_GeV,
        weights=(mask_detected * num_grid_cells_above_lose_threshold**2),
    )[0]
    assert np.all(A_square >= 0.0)

    A = np.histogram(
        energy_GeV,
        bins=energy_bin_edges_GeV,
        weights=(mask_detected * num_grid_cells_above_lose_threshold),
    )[0]
    assert np.all(A >= 0.0)

    effective_quantity_relative_uncertainty = divide_silent(
        numerator=np.sqrt(A_square), denominator=A, default=0.0
    )

    effective_quantity_absolute_uncertainty = (
        effective_quantity * effective_quantity_relative_uncertainty
    )

    return effective_quantity, effective_quantity_absolute_uncertainty
    """


"""
def divide_silent(numerator, denominator, default):
    valid = denominator != 0
    division = np.ones(shape=numerator.shape) * default
    division[valid] = numerator[valid] / denominator[valid]
    return division

"""
