import argparse
import os
import binning_utils
import sebastians_matplotlib_addons as sebplt
import atmospheric_cherenkov_response
from atmospheric_cherenkov_response import plot
import numpy as np


parser = argparse.ArgumentParser(
    prog="plot_allsky_cherenkov_density.py",
    description=("Make plots of Cherenkov statistics"),
)
parser.add_argument(
    "--out_path",
    metavar="PATH",
    type=str,
    help="Path of figure.",
)

args = parser.parse_args()

PLT = atmospheric_cherenkov_response.plot.config()
sebplt.matplotlib.rcParams.update(PLT["rcParams"])

_NUM_STEPS = 101
Nacc = np.floor(np.geomspace(1e3, 1e5, _NUM_STEPS)).astype(int)
Nacc = sorted(list(set(Nacc)))
Nacc = np.asarray(Nacc)
NUM_STEPS = len(Nacc)

prng = np.random.Generator(np.random.PCG64(100))

deltaQ_over_Q_disk_method = np.sqrt(Nacc) / Nacc

deltaQ_over_Q_grid_method_Nm_const = np.zeros(NUM_STEPS)
deltaQ_over_Q_grid_method_Nm_normal_25 = np.zeros(NUM_STEPS)
deltaQ_over_Q_grid_method_Nm_normal_50 = np.zeros(NUM_STEPS)
deltaQ_over_Q_grid_method_Nm_normal_100 = np.zeros(NUM_STEPS)


def integers_normal(prng, low, loc, scale, size, max_num_iterations=10000):
    n = np.zeros(size, dtype=int)

    i = 0
    while np.any(n < low):
        assert i < max_num_iterations
        mask = n < 1
        n[mask] = prng.normal(loc=loc, scale=scale, size=np.sum(mask)).astype(
            int
        )
        i += 1
    return n


def evaluate_dQoverQ(N_sm):
    return np.sqrt(np.sum(N_sm**2)) / np.sum(N_sm)


EXPECTATION_VALUE_N_SM = 50
for i in range(NUM_STEPS):
    N_sm_const = EXPECTATION_VALUE_N_SM * np.ones(Nacc[i])
    deltaQ_over_Q_grid_method_Nm_const[i] = evaluate_dQoverQ(N_sm=N_sm_const)

    N_sm_normal = integers_normal(
        prng=prng,
        low=1,
        loc=EXPECTATION_VALUE_N_SM,
        scale=(1 / 4) * EXPECTATION_VALUE_N_SM,
        size=Nacc[i],
    )
    deltaQ_over_Q_grid_method_Nm_normal_25[i] = evaluate_dQoverQ(
        N_sm=N_sm_normal
    )

    N_sm_normal = integers_normal(
        prng=prng,
        low=1,
        loc=EXPECTATION_VALUE_N_SM,
        scale=(1 / 2) * EXPECTATION_VALUE_N_SM,
        size=Nacc[i],
    )
    deltaQ_over_Q_grid_method_Nm_normal_50[i] = evaluate_dQoverQ(
        N_sm=N_sm_normal
    )

    N_sm_normal = integers_normal(
        prng=prng,
        low=1,
        loc=EXPECTATION_VALUE_N_SM,
        scale=EXPECTATION_VALUE_N_SM,
        size=Nacc[i],
    )
    deltaQ_over_Q_grid_method_Nm_normal_100[i] = evaluate_dQoverQ(
        N_sm=N_sm_normal
    )


fig = sebplt.figure({"rows": 1080, "cols": 1920, "fontsize": 1.875})
ax = sebplt.add_axes(fig=fig, span=(0.175, 0.2, 0.8, 0.75))
ax.plot(
    Nacc, deltaQ_over_Q_grid_method_Nm_const / deltaQ_over_Q_disk_method, "k--"
)
ax.plot(
    Nacc,
    deltaQ_over_Q_grid_method_Nm_normal_25 / deltaQ_over_Q_disk_method,
    "r--",
)
ax.plot(
    Nacc,
    deltaQ_over_Q_grid_method_Nm_normal_50 / deltaQ_over_Q_disk_method,
    "g--",
)
ax.plot(
    Nacc,
    deltaQ_over_Q_grid_method_Nm_normal_100 / deltaQ_over_Q_disk_method,
    "b--",
)
ax.semilogx()
ax.set_xlabel("number of accepted showers / 1")
ax.set_ylabel("relative uncertainty\n(disk method)/(grid method) / 1")
ax.set_ylim(0.95, 1.35)
ax.set_xlim(min(Nacc), max(Nacc))
fig.savefig(args.out_path)
sebplt.close(fig)


num_normal_shapes = 51
scale_over_locs = np.linspace(0, 2, num_normal_shapes)
Q_runc = np.zeros(num_normal_shapes)

for ns in range(num_normal_shapes):
    loc = EXPECTATION_VALUE_N_SM
    scale = scale_over_locs[ns] * loc

    dQ_Q = np.zeros(NUM_STEPS)
    for i in range(NUM_STEPS):
        N_sm_normal = integers_normal(
            prng=prng, low=1, loc=loc, scale=scale, size=Nacc[i]
        )
        dQ_Q[i] = evaluate_dQoverQ(N_sm=N_sm_normal)

    Q_runc[ns] = np.median(dQ_Q / deltaQ_over_Q_disk_method)


fig = sebplt.figure({"rows": 1080, "cols": 1920, "fontsize": 1.875})
ax = sebplt.add_axes(fig=fig, span=(0.175, 0.2, 0.8, 0.75))
ax.plot(scale_over_locs, Q_runc, "k-")
ax.set_xlabel(r"Var($N_\mathrm{S}$)/$\langle N_\mathrm{S} \rangle$  / 1")
ax.set_ylabel("relative uncertainty\n(disk method)/(grid method) / 1")
ax.set_ylim(0.95, 1.35)
ax.set_xlim(min(scale_over_locs), max(scale_over_locs))
fig.savefig(args.out_path)
sebplt.close(fig)
