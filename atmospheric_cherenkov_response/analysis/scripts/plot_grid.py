import svg_cartesian_plot as svgplt
import numpy as np
import os

fig = svgplt.Fig(cols=1280, rows=1280)

ax = svgplt.Ax(fig=fig)
ax["span"] = (0.0, 0.0, 1, 1)
ax["xlim"] = [-1, 1]
ax["ylim"] = [-1, 1]

NUM_HORIZONTAL = 48
NUM_VERTICAL = 24
FLUCHT_POINT = [0, 0]


def line_line_intersection(a_start_xy, a_stop_xy, b_start_xy, b_stop_xy):
    a_start_xy = np.asarray(a_start_xy)
    a_stop_xy = np.asarray(a_stop_xy)
    b_start_xy = np.asarray(b_start_xy)
    b_stop_xy = np.asarray(b_stop_xy)

    a_sup = a_start_xy
    a_dir = a_stop_xy - a_sup

    a_m = a_dir[1] / a_dir[0]
    a_b = a_sup[1] - a_m * a_sup[0]

    b_sup = b_start_xy
    b_dir = b_stop_xy - b_sup

    b_m = b_dir[1] / b_dir[0]
    b_b = b_sup[1] - b_m * b_sup[0]

    x_inter = (b_b - a_b) / (a_m - b_m)
    y_inter = a_m * x_inter + a_b

    return np.asarray([x_inter, y_inter])


Y_LIM = -0.01

x_ticks = np.linspace(
    -NUM_HORIZONTAL / 2, NUM_HORIZONTAL / 2, NUM_HORIZONTAL * 2
)
y_ticks = np.geomspace(Y_LIM, -0.99, NUM_VERTICAL)

xy_ylim_start = [-1, Y_LIM]
xy_ylim_stop = [1, Y_LIM]

"""
for xflucht in np.linspace(-NUM_HORIZONTAL/4, NUM_HORIZONTAL/4, NUM_HORIZONTAL):
    xy_flucht_start = [xflucht, -1]
    xy_flucht_final = FLUCHT_POINT

    xy_stop = line_line_intersection(
        a_start_xy=xy_flucht_start,
        a_stop_xy=xy_flucht_final,
        b_start_xy=xy_ylim_start,
        b_stop_xy=xy_ylim_stop,
    )

    svgplt.ax_add_path(ax=ax, xy=[xy_flucht_start, xy_stop],
        stroke=svgplt.color.css("black"),
    )
"""

xy_flucht_final = FLUCHT_POINT
for xx in range(len(x_ticks)):
    xy_flucht_start = [x_ticks[xx], -1]

    for yy in range(len(y_ticks) - 1):
        y_low = y_ticks[yy]
        y_upp = y_ticks[yy + 1]

        xy_hor_low_start = [-1, y_low]
        xy_hor_low_stop = [1, y_low]

        xy_hor_upp_start = [-1, y_upp]
        xy_hor_upp_stop = [1, y_upp]

        xy_start = line_line_intersection(
            a_start_xy=xy_flucht_start,
            a_stop_xy=xy_flucht_final,
            b_start_xy=xy_hor_low_start,
            b_stop_xy=xy_hor_low_stop,
        )
        xy_stop = line_line_intersection(
            a_start_xy=xy_flucht_start,
            a_stop_xy=xy_flucht_final,
            b_start_xy=xy_hor_upp_start,
            b_stop_xy=xy_hor_upp_stop,
        )

        if -1 <= xy_start[0] <= 1 or -1 <= xy_stop[0] <= 1:
            svgplt.ax_add_path(
                ax=ax,
                xy=[xy_start, xy_stop],
                stroke=svgplt.color.css("black"),
                stroke_width=np.abs(y_low) * 2,
            )


for ypos in y_ticks:
    xy_start = [-1, ypos]
    xy_stop = [+1, ypos]
    svgplt.ax_add_path(
        ax=ax,
        xy=[xy_start, xy_stop],
        stroke=svgplt.color.css("black"),
        stroke_width=np.abs(ypos) * 2,
    )


svgplt.fig_write(fig=fig, path="grid.svg")
