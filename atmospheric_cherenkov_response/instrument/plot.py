import sebastians_matplotlib_addons as sebplt
import plenopy
import numpy as np


def plot_response(path, toy_instrument, response, vmin=None, vmax=None):
    toy = toy_instrument
    image_response = np.sum(response, axis=0)
    time_response = np.sum(response, axis=1)

    f_m = toy["mirror"]["focal_length_m"]
    r_deg = toy["camera"]["field_of_view_half_angle_deg"]

    fig = sebplt.figure({"rows": 1280, "cols": 1280, "fontsize": 1})
    ax_img = sebplt.add_axes(fig, [0.15, 0.3, 0.7, 0.7])
    plenopy.plot.image.add2ax(
        ax=ax_img,
        I=image_response,
        px=np.rad2deg(np.arctan(toy["camera"]["pixel"]["x_m"] / f_m)),
        py=np.rad2deg(np.arctan(toy["camera"]["pixel"]["y_m"] / f_m)),
        colormap="magma_r",
        hexrotation=0,
        vmin=vmin,
        vmax=vmax,
        colorbar=True,
        norm=None,
    )
    oh = 1.05
    roh_deg = oh * r_deg
    ax_img.set_xlim([-roh_deg, roh_deg])
    ax_img.set_ylim([-roh_deg, roh_deg])
    ax_img.set_xlabel(r"$c_x\,/\,1^{\circ}$")
    ax_img.set_ylabel(r"$c_y\,/\,1^{\circ}$")
    ax_img.set_aspect("equal")

    ax_time = sebplt.add_axes(fig, [0.15, 0.1, 0.7, 0.1])
    sebplt.ax_add_histogram(
        ax=ax_time,
        bin_edges=1e9 * toy["camera"]["time_slices"]["edges_s"],
        bincounts=time_response,
        linestyle="-",
        linecolor="k",
        linealpha=1.0,
        draw_bin_walls=True,
    )
    ax_time.set_xlabel(r"time$\,/\,$ns")
    fig.savefig(path)
    sebplt.close(fig)
