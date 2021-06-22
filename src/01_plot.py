# -*- coding: utf-8 -*-
"""
Transport example using GSTools.

Plotting the plumes at t=15d and calculating the breakthrough curves at
the observation wells.

Authors: Alraune Zech and Sebastian Müller
"""
import os
import numpy as np
from ogs5py.reader import readtec_polyline
import meshio as mio
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter


def dashes(i=1, max_n=6, width=1):
    """Dashes for matplotlib."""
    return i * [width, width] + [max_n * 2 * width - 2 * i * width, width]


out_path = os.path.join("..", "results", "tracer_test_conn")
ow_file_name = "model_ply_ow{0}_t{0}_MASS_TRANSPORT.tec"
connectivity = ["mean", "low", "high"]
conn_count = len(connectivity)
ow_pos = [2.5, 5.5, 8.5]
ow_count = len(ow_pos)
ow_label = ", ".join([f"ow{j+1}" for j in range(ow_count)])

sub_factor = 2
fig = plt.figure(constrained_layout=False, figsize=[8, 9])
gs = fig.add_gridspec(conn_count + 1, ow_count * sub_factor, wspace=0)

ax_fld = {}
for i, conn in enumerate(connectivity):
    ax_shr = None if i == 0 else ax_fld["mean"]
    ax_fld[conn] = fig.add_subplot(gs[i, :-2], sharex=ax_shr)
    plt.setp(ax_fld[conn].get_xticklabels(), visible=(i == conn_count - 1))
    ax_fld[conn].margins(x=0)

ax = [fig.add_subplot(gs[-1, 0:sub_factor])]
for i in range(1, ow_count):
    ax.append(
        fig.add_subplot(
            gs[-1, sub_factor * i : sub_factor * (i + 1)],
            sharey=ax[0],
            sharex=ax[0],
        )
    )
    plt.setp(ax[-1].get_yticklabels(), visible=False)

cax1 = fig.add_subplot(gs[:-1, -2])
cax1.axis("off")
cax2 = fig.add_subplot(gs[:-1, -1])
cax2.axis("off")

fld_t = {}
fld_c = {}
bnd_t = [np.inf, -np.inf]
max_c = -np.inf
for conn in connectivity:
    # transmissivity
    fld = mio.read(os.path.join(out_path, f"trans_field_{conn}.vtk"))
    fld_t[conn] = np.log10(fld.cell_data["transmissivity"][0].reshape(40, 400))
    bnd_t[0] = min(bnd_t[0], fld_t[conn].min())
    bnd_t[1] = max(bnd_t[1], fld_t[conn].max())
    # concentration at t=15d (step 180)
    fld = mio.read(os.path.join(out_path, f"out_{conn}", "model0180.vtk"))
    x, z = fld.points[:, [0, 2]].T
    c = fld.point_data["CONCENTRATION1"][:, 0]
    fld_c[conn] = [x, z, c]
    max_c = max(max_c, c.max())

# low-cut for concentration at 0.001
lvl_c = np.linspace(0.001, max_c, 16)
kw_t = dict(vmin=bnd_t[0], vmax=bnd_t[1], interpolation="bicubic")

for conn in connectivity:
    im = ax_fld[conn].imshow(
        fld_t[conn], origin="lower", extent=[-1, 9, -5, -1], **kw_t
    )
    ax_fld[conn].set_ylabel("z / m")
    ax_fld[conn].set_title(f"Connected region: '{conn}' transmissivity")
    # plot source line
    lbl = "source" if conn == "mean" else None
    ax_fld[conn].plot(
        [0, 0], [-2, -4], color="k", alpha=0.7, linewidth=3, label=lbl
    )
    # plot observation wells
    for i, ow in enumerate(ow_pos, start=1):
        lbl = f"ow{i}"
        ax_fld[conn].axvline(
            ow, color="k", alpha=0.7, linewidth=2, dashes=dashes(i), label=lbl
        )
    x, z, c = fld_c[conn]
    tri = ax_fld[conn].tricontourf(
        x, z, c, levels=lvl_c, cmap="Reds", alpha=0.7
    )
ax_fld["high"].set_xlabel("x / m")
handles, labels = ax_fld["mean"].get_legend_handles_labels()
fig.legend(
    handles, labels, loc="upper right", ncol=2, title="Plumes for t=15d"
)

c1 = fig.colorbar(im, ax=cax1, use_gridspec=False, pad=0, fraction=0.25)
c2 = fig.colorbar(tri, ax=cax2, use_gridspec=False, pad=0, fraction=0.25)
c1.ax.set_ylabel(r"transmissivity / $\frac{m^2}{s}$")
c1.ax.yaxis.set_label_position("left")
c2.ax.set_ylabel(r"concentration")
c2.ax.yaxis.set_label_position("left")
c1.ax.margins(x=0)
c2.ax.margins(x=0)

# to actually determine the labels, draw the figure
fig.canvas.draw()
lbl_t = c1.ax.get_yticklabels()
lbl_t = [f"$10^{{{lbl.get_text()}}}$" for lbl in lbl_t]
c1.ax.set_yticklabels(lbl_t)
c2.ax.yaxis.set_major_formatter(FormatStrFormatter("%.3f"))

for i, conn in enumerate(connectivity):
    for j in range(ow_count):
        name = os.path.join(f"out_{conn}", ow_file_name.format(j + 1))
        out = readtec_polyline(single_file=os.path.join(out_path, name))
        t = out["TIME"] / 3600 / 24  # convert secs to days
        # mean concentration along observation well
        c = np.trapz(out["CONCENTRATION1"], x=out["DIST"]) / out["DIST"][:, -1]
        ax[j].plot(np.cumsum(c), t, label=f"{conn}")
        ax[j].set_xticks([])
        ax[j].invert_yaxis()
        ax[j].set_xlabel(f"breakthrough curves: ow{j+1}")

ax[-1].legend(loc="upper right")
for j in range(ow_count):
    ax[j].axhline(15, color="k", alpha=0.7, linestyle=":")

ax[0].set_ylabel("time / d")
fig.tight_layout()
fig.savefig(os.path.join("..", "results", "comparison.pdf"), dpi=300)
