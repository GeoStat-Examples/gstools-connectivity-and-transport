# -*- coding: utf-8 -*-
"""
Transport example using GSTools.

Plotting the plumes at t=15d.

Authors: Alraune Zech and Sebastian Müller
"""
import os
import numpy as np
import meshio as mio
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter


out_path = os.path.join("..", "results", "tracer_test_conn")
connectivity = ["mean", "low", "high"]
conn_count = len(connectivity)

fac = 3  # ratio between axis and colorbar
fig = plt.figure(constrained_layout=False, figsize=[6.4, 10])
gs = fig.add_gridspec(conn_count * fac + 2, 1)

ax_fld = {}
for i, conn in enumerate(connectivity):
    axs = ax_fld.get("mean", None)
    ax_fld[conn] = fig.add_subplot(gs[fac * i : fac * (i + 1), 0], sharex=axs)
    plt.setp(ax_fld[conn].get_xticklabels(), visible=(i == conn_count - 1))

cax1 = fig.add_subplot(gs[-2, 0])
cax2 = fig.add_subplot(gs[-1, 0])

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
    ax_fld[conn].plot([0, 0], [-2, -4], color="k", alpha=0.7, linewidth=3)
    tri = ax_fld[conn].tricontourf(
        *fld_c[conn], levels=lvl_c, cmap="Reds", alpha=0.7
    )
ax_fld["high"].set_xlabel("x / m")

c1 = fig.colorbar(im, orientation="horizontal", cax=cax1, shrink=0.2, pad=0.5)
c2 = fig.colorbar(tri, orientation="horizontal", cax=cax2, shrink=0.2, pad=0.5)
c1.ax.set_xlabel(r"Transmissivity / m s$^{-1}$")
c2.ax.set_xlabel(r"Concentration")

# to actually determine the labels, draw the figure
fig.canvas.draw()
lbl_t = [f"$10^{{{lbl.get_text()}}}$" for lbl in c1.ax.get_xticklabels()]
c1.ax.set_xticklabels(lbl_t)
c2.ax.xaxis.set_major_formatter(FormatStrFormatter("%.3f"))
fig.tight_layout()
fig.subplots_adjust(right=0.95)
fig.savefig(os.path.join("..", "results", "comparison_simple.pdf"), dpi=300)
