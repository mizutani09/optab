import argparse

import h5py
import matplotlib.pyplot as plt
import numpy as np


def read_table(h5_path):
    with h5py.File(h5_path, "r") as f:
        logt = f["axes/log10_temperature"][:]
        logrho = f["axes/log10_density"][:]
        logk_ross = f["kappa/log10_rosseland"][:]
        logk_pla = f["kappa/log10_planck"][:]
    return logt, logrho, logk_ross, logk_pla


def nearest_index(array, value):
    return int(np.argmin(np.abs(array - value)))


# def pick_color_limits(mean, logk, vmin=None, vmax=None, qmin=1.0, qmax=99.0):
#     if vmin is not None and vmax is not None:
#         return float(vmin), float(vmax)
#     vals = logk[np.isfinite(logk)]
#     cmin = np.percentile(vals, qmin) if vmin is None else vmin
#     cmax = np.percentile(vals, qmax) if vmax is None else vmax
#     if cmax <= cmin:
#         cmax = cmin + 1e-6
#     # Keep the same range style as plot_opacity_table_hdf5.py
#     cmin = -5
#     cmax = 4
#     return float(cmin), float(cmax)


def draw_colormap(ax, logt, logrho, logk, mean, vmin=None, vmax=None, qmin=1.0, qmax=99.0):
    t = 10.0 ** logt
    rho = 10.0 ** logrho
    # cmin, cmax = pick_color_limits(mean, logk, vmin=vmin, vmax=vmax, qmin=qmin, qmax=qmax)
    if mean == "ross":
        cmin, cmax = -6, 6
    else:
        cmin, cmax = -2, 7
    mesh = ax.pcolormesh(rho, t, logk.T, shading="auto", cmap="jet", vmin=cmin, vmax=cmax)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.tick_params(axis="both", which="major", labelsize=18)
    ax.set_xlabel(r"$\rho\,[\mathrm{g\,cm^{-3}}]$", fontsize=20)
    ax.set_ylabel(r"$T\,[\mathrm{K}]$", fontsize=20)
    # if mean == "ross":
    #     ax.set_title("Rosseland colormap", fontsize=14)
    # else:
    #     ax.set_title("Planck colormap", fontsize=14)
    ax.grid(True, alpha=0.3)
    return mesh


def draw_kappa_vs_t(ax, logt, logrho, logk, log_rho_targets, mean):
    t = 10.0 ** logt
    for target in log_rho_targets:
        j = nearest_index(logrho, target)
        label = r"$\rho=10^{" + f"{logrho[j]:.0f}" + r"}\,\mathrm{g\,cm^{-3}}$"
        ax.plot(t, logk[j, :], marker="o", ms=3.0, lw=1.4, label=label)

    ax.tick_params(axis="both", which="major", labelsize=18)
    ax.set_xscale("log")
    ax.set_ylim(-5, 4)
    ax.set_xlabel(r"$T\,[\mathrm{K}]$", fontsize=20)
    if mean == "ross":
        ylabel = r"$\log_{10}\,\kappa_\mathrm{R}\,[\mathrm{cm^2\,g^{-1}}]$"
        # ax.set_title("Rosseland at fixed density", fontsize=14)
    else:
        ylabel = r"$\log_{10}\,\kappa_\mathrm{P}\,[\mathrm{cm^2\,g^{-1}}]$"
        # ax.set_title("Planck at fixed density", fontsize=14)
    ax.set_ylabel(ylabel, fontsize=20)
    ax.grid(True, alpha=0.3)
    if mean == "ross":
        # ax.legend(fontsize=14)
        # make it in two row
        ax.legend(fontsize=12)


def main():
    parser = argparse.ArgumentParser(
        description="Create a 4-panel figure (ross/pla x colormap/vsT) from opacity table HDF5."
    )
    parser.add_argument("h5_file", help="Input opacity table HDF5 path")
    parser.add_argument(
        "--log-rho-targets",
        type=float,
        nargs="+",
        default=[-12, -11, -10, -9, -8],
        help="Target log10(rho[g/cm^3]) values for kappa(T) curves",
    )
    parser.add_argument(
        "--out",
        default="opacity_table_4panel.png",
        help="Output PNG filename",
    )
    parser.add_argument("--vmin", type=float, default=None, help="Fixed vmin for colormap")
    parser.add_argument("--vmax", type=float, default=None, help="Fixed vmax for colormap")
    parser.add_argument("--qmin", type=float, default=1.0, help="Lower percentile for auto color range")
    parser.add_argument("--qmax", type=float, default=99.0, help="Upper percentile for auto color range")
    args = parser.parse_args()

    logt, logrho, logk_ross, logk_pla = read_table(args.h5_file)

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    mesh_r = draw_colormap(
        axes[0, 0], logt, logrho, logk_ross, "ross",
        vmin=args.vmin, vmax=args.vmax, qmin=args.qmin, qmax=args.qmax
    )
    mesh_p = draw_colormap(
        axes[0, 1], logt, logrho, logk_pla, "pla",
        vmin=args.vmin, vmax=args.vmax, qmin=args.qmin, qmax=args.qmax
    )

    draw_kappa_vs_t(axes[1, 0], logt, logrho, logk_ross, args.log_rho_targets, "ross")
    draw_kappa_vs_t(axes[1, 1], logt, logrho, logk_pla, args.log_rho_targets, "pla")

    cbar_r = fig.colorbar(mesh_r, ax=axes[0, 0], fraction=0.046, pad=0.04)
    cbar_p = fig.colorbar(mesh_p, ax=axes[0, 1], fraction=0.046, pad=0.04)
    cbar_r.set_label(r"$\log_{10}\,\kappa_\mathrm{R}\,[\mathrm{cm^2\,g^{-1}}]$", fontsize=20)
    cbar_p.set_label(r"$\log_{10}\,\kappa_\mathrm{P}\,[\mathrm{cm^2\,g^{-1}}]$", fontsize=20)
    cbar_r.ax.tick_params(labelsize=18)
    cbar_p.ax.tick_params(labelsize=18)

    plt.tight_layout()
    plt.savefig(args.out, dpi=180)
    plt.close(fig)
    print(f"Input: {args.h5_file}")
    print(f"Output: {args.out}")


if __name__ == "__main__":
    main()
