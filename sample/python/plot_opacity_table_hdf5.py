import argparse

import h5py
import matplotlib.pyplot as plt
import numpy as np


def read_table(h5_path, mean):
    with h5py.File(h5_path, "r") as f:
        logt = f["axes/log10_temperature"][:]
        logrho = f["axes/log10_density"][:]
        if mean == "ross":
            logk = f["kappa/log10_rosseland"][:]
            title = "Rosseland-mean opacity"
        elif mean == "pla":
            logk = f["kappa/log10_planck"][:]
            title = "Planck-mean opacity"
        else:
            raise ValueError(f"Unknown mean: {mean}")
    return logt, logrho, logk, title


def nearest_index(array, value):
    return int(np.argmin(np.abs(array - value)))


def pick_color_limits(logk, vmin=None, vmax=None, qmin=1.0, qmax=99.0):
    if vmin is not None and vmax is not None:
        return vmin, vmax
    vals = logk[np.isfinite(logk)]
    cmin = np.percentile(vals, qmin) if vmin is None else vmin
    cmax = np.percentile(vals, qmax) if vmax is None else vmax
    if cmax <= cmin:
        cmax = cmin + 1e-6
    cmin = -5
    cmax = +4
    return float(cmin), float(cmax)


def plot_colormap(logt, logrho, logk, mean, out_png, vmin=None, vmax=None, qmin=1.0, qmax=99.0):
    t = 10.0 ** logt
    rho = 10.0 ** logrho

    cmin, cmax = pick_color_limits(logk, vmin=vmin, vmax=vmax, qmin=qmin, qmax=qmax)
    fig, ax = plt.subplots(figsize=(8, 6))
    mesh = ax.pcolormesh(rho, t, logk.T, shading="auto", cmap="jet", vmin=cmin, vmax=cmax)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.tick_params(axis="both", which="major", labelsize=18)
    ax.set_xlabel(r"$\rho\,[\mathrm{g\,cm^{-3}}]$", fontsize=20)
    ax.set_ylabel(r"$T\,[\mathrm{K}]$", fontsize=20)
    # ax.set_title(title, fontsize=18)
    cbar = plt.colorbar(mesh, ax=ax)
    if mean == "ross":
        clabel = r"$\log_{10}\,\kappa_\mathrm{R}\,[\mathrm{cm^2\,g^{-1}}]$"
    else:
        clabel = r"$\log_{10}\,\kappa_\mathrm{P}\,[\mathrm{cm^2\,g^{-1}}]$"
    cbar.set_label(clabel, fontsize=18)
    cbar.ax.tick_params(labelsize=14)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_png, dpi=160)
    plt.close(fig)


def plot_kappa_vs_t(logt, logrho, logk, log_rho_targets, out_png, mean):
    t = 10.0 ** logt

    fig, ax = plt.subplots(figsize=(8, 6))
    for target in log_rho_targets:
        j = nearest_index(logrho, target)
        label = (
            r"$\rho=10^{" + f"{logrho[j]:.0f}" + r"}\,\mathrm{g\,cm^{-3}}$"
        )
        ax.plot(t, logk[j, :], marker="o", ms=3.0, lw=1.4, label=label)

    ax.tick_params(axis="both", which="major", labelsize=18)
    ax.set_xscale("log")
    ax.set_ylim(-5, 4)
    ax.set_xlabel(r"$T\,[\mathrm{K}]$", fontsize=20)
    if mean == "ross":
        ylabel = r"$\log_{10}\,\kappa_\mathrm{R}\,[\mathrm{cm^2\,g^{-1}}]$"
    else:
        ylabel = r"$\log_{10}\,\kappa_\mathrm{P}\,[\mathrm{cm^2\,g^{-1}}]$"
    ax.set_ylabel(ylabel, fontsize=20)
    # ax.set_title(title, fontsize=18)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=14)
    plt.tight_layout()
    plt.savefig(out_png, dpi=160)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="Plot opacity table HDF5 produced by make_opacity_table_hdf5.py"
    )
    parser.add_argument("h5_file", help="Input opacity table HDF5 path")
    parser.add_argument(
        "--mean",
        choices=["ross", "pla"],
        default="ross",
        help="Opacity mean to plot (default: ross)",
    )
    parser.add_argument(
        "--log-rho-targets",
        type=float,
        nargs="+",
        # default=[-10, -8, -6, -4, -2],
        default=[-12, -11, -10, -9, -8],
        help="Target log10(rho[g/cm^3]) values for kappa(T) curves",
    )
    parser.add_argument(
        "--prefix",
        default="opacity_table",
        help="Prefix for output png names",
    )
    parser.add_argument(
        "--vmin",
        type=float,
        default=None,
        help="Fixed vmin for colormap (log10 kappa).",
    )
    parser.add_argument(
        "--vmax",
        type=float,
        default=None,
        help="Fixed vmax for colormap (log10 kappa).",
    )
    parser.add_argument(
        "--qmin",
        type=float,
        default=1.0,
        help="Lower percentile for auto color range (default: 1).",
    )
    parser.add_argument(
        "--qmax",
        type=float,
        default=99.0,
        help="Upper percentile for auto color range (default: 99).",
    )
    args = parser.parse_args()

    logt, logrho, logk, title = read_table(args.h5_file, args.mean)

    suffix = "ross" if args.mean == "ross" else "pla"
    out_colormap = f"{args.prefix}_{suffix}_colormap.png"
    out_curve = f"{args.prefix}_{suffix}_vsT.png"

    plot_colormap(
        logt,
        logrho,
        logk,
        args.mean,
        out_colormap,
        vmin=args.vmin,
        vmax=args.vmax,
        qmin=args.qmin,
        qmax=args.qmax,
    )
    plot_kappa_vs_t(
        logt,
        logrho,
        logk,
        args.log_rho_targets,
        out_curve,
        mean=args.mean,
    )

    print(f"Input: {args.h5_file}")
    print(f"Output: {out_colormap}")
    print(f"Output: {out_curve}")


if __name__ == "__main__":
    main()
