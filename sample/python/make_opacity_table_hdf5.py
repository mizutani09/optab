import argparse
import glob
import os

import h5py
import numpy as np


def read_layer(path):
    with h5py.File(path, "r") as f:
        temp = float(f["temp"][0])
        rho = float(f["rho"][0])
        ros = float(f["ros"][0]) / rho
        pla = (float(f["plac"][0]) + float(f["plal"][0])) / rho
    return temp, rho, ros, pla


def build_points(files):
    temp = []
    rho = []
    ros = []
    pla = []
    for file_path in files:
        t, r, kr, kp = read_layer(file_path)
        temp.append(t)
        rho.append(r)
        ros.append(kr)
        pla.append(kp)

    temp = np.asarray(temp, dtype=np.float64)
    rho = np.asarray(rho, dtype=np.float64)
    ros = np.asarray(ros, dtype=np.float64)
    pla = np.asarray(pla, dtype=np.float64)

    # Guard against log10(0). Very small values are still represented.
    ros = np.clip(ros, 1e-300, None)
    pla = np.clip(pla, 1e-300, None)
    return temp, rho, ros, pla


def compute_full_rectangle(logt_pts, logr_pts, round_decimals=10):
    tkey = np.round(logt_pts, round_decimals)
    unique_t = np.unique(tkey)
    rho_min_each_t = []
    rho_max_each_t = []
    for t in unique_t:
        m = tkey == t
        rho_min_each_t.append(np.min(logr_pts[m]))
        rho_max_each_t.append(np.max(logr_pts[m]))

    # Rectangle fully covered across all temperatures:
    # logT in [min(unique_t), max(unique_t)],
    # logRho in [max(min_rho(T)), min(max_rho(T))]
    logt_min = float(np.min(unique_t))
    logt_max = float(np.max(unique_t))
    logrho_min = float(np.max(rho_min_each_t))
    logrho_max = float(np.min(rho_max_each_t))
    if not np.isfinite(logrho_min) or not np.isfinite(logrho_max) or logrho_min >= logrho_max:
        raise RuntimeError("Failed to find a non-empty fully-contained (logT, logRho) rectangle.")
    return logt_min, logt_max, logrho_min, logrho_max


def idw_interpolate(logt_pts, logr_pts, values, logt_grid, logr_grid, power=2.0, eps=1e-12):
    tt, rr = np.meshgrid(logt_grid, logr_grid, indexing="xy")
    out = np.empty_like(tt)

    for j in range(tt.shape[0]):
        dt = logt_pts[:, None] - tt[j, :][None, :]
        dr = logr_pts[:, None] - rr[j, :][None, :]
        dist2 = dt * dt + dr * dr

        exact = dist2 < eps
        row = np.empty(tt.shape[1], dtype=np.float64)

        for i in range(tt.shape[1]):
            if np.any(exact[:, i]):
                row[i] = values[np.argmax(exact[:, i])]
                continue
            w = 1.0 / np.power(dist2[:, i], 0.5 * power)
            row[i] = np.sum(w * values) / np.sum(w)

        out[j, :] = row

    return out


def structured_interpolate(logt_pts, logr_pts, values, logt_grid, logr_grid, round_decimals=10):
    tkey = np.round(logt_pts, round_decimals)
    unique_t = np.unique(tkey)

    # First pass: for each source T, interpolate along rho -> common rho grid.
    vrho = np.empty((len(unique_t), len(logr_grid)), dtype=np.float64)
    for i, t in enumerate(unique_t):
        m = tkey == t
        x = logr_pts[m]
        y = values[m]
        order = np.argsort(x)
        x = x[order]
        y = y[order]
        # Remove duplicate rho samples if any.
        xu, idx = np.unique(x, return_index=True)
        yu = y[idx]
        vrho[i, :] = np.interp(logr_grid, xu, yu)

    # Second pass: for each rho on target grid, interpolate in T.
    out = np.empty((len(logr_grid), len(logt_grid)), dtype=np.float64)
    for j in range(len(logr_grid)):
        out[j, :] = np.interp(logt_grid, unique_t, vrho[:, j])
    return out


def main():
    parser = argparse.ArgumentParser(
        description="Pack OPTAB Rosseland/Planck means into a single HDF5 table."
    )
    parser.add_argument(
        "dir_path",
        help="Path containing output/mono_*.h5 (e.g., ../table_detail).",
    )
    parser.add_argument(
        "--out",
        default="opacity_table_athena.h5",
        help="Output HDF5 path. Default: opacity_table_athena.h5",
    )
    parser.add_argument(
        "--nlogt",
        type=int,
        default=128,
        help="Number of log10(T) grid points. Default: 128",
    )
    parser.add_argument(
        "--nlogrho",
        type=int,
        default=128,
        help="Number of log10(rho) grid points. Default: 128",
    )
    parser.add_argument(
        "--power",
        type=float,
        default=2.0,
        help="IDW interpolation power. Default: 2.0",
    )
    parser.add_argument(
        "--interp",
        choices=["structured", "idw"],
        default="structured",
        help=(
            "Interpolation method. "
            "'structured' uses 1D(rho)+1D(T) interpolation on source T-slices; "
            "'idw' uses inverse-distance weighting in 2D."
        ),
    )
    parser.add_argument(
        "--crop",
        choices=["full", "none"],
        default="full",
        help=(
            "Grid range selection in (logT,logRho). "
            "'full': use fully-contained rectangle. "
            "'none': use full min-max bounding box."
        ),
    )
    args = parser.parse_args()

    pattern = os.path.join(args.dir_path, "output", "mono_*.h5")
    files = sorted(glob.glob(pattern))
    if not files:
        raise FileNotFoundError(f"No files matched: {pattern}")

    temp, rho, ros, pla = build_points(files)
    logt_pts = np.log10(temp)
    logr_pts = np.log10(rho)
    logros_pts = np.log10(ros)
    logpla_pts = np.log10(pla)

    if args.crop == "full":
        logt_min, logt_max, logrho_min, logrho_max = compute_full_rectangle(logt_pts, logr_pts)
    else:
        logt_min = float(np.min(logt_pts))
        logt_max = float(np.max(logt_pts))
        logrho_min = float(np.min(logr_pts))
        logrho_max = float(np.max(logr_pts))

    logt_grid = np.linspace(logt_min, logt_max, args.nlogt)
    logr_grid = np.linspace(logrho_min, logrho_max, args.nlogrho)

    if args.interp == "structured":
        logros_table = structured_interpolate(logt_pts, logr_pts, logros_pts, logt_grid, logr_grid)
        logpla_table = structured_interpolate(logt_pts, logr_pts, logpla_pts, logt_grid, logr_grid)
        interp_desc = "structured(1D-rho + 1D-T) in (log10T, log10rho)"
    else:
        logros_table = idw_interpolate(logt_pts, logr_pts, logros_pts, logt_grid, logr_grid, power=args.power)
        logpla_table = idw_interpolate(logt_pts, logr_pts, logpla_pts, logt_grid, logr_grid, power=args.power)
        interp_desc = f"IDW(power={args.power}) in (log10T, log10rho)"

    # Table shape follows [n_logrho, n_logt].
    with h5py.File(args.out, "w") as f:
        axes = f.create_group("axes")
        axes.create_dataset("log10_temperature", data=logt_grid)
        axes.create_dataset("log10_density", data=logr_grid)
        axes["log10_temperature"].attrs["unit"] = "K"
        axes["log10_density"].attrs["unit"] = "g cm^-3"

        kappa = f.create_group("kappa")
        kappa.create_dataset("log10_rosseland", data=logros_table)
        kappa.create_dataset("log10_planck", data=logpla_table)
        kappa.create_dataset("rosseland", data=np.power(10.0, logros_table))
        kappa.create_dataset("planck", data=np.power(10.0, logpla_table))
        kappa["log10_rosseland"].attrs["unit"] = "log10(cm^2 g^-1)"
        kappa["log10_planck"].attrs["unit"] = "log10(cm^2 g^-1)"
        kappa["rosseland"].attrs["unit"] = "cm^2 g^-1"
        kappa["planck"].attrs["unit"] = "cm^2 g^-1"
        kappa.attrs["shape"] = "[n_logrho, n_logt]"

        points = f.create_group("source_points")
        points.create_dataset("temperature", data=temp)
        points.create_dataset("density", data=rho)
        points.create_dataset("log10_rosseland", data=logros_pts)
        points.create_dataset("log10_planck", data=logpla_pts)

        f.attrs["generator"] = "sample/python/make_opacity_table_hdf5.py"
        f.attrs["input_pattern"] = pattern
        f.attrs["interpolation"] = interp_desc
        f.attrs["crop_mode"] = args.crop
        f.attrs["log10T_min"] = logt_min
        f.attrs["log10T_max"] = logt_max
        f.attrs["log10Rho_min"] = logrho_min
        f.attrs["log10Rho_max"] = logrho_max
        f.attrs["note"] = (
            "Built from scattered OPTAB layers. "
            "Use kappa/log10_* datasets for Athena++ table conversion."
        )

    print(f"Read {len(files)} files.")
    print(f"Wrote: {args.out}")
    print(f"Grid: n_logrho={args.nlogrho}, n_logt={args.nlogt}")
    print(f"crop={args.crop}: log10T=[{logt_min:.3f}, {logt_max:.3f}], log10Rho=[{logrho_min:.3f}, {logrho_max:.3f}]")


if __name__ == "__main__":
    main()
