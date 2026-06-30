# %%
import os
import numpy as np
import matplotlib.pyplot as plt
import h5py
import argparse

# Constants
K_BOL = 1.3806503e-16

def read_mono_hdf5(file_path):
    try:
        with h5py.File(file_path, 'r') as f:
            return {name: f[name][0] for name in f.keys()}
    except Exception as e:
        raise RuntimeError(f"Failed to read {file_path}: {e}")


# %%
dir_path = '../table'
mean = 'ross'
# log_rho_targets = [-10, -8, -6, -4, -2]  # in log10(g/cm^3)
log_rho_targets = [-12, -11, -10, -9, -8]  # in log10(g/cm^3)
tol_dex = 0.2

files = [os.path.join(dir_path, 'output', f)
            for f in os.listdir(os.path.join(dir_path, 'output'))
            if f.startswith("mono_") and f.endswith(".h5")]
nmax = len(files)
if nmax == 0:
    raise ValueError("No HDF5 files found.")

print(f"Reading {nmax} layers ... ")

datasets = [read_mono_hdf5(file) for file in files]
tmp_tot = np.array([data['temp'] for data in datasets])         # T [K]
rho_tot = np.array([data['rho'] for data in datasets])          # rho [g/cm^3]
tmp2    = np.array([data['temp2'] for data in datasets])        # for pla2 title

# eps = 1e-30
# rho_tot = np.maximum(rho_tot, eps)


# log10(kappa / rho) [cm^2/g]
ros_tot  = np.log10(np.array([data['ros'] for data in datasets]) / rho_tot)
pla_tot  = np.log10(np.array([data['plac'] + data['plal']  for data in datasets]) / rho_tot)
pla2_tot = np.log10(np.array([data['plac2'] + data['plal2'] for data in datasets]) / rho_tot)

# Choose which opacity to plot
if mean == 'ross':
    logkappa = ros_tot
    title = 'Rosseland-mean opacity'
elif mean == 'pla':
    logkappa = pla_tot
    title = 'Planck-mean opacity'
elif mean == 'pla2':
    logkappa = pla2_tot
    title = f'Planck-mean opacity at'# T_rad={int(tmp2[0])} K'
    tmp = r'T_rad=' + f'{int(tmp2[0])} K'
    title += f' ({tmp})'
else:
    raise ValueError(f"Unknown mean type: {mean}")

# %%
# plt.rcParams.update({'font.size': 20})
fig, ax = plt.subplots(figsize=(8, 6))

log_rho_tot = np.log10(rho_tot)

for log_rho_target in log_rho_targets:
    # select points within tol_dex in log10(rho)
    mask = np.abs(log_rho_tot - log_rho_target) < tol_dex

    if not np.any(mask):
        print(f"[WARN] No points found near rho = {10**log_rho_target:.3e} g/cm^3")
        continue

    T_slice   = tmp_tot[mask]
    logk_slice = logkappa[mask]

    # sort by temperature
    idx = np.argsort(T_slice)
    T_slice   = T_slice[idx]
    logk_slice = logk_slice[idx]

    label = r'$\rho = $' + f'$10^{{{log_rho_target}}}$' + r' $\mathrm{g\, cm^{-3}}$'
    ax.plot(T_slice, logk_slice, marker='o', linestyle='-', label=label)


ax.tick_params(axis='both', which='major', labelsize=18)
ax.set_xscale('log')
ax.set_xlabel(r'$T\,\mathrm{[K]}$', fontsize=20)
ax.set_ylabel(r'$\log_{10}\,\kappa\,\mathrm{[cm^2\,g^{-1}]}$', fontsize=20)
ax.set_title(title, fontsize=20)
ax.grid(True)
ax.legend(fontsize=18)
plt.tight_layout()

plt.savefig('opac_vs_T_es_brems2.png', format='png')
plt.show()

# %%
print(10**(-0.462))