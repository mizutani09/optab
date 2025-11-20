# %%
import os
import numpy as np
import matplotlib.pyplot as plt
import h5py
# import argparse
from scipy.interpolate import griddata

# Constants
K_BOL = 1.3806503e-16


def read_mono_hdf5(file_path):
    try:
        with h5py.File(file_path, 'r') as f:
            return {name: f[name][0] for name in f.keys()}
    except Exception as e:
        raise RuntimeError(f"Failed to read {file_path}: {e}")


# def opac(dir_path, mean, syms):

# opac('../table', 'ross', 150)
dir_path = '../table'
mean = 'ross'
syms = 150
files = [os.path.join(dir_path, 'output', f) for f in os.listdir(os.path.join(dir_path, 'output')) if f.startswith("mono_") and f.endswith(".h5")]
nmax = len(files)
files.sort()

if nmax == 0:
    raise ValueError("No HDF5 files found.")

print(f"Reading {nmax} layers ... ")

datasets = [read_mono_hdf5(file) for file in files]
tmp2 = np.array([data['temp2'] for data in datasets])
tmp_tot = np.array([data['temp'] for data in datasets])
rho_tot = np.array([data['rho'] for data in datasets])
nden_tot = np.array([data['nden'] for data in datasets])
pre_tot = np.array([data['nden'] * K_BOL * data['temp'] for data in datasets])
ros_tot = np.log10(np.array([data['ros'] for data in datasets]) / rho_tot)
pla_tot = np.log10(np.array([data['plac'] + data['plal'] for data in datasets]) / rho_tot)
pla2_tot = np.log10(np.array([data['plac2'] + data['plal2'] for data in datasets]) / rho_tot)
data_dict = {
    'rho': rho_tot,
    'temp': tmp_tot,
    'pre': pre_tot,
}

# Title and label mapping
titles = {
    'ross': ('Rosseland-mean opacity', 'log $\kappa$ [cm$^2$/g]'),
    'pla': ('Planck-mean opacity', 'log $\kappa$ [cm$^2$/g]'),
    'pla2': (f'Planck-mean opacity at T_rad={int(tmp2[0])}K', 'log $\kappa$ [cm$^2$/g]')
}

if mean not in titles:
    raise ValueError(f"Unknown mean type: {mean}")
title, btitle = titles[mean]

# first make 2D coord for rho and temp
temp_coord = np.logspace(np.log10(np.min(tmp_tot)), np.log10(np.max(tmp_tot)), num=100)

rho_min = 1e-12
rho_max = 1e-8
rho_coord = np.logspace(np.log10(rho_min), np.log10(rho_max), num=100)

# create meshgrid
T_grid, R_grid = np.meshgrid(temp_coord, rho_coord)

# interpolate ros_tot onto the new grid

Z_grid = griddata((tmp_tot, rho_tot), ros_tot, (T_grid, R_grid), method='linear')
# print(Z_grid.shape)

# # make plot
# plt.figure(figsize=(10, 8))
# plt.contourf(T_grid, R_grid, Z_grid, levels=100, cmap='plasma')
# plt.yscale('log')
# plt.xscale('log')
# plt.xlabel('T [K]')
# plt.ylabel(r'$\rho$ [g/cm$^3$]')
# plt.title(f'Interpolated {title}')
# cbar = plt.colorbar(label=btitle)
# cbar.set_label(btitle, rotation=270, labelpad=15)
# plt.grid(True)
# plt.savefig(f'{file_basename}_interpolated.png', format='png')
# plt.show()

# make plot for each density
density_values = [1e-8, 1e-9, 1e-10, 1e-11, 1e-12]
fig, ax = plt.subplots(figsize=(10, 8))
for density in density_values:
    idx = np.argmin(np.abs(rho_coord - density))
    opacity_slice = np.power(10, Z_grid[idx, :])
    ax.plot(temp_coord, opacity_slice,
            label=r'$\rho$ = {:.1e} $\mathrm{{g}}\,\mathrm{{cm}}^{{-3}}$'.format(density))
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$T$ [$\mathrm{K}$]', fontsize=25)
    ax.set_ylabel(r'$\kappa_{\mathrm{R}}$ [$\mathrm{cm}^2\,\mathrm{g}^{-1}$]', fontsize=25)
    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.set_ylim(1e-5, 1e4)
    ax.set_xlim(8e2, 3e5)
    ax.grid(True)

ax.legend(fontsize=16)
figname = os.path.join(dir_path, 'output', f'{mean}_slices.png')
plt.tight_layout()
plt.savefig(figname, format='png', dpi=300)
plt.show()
