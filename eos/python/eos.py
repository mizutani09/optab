# %%
import h5py
import numpy as np
import matplotlib.pyplot as plt
import os
import subprocess
import matplotlib.colors as mcolors
from scipy.constants import physical_constants
import argparse

# Define the atomic mass unit in grams using scipy's physical constants
amu = physical_constants['atomic mass constant'][0] * 1000  # Convert kg to g

def read_dataset(file, name):
    """
    Read and return the dataset from an HDF5 file.
    """
    with h5py.File(file, 'r') as h5file:
        if name == 'n_layer':
            for key in h5file.keys():
                print(f"{key}, Shape: {h5file[key].shape}, Dtype: {h5file[key].dtype}")
        return h5file[name][()]

def plot_data(fname, variable, syms):
    """
    Plot the data based on the variable type ('rho' or 'mmw') using a color map.
    """
    # Read the datasets
    n_layer = read_dataset(fname, 'n_layer')
    temp = read_dataset(fname, 'temp')
    pres = read_dataset(fname, 'pres')
    rho = read_dataset(fname, 'rho')
    ndens = read_dataset(fname, 'ndens')
    
    print("temp shape:", temp.shape)
    print("pres shape:", pres.shape)
    print("rho shape:", rho.shape)
    print("ndens shape:", ndens.shape)
    
    # count identical temperatures
    unique_temps = np.unique(temp)
    print(f"Number of unique temperatures: {len(unique_temps)}")
    unique_pressures = np.unique(pres)
    print(f"Number of unique pressures: {len(unique_pressures)}")
    unique_densities = np.unique(rho)
    print(f"Number of unique densities: {len(unique_densities)}")
    # return
    
    # Calculate the variable to be plotted and its range
    if variable == 'rho':
        title = 'Mass Density'
        cb_title = 'log ρ [g cm⁻³]'
        v = np.log10(rho)
        v_min, v_max = -16, -2
    elif variable == 'mmw':
        title = 'Mean Molecular Weight'
        cb_title = 'Mean Molecular Weight'
        v = rho / ndens[:, 0] / amu
        v_min, v_max = np.min(v), np.max(v)

    # Set up the normalization and color map
    norm = mcolors.Normalize(vmin=v_min, vmax=v_max)
    cmap = plt.cm.jet

    # Create the plot
    fig, ax = plt.subplots()
    sc = ax.scatter(temp, pres, c=v, s=syms, marker='s', cmap=cmap, norm=norm)

    # Add titles and labels
    plt.title(title)
    cbar = plt.colorbar(sc, label=cb_title)
    cbar.set_label(cb_title, rotation=270, labelpad=15)

    # Set the scales and labels of the axes
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Temperature [K]')
    ax.set_ylabel('Pressure [Ba]')    
    
    # Show the plot
    png_filename = fname.replace('.h5', '.png')
    plt.savefig(png_filename, format='png', transparent=True, dpi=600)
    plt.show()

def main():
    """
    Main function to execute the script.
    """
    # plot_data('../FastChem/withline_101_101.h5', 'rho', syms=100)
    # return
    # Define the filename and variable to plot
    parser = argparse.ArgumentParser(description='Plot data from an HDF5 file.')
    parser.add_argument('filename', type=str, help='Name of the HDF5 file.')
    parser.add_argument('variable', type=str, choices=['rho', 'mmw'], help='The variable to plot (rho or mmw).')
    parser.add_argument('--syms', type=int, default=5, help='Symbol size for the scatter plot (default: 5).')
    args = parser.parse_args()
    
    # Call the plot function
    plot_data(args.filename, args.variable, args.syms)

# Execute the main function
if __name__ == "__main__":
    main()
