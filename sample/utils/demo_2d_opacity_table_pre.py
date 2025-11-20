#!/usr/bin/env python3
"""
Script to demonstrate how to use the 2D opacity lookup table
"""
# %%
import h5py
import numpy as np
import matplotlib.pyplot as plt


def load_2d_opacity_table(filename):
    """Load the 2D opacity lookup table from HDF5 file"""
    with h5py.File(filename, 'r') as f:
        data = {
            'log_temperature': f['log_temperature'][:],
            'log_pressure': f['log_pressure'][:],
            'planck_mean_opacity': f['planck_mean_opacity'][:],  # 2D
            'rosseland_mean_opacity': f['rosseland_mean_opacity'][:]  # 2D
        }

        # Print structure information
        print("2D Opacity table loaded:")
        log_temp_len = len(data['log_temperature'])
        log_press_len = len(data['log_pressure'])
        print(f"  Log temperature grid: {data['log_temperature'].shape} - "
              f"{log_temp_len} points")
        print(f"  Log pressure grid: {data['log_pressure'].shape} - "
              f"{log_press_len} points")
        print(f"  Planck opacity: {data['planck_mean_opacity'].shape} - "
              f"2D grid")
        print(f"  Rosseland opacity: {data['rosseland_mean_opacity'].shape} - "
              f"2D grid")

        # Get grid structure attribute
        grid_attr = f.attrs.get('grid_structure', '')
        if isinstance(grid_attr, bytes):
            grid_structure = grid_attr.decode('utf-8')
        else:
            grid_structure = str(grid_attr)
        print(f"  Grid structure: {grid_structure}")

    return data


def interpolate_2d_opacity(data, target_temp, target_press):
    """
    Interpolate opacity at given temperature and pressure using 2D grid
    This uses logarithmic interpolation for better accuracy
    """

    # Convert targets to log space
    log_target_temp = np.log10(target_temp)
    log_target_press = np.log10(target_press)

    # Find surrounding grid points
    log_temp_grid = data['log_temperature']
    log_press_grid = data['log_pressure']

    # Find temperature indices
    if log_target_temp <= log_temp_grid[0]:
        i1, i2 = 0, 0
        t_frac = 0.0
    elif log_target_temp >= log_temp_grid[-1]:
        i1, i2 = len(log_temp_grid)-1, len(log_temp_grid)-1
        t_frac = 0.0
    else:
        i2 = np.searchsorted(log_temp_grid, log_target_temp)
        i1 = i2 - 1
        temp_diff = log_temp_grid[i2] - log_temp_grid[i1]
        t_frac = (log_target_temp - log_temp_grid[i1]) / temp_diff
        print("i1, i2 = ", i1, i2)
        print("t_frac = ", t_frac)

    # Find pressure indices
    if log_target_press <= log_press_grid[0]:
        j1, j2 = 0, 0
        r_frac = 0.0
    elif log_target_press >= log_press_grid[-1]:
        j1, j2 = len(log_press_grid)-1, len(log_press_grid)-1
        r_frac = 0.0
    else:
        j2 = np.searchsorted(log_press_grid, log_target_press)
        j1 = j2 - 1
        press_diff = log_press_grid[j2] - log_press_grid[j1]
        r_frac = (log_target_press - log_press_grid[j1]) / press_diff
        print("j1, j2 = ", j1, j2)
        print("r_frac = ", r_frac)

    # Bilinear interpolation on 2D grid: opacity[log_pressure_idx, log_temp_idx]
    planck_grid = data['planck_mean_opacity']
    ross_grid = data['rosseland_mean_opacity']

    # Get corner values - NOTE: grid is [log_pressure, log_temperature]
    p11 = planck_grid[j1, i1]
    p12 = planck_grid[j2, i1]
    p21 = planck_grid[j1, i2]
    p22 = planck_grid[j2, i2]

    r11 = ross_grid[j1, i1]
    r12 = ross_grid[j2, i1]
    r21 = ross_grid[j1, i2]
    r22 = ross_grid[j2, i2]

    # Interpolate
    # planck_interp = ((1-t_frac)*(1-r_frac)*p11 + (1-t_frac)*r_frac*p12 +
    #                  t_frac*(1-r_frac)*p21 + t_frac*r_frac*p22)
    # ross_interp = ((1-t_frac)*(1-r_frac)*r11 + (1-t_frac)*r_frac*r12 +
    #                t_frac*(1-r_frac)*r21 + t_frac*r_frac*r22)
    
    x1norm = (len(log_temp_grid) - 1) / (log_temp_grid[-1] - log_temp_grid[0])
    x2norm = (len(log_press_grid) - 1) / (log_press_grid[-1] - log_press_grid[0])
    
    t = (log_target_temp - log_temp_grid[0]) * x1norm
    r = (log_target_press - log_press_grid[0]) * x2norm
    
    til = int(np.floor(t))
    trl = 1 + til - t
    ril = int(np.floor(r))
    rrl = 1 + ril - r
    
    t_frac = trl
    r_frac = rrl
    print("t_frac = ", t_frac)
    print("r_frac = ", r_frac)
    
    planck_interp = t_frac*r_frac*planck_grid[j1, i1] \
        + (1-t_frac)*r_frac*planck_grid[j1, i2] \
        + t_frac*(1-r_frac)*planck_grid[j2, i1] \
        + (1-t_frac)*(1-r_frac)*planck_grid[j2, i2]

    ross_interp = t_frac*r_frac*ross_grid[j1, i1] \
        + (1-t_frac)*r_frac*ross_grid[j1, i2] \
        + t_frac*(1-r_frac)*ross_grid[j2, i1] \
        + (1-t_frac)*(1-r_frac)*ross_grid[j2, i2]
    return planck_interp, ross_interp


def demonstrate_2d_table(data_folder, output_folder):
    """Demonstrate usage of the 2D opacity lookup table"""

    # Load the 2D lookup table
    table_file = f"{data_folder}/opacity_lookup_table_2d.h5"
    print(f"Loading 2D opacity lookup table from: {table_file}")

    try:
        data = load_2d_opacity_table(table_file)
    except FileNotFoundError:
        print(f"Error: Could not find {table_file}")
        print("Please run the create_opacity_table_2d program first.")
        return

    # Print basic statistics
    print("\nGrid information:")

    # Convert log values back to linear for display
    temp_min = 10**data['log_temperature'].min()
    temp_max = 10**data['log_temperature'].max()
    press_min = 10**data['log_pressure'].min()
    press_max = 10**data['log_pressure'].max()

    print(f"Temperature range: {temp_min:.2e} - {temp_max:.2e} K")
    print(f"Pressure range: {press_min:.2e} - {press_max:.2e} g/cm³")

    planck_min = data['planck_mean_opacity'].min()
    planck_max = data['planck_mean_opacity'].max()
    print(f"Planck opacity range: {planck_min:.2e} - {planck_max:.2e} cm²/g")

    ross_min = data['rosseland_mean_opacity'].min()
    ross_max = data['rosseland_mean_opacity'].max()
    print(f"Rosseland opacity range: {ross_min:.2e} - "
          f"{ross_max:.2e} cm²/g")

    # # Demonstrate fast 2D interpolation in log space
    # print("\n=== 2D Interpolation Examples (Log Space) ===")

    # test_points = [
    #     (1e4, 1e-12),   # Low temperature, low pressure
    #     (1e6, 1e-6),    # High temperature, high pressure
    #     (5623, 7e-11),  # Mid-range values
    # ]

    # for temp, press in test_points:
    #     planck, ross = interpolate_2d_opacity(data, temp, press)
    #     print(f"T = {temp:.0e} K, ρ = {press:.0e} g/cm³:")
    #     print(f"  Planck opacity = {planck:.6e} cm²/g")
    #     print(f"  Rosseland opacity = {ross:.6e} cm²/g")

    # Create visualization of the 2D grid
    # print("\nCreating 2D grid visualization...")

    # fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))
    # title_fontsize = 22
    # label_fontsize = 20
    # ticks_fontsize = 15

    # # Convert log values back to linear for display
    # temp_linear = 10**data['log_temperature']
    # press_linear = 10**data['log_pressure']

    # # Temperature vs Pressure grid points
    # T_grid, R_grid = np.meshgrid(temp_linear, press_linear,
    #                              indexing='ij')

    # # 1. Temperature-Pressure grid
    # ax1.scatter(T_grid.flatten(), R_grid.flatten(), s=20, alpha=0.7)
    # ax1.set_xscale('log')
    # ax1.set_yscale('log')
    # ax1.set_xlabel('Temperature [K]', fontsize=label_fontsize)
    # ax1.set_ylabel('Pressure [g/cm³]', fontsize=label_fontsize)
    # ax1.set_title('2D Grid Points', fontsize=title_fontsize)
    # ax1.tick_params(axis='both', which='major', labelsize=ticks_fontsize)
    # ax1.tick_params(axis='both', which='minor', labelsize=ticks_fontsize)
    # ax1.grid(True, alpha=0.3)

    # # 2. Planck opacity surface (transpose for display)
    # planck_display = data['planck_mean_opacity'].T
    # im1 = ax2.imshow(np.log10(planck_display),
    #                  extent=[data['log_pressure'].min(),
    #                          data['log_pressure'].max(),
    #                          data['log_temperature'].min(),
    #                          data['log_temperature'].max()],
    #                  aspect='auto', origin='lower', cmap='viridis')
    # ax2.set_xlabel(
    #     r'$\log_{10}\press\,[\mathrm{g}\,\mathrm{cm}^{-3}]$',
    #     fontsize=label_fontsize)
    # ax2.set_ylabel(
    #     r'$\log_{10}T\,[\mathrm{K}]$', fontsize=label_fontsize)
    # ax2.set_title('log₁₀(Planck Opacity)', fontsize=title_fontsize)
    # ax2.tick_params(axis='both', which='major', labelsize=ticks_fontsize)
    # ax2.tick_params(axis='both', which='minor', labelsize=ticks_fontsize)
    # im1.set_clim(vmin=np.log10(planck_min), vmax=np.log10(planck_max))
    # cb = ax2.figure.colorbar(im1, ax=ax2, orientation='vertical')
    # planck_label = r'$\log_{10}\kappa_{\mathrm{P}}\,[\mathrm{cm}^2\,\mathrm{g}^{-1}]$'
    # cb.set_label(planck_label, fontsize=label_fontsize)

    # # 3. Rosseland opacity surface (transpose for correct display)
    # ross_display = data['rosseland_mean_opacity'].T
    # im2 = ax3.imshow(np.log10(ross_display),
    #                  extent=[data['log_pressure'].min(),
    #                          data['log_pressure'].max(),
    #                          data['log_temperature'].min(),
    #                          data['log_temperature'].max()],
    #                  aspect='auto', origin='lower', cmap='plasma')
    # ax3.set_xlabel(
    #     r'$\log_{10}\press\,[\mathrm{g}\,\mathrm{cm}^{-3}]$',
    #     fontsize=label_fontsize)
    # ax3.set_ylabel(
    #     r'$\log_{10}T\,[\mathrm{K}]$', fontsize=label_fontsize)
    # ax3.set_title('log₁₀(Rosseland Opacity)', fontsize=title_fontsize)
    # ax3.tick_params(axis='both', which='major', labelsize=ticks_fontsize)
    # ax3.tick_params(axis='both', which='minor', labelsize=ticks_fontsize)
    # im2.set_clim(vmin=np.log10(ross_min), vmax=np.log10(ross_max))
    # cb = ax3.figure.colorbar(im2, ax=ax3, orientation='vertical')
    # ross_label = r'$\log_{10}\kappa_{\mathrm{R}}\,[\mathrm{cm}^2\,\mathrm{g}^{-1}]$'
    # cb.set_label(ross_label, fontsize=label_fontsize)

    # # 4. Opacity comparison at different temperatures (mid pressure slice)
    # mid_pressure_idx = len(data['log_pressure']) // 2
    # mid_press = 10**data["log_pressure"][mid_pressure_idx]
    # # planck_label = r'Planck ($\press$={{:.0e}} $\mathrm{g}\,\mathrm{cm}^{-3}$)'.format(mid_press)
    # # ross_label = r'Rosseland ($\press$={:.0e} $\mathrm{g}\,\mathrm{cm}^{-3}$)'.format(mid_press)
    # planck_label = f'Planck ($\\press$={mid_press:.0e} g/cm$^3$)'
    # ross_label = f'Rosseland ($\\press$={mid_press:.0e} g/cm$^3$)'
    # ax4.loglog(temp_linear,
    #            data['planck_mean_opacity'][mid_pressure_idx, :],
    #            'o-',
    #            label=planck_label,
    #            markersize=4)
    # ax4.loglog(temp_linear,
    #            data['rosseland_mean_opacity'][mid_pressure_idx, :],
    #            's-',
    #            label=ross_label,
    #            markersize=4)
    # # ax4.set_xlabel('Temperature [K]')
    # # ax4.set_ylabel('Opacity [cm²/g]')
    # # ax4.set_title('Opacity vs Temperature (Mid Pressure)')
    # ax4.set_xlabel('Temperature [K]', fontsize=label_fontsize)
    # ax4.set_ylabel('Opacity [cm²/g]', fontsize=label_fontsize)
    # # ax4.set_title('Opacity vs Temperature (Mid Pressure)', fontsize=title_fontsize)
    # ax4.tick_params(axis='both', which='major', labelsize=ticks_fontsize)
    # ax4.tick_params(axis='both', which='minor', labelsize=ticks_fontsize)
    # ax4.legend(fontsize=label_fontsize-6)
    # ax4.grid(True, alpha=0.3)

    # plt.tight_layout()

    # # Save the plot
    # output_plot = f"{output_folder}/opacity_lookup_2d_visualization.png"
    # plt.savefig(output_plot, dpi=150, bbox_inches='tight')
    # print(f"2D grid visualization saved as: {output_plot}")

    # plt.show()
    
    tmp_rho = 1e-10
    tmp_tempature = 1e5
    tmp_pressure = 1.662902000000000e+03
    planck_interp, ross_interp = interpolate_2d_opacity(data, tmp_tempature, tmp_pressure)
    print(f"\nInterpolated Planck opacity at T={tmp_tempature:.0e} K, "
          f"press={tmp_pressure:.0e} g/cm³: {planck_interp:.6e} cm²/g")
    print(f"Interpolated Rosseland opacity at T={tmp_tempature:.0e} K, "
          f"press={tmp_pressure:.0e} g/cm³: {ross_interp:.6e} cm²/g")


if __name__ == "__main__":
    print("=== 2D Opacity Lookup Table Demonstration ===")
    data_folder = "../table/output"
    output_folder = data_folder  # Assuming output is in the same folder
    demonstrate_2d_table(data_folder, output_folder)
    print("\nDemonstration completed!")
