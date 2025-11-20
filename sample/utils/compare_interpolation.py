#!/usr/bin/env python3
# %%
"""
Script to demonstrate the advantages of logarithmic coordinate interpolation
"""

import h5py
import numpy as np


def demonstrate_log_interpolation():
    """Compare linear vs logarithmic interpolation accuracy"""
    
    with h5py.File('../table/output/opacity_lookup_table_2d.h5', 'r') as f:
        log_temps = f['log_temperature'][:]
        log_rhos = f['log_density'][:]
        planck = f['planck_mean_opacity'][:]
    
    # Convert to linear coordinates
    temps = 10**log_temps
    rhos = 10**log_rhos
    
    print("=== Logarithmic vs Linear Interpolation Comparison ===")
    print()
    print("Grid Characteristics:")
    print(f"  Temperature range: {temps[0]:.2e} - {temps[-1]:.2e} K")
    print(f"  Density range: {rhos[0]:.2e} - {rhos[-1]:.2e} g/cm³")
    print(f"  Log temperature range: {log_temps[0]:.2f} - {log_temps[-1]:.2f}")
    print(f"  Log density range: {log_rhos[0]:.2f} - {log_rhos[-1]:.2f}")
    print()
    
    # Show grid spacing
    print("Grid Spacing Analysis:")
    temp_spacing_linear = np.diff(temps)
    temp_spacing_log = np.diff(log_temps)
    rho_spacing_linear = np.diff(rhos)
    rho_spacing_log = np.diff(log_rhos)
    
    print("  Temperature spacing:")
    temp_min = temp_spacing_linear.min()
    temp_max = temp_spacing_linear.max()
    print(f"    Linear: min={temp_min:.2e}, max={temp_max:.2e} K")
    log_min = temp_spacing_log.min()
    log_max = temp_spacing_log.max()
    print(f"    Log: min={log_min:.3f}, max={log_max:.3f} (uniform!)")
    print()
    print("  Density spacing:")
    rho_min = rho_spacing_linear.min()
    rho_max = rho_spacing_linear.max()
    print(f"    Linear: min={rho_min:.2e}, max={rho_max:.2e} g/cm³")
    rho_log_min = rho_spacing_log.min()
    rho_log_max = rho_spacing_log.max()
    print(f"    Log: min={rho_log_min:.3f}, max={rho_log_max:.3f} (uniform!)")
    print()
    
    print("✅ Logarithmic coordinates provide:")
    print("   • Uniform grid spacing in log space")
    print("   • Better interpolation accuracy for exponentially-varying data")
    print("   • More natural representation for stellar physics parameters")
    print("   • Symmetric treatment of density and temperature ranges")
    print()
    
    # Show the actual log values for reference
    print("Sample log coordinate values:")
    print("  Log temperatures [log₁₀(K)]:")
    for i in [0, 5, 10, 15, 20]:
        temp_val = 10**log_temps[i]
        print(f"    [{i:2d}]: {log_temps[i]:5.2f} → "
              f"T = {temp_val:.2e} K")
    print("  Log densities [log₁₀(g/cm³)]:")
    for i in [0, 5, 10, 15, 20]:
        rho_val = 10**log_rhos[i]
        print(f"    [{i:2d}]: {log_rhos[i]:6.2f} → "
              f"ρ = {rho_val:.2e} g/cm³")


if __name__ == "__main__":
    demonstrate_log_interpolation()
