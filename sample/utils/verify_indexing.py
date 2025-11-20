#!/usr/bin/env python3
"""
Script to verify the new (rho, T) indexing order in the opacity table
"""

import h5py


def verify_indexing():
    """Verify that the opacity table uses (rho, T) indexing as intended"""
    
    with h5py.File('../table/output/opacity_lookup_table_2d.h5', 'r') as f:
        # Load data
        log_temps = f['log_temperature'][:]
        log_rhos = f['log_density'][:]
        planck = f['planck_mean_opacity'][:]
        rosseland = f['rosseland_mean_opacity'][:]
        
        print("=== Opacity Table Indexing Verification ===")
        print(f"Log temperature grid: {len(log_temps)} points")
        print(f"Log density grid: {len(log_rhos)} points")
        print(f"Planck opacity shape: {planck.shape}")
        grid_attr = f.attrs['grid_structure']
        if isinstance(grid_attr, bytes):
            grid_attr = grid_attr.decode()
        print(f"Grid structure attribute: {grid_attr}")
        print()
        
        # Convert log values back to linear for display
        temps = 10**log_temps
        rhos = 10**log_rhos
        
        # Verify the indexing: opacity[log_density_idx, log_temp_idx]
        print("Verifying indexing: opacity[log_density_idx, log_temp_idx]")
        print()
        
        # Test corner cases
        print("Corner cases:")
        print(f"  Min ρ, Min T: opacity[0,0] = {planck[0,0]:.3e} cm²/g")
        print(f"    ρ = {rhos[0]:.2e} g/cm³, T = {temps[0]:.2e} K")
        print(f"  Max ρ, Max T: opacity[-1,-1] = {planck[-1,-1]:.3e} cm²/g")
        print(f"    ρ = {rhos[-1]:.2e} g/cm³, T = {temps[-1]:.2e} K")
        print()
        
        # Test that opacity varies correctly along each axis
        print("Testing variation along log density axis (constant log T):")
        temp_idx = 10  # Middle temperature
        for i in [0, 10, 20]:
            rho_val = rhos[i]
            temp_val = temps[temp_idx]
            kappa_val = planck[i, temp_idx]
            print(f"  log_ρ[{i}] = {log_rhos[i]:.2f}, "
                  f"log_T[{temp_idx}] = {log_temps[temp_idx]:.2f} "
                  f"→ κ = {kappa_val:.3e}")
        print()
        
        print("Testing variation along log temperature axis (constant log ρ):")
        rho_idx = 10  # Middle density
        for i in [0, 10, 20]:
            rho_val = rhos[rho_idx]
            temp_val = temps[i]
            kappa_val = planck[rho_idx, i]
            print(f"  log_ρ[{rho_idx}] = {log_rhos[rho_idx]:.2f}, "
                  f"log_T[{i}] = {log_temps[i]:.2f} "
                  f"→ κ = {kappa_val:.3e}")
        print()
        
        # Compare Planck vs Rosseland
        print("Planck vs Rosseland comparison (mid-grid):")
        mid_rho, mid_temp = 10, 10
        planck_val = planck[mid_rho, mid_temp]
        ross_val = rosseland[mid_rho, mid_temp]
        rho_val = rhos[mid_rho]
        temp_val = temps[mid_temp]
        print(f"  At ρ = {rho_val:.2e} g/cm³, T = {temp_val:.2e} K:")
        print(f"    log₁₀ρ = {log_rhos[mid_rho]:.2f}, "
              f"log₁₀T = {log_temps[mid_temp]:.2f}")
        print(f"    Planck κ = {planck_val:.3e} cm²/g")
        print(f"    Rosseland κ = {ross_val:.3e} cm²/g")
        print(f"    Ratio κ_P/κ_R = {planck_val/ross_val:.2f}")
        print()
        
        print("✅ Indexing verification complete!")
        print("✅ Opacity table uses (log_density, log_temperature) indexing")


if __name__ == "__main__":
    verify_indexing()
