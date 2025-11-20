# Opacity Table Quick Reference

## File Format: `opacity_lookup_table_2d.h5` ⭐ RECOMMENDED
```
/log_temperature[21]         # log10(K) - Log temperature grid
/log_density[21]             # log10(g/cm³) - Log density grid  
/planck_mean_opacity[21,21]  # cm²/g - opacity[log_rho_idx, log_temp_idx]
/rosseland_mean_opacity[21,21] # cm²/g - opacity[log_rho_idx, log_temp_idx]
```

## Data Ranges
- **Log Temperature**: log₁₀(3.16×10²) to log₁₀(3.16×10⁷) (21 points, linear in log space)
- **Log Density**: log₁₀(8.89×10⁻¹⁵) to log₁₀(8.89×10⁻²) g/cm³ (21 points per temperature)
- **Total data points**: 441 (21×21 grid)

## Quick Usage

### Python
```python
import h5py
with h5py.File('opacity_lookup_table_2d.h5', 'r') as f:
    log_T_grid = f['log_temperature'][:]      # [21] log temperature points
    log_rho_grid = f['log_density'][:]        # [21] log density points
    kappa_P = f['planck_mean_opacity'][:]     # [21,21] Planck opacity
    kappa_R = f['rosseland_mean_opacity'][:] # [21,21] Rosseland opacity

# Access: kappa_P[j,i] = opacity at 10^log_rho_grid[j], 10^log_T_grid[i]
```

### C/C++
```cpp
#include <hdf5.h>
double log_temperature[21], log_density[21];
double planck_opacity[21][21], rosseland_opacity[21][21];

hid_t file = H5Fopen("opacity_lookup_table_2d.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
// Read datasets...
// Access: planck_opacity[log_density_idx][log_temp_idx]
```

### Fortran
```fortran
real(8) :: log_temperature(21), log_density(21)
real(8) :: planck_opacity(21,21), rosseland_opacity(21,21)
! Read from HDF5...
! Access: planck_opacity(j,i) for log_density(j), log_temperature(i)
```

## Interpolation Template (Python)
```python
def interpolate_opacity(target_T, target_rho):
    # Convert to log space
    log_target_T = np.log10(target_T)
    log_target_rho = np.log10(target_rho)
    
    # Find grid positions
    i = np.searchsorted(log_T_grid, log_target_T)
    j = np.searchsorted(log_rho_grid, log_target_rho)
    
    # Bilinear interpolation in log space (opacity[log_rho_idx, log_temp_idx])
    t_frac = (log_target_T - log_T_grid[i-1]) / (log_T_grid[i] - log_T_grid[i-1])
    r_frac = (log_target_rho - log_rho_grid[j-1]) / (log_rho_grid[j] - log_rho_grid[j-1])
    
    opacity = (1-t_frac)*(1-r_frac)*kappa[j-1,i-1] + \
              (1-t_frac)*r_frac*kappa[j,i-1] + \
              t_frac*(1-r_frac)*kappa[j-1,i] + \
              t_frac*r_frac*kappa[j,i]
    return opacity
```

## File Size
- **2D format**: 16 KB (884 values)

## Demo Script
- `demo_2d_opacity_table.py` - 2D format examples and interpolation
