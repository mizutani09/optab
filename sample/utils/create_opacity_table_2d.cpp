// Program to create a 2D lookup table HDF5 file from opacity calculation results
// Creates efficient 2D grid structure: unique temperature array, unique density array,
// and 2D opacity arrays indexed by [temp_index][density_index]

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <hdf5.h>
#include <vector>
#include <iostream>
#include <dirent.h>
#include <string>
#include <algorithm>
#include <map>
#include <set>

struct OpacityData {
    double temperature;
    double density;
    double planck_mean;
    double rosseland_mean;
};

// Function to find all H5 files in a directory
std::vector<std::string> find_h5_files(const char* folder_path) {
    std::vector<std::string> h5_files;

    DIR* dir = opendir(folder_path);
    if (dir == NULL) {
        printf("Error: Cannot open directory %s\n", folder_path);
        return h5_files;
    }

    struct dirent* entry;
    while ((entry = readdir(dir)) != NULL) {
        std::string filename(entry->d_name);

        // Check if file has .h5 extension and exclude the lookup table itself
        if (filename.length() > 3 &&
            filename.substr(filename.length() - 3) == ".h5" &&
            filename != "opacity_lookup_table.h5" &&
            filename != "opacity_lookup_table_2d.h5") {

            std::string full_path = std::string(folder_path) + "/" + filename;
            h5_files.push_back(full_path);
        }
    }

    closedir(dir);

    // Sort the files
    std::sort(h5_files.begin(), h5_files.end());

    return h5_files;
}

// Function to extract data from a single HDF5 file
bool extract_data_from_file(const char* filepath, OpacityData& data) {
    hid_t file_id = H5Fopen(filepath, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id < 0) {
        fprintf(stderr, "Error opening file: %s\n", filepath);
        return false;
    }

    bool success = true;

    // Extract temperature
    if (H5Lexists(file_id, "/temp", H5P_DEFAULT) > 0) {
        hid_t dataset_id = H5Dopen(file_id, "/temp", H5P_DEFAULT);
        if (dataset_id >= 0) {
            if (H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data.temperature) < 0) {
                success = false;
            }
            H5Dclose(dataset_id);
        } else {
            success = false;
        }
    } else {
        success = false;
    }

    // Extract density
    if (success && H5Lexists(file_id, "/rho", H5P_DEFAULT) > 0) {
        hid_t dataset_id = H5Dopen(file_id, "/rho", H5P_DEFAULT);
        if (dataset_id >= 0) {
            if (H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data.density) < 0) {
                success = false;
            }
            H5Dclose(dataset_id);
        } else {
            success = false;
        }
    } else {
        success = false;
    }

    // Extract Planck mean opacity
    if (success && H5Lexists(file_id, "/pla", H5P_DEFAULT) > 0) {
        hid_t dataset_id = H5Dopen(file_id, "/pla", H5P_DEFAULT);
        if (dataset_id >= 0) {
            if (H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data.planck_mean) < 0) {
                success = false;
            }
            H5Dclose(dataset_id);
        } else {
            success = false;
        }
    } else {
        success = false;
    }

    // Extract Rosseland mean opacity
    if (success && H5Lexists(file_id, "/ros", H5P_DEFAULT) > 0) {
        hid_t dataset_id = H5Dopen(file_id, "/ros", H5P_DEFAULT);
        if (dataset_id >= 0) {
            if (H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data.rosseland_mean) < 0) {
                success = false;
            }
            H5Dclose(dataset_id);
        } else {
            success = false;
        }
    } else {
        success = false;
    }

    H5Fclose(file_id);
    return success;
}

// Function to create the 2D opacity lookup table HDF5 file
bool create_2d_opacity_table(const std::vector<OpacityData>& data, const char* output_filename) {
    // Organize data to understand the actual structure
    std::map<double, std::vector<OpacityData>> temp_groups;

    // Group data by temperature
    for (const auto& d : data) {
        temp_groups[d.temperature].push_back(d);
    }

    printf("Analyzing data structure:\n");
    printf("Found %zu temperature groups\n", temp_groups.size());

    // Check if all temperature groups have the same number of density points
    size_t expected_density_count = temp_groups.begin()->second.size();
    bool uniform_structure = true;

    for (const auto& group : temp_groups) {
        if (group.second.size() != expected_density_count) {
            uniform_structure = false;
            break;
        }
    }

    if (!uniform_structure) {
        printf("Warning: Non-uniform grid structure detected. Using 1D storage format.\n");

        // Fall back to 1D storage (more memory efficient for irregular grids)
        hid_t file_id = H5Fcreate(output_filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        if (file_id < 0) {
            fprintf(stderr, "Error creating output file: %s\n", output_filename);
            return false;
        }

        hsize_t npoints = data.size();
        hsize_t dims[1] = {npoints};
        hid_t space_id = H5Screate_simple(1, dims, NULL);

        // Prepare data arrays
        std::vector<double> temperatures(npoints);
        std::vector<double> densities(npoints);
        std::vector<double> planck_means(npoints);
        std::vector<double> rosseland_means(npoints);

        for (size_t i = 0; i < npoints; i++) {
            temperatures[i] = data[i].temperature;
            densities[i] = data[i].density;
            planck_means[i] = data[i].planck_mean;
            rosseland_means[i] = data[i].rosseland_mean;
        }

        // Create and write datasets with 1D structure
        hid_t temp_dataset = H5Dcreate2(file_id, "/temperature", H5T_NATIVE_DOUBLE, space_id,
                                        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if (temp_dataset >= 0) {
            H5Dwrite(temp_dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temperatures.data());
            H5Dclose(temp_dataset);
        }

        hid_t rho_dataset = H5Dcreate2(file_id, "/density", H5T_NATIVE_DOUBLE, space_id,
                                       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if (rho_dataset >= 0) {
            H5Dwrite(rho_dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, densities.data());
            H5Dclose(rho_dataset);
        }

        hid_t planck_dataset = H5Dcreate2(file_id, "/planck_mean_opacity", H5T_NATIVE_DOUBLE, space_id,
                                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if (planck_dataset >= 0) {
            H5Dwrite(planck_dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, planck_means.data());
            H5Dclose(planck_dataset);
        }

        hid_t ross_dataset = H5Dcreate2(file_id, "/rosseland_mean_opacity", H5T_NATIVE_DOUBLE, space_id,
                                        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if (ross_dataset >= 0) {
            H5Dwrite(ross_dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, rosseland_means.data());
            H5Dclose(ross_dataset);
        }

        // Add structure information
        hid_t attr_space = H5Screate(H5S_SCALAR);
        hid_t str_type = H5Tcopy(H5T_C_S1);
        H5Tset_size(str_type, H5T_VARIABLE);

        hid_t structure_attr = H5Acreate2(file_id, "data_structure", str_type, attr_space, H5P_DEFAULT, H5P_DEFAULT);
        const char* structure_desc = "1D arrays: irregular grid, memory optimized";
        H5Awrite(structure_attr, str_type, &structure_desc);
        H5Aclose(structure_attr);

        H5Tclose(str_type);
        H5Sclose(attr_space);
        H5Sclose(space_id);
        H5Fclose(file_id);

        printf("Created 1D lookup table (irregular grid, memory optimized)\n");
        return true;
    }

    // Create true 2D grid structure
    printf("Creating true 2D grid structure:\n");

    // Extract unique temperatures and densities (assumes uniform structure)
    std::vector<double> temperatures;
    std::vector<double> densities;
    std::vector<double> log_temperatures;
    std::vector<double> log_densities;

    for (const auto& group : temp_groups) {
        temperatures.push_back(group.first);
        log_temperatures.push_back(log10(group.first));
    }

    // Get densities from the first temperature group (should be same for all)
    auto first_group = temp_groups.begin()->second;
    std::sort(first_group.begin(), first_group.end(),
              [](const OpacityData& a, const OpacityData& b) { return a.density < b.density; });

    for (const auto& d : first_group) {
        densities.push_back(d.density);
        log_densities.push_back(log10(d.density));
    }

    size_t n_temp = temperatures.size();
    size_t n_rho = densities.size();

    printf("Grid: %zu temperatures x %zu densities = %zu points\n", n_temp, n_rho, n_temp * n_rho);

    // Create 2D arrays for opacities (density first, then temperature)
    std::vector<std::vector<double>> planck_2d(n_rho, std::vector<double>(n_temp, 0.0));
    std::vector<std::vector<double>> rosseland_2d(n_rho, std::vector<double>(n_temp, 0.0));

    // Fill the 2D arrays [density_idx][temp_idx]
    for (size_t i = 0; i < n_temp; i++) {
        double temp = temperatures[i];
        auto& temp_data = temp_groups[temp];

        // Sort by density
        std::sort(temp_data.begin(), temp_data.end(),
                  [](const OpacityData& a, const OpacityData& b) { return a.density < b.density; });

        for (size_t j = 0; j < n_rho; j++) {
            if (j < temp_data.size()) {
                planck_2d[j][i] = temp_data[j].planck_mean;      // [density_idx][temp_idx]
                rosseland_2d[j][i] = temp_data[j].rosseland_mean; // [density_idx][temp_idx]
            }
        }
    }

    // Create new HDF5 file
    hid_t file_id = H5Fcreate(output_filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file_id < 0) {
        fprintf(stderr, "Error creating output file: %s\n", output_filename);
        return false;
    }

    // Create 1D dataspaces for temperature and density arrays
    hsize_t temp_dims[1] = {n_temp};
    hsize_t rho_dims[1] = {n_rho};
    hid_t temp_space = H5Screate_simple(1, temp_dims, NULL);
    hid_t rho_space = H5Screate_simple(1, rho_dims, NULL);

    // Create 2D dataspace for opacity arrays (density × temperature)
    hsize_t opacity_dims[2] = {n_rho, n_temp};
    hid_t opacity_space = H5Screate_simple(2, opacity_dims, NULL);

    // Create and write temperature dataset (log10 values)
    hid_t temp_dataset = H5Dcreate2(file_id, "/log_temperature", H5T_NATIVE_DOUBLE, temp_space,
                                    H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (temp_dataset >= 0) {
        H5Dwrite(temp_dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, log_temperatures.data());
        H5Dclose(temp_dataset);
    }

    // Create and write density dataset (log10 values)
    hid_t rho_dataset = H5Dcreate2(file_id, "/log_density", H5T_NATIVE_DOUBLE, rho_space,
                                   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (rho_dataset >= 0) {
        H5Dwrite(rho_dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, log_densities.data());
        H5Dclose(rho_dataset);
    }

    // Flatten 2D arrays for HDF5 storage (row-major order: density × temperature)
    std::vector<double> planck_flat(n_rho * n_temp);
    std::vector<double> rosseland_flat(n_rho * n_temp);

    for (size_t i = 0; i < n_rho; i++) {
        for (size_t j = 0; j < n_temp; j++) {
            planck_flat[i * n_temp + j] = planck_2d[i][j];
            rosseland_flat[i * n_temp + j] = rosseland_2d[i][j];
        }
    }

    // Create and write Planck mean opacity dataset (2D)
    hid_t planck_dataset = H5Dcreate2(file_id, "/planck_mean_opacity", H5T_NATIVE_DOUBLE, opacity_space, 
                                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (planck_dataset >= 0) {
        H5Dwrite(planck_dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, planck_flat.data());
        H5Dclose(planck_dataset);
    }

    // Create and write Rosseland mean opacity dataset (2D)
    hid_t ross_dataset = H5Dcreate2(file_id, "/rosseland_mean_opacity", H5T_NATIVE_DOUBLE, opacity_space, 
                                    H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (ross_dataset >= 0) {
        H5Dwrite(ross_dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, rosseland_flat.data());
        H5Dclose(ross_dataset);
    }

    // Add attributes for units and description
    hid_t attr_space = H5Screate(H5S_SCALAR);
    hid_t str_type = H5Tcopy(H5T_C_S1);
    H5Tset_size(str_type, H5T_VARIABLE);

    // Log Temperature units and description
    if (H5Lexists(file_id, "/log_temperature", H5P_DEFAULT) > 0) {
        hid_t temp_id = H5Dopen(file_id, "/log_temperature", H5P_DEFAULT);
        hid_t attr = H5Acreate2(temp_id, "units", str_type, attr_space, H5P_DEFAULT, H5P_DEFAULT);
        const char* temp_units = "log10(K)";
        H5Awrite(attr, str_type, &temp_units);
        H5Aclose(attr);

        attr = H5Acreate2(temp_id, "description", str_type, attr_space, H5P_DEFAULT, H5P_DEFAULT);
        const char* temp_desc = "Log10 of temperature grid points";
        H5Awrite(attr, str_type, &temp_desc);
        H5Aclose(attr);
        H5Dclose(temp_id);
    }

    // Log Density units and description
    if (H5Lexists(file_id, "/log_density", H5P_DEFAULT) > 0) {
        hid_t rho_id = H5Dopen(file_id, "/log_density", H5P_DEFAULT);
        hid_t attr = H5Acreate2(rho_id, "units", str_type, attr_space, H5P_DEFAULT, H5P_DEFAULT);
        const char* rho_units = "log10(g/cm^3)";
        H5Awrite(attr, str_type, &rho_units);
        H5Aclose(attr);

        attr = H5Acreate2(rho_id, "description", str_type, attr_space, H5P_DEFAULT, H5P_DEFAULT);
        const char* rho_desc = "Log10 of density grid points";
        H5Awrite(attr, str_type, &rho_desc);
        H5Aclose(attr);
        H5Dclose(rho_id);
    }

    // Opacity units and descriptions
    const char* opacity_units = "cm^2/g";
    if (H5Lexists(file_id, "/planck_mean_opacity", H5P_DEFAULT) > 0) {
        hid_t planck_id = H5Dopen(file_id, "/planck_mean_opacity", H5P_DEFAULT);
        hid_t attr = H5Acreate2(planck_id, "units", str_type, attr_space, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attr, str_type, &opacity_units);
        H5Aclose(attr);

        attr = H5Acreate2(planck_id, "description", str_type, attr_space, H5P_DEFAULT, H5P_DEFAULT);
        const char* planck_desc = "Planck mean opacity [log_density, log_temperature]";
        H5Awrite(attr, str_type, &planck_desc);
        H5Aclose(attr);
        H5Dclose(planck_id);
    }

    if (H5Lexists(file_id, "/rosseland_mean_opacity", H5P_DEFAULT) > 0) {
        hid_t ross_id = H5Dopen(file_id, "/rosseland_mean_opacity", H5P_DEFAULT);
        hid_t attr = H5Acreate2(ross_id, "units", str_type, attr_space, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attr, str_type, &opacity_units);
        H5Aclose(attr);

        attr = H5Acreate2(ross_id, "description", str_type, attr_space, H5P_DEFAULT, H5P_DEFAULT);
        const char* ross_desc = "Rosseland mean opacity [log_density, log_temperature]";
        H5Awrite(attr, str_type, &ross_desc);
        H5Aclose(attr);
        H5Dclose(ross_id);
    }

    // Add global attributes
    hid_t grid_attr = H5Acreate2(file_id, "grid_structure", str_type, attr_space, H5P_DEFAULT, H5P_DEFAULT);
    const char* grid_desc = "2D grid: opacity[i,j] corresponds to log_density[i] and log_temperature[j]";
    H5Awrite(grid_attr, str_type, &grid_desc);
    H5Aclose(grid_attr);

    H5Tclose(str_type);
    H5Sclose(attr_space);
    H5Sclose(temp_space);
    H5Sclose(rho_space);
    H5Sclose(opacity_space);
    H5Fclose(file_id);

    return true;
}

int main() {
    printf("=== 2D Opacity Lookup Table Creator ===\n\n");

    // variable for the path
    const char* folder_path = "../table/output";

    // Find all H5 files
    std::vector<std::string> h5_files = find_h5_files(folder_path);

    if (h5_files.empty()) {
        printf("No HDF5 files found in %s\n", folder_path);
        return 1;
    }

    printf("Found %zu HDF5 files\n", h5_files.size());
    printf("Extracting data from all files...\n");

    // Extract data from all files
    std::vector<OpacityData> opacity_data;
    int successful_reads = 0;

    for (const auto& filepath : h5_files) {
        OpacityData data;
        if (extract_data_from_file(filepath.c_str(), data)) {
            opacity_data.push_back(data);
            successful_reads++;

            // Progress indicator
            if (successful_reads % 50 == 0) {
                printf("Processed %d files...\n", successful_reads);
            }
        } else {
            printf("Warning: Failed to extract data from %s\n", filepath.c_str());
        }
    }

    printf("Successfully extracted data from %d files\n", successful_reads);

    if (opacity_data.empty()) {
        printf("No valid data extracted. Exiting.\n");
        return 1;
    }

    // Create the 2D lookup table
    std::string output_filename_pre = std::string(folder_path) + "/opacity_lookup_table_2d.h5";
    const char* output_filename = output_filename_pre.c_str();

    printf("Creating 2D opacity lookup table: %s\n", output_filename);

    if (create_2d_opacity_table(opacity_data, output_filename)) {
        printf("Successfully created opacity lookup table\n");

        // Print some statistics
        std::set<double> temp_set, rho_set;
        for (const auto& d : opacity_data) {
            temp_set.insert(d.temperature);
            rho_set.insert(d.density);
        }

        std::vector<double> temps(temp_set.begin(), temp_set.end());
        std::vector<double> rhos(rho_set.begin(), rho_set.end());

        printf("\nFinal structure:\n");
        printf("Temperature points: %zu\n", temps.size());
        printf("Total data points: %zu\n", opacity_data.size());
        printf("Data ranges:\n");
        printf("Temperature: %.3e - %.3e K\n", temps.front(), temps.back());

        if (opacity_data.size() == temps.size() * 21) {
            printf("Density points per temperature: 21\n");
            printf("Structure: Regular 2D grid (%zu x 21)\n", temps.size());
        } else {
            printf("Structure: Optimized 1D storage for irregular grid\n");
        }

    } else {
        printf("Error creating opacity lookup table\n");
        return 1;
    }

    return 0;
}
