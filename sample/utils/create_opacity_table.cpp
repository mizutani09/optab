// Program to create a lookup table HDF5 file from opacity calculation results
// Combines temperature, density, Planck mean, and Rosseland mean opacity data

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <hdf5.h>
#include <vector>
#include <iostream>
#include <dirent.h>
#include <string>
#include <algorithm>

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

        // Check if file has .h5 extension
        if (filename.length() > 3 &&
            filename.substr(filename.length() - 3) == ".h5") {

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

// Function to create the opacity lookup table HDF5 file
bool create_opacity_table(const std::vector<OpacityData>& data, const char* output_filename) {
    // Create new HDF5 file
    hid_t file_id = H5Fcreate(output_filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file_id < 0) {
        fprintf(stderr, "Error creating output file: %s\n", output_filename);
        return false;
    }

    hsize_t npoints = data.size();
    hsize_t dims[1] = {npoints};

    // Create dataspace
    hid_t space_id = H5Screate_simple(1, dims, NULL);
    if (space_id < 0) {
        H5Fclose(file_id);
        return false;
    }

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

    // Create and write temperature dataset
    hid_t temp_dataset = H5Dcreate2(file_id, "/temperature", H5T_NATIVE_DOUBLE, space_id,
                                    H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (temp_dataset >= 0) {
        H5Dwrite(temp_dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temperatures.data());
        H5Dclose(temp_dataset);
    }

    // Create and write density dataset
    hid_t rho_dataset = H5Dcreate2(file_id, "/density", H5T_NATIVE_DOUBLE, space_id,
                                   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (rho_dataset >= 0) {
        H5Dwrite(rho_dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, densities.data());
        H5Dclose(rho_dataset);
    }

    // Create and write Planck mean opacity dataset
    hid_t planck_dataset = H5Dcreate2(file_id, "/planck_mean_opacity", H5T_NATIVE_DOUBLE, space_id, 
                                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (planck_dataset >= 0) {
        H5Dwrite(planck_dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, planck_means.data());
        H5Dclose(planck_dataset);
    }

    // Create and write Rosseland mean opacity dataset
    hid_t ross_dataset = H5Dcreate2(file_id, "/rosseland_mean_opacity", H5T_NATIVE_DOUBLE, space_id, 
                                    H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (ross_dataset >= 0) {
        H5Dwrite(ross_dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, rosseland_means.data());
        H5Dclose(ross_dataset);
    }

    // Add attributes for units and description
    hid_t attr_space = H5Screate(H5S_SCALAR);
    hid_t str_type = H5Tcopy(H5T_C_S1);
    H5Tset_size(str_type, H5T_VARIABLE);

    // Temperature units
    if (H5Lexists(file_id, "/temperature", H5P_DEFAULT) > 0) {
        hid_t temp_id = H5Dopen(file_id, "/temperature", H5P_DEFAULT);
        hid_t attr = H5Acreate2(temp_id, "units", str_type, attr_space, H5P_DEFAULT, H5P_DEFAULT);
        const char* temp_units = "K";
        H5Awrite(attr, str_type, &temp_units);
        H5Aclose(attr);
        H5Dclose(temp_id);
    }

    // Density units
    if (H5Lexists(file_id, "/density", H5P_DEFAULT) > 0) {
        hid_t rho_id = H5Dopen(file_id, "/density", H5P_DEFAULT);
        hid_t attr = H5Acreate2(rho_id, "units", str_type, attr_space, H5P_DEFAULT, H5P_DEFAULT);
        const char* rho_units = "g/cm^3";
        H5Awrite(attr, str_type, &rho_units);
        H5Aclose(attr);
        H5Dclose(rho_id);
    }

    // Opacity units
    const char* opacity_units = "cm^2/g";
    if (H5Lexists(file_id, "/planck_mean_opacity", H5P_DEFAULT) > 0) {
        hid_t planck_id = H5Dopen(file_id, "/planck_mean_opacity", H5P_DEFAULT);
        hid_t attr = H5Acreate2(planck_id, "units", str_type, attr_space, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attr, str_type, &opacity_units);
        H5Aclose(attr);
        H5Dclose(planck_id);
    }

    if (H5Lexists(file_id, "/rosseland_mean_opacity", H5P_DEFAULT) > 0) {
        hid_t ross_id = H5Dopen(file_id, "/rosseland_mean_opacity", H5P_DEFAULT);
        hid_t attr = H5Acreate2(ross_id, "units", str_type, attr_space, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attr, str_type, &opacity_units);
        H5Aclose(attr);
        H5Dclose(ross_id);
    }

    H5Tclose(str_type);
    H5Sclose(attr_space);
    H5Sclose(space_id);
    H5Fclose(file_id);

    return true;
}

int main() {
    printf("=== Opacity Lookup Table Creator ===\n\n");

    // Find all H5 files
    std::vector<std::string> h5_files = find_h5_files("../table/output");

    if (h5_files.empty()) {
        printf("No HDF5 files found in ../table/output\n");
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

    // Create the lookup table
    const char* output_filename = "../table/output/opacity_lookup_table.h5";
    printf("Creating opacity lookup table: %s\n", output_filename);

    if (create_opacity_table(opacity_data, output_filename)) {
        printf("Successfully created opacity lookup table with %zu data points\n", opacity_data.size());

        // Print some statistics
        double min_temp = opacity_data[0].temperature, max_temp = opacity_data[0].temperature;
        double min_rho = opacity_data[0].density, max_rho = opacity_data[0].density;

        for (const auto& data : opacity_data) {
            if (data.temperature < min_temp) min_temp = data.temperature;
            if (data.temperature > max_temp) max_temp = data.temperature;
            if (data.density < min_rho) min_rho = data.density;
            if (data.density > max_rho) max_rho = data.density;
        }

        printf("\nData ranges:\n");
        printf("Temperature: %.3e - %.3e K\n", min_temp, max_temp);
        printf("Density: %.3e - %.3e g/cm³\n", min_rho, max_rho);

        // Show first few entries
        printf("\nFirst 5 data points:\n");
        printf("%-12s %-12s %-12s %-12s\n", "Temp [K]", "Rho [g/cm³]", "Planck", "Rosseland");
        for (int i = 0; i < std::min(5, (int)opacity_data.size()); i++) {
            printf("%-12.3e %-12.3e %-12.3e %-12.3e\n",
                   opacity_data[i].temperature, opacity_data[i].density,
                   opacity_data[i].planck_mean, opacity_data[i].rosseland_mean);
        }

    } else {
        printf("Error creating opacity lookup table\n");
        return 1;
    }

    return 0;
}
