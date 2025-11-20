// Program to extract opacity means, temperature, density, and calculate pressure from HDF5 files

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <hdf5.h>
#include <vector>
#include <iostream>
#include <dirent.h>
#include <string>
#include <algorithm>

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

    printf("Found %zu H5 files in %s\n", h5_files.size(), folder_path);

    return h5_files;
}

// Function to extract temperature, density, and opacity means
void extract_all_data(hid_t file_id, const char* filename) {
    printf("=== Processing %s ===\n", filename);

    // Extract temperature
    if (H5Lexists(file_id, "/temp", H5P_DEFAULT) > 0) {
        hid_t dataset_id = H5Dopen(file_id, "/temp", H5P_DEFAULT);
        if (dataset_id >= 0) {
            double temp;
            if (H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &temp) >= 0) {
                printf("Temperature: %g K\n", temp);
            }
            H5Dclose(dataset_id);
        }
    }

    // Extract density
    if (H5Lexists(file_id, "/rho", H5P_DEFAULT) > 0) {
        hid_t dataset_id = H5Dopen(file_id, "/rho", H5P_DEFAULT);
        if (dataset_id >= 0) {
            double rho;
            if (H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &rho) >= 0) {
                printf("Density: %g g/cm³\n", rho);
            }
            H5Dclose(dataset_id);
        }
    }

    // Extract opacity means
    const char* opacity_datasets[] = {"/pla", "/ros"};
    const char* opacity_names[] = {"Planck mean opacity", "Rosseland mean opacity"};

    for (int i = 0; i < 2; i++) {
        const char* dataset_name = opacity_datasets[i];

        if (H5Lexists(file_id, dataset_name, H5P_DEFAULT) > 0) {
            hid_t dataset_id = H5Dopen(file_id, dataset_name, H5P_DEFAULT);
            if (dataset_id >= 0) {
                double opacity;
                if (H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &opacity) >= 0) {
                    printf("%s: %g cm²/g\n", opacity_names[i], opacity);
                }
                H5Dclose(dataset_id);
            }
        }
    }
    printf("\n");
}

// Function to process a single H5 file
void process_h5_file(const char* filepath) {
    hid_t file_id = H5Fopen(filepath, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id < 0) {
        fprintf(stderr, "Error opening file: %s\n", filepath);
        return;
    }

    // Extract all data: temperature, density, pressure, and opacity means
    extract_all_data(file_id, filepath);

    H5Fclose(file_id);
}

int main() {
    printf("=== Stellar Opacity Data Extractor ===\n\n");

    // // Process single file example
    // const char* single_file = "../table/output/mono_00314.h5";
    // printf("Processing single file example:\n");
    // process_h5_file(single_file);

    // Process all files in directory
    printf("\nProcessing first 5 files from directory:\n");
    std::vector<std::string> h5_files = find_h5_files("../table/output");

    if (!h5_files.empty()) {
        // Process only the first 5 files to avoid too much output
        int files_to_process = std::min(5, (int)h5_files.size());
        for (int i = 0; i < files_to_process; i++) {
            process_h5_file(h5_files[i].c_str());
        }

        if (h5_files.size() > 5) {
            printf("... and %zu more files (not shown)\n", h5_files.size() - 5);
        }
    }

    printf("Data extraction completed.\n");
    return 0;
}
