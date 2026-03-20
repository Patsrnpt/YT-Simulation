import os
import numpy as np
from library.calculate_quantities import calculate_original_resolution

# function to create the filter without redshifting
def create_multiple_filter_files(output_dir, wl_initial, wl_final, num_filters = 20, jwst_filter_path = "F200W_filter.txt"):
    """"
    Create multiple filter files. Each filter covers the full wavelength range,
    but transmission = 1 only in its specific bin, 0 elsewhere.

    Parameters:
    ----------
    output_dir (str): location to save these filter files.
    wl_initial (float): initial wavelength of the full wavelength range (in microns).
    wl_final (float): final wavelength of the full wavelength range (in microns).
    num_filters (int): number of filters you want to create.
    jwst_filter_path (str): location of your JWST filter. This will be used for calculating resolution, which will be matched with filter's resolution created using this function. 
    This means that this function requires you to download JWST filter from this link: https://jwst-docs.stsci.edu/jwst-near-infrared-camera/nircam-instrumentation/nircam-filters#gsc.tab=0

    Outputs:
    -------
    bin_edge (np.ndarray): wavelength boundaries for each filter bin  
    jwst_filter_paths (str): file paths for each created filter file
    """

    # if you choose to create store this filter, this part will create output directory
    os.makedirs(output_dir, exist_ok = True)

    # store values of JWST filter data
    jwst_data = np.loadtxt(jwst_filter_path, skiprows = 1)
    jwst_wavelength = jwst_data[:, 0]  # wavelength in microns

    # calculate resolution of the jwst filter
    jwst_resolution = calculate_original_resolution(jwst_wavelength) # this resolution will be used as a resolution of each filter

    # define balmer break wavelength
    balmer_wl = 0.3645 # in microns

    # check if balmer break is in the wavelength range
    if (wl_initial < balmer_wl) and (balmer_wl < wl_final):

        # convert to log space
        log_wl_initial = np.log10(wl_initial)
        log_wl_final = np.log10(wl_final)
        log_balmer_wl = np.log10(balmer_wl)

        # calculate fraction in log space
        f1 = (log_balmer_wl - log_wl_initial)/(log_wl_final - log_wl_initial)
        f2 = (log_wl_final - log_balmer_wl)/(log_wl_final - log_wl_initial)

        # calculate number of filters in each region (flooring)
        N1 = int(num_filters * f1)
        N2 = int(num_filters * f2)

        # create wavelength grid for each region
        wavelength_filter_grid_1 = np.logspace(log_wl_initial, log_balmer_wl, N1 + 1)
        wavelength_filter_grid_2 = np.logspace(log_balmer_wl, log_wl_final, N2 + 1)

        # combine grid and remove duplicate balmer point
        wavelength_filter_grid = np.concatenate((wavelength_filter_grid_1[:-1], wavelength_filter_grid_2))

    else:
        
        # wavelength grid of multiple filters to calculate spacing between each data point
        wavelength_filter_grid = np.logspace(np.log10(wl_initial), np.log10(wl_final), num_filters + 1) # this variable will tell you where the bin edge is

    # calculate number of data points in each filter to match with jwst resolution using the first filter
    central_wavelength = (wavelength_filter_grid[0] + wavelength_filter_grid[1])/2
    delta_wavelength = wavelength_filter_grid[1] - wavelength_filter_grid[0]
    num_points = 1 + (delta_wavelength * jwst_resolution)/central_wavelength
    num_points = int(np.round(num_points))
    
    # calculate spacing between data points
    wavelength_spacing = (wavelength_filter_grid[1] - wavelength_filter_grid[0])/(num_points - 1)

    # create wavelength array using linspace to include last point exactly
    wl_array = np.linspace(wl_initial, wl_final, int((wl_final - wl_initial)/wavelength_spacing) + 1) # x-axis of all filters

    # empty array of filters
    filter_paths = []

    # for loop to create filters
    for i in range(len(wavelength_filter_grid) - 1):
        # empty array of transmission
        transmission = np.zeros_like(wl_array)

        # consider range of each filter
        if i == len(wavelength_filter_grid) - 2:
            # include last point for the last bin
            mask = (wl_array >= wavelength_filter_grid[i]) & (wl_array <= wavelength_filter_grid[i + 1])
        else:
            mask = (wl_array >= wavelength_filter_grid[i]) & (wl_array < wavelength_filter_grid[i + 1])

        # transmission in the consider range = 1
        transmission[mask] = 1

        # skip empty filters (in case flooring produced 0 points)
        if np.any(transmission):
            # save filter files
            file_name = f"filter_{i+1}.txt"
            file_path = os.path.join(output_dir, file_name)

            # save as two columns: wavelength and transmission
            data_to_save = np.column_stack((wl_array, transmission))
            np.savetxt(file_path, data_to_save, header = "Wavelength[Microns] Transmission", fmt = "%.6f")
            filter_paths.append(file_path)

    print(f"Successfully created {len(filter_paths)} filters in {output_dir}")

    return wavelength_filter_grid, filter_paths