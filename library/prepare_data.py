import os
import numpy as np
import pandas as pd
from yt import load
from library.filter_tools import create_multiple_filter_files

# function to prepare simulation data
def prepare_simulation_data(input_path, cell_fields = None, epf = None, 
                            filter_path = None, z = None, 
                            filter_dir = None, wl_initial = 0.6, 
                            wl_final = 2.3, num_bins = 20, 
                            jwst_filter_file = "F200W_filter.txt"):
    """"
    Prepare simulation dataset, star particle positions, and load filter data. 

    Parameters:
    ----------
    - input_path (str): path to the simulation dataset file
    - cell_fields (list or None): list of cell fields to load
    - epf (list or None): extra particle fields to load
    - filter_path (str or None): 'F200W_filter.txt' to load JWST filter, None to create filters automatically
    - z (float): redshift for wavelength shifting
    - filter_dir (str): directory for filter bins (default: 'filter_bins' in current directory)
    - wl_initial (float): starting wavelength in microns
    - wl_final (float): ending wavelength in microns
    - num_bins (int): number of filters to create
    - jwst_filter_file (str): JWST filter for resolution reference
    """

    # set default filter directory to current working directory if none provided
    if filter_dir is None:
        filter_dir = os.path.join(os.getcwd(), "filter_bins")

    # default fields
    if cell_fields is None:
        cell_fields = [
            "Density", "x-velocity", "y_velocity", "z-velocity", "Pressure",
            "Metallicity", "xHI", "xHII", "xHeII", "xHeIII"
        ]

    if epf is None:
        epf = [
            ("particle_family", "b"),
            ("particle_tag", "b"),
            ("particle_birth_epoch", "d"),
            ("particle_metallicity", "d"),
        ]

    # load dataset
    ds = load(input_path, fields = cell_fields, extra_particle_fields = epf, default_species_fields = "ionized")
    ad = ds.all_data()

    # star particle positions
    x_pos = np.array(ad["star", "particle_position_x"])
    y_pos = np.array(ad["star", "particle_position_y"])
    z_pos = np.array(ad["star", "particle_position_z"])
    x_center, y_center, z_center = np.mean(x_pos), np.mean(y_pos), np.mean(z_pos)
    ctr_at_code = np.array([x_center, y_center, z_center])

    # center positions and convert to pc
    x_pos -= x_center
    y_pos -= y_center
    z_pos -= z_center
    pop2_xyz = np.array(ds.arr(np.vstack([x_pos, y_pos, z_pos]), "code_length").to("pc")).T

    wavelength_filter_shifted = []
    output_filter = []

    # create the filter for users if they don't have it
    if filter_path is None:
        print(f"no filter file provided... creating filters automatically in {filter_dir}")

        # uses the function from your other file
        bin_edge, filter_files = create_multiple_filter_files(
            filter_dir, wl_initial, wl_final, num_bins, jwst_filter_file
        )

        for i, f in enumerate(filter_files):
            filter_data = np.loadtxt(f, skiprows = 1)
            df_filter = pd.DataFrame(filter_data, columns = ["Wavelength [Microns]", "Output"])
            
            # convert microns to Angstroms and shift by redshift
            shifted_wl = df_filter["Wavelength [Microns]"] * 1e4 / (1 + z)
            
            wavelength_filter_shifted.append(shifted_wl.values)
            output_filter.append(df_filter["Output"].values)
            
            print(f"loaded filter {i+1}: {os.path.basename(f)}")
        
        wavelength_filter_shifted = np.array(wavelength_filter_shifted)
        output_filter = np.array(output_filter)

    # use existing filter file
    else:
        if not os.path.isfile(filter_path):
            raise FileNotFoundError(f"filter file '{filter_path}' does not exist.")

        print(f"loading filter data from file: {filter_path}")
        filter_data = np.loadtxt(filter_path, skiprows = 1)
        df_filter = pd.DataFrame(filter_data, columns = ["Wavelength [Microns]", "Output"])
        wavelength_filter_shifted = df_filter["Wavelength [Microns]"].values * 1e4 / (1 + z)
        output_filter = df_filter["Output"].values

    plt_wdth = 400 

    return ds, ad, ctr_at_code, pop2_xyz, wavelength_filter_shifted, output_filter, plt_wdth