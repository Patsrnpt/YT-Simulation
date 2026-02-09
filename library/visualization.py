# function to visualize projection plot
def create_projection_plot(ds, field_name, ctr_at_code, z, distance_pc, plt_wdth = 400, plot_units = "flux", axis_units = "arcsec", filter_path = None):
    """"
    Create a projection plot of a specific flux field.
    """

    # create projection to get data
    proj = yt.ProjectionPlot(ds, "z", ("gas", field_name), width = (plt_wdth, "pc"), center = ctr_at_code)
    p_frb = proj.frb
    flux_data = np.array(p_frb["gas", field_name])

    # handle filter width
    delta_wl = 1.988e4 
    if filter_path:
        delta_wl = calculate_filter_width(filter_path)

    # unit conversion logic
    if plot_units == "flux" or plot_units == None:
        plot_data = np.log10(flux_data)
        colorbar_label = "log10(flux [erg s$^{-1}$ cm$^{-2}$])"
        cmap = "viridis"
    elif plot_units == "jy_arcsec2":
        flux_per_A = flux_data / delta_wl
        sb = flux_to_surface_brightness(flux_per_A, z = z)
        plot_data = np.log10(sb) + 9
        colorbar_label = "log10(nJy arcsec$^{-2}$)"
        cmap = "viridis"
    elif plot_units == "magnitude_arcsec2":
        flux_per_A = flux_data / delta_wl
        sb = flux_to_surface_brightness(flux_per_A, z = z)
        plot_data = flux_to_magnitude(sb, z = z)
        colorbar_label = "magnitude arcsec$^{-2}$"
        cmap = "viridis_r"

    # axis extent
    if axis_units == "arcsec":
        pc_to_cm = 3.086e18
        ang_rad = ((plt_wdth * pc_to_cm) / 2) / (distance_pc * pc_to_cm) * (180 / np.pi) * 3600
        extent = [-ang_rad, ang_rad, -ang_rad, ang_rad]
    elif axis_units == "pc":
        extent = [-plt_wdth / 2, plt_wdth / 2, -plt_wdth / 2, plt_wdth / 2]

    fig, ax = plt.subplots(figsize = (10, 8))
    im = ax.imshow(plot_data, extent = extent, cmap = cmap, origin = "lower")
    ax.set_xlabel(axis_units)
    ax.set_ylabel(axis_units)
    plot_title = f"projection plot: {field_name} ({plot_units})"
    ax.set_title(plot_title, fontsize = 14, pad = 20)
    plt.colorbar(im, label = colorbar_label)
    plt.show()

    return None

# function to visualize phase plot
def create_phase_plot(ad, x_field = "density", y_field = "temperature", z_field = "flux_total", x_bins = None, y_bins = None, weight_field = None):
    """"
    Create a phase plot showing distribution of flux across two variables.
    """

    plot = yt.PhasePlot(ad, x_field, y_field, [z_field], weight_field = weight_field)
    
    if x_bins:
        plot.set_xlim(x_bins[0], x_bins[1])
    if y_bins:
        plot.set_ylim(y_bins[0], y_bins[1])
        
    plot.set_title(z_field, f"phase plot: {z_field} distribution")
    plot.set_cmap(z_field, "viridis")
    plot.show()

    return None

# function to calculate related variable to the spectrum
def create_spectrum_plot(ds, filter_list, z, filter_dir, flux_type = "total", plot_units = "flux"):
    """"
    Create a spectrum plot from multiple filter bins.
    """

    centers_wl = []
    y_values = []

    for f_num in filter_list:
        if flux_type == "total":
            field = f"flux_total_filter_{f_num}"
        else:
            field = f"flux_{flux_type}_filter_{f_num}"

        if ("gas", field) not in ds.derived_field_list:
            continue
            
        total_flux = float(ds.all_data()[("gas", field)].sum())
        
        file_path = os.path.join(filter_dir, f"filter_{f_num}.txt")
        if not os.path.exists(file_path):
            continue
            
        f_data = np.loadtxt(file_path, skiprows = 1)
        mask = f_data[:, 1] > 0
        wls = f_data[:, 0][mask]
        
        center = (wls.max() + wls.min()) / 2
        delta_wl = (wls.max() - wls.min()) * 1e4
        
        centers_wl.append(center)
        
        if plot_units == "flux" or plot_units == None:
            y_values.append(total_flux)
        else:
            flux_per_A = total_flux / delta_wl
            sb = flux_to_surface_brightness(flux_per_A, z = z)
            
            if plot_units == "jy_arcsec2":
                y_values.append(sb * 1e9) 
            elif plot_units == "magnitude_arcsec2":
                y_values.append(flux_to_magnitude(sb, z = z))

    # sort the data by wavelength for a proper spectrum line
    sorted_indices = np.argsort(centers_wl)
    centers_wl = np.array(centers_wl)[sorted_indices]
    y_values = np.array(y_values)[sorted_indices]

    return centers_wl, y_values, plot_units