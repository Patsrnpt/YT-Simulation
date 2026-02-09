# function to add flux fields to yt package
def add_flux_fields(ds, interp_funcs, min_temp, max_temp, min_dense, max_dense, he_h_ratio = 0.1):
    """"
    Add interpolated flux fields (bound-free, two-photon, free-free, total) to the yt dataset.
    Automatically detects and processes single or multiple filters.
    """
    
    # check if we have multiple filters by checking if keys start with "filter_"
    has_multiple_filters = any(key.startswith("filter_") for key in interp_funcs.keys())
    
    if has_multiple_filters:
        # multiple filters case
        filter_keys = [key for key in interp_funcs.keys() if key.startswith("filter_")]
        filter_keys.sort()
        n_filters = len(filter_keys)
        
        print(f"processing {n_filters} filters...")
        
        filter_list = []
        
        # process each filter
        for filter_key in filter_keys:
            filter_number = filter_key.replace("filter_", "")
            filter_list.append(filter_number)
            
            print(f"\nadding fields for filter {filter_number}...")
            
            # get interpolation functions for this specific filter
            current_interp_funcs = interp_funcs[filter_key]
            
            # add individual component fields for this filter
            for comp in ["contH", "cont2p", "contff"]:
                if comp in current_interp_funcs:
                    
                    def _flux_field(field, data, component = comp, interp = current_interp_funcs[comp]):
                        # clip temperature and density to valid ranges
                        temperature = np.clip(data["gas", "temperature"].value, min_temp, max_temp)
                        number_density = np.clip(data["gas", "number_density"].value, min_dense, max_dense)
                        
                        # create mask for valid density values
                        mask = (number_density <= 1e10)
                        
                        # get hydrogen and electron densities
                        nH_tot = data["gas", "H_nuclei_density"].value
                        xHII = data["ramses", "xHII"].value
                        xHeII = data["ramses", "xHeII"].value
                        xHeIII = data["ramses", "xHeIII"].value
                        
                        n_p = xHII * nH_tot 
                        n_e = nH_tot * (xHII + he_h_ratio * (xHeII + 2 * xHeIII)) 
                        
                        # prepare points for interpolation
                        points = np.column_stack((temperature.flatten(), number_density.flatten()))
                        
                        # interpolate flux
                        flux = interp(points)
                        flux = flux.reshape(n_p.shape)
                        
                        # calculate final flux: flux * n_p * n_e
                        flux = flux * n_p * n_e
                        
                        # set flux to 0 where density is invalid
                        flux[~mask] = 0.0
                        
                        return flux.reshape(temperature.shape)
                    
                    # add field to dataset
                    field_name = f"flux_{comp}_filter_{filter_number}"
                    ds.add_field(
                        ("gas", field_name),
                        function = _flux_field,
                        sampling_type = "cell",
                        units = "",
                        force_override = True
                    )
                    
                    print(f"  added: {field_name}")
            
            # add total flux field for this filter
            def _flux_total_filter(field, data, fnum = filter_number):
                return (data["gas", f"flux_contH_filter_{fnum}"] +
                        data["gas", f"flux_cont2p_filter_{fnum}"] +
                        data["gas", f"flux_contff_filter_{fnum}"])
            
            ds.add_field(
                ("gas", f"flux_total_filter_{filter_number}"),
                function = _flux_total_filter,
                sampling_type = "cell",
                units = "",
                force_override = True
            )
            
            print(f"  added: flux_total_filter_{filter_number}")
        
        print(f"\nsuccessfully added fields for {n_filters} filters")
        return ds, filter_list
        
    else:
        # single filter case
        print("processing single filter...")
        filter_list = ["01"] 
        
        for comp in ["contH", "cont2p", "contff"]:
            if comp in interp_funcs:
                
                def _flux_field(field, data, component = comp, interp = interp_funcs[comp]):
                    temperature = np.clip(data["gas", "temperature"].value, min_temp, max_temp)
                    number_density = np.clip(data["gas", "number_density"].value, min_dense, max_dense)
                    mask = (number_density <= 1e10)
                    
                    nH_tot = data["gas", "H_nuclei_density"].value
                    xHII = data["ramses", "xHII"].value
                    xHeII = data["ramses", "xHeII"].value
                    xHeIII = data["ramses", "xHeIII"].value
                    
                    n_p = xHII * nH_tot 
                    n_e = nH_tot * (xHII + he_h_ratio * (xHeII + 2 * xHeIII)) 
                    
                    points = np.column_stack((temperature.flatten(), number_density.flatten()))
                    flux = interp(points)
                    flux = flux.reshape(n_p.shape)
                    flux = flux * n_p * n_e
                    flux[~mask] = 0.0
                    
                    return flux.reshape(temperature.shape)
                
                field_name = f"flux_{comp}"
                ds.add_field(
                    ("gas", field_name),
                    function = _flux_field,
                    sampling_type = "cell",
                    units = "",
                    force_override = True
                )
                print(f"added: {field_name}")
        
        def _flux_total_field(field, data):
            return (data["gas", "flux_contH"] +
                    data["gas", "flux_total_2p"] + # Note: fixed potential typo from previous version
                    data["gas", "flux_contff"])
        
        ds.add_field(
            ("gas", "flux_total"),
            function = _flux_total_field,
            sampling_type = "cell",
            units = "",
            force_override = True
        )
        
        print("added: flux_total")
        return ds, filter_list