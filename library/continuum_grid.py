# function to calculate continuum grid
def compute_continuum_grid(min_temp, max_temp, num_temp_grid, min_dense, max_dense, num_dense_grid, min_wl, max_wl, num_wl_grid, filter_wl = None, filter_output = None, save_dir = None):
    """"
    Compute continuum flux (bound-free, two-photon, free-free) over a grid of temperature and number density.
    This function will also create variables like df01, df02 in your notebook for easy access.
    """

    # build the grid for temperature, density, and wavelength
    temperature_grid = np.logspace(np.log10(min_temp), np.log10(max_temp), num_temp_grid)
    density_grid = np.logspace(np.log10(min_dense), np.log10(max_dense), num_dense_grid)
    wl = np.logspace(np.log10(min_wl), np.log10(max_wl), num_wl_grid)

    # create a continuum using Pyneb package
    C = Continuum()
    
    # check if we have multiple filters
    filter_wl_array = np.asarray(filter_wl) if filter_wl is not None else None
    filter_output_array = np.asarray(filter_output) if filter_output is not None else None
    
    # create directory for saving dataframes
    if save_dir is None:
        save_dir = "df_filter"
    os.makedirs(save_dir, exist_ok = True)
    
    all_df_results = {}
    all_interp_funcs = {}

    # logic for multiple filters
    if filter_wl_array is not None and filter_wl_array.ndim == 2:
        n_filters = filter_wl_array.shape[0]
        print(f"processing {n_filters} filters...")
        
        for filter_idx in range(n_filters):
            filter_number = filter_idx + 1 
            current_filter_wl = filter_wl_array[filter_idx]
            current_filter_output = filter_output_array[filter_idx]
            
            mask = (current_filter_wl >= min_wl) & (current_filter_wl <= max_wl) & (current_filter_output > 0)
            calc_filter_wl = current_filter_wl[mask]
            calc_filter_output = current_filter_output[mask]
            
            if len(calc_filter_wl) == 0:
                df_results = pd.DataFrame()
            else:
                filter_min = calc_filter_wl.min()
                filter_max = calc_filter_wl.max()
                results = []
                
                for tem in temperature_grid:
                    for den in density_grid:
                        contH = C.get_continuum(tem = tem, den = den, wl = wl, HI_label = None, cont_HI = True, cont_HeI = True, cont_HeII = True, cont_2p = False, cont_ff = False)
                        cont2p = C.get_continuum(tem = tem, den = den, wl = wl, HI_label = None, cont_HI = False, cont_HeI = False, cont_HeII = False, cont_2p = True, cont_ff = False)
                        contff = C.get_continuum(tem = tem, den = den, wl = wl, HI_label = None, cont_HI = False, cont_HeI = False, cont_HeII = False, cont_2p = False, cont_ff = True)
                        
                        mask_wl = (wl >= max(filter_min, wl.min())) & (wl <= min(filter_max, wl.max()))
                        
                        if mask_wl.sum() == 0:
                            flux_avg_H = flux_avg_2p = flux_avg_ff = np.nan
                        else:
                            fluxH_interp = interp1d(wl[mask_wl], contH[mask_wl] * wl[mask_wl], fill_value = "extrapolate")
                            flux2p_interp = interp1d(wl[mask_wl], cont2p[mask_wl] * wl[mask_wl], fill_value = "extrapolate")
                            fluxff_interp = interp1d(wl[mask_wl], contff[mask_wl] * wl[mask_wl], fill_value = "extrapolate")
                            
                            resultH = fluxH_interp(calc_filter_wl) * calc_filter_output
                            result2p = flux2p_interp(calc_filter_wl) * calc_filter_output
                            resultff = fluxff_interp(calc_filter_wl) * calc_filter_output
                            integral_filter = np.trapz(calc_filter_output * calc_filter_wl, calc_filter_wl)
                            
                            flux_avg_H = np.trapz(resultH, calc_filter_wl) / integral_filter if integral_filter > 0 else np.nan
                            flux_avg_2p = np.trapz(result2p, calc_filter_wl) / integral_filter if integral_filter > 0 else np.nan
                            flux_avg_ff = np.trapz(resultff, calc_filter_wl) / integral_filter if integral_filter > 0 else np.nan
                        
                        results.append({
                            "Temperature": tem,
                            "Density": den,
                            "Average Specific Flux ContH": flux_avg_H,
                            "Average Specific Flux Cont2p": flux_avg_2p,
                            "Average Specific Flux Contff": flux_avg_ff
                        })
                df_results = pd.DataFrame(results)

            filter_name = f"filter_{filter_number:02d}"
            all_df_results[filter_name] = df_results
            
            # create global variable for notebook access
            var_name = f"df{filter_number:02d}"
            globals()[var_name] = df_results
            
            df_results_clean = df_results.dropna()
            if len(df_results_clean) > 0:
                points = np.column_stack((df_results_clean["Temperature"], df_results_clean["Density"]))
                all_interp_funcs[filter_name] = {
                    "contH": LinearNDInterpolator(points, df_results_clean["Average Specific Flux ContH"]),
                    "cont2p": LinearNDInterpolator(points, df_results_clean["Average Specific Flux Cont2p"]),
                    "contff": LinearNDInterpolator(points, df_results_clean["Average Specific Flux Contff"]),
                }
            
            filename = f"df_{filter_number:02d}.txt"
            df_results.to_csv(os.path.join(save_dir, filename), sep = '\t', index = False)
        
        return all_df_results, all_interp_funcs

    # logic for single filter
    else:
        filter_wl_1d = filter_wl_array.flatten()
        filter_output_1d = filter_output_array.flatten()
        
        mask = (filter_wl_1d >= min_wl) & (filter_wl_1d <= max_wl) & (filter_output_1d > 0)
        calc_filter_wl = filter_wl_1d[mask]
        calc_filter_output = filter_output_1d[mask]
        
        filter_min = calc_filter_wl.min()
        filter_max = calc_filter_wl.max()
        results = []
        
        for tem in temperature_grid:
            for den in density_grid:
                contH = C.get_continuum(tem = tem, den = den, wl = wl, HI_label = None, cont_HI = True, cont_HeI = True, cont_HeII = True, cont_2p = False, cont_ff = False)
                cont2p = C.get_continuum(tem = tem, den = den, wl = wl, HI_label = None, cont_HI = False, cont_HeI = False, cont_HeII = False, cont_2p = True, cont_ff = False)
                contff = C.get_continuum(tem = tem, den = den, wl = wl, HI_label = None, cont_HI = False, cont_HeI = False, cont_HeII = False, cont_2p = False, cont_ff = True)
                
                mask_wl = (wl >= max(filter_min, wl.min())) & (wl <= min(filter_max, wl.max()))
                fluxH_interp = interp1d(wl[mask_wl], contH[mask_wl] * wl[mask_wl], fill_value = "extrapolate")
                flux2p_interp = interp1d(wl[mask_wl], cont2p[mask_wl] * wl[mask_wl], fill_value = "extrapolate")
                fluxff_interp = interp1d(wl[mask_wl], contff[mask_wl] * wl[mask_wl], fill_value = "extrapolate")
                
                resultH = fluxH_interp(calc_filter_wl) * calc_filter_output
                result2p = flux2p_interp(calc_filter_wl) * calc_filter_output
                resultff = fluxff_interp(calc_filter_wl) * calc_filter_output
                integral_filter = np.trapz(calc_filter_output * calc_filter_wl, calc_filter_wl)
                
                flux_avg_H = np.trapz(resultH, calc_filter_wl) / integral_filter
                flux_avg_2p = np.trapz(result2p, calc_filter_wl) / integral_filter
                flux_avg_ff = np.trapz(resultff, calc_filter_wl) / integral_filter
                
                results.append({
                    "Temperature": tem,
                    "Density": den,
                    "Average Specific Flux ContH": flux_avg_H,
                    "Average Specific Flux Cont2p": flux_avg_2p,
                    "Average Specific Flux Contff": flux_avg_ff
                })
        
        df_results = pd.DataFrame(results)
        globals()["df01"] = df_results
        
        filename = f"df_single_filter.txt"
        df_results.to_csv(os.path.join(save_dir, filename), sep = '\t', index = False)
        
        points = np.column_stack((df_results["Temperature"], df_results["Density"]))
        interp_funcs = {
            "contH": LinearNDInterpolator(points, df_results["Average Specific Flux ContH"]),
            "cont2p": LinearNDInterpolator(points, df_results["Average Specific Flux Cont2p"]),
            "contff": LinearNDInterpolator(points, df_results["Average Specific Flux Contff"]),
        }
        return df_results, interp_funcs