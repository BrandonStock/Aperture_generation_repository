import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import iqr
import pandas as pd
import os
import sys
func_path = os.path.join(os.path.dirname(__file__), '../required_functions')
sys.path.append(func_path)
from req_functions import RMS_COR
    
##path to pre processed surface data
path_to_surface_data='../Pre_processed_data'
## Load surface data
upper=np.load(path_to_surface_data+'/upper.npy')
lower=np.load(path_to_surface_data+'/lower.npy')

## Calculate H for parallel traces
find_H_parallel = True
if find_H_parallel:
    x_length = upper.shape[1]
    y_length = upper.shape[0]

    ## Initialize dictionaries to store results
    results = {
        "H_lower_x": [],
        "int_lower_x": [],
        "H_lower_y": [],
        "int_lower_y": [],
        "H_upper_x": [],
        "int_upper_x": [],
        "H_upper_y": [],
        "int_upper_y": []
    }

    ## Function to process rows and columns
    def process_traces(data, prefix):
        for i in range(y_length):
            t = data[i, :].reshape(x_length, 1)
            h, intercept, std = RMS_COR(t)
            results[f"H_{prefix}_x"].append(h)
            results[f"int_{prefix}_x"].append(intercept)

        for i in range(x_length):
            t = data[:, i].reshape(y_length, 1)
            h, intercept, std = RMS_COR(t)
            results[f"H_{prefix}_y"].append(h)
            results[f"int_{prefix}_y"].append(intercept)

    ## Process lower and upper data
    process_traces(lower, "lower")
    process_traces(upper, "upper")

    ## Save results to disk
    for key, value in results.items():
        np.save(f'calc_H_RMS_COR/{key}.npy', value)
        
    
    # ## Lower H and scaling
    Hs_lower_all_array=results.get("H_lower_x")+results.get("H_lower_y")
    Hs_lower_all_array_pc75=np.percentile(Hs_lower_all_array,75)
    np.save('calc_H_RMS_COR/H_lower.npy',Hs_lower_all_array_pc75)

    int_lower_all_array=results.get("int_lower_x")+results.get("int_lower_y")
    int_lower_all_array_pc75=np.percentile(int_lower_all_array,75)
    np.save('calc_H_RMS_COR/Sp_lower.npy',int_lower_all_array_pc75)
    
    # ## Upper H and scaling
    Hs_upper_all_array=results.get("H_upper_x")+results.get("H_upper_y")
    Hs_upper_all_array_pc75=np.percentile(Hs_upper_all_array,75)
    np.save('calc_H_RMS_COR/H_upper.npy',Hs_upper_all_array_pc75)
    
    int_upper_all_array=results.get("int_upper_x")+results.get("int_upper_y")
    int_upper_all_array_pc75=np.percentile(int_upper_all_array,75)
    np.save('calc_H_RMS_COR/Sp_upper.npy',int_upper_all_array_pc75)

    # ## get anisotropy data
    # ## lower surface
    lower_aniso=np.asarray(results.get("H_lower_x"))/np.asarray(results.get("H_lower_y"))
    lower_aniso_pc75=np.percentile(lower_aniso,75)
    np.save('calc_H_RMS_COR/lower_aniso_pc75.npy',lower_aniso_pc75)

    # ## upper surface
    upper_aniso=np.asarray(results.get("H_upper_x"))/np.asarray(results.get("H_upper_y"))
    upper_aniso_pc75=np.percentile(upper_aniso,75)
    np.save('calc_H_RMS_COR/upper_aniso_pc75.npy',upper_aniso_pc75)

