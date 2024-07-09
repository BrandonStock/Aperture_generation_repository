#!/bin/bash

## Set working directory
cd "$(dirname "${BASH_SOURCE[0]}")"

# Calculate Hurst exponent and scaling parameter
run_calculate_H=true
if [ "$run_calculate_H" = true ]; then
    ## Check which calculation method should be run based on the configuration file
    selected_method=$(python3 -c "
config = {}
with open('input_variables.py') as f:
    exec(f.read(), config)
methods = ['calc_H_PSD_filtered', 'calc_H_PSD_unfiltered', 'calc_H_PSD_full_surf_filtered', 'calc_H_PSD_full_surf_unfiltered', 'calc_H_RMS_COR','calc_H_RMS_COR_comb_profs']
selected = None
for method in methods:
    if config.get(method, False):
        selected = method
        break
print(selected)
")
    
    if [ "$selected_method" ]; then
        method_script="${selected_method}.py"
        cd aperture_generation/Find_H_and_int
        start=$(date +%s)
        python3 "$method_script"
        end=$(date +%s)
        cd ../..
        echo Finished running $method_script in $(($end - $start)) seconds.
    else
        echo No valid calculation method is set to True in input_variables.py
    fi
else
    echo Run_calculate_H is false
fi

# Number swapping algorithm for correlation
run_number_swap=true
if [ "$run_number_swap" = true ]; then
    cd aperture_generation/multi_swap_and_cor
    # Run number_swap_corr_min_to_max_PSDR.py
    start=$(date +%s)
    python3 number_swap_corr_min_to_max_PSDR.pyS
    end=$(date +%s)
    cd ../..
    echo Finished running number_swap_corr_min_to_max_PSDR in $(($end - $start)) seconds.
else
    echo Run_number_swap is false
fi

# Generate apertures
run_gen_aps=true
if [ "$run_gen_aps" = true ]; then
    cd aperture_generation/generate_surfs_and_aps # For testing
    # Run gen_aps.py
    start=$(date +%s)
    python3 gen_aps.py
    end=$(date +%s)
    cd ../..
    echo Finished running gen_aps in $(($end - $start)) seconds.
else
    echo Run_gen_aps is false
fi

