#!/bin/bash

#import os
#path_to_H = os.path.join(os.path.dirname(__file__), 'aperture_generation/Find_H_and_int')
cd "$(dirname "${BASH_SOURCE[0]}")"

# path to calculate H
cd aperture_generation/Find_H_and_int

# run calculate_H.py
start=`date +%s`
python3 calculate_H.py
end=`date +%s`
echo Finished running calculate_H in `expr $end - $start` seconds.

# path to swapping algorithm
cd ../multi_swap_and_cor

# run number_swap_corr_min_to_max_PSDR.py
start=`date +%s`
python3 number_swap_corr_min_to_max_PSDR.py
end=`date +%s`
echo Finished running number_swap_corr_min_to_max_PSDR in `expr $end - $start` seconds.

# path to generate aps
cd ../generate_surfs_and_aps

# run number_swap_corr_min_to_max_PSDR.py
start=`date +%s`
python3 gen_aps.py
end=`date +%s`
echo Finished running gen_aps in `expr $end - $start` seconds.
