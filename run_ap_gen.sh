#!/bin/bash

## set working directory
cd "$(dirname "${BASH_SOURCE[0]}")"

# calculate Hurst exponent and scaling parameter
run_calculate_H=false
if $run_calculate_H; then
	cd aperture_generation/Find_H_and_int
	# run calculate_H.py
	start=`date +%s`
	python3 calculate_H.py
	end=`date +%s`
	cd ../..
	echo Finished running calculate_H in `expr $end - $start` seconds.
else
	echo Run_calculate_H is false
fi

# number swapping algorithm for correlation
run_number_swap=false
if $run_number_swap; then
	cd aperture_generation/multi_swap_and_cor
	# run number_swap_corr_min_to_max_PSDR.py
	start=`date +%s`
	python3 number_swap_corr_min_to_max_PSDR.py
	end=`date +%s`
	cd ../..
	echo Finished running number_swap_corr_min_to_max_PSDR in `expr $end - $start` seconds.
else
	echo Run_number_swap is false
fi

# generate apertures
run_gen_aps=true
if $run_gen_aps; then
	cd aperture_generation/generate_surfs_and_aps # for testing
	# run gen_aps.py
	start=`date +%s`
	python3 gen_aps.py
	end=`date +%s`
	cd ../..
	echo Finished running gen_aps in `expr $end - $start` seconds.
else
	echo Run_gen_aps is false
fi
