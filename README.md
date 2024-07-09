# aperture_generation

A python code for generating synthetic rock fracture apertures using self-affine conepts and information obtained from upper and lower surface scans of natural rock fractures

***
## To run code
Pre processed surfaces must be named upper.npy and lower.npy (.npy file format) and placed in /aperture_generation/Pre_processed_data directory.

input_variables.py file is used to define number of realisations, resolution, number of processors and the method for extracting Hurst exponent.

run_ap_gen.sh executes all the required scripts. Turn options to false if not required, for example number_swapping only needs to be run once for one set of input fracture surfaces.

generated apertures and surfaces are found in /aperture_generation/generate_surfs_and_aps directory

