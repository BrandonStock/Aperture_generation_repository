# aperture_generation

***

A python code for generating synthetic rock fracture aperture based on the upper and lower surface scans of natural rock fractures

***
## To run

Pre processed surfaces must be called upper and lower and saved in .npy format and placed in directory /aperture_generation/Pre_processed_data. 
Aperture generation code, key input parameters for executing number swapping algorithm and generating apertures can be altered in input_variables.py
run_ap_gen.sh executes all the required codes, can turn some key methods to false if number swapping or calculating Hurst has already been run.
number swapping is only required once for one fracture, if new surface scans are used it should be run again
generated apertures and surfaces found in /aperture_generation/generate_surfs_and_aps/.
