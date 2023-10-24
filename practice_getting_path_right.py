import numpy as np
import os
# dirname = os.path.dirname(__file__)
func_path = os.path.join(os.path.dirname(__file__), 'aperture_generation/Pre_processed_data')

# filename = os.path.join(dirname, 'relative/path/to/file/you/want')

##path to pre processed surface data
# path_to_surface_data='/home/brandon/Documents/TF10/Task 10.2.2 Delivery v2/ap_gen_upscaling/Pre_processed_data'
# path_to_surface_data='/home/brandon/Documents/TF10/Task 10.2.2 Delivery v2/Aperture_generation_repository/aperture_generation/Pre_processed_data'

## Load surface data
# shift=0.16 # this is calculated during pre processsing (upper_lower_surfaces.py), fitting with pressure film
# shift=0.16 
upper=np.load(func_path+'/upper.npy')
lower=np.load(func_path+'/lower.npy')

# np.save('upper.npy',upper)
# np.save('lower.npy',lower)


