import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as stats
from scipy.stats.mstats import gmean
from scipy import signal
import time
import os
import sys
func_path = os.path.join(os.path.dirname(__file__), '../required_functions')
sys.path.append(func_path)
from req_functions import psd_data
from req_functions import solve_R
from req_functions import WavenumberGrid_flattened
from req_functions import getAandBforRoughsurf
from req_functions import makeFracSurf_updated 
from req_functions import correlationValues 
from req_functions import calculate_scaling 
from req_functions import calculate_scaling_gen_surf
from req_functions import calculate_scaling_testing


path_to_surface_data='../Pre_processed_data'

## Load surface data
upper=np.load(path_to_surface_data+'/upper.npy')
lower=np.load(path_to_surface_data+'/lower.npy')
aperture_data=upper-lower

## path to correlated arrays
path_to_corr_arrays='../multi_swap_and_cor'
## load correlated arrays
all_B=np.load(path_to_corr_arrays+'/B.npy') 
all_A=np.load(path_to_corr_arrays+'/A.npy')

## add the saved arrays rather than calculating in here
R_unique=np.load(path_to_corr_arrays+'/R_unique.npy')
R_values=np.load(path_to_corr_arrays+'/corrected_corr.npy')
reordered_R=np.asarray(solve_R(R_unique,R_values)) 
min_r_index=np.where(reordered_R == reordered_R.min())

## path to H and scaling values
path_to_H='../Find_H_and_int/H_data'

## scaling
int_lower_params=pd.read_csv(path_to_H+'/int_lower_all.csv')
int_upper_params=pd.read_csv(path_to_H+'/int_upper_all.csv')
int_lower_params=int_lower_params['Values'].tolist()
int_upper_params=int_upper_params['Values'].tolist()
# 75th perc
int_lower_params=int_lower_params[2:3][0]
int_upper_params=int_upper_params[2:3][0]

## Hurst
H_lower_params=pd.read_csv(path_to_H+'/H_lower_all.csv')
H_upper_params=pd.read_csv(path_to_H+'/H_upper_all.csv')
H_lower_params=H_lower_params['Values'].tolist()
H_upper_params=H_upper_params['Values'].tolist()
# 75th perc
H_lower_params=H_lower_params[2:3][0]
H_upper_params=H_upper_params[2:3][0]

## anisotropy
aniso_lower_params=pd.read_csv(path_to_H+'/aniso_lower_all.csv') 
aniso_upper_params=pd.read_csv(path_to_H+'/aniso_upper_all.csv')
aniso_lower_params=aniso_lower_params['Values'].tolist()
aniso_upper_params=aniso_upper_params['Values'].tolist()
## med
aniso_lower_params=aniso_lower_params[0:1][0]
aniso_upper_params=aniso_upper_params[0:1][0]

## Generate surfaces and apertures
seed=np.arange(1,3,1) #iterater
all_aps=[]
up_surf=[]
low_surf=[]
shift1=[]

for i in range(0,len(seed)):
      N=8
      seed_for_rearrange=seed[i]
      print (i)
      S = np.arange((2*2**N+1)**2) 
      # S = np.arange(10000) 
      np.random.seed(seed_for_rearrange)
      np.random.shuffle(S)
      # print(S)
      all_B=all_B[0::,S]
      all_A=all_A[0::,S]
      #vairbales for K grid generation and subsquent size of generated surface (number of cells)
      nval=2*2**N+1
      K_cutoff=int(min_r_index[0])-1 # plus 1 added, remove to get to original
      cutoff_length=int(nval/2)+1
      ## calculate wavenumbers and put into flattened array for correct position
      K,K_c=WavenumberGrid_flattened(nval=nval,cutoff_length=cutoff_length,R_values=R_unique,K_cutoff=K_cutoff)
      K_c=K_c.reshape(nval,nval)
      # calculated A and B values, which in this case are phase1 and phase2 over the full scale
      N=nval
      phase1,phase2=getAandBforRoughsurf(reordered_R=reordered_R,R_unique=R_unique,K=K,all_A=all_A,all_B=all_B)
      phase1=phase1.reshape(N,N)
      phase2=phase2.reshape(N,N)
      average_intercept_upper=int_upper_params
      average_intercept_lower=int_lower_params
     
      #generate upper surface ##########################################
      aniso_upper = aniso_upper_params 
      H_upper = H_upper_params 
      upper_surf=makeFracSurf_updated(N=N,H=H_upper,anisotropy=aniso_upper,phase1=phase1,wavenumber_grid=K_c) #phase2=phase2
      # scaling parameters
      target_rms=average_intercept_upper 
      rms_heights_upper_surf=calculate_scaling_gen_surf(upper,upper_surf) ## correct
      # print (average_intercept_upper, rms_heights_upper_surf)
      # rms_heights_upper_surf=calculate_scaling(upper_surf) ## for testing
      upper_surf *= target_rms/rms_heights_upper_surf
      up_surf.append(upper_surf.flatten())

      #generate lower surface #########################################
      aniso_lower =aniso_lower_params
      H_lower=H_lower_params 
      # make surface
      lower_surf=makeFracSurf_updated(N=N,H=H_lower,anisotropy=aniso_lower,phase1=phase2,wavenumber_grid=K_c)
      target_rms=average_intercept_lower     
      rms_heights_lower_surf=calculate_scaling_gen_surf(lower,lower_surf) ## correct
      # print (average_intercept_lower, rms_heights_lower_surf)
      # rms_heights_lower_surf=calculate_scaling(lower_surf) ## for testing
      lower_surf *= target_rms/rms_heights_lower_surf
      low_surf.append(lower_surf.flatten())
      aperture=upper_surf-lower_surf
      ## shift to remove zero    
      shift=np.min(aperture)
      shift1.append(shift)
      aperture=upper_surf-(lower_surf+shift)
     
      # # #test for shfting to same contact area
      shift_to_same_contact = False
      if shift_to_same_contact == True:
          n_zero=np.where(aperture<=0)
          n_zero=n_zero[0].size
          percent_contact_gen=(n_zero/aperture.size)*100
          print ('before, gen',percent_contact_gen)
          n_zero=np.where(aperture_data<=0)
          n_zero=n_zero[0].size
          percent_contact_data=(n_zero/aperture_data.size)*100
          print ('before,data',percent_contact_data)
          shift=np.linspace(0.1,-1.3,1201)
          # A=np.zeros((N,N))
          ##A=np.zeros((2**N*2+1,2**N*2+1))
          for i in range(0,len(shift)):
              # print (shift[i])
              aperture=upper_surf-lower_surf-shift[i]
              # print (np.min(apX))
              n_zero=np.where(aperture<=0)
              n_zero=n_zero[0].size
              percent_con_gen=(n_zero/aperture.size)*100
              # print (percent_con_gen)
              if np.isclose(percent_con_gen,percent_contact_data,rtol=0.05):
                  aperture=aperture
                  break
              # aperture=A
          n_zero=np.where(aperture<=0)
          n_zero=n_zero[0].size
          p=(n_zero/aperture.size)*100
          print ('after, gen',p)

      all_aps.append(aperture.flatten())
             
np.save('gen_apertures.npy',all_aps)
np.save('gen_upper_surf.npy',up_surf)
np.save('gen_lower_surf.npy',low_surf)


# calculate_scaling_gen_surf = True
# if calculate_scaling_gen_surf == True:
#     ## run this with correct testing, test_intercept output should be 1, showing the same resolution between real and gen, 
#     ## and test output should be the desired resolution, same as target rms
#     upper_surf_for_testing=up_surf[0].reshape(513,513)
#     test_intercept=calculate_scaling_testing(upper_surf_for_testing)
#     test=target_rms/test_intercept 

# calculate_scaling = False
# if calculate_scaling == True:
#     ## This test needs to be run with scaling for testing, checks rescaling is working correctly assuming same resolution
#     check_for_correct_scaling=calculate_scaling(upper_surf_for_testing) # should be the same as target_rms

