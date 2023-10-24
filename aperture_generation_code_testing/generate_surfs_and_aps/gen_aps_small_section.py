import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as stats
from scipy.stats.mstats import gmean
from scipy import signal
import time
import sys
sys.path.append('/home/brandon/Documents/TF10/Task 10.2.2 Delivery v2/ap_gen_upscaling/required_functions')
from req_functions import psd_data
from req_functions import solve_R
from req_functions import WavenumberGrid_flattened
from req_functions import getAandBforRoughsurf
from req_functions import makeFracSurf_updated 
from req_functions import correlationValues 
from req_functions import calculate_scaling 
from req_functions import calculate_scaling_gen_surf
from req_functions import calculate_scaling_testing
# from req_functions import RMS_COR

## this this function for when generating apertures using best fit line with varying degrees of freedom
def correlationValues(lower,upper):
    """
    Finds the correlation for each wavenumber from the PSDR of the real data
    """    
    aper=upper-lower
    # aper[aper>0.2773] = 0 #np.nan 5% added for testing against filtered data in psd
    # aper[aper>0.889] = 0 #np.nan 1%

    # aper[aper<0] = 0

    kvals_lower,abins_lower=psd_data(lower)
    kvals_upper,abins_upper=psd_data(upper)
    kvals_aper,abins_aper=psd_data(aper)
    PR_r=abins_aper/(abins_upper+abins_lower)
    
    lamda=((2*np.pi)/kvals_lower)
    wl=(1/kvals_lower)*len(lower[0,:])/10 
    
    logx=np.log(wl)
    logy=np.log(PR_r)
    
    R_values=[]
    for i in range(1,21):
        coeffs=np.polyfit(logx,logy,deg=i)
        poly=np.poly1d(coeffs)
        yfit = lambda wl: np.exp(poly(np.log(wl)))
        test_corr=np.corrcoef(yfit(wl),PR_r)[0,1]**2
        R_values.append(test_corr)
    R_values=np.asarray(R_values)
    deg=25 #np.where(R_values>=R_values[-1])[0][0]  
    
    coeffs=np.polyfit(logx,logy,deg)
    poly=np.poly1d(coeffs)
    yfit = lambda wl: np.exp(poly(np.log(wl)))
    best_fit=yfit(wl)
    plt.plot(kvals_aper, PR_r,label='PSDR')
    plt.plot(kvals_aper,yfit(wl), label='best fit')
    plt.xscale('log')
    plt.ylabel("$PSDR(k)$",fontsize=12)
    plt.xlabel("wavenumber",fontsize=12)
    plt.legend(fontsize=12)
    plt.tick_params(labelsize=12)
    # plt.savefig('large_psdr_best_fit.png', dpi=600, bbox_inches='tight')
    plt.show()
    return 1-best_fit

##path to pre processed surface data
path_to_surface_data='/home/brandon/Documents/TF10/Task 10.2.2 Delivery v2/ap_gen_upscaling/Pre_processed_data'

## Load surface data
# shift=0.16 # this is calculated during pre processsing (upper_lower_surfaces.py), fitting with pressure film
shift=0.16 
upper=np.load(path_to_surface_data+'/upper_surface_uncorrected_res_0_1mm.npy')[15:1985,15:1985]-shift
lower=np.load(path_to_surface_data+'/lower_surface_uncorrected_res_0_1mm.npy')[15:1985,15:1985]
upper=upper[400:1400,400:1400]
lower=lower[400:1400,400:1400]
aperture_data=upper-lower

## correct apertur to remove largest 1%, excluding isolated large apertures and removing negatives
# aperture_data=aperture_data.flatten()

## path to correlated arrays
path_to_corr_arrays='/home/brandon/Documents/TF10/Task 10.2.2 Delivery v2/ap_gen_upscaling/multi_swap_and_cor'
## load correlated arrays
all_B=np.load(path_to_corr_arrays+'/B_all_8N_len_256_deg3_small_section.npy') 
all_A=np.load(path_to_corr_arrays+'/A_all_8N_len_256_deg3_small_section.npy')

## add the saved arrays rather than calculating in here
R_unique=np.load(path_to_corr_arrays+'/R_unique_full_deg3_small_section.npy')
R_values=np.load(path_to_corr_arrays+'/corrected_corr_full_deg3_small_section.npy')
reordered_R=np.asarray(solve_R(R_unique,R_values)) # this only becomes relevant when number of
min_r_index=np.where(reordered_R == reordered_R.min())

# ## path to correlated arrays, testing whether it is purely just different correlation that is making the big difference
# path_to_corr_arrays='/home/brandon/Documents/TF10/Task 10.2.2 Delivery v2/ap_gen_upscaling/multi_swap_and_cor'
# ## load correlated arrays
# all_B=np.load(path_to_corr_arrays+'/B_all_8N_len_256_deg3.npy') 
# all_A=np.load(path_to_corr_arrays+'/A_all_8N_len_256_deg3.npy')

# ## add the saved arrays rather than calculating in here
# R_unique=np.load(path_to_corr_arrays+'/R_unique_full_deg3.npy')
# R_values=np.load(path_to_corr_arrays+'/corrected_corr_full_deg3.npy')
# reordered_R=np.asarray(solve_R(R_unique,R_values)) # this only becomes relevant when number of
# min_r_index=np.where(reordered_R == reordered_R.min())


## path to H and scaling values
path_to_H='/home/brandon/Documents/TF10/Task 10.2.2 Delivery v2/ap_gen_upscaling/Find_H_and_int/H_data'

## scaling
int_lower_params=pd.read_csv(path_to_H+'/int_lower_all_1_params_small_section.csv')
int_upper_params=pd.read_csv(path_to_H+'/int_upper_all_1_params_small_section.csv')
int_lower_params=int_lower_params['Values'].tolist()
int_upper_params=int_upper_params['Values'].tolist()
# med
# int_lower_params=int_lower_params[0:1][0]
# int_upper_params=int_upper_params[0:1][0]
# int_lower_params=0.0355
# int_upper_params=0.0396
# 25th perc
# int_lower_params=int_lower_params[1:2][0]
# int_upper_params=int_upper_params[1:2][0]
# 75th perc
int_lower_params=int_lower_params[2:3][0]
int_upper_params=int_upper_params[2:3][0]
# int_lower_params=0.041
# int_upper_params=0.04
## bottom cap
# int_lower_params_bottom_cap=int_lower_params[3:4][0]
# int_upper_params_bottom_cap=int_upper_params[3:4][0]
## top cap
# int_lower_params=int_lower_params[4:5][0]
# int_upper_params=int_upper_params[4:5][0]
## add required values
# int_lower_params=0.03 #int_lower_params[3:4][0]
# int_upper_params=0.03 #int_upper_params[3:4][0]

## Hurst
H_lower_params=pd.read_csv(path_to_H+'/H_lower_all_1_params_small_section.csv')
H_upper_params=pd.read_csv(path_to_H+'/H_upper_all_1_params_small_section.csv')
H_lower_params=H_lower_params['Values'].tolist()
H_upper_params=H_upper_params['Values'].tolist()
## med
# H_lower_params=H_lower_params[0:1][0]
# H_upper_params=H_upper_params[0:1][0]
# H_lower_params=0.85
# H_upper_params=0.851
# # 25th perc
# H_lower_params=H_lower_params[1:2][0]
# H_upper_params=H_upper_params[1:2][0]
# 75th perc
H_lower_params=H_lower_params[2:3][0]
H_upper_params=H_upper_params[2:3][0]
# H_lower_params=0.892
# H_upper_params=0.899
## bottom cap
# H_lower_params_bottom_cap=H_lower_params[3:4][0]
# H_upper_params_bottom_cap=H_upper_params[3:4][0]
## top cap
# H_lower_params_top_cap=H_lower_params[4:5][0]
# H_upper_params_top_cap=H_upper_params[4:5][0]
## add required values
# H_lower_params=0.95
# H_upper_params=0.95

## anisotropy
aniso_lower_params=pd.read_csv(path_to_H+'/aniso_lower_params_1_small_section.csv') ## made a mistake and this data has not been saved
aniso_upper_params=pd.read_csv(path_to_H+'/aniso_upper_params_1_small_section.csv')
aniso_lower_params=aniso_lower_params['Values'].tolist()
aniso_upper_params=aniso_upper_params['Values'].tolist()
## med
aniso_lower_params=aniso_lower_params[0:1][0]
aniso_upper_params=aniso_upper_params[0:1][0]
# 25th perc
# aniso_lower_params=aniso_lower_params[1:2][0]
# aniso_upper_params=aniso_upper_params[1:2][0]
# 75th perc
# aniso_lower_params=aniso_lower_params[2:3][0]
# aniso_upper_params=aniso_upper_params[2:3][0]



## Generate surfaces and apertures
# seed=np.arange(1,101,1) #iterater
seed=np.arange(1,101,1) #iterater
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
     
             
np.save('aps_de3_small_section_nc.npy',all_aps)
np.save('up_surf_de3_small_section_nc.npy',up_surf)
np.save('low_surf_de3_small_section_nc.npy',low_surf)

# np.save('aps_de3_small_section_cor_check.npy',all_aps)
# np.save('up_surf_de3_small_section_cor_check.npy',up_surf)
# np.save('low_surf_de3_small_section_cor_check.npy',low_surf)
# np.save('shift1_75th_input_deg25.npy',shift1)


calculate_scaling_gen_surf = True
if calculate_scaling_gen_surf == True:
    ## run this with correct testing, test_intercept output should be 1, showing the same resolution between real and gen, 
    ## and test output should be the desired resolution, same as target rms
    upper_surf_for_testing=up_surf[0].reshape(513,513)
    test_intercept=calculate_scaling_testing(upper_surf_for_testing)
    test=target_rms/test_intercept 

calculate_scaling = False
if calculate_scaling == True:
    ## This test needs to be run with scaling for testing, checks rescaling is working correctly assuming same resolution
    check_for_correct_scaling=calculate_scaling(upper_surf_for_testing) # should be the same as target_rms


# from req_functions import RMS_COR
# x_length=len(up_surf[1,:])
# y_length=len(up_surf[:,1])   
# # H and intercept
# H_x=[]
# int_x=[]
# H_y=[]
# int_y=[]

# for i in range(0,y_length):
#     t=surface_gen[i,:].reshape(x_length,1)
#     h,intercept,std=RMS_COR(t)
#     H_x.append(h)
#     intercept=intercept*(len(surface_real)**2/len(surface_gen)**2)
#     int_x.append(intercept) 
# # y

# t=up_surf[0].reshape(513,1)
# h,intercept,std=RMS_COR(t)
# # H_x.append(h)
# intercept=intercept*(len(upper)**2/len(up_surf)**2)
# # int_x.append(intercept) 


