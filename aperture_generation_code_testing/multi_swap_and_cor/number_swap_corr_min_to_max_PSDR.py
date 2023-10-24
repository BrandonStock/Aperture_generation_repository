import numpy as np
import time
import random
import scipy.stats as stats
import multiprocessing as mp
import os
from scipy.stats import mielke
# from sklearn.neighbors import KernelDensity
import matplotlib.pyplot as plt
import sys
sys.path.append('/home/brandon/Documents/TF10/Task 10.2.2 Delivery v2/ap_gen_upscaling/required_functions')
# sys.path.append('./ap_gen_upscaling/required_functions')
from req_functions import correlationValues
from req_functions import psd_data
from req_functions import solve_R

##path to pre processed surface data
path_to_surface_data='/home/brandon/Documents/TF10/Task 10.2.2 Delivery v2/ap_gen_upscaling/Pre_processed_data'

## Load surface data
# shift=0.16 # this is calculated during pre processsing (upper_lower_surfaces.py), fitting with pressure film
shift=0.16 #-1.7611702612734348 # shifted so no minimum
upper=np.load(path_to_surface_data+'/upper_surface_uncorrected_res_0_1mm.npy')[15:1985,15:1985]-shift
lower=np.load(path_to_surface_data+'/lower_surface_uncorrected_res_0_1mm.npy')[15:1985,15:1985]
upper=upper[400:1400,400:1400]
lower=lower[400:1400,400:1400]

# upper=upper[0:1000,0:1000]
# lower=lower[0:1000,0:1000]

aper=upper-lower
aper[aper<0] = 0 #np.nan

# ## % cropped aperture
aperture5=aper
# aperture5[aperture5>0.2773] = 0 #np.nan 5%
# aperture5[aperture5>0.297] = 0 #np.nan #4%
# aperture5[aperture5>0.34] = 0 #np.nan 3%
# aperture5[aperture5>0.48] = 0 #np.nan 2%
# aperture5[aperture5>0.889] = 0 #np.nan 1%


## get PSD values
kvals_lower,abins_lower=psd_data(lower)
kvals_upper,abins_upper=psd_data(upper)
kvals_aper,abins_aper=psd_data(aperture5)
PR_r=abins_aper/(abins_upper+abins_lower)

lamda=(1/kvals_lower)*int(len(upper)/2)-1# 984 #value=half length -1

plt.plot(kvals_lower,abins_lower)
plt.plot(kvals_lower,abins_upper)
plt.plot(kvals_lower,abins_aper)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('k')
plt.ylabel('PSD')
# plt.savefig('PSD_filtered_ap_5%.png', dpi=600, bbox_inches='tight')
plt.show()

ind_max=np.where(PR_r==PR_r[0:-1].max())
test_line=np.linspace(0,PR_r[ind_max[0][0]],ind_max[0][0])

plt.plot(kvals_lower[0:-1],PR_r[0:-1])
plt.plot(kvals_lower[0:ind_max[0][0]],test_line)

plt.xscale('log')
plt.xlabel('k')
plt.ylabel('PSDR(k)')
# plt.yscale('log')
# plt.savefig('PSDR_filtered_ap_5%.png', dpi=600, bbox_inches='tight')
plt.show()


## doing matching from min to max
lamda=((2*np.pi)/kvals_lower)
wl=(1/kvals_lower)*len(lower[0,:])/10 

logx=np.log(wl)
logy=np.log(PR_r)
min_r_index=np.where(PR_r == PR_r.min())[0][0]
max_r_index=np.where(PR_r == PR_r[0:-1].max())[0][0]

PR_new=PR_r[0:max_r_index+1]
wl=wl[0:max_r_index+1]
logx=np.log(wl)
logy=np.log(PR_new)

deg=4
coeffs=np.polyfit(logx,logy,deg=deg)
poly=np.poly1d(coeffs)
yfit = lambda wl: np.exp(poly(np.log(wl)))
# test_corr=np.corrcoef(yfit(wl),PR_r)[0,1]**2
# R_values.append(test_corr)
# R_values=np.asarray(R_values)
#np.where(R_values>=R_values[-1])[0][0]  ## Change 0.93 to different values for different cases to get best fit

coeffs=np.polyfit(logx,logy,deg)
poly=np.poly1d(coeffs)
yfit = lambda wl: np.exp(poly(np.log(wl)))

plt.plot(kvals_aper[0:max_r_index+1],PR_new)
plt.plot(kvals_aper[0:max_r_index+1],yfit(wl), label='polyfit O(3)')
plt.xscale('log')
plt.show()

polyfit=yfit(wl)
# test=np.round(wl,7)
# R_max_to_min_2=corr_deg_2[0:min_r_index_2[0][0]]
line_2=np.full(int(len(upper)/2-len(polyfit)),polyfit[-1])
corrected_corr=np.append(polyfit,line_2)

# plt.plot(kvals_aper,PR_r)
# plt.plot(kvals_aper,corrected_corr)
# plt.xscale('log')
# plt.show()

import scipy.interpolate as interp
arr1_interp = interp.interp1d(np.arange(corrected_corr.size),corrected_corr)
corrected_corr = arr1_interp(np.linspace(0,corrected_corr.size-1,256))

# ## change to corr
corrected_corr=1-corrected_corr
# # plt.plot(kvals_aper,PR_r)
plt.plot(corrected_corr)
plt.xscale('log')
plt.show()

R_unique=np.unique(corrected_corr)

np.save('R_unique_full_deg3_small_section.npy',R_unique)
np.save('corrected_corr_full_deg3_small_section.npy',corrected_corr)
# 

# Input parameters
R=R_unique
N=8
num_proc=9
seed1=1
seed2=2 
n=((2*2**N+1)**2)

## input random number sequences A and B
np.random.seed(seed1)
A=np.random.randn(n)*np.pi*2 ## change to randn for normal distribution 

np.random.seed(seed2)
B=np.random.randn(n)*np.pi*2

all_B=[]
all_A=[]
for k in range(len(R)):
    print (k)
        
    # these values can be customised
    R=R_unique[k]
    # R=0.4
    maxIters=100000
    minSwaps=1000
    tol=0.001
    
    # length of the array
    n=((2*2**N+1)**2)
    n = n + (n%num_proc) # always keep this, makes sure we can reshape array to divide it over processes
    
    # function that checks whether to swap based on chosen pairs, withoot calculating correlation over full array256
    def swap_check(A,B,t): # A and B are the pairs to be potentially swapped, t is target direction
        if A[0] == A[1] or B[0] == B[1]:
            return False
        before = stats.pearsonr(A,B)[0]#.real #np.mean(A/B) #
        B_swap = np.array([B[1],B[0]])
        after =stats.pearsonr(A,B_swap)[0]#.real # np.mean(A/B_swap) #
        return np.sign(after-before) == t
    
    def parallel_swap(A,B):
        #print('process', os.getpid(), 'starting...')
        corr=stats.pearsonr(A,B)[0]
        target_direction=np.sign(R-corr)
    
        iterations=0
        swaps=0
        l = len(A)
    
        mod=1000 #determines how often correlation is evaluated (adapts automatically if necessary) shhould be in proportion to array length
    
        # swaps until value within tolerance of target correlation or maximum iteration number reached
        # while iterations < maxIters and not np.isclose(R,corr,rtol=0.05):
        while not np.isclose(corr,R,rtol=tol):
            # randomly select 2 position to swap
            p1 = random.randint(0,l-1)
            p2 = random.randint(0,l-1)
            A_swap = np.array([A[p1],A[p2]])
            B_swap = np.array([B[p1],B[p2]])
            
            # swap if correlation is improved
            if  swap_check(A_swap,B_swap,target_direction):
                B[p1] = B_swap[1]
                B[p2] = B_swap[0]
                # corr = stats.pearsonr(A,B)[0]
                if swaps%mod == 0:
                    # print (mod)
                    corr_new=stats.pearsonr(A,B)[0]#.real # original
                    # corr_new=stats.pearsonr(A.real,B.real)[0] # added for correlating the fft phase
                    # if getting close to target too fast reduce 'mod'
                    if abs(R-corr)/abs(corr-corr_new) < 2:
                        mod = int(mod/2)
                        if (mod == 0): mod = 1 ## add corrections so mod can never get to zero (needs doing)
                    corr = corr_new
                    target_direction = np.sign(R-corr)
                    # print (corr)
                swaps += 1
                
            iterations += 1
            if iterations >= maxIters and swaps >= minSwaps:
                break
        
        #print('process', os.getpid(), 'finished!\nIn', iterations, 'iterations, with correlation', corr, ', through', swaps, 'swaps')
        return np.array([np.array(A),np.array(B)])
    
    if __name__ == '__main__' :
    
        #print('target correlation:', R)
    
        # for testing
        t1 = time.time()
    
        corr=stats.pearsonr(A,B)[0]
        #print('initial correlation:', corr)
    
        with mp.Pool(processes = num_proc) as pool:
            
            while True:
                # prepare arrays for multiprocessing
                A = np.reshape(A,[num_proc,-1])
                B = np.reshape(B,[num_proc,-1])
                # collect results in C
                C = np.array([])
                C = np.array(pool.starmap(parallel_swap,zip(A,B)))
                # reshape C back into A and B in correct order
                D = C[0]
                for i in range(1,num_proc): ## was(1,4)
                    D = np.concatenate((D,C[i]), axis = 1)
                A = D[0]
                B = D[1]
                # print(A)
                # print(B)
                corr=stats.pearsonr(A,B)[0]                
                # print('current correlation:', corr)
                # check for completion criterium
                if np.isclose(corr,R,rtol=tol):
                    break
                # otherwise shuffle arrays and start over (so that each process gets a different selection of the array)
                else:
                            
                    S = np.arange(n)
                    np.random.shuffle(S)
                    A = A[S]
                    B = B[S]
    
        t2 = time.time()
    
        print('time (s):', t2-t1)
        print ('target correlation', R)
        print('final correlation:', corr)
        all_B.append(B)
        all_A.append(A)
        
np.save('B_all_8N_len_256_deg3_small_section.npy',all_B)
np.save('A_all_8N_len_256_deg3_small_section.npy', all_A) 





