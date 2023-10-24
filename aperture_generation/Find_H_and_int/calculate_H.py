import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import iqr
import pandas as pd
import os
import sys
# sys.path.append('/home/brandon/Documents/TF10/Task 10.2.2 Delivery v2/ap_gen_upscaling/required_functions')
func_path = os.path.join(os.path.dirname(__file__), '../required_functions')
sys.path.append(func_path)
from req_functions import RMS_COR
    
##path to pre processed surface data
path_to_surface_data='../Pre_processed_data'
## Load surface data
upper=np.load(path_to_surface_data+'/upper.npy')
lower=np.load(path_to_surface_data+'/lower.npy')

## calculate H for parallel traces
find_H_parallel=True
if find_H_parallel == True:
        
    x_length=len(upper[1,:])
    y_length=len(upper[:,1])
    
    ## Lower H and intercept
    H_lower_x=[]
    int_lower_x=[]
    H_lower_y=[]
    int_lower_y=[]
    
    H_upper_x=[]
    int_upper_x=[]
    H_upper_y=[]
    int_upper_y=[]
    
    # x
    for i in range(0,y_length):
        # print (i)
        t=lower[i,:].reshape(x_length,1)
        h,intercept,std=RMS_COR(t)
        H_lower_x.append(h)
        int_lower_x.append(intercept) 
    # y
    for i in range(0,x_length):
        t=lower[:,i].reshape(y_length,1)
        h,intercept,std=RMS_COR(t)
        H_lower_y.append(h)
        int_lower_y.append(intercept)
        
    ## upper H and intercept
    # x
    for i in range(0,y_length):
        # print (i)
        t=upper[i,:].reshape(x_length,1)
        h,intercept,std=RMS_COR(t)
        H_upper_x.append(h)
        int_upper_x.append(intercept) 
    # y
    for i in range(0,x_length):
        t=upper[:,i].reshape(y_length,1)
        h,intercept,std=RMS_COR(t)
        H_upper_y.append(h)
        int_upper_y.append(intercept)
    
    np.save('H_data/H_lower_x.npy',H_lower_x)
    np.save('H_data/H_lower_y.npy',H_lower_y)
    np.save('H_data/H_upper_x.npy',H_upper_x)
    np.save('H_data/H_upper_y.npy',H_upper_y)

    np.save('H_data/int_lower_x.npy',int_lower_x)
    np.save('H_data/int_lower_y.npy',int_lower_y)
    np.save('H_data/int_upper_x.npy',int_upper_x)
    np.save('H_data/int_upper_y.npy',int_upper_y)
    
    ## Lower H and scaling
    Hs_lower_all_array=H_lower_x+H_lower_y
    median=np.median(Hs_lower_all_array)
    pc25=np.percentile(Hs_lower_all_array,25)
    pc75=np.percentile(Hs_lower_all_array,75)
    d={'Index_Title':['median','pc25','pc75'],'Values':[median,pc25,pc75]}
    df = pd.DataFrame(d).set_index('Index_Title')
    df.to_csv('H_data/H_lower_all.csv')

    int_lower_all_array=int_lower_x+int_lower_y
    median=np.median(int_lower_all_array)
    pc25=np.percentile(int_lower_all_array,25)
    pc75=np.percentile(int_lower_all_array,75)
    d={'Index_Title':['median','pc25','pc75'],'Values':[median,pc25,pc75]}
    df = pd.DataFrame(d).set_index('Index_Title')
    df.to_csv('H_data/int_lower_all.csv')

    
    ## Upper H and scaling
    Hs_upper_all_array=H_upper_x+H_upper_y
    median=np.median(Hs_upper_all_array)
    pc25=np.percentile(Hs_upper_all_array,25)
    pc75=np.percentile(Hs_upper_all_array,75)
    d={'Index_Title':['median','pc25','pc75'],'Values':[median,pc25,pc75]}
    df = pd.DataFrame(d).set_index('Index_Title')
    df.to_csv('H_data/H_upper_all.csv')

    int_upper_all_array=int_upper_x+int_upper_y
    median=np.median(int_upper_all_array)
    pc25=np.percentile(int_upper_all_array,25)
    pc75=np.percentile(int_upper_all_array,75)
    d={'Index_Title':['median','pc25','pc75'],'Values':[median,pc25,pc75]}
    df = pd.DataFrame(d).set_index('Index_Title')
    df.to_csv('H_data/int_upper_all.csv')
    
    
    ## get anisotropy data
    ## lower surface
    lower_aniso=np.asarray(H_lower_x)/np.asarray(H_lower_y)   
    np.save('H_data/lower_aniso.npy',lower_aniso)
    
    median=np.median(lower_aniso)
    pc25=np.percentile(lower_aniso,25)
    pc75=np.percentile(lower_aniso,75)
    d={'Index_Title':['median','pc25','pc75'],'Values':[median,pc25,pc75]}
    df = pd.DataFrame(d).set_index('Index_Title')
    df.to_csv('H_data/aniso_lower_all.csv')


    ## upper surface
    upper_aniso=np.asarray(H_upper_x)/np.asarray(H_upper_y)   
    np.save('H_data/upper_aniso.npy',lower_aniso)

    median=np.median(upper_aniso)
    pc25=np.percentile(upper_aniso,25)
    pc75=np.percentile(upper_aniso,75)
    d={'Index_Title':['median','pc25','pc75'],'Values':[median,pc25,pc75]}
    df = pd.DataFrame(d).set_index('Index_Title')
    df.to_csv('H_data/aniso_upper_all.csv')
