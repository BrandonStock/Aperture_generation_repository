import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import iqr
import pandas as pd
import sys
sys.path.append('/home/brandon/Documents/TF10/Task 10.2.2 Delivery v2/ap_gen_upscaling/required_functions')
from req_functions import RMS_COR

def make_labels(ax, boxplot):
    # Grab the relevant Line2D instances from the boxplot dictionary
    iqr = boxplot['boxes'][0]
    caps = boxplot['caps']
    med = boxplot['medians'][0]
    # fly = boxplot['fliers'][0]
    # The x position of the median line
    xpos = med.get_xdata()
    # Lets make the text have a horizontal offset which is some 
    # fraction of the width of the box
    xoff = 0.10 * (xpos[1] - xpos[0])
    # The x position of the labels
    xlabel = xpos[1] + xoff
    # The median is the y-position of the median line
    median = med.get_ydata()[1]
    # The 25th and 75th percentiles are found from the
    # top and bottom (max and min) of the box
    pc25 = iqr.get_ydata().min()
    pc75 = iqr.get_ydata().max()
    # The caps give the vertical position of the ends of the whiskers
    capbottom = caps[0].get_ydata()[0]
    captop = caps[1].get_ydata()[0]
    # Make some labels on the figure using the values derived above
    ax.text(xlabel, median,'Median = {:6.3g}'.format(median), va='center')#, fontsize=8)
    ax.text(xlabel, pc25,'25th percentile = {:6.3g}'.format(pc25), va='center')
    ax.text(xlabel, pc75,'75th percentile = {:6.3g}'.format(pc75), va='center')
    ax.text(xlabel, capbottom,'Bottom cap = {:6.3g}'.format(capbottom), va='center')
    ax.text(xlabel, captop,'Top cap = {:6.3g}'.format(captop), va='center')
    
    
##path to pre processed surface data
path_to_surface_data='/home/brandon/Documents/TF10/Task 10.2.2 Delivery v2/ap_gen_upscaling/Pre_processed_data'

## Load surface data
upper=np.load(path_to_surface_data+'/upper_surface_uncorrected_res_0_1mm.npy')[15:1985,15:1985]
lower=np.load(path_to_surface_data+'/lower_surface_uncorrected_res_0_1mm.npy')[15:1985,15:1985]
upper=upper[400:1400,400:1400]
lower=lower[400:1400,400:1400]

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
    
    np.save('H_data/H_lower_x_1_small_section.npy',H_lower_x)
    np.save('H_data/H_lower_y_1_small_section.npy',H_lower_y)
    np.save('H_data/H_upper_x_1_small_section.npy',H_upper_x)
    np.save('H_data/H_upper_y_1_small_section.npy',H_upper_y)

    np.save('H_data/int_lower_x_1_small_section.npy',int_lower_x)
    np.save('H_data/int_lower_y_1_small_section.npy',int_lower_y)
    np.save('H_data/int_upper_x_1_small_section.npy',int_upper_x)
    np.save('H_data/int_upper_y_1_small_section.npy',int_upper_y)
    
    ## plot boxplot of all lower
    Hs_lower_all_array=H_lower_x+H_lower_y
    red_diamond = dict(markerfacecolor='r', marker='D')
    fig3, ax3 = plt.subplots()
    ax3.set_title('lower all, 1 sections')#+str(title))
    ax3.set_ylabel('Hurst exponent')
    # Create the boxplot and store the resulting python dictionary
    my_boxes = ax3.boxplot(Hs_lower_all_array,showfliers=True)#, flierprops=red_diamond)
    iqr = my_boxes['boxes'][0]
    caps = my_boxes['caps']
    med = my_boxes['medians'][0]
    median = med.get_ydata()[1]
    pc25 = iqr.get_ydata().min() ## 25th perentile 
    pc75 = iqr.get_ydata().max() ## 75th 
    capbottom = caps[0].get_ydata()[0]
    captop = caps[1].get_ydata()[0]
    data=[median,pc25,pc75,capbottom,captop]
    d={'Index_Title':['median','pc25','pc75','capbottom','captop'],'Values':[median,pc25,pc75,capbottom,captop]}
    df = pd.DataFrame(d).set_index('Index_Title')
    df.to_csv('H_data/H_lower_all_1_params_small_section.csv')
    # Call the function to make labels
    make_labels(ax3, my_boxes)
    plt.savefig('figs/H_lower_all_1_small_section.png', dpi=600, bbox_inches='tight')
    plt.show()

    int_lower_all_array=int_lower_x+int_lower_y
    red_diamond = dict(markerfacecolor='r', marker='D')
    fig3, ax3 = plt.subplots()
    ax3.set_title('lower all, 1 sections')#+str(title))
    ax3.set_ylabel('Intercept')
    # Create the boxplot and store the resulting python dictionary
    my_boxes = ax3.boxplot(int_lower_all_array,showfliers=True)#, flierprops=red_diamond)
    iqr = my_boxes['boxes'][0]
    caps = my_boxes['caps']
    med = my_boxes['medians'][0]
    median = med.get_ydata()[1]
    pc25 = iqr.get_ydata().min() ## 25th perentile 
    pc75 = iqr.get_ydata().max() ## 75th 
    capbottom = caps[0].get_ydata()[0]
    captop = caps[1].get_ydata()[0]
    data=[median,pc25,pc75,capbottom,captop]
    d={'Index_Title':['median','pc25','pc75','capbottom','captop'],'Values':[median,pc25,pc75,capbottom,captop]}
    df = pd.DataFrame(d).set_index('Index_Title')
    df.to_csv('H_data/int_lower_all_1_params_small_section.csv')
    # Call the function to make labels
    make_labels(ax3, my_boxes)
    plt.savefig('figs/scaling_lower_all_1_small_section.png', dpi=600, bbox_inches='tight')
    plt.show()
    
    ## plot boxplot of all upper
    Hs_upper_all_array=H_upper_x+H_upper_y
    red_diamond = dict(markerfacecolor='r', marker='D')
    fig3, ax3 = plt.subplots()
    ax3.set_title('upper all, 1 sections')#+str(title))
    ax3.set_ylabel('Hurst exponent')
    # Create the boxplot and store the resulting python dictionary
    my_boxes = ax3.boxplot(Hs_upper_all_array,showfliers=True)#, flierprops=red_diamond)
    iqr = my_boxes['boxes'][0]
    caps = my_boxes['caps']
    med = my_boxes['medians'][0]
    median = med.get_ydata()[1]
    pc25 = iqr.get_ydata().min() ## 25th perentile 
    pc75 = iqr.get_ydata().max() ## 75th 
    capbottom = caps[0].get_ydata()[0]
    captop = caps[1].get_ydata()[0]
    data=[median,pc25,pc75,capbottom,captop]
    d={'Index_Title':['median','pc25','pc75','capbottom','captop'],'Values':[median,pc25,pc75,capbottom,captop]}
    df = pd.DataFrame(d).set_index('Index_Title')
    df.to_csv('H_data/H_upper_all_1_params_small_section.csv')
    # Call the function to make labels
    make_labels(ax3, my_boxes)
    plt.savefig('figs/H_upper_all_1_small_section.png', dpi=600, bbox_inches='tight')
    plt.show()

    int_upper_all_array=int_upper_x+int_upper_y
    red_diamond = dict(markerfacecolor='r', marker='D')
    fig3, ax3 = plt.subplots()
    ax3.set_title('upper all, 1 sections')#+str(title))
    ax3.set_ylabel('Intercept')
    # Create the boxplot and store the resulting python dictionary
    my_boxes = ax3.boxplot(int_upper_all_array,showfliers=True)#, flierprops=red_diamond)
    iqr = my_boxes['boxes'][0]
    caps = my_boxes['caps']
    med = my_boxes['medians'][0]
    median = med.get_ydata()[1]
    pc25 = iqr.get_ydata().min() ## 25th perentile 
    pc75 = iqr.get_ydata().max() ## 75th 
    capbottom = caps[0].get_ydata()[0]
    captop = caps[1].get_ydata()[0]
    data=[median,pc25,pc75,capbottom,captop]
    d={'Index_Title':['median','pc25','pc75','capbottom','captop'],'Values':[median,pc25,pc75,capbottom,captop]}
    df = pd.DataFrame(d).set_index('Index_Title')
    df.to_csv('H_data/int_upper_all_1_params_small_section.csv')
    # Call the function to make labels
    make_labels(ax3, my_boxes)
    plt.savefig('figs/scaling_upper_all_1_small_section.png', dpi=600, bbox_inches='tight')
    plt.show()
   
    ## get anisotropy data
    ## lower surface
    lower_aniso=np.asarray(H_lower_x)/np.asarray(H_lower_y)   
    np.save('H_data/lower_aniso_1_small_section.npy',lower_aniso)

    ## boxplot of aniso lower   
    red_diamond = dict(markerfacecolor='r', marker='D')
    fig3, ax3 = plt.subplots()
    ax3.set_title('Lower surface')#+str(title))
    ax3.set_ylabel('Anisotropy')
    # Create the boxplot and store the resulting python dictionary
    my_boxes = ax3.boxplot(lower_aniso,showfliers=True)#, flierprops=red_diamond)
    iqr = my_boxes['boxes'][0]
    caps = my_boxes['caps']
    med = my_boxes['medians'][0]
    median = med.get_ydata()[1]
    pc25 = iqr.get_ydata().min() ## 25th perentile 
    pc75 = iqr.get_ydata().max() ## 75th 
    capbottom = caps[0].get_ydata()[0]
    captop = caps[1].get_ydata()[0]
    data=[median,pc25,pc75,capbottom,captop]
    d={'Index_Title':['median','pc25','pc75','capbottom','captop'],'Values':[median,pc25,pc75,capbottom,captop]}
    df = pd.DataFrame(d).set_index('Index_Title')
    df.to_csv('H_data/aniso_lower_params_1_small_section.csv')    
    # Call the function to make labels
    make_labels(ax3, my_boxes)
    plt.savefig('figs/lower_aniso_1_small_section.png', dpi=600, bbox_inches='tight')
    plt.show()
    
    ## upper surface
    upper_aniso=np.asarray(H_upper_x)/np.asarray(H_upper_y)   
    np.save('H_data/upper_aniso_1_small_section.npy',lower_aniso)

    ## boxplot of aniso lower   
    red_diamond = dict(markerfacecolor='r', marker='D')
    fig3, ax3 = plt.subplots()
    ax3.set_title('Upper surface')#+str(title))
    ax3.set_ylabel('Anisotropy')
    # Create the boxplot and store the resulting python dictionary
    my_boxes = ax3.boxplot(upper_aniso,showfliers=True)#, flierprops=red_diamond)
    iqr = my_boxes['boxes'][0]
    caps = my_boxes['caps']
    med = my_boxes['medians'][0]
    median = med.get_ydata()[1]
    pc25 = iqr.get_ydata().min() ## 25th perentile 
    pc75 = iqr.get_ydata().max() ## 75th 
    capbottom = caps[0].get_ydata()[0]
    captop = caps[1].get_ydata()[0]
    data=[median,pc25,pc75,capbottom,captop]
    d={'Index_Title':['median','pc25','pc75','capbottom','captop'],'Values':[median,pc25,pc75,capbottom,captop]}
    df = pd.DataFrame(d).set_index('Index_Title')
    df.to_csv('H_data/aniso_upper_params_1_small_section.csv')    
    # Call the function to make labels
    make_labels(ax3, my_boxes)
    plt.savefig('figs/upper_aniso_1_small_section.png', dpi=600, bbox_inches='tight')
    plt.show()