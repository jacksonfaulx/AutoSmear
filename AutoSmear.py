#!/usr/bin/env python3
import sys
import os
import matplotlib.pyplot as plt # type: ignore
import pandas as pd # type: ignore
import numpy as np # type: ignore
from scipy.signal import find_peaks, argrelmin, savgol_filter # type: ignore
import math

def AutoSmear(file,folder,length=10000):
    ''' 
    Gates are placed +/- the `window` size from the target fragment peak. 
    The gates are then further refined using curve-smoothing
    and relative minima to find the true boundaries of the peak.

    Calculated RNA integrity along with the peak boundaries
    along with grid of plots visualizing the gating for each capillary
    are output as a .csv file to 'smear_results' folder
    '''

    #a window used to look for the main peak is established. After lots of testing, 1200 seems to provide the most accurate results
    window=1200
    #data is read in
    try:
        data=pd.read_csv(file)
    except Exception:
        print("Error: could not read file. Please make sure the input file is in csv format and in the current directory")
    #the lengths of columns and the names of each sample lane are defined
    cols=len(data.columns)
    names=list(data.columns.values)
    names_split=[n.split(": ") for n in names[:-1]]
    #input check to make sure that the data was correctly exported from ProSize
    if names_split[0][0]!='Size (nt)':
        if names_split[0][0]!='Size (bp)':
            raise Exception("Error: x-axis values must be size (nt). Please export electropherogram data from ProSize with size (nt) as the x-axis")
    #the results table is then set up with a list of sample names and their associated capillary
    IDs=[i[1] for i in names_split[1:]]
    lanes=[l[0] for l in names_split[1:]]
    smear_table=pd.DataFrame(lanes,columns=['Lane'])
    smear_table=smear_table.assign(ID=IDs)
    #x values (size in nt) are established
    xvals=data.iloc[:,0]
    #A few empty lists are established to add data to as the algorithm progresses
    percents=[]
    left_gates=[]
    right_gates=[]
    #If there is no peak found for the first lane, the y value at x=length is set
    pp=xvals[min(range(len(xvals)),key=lambda i: abs(xvals[i]-length))]
    previous_peak=list(xvals).index(pp)
    #I then set up the grid dimensions for the peak calling visualization graphs
    sqrt = math.sqrt(cols-2)
    closest = [math.floor(sqrt)**2,math.ceil(sqrt)**2]
    nearsq=closest[min(range(len(closest)),key=lambda i: abs(closest[i]-(cols-2)))]
    #This code constructs the grid of gating plots 
    if nearsq>(cols-2):
        sqr=math.sqrt(nearsq)
        sqy=math.sqrt(nearsq)
        sqx=math.sqrt(nearsq)
        diff=nearsq-(cols-2)
        while diff>sqr:
            sqy-=1
            diff-=sqr
        fig,axs=plt.subplots(int(sqy),int(sqx))
    elif nearsq<(cols-2):
        sqr=math.sqrt(nearsq)
        sqy=math.sqrt(nearsq)+1
        sqx=math.sqrt(nearsq)
        diff=(cols-2)-nearsq
        while diff>sqr:
            sqy+=1
            diff-=sqr
        fig,axs=plt.subplots(int(sqy),int(sqx))
    elif nearsq==(cols-2):
        sqy=math.sqrt(nearsq)
        sqx=math.sqrt(nearsq)
        fig,axs=plt.subplots(int(sqy),int(sqx))
    grid_ax=[0,0]
    start_peak=xvals[min(range(len(xvals)),key=lambda i: abs(xvals[i]-10000))]
    previous_peak=list(xvals).index(start_peak)
    #I then iterate over each sample run (capillary)
    for i in range(1,cols-1):
        #RFU (y-values) are established
        yvals=data.iloc[:,i]
        #peak calling is performed, and the peak with the largest size (x-value) is isolated. the minimum peak height of 50 is used (same as ProSize)
        peaks,_=find_peaks(yvals,height=50)
        peak_list=[]
        #peaks that are between 9000 and 15000 bp length are added to the peaks list
        for p in peaks:
            if xvals[p]>(length-1000) and xvals[p]<(length+5000):
                peak_list.append(p)
        #if no peaks are detected, then the peak determined from the previous capillary is used to as the peak
        if len(peak_list)==0:
            peak=previous_peak
        #else, the peak with the highest x-value (farthest right) in the list of called peaks is used as the peak of interest
        else:
            peak_ys=[yvals[i] for i in peak_list]
            max_idx=peak_ys.index(max(peak_ys))
            peak=peak_list[max_idx]
            previous_peak=peak
            #If there are "twin peaks" the rightmost one is selected
            for p in peak_list:
                #twin peaks are defined as two peaks close together that are at least 60% the height of the taller one
                if (yvals[p]/yvals[peak])>0.6 and p!=peak:
                    peak=max(peak,p)
        #A window of values around the peak is defined to perform curve smoothing and local minimum analysis
        left_gate=xvals[peak]-window
        right_gate=xvals[peak]+window
        #I find the value closest to the window x-value as well as its index
        lg=xvals[min(range(len(xvals)),key=lambda i: abs(xvals[i]-left_gate))]
        rg=xvals[min(range(len(xvals)),key=lambda i: abs(xvals[i]-right_gate))]
        lg_idx=list(xvals).index(lg)
        rg_idx=list(xvals).index(rg)
        #The points between each gate and the peak are subset from the complete data (left and right sides) in order to perform curve smoothing
        x_left=np.array(xvals[lg_idx:peak])
        y_left=np.array(yvals[lg_idx:peak])
        x_right=np.array(xvals[peak+1:rg_idx+1])
        y_right=np.array(yvals[peak+1:rg_idx+1])
        #the curves are smoothed and fitted to x^2 to easily identify the edges of the peak
        y_left=savgol_filter(y_left,int(2*len(y_left)/3),2)
        y_right=savgol_filter(y_right,int(2*len(y_right)/3),2)
        #The relative minima of the smoothed curve are found and identified as the boundaries of the peak, with the minima closest to the peak
        #used as a marker for where the gates will be placed
        #left
        mins=list(argrelmin(y_left))
        left_midx=mins[-1]
        #if there are no minima found, then gates are placed 1000 nt from the peak
        if len(left_midx)==0:
            nomins=x_left[min(range(len(x_left)),key=lambda i: abs(x_left[i]-(x_left[-1]-1000)))]
            left_midx=list(x_left).index(nomins)
        else:
            left_midx=left_midx[-1]
        #right
        rmins=list(argrelmin(y_right))
        right_midx=rmins[0]
        #if there are no minima found, then gates are placed 1000 nt from the peak
        if len(right_midx)==0:
            nomins=x_right[min(range(len(x_right)),key=lambda i: abs(x_right[i]-(x_right[0]+1000)))]
            right_midx=list(x_right).index(nomins)
        else:
            right_midx=right_midx[0]
        #The indices of the two gates in the complete data are found
        lidx=list(xvals).index(x_left[left_midx])
        ridx=list(xvals).index(x_right[right_midx])
        left_gates.append(round(xvals[lidx]))
        right_gates.append(round(xvals[ridx]))
        #the integration boundaries are then set to leave out the control spike and the noise at the end of the data
        left_bound=xvals[min(range(len(xvals)),key=lambda i: abs(xvals[i]-125))]
        right_bound=xvals[min(range(len(xvals)),key=lambda i: abs(xvals[i]-13800))]
        lb_idx=list(xvals).index(left_bound)
        rb_idx=list(xvals).index(right_bound)
        #The total area and area within the gates are then calculated to give a percentage
        total_area=np.trapezoid(y=yvals[lb_idx:rb_idx+1],x=xvals[lb_idx:rb_idx+1])
        gated_area=np.trapezoid(y=yvals[lidx:ridx+1],x=xvals[lidx:ridx+1])
        percent=(gated_area/total_area)*100
        percents.append(round(percent,2))
        #The data and gates are then plotted in a grid
        plot_x=xvals[min(range(len(xvals)),key=lambda i: abs(xvals[i]-15000))]
        cutoff=list(xvals).index(plot_x)
        axs[grid_ax[0],grid_ax[1]].plot(xvals[:cutoff],yvals[:cutoff])
        axs[grid_ax[0],grid_ax[1]].plot(xvals[peak],yvals[peak],'x')
        axs[grid_ax[0],grid_ax[1]].axvline(x=x_right[right_midx],color='r',linestyle='dashed')
        axs[grid_ax[0],grid_ax[1]].axvline(x=x_left[left_midx],color='r',linestyle='dashed')
        axs[grid_ax[0],grid_ax[1]].set_title(names[i])
        grid_ax[1]+=1
        if grid_ax[1]==int(sqx):
            grid_ax[1]=0
            grid_ax[0]+=1
    if cols<90:
        plt.subplots_adjust(bottom=0.1,right=3,top=3)
    else:
        plt.subplots_adjust(bottom=0.1,right=10,top=10)
    grid_total=int(sqy*sqx)
    gdiff=grid_total-(cols-2)
    if gdiff>0:
        for num in range(int(gdiff)):
            fig.delaxes(axs[int(sqy)-1][(int(sqx)-1)-int(num)])
    #following the above analysis for each run, the results table is updated and written into a csv file in a new folder called "smear_analysis"
    smear_table=smear_table.assign(Left_Gate=left_gates)
    smear_table=smear_table.assign(Right_Gate=right_gates)
    smear_table=smear_table.assign(Percent_Total=percents)
    os.makedirs('smear_results',exist_ok=True)
    smear_table.to_csv('smear_results/'+folder+'_smear_results.csv',index=False)
    #the graphs of peak calling for each sample are written to a jpg in the same folder
    plt.savefig("smear_results/"+folder+"_peaks.jpg",bbox_inches='tight')
    return print(folder+' smear analysis complete')


#for each folder, the program looks for the electropherogram data file
os.chdir("..")
#for each data folder exported from ProSize, the electropherogram file is located and analyzed by the function above
for folder in os.listdir():
    input_folder=folder
    sub_path=os.path.join(os.getcwd(),folder)
    if os.path.isdir(sub_path):
        for file in os.listdir(sub_path):
            if file.endswith('Electropherogram.csv'):
                #The above function is then executed on the data and the results are written to a new folder called smear_results
                input_file=os.path.join(sub_path,file)
                AutoSmear(input_file,input_folder)
print("All smear analyses complete")

