#!bin/bash python3

from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import re
import os

#samplesheet
# - construct
# - lanes
# - skip lanes

#look into the formula and half life functions to make sure they match up with stack exchange version


#I define the half life formula as a funcion
def formula(x,a,b,k):
    return (a-b)*np.exp(-k*x)+b


#This function calculates half life given the raw integrity data and timepoints
def half_life(x,y,factor):
    ynorm=np.array([a*factor for a in y])
    p0=(1.,1.e-5,1.)
    #opt,pcov=curve_fit(formula,x,ynorm,bounds=(0,[100.,1.]))
    opt,pcov=curve_fit(formula,x,ynorm)
    a,b,k=opt
    #print(k,b)
    y2=formula(x,a,b,k)
    residuals=ynorm-formula(x,*opt)
    ss_res=np.sum(residuals**2)
    ss_tot=np.sum((ynorm-np.mean(ynorm))**2)
    R_sq=1-(ss_res/ss_tot)
    hl=np.log(2)/k
    print(hl)
    return x,y2,ynorm,hl,R_sq

def AutoHL(samplesheet):
    '''
    AutoHL reads a samplesheet containing paths to smear tables generated by AutoSmear, experimental groups and their corresponding lanes, 
    and any lanes that need to be skipped due to low quality readouts. The function is able to derive associated timepoints based on naming conventions,
    with the timepoint following the name of the sample and an underscore (ex. pBPD16_4hr, will read that it is the 4 hr timepoint). 

    A folder called `half_life_analysis` is created, containing individual half-life curves and normalized integrity tablesfor each of the experimental 
    groups/constructs specified in the samplesheet. Additionally, a plot containing half-life curves of all constructs is created along woth a table containing the 
    calculated half-life for all of the contructs specified in the samplesheet

    '''
    #read in samplesheet
    sample_sheet=pd.read_csv(samplesheet,sep='\t')
    #set up a dictionary that will contain experimental groups as keys, with their corresponding rows from each of the files specified in the samplesheet
    data_dict={}
    os.chdir("..")
    #set up the half-life results table
    hl_table=pd.DataFrame(columns=['Construct','Half_Life'])
    #get the list of files to be analyzed
    file_set=sample_sheet["File"].tolist()
    #for each file, read the smear table in and find all of the samples that are to be analyzed in this file
    for file in file_set:
        smear_table=pd.read_csv(file)
        #the sample groups contained within this file are subset from the rest of the samplesheet
        samples=sample_sheet[sample_sheet["File"]==file]
        #for each sample in the file, create/add to a dictionary entry that will store the selected rows minus the skipped rows in 
        #the samplesheet
        for s in samples['Construct'].tolist():
            #for each sample group in the subset, I find their corresponding lanes as well as lanes that will be skipped
            lanes=samples[samples["Construct"]==s]['Lanes'].to_string(index=False)
            lane_list=lanes.split(',')
            skip=samples[samples["Construct"]==s]['Skip Lanes'].to_string(index=False)
            skip_list=skip.split(',')
            skip_list=[x.replace(' ','') for x in skip_list]
            sample_lanes=[]
            #for any lanes designated with a range (eg. A1-A8), this format is processed below
            for l in lane_list:
                l=l.replace(' ','')
                if '-' in l:
                    fl=l.split('-')
                    letter=fl[0][0]
                    inds=[re.findall(r'\d+', n) for n in fl]
                    for x in range(int(inds[0][0]),int(inds[1][0])+1):
                        addlane=letter+str(x)
                        sample_lanes.append(addlane)
                else:
                    sample_lanes.append(l)
            #lanes contained within the skip list are removed from the data
            for y in skip_list:
                if y in sample_lanes:
                    sample_lanes.remove(y)
            #If the sample group is already present in the dictionary, the table rows are appended to the existing table
            if s in data_dict.keys():
                data_dict[s]=pd.concat([data_dict[s], smear_table[smear_table["Lane"].isin(sample_lanes)]], ignore_index=True, axis=0)
            #Else, a new sample group key is created and the table rows are added to the dictionary
            else:
                data_dict[s]=smear_table[smear_table["Lane"].isin(sample_lanes)]
    #I then loop through data for each sample group and extract the timepoint from the sample name
    for key,tab in data_dict.items():
        timescale=[]
        IDs=tab["ID"].tolist()
        for a in IDs:
            #will find the timepoint as long as it is the last number found in the sample ID
            hour=re.findall(r'\d+', a)[-1]
            timescale.append(hour)
        #a new timepoint co0lumn is created
        data_dict[key]=tab.assign(timepoint=[eval(i) for i in timescale])
    #a list that will contain half-life values for each of the sample groups is initialized
    hls=[]
    #find the summary statistics
    plot_data=[]
    #for each of the sample groups, I generate a normalized integrity table containing averages across replicates, normalized average integrity, number of
    #replicates, and standard deviation across replicated for each timepoint
    for key,val in data_dict.items():
        norm_data=pd.DataFrame(columns=['Timepoint','Replicates','Average_Integrity','Normalized_Integrity','Standard_Deviation'])
        time=[]
        integrity=[]
        stdevs=[]
        replicates=[]
        #loop through all unique timepoints
        for num in list(set(val["timepoint"].tolist())):
            time.append(num)
            reps=val.loc[val['timepoint'] == num]
            replicates.append(len(reps))
            integs=reps["Percent_Total"].tolist()
            average_int=sum(integs)/len(integs)
            variance=sum([((x - average_int) ** 2) for x in integs]) / len(integs) 
            stdev_int= variance ** 0.5
            stdevs.append(stdev_int)
            integrity.append(average_int)
            if num==0:
                norm_factor=100/average_int
        #I then write the table to a csv file and perform half-life calcutations using average integrity and timepoint values
        norm_data=norm_data.assign(Timepoint=time)
        norm_data=norm_data.assign(Replicates=replicates)
        norm_data=norm_data.assign(Average_Integrity=[round(x,2) for x in integrity])
        norm_data=norm_data.assign(Standard_Deviation=stdevs)
        model=half_life(np.array(time),np.array(integrity),norm_factor)
        norm_data=norm_data.assign(Normalized_Integrity=[round(x,2) for x in model[2]])
        os.makedirs('half_life_results',exist_ok=True)
        norm_data.to_csv('half_life_results/'+key+'_normalized_data.csv',index=False)
        #I append the half-life values to the list and add the half-life curve points to a list that will help with plotting all of the sample groups together
        hls.append(model[3])
        plot_data.append([model[0],model[1],model[2]])
        #I then write the individual plots
        plt.plot(model[0],model[1],label=key)
        plt.scatter(model[0],model[2],marker='o')
        plt.title("Half-Life of "+str(key)+" at 10mM Mg")
        plt.legend()
        plt.xlabel('Time (hr)')
        plt.ylabel('Normalized RNA Integrity (%)')
        plt.savefig('half_life_results/'+str(key)+"_half_life.jpg",bbox_inches='tight')
        plt.clf()
        #plt.errorbar(model[0], model[2], stdevs, linestyle='None', capsize=3)
        #plt.scatter(manual_half[0],manual_half[2],color='#808284',marker='o')
    #I then create a table containing half-life values for all experimental groups and write it to a csv
    hl_table=hl_table.assign(Half_Life=hls)
    hl_table=hl_table.assign(Construct=data_dict.keys())
    hl_table.to_csv('half_life_results/half_life_table.csv',index=False)
    #I then loop through the plot data for each sample group and plot them together
    for p in plot_data:
        print("pl;ot")
        plt.plot(p[0],p[1],label=key)
        plt.scatter(p[0],p[2],marker='o')
    plt.title("Half-Life of Constructs at 10mM Mg")
    plt.legend()
    plt.xlabel('Time (hr)')
    plt.ylabel('Normalized RNA Integrity (%)')
    #The plot is then saved as a png
    plt.savefig("half_life_results/half_life.jpg",bbox_inches='tight')
    
    
    return "Half-Life Analysis Complete"

#the samplesheet found within the current folder is used as an input for the AutoHL function
for file in os.listdir():
    if file.endswith("samplesheet.tsv"):
        ss_path=os.path.join(os.getcwd(),file)
        AutoHL(ss_path)
print("Half-Life Analysis Complete")