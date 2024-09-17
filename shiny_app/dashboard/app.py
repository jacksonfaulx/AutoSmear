import sys
sys.path.append('/opt/anaconda3_24/lib/python3.11/site-packages/')
import seaborn as sns
from faicons import icon_svg
from shiny import App, reactive, render, ui, req, Outputs, Session
import os
import re
import io
import matplotlib.pyplot as plt # type: ignore
import pandas as pd # type: ignore
import numpy as np # type: ignore
from pathlib import Path
from scipy.signal import find_peaks, argrelmin, savgol_filter # type: ignore
from scipy.optimize import curve_fit
import math


app_dir = Path(__file__).parent
app_ui = ui.page_fluid(
    ui.page_navbar(
        ui.nav_spacer(),
        ui.nav_panel(ui.h4("AutoSmear"),
                ui.layout_sidebar(
                    ui.sidebar(
                        ui.h2("Automated Smear Analysis"),
                        ui.p("This tool is designed to automatically define the edges of the target peak to calculate the percent RNA integrity from fragment analyzer electropherogram data exported from ProSize."),
                        ui.h2("Data Input"),
                        ui.input_file("input_csv",ui.h4("Choose Electropherogram File (CSV)"),accept=[".csv"],multiple=False),
                        ui.p("Upload an Electropherogram.csv file exported from ProSize, with size (nt) as the x-axis value"),
                        ui.input_text("exp_name",ui.h4("Experiment Name"),placeholder="ex. EXP_01"),
                        ui.p("Type in a name for the experiment. This name will be added to the file names of the downloadable plots and tables"),
                        ui.input_action_button('submit','Submit'),
                        width='20%'
                        ),
                ui.layout_columns(
                    ui.card(
                        ui.card_header("Smear Table"),
                        ui.output_data_frame("smear_table"),
                        ui.card_footer(
                            ui.download_button("smear_download","Download Smear Table"),
                            ),
                        full_screen=True,
                        ),
                    ui.card(
                        ui.card_header("Gating Plots"),
                        #ui.output_plot("gating"),
                        ui.output_image('plots',fill=True),
                        ui.card_footer(
                            ui.download_button("plot_download","Download Smear Plots"),
                            ),
                        full_screen=True,
                        )
                    )
                )
            ),
        ui.nav_spacer(),
        ui.nav_panel(ui.h4("AutoHL"),
                ui.layout_sidebar(
                    ui.sidebar(
                        ui.h2("Automated Half-Life Calculation"),
                        ui.p("This tool is designed to fit the integrity values calculated using ProSize or AutoSmear to a one-phase decay equation, plotting the resulting curve and calculating the half-life of a construct based on timepoints.."),
                        #ui.p("The script takes a samplesheet as an input, indicating the paths of relevant smear analysis files to be analyzed, along with details regarding sample groups and their corresponding lanes in the smear table. If multiple constructs are present within the data, a combined plot and data table will be generated as an output for comparison"),
                        ui.h2("Data Input"),
                        ui.input_file("input_hl",ui.h4("Choose Smear Results File (CSV)"),accept=[".csv"],multiple=True),
                        ui.p("Upload one or more smear results files (.csv format), obtained from AutoSmear or ProSize"),
                        ui.input_file("input_ss",ui.h4("Upload Sample Sheet"),accept=[".csv",".tsv",".txt"],multiple=False),
                        ui.p("Upload a samplesheet formatted as detailed below. An example template can be downloaded below"),
                        ui.tags.ul(
                            ui.tags.li("experiment_name"),
                                ui.tags.ul(
                                 ui.tags.li("only need to specify once in the first row, name is used for combined outputs"),
                             ),
                            ui.tags.li("File"),
                                ui.tags.ul(
                                 ui.tags.li("The full name of a smear results file (NOT full path)"),
                             ),
                            ui.tags.li("Construct"),
                                ui.tags.ul(
                                 ui.tags.li("name of the sample being tested across timepoints"),
                             ),
                            ui.tags.li("Lanes"),
                                ui.tags.ul(
                                 ui.tags.li("names of all the lanes associated with the Construct. Lanes must be separated by commas with no spaces in between"),
                             ),
                            ui.tags.li("Skip_Lanes"),
                                ui.tags.ul(
                                 ui.tags.li("lanes containing bad quality data can be listed here and will be skipped in analysis. Same format as 'Lanes' column"),
                             ),
                        ),
                        ui.download_link("example_sheet","Download Samplesheet Template"),
                        ui.input_text("exp_hl",ui.h4("Experiment Name"),placeholder="ex. EXP_01"),
                        ui.p("Type in a name for the experiment. This name will be added to the file names of the downloadable plots and tables"),
                        ui.input_action_button('submit_hl','Submit'),
                        ui.output_ui("sample_select"),
                        width='20%'
                        ),
                ui.layout_columns(
                    ui.card(
                        ui.card_header("Normalized Smear Summary Table"),
                        ui.output_data_frame("norm_smear"),
                        ui.card_footer(
                            ui.download_button("norm_smear_download","Download Summary Table"),
                            ),
                        full_screen=True,
                        ),
                    ui.card(
                        ui.card_header("Construct Half-Life Plots"),
                        #need to create a slideshow format to go through all of the individual plots for each sample
                        #also need to brainstorm how to download the individual sample group plots
                        ui.output_plot('hl_plots',fill=True),
                        ui.card_footer(
                            ui.download_button("hl_plots_download","Download Half-Life Plots"),
                            ),
                        full_screen=True,
                        )
                    ),
                ui.layout_columns(
                    ui.card(
                        ui.card_header("Half-Life Table"),
                        ui.output_data_frame('hl_table'),
                        ui.card_footer(
                            ui.download_button("hl_table_download","Download Half-Life Table"),
                            ),
                        full_screen=True,
                        ),
                    ui.card(
                        ui.card_header("Combined Half-Life Plot"),
                        ui.output_plot('combined_plot',fill=True),
                        ui.card_footer(
                            ui.download_button("combined_plot_download","Download Combined Half-Life Plot"),
                            ),
                        full_screen=True,
                        )
                    )
                )
            ),
        ui.nav_spacer(),
        title=ui.div(ui.h1("Automated Smear Analysis"),
                     ui.h1("and Half-Life Calculator")
        )
    )
)


def server(input, output, session):

# AutoSmear - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def AutoSmear(data,length=10000):
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
            total_area=np.trapz(y=yvals[lb_idx:rb_idx+1],x=xvals[lb_idx:rb_idx+1])
            gated_area=np.trapz(y=yvals[lidx:ridx+1],x=xvals[lidx:ridx+1])
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
        #the graphs of peak calling for each sample are written to a jpg in the same folder
        plt.savefig("peaks_plot.png",bbox_inches='tight')
        return [smear_table,fig]
    

    # REACTIVES

    @reactive.calc
    def parsed_file():
        req(input.input_csv(),input.submit())
        file: list[FileInfo] | None = input.input_csv()
        return pd.read_csv( # pyright: ignore[reportUnknownMemberType]
            file[0]["datapath"])

    @reactive.calc
    def smear_results():
        req(input.input_csv(),input.exp_name(),input.submit())
        df=parsed_file()
        smear_data=AutoSmear(df)
        return smear_data


    @render.data_frame
    def smear_table():
        req(input.input_csv(),input.exp_name(),input.submit())
        return render.DataGrid(smear_results()[0],width='100%')

    
    @render.download(filename=lambda: f'{input.exp_name()}_smear_table.csv')
    def smear_download():
        req(input.input_csv(),input.exp_name(),input.submit())
        with io.BytesIO() as buf:
            smear_results()[0].to_csv(buf,index=False)
            yield buf.getvalue()

    @render.download(filename=lambda: f'{input.exp_name()}_peaks.png')
    def plot_download():
        req(input.input_csv(),input.exp_name(),input.submit())
        filename=str(input.exp_name())+'_peaks.png'
        path=os.path.join(os.path.dirname(__file__),"peaks_plot.png")
        return path

    @render.image
    def plots():
        req(input.input_csv(),input.exp_name(),input.submit())
        fullpath=os.path.join(os.path.dirname(__file__),"peaks_plot.png")
        img: ImgData = {"src": fullpath}
        return img
    
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#AutoHL - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# FUNCTIONS

    def formula(x,a,b,k):
        return (a-b)*np.exp(-k*x)+b

    def half_life(x,y,factor):
        ynorm=np.array([a*factor for a in y])
        #p0=(1.,1.e-5,1.)
        opt,pcov=curve_fit(formula,x,ynorm,bounds=([0,0,-np.inf],[100.,1.,np.inf]))
        #opt,pcov=curve_fit(formula,x,ynorm)
        a,b,k=opt
        #print(k,b)
        y2=formula(x,a,b,k)
        residuals=ynorm-formula(x,*opt)
        ss_res=np.sum(residuals**2)
        ss_tot=np.sum((ynorm-np.mean(ynorm))**2)
        R_sq=1-(ss_res/ss_tot)
        hl=np.log(2)/k
        #print(hl)
        return x,y2,ynorm,hl,R_sq
    
    def calculate_half_life(input_files,file_names,ss):
        #read in samplesheet
        sample_sheet=pd.read_csv(ss,sep='\t')
        data_dict={}
        #set up half life results table
        hl_table=pd.DataFrame(columns=['Construct','Half_Life'])
        #get the list of file to be analyzed
        #file_set=sample_sheet["File"].tolist()
        #print(file_set)
        #for each file, read the smear table in and find all of the samples that are to be analyzed in this file
        for x in range(len(input_files)):
            file=file_names[x]
            smear_table=pd.read_csv(input_files[x])
            samples=sample_sheet[sample_sheet["File"]==file]
            #print(samples)
            #for each sample in the file, create/add to a dictionary entry that will store the selected rows minus the skipped rows in 
            #the samplesheet
            #print(samples)
            for s in samples['Construct'].tolist():
                lanes=samples[samples["Construct"]==s]['Lanes'].to_string(index=False)
                lane_list=lanes.split(',')
                skip=samples[samples["Construct"]==s]['Skip_Lanes'].to_string(index=False)
                skip_list=skip.split(',')
                skip_list=[x.replace(' ','') for x in skip_list]
                #print(lane_list)
                sample_lanes=[]
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
                #print(sample_lanes)
                for y in skip_list:
                    if y in sample_lanes:
                        sample_lanes.remove(y)
                if s in data_dict.keys():
                    data_dict[s]=pd.concat([data_dict[s], smear_table[smear_table["Lane"].isin(sample_lanes)]], ignore_index=True, axis=0)
                else:
                    data_dict[s]=smear_table[smear_table["Lane"].isin(sample_lanes)]
    #for each sample group, want to make an individual half-life graph, including error bars for each timepoint 
        # Also want to return a table of normalized integrity values along with the number of replicates that werernt skipped for each timepoint
        #loop through data for each sample group and extract the timepoint from the sample name
        sample_list=[]
        for key,tab in data_dict.items():
            sample_list.append(key)
            timescale=[]
            IDs=tab["ID"].tolist()
            for a in IDs:
                #will find the timepoint as long as it is the last number found in the sample ID
                hour=re.findall(r'\d+', a)[-1]
                timescale.append(hour)
            data_dict[key]=tab.assign(timepoint=[eval(i) for i in timescale])
        #print(data_dict)
        hls=[]
        #find the summary statistics
        plot_data=[]
        st_dict={}
        plt_dict={}
        for key,val in data_dict.items():
            norm_data=pd.DataFrame(columns=['Timepoint','Replicates','Average_Integrity','Normalized_Integrity','Standard_Deviation'])
            time=[]
            integrity=[]
            stdevs=[]
            replicates=[]
            #loop through all unique timepoints
            for num in list(set(val["timepoint"].tolist())):
                #print(key,num)
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
            #print(time)
            #print(replicates)
            #print(integrity)
            #print(stdevs)
            norm_data=norm_data.assign(Timepoint=time)
            norm_data=norm_data.assign(Replicates=replicates)
            norm_data=norm_data.assign(Average_Integrity=integrity)
            norm_data=norm_data.assign(Standard_Deviation=[round(x,3) for x in stdevs])
            model=half_life(np.array(time),np.array(integrity),norm_factor)
            norm_data=norm_data.assign(Normalized_Integrity=[round(x,2) for x in model[2]])
            #print(norm_data)
            hls.append(model[3])
            plot_data.append([model[0],model[1],model[2],key])
            ind_plot=plt.figure()
            plt.plot(model[0],model[1],label=key)
            #plt.plot(manual_half[0],manual_half[1],'#808284')
            plt.scatter(model[0],model[2],marker='o')
            plt.title("Half-Life of "+str(key)+" at 10mM Mg")
            plt.legend()
            plt.xlabel('Time (hr)')
            plt.ylabel('Normalized RNA Integrity (%)')
            plt_dict[key]=ind_plot
            st_dict[key]=norm_data
            #plt.errorbar(model[0], model[2], stdevs, linestyle='None', capsize=3)
            #plt.scatter(manual_half[0],manual_half[2],color='#808284',marker='o')
        #plt.legend(['Automated','Manual (Mark)'])
        #Need a way to first plot individual plots, then add plot object to a list for graphing against other sample groups
        hl_table=hl_table.assign(Half_Life=[round(x,2) for x in hls])
        hl_table=hl_table.assign(Construct=data_dict.keys())
        hl_table.to_csv('half_life_table.csv',index=False)
        combined_plot=plt.figure()
        cm = plt.get_cmap('tab20')
        for c,p in enumerate(plot_data):
            plt.plot(p[0],p[1],label=p[3],color=cm(c))
            #plt.plot(manual_half[0],manual_half[1],'#808284')
            plt.scatter(p[0],p[2],marker='o',color=cm(c))
        plt.title("Half-Life of Constructs at 10mM Mg")
        plt.legend()
        plt.xlabel('Time (hr)')
        plt.ylabel('Normalized RNA Integrity (%)')
        
        return st_dict, hl_table, combined_plot, plt_dict, sample_list, plot_data

# REACTIVES

    @reactive.calc
    def hl_results():
        req(input.input_hl(),input.exp_hl(),input.submit_hl(),input.input_ss())
        inputs=input.input_hl()
        path_list=[]
        name_list=[]
        for dict in inputs:
            name_list.append(dict['name'])
            path_list.append(dict['datapath'])
        hl_data=calculate_half_life(path_list,name_list,input.input_ss()[0]['datapath'])
        #print(hl_data[0])
        return hl_data
    
    @render.ui
    def sample_select():
        req(input.input_hl(),input.exp_hl(),input.submit_hl(),input.input_ss(),hl_results())
        return ui.input_select('samples','Choose A Sample to Visualize',
                               hl_results()[4])


    @render.data_frame
    def norm_smear():
        req(input.input_hl(),input.exp_hl(),input.submit_hl(),input.input_ss(),input.samples())
        return render.DataGrid(hl_results()[0][input.samples()],width='100%')
    
    @render.data_frame
    def hl_table():
        req(input.input_hl(),input.exp_hl(),input.submit_hl(),input.input_ss())
        return render.DataGrid(hl_results()[1],width='100%')

    
    @render.plot
    def combined_plot():
        req(input.input_hl(),input.exp_hl(),input.submit_hl(),input.input_ss())
        return hl_results()[2]
    
    @render.plot
    def hl_plots():
        req(input.input_hl(),input.exp_hl(),input.submit_hl(),input.input_ss(),input.samples())
        return hl_results()[3][input.samples()]

    @render.download(filename=lambda: f'{input.samples()}_smear_summary_table.csv')
    def norm_smear_download():
        req(input.input_hl(),input.exp_hl(),input.submit_hl(),input.input_ss(),input.samples())
        with io.BytesIO() as buf:
            hl_results()[0][input.samples()].to_csv(buf,index=False)
            yield buf.getvalue()
    
    @render.download(filename=lambda: f'{input.exp_hl()}_half_life_table.csv')
    def hl_table_download():
        req(input.input_hl(),input.exp_hl(),input.submit_hl(),input.input_ss())
        with io.BytesIO() as buf:
            hl_results()[1].to_csv(buf,index=False)
            yield buf.getvalue()
    
    @render.download(filename=lambda: f'{input.exp_hl()}_combined_half_life.png')
    def combined_plot_download():
        req(input.input_hl(),input.exp_hl(),input.submit_hl(),input.input_ss())
        plt.figure()
        plot_data=hl_results()[5]
        #plot=hl_results()[2]]
        cm = plt.get_cmap('tab20')
        for c,p in enumerate(plot_data):
            plt.plot(p[0],p[1],label=p[3],color=cm(c))
            plt.scatter(p[0],p[2],marker='o',color=cm(c))
        plt.title("Half-Life of Constructs at 10mM Mg")
        plt.legend()
        plt.xlabel('Time (hr)')
        plt.ylabel('Normalized RNA Integrity (%)')
        with io.BytesIO() as buf:
            plt.savefig(buf,format='png')
            yield buf.getvalue()

    @render.download(filename=lambda: f'{input.samples()}_half_life.png')
    def hl_plots_download():
        req(input.input_hl(),input.exp_hl(),input.submit_hl(),input.input_ss(),input.samples())
        plt.figure()
        plot_data=hl_results()[5]
        #plot=hl_results()[2]
        try:
            for p in plot_data:
                if p[3]==input.samples():
                    hl_plot_data=p
            plt.plot(hl_plot_data[0],hl_plot_data[1],label=hl_plot_data[3])
            plt.scatter(hl_plot_data[0],hl_plot_data[2],marker='o')
            plt.title("Half-Life of "+str(hl_plot_data[3])+" at 10mM Mg")
            plt.legend()
            plt.xlabel('Time (hr)')
            plt.ylabel('Normalized RNA Integrity (%)')
        except:
            print("Error: Sample not found in plot data, please try another sample")
        with io.BytesIO() as buf:
            plt.savefig(buf,format='png')
            yield buf.getvalue()
    
    @render.download(filename='samplesheet_template.csv')
    def example_sheet():
        if input.exp_hl():
            exp_name=input.exp_hl()
        else:
            exp_name="experiment name"
        ss_df=pd.DataFrame(columns=["experiment_name","File","Construct","Lanes","Skip_Lanes"])
        ss_df=ss_df.assign(experiment_name=[exp_name])
        ss_df=ss_df.assign(File=["file_name_smear_results.csv"])
        ss_df=ss_df.assign(Construct=["sample_name (eg. pBPD16)"])
        ss_df=ss_df.assign(Lanes=["A1-A12,B1"])
        ss_df=ss_df.assign(Skip_Lanes=["A1-A3,A9"])
        with io.BytesIO() as buf:
            ss_df.to_csv(buf,index=False)
            yield buf.getvalue()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

app = App(app_ui, server)
