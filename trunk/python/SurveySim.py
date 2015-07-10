#!/usr/bin/env python

#Python prelims
import os
import sys, getopt
import datetime

thisdir=os.getcwd()+'/'
if (os.getenv("SURVEYSIMPATH") != None):
    codedir=os.getenv("SURVEYSIMPATH")+'/'
else:
    codedir='/usr/local/surveysim/'
pydir=codedir+'python'
#ensure this is in the path
sys.path.append(pydir)

import time
import pyfits as fits
import matplotlib.pyplot as plt
import numpy as np
import img_scale
from pylab import *
import math
from functools import partial
import Tkinter as tk
from Tkinter import *
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from ModelFile import *
from OutputFile import *
from filters import read_filters

filter_id,filter_choices=read_filters();

mod=ModelFile()

#=================================================================================
#SurveySim settings
#---------------------------------------------------------------------------------
print("Welcome to SurveySim");
print("a MCMC-based galaxy evolution fitter and simulator");

#GUI top frame
fields_files='Survey data','SEDs'

#initialize datafiles
seddir=codedir+'templates/'
sedfile=seddir+'Kirkpatrick2015_templates.fits'
dmodelfile=codedir+'model/default_model.fits'
modelfile=thisdir+'model.fits'
obsdir=codedir+'obs/'
obsfile=obsdir+'wise.fits'
outdir=os.getcwd()+'/OUTPUT/'
outfile=os.getcwd()+'/output.fits'
fitcode='/usr/local/bin/SurveySim'

try:
    opts, args = getopt.getopt(sys.argv[1:],"hi:o:",["ifile=","ofile="])
except getopt.GetoptError:
    print 'SurveySim.py -i <modelfile> -o <outputfile>'
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print 'SurveySim.py -i <modelfile> -o <outputfile>'
        sys.exit()
    elif opt in ("-i", "--ifile"):
        dmodelfile = arg
    elif opt in ("-o", "--ofile"):
        outfile = arg
print 'Input file is "', dmodelfile
print 'Output file is "', outfile

mod.load(dmodelfile)

if(not os.path.exists(outdir)):
    os.makedirs(outdir)

lfform_tmp=mod.survey['lfForm']
if(lfform_tmp == 0):
    lfform='DoublePowerLaw'
if(lfform_tmp == 1):
    lfform='ModifiedSchecter'
if(lfform_tmp == 2):
    lfform='Schecter'
if(lfform_tmp == 3):
    lfform='DoubleExponential'

#initialize luminosity function parameters (GUI middle frame)
fields_lf=mod.params.keys()
npar=len(fields_lf)
value_initial,value_min,value_max,value_fix=[],[],[],[]
for lffield in fields_lf:
    value_initial.append(mod.params[lffield].value)
    value_min.append(mod.params[lffield].pmin)
    value_max.append(mod.params[lffield].pmax)
    value_fix.append(mod.params[lffield].fixed)

#initialize survey parameters (GUI bottom frame)
axes=mod.axis1,mod.axis2

area=mod.survey['area']
band=[getFilterName(mod.filters[0].fid),getFilterName(mod.filters[1].fid),getFilterName(mod.filters[2].fid)]
band_units=[mod.filters[0].unit,mod.filters[1].unit,mod.filters[2].unit]
flim=[mod.filters[0].limit,mod.filters[1].limit,mod.filters[2].limit]
ferr=[mod.filters[0].err,mod.filters[1].err,mod.filters[2].err]

f_id=[0,0,0] #placeholder for the filter ids
fields_bands=band[0],band[1],band[2]

#initialize MCMC settings parameters
#fields3='zmin','zmax','dz','Runs','Nchain','verbosity','Tmax','conv_rmax','conv_step','conv_ci'
fields3='zmin','zmax','dz','Runs','Nchain','verbosity','Tmax','conv_ci'

zmin1=mod.redshift['min'] #Simulation minimum redshift
zmax1=mod.redshift['max'] #Simulation maximum redshift
dz1=mod.redshift['delta'] #Redshift Bin Width

runs1=mod.settings['runs'] #Number of runs
nchain1=mod.settings['nchain'] #Chain Number
mesprint1=mod.settings['verbosity'] #Print Debug MSGs (0=silent,1=critical, 2=info,3=debug)')

tmax1=mod.annealing['temp'] #Starting Anneal Temperature
#ann_pct=mod.annealing['ideal_pct'] #Ideal acceptance Percentage
#ann_rng=mod.annealing['range'] #Range to maintain acceptance
#burn_ste=mod.annealing['burn_step'] #Steps btw anneal calls in burn-in
#burnvrun=10 #Ratio of normal to burn-in steps

#conv_rma1=mod.convergence['r_max'] #Convergence Rmax Criterion
#conv_ste1=mod.convergence['step'] #Steps btw convergence checks
conv_con1=mod.convergence['CI'] #Convergence confidence interval

#-----------------------------------------------------------------------------------------

#==========================================================================================
#the GUI definition
#-------------------------------------------------------------------------------------------
class SurveySimGUI:
    def __init__(self, master):
        #local variables to hold the entries in the GUI
        self.v_fixed=[IntVar(),IntVar(),IntVar(),IntVar(),IntVar(),IntVar(),IntVar(),IntVar(),IntVar(),IntVar(),IntVar(),IntVar(),IntVar()] 
        self.v_min=[DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar()] 
        self.v_max=[DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar()] 
        self.v_init=[DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar()]
        self.obsfile_set=StringVar()
        self.lfform_set=StringVar()
        self.sedfile_set=StringVar()
        self.colsel=StringVar()
        self.area=DoubleVar()
        self.fitaxes=[StringVar(),StringVar()]
        self.limits=[DoubleVar(),DoubleVar(),DoubleVar()]
        self.units=[StringVar(),StringVar(),StringVar()]
        self.bands=[StringVar(),StringVar(),StringVar()]
        self.defband=[StringVar(),StringVar(),StringVar()]
        self.settings_on='no' #a switch to say whether or not the SettingsWindow was used

        self.master=master
        master.title("SurveySim")
        self.labelframe_top = LabelFrame(master, text="Data files",bg='grey')
        self.labelframe_top.grid(columnspan=3,sticky=W+E)

        self.labelframe_lf = LabelFrame(master, text="Luminosity Function Parameters",bg='light blue') 
        self.labelframe_lf.grid(column=0,row=1,rowspan=2)
        self.label_lf=Label(self.labelframe_lf,text="Initial/Min/Max/Fix",bg='light blue')
        self.label_lf.grid(in_=self.labelframe_lf,row=2,column=1)

        self.labelframe_survey = LabelFrame(master, text="Survey fitting properties",bg='pink') 
        self.labelframe_survey.grid(column=1,row=1,sticky=W+N+S)
        self.label_survey=Label(self.labelframe_survey,text="Filter/Limit/Units",bg='pink') 
        self.label_survey.grid(in_=self.labelframe_survey,row=1,column=1)

        self.labelframe_seds=LabelFrame(master,text="SEDs",bg='green')
        self.labelframe_seds.grid(column=2,row=1)

#the GUI buttons
        self.labelframe_buttons=LabelFrame(master)
        self.labelframe_buttons.grid(column=1,columnspan=2,row=2,sticky=S+E+W)
        row=Frame(self.labelframe_buttons)
        row.pack(side=LEFT)
        self.settings_button = Button(row,text='Settings',command=self.settings)
        self.settings_button.pack(side=LEFT,padx=5,pady=5)
        self.update_button = Button(self.labelframe_buttons, text='Update',command=self.update_mfile)
        self.update_button.pack(side=LEFT,padx=5,pady=5)
        self.run_button = Button(self.labelframe_buttons, text='Run',command=self.runcode)
        self.run_button.pack(side=LEFT,padx=5,pady=5)
        self.results_button = Button(self.labelframe_buttons, text='Show results',command=self.showresults)
        self.results_button.pack(side=LEFT,padx=5,pady=5)
        self.mcmc_button = Button(self.labelframe_buttons, text='MCMC diagnostics',command=self.mcmcdiag)
        self.mcmc_button.pack(side=LEFT,padx=5,pady=5)
        self.quit_button = Button(self.labelframe_buttons, text='Quit', command=self.quit)
        self.quit_button.pack(side=LEFT, padx=5, pady=5)

# Data files frame
        ind=0;
        for field0 in fields_files:
            row=Frame(self.labelframe_top,bg='grey')
            lab=Label(row,text=field0,bg='grey')
            lab.pack(side=LEFT)
            row.pack(side=TOP,padx=1,pady=5)
            if (ind == 0):
                obsfiles=[]
                for file in os.listdir(obsdir):
                    if file.endswith(".fits"):
                        obsfiles.append(obsdir+file)
                file1=Entry(row,textvariable=obsfile)
                option1=OptionMenu(row,self.obsfile_set,*obsfiles)
                option1.pack(side='left',padx=0,pady=0)
            if (ind == 1):
                sedfiles=[]
                for file in os.listdir(seddir):
                    if file.endswith(".fits"):
                        sedfiles.append(seddir+file)
                file1=Entry(row,textvariable=sedfile)
                option2=OptionMenu(row,self.sedfile_set,*sedfiles)
                option2.pack(side='left',padx=0,pady=0)
            self.obsfile_set.set(obsfile)
            self.sedfile_set.set(sedfile)
            ind=ind+1

#  LF frame
        ind=0;
        lab=Label(self.labelframe_lf,text='LF form',bg='light blue')
        lab.grid(in_=self.labelframe_lf,row=ind+1,column=0)
        lfforms=['DoublePowerLaw','ModifiedSchecter','Schecter']
        lf_option=OptionMenu(self.labelframe_lf,self.lfform_set,*lfforms)
        self.lfform_set.set(lfform)
        lf_option.grid(in_=self.labelframe_lf,row=ind+1,column=1)
 
        for field in fields_lf:
            lab = Label(self.labelframe_lf, text=field,bg='light blue')
            lab.grid(in_=self.labelframe_lf,row=ind+3,column=0)

            ent3 = Entry(self.labelframe_lf,textvariable=self.v_init[ind],width=5)
            self.v_init[ind].set(mod.params[field].value)
            ent3.grid(in_=self.labelframe_lf,row=ind+3,column=1)

            ent2 = Entry(self.labelframe_lf,textvariable=self.v_min[ind],width=5)
            self.v_min[ind].set(value_min[ind])
            ent2.grid(in_=self.labelframe_lf,row=ind+3,column=2)

            ent1 = Entry(self.labelframe_lf,textvariable=self.v_max[ind],width=5)
            self.v_max[ind].set(value_max[ind])
            ent1.grid(in_=self.labelframe_lf,row=ind+3,column=3)

            ent0 = Entry(self.labelframe_lf,textvariable=self.v_fixed[ind],width=1)
            self.v_fixed[ind].set(value_fix[ind])
            ent0.grid(in_=self.labelframe_lf,row=ind+3,column=4)
            ind=ind+1;
        
# Survey frame
        ind=0;
        for field in fields_bands:
            option1=OptionMenu(self.labelframe_survey,self.bands[ind],*filter_choices)
            self.bands[ind].set(band[ind])
            option1.grid(in_=self.labelframe_survey,row=ind+2,column=0)
            ent0_1=Entry(self.labelframe_survey,textvar=self.limits[ind],width=5)
            ent1_1=Entry(self.labelframe_survey,textvar=self.units[ind],width=5)
            self.limits[ind].set(flim[ind])
            self.units[ind].set(band_units[ind])
            ent0_1.grid(in_=self.labelframe_survey,row=ind+2,column=1)
            ent1_1.grid(in_=self.labelframe_survey,row=ind+2,column=2)
            ind=ind+1;

        lab = Label(self.labelframe_survey, width=10, text='Color cut:', anchor='w',bg='pink')
        lab.grid(in_=self.labelframe_survey,row=ind+2,column=0)
        ent2_0=Entry(self.labelframe_survey,textvar=self.colsel,width=8)
        self.colsel.set('None')
        ent2_0.grid(in_=self.labelframe_survey,row=ind+2,column=1,pady=2)
        self.info_button = Button(self.labelframe_survey, text='?',command=self.colsel_info)
        self.info_button.grid(row=ind+2,column=2)

        lab = Label(self.labelframe_survey, width=10, text='Area[sq.deg.]=', anchor='w',bg='pink')
        lab.grid(in_=self.labelframe_survey,row=ind+3,column=0,pady=2)
        ent2_0=Entry(self.labelframe_survey,textvar=self.area,width=8)
        self.area.set(area)
        ent2_0.grid(in_=self.labelframe_survey,row=ind+3,column=1,pady=2)

        lab = Label(self.labelframe_survey, width=6, text='AXIS1=', anchor='w',bg='pink')
        lab.grid(in_=self.labelframe_survey,row=ind+4,column=0,pady=2)
        ent2_1=Entry(self.labelframe_survey,textvar=self.fitaxes[0],width=8)
        self.fitaxes[0].set(axes[0])
        ent2_1.grid(in_=self.labelframe_survey,row=ind+4,column=1,pady=2)

        lab = Label(self.labelframe_survey, width=6, text='AXIS2=', anchor='w',bg='pink')
        lab.grid(in_=self.labelframe_survey,row=ind+5,column=0,pady=2)
        ent3_1=Entry(self.labelframe_survey,textvar=self.fitaxes[1],width=8)
        self.fitaxes[1].set(axes[1])
        ent3_1.grid(in_=self.labelframe_survey,row=ind+5,column=1,pady=2)

#Plot_SED_templates frame
        fig=plt.Figure(figsize=(4,3)) 
        x = np.arange(0, 2*np.pi, 0.01)        # x-array
        canvas = FigureCanvasTkAgg(fig, master=self.labelframe_seds)
        canvas.get_tk_widget().grid(in_=self.labelframe_seds,column=0,row=0,sticky=N+S)
        ax = fig.add_subplot(111)
        sfile=fits.open(self.sedfile_set.get())
        shdr=sfile[1].header
        nlam=shdr["NAXIS1"]
        ntmp=shdr["TFIELDS"]-1 #number of templates -1 is because the first field is the lambda array
        seds=sfile[1].data #the first extension contains the data
        lam=seds['lambda'] #seds[0,0,0,]
        for i in range(1,ntmp):
            fnu=seds['lnu'+str(i)]
            ax.plot(lam,fnu)

        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('wavelength [um]')
        ax.set_ylabel('Fnu [W/Hz]')
        show()
        sfile.close()

    def colsel_info(self):
        self.colselinfo = ColInfoWindow()

    def settings(self):
        self.settings = SettingsWindow()
        self.settings_on='yes'

    def showresults(self): #read-in the results from output.fits and show plots
        print("Showing results....")
        output=OutputFile(outfile)
        output.show(block=False)
        return

    def mcmcdiag(self): #show chain behavior
        print("MCMC diagnostics....")
        output=OutputFile(outfile)
        if(output.fit()):
            output.MCMC.showFit(block=False)
            output.MCMC.showChains(block=False)
        return

    def update_mfile(self):
        if(self.lfform_set.get() == "DoublePowerLaw"):
            mod.survey['lfForm']=0
        if(self.lfform_set.get() == "ModifiedSchecter"):
            mod.survey['lfForm']=1
        if(self.lfform_set.get() == "Schecter"):
            mod.survey['lfForm']=2
        mod.survey['area']=self.area.get()
        ind=0
        lfform=self.lfform_set.get()
        for field in fields_lf:
            mod.params[field].value=self.v_init[ind].get()
            mod.params[field].pmin=self.v_min[ind].get()
            mod.params[field].pmax=self.v_max[ind].get()
            mod.params[field].fixed=self.v_fixed[ind].get()
            ind=ind+1

        #mod.zbc=self.zbc.get()
        mod.zbc=2.0
        mod.area=self.area.get()
        mod.axis1=self.fitaxes[0].get()
        mod.axis2=self.fitaxes[1].get() 
        ind=0
        for field in fields_bands:
            mod.filters[ind].setID(self.bands[ind].get())
            mod.filters[ind].fid=mod.filters[ind].fid
            mod.filters[ind].limit=self.limits[ind].get()
            mod.filters[ind].error=mod.filters[ind].limit/3.0 #assume 3sigma det.limit, we do not use this at present
            mod.filters[ind].unit=self.units[ind].get()
            ind=ind+1

        #Update MCMC settings
        if(self.settings_on == 'yes'): #to ensure that this is run only if the settings were touched
            mod.redshift['min']=self.settings.s_set[0].get()
            mod.redshift['max']=self.settings.s_set[1].get()
            mod.redshift['delta']=self.settings.s_set[2].get()

            mod.settings['runs']=self.settings.s_set[3].get()
            mod.settings['nchain']=self.settings.s_set[4].get()

            mod.settings['verbosity']=self.settings.s_set[5].get()
            mod.annealing['temp']=self.settings.s_set[6].get() 
            #mod.annealing['learningRate']=need to add this!
            
            #mod.convergence['r_max']=self.settings.s_set[7].get()
            #mod.convergence['step']=self.settings.s_set[8].get()
            mod.convergence['CI']=self.settings.s_set[7].get()

        if (band[0] != 'Band1'):
            mod.filters[0].filterIDs=filter_choices.index(band[0])
            mod.filters[1].filterIDs=filter_choices.index(band[1])
            mod.filters[2].filterIDs=filter_choices.index(band[2])
    
        mod.write(modelfile)

    def plot_seds(self):
        fig=plt.Figure(figsize=(3.0,2.0)) #,dpi=100, facecolor='w')
        x = np.arange(0, 2*np.pi, 0.01)        # x-array
        canvas = FigureCanvasTkAgg(fig, master=self.labelframe_seds)
        canvas.get_tk_widget().grid(in_=self.labelframe_seds,column=0,row=0,sticky=N+S)
        ax = fig.add_subplot(111)
        sedfile=self.sedfile_set.get()
        sfile=fits.open(sedfile)
        shdr=sfile[0].header
        seds=sfile[0].data #the first extension contains the data
        lam=seds[0,0,0,]
        fnu1=seds[0,0,1,]
        for i in range(1,ntmp):
            fnu=seds[0,0,i,]
            ax.plot(lam,fnu)

        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('wavelength [um]')
        ax.set_ylabel('Fnu [W/Hz]')
        show()
        sfile.close()
        return


    def runcode(self):
        print("Runnning fitter code....")
        obsfile=self.obsfile_set.get()
        sedfile=self.sedfile_set.get()
#        os.system(fitcode+' '+obsfile+' '+modelfile+' '+sedfile+' '+outfile)
        os.system(fitcode+' '+modelfile+' '+sedfile+' '+obsfile+' -o '+outfile)

    def quit(self):
        self.master.destroy()

class ColInfoWindow(Frame):     
    def __init__(self):
        popup =tk.Frame.__init__(self)
        popup = Toplevel(self)
        about_message='E.g.: mag1-mag2>2 [AB]  or    colF1F2>2          [Fnu spectral index]'
        msg = Message(popup, text=about_message,width=160)
        msg.pack()
        button = Button(popup, text="Okay", command=popup.destroy)
        button.pack()

class SettingsWindow(Frame):     
    def __init__(self):
        new =tk.Frame.__init__(self)
        new = Toplevel(self)
        new.title("MCMC Settings")

        #variables to hold MCMC settings values
        self.s_set=[DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar()]

        ind=0;
        for field3_i in fields3:
            row = Frame(new)
            lab = Label(row, width=10, text=field3_i,anchor='w')
            row.pack(side=TOP, padx=2, pady=5)
            lab.pack(side=LEFT)
            e = Entry(row,textvariable=self.s_set[ind])
            e.pack(side=RIGHT)
            if(ind == 0): 
                self.s_set[ind].set(zmin1)
            if(ind == 1):
                self.s_set[ind].set(zmax1)
            if(ind == 2):
                self.s_set[ind].set(dz1)
            if(ind == 3):
                self.s_set[ind].set(runs1)
            if(ind == 4):
                self.s_set[ind].set(nchain1)
            if(ind == 5):
                self.s_set[ind].set(mesprint1)
            if(ind == 6):
                self.s_set[ind].set(tmax1)
 #           if(ind == 7):
 #               self.s_set[ind].set(conv_rma1)
 #           if(ind == 8):
 #               self.s_set[ind].set(conv_ste1)
            if(ind == 7):
                self.s_set[ind].set(conv_con1)
            ind=ind+1;

def main():
    root=Tk()   
    mygui=SurveySimGUI(root)
    root.mainloop() 

if __name__ == '__main__':
    main()

