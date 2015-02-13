#!/usr/bin/env python

#Python prelims
import os
import sys

codedir=os.getcwd()+'/trunk/';
pydir=os.getcwd()+'/Python/';
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

#=================================================================================
#SurveySim settings
#---------------------------------------------------------------------------------
print("Welcome to SurveySim");
print("a MCMC-based galaxy evolution fitter and simulator");

#GUI top frame
fields_files='Survey data','SEDs'

#initialize datafiles
seddir=codedir+'templates/'
sedfile=seddir+'my_templates.fits'
dmodelfile=codedir+'model/default_model.fits'
modelfile=codedir+'model/model.fits'
obsdir=codedir+'obs/'
obsfile=obsdir+'wise.fits'
outfile=codedir+'output/output.fits'
fitcode='fitter'

#initialize luminosity function parameters (GUI middle frame)
fields_lf='Phi_0','Lstar_0','alpha','beta','p','q','zcut','c_evol'
value_initial=[-2.20,10.14,0.50,3.00,-4.50,4.50,2.00,0.0]
value_min=[-3.00,9.00,0.00,0.00,-8.00,3.00,0.00,0.0]
value_max=[-1.00,11.00,2.00,5.00,-1.00,6.00,0.00,1.0]
value_fix=[1,1,1,1,0,0,1,1]

#initialize survey parameters (GUI bottom frame)
area=[4.0,4.0,4.0]
band=['Band1','Band2','Band3']
flim=[25.0,20.0,15.0]
f_id=[0,0,0] #placeholder for the filter ids
fields_bands=band[0],band[1],band[2]

#initialize MCMC settings parameters
fields3='zmin','zmax','dz','Runs','Nchain','Tmax','Verbose(1/0)'
zmin=0.01 #Simulation minimum redshift
zmax=5.0 #Simulation maximum redshift
dz=0.05 #Redshift Bin Width
runs=1.e4 #Number of runs
nchain=5 #Chain Number
tmax=100.0 #Starting Anneal Temperature
ann_pct=0.25 #Ideal acceptance Percentage
ann_rng=0.05 #Range to maintain acceptance
conv_con=0.05 #Convergence confidence interval
conv_rma=1.05 #Convergence Rmax Criterion
conv_ste=20 #Steps btw convergence checks
burn_ste=10 #Steps btw anneal calls in burn-in
burnvrun=10 #Ratio of normal to burn-in steps
mesprint=1 #Print Debug MSGs (1=verbose,0=silent)')

#read-in values if model.fits exists
if os.path.isfile(modelfile): # and os.access(modelfile,os.R.OK):
    mfile=fits.open(modelfile)
    mhdr=mfile[0].header
    value_initial[0]=mhdr['PHI0']
    band[0]=mhdr['Band_1']
    band[1]=mhdr['Band_2']
    band[2]=mhdr['Band_3']
    area[0]=mhdr['AREA']
    mfile.close()

#-----------------------------------------------------------------------------------------
#read-in filter transmission curves/ use the FSPS filter set
with open (codedir+'filters/FILTER_LIST','r') as f: 
    flines=f.readlines()

filter_id,filter_choices=[],[]

for fline in flines:
    if fline.strip():
        ncol=len(fline.split())
        if 3<ncol<6:
            tmp1,tmp2,tmp3,tmp4=fline.split()[0:4]
            if (tmp4 != '(from') and (tmp2 != '2MASS') and (tmp2 != 'Steidel') and (tmp2 != 'Stromgren') and (tmp2 != 'Idealized') and (tmp2 != 'SCUBA') and (tmp2 != 'JWST') :
                filter_id.append(int(tmp1));filter_choices.append(tmp2+'_'+tmp3+'_'+tmp4)
            else:
                filter_id.append(int(tmp1));filter_choices.append(tmp2+'_'+tmp3)
        if ncol >= 6:
            tmp1,tmp2,tmp3,tmp4,tmp5=fline.split()[0:5]
            if (tmp4 != '(from') and (tmp3 != 'IRAC') and (tmp3 != 'WFC3') and (tmp2 != '2MASS') and (tmp2 != 'Steidel') and (tmp2 != 'Stromgren') and (tmp2 != 'Idealized') and (tmp2 != 'SCUBA') and (tmp2 != 'JWST') :
                filter_id.append(int(tmp1));filter_choices.append(tmp2+'_'+tmp3+'_'+tmp4)
            else:
                if (tmp3 != 'IRAC') and (tmp3 != 'WFC3'):
                    filter_id.append(int(tmp1));filter_choices.append(tmp2+'_'+tmp3)
            if (tmp3 == 'IRAC') or (tmp3 == 'WFC3'):
                filter_id.append(int(tmp1));filter_choices.append(tmp2+'_'+tmp3+'_'+tmp5)
        if ncol == 3:
            tmp1,tmp2,tmp3=fline.split()[0:3]
            filter_id.append(int(tmp1));filter_choices.append(tmp2+'_'+tmp3)
#-----------------------------------------------------------------------------------------

def update_mfile(modelfile,f_id):
    hdulist=fits.open(modelfile)
    
#read-in filter transmission curves for selected bands
    with open (codedir+'filters/allfilters.dat','r') as f: 
        flines=f.readlines()

    fcount=-1
    lam1,lam2,lam3,trans1,trans2,trans3=[],[],[],[],[],[]
    for fline in flines:
        if fline.strip():
            ncol=len(fline.split())
            if ncol > 2:
                fcount=fcount+1
            if ncol == 2:
                if fcount == f_id[0]: 
                    tmp1,tmp2=fline.split()[0:2]
                    lam1.append(float(tmp1)),trans1.append(float(tmp2))
                if fcount == f_id[1]:
                    tmp1,tmp2=fline.split()[0:2]
                    lam2.append(float(tmp1)),trans2.append(float(tmp2))
                if fcount == f_id[2]:
                    tmp1,tmp2=fline.split()[0:2]
                    lam3.append(float(tmp1)),trans3.append(float(tmp2))
                        
    col1=fits.Column(name='lambda1',format='FLOAT',array=lam1)
    col2=fits.Column(name='transmission1',format='FLOAT',array=trans1)
    col3=fits.Column(name='lambda2',format='FLOAT',array=lam2)
    col4=fits.Column(name='transmission2',format='FLOAT',array=trans2)
    col5=fits.Column(name='lambda3',format='FLOAT',array=lam3)
    col6=fits.Column(name='transmission3',format='FLOAT',array=trans3)

    cols=fits.ColDefs([col1,col2,col3,col4,col5,col6])
    tbhdu=fits.new_table(cols)

    hdr=hdulist[0].header #the header associated with extension=0

    #create/update luminosity function parameters in model file header
    hdr.set('PHI0',value_initial[0],'Luminosity Function Normalization')
    hdr.set('PHI0_FIX',value_fix[0],'Fix Phi0 (Y=1/N=0)')
    hdr.set('PHI0_MIN',value_min[0],'Minimum Phi0 value')
    hdr.set('PHI0_MAX',value_max[0],'Maximum Phi0 value')
    
    hdr.set('L0',value_initial[1],'Luminosity Function knee')
    hdr.set('L0_FIX',value_fix[1],'Fix L0 (Y=1/N=0)')
    hdr.set('L0_MIN',value_min[1],'Minimum L0 value')
    hdr.set('L0_MAX',value_max[1],'Maximum L0 value')
 
    hdr.set('ALPHA',value_initial[2],'Luminosity Function upper slope')
    hdr.set('ALPHA_FI',value_fix[2],'Fix alpha (Y=1/N=0)')
    hdr.set('ALPHA_MI',value_min[2],'Minimum alpha value')
    hdr.set('ALPHA_MA',value_max[2],'Maximum alpha value')

    hdr.set('BETA',value_initial[3],'Luminosity Function lower slope')
    hdr.set('BETA_FIX',value_fix[3],'Fix beta (Y=1/N=0)')
    hdr.set('BETA_MIN',value_min[3],'Minimum beta value')
    hdr.set('BETA_MAX',value_max[3],'Maximum beta value')

    hdr.set('P',value_initial[4],'LF density evolution parameter')
    hdr.set('P_FIX',value_fix[4],'Fix p (Y=1/N=0)')
    hdr.set('P_MIN',value_min[4],'Minimum p value')
    hdr.set('P_MAX',value_max[4],'Maximum p value')

    hdr.set('Q',value_initial[5],'LF luminosity evolution parameter')
    hdr.set('Q_FIX',value_fix[5],'Fix q (Y=1/N=0)')
    hdr.set('Q_MIN',value_min[5],'Minimum q value')
    hdr.set('Q_MAX',value_max[5],'Maximum q value')
    
    hdr.set('ZCUT',value_initial[6],'LF evolution redshift cutoff')
    hdr.set('ZCUT_FIX',value_fix[6],'Fix zcut (Y=1/N=0)')
    hdr.set('ZCUT_MIN',value_min[6],'Minimum zcut value')
    hdr.set('ZCUT_MAX',value_max[6],'Maximum zcut value')

    hdr.set('CEXP',value_initial[6],'Color evolution term')
    hdr.set('CEXP_FIX',value_fix[6],'Fix cexp (Y=1/N=0)')
    hdr.set('CEXP_MIN',value_min[6],'Minimum cexp value')
    hdr.set('CEXP_MAX',value_max[6],'Maximum cexp value')

#====================================================================
# Survey properties 
#--------------------------------------------------------------------   
    hdr.set('AREA',area[0],'Observed Solid Angle of survey')
    hdr.set('Axis1','ColorF1F2','Diagnostic x-axis')
    hdr.set('Axis2','Flux1','Diagnostic y-axis')
    hdr.set('Band_1',band[0],'1st filter name')
    hdr.set('Band_2',band[1],'2nd filter name')
    hdr.set('Band_3',band[2],'3rd filter name')
    hdr.set('Fluxlim1',flim[0],'1st band flux limit')
    hdr.set('Fluxlim2',flim[1],'2nd band flux limit')
    hdr.set('Fluxlim3',flim[2],'3rd band flux limit')

#====================================================================
# Code settings
#---------------------------------------------------------------------  
    hdr.set('ZMIN',zmin,'Simulation minimum redshift')
    hdr.set('ZMAX',zmax,'Simulation maximum redshift')
    hdr.set('DZ',dz,'Redshift Bin Width')
    hdr.set('RUNS',runs,'Number of runs')

    hdr.set('NCHAIN',nchain,'Chain Number')
    hdr.set('TMAX',tmax,'Starting Anneal Temperature')
    hdr.set('ANN_PCT',ann_pct,'Ideal acceptance Percentage')
    hdr.set('ANN_RNG',ann_rng,'Range to maintain acceptance')
    hdr.set('CONV_CON',conv_con,'Convergence confidence interval') #was 0.05 but having trouble converging so set to 0.5 just for now
    hdr.set('CONV_RMA',conv_rma,'Convergence Rmax Criterion')
    hdr.set('CONV_STE',conv_ste,'Steps btw convergence checks')
    hdr.set('BURN_STE',burn_ste,'Steps btw anneal calls in burn-in')
    hdr.set('BURNVRUN',burnvrun,'Ratio of normal to burn-in steps')
    hdr.set('PRINT',mesprint,'Print Debug MSGs (1=verbose,0=silent)')

    hdr['HISTORY']='Last updated on: '+time.strftime("%c") #get current date+time
    hdulist.close

    thdulist=fits.HDUList([hdulist[0],tbhdu])
    if(os.path.isfile(modelfile)):
           print 'Replacing existing model file....'
           os.remove(modelfile)
           thdulist.writeto(modelfile);
    else:
           print 'Writing a new model file....'
           thdulist.writeto(modelfile);
    return

#==========================================================================================
#the GUI definition
#-------------------------------------------------------------------------------------------
class SurveySimGUI:
    def __init__(self, master):
        #local variables to hold the entries in the GUI
        self.v_fixed=[DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar()]
        self.v_min=[DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar()]
        self.v_max=[DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar()]
        self.v_init=[DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar()]
        self.obsfile_set=StringVar()
        self.sedfile_set=StringVar()
        self.f1=DoubleVar()
        self.f2=DoubleVar()
        self.f3=DoubleVar()
        self.db1=StringVar()
        self.db2=StringVar()
        self.db3=StringVar()
        self.defband=[StringVar(),StringVar(),StringVar()]

#variables to hold MCMC settings values
        self.s_set=[DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),IntVar()]
        self.settings_on='no'

        self.master=master
        master.title("SurveySim")
        self.labelframe_top = LabelFrame(master, text="Data files",width=50);
        self.labelframe_top.pack(fill="both",expand="yes");
        self.label_top=Label(self.labelframe_top,text=" ",width=50);
        self.label_top.pack();

        self.labelframe_middle = LabelFrame(master, text="Luminosity Function Parameters",width=50);
        self.labelframe_middle.pack(fill="both", expand="yes");
        self.label_middle=Label(self.labelframe_middle,text="Initial value/Minimum/Maximum/Fixed",width=50);
        self.label_middle.pack();

        self.labelframe_bottom = LabelFrame(master, text="Survey properties",width=50);
        self.labelframe_bottom.pack(fill="both", expand="yes");
        self.label=Label(self.labelframe_bottom,text="Filter/Area [sqdeg]/Flux limit [mJy]",width=50);
        self.label.pack();

#the GUI buttons
        self.settings_button = Button(master, text='Settings',command=self.settings)
        self.settings_button.pack(side=LEFT,padx=5,pady=5)
        self.update_button = Button(master, text='Update',command=self.callback)
        self.update_button.pack(side=LEFT,padx=5,pady=5)
        self.run_button = Button(master, text='Run',command=self.runcode)
        self.run_button.pack(side=LEFT,padx=5,pady=5)
        self.results_button = Button(master, text='Show results',command=self.showresults)
        self.results_button.pack(side=LEFT,padx=5,pady=5)
        self.mcmc_button = Button(master, text='MCMC diagnostics',command=self.mcmcdiag)
        self.mcmc_button.pack(side=LEFT,padx=5,pady=5)
        self.quit_button = Button(master, text='Quit', command=self.quit)
        self.quit_button.pack(side=LEFT, padx=5, pady=5)

# Top frame
        ind=0;
        for field0 in fields_files:
            row=Frame(self.labelframe_top)
            lab=Label(row,width=10,text=field0,anchor='w')
            row.pack(side=TOP,padx=2,pady=5)
            lab.pack(side=LEFT)
            if (ind == 0):
                obsfiles=[]
                for file in os.listdir(obsdir):
                    if file.endswith(".fits"):
                        obsfiles.append(obsdir+file)
                file1=Entry(row,textvariable=obsfile)
                option1=OptionMenu(row,self.obsfile_set,*obsfiles)
                option1.pack(side='left',padx=10,pady=0)
            if (ind == 1):
                sedfiles=[]
                for file in os.listdir(seddir):
                    if file.endswith(".fits"):
                        sedfiles.append(seddir+file)
                file1=Entry(row,textvariable=sedfile)
                option2=OptionMenu(row,self.sedfile_set,*sedfiles)
                option2.pack(side='left',padx=10,pady=0)
            self.obsfile_set.set(obsfile)
            self.sedfile_set.set(sedfile)
            ind=ind+1

#  Middle frame
        ind=0;
        for field in fields_lf:
            row = Frame(self.labelframe_middle)
            lab = Label(row, width=10, text=field, anchor='w')
            row.pack(side=TOP, padx=2, pady=5)
            lab.pack(side=LEFT)

            ent0 = Entry(row,textvariable=self.v_fixed[ind])
            ent0.pack(side=RIGHT)
            self.v_fixed[ind].set(value_fix[ind])
            ent0.pack(side=RIGHT) 

            ent1 = Entry(row,textvariable=self.v_max[ind])
            ent1.pack(side=RIGHT)
            self.v_max[ind].set(value_max[ind])

            ent2 = Entry(row,textvariable=self.v_min[ind])
            ent2.pack(side=RIGHT)
            self.v_min[ind].set(value_min[ind])

            ent3 = Entry(row,textvariable=self.v_init[ind])
            ent3.pack(side=RIGHT)
            self.v_init[ind].set(value_initial[ind])
            ind=ind+1;
        
# Bottom frame
        ind=0;
        for field in fields_bands:
            row = Frame(self.labelframe_bottom)
            lab = Label(row, width=10, anchor='w')
            row.pack(side=TOP, padx=2, pady=5)
            lab.pack(side=LEFT)
            if(ind == 0):
                ent0_1 = Entry(row,textvar=self.f1)
                option1=OptionMenu(row,self.db1,*filter_choices)
                option1.pack(side='left',padx=10,pady=0)
            if(ind == 1):
                ent0_1 = Entry(row,textvar=self.f2) 
                option2=OptionMenu(row,self.db2,*filter_choices)
                option2.pack(side='left',padx=10,pady=0)
            if(ind == 2):
                ent0_1 = Entry(row,textvar=self.f3)
                option3=OptionMenu(row,self.db3,*filter_choices)
                option3.pack(side='left',padx=10,pady=0)
            ent0_1.pack(side=RIGHT) 
            self.db1.set(band[0])
            self.db2.set(band[1])
            self.db3.set(band[2])
            self.f1.set(flim[0])
            self.f2.set(flim[1])
            self.f3.set(flim[2])
            ent2_1 = Entry(row)
            ent2_1.insert(10,area[ind])
            ent2_1.pack(side=RIGHT) 
            ind=ind+1;

    def callback(self):
        ind=0
        for field in fields_lf:
            value_initial[ind]=self.v_init[ind].get()
            value_min[ind]=self.v_min[ind].get()
            value_max[ind]=self.v_max[ind].get()
            value_fix[ind]=self.v_fixed[ind].get()
            ind=ind+1
        obsfile=self.obsfile_set.get()
        sedfile=self.sedfile_set.get()
        band[0]=self.db1.get()
        band[1]=self.db2.get()
        band[2]=self.db3.get()

#update MCMC settings
        if(self.settings_on == 'yes'): #to ensure that this is run only if the settings were touched
            zmin=self.settings.s_set[0].get()
            zmax=self.settings.s_set[1].get()
            dz=self.settings.s_set[2].get()
            runs=self.settings.s_set[3].get()
            nchain=self.settings.s_set[4].get()
            tmax=self.settings.s_set[5].get()
            mesprint=self.settings.s_set[6].get()

        print tmax #just checking
        if (band[0] != 'Band1'):
            f_id=[filter_choices.index(band[0]),filter_choices.index(band[1]),filter_choices.index(band[2])]
        flim[0]=self.f1.get()
        flim[1]=self.f2.get()
        flim[2]=self.f3.get()

        update_mfile(modelfile,f_id)
        return obsfile,sedfile

    def settings(self):
        self.settings = SettingsWindow()
        self.settings_on='yes'

    def showresults(self): #read-in the results from output.fits and show plots
        print("Showing results....")
        hdulist=fits.open(outfile)
#all the extensions headers here are identical -- apply to ccd images
        hdr=hdulist[0].header #the header associated with extension=0
        model_ccd=hdulist[1].data #model color-color distribution
        obs_ccd=hdulist[2].data #observations color-color distribution
        res_ccd=hdulist[0].data #residual color-color distribution
    
        params=hdulist[3].data
        chain=hdulist[4].data
        conv=hdulist[5].data
        res=hdulist[6].data

        counts=hdulist[7].data
    #dnds1=counts.
#read-in the simulated source list
        srcs_hdr=hdulist[3].header
        srcs_data=hdulist[3].data
    #f1=src
    #sim_srcs_table=hdulist[3]# the table with simulated source properties
    #sim_srcs_data=hdulist[3].data
    #the field names are f1,f2,f3,z,m,lum,c7,p0,p1
        sim_srcs_zs=srcs_data.field('z')
        sim_srcs_f1=srcs_data.field('f1')

#seems like these might be in Jy rather than mJy is that right?
        sim_srcs_f1=sim_srcs_f1*1000.0 #turn into mJy
        fbins=np.logspace(0,3,num=20,base=10)
        ds=(fbins[1]-fbins[0])*fbins
        counts1,bin_edges=np.histogram(sim_srcs_f1,bins=fbins)
        ds1=ds[0:18]
        counts=counts1 #[:,0]/ds1[:,0] #*pow(fbins,2.5)
        for i in range(counts1.shape[0]):
            counts[i]=(counts1[i]/ds[i])*(bin_edges[i]**2.5)
            
 #       print bin_edges
 #       print counts

        chain=hdulist[4].data
        p_c1=chain.field('P0')
        q_c1=chain.field('Q0')

        img1=np.zeros((model_ccd.shape[0],model_ccd.shape[1]),dtype=float)
        img2=np.zeros((model_ccd.shape[0],model_ccd.shape[1]),dtype=float)
        img3=np.zeros((model_ccd.shape[0],model_ccd.shape[1]),dtype=float)
            
        img1[:,:]=img_scale.linear(model_ccd,scale_min=0,scale_max=100.0)
        img2[:,:]=img_scale.linear(obs_ccd,scale_min=0,scale_max=100.0)
        img3[:,:]=img_scale.linear(res_ccd,scale_min=0,scale_max=100.0)

    #make actual figure
        fig=plt.figure(1)
        a1=fig.add_subplot(2,3,1)
        a1.plot(p_c1,q_c1,'k.')
        a1.set_xlabel('P')
        a1.axis([value_min[4],value_max[4],value_min[5],value_max[5]])
        a1.set_ylabel('Q')

        a2=fig.add_subplot(2,3,2)
        a2.set_title('Redshifts')
        zbins=range(25)
        zbins=np.divide(zbins,5.0)
        hist(sim_srcs_zs,bins=zbins,histtype='step',color='black')
        a2.set_xlabel('z')
        a2.set_ylabel('N(z)')
        a2.set_xlim(0, 5)

        a3=fig.add_subplot(2,3,3)
        a3.set_title('Counts')
        a3.set_xscale('log')
        a3.set_yscale('log')
    
        nbins=size(counts)
        bins=bin_edges[0:nbins]
        plt.plot(bins,counts,'+')
        a3.set_xlabel('S [mJy]')
        a3.set_ylabel('(dN/dS)*S^2.5')
        a3.set_xlim(20,1000)
        a3.set_ylim(200,8000000)

        a4=fig.add_subplot(2,3,4) 
        imgplot=plt.imshow(img1)
        a4.set_title('Model')

        a5=fig.add_subplot(2,3,5)
        imgplot=plt.imshow(img2)
        a5.set_title('Observations')

        a6=fig.add_subplot(2,3,6)
        imgplot=plt.imshow(img3)
        a6.set_title('Residual')

        plt.tight_layout(1)
        plt.show(1)
        return

    def mcmcdiag(self): #show chain behavior
        print("MCMC diagnostics....")
        hdulist=fits.open(outfile)
#all the extensions headers here are identical -- apply to ccd images
        hdr=hdulist[0].header #the header associated with extension=0
        chain=hdulist[4].data
        chi2_c1=chain.field('CHISQ0')
        chi2_c2=chain.field('CHISQ1')
        chi2_c3=chain.field('CHISQ2')
        chi2_c4=chain.field('CHISQ3')
        chi2_c5=chain.field('CHISQ4')
    
        csize=size(chi2_c1)
        a=np.arange(csize)
        plt.figure(2)
        plt.plot(a,chi2_c1,'r-',lw=2)
        plt.xlabel('step')
        plt.ylabel('chi2')
        plt.plot(a,chi2_c2,'k-',lw=0.7,color='0.5')
        plt.plot(a,chi2_c3,'k-',lw=0.7,color='0.5')
        plt.plot(a,chi2_c4,'k-',lw=0.7,color='0.5')
        plt.plot(a,chi2_c5,'k-',lw=0.7,color='0.5')
        plt.show(2)
        return

    def runcode(self):
        print("Runnning fitter code....")
        obsfile=self.obsfile_set.get()
        sedfile=self.sedfile_set.get()
        os.system(fitcode+' '+obsfile+' '+modelfile+' '+sedfile+' '+outfile)

    def quit(self):
        self.master.destroy()

class SettingsWindow(Frame):     
    def __init__(self):
        new =tk.Frame.__init__(self)
        new = Toplevel(self)
        new.title("MCMC Settings")

        #variables to hold MCMC settings values
        self.s_set=[DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),IntVar()]

        ind=0;
        for field3_i in fields3:
            row = Frame(new)
            lab = Label(row, width=10, text=field3_i,anchor='w')
            row.pack(side=TOP, padx=2, pady=5)
            lab.pack(side=LEFT)
            e = Entry(row,textvariable=self.s_set[ind])
            e.pack(side=RIGHT)
            if(ind == 0): 
                self.s_set[ind].set(zmin)
            if(ind == 1):
                self.s_set[ind].set(zmax)
            if(ind == 2):
                self.s_set[ind].set(dz)
            if(ind == 3):
                self.s_set[ind].set(runs)
            if(ind == 4):
                self.s_set[ind].set(nchain)
            if(ind == 5):
                self.s_set[ind].set(tmax)
            if(ind == 6):
                self.s_set[ind].set(mesprint)
            ind=ind+1;


def main():
    root=Tk()   
    mygui=SurveySimGUI(root)
    root.mainloop() 

if __name__ == '__main__':
    main()

