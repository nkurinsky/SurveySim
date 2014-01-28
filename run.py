#!/usr/bin/python

print("Welcome to SurveySim");
print("a MCMC-based galaxy evolution fitter and simulator");

import time
from astropy.io import fits
import matplotlib.pyplot as plt
import os
import numpy as np
import img_scale
from pylab import *
import math

codedir=os.getcwd();
#if want to be able to run this from wherever, should simply change the codedir 
sedfile=codedir+'/Widget/sf_templates.fits'
modelfile=codedir+'/Widget/model.fits'
obsfile=codedir+'/Widget/observation.fits'
outfile=codedir+'/Widget/output.fits'

def runcode():
    print("Runnning fitter code....")
    os.system(codedir+"/Fitting/v0/fitter"+' '+obsfile+' '+modelfile+' '+sedfile+' '+outfile)

def showresults(): #read-in the results from output.fits and show plots
    print("Showing results....")
    hdulist=fits.open(outfile)
#all the extensions headers here are identical -- apply to ccd images
    hdr=hdulist[0].header #the header associated with extension=0
    model_ccd=hdulist[1].data #model color-color distribution
    obs_ccd=hdulist[2].data #observations color-color distribution
    res_ccd=hdulist[0].data #residual color-color distribution
    sim_srcs_table=hdulist[3]# the table with simulated source properties
    sim_srcs_table_data=hdulist[3].data
    #the field names are f1,f2,f3,z,m,lum,c7,p0,p1
    sim_srcs_zs=sim_srcs_table_data.field('z')
    sim_srcs_f1=sim_srcs_table_data.field('f1')
#    mctable=hdulist[3]#.columns #MCMC chain
#    mctable.columns
#    mctable.shape
    #cols=hdulist[4].columns
#    partable=hdulist[4].data #parameters table
#    print mctable.shape,partable.shape
#    print model_ccd.dtype.name
    img1=np.zeros((model_ccd.shape[0],model_ccd.shape[1]),dtype=float)
    img2=np.zeros((model_ccd.shape[0],model_ccd.shape[1]),dtype=float)
    img3=np.zeros((model_ccd.shape[0],model_ccd.shape[1]),dtype=float)

    img1[:,:]=img_scale.linear(model_ccd,scale_min=0,scale_max=100.0)
    img2[:,:]=img_scale.linear(obs_ccd,scale_min=0,scale_max=100.0)
    img3[:,:]=img_scale.linear(res_ccd,scale_min=0,scale_max=100.0)

    #make actual figure
    fig=plt.figure()
    a1=fig.add_subplot(2,3,1)
    a1.set_title('Redshifts')
    zbins=range(25)
    zbins=np.divide(zbins,5.0)
    hist(sim_srcs_zs,bins=zbins,histtype='step',color='black')

    a1.set_xlabel('z')
    a1.set_ylabel('N(z)')
    a1.set_xlim(0, 5)

    a2=fig.add_subplot(2,3,2)
    a2.set_title('250um counts')
    a2.set_xscale('log')
    a2.set_yscale('log')
    a2.set_xlabel('S [mJy]')
    a2.set_ylabel('(dN/dS)*S^2.5')
    a2.set_xlim(20,1000)
    a2.set_ylim(200,50000)
    cfile=open(codedir+'/Widget/counts_clements10.dat','r')
#    flux250=np.empty([15],dtype=float)
    flux250=[]
    dnds250=[]
    dnds_err250=[]
    ind=0
    for line in cfile:
        if((ind >= 2) and (ind <=17)):
            line=line.strip()
            columns=line.split()
#            flux250[ind-2]=columns[0]
#            dnds250[ind-2]=columns[5]
#            dnds_err250[ind-2]=columns[6]
            flux250.append(columns[0])
            dnds250.append(columns[5])
            dnds_err250.append(columns[6])
            #print flux[ind],dnds[ind]
        ind=ind+1
#    counts=cfile.read()
#    print counts.shape
    flux250=np.asarray(flux250,dtype=float)
    dnds250=np.asarray(dnds250,dtype=float)
    dnds_err250=np.asarray(dnds_err250,dtype=float)
    plt.errorbar(flux250,dnds250,dnds_err250,capsize=0,ls='none',marker='o',color='black',elinewidth=2)
    cfile.close()

    fbins=flux250
    dnds_sim=fbins*0.0
    ind=0
    for fbin in fbins:
        if(ind < 15):
            a=find((sim_srcs_f1 <=fbins[ind]))
            tmp=sim_srcs_f1[a]
            b=find(tmp > fbins[ind+1])
            tmp=fbins[ind]
            dnds_sim[ind]=(tmp**(2.5))*len(b)/area[0]
            #print fbins[ind],len(b),dnds_sim[ind]
        ind=ind+1
    print dnds_sim
#    plt.errorbar(flux250,dnds_sim,dnds_sim*0.0,color='black')
    plt.plot(flux250,dnds_sim,'-',color='black')
#    hist(sim_srcs_f1,bins=fbins,histtype='step',color='black')


    a3=fig.add_subplot(2,3,3)
    a3.set_xlabel('x')
    a3.set_ylabel('y')
    #plot(x,y)
    a3.set_xlim(0,2)
    a3.set_ylim(-2, 2)
    a3.set_aspect('equal', 'datalim')

    a4=fig.add_subplot(2,3,4) 
    imgplot=plt.imshow(img1)
    a4.set_title('Model')

    a5=fig.add_subplot(2,3,5)
    imgplot=plt.imshow(img2)
    a5.set_title('Observations')

    a6=fig.add_subplot(2,3,6)
    imgplot=plt.imshow(img3)
    a6.set_title('Residual')

    plt.show()
    return

defaults=raw_input("Do you wish to see/edit the default file settings (y/n)?");
if (defaults == 'y'):
    print("The SED templates are defined in:");
    print sedfile;
    change=raw_input("Accept (y/n)?")
    if(change == 'n'):
        sedfile=raw_input("New SED template file:");
    print("The observations to be fit are defined in:");
    print obsfile;
    change=raw_input("Accept (y/n)?")
    if(change == 'n'):
        obsfile=raw_input("New observations file:");
    print("The output is defined in:");
    print outfile;
    change=raw_input("Accept (y/n)?")
    if(change == 'n'):
        outfile=raw_input("New output file:");

mdefaults=raw_input("Do you wish to see/edit the default model settings (y/n)?");

#initialize luminosity function parameters
fields='Phi_0','Lstar_0','alpha','beta','p','q','zcut','c_evol'
value_initial=[-2.20,10.14,0.50,3.00,-4.50,4.50,2.00,0.0]
value_min=[-3.00,9.00,0.00,0.00,-8.00,3.00,0.00,0.0]
value_max=[-1.00,11.00,2.00,5.00,-1.00,6.00,0.00,1.0]
value_fix=[1,1,1,1,0,0,1,1]
#not sure why we need this here, I feel that the sigma of the proposal distribution shouldn't be user controllable
value_dp=[0.0,0.0,0.0,0.0,0.3,0.1,0.0,0.0]

#initialize survey parameters
area=[4.0,4.0,4.0]
wavelength=[250,350,500]
flim=[25.0,20.0,15.0]

def update_mfile(modelfile):
    hdu1=fits.open(modelfile,mode='update') 
    hdr=hdu1[0].header #the header associated with extension=0

    #create/update luminosity function parameters in model file header
    hdr.set('PHI0',value_initial[0],'Luminosity Function Normalization')
    hdr.set('PHI0_FIX',value_fix[0],'Fix Phi0 (Y=1/N=0)')
    hdr.set('PHI0_MIN',value_min[0],'Minimum Phi0 value')
    hdr.set('PHI0_MAX',value_max[0],'Maximum Phi0 value')
    hdr.set('PHI0_DP',value_dp[0],'sigma phi0')
    
    hdr.set('L0',value_initial[1],'Luminosity Function knee')
    hdr.set('L0_FIX',value_fix[1],'Fix L0 (Y=1/N=0)')
    hdr.set('L0_MIN',value_min[1],'Minimum L0 value')
    hdr.set('L0_MAX',value_max[1],'Maximum L0 value')
    hdr.set('L0_DP',value_dp[0],'sigma L0')
 
    hdr.set('ALPHA',value_initial[2],'Luminosity Function upper slope')
    hdr.set('ALPHA_FI',value_fix[2],'Fix alpha (Y=1/N=0)')
    hdr.set('ALPHA_MI',value_min[2],'Minimum alpha value')
    hdr.set('ALPHA_MA',value_max[2],'Maximum alpha value')
    hdr.set('ALPHA_DP',value_dp[2],'sigma alpha')

    hdr.set('BETA',value_initial[3],'Luminosity Function lower slope')
    hdr.set('BETA_FIX',value_fix[3],'Fix beta (Y=1/N=0)')
    hdr.set('BETA_MIN',value_min[3],'Minimum beta value')
    hdr.set('BETA_MAX',value_max[3],'Maximum beta value')
    hdr.set('BETA_DP',value_dp[3],'sigma beta')

    hdr.set('P',value_initial[4],'LF density evolution parameter')
    hdr.set('P_FIX',value_fix[4],'Fix p (Y=1/N=0)')
    hdr.set('P_MIN',value_min[4],'Minimum p value')
    hdr.set('P_MAX',value_max[4],'Maximum p value')
    hdr.set('P_DP',value_dp[4],'sigma p')

    hdr.set('Q',value_initial[5],'LF luminosity evolution parameter')
    hdr.set('Q_FIX',value_fix[5],'Fix q (Y=1/N=0)')
    hdr.set('Q_MIN',value_min[5],'Minimum q value')
    hdr.set('Q_MAX',value_max[5],'Maximum q value')
    hdr.set('Q_DP',value_dp[5],'sigma q')
    
    hdr.set('ZCUT',value_initial[6],'LF evolution redshift cutoff')
    hdr.set('ZCUT_FIX',value_fix[6],'Fix zcut (Y=1/N=0)')
    hdr.set('ZCUT_MIN',value_min[6],'Minimum zcut value')
    hdr.set('ZCUT_MAX',value_max[6],'Maximum zcut value')
    hdr.set('ZCUT_DP',value_max[6],'sigma zcut')

    hdr.set('CEXP',value_initial[6],'Color evolution term')
    hdr.set('CEXP_FIX',value_fix[6],'Fix cexp (Y=1/N=0)')
    hdr.set('CEXP_MIN',value_min[6],'Minimum cexp value')
    hdr.set('CEXP_MAX',value_max[6],'Maximum cexp value')
    hdr.set('CEXP_DP',value_dp[6],'sigma cexp')
    
    hdr.set('AREA',area[0],'Observed Solid Angle of survey')
  
    #these aren't currently user controllable in widget, but can change them right here at own risk  
    hdr.set('ZMIN',0.0,'Simulation minimum redshift')
    hdr.set('ZMAX',5.0,'Simulation maximum redshift')
    hdr.set('DZ',0.1,'Redshift Bin Width')
    hdr.set('RUNS',1.e3,'Number of runs')

    hdr.set('NCHAIN',5,'Chain Number')
    hdr.set('TMAX',20.0,'Starting Anneal Temperature')
    hdr.set('ANN_PCT',0.25,'Ideal acceptance Percentage')
    hdr.set('ANN_RNG',0.05,'Range from ideal within which to maintainacceptance')
    hdr.set('CONV_CON',0.5,'Convergence confidence interval') #was 0.05 but having trouble converging so set to 0.5 just for now
    hdr.set('CONV_RMA',1.05,'Convergence Rmax Criterion')
    hdr.set('CONV_STE',20,'Iterations between convergence checks')
    hdr.set('BURN_STE',10,'Iterations between anneal calls in burn-in')
    hdr.set('BURNVRUN',10,'Ratio of normal to burn-in steps')
    hdr.set('PRINT',1,'Whether to Print Debug MSGs (1=verbose,0=silent)')

    hdr['HISTORY']='Last updated on: '+time.strftime("%c") #get current date+time
    hdu1.close

#this is all to do with the GUI where the user is allowed to change the model parameter values etc
if (mdefaults == 'y'):
    import Tkinter
    from functools import partial
    from Tkinter import *
    root=Tk();
    root.title("SurveySim")
    labelframe = LabelFrame(root, text="Luminosity Function Parameters",width=50);
    labelframe.pack(fill="both", expand="yes");
    label=Label(labelframe,text="Initial value/Minimum/Maximum/Fixed",width=50);
    label.pack();

    v_fixed=[DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar()]
    v_min=[DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar()]
    v_max=[DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar()]
    v_init=[DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar(),DoubleVar()]

    def makeform(labelframe, fields):
        entries = []
        ind=0;
        for field in fields:
            row = Frame(labelframe)
            lab = Label(row, width=10, text=field, anchor='w')
            row.pack(side=TOP, padx=2, pady=5)
            lab.pack(side=LEFT)

            ent0 = Entry(row,textvariable=v_fixed[ind])
            ent0.pack(side=RIGHT)
            v_fixed[ind].set(value_fix[ind])
            ent0.pack(side=RIGHT) 
            entries.append((field, ent0))

            ent1 = Entry(row,textvariable=v_max[ind])
            ent1.pack(side=RIGHT)
            v_max[ind].set(value_max[ind])
            entries.append((field, ent1))

            ent2 = Entry(row,textvariable=v_min[ind])
            ent2.pack(side=RIGHT)
            v_min[ind].set(value_min[ind])
            entries.append((field, ent2))

            ent3 = Entry(row,textvariable=v_init[ind])
            ent3.pack(side=RIGHT)
            v_init[ind].set(value_initial[ind])
            entries.append((field, ent3))

            ind=ind+1;
        return entries
        
    labelframe2 = LabelFrame(root, text="Survey properties",width=50);
    labelframe2.pack(fill="both", expand="yes");
    label=Label(labelframe2,text="Area [sqdeg]/wavelength [um]/Flux limit [mJy]",width=50);
    label.pack();
    fields2='Band 1','Band 2','Band 3'

    def fetch(entries):
        for entry in entries:
            field = entry[0]
            text  = entry[1].get()
            print('%s: "%s"' % (field, text)) 

    f1=DoubleVar()
    f2=DoubleVar()
    f3=DoubleVar()
    def makeform2(labelframe2, fields2):
        entries2 = []
        ind=0;
        for field in fields2:
            row = Frame(labelframe2)
            lab = Label(row, width=10, text=field, anchor='w')
            row.pack(side=TOP, padx=2, pady=5)
            lab.pack(side=LEFT)
            if(ind == 0):
                ent0_1 = Entry(row,textvar=f1)
            if(ind == 1):
                ent0_1 = Entry(row,textvar=f2)
            if(ind == 2):
                ent0_1 = Entry(row,textvar=f3)
            ent0_1.pack(side=RIGHT) 
            f1.set(flim[0])
            f2.set(flim[1])
            f3.set(flim[2])

            entries2.append((field, ent0_1))
            ent1_1 = Entry(row)
            ent1_1.insert(10,wavelength[ind])
            ent1_1.pack(side=RIGHT) #, expand=YES, fill=X)
            entries2.append((field, ent1_1))
            ent2_1 = Entry(row)
            ent2_1.insert(10,area[ind])
            ent2_1.pack(side=RIGHT) #, expand=YES, fill=X)
            entries2.append((field, ent2_1))
            ind=ind+1;
        return entries2

    def callback():
        print("Updating parameter values...")
        ind=0
        for field in fields:
            value_initial[ind]=v_init[ind].get()
            value_min[ind]=v_min[ind].get()
            value_max[ind]=v_max[ind].get()
            value_fix[ind]=v_fixed[ind].get()
            ind=ind+1
        flim[0]=f1.get()
        flim[1]=f2.get()
        flim[2]=f3.get()
        update_mfile(modelfile)

    if __name__ == '__main__':
        ents = makeform(labelframe, fields)
        root.bind('<Return>', partial(fetch, ents))
        ents2 = makeform2(labelframe2, fields2)
        root.bind('<Return>', partial(fetch, ents2))
        b1 = Button(labelframe2, text='Update',command=callback)
        b1.pack(side=LEFT,padx=5,pady=5)
        b2 = Button(labelframe2, text='Run',command=runcode)
        b2.pack(side=LEFT,padx=5,pady=5)
        b3 = Button(labelframe2, text='Show results',command=showresults)
        b3.pack(side=LEFT,padx=5,pady=5)
        b4 = Button(labelframe2, text='Quit', command=root.quit)
        b4.pack(side=LEFT, padx=5, pady=5)
        root.mainloop() 

#in case don't want to see GUI, use command line to run the code and display results
if (mdefaults == 'n'):
    torun=raw_input("Run code (y/n)?");
    if (torun == 'y'):
        runcode()
