#!/usr/bin/env python

import pyfits

flux_limits=[0.001,0.001,0.001]
flux_errors=[0.000001,0.000001,0.000001]

value_initial=[-2.02,10.34,0.50,3.00,-3.75,4.88,3.49,0.0]
value_min=[-3.00,9.00,0.00,0.00,-8.00,0.00,0.00,0.0]
value_max=[-1.00,11.00,2.00,5.00,0.00,8.00,4.00,1.0]
value_fix=[0,0,1,1,0,0,0,0]

zmin=0.01 #Simulation minimum redshift
zmax=5.0 #Simulation maximum redshift
dz=0.05 #Redshift Bin Width
runs=1.e4 #Number of runs
nchain=5 #Chain Number
tmax=10.0 #Starting Anneal Temperature
ann_pct=0.25 #Ideal acceptance Percentage
ann_rng=0.05 #Range to maintain acceptance
conv_con=0.01 #Convergence confidence interval
conv_rma=1.01 #Convergence Rmax Criterion
conv_ste=20 #Steps btw convergence checks
burn_ste=10 #Steps btw anneal calls in burn-in
burnvrun=10 #Ratio of normal to burn-in steps
mesprint=1 #Print Debug MSGs (0=silent,1=critical, 2=info,3=debug)')

hdus=pyfits.open("jwst_model.fits",mode='update')
hdr=hdus[0].header

hdr.set('limit1',flux_limits[0],'flux/magnitude limit')
hdr.set('limit2',flux_limits[1],'flux/magnitude limit')
hdr.set('limit3',flux_limits[2],'flux/magnitude limit')
hdr.set('error1',flux_errors[0],'flux/magnitude error')
hdr.set('error2',flux_errors[1],'flux/magnitude error')
hdr.set('error3',flux_errors[2],'flux/magnitude error')

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
hdr.set('CONV_CON',conv_con,'Convergence confidence interval') 
hdr.set('CONV_RMA',conv_rma,'Convergence Rmax Criterion')
hdr.set('CONV_STE',conv_ste,'Steps btw convergence checks')
hdr.set('BURN_STE',burn_ste,'Steps btw anneal calls in burn-in')
hdr.set('BURNVRUN',burnvrun,'Ratio of normal to burn-in steps')
hdr.set('PRINT',mesprint,'Print level (0=silent,3=debug)')

hdus.flush()
hdus.close()

print(hdr)
