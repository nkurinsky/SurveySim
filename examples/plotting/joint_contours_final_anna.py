#!/usr/bin/env python

import numpy
from astropy.io import fits
from scipy.stats import gaussian_kde
import pyfits


type_model=['A','B','C','D','E','F']
comp=[7,7,11,11,12,12]
cmaps=['Paired','Paired']

ymin_par=range(12)
ymax_par=range(12)

L0_mean=[0.0 for x in range(len(type_model))]
PHI0_mean=[0.0 for x in range(len(type_model))]
ZBP_mean=[0.0 for x in range(len(type_model))]
ZBQ_mean=[0.0 for x in range(len(type_model))]
P_mean=[0.0 for x in range(len(type_model))]
Q_mean=[0.0 for x in range(len(type_model))]
P2_mean=[0.0 for x in range(len(type_model))]
Q2_mean=[0.0 for x in range(len(type_model))]
fa0_mean=[0.0 for x in range(len(type_model))]
t1_mean=[0.0 for x in range(len(type_model))]
t2_mean=[0.0 for x in range(len(type_model))]
zbt_mean=[0.0 for x in range(len(type_model))]
fcomp_mean=[0.0 for x in range(len(type_model))]

L0_min=[0.0 for x in range(len(type_model))]
PHI0_min=[0.0 for x in range(len(type_model))]
ZBP_min=[0.0 for x in range(len(type_model))]
ZBQ_min=[0.0 for x in range(len(type_model))]
P_min=[0.0 for x in range(len(type_model))]
Q_min=[0.0 for x in range(len(type_model))]
P2_min=[0.0 for x in range(len(type_model))]
Q2_min=[0.0 for x in range(len(type_model))]
fa0_min=[0.0 for x in range(len(type_model))]
t1_min=[0.0 for x in range(len(type_model))]
t2_min=[0.0 for x in range(len(type_model))]
zbt_min=[0.0 for x in range(len(type_model))]
fcomp_min=[0.0 for x in range(len(type_model))]

L0_max=[0.0 for x in range(len(type_model))]
PHI0_max=[0.0 for x in range(len(type_model))]
ZBP_max=[0.0 for x in range(len(type_model))]
ZBQ_max=[0.0 for x in range(len(type_model))]
P_max=[0.0 for x in range(len(type_model))]
Q_max=[0.0 for x in range(len(type_model))]
P2_max=[0.0 for x in range(len(type_model))]
Q2_max=[0.0 for x in range(len(type_model))]
fa0_max=[0.0 for x in range(len(type_model))]
t1_max=[0.0 for x in range(len(type_model))]
t2_max=[0.0 for x in range(len(type_model))]
zbt_max=[0.0 for x in range(len(type_model))]
fcomp_max=[0.0 for x in range(len(type_model))]

toshow1=['L0' for x in range(12)]
toshow2=["" for x in range(12)]

outputdir='Final_output_2/'

xmin_par=[8.45 for x in range(12)]
xmax_par=[10.45 for x in range(12)]
xmin_par=[9.7 for x in range(12)]
xmax_par=[11.2 for x in range(12)]

toshow2[0]='PHI0'
ymin_par[0]=-3.70#-3.55#-3.70
ymax_par[0]=-1.90

toshow2[1]='ZBP'
ymin_par[1]=0.1
ymax_par[1]=4.2

toshow2[2]='ZBQ'
ymin_par[2]=-0.1
ymax_par[2]=4.2

toshow2[3]='P'
ymin_par[3]=-7.0
ymax_par[3]=7.9

toshow2[4]='Q'
ymin_par[4]=-4.1
ymax_par[4]=7.6

toshow2[5]='P2'
ymin_par[5]=-7.4
ymax_par[5]=7.7

toshow2[6]='Q2'
ymin_par[6]=-7.6
ymax_par[6]=8.4

toshow2[7]='t1'
ymin_par[7]=-2.0#-7.3
ymax_par[7]=8.2

toshow2[8]='t2'
ymin_par[8]=-7.3
ymax_par[8]=8.3

toshow2[9]='fa0'
ymin_par[9]=0.05
ymax_par[9]=0.60

toshow2[10]='zbt'
ymin_par[10]=0.9
ymax_par[10]=4.2

toshow2[11]='fcomp'
ymin_par[11]=-0.2
ymax_par[11]=1.2

import pyfits as fits
from matplotlib import gridspec
import matplotlib.cm as cm
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.figure as fig
import numpy as np
from numpy import std
from scipy import stats
from scipy.stats import norm
from math import ceil
import os
import Tkinter as tk
from matplotlib.colors import Normalize

#Prelims to make prettier colors
import seaborn as sns
sns.set(style='white',font='serif',font_scale=1.5,color_codes='True')
#sns.despine()
#palette=sns.color_palette()

plt.rc('text', usetex=True)

params = {'text.latex.preamble' : [r'\usepackage{siunitx}', r'\usepackage{sfmath}']}
plt.rcParams.update(params)

filename=outputdir+type_model[4]+'_mips'+'_model.fits'

hdus=fits.open(filename)
phdr=hdus[0].header
nchains=phdr['NCHAIN']
conv_step=phdr['CONV_STE']
keywords=phdr.keys()
matching = np.array([s for s in keywords if "_FI" in s])
fixflag=np.array([phdr[i] for i in matching])
gpts=(fixflag < 1)
tmp=matching[gpts]
lpars = ["" for x in range(len(tmp))]
npars=12#len(toshow1)+len(toshow2)

rcut=1.2
chi2_cut=1.0

#for i_j in range (len(type_model)):
for i_j in range (4,6):
    print ' '
    print 'Model',type_model[i_j]
    fits_filename_output='Output_ellipses_2/ellipses_'+type_model[i_j]+'.fits'

    if(os.path.isfile(fits_filename_output)):
        os.remove(fits_filename_output)
    if(type_model[i_j] == 'E'):
        xmin_par=[10.65 for x in range(12)]#10.65
        xmax_par=[11.05 for x in range(12)]#11.05
        ymin_par[0]=-3.45
        ymax_par[0]=-3.05
        #xmin_par=[9.8 for x in range(12)]
        #xmax_par=[11.15 for x in range(12)]
        #ymin_par[0]=-3.55
        #ymax_par[0]=-2.05
    if(type_model[i_j] == 'F'):
        xmin_par=[9.95 for x in range(12)]
        xmax_par=[10.25 for x in range(12)]
        ymin_par[0]=-2.45
        ymax_par[0]=-2.15    
    models=[type_model[i_j]+'_mips',type_model[i_j]+'_spire']
    #print models
    num_comp=comp[i_j]
    for i_pam in range (num_comp):
        #toshow='toshow'+str(i_pam)
        xmin=xmin_par[i_pam]
        xmax=xmax_par[i_pam]
        ymin=ymin_par[i_pam]
        ymax=ymax_par[i_pam]

        # Perform a kernel density estimate (KDE) on the data
        x, y = np.mgrid[xmin:xmax:50j, ymin:ymax:50j]
        positions = np.vstack([x.ravel(), y.ravel()])
            
        for im in range(0,len(models)):
            filename=outputdir+models[im]+'_output.fits'
            hdus=fits.open(filename)
            phdr=hdus[0].header
            converge=hdus[5].data
            half_len=-1
            r=np.ndarray((npars,len(converge)))
            for i_p in range(npars):
                rvalue='R'+str(i_p)
                r[i_p:]=converge[rvalue]
                half_len = -1
            for iconv in range(0,len(converge)):
                if(half_len == -1):
                    tmp=np.array([row[iconv] for row in r])
                    x=tmp[np.where(tmp <= rcut)]
                    if(len(x) >= npars):
                        half_len=iconv                    
            if(im == 0):
                half_len1=half_len*conv_step
                chain1=hdus[4].data
            if(im == 1):
                half_len2=half_len*conv_step
                chain2=hdus[4].data
    
        for i_pam in range(npars-1):
            for im in range(0,len(models)):
                if(im == 0):
                    chain=chain1
                    half_len=half_len1
                if(im == 1):
                    chain=chain2
                    half_len=half_len2
                full_len=len(chain['ACPT0'])
                accepted=np.empty([full_len,nchains])
                chi2=np.empty([full_len,nchains])
                par1=np.empty([full_len,nchains])
                par2=np.empty([full_len,nchains])
                
                for ichain in range(0,nchains):
                    accepted[:,ichain]=chain['ACPT'+str(ichain)]
                    chi2[:,ichain]=chain['CHISQ'+str(ichain)]
                    #print 'median chi2:', np.median(chi2[half_len:full_len,chain])
                    chi2_cut=np.median(chi2[half_len:full_len,ichain])
                    par1[:,ichain]=chain[toshow1[i_pam]+str(ichain)]
                    par2[:,ichain]=chain[toshow2[i_pam]+str(ichain)]
                    if(ichain == 0):
                        par1_total=par1[half_len:full_len,0]
                        par2_total=par2[half_len:full_len,0]
                        #exclude poor chi2 solutions
                        #gpts=(chi2[half_len:full_len,0] < chi2_cut)
                        gpts=(accepted[half_len:full_len,0] == 1)
                        par1_total=par1_total[gpts]
                        par2_total=par2_total[gpts]
                    if(ichain > 0):
                        #gpts=(chi2[half_len:full_len,chain] < chi2_cut)
                        gpts=(accepted[half_len:full_len,ichain] == 1)
                        tmp_par1=par1[half_len:full_len,ichain]
                        tmp_par2=par2[half_len:full_len,ichain]
                        tmp_par1=tmp_par1[gpts]
                        tmp_par2=tmp_par2[gpts]
                        par1_total=np.append(par1_total,tmp_par1)
                        par2_total=np.append(par2_total,tmp_par2)
                par1=par1_total
                par2=par2_total
    
                if(im == 0):
                    #xmin=min(par1)-0.2
                    #xmax=max(par1)+0.2
                    #ymin=min(par2)-0.2
                    #ymax=max(par2)+0.2
                    x, y = np.mgrid[xmin:xmax:50j, ymin:ymax:50j]
                    positions = np.vstack([x.ravel(), y.ravel()])
    
                values = np.vstack([par1, par2])
                kernel = stats.gaussian_kde(values)
                if(im ==0):
                    f_m1 = np.reshape(kernel(positions).T, x.shape)
                    f_tmp=f_m1
                if(im == 1):
                    f_m2 = np.reshape(kernel(positions).T,x.shape)
                    f_tmp=f_m2
                conf68_level=0
                conf95_level=0
                f_tmp=f_tmp/np.sum(f_tmp)
                for i in range(0,50):
                    test_cut=i*0.0003 #need to play with the steps here a bit to ensure good sampling of f
                    gpts=(f_tmp > test_cut)
                    frac=np.sum(f_tmp[gpts])/np.sum(f_tmp)
                #print gpts
                    if((frac <= 0.95) and (conf95_level == 0)):
                        print '95% confidence: ',frac,test_cut
                        conf95_level=test_cut
                    if((frac <= 0.68) and (conf68_level == 0)):
                        print '68% confidence: ',frac,test_cut
                        conf68_level=test_cut
    
                if(im == 0):
                    g=sns.JointGrid(par1,par2,space=0,xlim=[xmin,xmax],ylim=[ymin,ymax])
                    f_m1=f_m1/np.sum(f_m1)
                    dist_par1=f_m1.sum(axis=1)
                    dist_par2=f_m1.sum(axis=0)
                    #transpose arrays for y-margin plot
                    dist_par2_t=dist_par2.T
                    y_t=y.T
                    g.ax_marg_x.plot(x,dist_par1,color='0.3',linewidth=1)
                    g.ax_marg_y.plot(dist_par2_t,y_t,color='0.3',linewidth=1)
                if(im==1):
                    f_m2=f_m2/np.sum(f_m2)
                    dist_par1=f_m2.sum(axis=1)
                    dist_par2=f_m2.sum(axis=0)
                    #transpose arrays for y-margin plot
                    dist_par2_t=dist_par2.T
                    y_t=y.T
                    g.ax_marg_x.plot(x,dist_par1,color='0.3',linewidth=1)
                    g.ax_marg_y.plot(dist_par2_t,y_t,color='0.3',linewidth=1)
               # sns.kdeplot(par1, ax=g.ax_marg_x, legend=False,color='0.3',linewidth=1)
                #sns.kdeplot(par2, ax=g.ax_marg_y, legend=False,vertical=True,color='0.3',linewidth=1)
                cset = g.ax_joint.contourf(x,y,f_tmp,levels=[conf68_level,conf95_level],cmap=plt.cm.bone,alpha=0.4)
                cset = g.ax_joint.contourf(x,y,f_tmp,levels=[conf95_level,conf68_level],cmap=plt.cm.bone,alpha=0.1)
    
            f=(f_m1/np.sum(f_m1))*(f_m2/np.sum(f_m2))
            
            f=f/np.sum(f) #the Joint normalized probability distribution
    
            joint_dist_par1=f.sum(axis=1)
            joint_dist_par2=f.sum(axis=0)
    
            #transpose arrays for y-margin plot
            joint_dist_par2_t=joint_dist_par2.T
            y_t=y.T
    
            g.ax_marg_x.plot(x,joint_dist_par1,color='red',linewidth=3)
            g.ax_marg_y.plot(joint_dist_par2_t,y_t,color='red',linewidth=3)
    
            conf68_level=0
            conf95_level=0
            for i in range(0,50):
                test_cut=i*0.0001 #need to play with the steps here a bit to ensure good sampling of f
                gpts=(f > test_cut)
                frac=np.sum(f[gpts])/np.sum(f)
                #print gpts
                if((frac <= 0.95) and (conf95_level == 0)):
                    print '95% confidence: ',frac,test_cut
                    conf95_level=test_cut
                if((frac <= 0.68) and (conf68_level == 0)):
                    print '68% confidence: ',frac,test_cut
                    conf68_level=test_cut
        
            # Draw contour lines
            cset = g.ax_joint.contour(x,y,f,levels=[conf68_level,conf95_level],colors=('red','red'))
    
            ell1 = cset.collections[0].get_paths()[0]
            ell11 = ell1.vertices
            ell2 = cset.collections[1].get_paths()[0]
            ell22 = ell2.vertices
            par1_68conf=ell11[0:86,0]
            par2_68conf=ell11[0:86,1]
            
            pyfits.append(fits_filename_output, ell11)
                
            
            fmt = {}
            strs = [ r'68\%', r'95\%' ]
            for l,s in zip( cset.levels, strs ):
                fmt[l] = s
            g.ax_joint.clabel(cset,cset.levels[::],fmt=fmt,inline=1, fontsize=10)
    

        #help(g.ax_joint.contour)

            mean_par1=x[np.where(joint_dist_par1 == np.max(joint_dist_par1)),0]
            mean_par2=y[0,np.where(joint_dist_par2 == np.max(joint_dist_par2))]

            print toshow1[i_pam],'=',mean_par1[0][0],'+/-',(np.max(par1_68conf)-mean_par1[0][0]),(mean_par1[0][0]-np.min(par1_68conf))
            print toshow2[i_pam],'=',mean_par2[0][0],'+/-',(np.max(par2_68conf)-mean_par2[0][0]),(mean_par2[0][0]-np.min(par2_68conf))

            xlabel=r'log(L$^*_0$)'
        
            if(toshow2[i_pam] == 'PHI0'):
                ylabel=r'$log(\Phi^*_0)$'
                PHI0_mean[i_j]=mean_par2[0][0]
                PHI0_min[i_j]=(mean_par2[0][0]-np.min(par2_68conf))
                PHI0_max[i_j]=(np.max(par2_68conf)-mean_par2[0][0])
                L0_mean[i_j]=mean_par1[0][0]
                L0_min[i_j]=(mean_par1[0][0]-np.min(par1_68conf))
                L0_max[i_j]=(np.max(par1_68conf)-mean_par1[0][0])
            if(toshow2[i_pam] == 'P'):
                ylabel=r'$p_1$'
                P_mean[i_j]=mean_par2[0][0]
                P_min[i_j]=(mean_par2[0][0]-np.min(par2_68conf))
                P_max[i_j]=(np.max(par2_68conf)-mean_par2[0][0])
            if(toshow2[i_pam] == 'Q'):
                ylabel=r'$q_1$'
                Q_mean[i_j]=mean_par2[0][0]
                Q_min[i_j]=(mean_par2[0][0]-np.min(par2_68conf))
                Q_max[i_j]=(np.max(par2_68conf)-mean_par2[0][0])
            if(toshow2[i_pam] == 'P2'):
                ylabel=r'$p_2$'
                P2_mean[i_j]=mean_par2[0][0]
                P2_min[i_j]=(mean_par2[0][0]-np.min(par2_68conf))
                P2_max[i_j]=(np.max(par2_68conf)-mean_par2[0][0])
            if(toshow2[i_pam] == 'Q2'):
                ylabel=r'$q_2$'
                Q2_mean[i_j]=mean_par2[0][0]
                Q2_min[i_j]=(mean_par2[0][0]-np.min(par2_68conf))
                Q2_max[i_j]=(np.max(par2_68conf)-mean_par2[0][0])
            if(toshow2[i_pam] == 'ZBP'):
                ylabel=r'$z_{break,p}$'
                ZBP_mean[i_j]=mean_par2[0][0]
                ZBP_min[i_j]=(mean_par2[0][0]-np.min(par2_68conf))
                ZBP_max[i_j]=(np.max(par2_68conf)-mean_par2[0][0])
            if(toshow2[i_pam] == 'ZBQ'):
                ylabel=r'$z_{break,q}$'
                ZBQ_mean[i_j]=mean_par2[0][0]
                ZBQ_min[i_j]=(mean_par2[0][0]-np.min(par2_68conf))
                ZBQ_max[i_j]=(np.max(par2_68conf)-mean_par2[0][0])
            if(toshow2[i_pam] == 't1'):
                ylabel=r'$t_{1}$'
                t1_mean[i_j]=mean_par2[0][0]
                t1_min[i_j]=(mean_par2[0][0]-np.min(par2_68conf))
                t1_max[i_j]=(np.max(par2_68conf)-mean_par2[0][0])
            if(toshow2[i_pam] == 't2'):
                ylabel=r'$t_{2}$'
                t2_mean[i_j]=mean_par2[0][0]
                t2_min[i_j]=(mean_par2[0][0]-np.min(par2_68conf))
                t2_max[i_j]=(np.max(par2_68conf)-mean_par2[0][0])
            if(toshow2[i_pam]=='fa0'):
                ylabel=r'f_{agn,0}'
                fa0_mean[i_j]=mean_par2[0][0]
                fa0_min[i_j]=(mean_par2[0][0]-np.min(par2_68conf))
                fa0_max[i_j]=(np.max(par2_68conf)-mean_par2[0][0])
            if(toshow2[i_pam]=='fcomp'):
                ylabel=r'f$_{comp}$'
                fcomp_mean[i_j]=mean_par2[0][0]
                fcomp_min[i_j]=(mean_par2[0][0]-np.min(par2_68conf))
                fcomp_max[i_j]=(np.max(par2_68conf)-mean_par2[0][0])
            if(toshow2[i_pam] == 'zbt'):
                ylabel=r'$z_{break,t}$'
                zbt_mean[i_j]=mean_par2[0][0]
                zbt_min[i_j]=(mean_par2[0][0]-np.min(par2_68conf))
                zbt_max[i_j]=(np.max(par2_68conf)-mean_par2[0][0])
                

            g.set_axis_labels(xlabel,ylabel)
            g.ax_joint.plot(mean_par1,mean_par2,'*',color='red')
            g = g.annotate(stats.pearsonr)
            plt.subplots_adjust(left=0.2, right=0.9, top=0.9, bottom=0.2)
    
#plt.show()

print ' '
print 'parameter_array=[',L0_mean,',$'
print PHI0_mean,',$'
print ZBP_mean,',$'
print ZBQ_mean,',$'
print P_mean,',$'
print Q_mean,',$'
print P2_mean,',$'
print Q2_mean,',$'
print fa0_mean,',$'
print t1_mean,',$'
print t2_mean,',$'
print zbt_mean,',$'
print fcomp_mean,']'

matrix=[L0_mean,PHI0_mean,ZBP_mean,ZBQ_mean,P_mean,Q_mean,P2_mean,Q2_mean,fa0_mean,t1_mean,t2_mean,zbt_mean,fcomp_mean]
np.savetxt('Output_ellipses_2/best_fit_array.txt', matrix)


print '      '
print "A & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ \\\\ & / & / & / & / & / & & & \\\\"% (L0_mean[0],L0_max[0],L0_min[0],PHI0_mean[0],PHI0_max[0],PHI0_min[0],ZBP_mean[0],ZBP_max[0],ZBP_min[0],ZBQ_mean[0],ZBQ_max[0],ZBQ_min[0],P_mean[0],P_max[0],P_min[0],Q_mean[0],Q_max[0],Q_min[0],P2_mean[0],P2_max[0],P2_min[0],Q2_mean[0],Q2_max[0],Q2_min[0]) 
print "B & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ \\\\ & / & / & / & / & / & & & \\\\"% (L0_mean[1],L0_max[1],L0_min[1],PHI0_mean[1],PHI0_max[1],PHI0_min[1],ZBP_mean[1],ZBP_max[1],ZBP_min[1],ZBQ_mean[1],ZBQ_max[1],ZBQ_min[1],P_mean[1],P_max[1],P_min[1],Q_mean[1],Q_max[1],Q_min[1],P2_mean[1],P2_max[1],P2_min[1],Q2_mean[1],Q2_max[1],Q2_min[1]) 
print "C & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ \\\\ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & / & & & \\\\"% (L0_mean[2],L0_max[2],L0_min[2],PHI0_mean[2],PHI0_max[2],PHI0_min[2],ZBP_mean[2],ZBP_max[2],ZBP_min[2],ZBQ_mean[2],ZBQ_max[2],ZBQ_min[2],P_mean[2],P_max[2],P_min[2],Q_mean[2],Q_max[2],Q_min[2],P2_mean[2],P2_max[2],P2_min[2],Q2_mean[2],Q2_max[2],Q2_min[2],fa0_mean[2],fa0_max[2],fa0_min[2],t1_mean[2],t1_max[2],t1_min[2],t2_mean[2],t2_max[2],t2_min[2],zbt_mean[2],zbt_max[2],zbt_min[2]) 
print "D & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ \\\\ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & / & & & \\\\"% (L0_mean[3],L0_max[3],L0_min[3],PHI0_mean[3],PHI0_max[3],PHI0_min[3],ZBP_mean[3],ZBP_max[3],ZBP_min[3],ZBQ_mean[3],ZBQ_max[3],ZBQ_min[3],P_mean[3],P_max[3],P_min[3],Q_mean[3],Q_max[3],Q_min[3],P2_mean[3],P2_max[3],P2_min[3],Q2_mean[3],Q2_max[3],Q2_min[3],fa0_mean[3],fa0_max[3],fa0_min[3],t1_mean[3],t1_max[3],t1_min[3],t2_mean[3],t2_max[3],t2_min[3],zbt_mean[3],zbt_max[3],zbt_min[3]) 
print "E & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ \\\\ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & & & \\\\"% (L0_mean[4],L0_max[4],L0_min[4],PHI0_mean[4],PHI0_max[4],PHI0_min[4],ZBP_mean[4],ZBP_max[4],ZBP_min[4],ZBQ_mean[4],ZBQ_max[4],ZBQ_min[4],P_mean[4],P_max[4],P_min[4],Q_mean[4],Q_max[4],Q_min[4],P2_mean[4],P2_max[4],P2_min[4],Q2_mean[4],Q2_max[4],Q2_min[4],fa0_mean[4],fa0_max[4],fa0_min[4],t1_mean[4],t1_max[4],t1_min[4],t2_mean[4],t2_max[4],t2_min[4],zbt_mean[4],zbt_max[4],zbt_min[4],fcomp_mean[4],fcomp_max[4],fcomp_min[4]) 
print "F & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ \\\\ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & & & \\\\"% (L0_mean[5],L0_max[5],L0_min[5],PHI0_mean[5],PHI0_max[5],PHI0_min[5],ZBP_mean[5],ZBP_max[5],ZBP_min[5],ZBQ_mean[5],ZBQ_max[5],ZBQ_min[5],P_mean[5],P_max[5],P_min[5],Q_mean[5],Q_max[5],Q_min[5],P2_mean[5],P2_max[5],P2_min[5],Q2_mean[5],Q2_max[5],Q2_min[5],fa0_mean[5],fa0_max[5],fa0_min[5],t1_mean[5],t1_max[5],t1_min[5],t2_mean[5],t2_max[5],t2_min[5],zbt_mean[5],zbt_max[5],zbt_min[5],fcomp_mean[5],fcomp_max[5],fcomp_min[5]) 

#quit ()
plt.show()
#quit()

