#!/usr/bin/env python

import numpy
import pyfits
from scipy.stats import gaussian_kde


type_model=['A','B','C','D','E','F']
type_model=['E','F']
comp=[7,7,11,11,12,12]
comp=[10,10]
cmaps=['Paired','Paired']

ymin_par=range(10)
ymax_par=range(10)

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

toshow1=['L0' for x in range(10)]
toshow2=["" for x in range(10)]

outputdir='Final_output_2/'

xmin_par=[8.45 for x in range(10)]
xmax_par=[10.45 for x in range(10)]
xmin_par=[9.7 for x in range(10)]
xmax_par=[11.2 for x in range(10)]

toshow2[0]='PHI0'
ymin_par[0]=-3.70#-3.55#-3.70
ymax_par[0]=-1.90

#toshow2[1]='ZBP'
#ymin_par[1]=0.1
#ymax_par[1]=4.2

#toshow2[2]='ZBQ'
#ymin_par[2]=-0.1
#ymax_par[2]=4.2

toshow2[1]='P'
ymin_par[1]=-6.5
ymax_par[1]=1.7
ymin_par[1]=-7.0
ymax_par[1]=7.9

toshow2[2]='Q'
ymin_par[2]=-0.1
ymax_par[2]=7.5
ymin_par[2]=-4.1
ymax_par[2]=7.6

toshow2[3]='P2'
ymin_par[3]=-7.2
ymax_par[3]=1.8
ymin_par[3]=-7.4
ymax_par[3]=7.7

toshow2[4]='Q2'
ymin_par[4]=-0.5
ymax_par[4]=7.2
ymin_par[4]=-7.6
ymax_par[4]=8.4

toshow2[5]='t1'
ymin_par[5]=-3.1
ymax_par[5]=1.6
ymin_par[5]=-7.3
ymax_par[5]=8.2

toshow2[6]='t2'
ymin_par[6]=-3.0
ymax_par[6]=9.6
ymin_par[6]=-7.3
ymax_par[6]=8.3

toshow2[7]='fa0'
ymin_par[7]=-0.3
ymax_par[7]=6.0
ymin_par[7]=0.05
ymax_par[7]=0.60

toshow2[8]='zbt'
ymin_par[8]=0.1
ymax_par[8]=4.2
ymin_par[8]=0.9
ymax_par[8]=4.2

toshow2[9]='fcomp'
ymin_par[9]=-0.2
ymax_par[9]=1.2

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

#from OutputFile import *

#Prelims to make prettier colors
import seaborn as sns
sns.set(style='white',font='serif',font_scale=1.5,color_codes='True')
#sns.despine()
#palette=sns.color_palette()

plt.rc('text', usetex=True)

params = {'text.latex.preamble' : [r'\usepackage{siunitx}', r'\usepackage{sfmath}']}
plt.rcParams.update(params)

rcut=1.1
nchains=10
chi2_cut=1.0

for i_j in range (len(type_model)):
    print ' '
    print 'Model',type_model[i_j]
    if(type_model[i_j] == 'E'):
        xmin_par=[10.65 for x in range(10)]
        xmax_par=[11.05 for x in range(10)]
        ymin_par[0]=-3.45
        ymax_par[0]=-3.05
    if(type_model[i_j] == 'F'):
        xmin_par=[9.90 for x in range(10)]
        xmax_par=[10.30 for x in range(10)]
        ymin_par[0]=-2.5
        ymax_par[0]=-2.10    
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
            r0=converge['R0']
            r1=converge['R1']
            r2=converge['R2']
            r3=converge['R3']
            r4=converge['R4']
            r5=converge['R5']
            r6=converge['R6']
            r7=converge['R7']
            r8=converge['R8']
            r9=converge['R9']
            r10=converge['R10']
            #r11=converge['R11']
            #r12=converge['R12']

            gpts_r=(r0 > 0)
            conv_len=len(r0[gpts_r])
            half_len = -1
            for iconv in range(0,conv_len):
                if(half_len == -1):
                    if((r0[iconv] <= rcut) & (r1[iconv] <= rcut) & (r2[iconv] <= rcut) & (r3[iconv] <= rcut) & (r4[iconv] <= rcut) & (r5[iconv] <= rcut) & (r6[iconv] <= rcut) & (r7[iconv] <= rcut) & (r8[iconv] <= rcut) & (r9[iconv] <= rcut) & (r10[iconv] <= rcut)):# & (r11[iconv] <= rcut) & (r12[iconv] <= rcut)):
                       half_len=iconv

            half_len=half_len*20
            chain1=hdus[4].data
            full_len=len(chain1['ACPT0'])
            accepted=np.empty([full_len,nchains])
            chi2=np.empty([full_len,nchains])
            par1=np.empty([full_len,nchains])
            par2=np.empty([full_len,nchains])
            for chain in range(0,nchains):
                accepted[:,chain]=chain1['ACPT'+str(chain)]
                chi2[:,chain]=chain1['CHISQ'+str(chain)]
                #print 'median chi2:', np.median(chi2[half_len:full_len,chain])
                chi2_cut=np.median(chi2[half_len:full_len,chain])
                par1[:,chain]=chain1[toshow1[i_pam]+str(chain)]
                par2[:,chain]=chain1[toshow2[i_pam]+str(chain)]
                if(chain == 0):
                    par1_total=par1[half_len:full_len,0]
                    par2_total=par2[half_len:full_len,0]
                    #exclude poor chi2 solutions
                    #gpts=(chi2[half_len:full_len,0] < chi2_cut)
                    gpts=(accepted[half_len:full_len,0] == 1)
                    par1_total=par1_total[gpts]
                    par2_total=par2_total[gpts]
                if(chain > 0):
                    gpts=(chi2[half_len:full_len,chain] < chi2_cut)
                    gpts=(accepted[half_len:full_len,chain] == 1)
                    tmp_par1=par1[half_len:full_len,chain]
                    tmp_par2=par2[half_len:full_len,chain]
                    tmp_par1=tmp_par1[gpts]
                    tmp_par2=tmp_par2[gpts]
                    par1_total=np.append(par1_total,tmp_par1)
                    par2_total=np.append(par2_total,tmp_par2)

            par1=par1_total
            par2=par2_total

            if(im == 0):
                xmin=min(par1)-0.2
                xmax=max(par1)+0.2
                ymin=min(par2)-0.2
                ymax=max(par2)+0.2
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
                test_cut=i*0.00003 #need to play with the steps here a bit to ensure good sampling of f
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

        fits_filename_output='Output_ellipses_2/ellipses_'+type_model[i_j]+'.fits'
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


print '      '
print "A & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ \\\\ & / & / & / & / & / & & & \\\\"% (L0_mean[0],L0_max[0],L0_min[0],PHI0_mean[0],PHI0_max[0],PHI0_min[0],ZBP_mean[0],ZBP_max[0],ZBP_min[0],ZBQ_mean[0],ZBQ_max[0],ZBQ_min[0],P_mean[0],P_max[0],P_min[0],Q_mean[0],Q_max[0],Q_min[0],P2_mean[0],P2_max[0],P2_min[0],Q2_mean[0],Q2_max[0],Q2_min[0]) 
print "B & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ \\\\ & / & / & / & / & / & & & \\\\"% (L0_mean[1],L0_max[1],L0_min[1],PHI0_mean[1],PHI0_max[1],PHI0_min[1],ZBP_mean[1],ZBP_max[1],ZBP_min[1],ZBQ_mean[1],ZBQ_max[1],ZBQ_min[1],P_mean[1],P_max[1],P_min[1],Q_mean[1],Q_max[1],Q_min[1],P2_mean[1],P2_max[1],P2_min[1],Q2_mean[1],Q2_max[1],Q2_min[1]) 
print "C & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ \\\\ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & / & & & \\\\"% (L0_mean[2],L0_max[2],L0_min[2],PHI0_mean[2],PHI0_max[2],PHI0_min[2],ZBP_mean[2],ZBP_max[2],ZBP_min[2],ZBQ_mean[2],ZBQ_max[2],ZBQ_min[2],P_mean[2],P_max[2],P_min[2],Q_mean[2],Q_max[2],Q_min[2],P2_mean[2],P2_max[2],P2_min[2],Q2_mean[2],Q2_max[2],Q2_min[2],fa0_mean[2],fa0_max[2],fa0_min[2],t1_mean[2],t1_max[2],t1_min[2],t2_mean[2],t2_max[2],t2_min[2],zbt_mean[2],zbt_max[2],zbt_min[2]) 
print "D & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ \\\\ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & / & & & \\\\"% (L0_mean[3],L0_max[3],L0_min[3],PHI0_mean[3],PHI0_max[3],PHI0_min[3],ZBP_mean[3],ZBP_max[3],ZBP_min[3],ZBQ_mean[3],ZBQ_max[3],ZBQ_min[3],P_mean[3],P_max[3],P_min[3],Q_mean[3],Q_max[3],Q_min[3],P2_mean[3],P2_max[3],P2_min[3],Q2_mean[3],Q2_max[3],Q2_min[3],fa0_mean[3],fa0_max[3],fa0_min[3],t1_mean[3],t1_max[3],t1_min[3],t2_mean[3],t2_max[3],t2_min[3],zbt_mean[3],zbt_max[3],zbt_min[3]) 
print "E & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ \\\\ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & & & \\\\"% (L0_mean[4],L0_max[4],L0_min[4],PHI0_mean[4],PHI0_max[4],PHI0_min[4],ZBP_mean[4],ZBP_max[4],ZBP_min[4],ZBQ_mean[4],ZBQ_max[4],ZBQ_min[4],P_mean[4],P_max[4],P_min[4],Q_mean[4],Q_max[4],Q_min[4],P2_mean[4],P2_max[4],P2_min[4],Q2_mean[4],Q2_max[4],Q2_min[4],fa0_mean[4],fa0_max[4],fa0_min[4],t1_mean[4],t1_max[4],t1_min[4],t2_mean[4],t2_max[4],t2_min[4],zbt_mean[4],zbt_max[4],zbt_min[4],fcomp_mean[4],fcomp_max[4],fcomp_min[4]) 
print "F & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ \\\\ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & $%6.2f^{+%6.2f}_{-%6.2f}$ & & & \\\\"% (L0_mean[5],L0_max[5],L0_min[5],PHI0_mean[5],PHI0_max[5],PHI0_min[5],ZBP_mean[5],ZBP_max[5],ZBP_min[5],ZBQ_mean[5],ZBQ_max[5],ZBQ_min[5],P_mean[5],P_max[5],P_min[5],Q_mean[5],Q_max[5],Q_min[5],P2_mean[5],P2_max[5],P2_min[5],Q2_mean[5],Q2_max[5],Q2_min[5],fa0_mean[5],fa0_max[5],fa0_min[5],t1_mean[5],t1_max[5],t1_min[5],t2_mean[5],t2_max[5],t2_min[5],zbt_mean[5],zbt_max[5],zbt_min[5],fcomp_mean[5],fcomp_max[5],fcomp_min[5]) 


plt.show()
#quit()

