#!/usr/bin/env python


type_model=['A','B','C','D','E','F']
comp=[4,4,6,6,7,7]
cmaps=['Paired','Paired']

xmin_par=range(7)
xmax_par=range(7)
ymin_par=range(7)
ymax_par=range(7)

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

toshow1=["" for x in range(7)]
toshow2=["" for x in range(7)]

outputdir=''

toshow1[0]='L0'
toshow2[0]='PHI0'
xmin_par[0]=8.45
xmax_par[0]=10.45
ymin_par[0]=-3.24
ymax_par[0]=-0.74

toshow1[1]='ZBP'
toshow2[1]='ZBQ'
xmin_par[1]=0.5
xmax_par[1]=2.5
ymin_par[1]=1.4
ymax_par[1]=2.1

toshow1[2]='P'
toshow2[2]='Q'
xmin_par[2]=-2.5
xmax_par[2]=1.5
ymin_par[2]=1.55
ymax_par[2]=5.55

toshow1[3]='P2'
toshow2[3]='Q2'
xmin_par[3]=-3.75
xmax_par[3]=-2.75
ymin_par[3]=0.8
ymax_par[3]=1.6

toshow1[4]='t1'
toshow2[4]='t2'
xmin_par[4]=-2.6
xmax_par[4]=-0.01
ymin_par[4]=1.5
ymax_par[4]=7.5

toshow1[5]='fa0'
toshow2[5]='zbt'
xmin_par[5]=0.01
xmax_par[5]=0.99
ymin_par[5]=0.5
ymax_par[5]=4.5

toshow1[6]='fcomp'
toshow2[6]='fa0'
xmin_par[6]=0.01
xmax_par[6]=0.99
ymin_par[6]=0.01
ymax_par[6]=0.99


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

print 'spacca'
for i_j in range (len(type_model)):
    print ' '
    print 'Model',type_model[i_j]
    models=[type_model[i_j],type_model[i_j]+'-C']
    #print models
    num_comp=comp[i_j]
    for i_pam in range (num_comp):
        toshow='toshow'+str(i_pam)
        xmin=xmin_par[i_pam]
        xmax=xmax_par[i_pam]
        ymin=ymin_par[i_pam]
        ymax=ymax_par[i_pam]

        # Perform a kernel density estimate (KDE) on the data
        x, y = np.mgrid[xmin:xmax:50j, ymin:ymax:50j]
        positions = np.vstack([x.ravel(), y.ravel()])
        
        for im in range(0,len(models)):
            filename=outputdir+models[im]+'_output.fits'
            #print 'For model: ',models[im]
            hdus=fits.open(filename)
            phdr=hdus[0].header
            chain1=hdus[4].data
            half_el_col=len(chain1['ACPT0'])/2
            par1=chain1[toshow1[i_pam]+'0']
            par1=np.append(par1,chain1[toshow1[i_pam]+'1'])
            par1=np.append(par1,chain1[toshow1[i_pam]+'2'])
            par1=np.append(par1,chain1[toshow1[i_pam]+'3'])
            par1=np.append(par1,chain1[toshow1[i_pam]+'4'])

            par2=chain1[toshow2[i_pam]+'0']
            par2=np.append(par2,chain1[toshow2[i_pam]+'1'])
            par2=np.append(par2,chain1[toshow2[i_pam]+'2'])
            par2=np.append(par2,chain1[toshow2[i_pam]+'3'])
            par2=np.append(par2,chain1[toshow2[i_pam]+'4'])

            ac=chain1['ACPT0']
            ac=np.append(ac,chain1['ACPT1'])
            ac=np.append(ac,chain1['ACPT2'])
            ac=np.append(ac,chain1['ACPT3'])
            ac=np.append(ac,chain1['ACPT4'])

            chisq=chain1['CHISQ0']
            chisq=np.append(chisq,chain1['CHISQ1'])
            chisq=np.append(chisq,chain1['CHISQ2'])
            chisq=np.append(chisq,chain1['CHISQ3'])
            chisq=np.append(chisq,chain1['CHISQ4'])

##            par1=chain1[toshow[0]+'0'][half_el_col:]
##            par1=np.append(par1,chain1[toshow[0]+'1'][half_el_col:])
##            par1=np.append(par1,chain1[toshow[0]+'2'][half_el_col:])
##            par1=np.append(par1,chain1[toshow[0]+'3'][half_el_col:])
##            par1=np.append(par1,chain1[toshow[0]+'4'][half_el_col:])
##
##            par2=chain1[toshow[1]+'0'][half_el_col:]
##            par2=np.append(par2,chain1[toshow[1]+'1'][half_el_col:])
##            par2=np.append(par2,chain1[toshow[1]+'2'][half_el_col:])
##            par2=np.append(par2,chain1[toshow[1]+'3'][half_el_col:])
##            par2=np.append(par2,chain1[toshow[1]+'4'][half_el_col:])
##
##            ac=chain1['ACPT0'][half_el_col:]
##            ac=np.append(ac,chain1['ACPT1'][half_el_col:])
##            ac=np.append(ac,chain1['ACPT2'][half_el_col:])
##            ac=np.append(ac,chain1['ACPT3'][half_el_col:])
##            ac=np.append(ac,chain1['ACPT4'][half_el_col:])
##
##            chisq=chain1['CHISQ0'][half_el_col:]
##            chisq=np.append(chisq,chain1['CHISQ1'][half_el_col:])
##            chisq=np.append(chisq,chain1['CHISQ2'][half_el_col:])
##            chisq=np.append(chisq,chain1['CHISQ3'][half_el_col:])
##            chisq=np.append(chisq,chain1['CHISQ4'][half_el_col:])
##    
##            par1=par1[(np.where(chisq < 100.0)) and (np.where(ac == 1.0))]
##            par2=par2[(np.where(chisq < 100.0)) and (np.where(ac == 1.0))]


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
            #for p1-q1=use 0.01
                test_cut=i*0.0001 #need to play with the steps here a bit to ensure good sampling of f
                gpts=(f_tmp > test_cut)
                frac=np.sum(f_tmp[gpts])/np.sum(f_tmp)
            #print gpts
                if((frac <= 0.95) and (conf95_level == 0)):
                    #print '95% confidence: ',frac,test_cut
                    conf95_level=test_cut
                if((frac <= 0.68) and (conf68_level == 0)):
                    #print '68% confidence: ',frac,test_cut
                    conf68_level=test_cut

            if(im == 0):
                g=sns.JointGrid(par1,par2,space=0,xlim=[xmin,xmax],ylim=[ymin,ymax])
            sns.kdeplot(par1, ax=g.ax_marg_x, legend=False,color='0.3',linewidth=1)
            sns.kdeplot(par2, ax=g.ax_marg_y, legend=False,vertical=True,color='0.3',linewidth=1)
            cset = g.ax_joint.contourf(x,y,f_tmp,levels=[conf68_level,conf95_level],cmap=plt.cm.bone,alpha=0.4)
            cset = g.ax_joint.contourf(x,y,f_tmp,levels=[conf95_level,conf68_level],cmap=plt.cm.bone,alpha=0.1)

        f=(f_m1/np.sum(f_m1))*(f_m2/np.sum(f_m2))
        f=f/np.sum(f) #normalize

        joint_dist_par1=20*f.sum(axis=1)
        joint_dist_par2=20*f.sum(axis=0)

        g.ax_marg_x.plot(x,joint_dist_par1,color='red')

        conf68_level=0
        conf95_level=0
        for i in range(0,50):
            test_cut=i*0.0001 #need to play with the steps here a bit to ensure good sampling of f
            gpts=(f > test_cut)
            frac=np.sum(f[gpts])/np.sum(f)
            #print gpts
            if((frac <= 0.95) and (conf95_level == 0)):
                #print '95% confidence: ',frac,test_cut
                conf95_level=test_cut
            if((frac <= 0.68) and (conf68_level == 0)):
                #print '68% confidence: ',frac,test_cut
                conf68_level=test_cut

        # Draw contour lines
        cset = g.ax_joint.contour(x,y,f,levels=[conf68_level,conf95_level],colors=('red','red'))
        fmt = {}
        strs = [ r'68\%', r'95\%' ]
        for l,s in zip( cset.levels, strs ):
            fmt[l] = s
        g.ax_joint.clabel(cset,cset.levels[::],fmt=fmt,inline=1, fontsize=10)

        mean_par1=x[np.where(joint_dist_par1 == np.max(joint_dist_par1)),0]
        mean_par2=y[0,np.where(joint_dist_par2 == np.max(joint_dist_par2))]

        #print mean_par1,mean_par2

        print toshow1[i_pam],'=',mean_par1[0][0]
        print toshow2[i_pam],'=',mean_par2[0][0]

        if(toshow1[i_pam] == 'L0'):
            xlabel=r'log(L$^*_0$)'
            L0_mean[i_j]=mean_par1[0][0]
        if(toshow2[i_pam] == 'PHI0'):
            ylabel=r'$log(\Phi^*_0)$'
            PHI0_mean[i_j]=mean_par2[0][0]
        if(toshow1[i_pam] == 'P'):
            xlabel=r'$p_1$'
            P_mean[i_j]=mean_par1[0][0]
        if(toshow2[i_pam] == 'Q'):
            ylabel=r'$q_1$'
            Q_mean[i_j]=mean_par2[0][0]
        if(toshow1[i_pam] == 'P2'):
            xlabel=r'$p_2$'
            P2_mean[i_j]=mean_par1[0][0]
        if(toshow2[i_pam] == 'Q2'):
            ylabel=r'$q_2$'
            Q2_mean[i_j]=mean_par2[0][0]
        if(toshow1[i_pam] == 'ZBP'):
            xlabel=r'$z_{break,p}$'
            ZBP_mean[i_j]=mean_par1[0][0]
        if(toshow2[i_pam] == 'ZBQ'):
            ylabel=r'$z_{break,q}$'
            ZBQ_mean[i_j]=mean_par2[0][0]
        if(toshow1[i_pam] == 't1'):
            xlabel=r'$t_{1}$'
            t1_mean[i_j]=mean_par1[0][0]
        if(toshow2[i_pam] == 't2'):
            ylabel=r'$t_{2}$'
            t2_mean[i_j]=mean_par2[0][0]
        if(toshow1[i_pam]=='fa0'):
            xlabel=r'f_{agn,0}'
            fa0_mean[i_j]=mean_par1[0][0]
        if(toshow1[i_pam]=='fcomp'):
            xlabel=r'f$_{comp}$'
            fcomp_mean[i_j]=mean_par1[0][0]
        if(toshow2[i_pam] == 'zbt'):
            ylabel=r'$z_{break,t}$'
            zbt_mean[i_j]=mean_par2[0][0]
        if(toshow2[i_pam]=='fa0'):
            ylabel=r'f_{agn,0}'    
            
        g.set_axis_labels(xlabel,ylabel)
        #print xlabel, ylabel

        g.ax_joint.plot(mean_par1,mean_par2,'*',color='red')
        g = g.annotate(stats.pearsonr)
        #ax=g.annotate(np.str('Mean')+xlabel)#,xytext=(0.8*xmax,0.8*ymax))#,mean_par1)
        plt.subplots_adjust(left=0.2, right=0.9, top=0.9, bottom=0.2)

print ' '
print 'Parameter arrays:'
print ' '
print L0_mean
print PHI0_mean
print ZBP_mean
print ZBQ_mean
print P_mean
print Q_mean
print P2_mean
print Q2_mean
print fa0_mean
print t1_mean
print t2_mean
print zbt_mean
print fcomp_mean




plt.show()
