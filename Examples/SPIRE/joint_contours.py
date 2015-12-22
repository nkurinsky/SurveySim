#!/usr/bin/env python

outputdir=''
models=['A_larger','A-C_larger']

toshow=['L0','PHI0']
xmin=8.45
xmax=10.45
ymin=-3.24
ymax=-0.74

#toshow=['ZBP','ZBQ']
#xmin=0.5
#xmax=2.5
#ymin=1.4
#ymax=2.1
#cmaps=['Paired','Paired']

#toshow=['P','Q']
#xmin=-2.5
#xmax=1.5
#ymin=1.55
#ymax=5.55

#toshow=['P2','Q2']
#xmin=-3.75
#xmax=-2.75
#ymin=0.8
#ymax=1.6

#toshow=['t1','t2']
#xmin=-2.6
#xmax=-0.01
#ymin=1.5
#ymax=7.5
#cmaps=['Paired','Paired']

#toshow=['fa0','fcomp']
#cmaps=['Paired','Paired']
#xmin=0.01
#xmax=1.0
#ymin=0.01
#ymax=0.99

#toshow=['fa0','zbt']
#cmaps=['Paired','Paired']
#xmin=0.01
#xmax=1.0
#ymin=0.5
#ymax=4.5

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

# Perform a kernel density estimate (KDE) on the data
x, y = np.mgrid[xmin:xmax:50j, ymin:ymax:50j]
positions = np.vstack([x.ravel(), y.ravel()])

for im in range(0,len(models)):
    filename=outputdir+models[im]+'_output.fits'
    print 'For model: ',models[im]
    hdus=fits.open(filename)
    phdr=hdus[0].header
    chain1=hdus[4].data
    half_el_col=len(chain1['ACPT0'])/2
    par1=chain1[toshow[0]+'0']
    par1=np.append(par1,chain1[toshow[0]+'1'])
    par1=np.append(par1,chain1[toshow[0]+'2'])
    par1=np.append(par1,chain1[toshow[0]+'3'])
    par1=np.append(par1,chain1[toshow[0]+'4'])

    par2=chain1[toshow[1]+'0']
    par2=np.append(par2,chain1[toshow[1]+'1'])
    par2=np.append(par2,chain1[toshow[1]+'2'])
    par2=np.append(par2,chain1[toshow[1]+'3'])
    par2=np.append(par2,chain1[toshow[1]+'4'])


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
            print '95% confidence: ',frac,test_cut
            conf95_level=test_cut
        if((frac <= 0.68) and (conf68_level == 0)):
            print '68% confidence: ',frac,test_cut
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
        print '95% confidence: ',frac,test_cut
        conf95_level=test_cut
    if((frac <= 0.68) and (conf68_level == 0)):
        print '68% confidence: ',frac,test_cut
        conf68_level=test_cut

# Draw contour lines
cset = g.ax_joint.contour(x,y,f,levels=[conf68_level,conf95_level],colors=('red','red'))
fmt = {}
strs = [ r'68\%', r'95\%' ]
for l,s in zip( cset.levels, strs ):
    fmt[l] = s
g.ax_joint.clabel(cset,cset.levels[::],fmt=fmt,inline=1, fontsize=10)

if(toshow[0] == 'L0'):
    xlabel=r'log(L$^*_0$)'
if(toshow[1] == 'PHI0'):
    ylabel=r'$log(\Phi^*_0)$'
if(toshow[0] == 'P'):
    xlabel=r'$p_1$'
if(toshow[1] == 'Q'):
    ylabel=r'$q_1$'
if(toshow[0] == 'P2'):
    xlabel=r'$p_2$'
if(toshow[1] == 'Q2'):
    ylabel=r'$q_2$'
if(toshow[0] == 'ZBP'):
    xlabel=r'$z_{break,p}$'
if(toshow[1] == 'ZBQ'):
    ylabel=r'$z_{break,q}$'
if(toshow[0] == 't1'):
    xlabel=r'$t_{1}$'
if(toshow[1] == 't2'):
    ylabel=r'$t_{2}$'
if(toshow[0] == 't1'):
    xlabel=r'$t_{1}$'
if(toshow[1] == 't2'):
    ylabel=r'$t_{2}$'
if(toshow[0]=='fa0'):
    xlabel=r'f_{agn,0}'
if(toshow[1]=='fcomp'):
    ylabel=r'f$_{comp}$'
if(toshow[1] == 'zbt'):
    ylabel=r'$z_{break,t}$'    
g.set_axis_labels(xlabel,ylabel)

mean_par1=x[np.where(joint_dist_par1 == np.max(joint_dist_par1)),0]
mean_par2=y[0,np.where(joint_dist_par2 == np.max(joint_dist_par2))]

print mean_par1,mean_par2
g.ax_joint.plot(mean_par1,mean_par2,'*',color='red')
g = g.annotate(stats.pearsonr)
#ax=g.annotate(np.str('Mean')+xlabel)#,xytext=(0.8*xmax,0.8*ymax))#,mean_par1)
plt.show()
