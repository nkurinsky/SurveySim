#!/usr/bin/env python

import sys
from astropy.io import fits
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
 
import seaborn as sns
sns.set(style='white',font='serif',font_scale=1.5,color_codes='True')

if(len(sys.argv) < 3):
    print "Calling Sequence: "+sys.argv[0]+" outputfile1 outputfile2 [outputfile3 ...]"
    quit()
else:
    ofile1=sys.argv[1]
    ofile2=sys.argv[2]

basename1=ofile1.split('_output.fits')
basename2=ofile2.split('_output.fits')
basename1=basename1[0]
basename2=basename2[0]

ofiles=[ofile1,ofile2]

toshow1=raw_input("Enter the name of a parameter: ")
toshow2=raw_input("Enter the name of another parameter: ")

toshow1=[toshow1]
toshow2=[toshow2]
print 'Displaying joint contours for parameters: ', toshow1, 'and', toshow2

filename=basename1+'_model.fits'
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
npars=len(toshow1)+len(toshow2)

rcut=1.2
chi2_cut=1.0

im=0
for filename in ofiles:
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
    im=im+1

for i_pam in range(npars-1):
    for im in range(0,len(ofiles)):
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
            test_cut=i*0.0003 #need to play with the steps here a bit to ensure good sampling of f
            gpts=(f_tmp > test_cut)
            frac=np.sum(f_tmp[gpts])/np.sum(f_tmp)
            #print gpts
            if((frac <= 0.95) and (conf95_level == 0)):
                print 'Test 95% confidence: ',frac,test_cut
                conf95_level=test_cut
            if((frac <= 0.68) and (conf68_level == 0)):
                print 'Test 68% confidence: ',frac,test_cut
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
            output1_par1=x[np.where(dist_par1 == np.max(dist_par1)),0]
            output1_par2=y[0,np.where(dist_par2 == np.max(dist_par2))]
            cset1 = g.ax_joint.contourf(x,y,f_tmp,levels=[conf68_level,conf68_level],cmap=plt.cm.bone,alpha=0.2)
            cset2 = g.ax_joint.contourf(x,y,f_tmp,levels=[conf95_level,conf95_level],cmap=plt.cm.bone,alpha=0.05)
            tmp_ell1 = cset1.collections[0].get_paths()[0]
            tmp_ell11 = tmp_ell1.vertices
            tmp_par1_68conf=tmp_ell11[0:86,0]
            tmp_par2_68conf=tmp_ell11[0:86,1]
            print '============================================================'
            print 'Outputfile#1: best-fit values w/ 68% confidence levels'
            print '------------------------------------------------------------'
            print toshow1[0],'=',output1_par1[0][0],'+/-',(np.max(tmp_par1_68conf)-output1_par1[0][0]),(output1_par1[0][0]-np.min(tmp_par1_68conf))
            print toshow2[0],'=',output1_par2[0][0],'+/-',(np.max(tmp_par2_68conf)-output1_par2[0][0]),(output1_par2[0][0]-np.min(tmp_par2_68conf))
        if(im==1):
            f_m2=f_m2/np.sum(f_m2)
            dist_par1=f_m2.sum(axis=1)
            dist_par2=f_m2.sum(axis=0)
            #transpose arrays for y-margin plot
            dist_par2_t=dist_par2.T
            y_t=y.T
            g.ax_marg_x.plot(x,dist_par1,color='0.3',linewidth=1)
            g.ax_marg_y.plot(dist_par2_t,y_t,color='0.3',linewidth=1)
            output2_par1=x[np.where(dist_par1 == np.max(dist_par1)),0]
            output2_par2=y[0,np.where(dist_par2 == np.max(dist_par2))]
            cset1 = g.ax_joint.contourf(x,y,f_tmp,levels=[conf68_level,conf68_level],cmap=plt.cm.bone,alpha=0.2)
            cset2 = g.ax_joint.contourf(x,y,f_tmp,levels=[conf95_level,conf95_level],cmap=plt.cm.bone,alpha=0.05)
            tmp_ell1 = cset1.collections[0].get_paths()[0]
            tmp_ell11 = tmp_ell1.vertices
            tmp_par1_68conf=tmp_ell11[0:86,0]
            tmp_par2_68conf=tmp_ell11[0:86,1]
            print '============================================================'
            print 'Outputfile#2: best-fit values w/ 68% confidence levels'
            print '------------------------------------------------------------'
            print toshow1[0],'=',output2_par1[0][0],'+/-',(np.max(tmp_par1_68conf)-output2_par1[0][0]),(output2_par1[0][0]-np.min(tmp_par1_68conf))
            print toshow2[0],'=',output2_par2[0][0],'+/-',(np.max(tmp_par2_68conf)-output2_par2[0][0]),(output2_par2[0][0]-np.min(tmp_par2_68conf))

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
        if((frac <= 0.95) and (conf95_level == 0)):
            print 'Test 95% confidence: ',frac,test_cut
            conf95_level=test_cut
        if((frac <= 0.68) and (conf68_level == 0)):
            print 'Test 68% confidence: ',frac,test_cut
            conf68_level=test_cut
    
    # Draw contour lines
    cset = g.ax_joint.contour(x,y,f,levels=[conf68_level,conf95_level],colors=('red','red'))
    ell1 = cset.collections[0].get_paths()[0]
    ell11 = ell1.vertices
    ell2 = cset.collections[1].get_paths()[0]
    ell22 = ell2.vertices
    par1_68conf=ell11[0:86,0]
    par2_68conf=ell11[0:86,1]

    fmt = {}
    strs = [ r'68\%', r'95\%' ]
    for l,s in zip( cset.levels, strs ):
        fmt[l] = s
    g.ax_joint.clabel(cset,cset.levels[::],fmt=fmt,inline=1, fontsize=10)

    mean_par1=x[np.where(joint_dist_par1 == np.max(joint_dist_par1)),0]
    mean_par2=y[0,np.where(joint_dist_par2 == np.max(joint_dist_par2))]

    print '============================================================'
    print 'JOINT best-fit values w/ 68% confidence levels'
    print '------------------------------------------------------------'
    print toshow1[i_pam],'=',mean_par1[0][0],'+/-',(np.max(par1_68conf)-mean_par1[0][0]),(mean_par1[0][0]-np.min(par1_68conf))
    print toshow2[i_pam],'=',mean_par2[0][0],'+/-',(np.max(par2_68conf)-mean_par2[0][0]),(mean_par2[0][0]-np.min(par2_68conf))
    
    xlabel=toshow1[0]
    ylabel=toshow2[0]
    g.set_axis_labels(xlabel,ylabel)
    g.ax_joint.plot(mean_par1,mean_par2,'*',color='red')
    g = g.annotate(stats.pearsonr,loc=3)
    plt.subplots_adjust(left=0.2, right=0.9, top=0.9, bottom=0.2)
 
plt.show()
