# Python class for reading and writing the SurveySim model.fits file
# Created by Noah Kurinsky, 5/9/2015

import os
import sys
import datetime
import time
sys.path.append("../Python/")
from filters import *
import pyfits as fits

class filterOptions:

    def __init__(self,fid,limit,err,unit):
        self.fid=fid
        self.limit=limit
        self.err=err
        self.unit=unit

    def __str__(self):
        return self.name()+", Limit="+str(self.limit)+" "+self.unit+", Error="+str(self.err)+" "+self.unit

    def setID(self,name):
        res=getFilterID(name)
        self.fid=res[0]

    def name(self):
        return getFilterName(self.fid)

    def writeKeys(self,hdr,number):
        hdr.set('Filter'+str(number),getFilterName(self.fid),'Filter '+str(number)+' name')
        hdr.set('units'+str(number),self.unit,'[mJy/ABmag]')
        hdr.set('limit'+str(number),self.limit,'flux/magnitude limit')
        hdr.set('error'+str(number),self.err,'flux/magnitude error')

class parameter:

    def __init__(self,value,pmin,pmax,fixed,desc):
        self.value=value
        self.pmin=pmin
        self.pmax=pmax
        self.fixed=fixed
        self.desc=desc

    def __str__(self):
        if(self.fixed == 1):
            return "("+self.desc+") Value="+str(self.value)+", Fixed"
        else:
            return "("+self.desc+") Value="+str(self.value)+", Range=("+str(self.pmin)+","+str(self.pmax)+")"

    def writeKeys(self,hdr,key):
        hdr.set(key,self.value,self.desc)
        keys=[]
        suffs=['FIX','MIN','MAX']
        for i in range(0,3):
            keys.append(key+'_'+suffs[i])
            if(len(keys[i]) > 7):
                keys[i]=keys[i][0:7]
        hdr.set(keys[0],self.fixed,'Fix '+str(key)+' (Y=1/N=0)')
        hdr.set(keys[1],self.pmin,'Minimum '+str(key)+' value')
        hdr.set(keys[2],self.pmax,'Maximum '+str(key)+' value')

def keyPrint(dict):
    for key,val in dict.iteritems():
        print '\t',key,'-\t',val

class ModelFile:    

    def __init__(self):
        self.filename="model.fits"
        self.redshift={ 
            'min':0.01,
            'max':5.00,
            'delta':0.05}
        self.filters=[filterOptions(98,9.3,3.1,'mJy'),
                      filterOptions(99,9.6,3.2,'mJy'),
                      filterOptions(100,13.5,4.5,'mJy')]
        self.survey={
            'area': 4.0,
            'lfForm':0}
        self.axis1='ColorF1F2'
        self.axis2='Flux1'
        self.colsel='None'
        self.convergence={
            'r_max':1.01,
            'step':20,
            'CI':0.01}
        self.annealing={
            'ideal_pct':0.25,
            'range':0.05,
            'burn_step':10,
            'temp':10.0}
        self.settings={
            'runs':1e4,
            'nchain':5,
            'burnRatio':10,
            'verbosity':2}

        self.params={
            #values based on Gruppioni+2013
            'Phi0':parameter(-2.29,-3.0,-1.0,1,"Log Normalization"),
            'L0':parameter(10.12,9.0,11.0,1,"Log Luminosity Knee"),
            'Alpha':parameter(1.15,0.0,1.0,1,"Primary Slope"),
            'Beta':parameter(0.52,2.5,3.5,1,"Secondary Slope"),
            'P1':parameter(-0.57,-8.0,0.0,0,"Low Z Density Evolution"),
            'P2':parameter(-3.92,-8.0,0.0,0,"High Z Density Evolution"),
            'Q1':parameter(3.55,0.0,8.0,0,"Low Z Luminosity Evolution"),
            'Q2':parameter(1.62,0.0,8.0,0,"High Z Luminosity Evolution"),
            'zcut1':parameter(1.1,0.0,5.0,1,"Low Z Cutoff"),
            'zcut2':parameter(1.85,0.0,5.0,1,"High Z Cutoff"),
            'cexp':parameter(0.0,-3.0,3.0,1,"SED Redshift Evolution")
        }

    def info(self):
        print "Model File Info:"
        print "  Name:",self.filename
        print "  Survey"
        keyPrint(self.survey)
        keyPrint(self.redshift)
        print "  Filters"
        for i in range (0,len(self.filters)) : print "\t",self.filters[i]
        print "  Settings"
        keyPrint(self.settings)
        print "  Convergence"
        keyPrint(self.convergence)
        print "  Annealing"
        keyPrint(self.annealing)
        print "  Parameters"
        keyPrint(self.params)
        print ""

    def filterIDs(self):
        return [self.filters[0].fid,
                self.filters[1].fid,
                self.filters[2].fid]

    def setLF(type):
        if(type == 'DoublePowerLaw'):
            self.survey['lfForm']=0
        elif (type == "ModifiedSchecter"):
            self.survey['lfForm']=1
        elif (type == "Schechter"):
            self.survey['lfForm']=2
        else:
            print "Invalid LF Form; options are \"DoublePowerLaw\", \"ModifiedSchechter\", and \"Schechter\""
        
    def load(self,filename):
        self.filename=filename

    def write(self,filename):
        self.filename=filename
        self.update()

    def update(self):

        #creating filter table
        lam1,lam2,lam3,trans1,trans2,trans3=fill_filters(self.filterIDs());
        tbhdu=fits.BinTableHDU.from_columns([
            fits.Column(name='lambda1',format='FLOAT',array=lam1),
            fits.Column(name='transmission1',format='FLOAT',array=trans1),
            fits.Column(name='lambda2',format='FLOAT',array=lam2),
            fits.Column(name='transmission2',format='FLOAT',array=trans2),
            fits.Column(name='lambda3',format='FLOAT',array=lam3),
            fits.Column(name='transmission3',format='FLOAT',array=trans3)])
        tbhdu.header['LSCALE']=(-10,"Wavelength of Filter Lambda")

        hdr=fits.Header() #the header associated with extension=0

        #create/update luminosity function parameters in model file header
        hdr.set('DATE',datetime.date.today().strftime("%B %d, %Y"),'Date of creation')
        hdr.set('LF_FORM',self.survey['lfForm'],'The functional form of the LF')
        for n,p in self.params.iteritems() : p.writeKeys(hdr,n)

        #====================================================================
        # Survey properties 
        #--------------------------------------------------------------------   
        hdr.set('AREA',self.survey['area'],'Solid Angle of survey [sq.deg.]')
        hdr.set('COLSEL',self.colsel,'Survey color cut')
        hdr.set('AXIS1',self.axis1,'1st axis to be fit')
        hdr.set('AXIS2',self.axis2,'2nd axis to be fit')

        for i in range(0,len(self.filters)): self.filters[i].writeKeys(hdr,i)
        
        #====================================================================
        # Code settings
        #---------------------------------------------------------------------  
        hdr.set('ZMIN',self.redshift['min'],'Simulation minimum redshift')
        hdr.set('ZMAX',self.redshift['max'],'Simulation maximum redshift')
        hdr.set('DZ',self.redshift['delta'],'Redshift Bin Width')
        hdr.set('RUNS',self.settings['runs'],'Number of runs')

        hdr.set('NCHAIN',self.settings['nchain'],'Chain Number')

        hdr.set('TMAX',self.annealing['temp'],'Starting Anneal Temperature')
        hdr.set('ANN_PCT',self.annealing['ideal_pct'],'Ideal acceptance Percentage')
        hdr.set('ANN_RNG',self.annealing['range'],'Range to maintain acceptance')

        hdr.set('CONV_CON',self.convergence['CI'],'Convergence confidence interval') 
        hdr.set('CONV_RMA',self.convergence['r_max'],'Convergence Rmax Criterion')
        hdr.set('CONV_STE',self.convergence['step'],'Steps btw convergence checks')

        hdr.set('BURN_STE',self.annealing['burn_step'],'Steps btw anneal calls in burn-in')
        hdr.set('BURNVRUN',self.settings['burnRatio'],'Ratio of normal to burn-in steps')
        hdr.set('PRINT',self.settings['verbosity'],'Print level (0=silent,3=debug)')

        hdr['HISTORY']='Created on: '+time.strftime("%c") #get current date+time

        thdulist=fits.HDUList([fits.PrimaryHDU(header=hdr),tbhdu])
        if(os.path.isfile(self.filename)):
            print 'Replacing existing model file'
            os.remove(self.filename)
            thdulist.writeto(self.filename);
        else:
            print 'Writing a new model file'
            thdulist.writeto(self.filename);
