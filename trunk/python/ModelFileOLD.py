# Python class for reading and writing the SurveySim model.fits file
# Created by Noah Kurinsky, 5/9/2015

import os
import sys
import datetime
import time
from filters import *
from astropy.io import fits
        
class filterOptions:

    def __init__(self,fid,limit,err,serr,unit):
        self.fid=fid
        self.limit=limit
        self.err=err
        self.serr=serr
        self.unit=unit

    def __str__(self):
        return self.name()+", Limit="+str(self.limit)+" "+self.unit+", Error="+str(self.err)+", Skew Error="+str(self.serr)+", "+self.unit

    def loadKeys(self,hdr,number):
        name=hdr['Filter'+str(number)]
        self.fid=getFilterID(name)[0]
        self.unit=hdr['units'+str(number)]
        self.limit=hdr['limit'+str(number)]
        self.err=hdr['error'+str(number)]
        self.serr=hdr['serror'+str(number)]

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
        hdr.set('serror'+str(number),self.serr,'flux/magnitude skew error')

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

    def loadKeys(self,hdr,key):
        self.value=hdr[key]
        self.desc=hdr.comments[key]
        keys=[]
        suffs=['FIX','MIN','MAX']
        for i in range(0,3):
            keys.append(key+'_'+suffs[i])
            if(len(keys[i]) > 8):
                keys[i]=keys[i][0:8]
        self.fixed=hdr[keys[0]]
        self.pmin=hdr[keys[1]]
        self.pmax=hdr[keys[2]]
            
    def writeKeys(self,hdr,key):
        hdr.set(key,self.value,self.desc)
        keys=[]
        suffs=['FIX','MIN','MAX']
        for i in range(0,3):
            keys.append(key+'_'+suffs[i])
            if(len(keys[i]) > 8):
                keys[i]=keys[i][0:8]
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
        self.filters=[filterOptions(98,9.3,3.1,0,'mJy'),
                      filterOptions(99,9.6,3.2,0,'mJy'),
                      filterOptions(100,13.5,4.5,0,'mJy')]
        self.survey={
            'area': 4.0,
            'lfForm':0}
        self.axis1='ColorF1F2'
        self.axis2='Flux1'
        self.colsel='None'
        self.convergence={
            'step':20,
            'CI':0.9}
        self.annealing={
            'ideal_pct':0.25,
            'range':0.05,
            'burn_step':10,
            'temp':0.01,
            'learningRate':1.0}
        self.settings={
            'runs':1e4,
            'nchain':5,
            'burnRatio':10,
            'verbosity':2}

#flags whether or not to include Composite and/or Cold Starbursts in model
        self.comp=0
        self.cold=0

        self.params={
            #values based on Gruppioni+2013 plus some additions (fa0)
            'Phi0':parameter(-2.29,-3.0,-1.0,1,"Log Normalization"),
            'L0':parameter(10.12,9.0,11.0,1,"Log Luminosity Knee"),
            'Alpha':parameter(3.0,2.8,3.2,1,"Primary Slope"),
            'Beta':parameter(0.52,0.4,0.6,1,"Secondary Slope"),
            'P':parameter(-0.57,-5.0,0.0,0,"Low Z Density Evolution"),
            'Q':parameter(3.55,0.0,5.0,0,"Low Z Luminosity Evolution"),
            'P2':parameter(-3.92,-5.0,0.0,0,"High Z Density Evolution"),
            'Q2':parameter(1.62,0.0,5.0,0,"High Z Luminosity Evolution"),
            'zbp':parameter(1.1,0.0,5.0,1,"P Z Break"),
            'zbq':parameter(1.85,0.0,5.0,1,"Q Z Break"),
            'fa0':parameter(0.08,0,0.2,1,"AGN fraction at logL=12,z=0"),
            't1':parameter(-0.5,-1,1,1,"T1"),
            't2':parameter(1.3,0.5,2,1,"T2"),
            'zbt':parameter(1,0.5,2,1,"T Z Break"),
            'cexp':parameter(0.0,-3.0,3.0,1,"SED Redshift Evolution"),
            'zbc':parameter(2.0,0.0,5.0,1,"SED Redshift Evolution"),
            'fcomp':parameter(0.5,0.0,1.0,1,"Fraction of AGN Composites"),
            'fcold':parameter(0.5,0.0,1.0,1,"Fraction of Cold SFG")
        }

    def info(self):
        print "Model File Info:"
        print "  Name:",self.filename
        print "  Survey"
        keyPrint(self.survey)
        keyPrint(self.redshift)
        print "  Filters"
        for i in range (0,len(self.filters)) : print "\t",self.filters[i]
        print "  Axes"
        print "\t",self.axis1
        print "\t",self.axis2
        print "\t",self.colsel
        print "  Settings"
        keyPrint(self.settings)
        print "  Convergence"
        keyPrint(self.convergence)
        print "  Annealing"
        keyPrint(self.annealing)
        print "  Parameters"
        keyPrint(self.params)
        print ""

    def fixParams(self,params):
        for param in params:
            self.params[param].fixed=1

    def fitParams(self,params):
        for param in params:
            self.params[param].fixed=0

    def fitAllParams(self):
        self.fitParams(self.params.keys())

    def fixAllParams(self):
        self.fixParams(self.params.keys())

    def filterIDs(self):
        return [self.filters[0].fid,
                self.filters[1].fid,
                self.filters[2].fid]

    def setLF(self,type):
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
        hdus=fits.open(filename)
        phdr=hdus[0].header

#structure needed to be run once if new parameters are to be added to model.fits
#        phdr['comp']=0
#        phdr['cold']=0

        #load LF parameters
        for n,p in self.params.iteritems() : p.loadKeys(phdr,n)

        self.survey['lfForm']=phdr['LF_FORM']
        
        self.survey['area']=phdr['AREA']
        self.colsel=phdr['COLSEL']
        self.axis1=phdr['AXIS1']
        self.axis2=phdr['AXIS2']
        for i in range(0,len(self.filters)): self.filters[i].loadKeys(phdr,i+1)
        
        self.redshift['min']=phdr['ZMIN']
        self.redshift['max']=phdr['ZMAX']
        self.redshift['delta']=phdr['DZ']

        self.settings['runs']=phdr['RUNS']
        self.settings['burnRatio']=phdr['BURNVRUN']
        self.settings['nchain']=phdr['NCHAIN']
        self.settings['verbosity']=phdr['PRINT']

        self.annealing['temp']=phdr['TEMP']
        self.annealing['learningRate']=phdr['LRATE']
        self.annealing['ideal_pct']=phdr['ANN_PCT']
        self.annealing['range']=phdr['ANN_RNG']
        self.annealing['burn_step']=phdr['BURN_STE']

        self.convergence['CI']=phdr['CONV_CON']
        self.convergence['step']=phdr['CONV_STE']
        

    def write(self,filename):
        self.filename=filename
        self.update()

    def update(self):
        #creating filter table
        lam1,lam2,lam3,trans1,trans2,trans3=fill_filters(self.filterIDs());
        col1=fits.Column(name='lambda1',format='FLOAT',array=lam1)
        col2=fits.Column(name='transmission1',format='FLOAT',array=trans1)
        col3=fits.Column(name='lambda2',format='FLOAT',array=lam2)
        col4=fits.Column(name='transmission2',format='FLOAT',array=trans2)
        col5=fits.Column(name='lambda3',format='FLOAT',array=lam3)
        col6=fits.Column(name='transmission3',format='FLOAT',array=trans3)

        cols=fits.ColDefs([col1,col2,col3,col4,col5,col6])
        tbhdu=fits.TableHDU.from_columns(cols)

        tbhdu.header['LSCALE']=(-10,"Wavelength of Filter Lambda")

        hdr=fits.Header() #the header associated with extension=0

        #create/update luminosity function parameters in model file header
        hdr.set('DATE',datetime.date.today().strftime("%B %d, %Y"),'Date of creation')

        #====================================================================
        # LF properties 
        #--------------------------------------------------------------------   
        hdr.set('LF_FORM',self.survey['lfForm'],'The functional form of the LF')

        for n,p in self.params.iteritems() : p.writeKeys(hdr,n)

        #====================================================================
        # Survey properties 
        #--------------------------------------------------------------------   
        hdr.set('AREA',self.survey['area'],'Solid Angle of survey [sq.deg.]')
        hdr.set('COLSEL',self.colsel,'Survey color cut')
        hdr.set('AXIS1',self.axis1,'1st axis to be fit')
        hdr.set('AXIS2',self.axis2,'2nd axis to be fit')

        for i in range(0,len(self.filters)): self.filters[i].writeKeys(hdr,i+1)
        
        #====================================================================
        # Code settings
        #---------------------------------------------------------------------  
        hdr.set('ZMIN',self.redshift['min'],'Simulation minimum redshift')
        hdr.set('ZMAX',self.redshift['max'],'Simulation maximum redshift')
        hdr.set('DZ',self.redshift['delta'],'Redshift Bin Width')
        hdr.set('RUNS',self.settings['runs'],'Number of runs')

        hdr.set('NCHAIN',self.settings['nchain'],'Chain Number')

        hdr.set('TEMP',self.annealing['temp'],'Starting Anneal Temperature')
        hdr.set('LRATE',self.annealing['learningRate'],'Annealing Learning Rate')
        hdr.set('ANN_PCT',self.annealing['ideal_pct'],'Ideal acceptance Percentage')
        hdr.set('ANN_RNG',self.annealing['range'],'Range to maintain acceptance')

        hdr.set('CONV_CON',self.convergence['CI'],'Convergence confidence interval (fraction)') 
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