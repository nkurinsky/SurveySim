#!/usr/bin/env python

#Python prelims
import sys, getopt
from astropy.io import fits

def main(argv):
   obsfile = ''
   outputfile = ''
   try:
      opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
   except getopt.GetoptError:
      print 'ObsFile.py -i <obsfile> -o <outputfile>'
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print 'ObsFile.py -i <obsfile> -o <outputfile>'
         sys.exit()
      elif opt in ("-i", "--ifile"):
         obsfile = arg
      elif opt in ("-o", "--ofile"):
         outputfile = arg
   print 'Observations file is ', obsfile
   print 'Output file is ', outputfile

   hdulist = fits.open(obsfile,mode='update')
   hdulist.info()
   prihdr = hdulist[0].header
#check if FHDU keyword exists -- this tells it which table extension to look into
   a='FHDU' in prihdr
   if(a == True):
       print prihdr['FHDU']
   if(a == False):
       prihdr.set('FHDU',1)
   extnum=prihdr['FHDU']
   thdr=hdulist[extnum].header
   a='F3COL' in thdr  # Check for existence
   if(a == True):
       print thdr['F1COL']
   if(a == False):
       f1 = raw_input("Name of column you want to be F1: ")
       ef1 = raw_input("Name of column you want to be EF1: ")
       f2 = raw_input("Name of column you want to be F2: ")
       ef2 = raw_input("Name of column you want to be EF2: ")
       f3 = raw_input("Name of column you want to be F3: ")
       ef3 = raw_input("Name of column you want to be EF3: ")
       thdr.set('F1COL', f1)
       thdr.set('F2COL',f2)
       thdr.set('F3COL', f3)
       thdr.set('EF1COL',ef1)
       thdr.set('EF2COL', ef2)
       thdr.set('EF3COL',ef3)

# uncomment if need to fix something that already has the above keywords set
#   thdr.set('F3COL','F350')
   thdr.set('EF1COL', 'e_F24')
#   thdr.set('EF2COL','et_F250')
#   thdr.set('EF3COL','et_F350')

#   thdr.set('EF1COL','et_F250')
#   thdr.set('EF2COL','et_F350')
#   thdr.set('EF3COL','et_F500')

   hdulist.flush()
   hdulist.close()

if __name__ == "__main__":
   main(sys.argv[1:])
