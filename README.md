Last updated: June 23, 2015 by Anna Sajina

SurveySim is an MCMC-based code to fit the luminosity function (LF) evolution based on infrared surveys or simulate such surveys by adopting a particular form of the LF evolution. This code was developed by Noah Kurinsky and Anna Sajina primarily as part of Noah Kurinsky's Tufts undergraduate senior thesis. Details on the code structure and usage are found in trunk/manual/UsersGuide.pdf. Please cite Kurinsky et al. 2015 (in prep) whenever making use of SurveySim or parts therein. 

In this main directory:
 
 SurveySim.py: opens a user-friendly GUI interface to SurveySim
 run.sh:       is the quickest, command line, means of running SurveySim based
               on pre-built model files (see trunk/model for examples)

 Examples      contains pre-built model files as well as their associated output files and figures. these contain the examples discussed in Kurinsky et al. 2015
 OUTPUT        the default directory where output fits files as well as any other figures are generated
 trunk/python  contains various useful python routines such as for generating the observations fits files, the model files, examining the model and output files, viewing filter profiles, summarizing output results etc.

The actual code is contained within the "trunk" directory. Please refer to the README file therein for compilation instructions.
