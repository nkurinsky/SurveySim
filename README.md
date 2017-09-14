Last updated: September 14, 2017 by Noah Kurinsky

SurveySim is an MCMC-based code to fit the luminosity function (LF) evolution based on infrared surveys or simulate such surveys by adopting a particular form of the LF evolution. This code was developed by Noah Kurinsky and Anna Sajina primarily as part of Noah Kurinsky's Tufts undergraduate senior thesis. Details on the code structure and usage are found in docs/UsersGuide.pdf. Please cite Kurinsky et al. 2017 (to appear in ApJ) whenever making use of SurveySim or parts therein. 

In this main directory:
 
* examples      contains pre-built model files as well as their associated output files and figures. these contain the examples discussed in Kurinsky et al. 2017
* templates     contains the SED template we provide with SurveySim
* tags          contains frozen release tar files for official SurveySim releases
* src           contains C++ code for compilation and makefiles
* python        contains the python interface for the SurveySim executable which builds modelfiles and helps parse output

Please refer to the docs/INSTALL file or docs/UsersGuide.pdf for compilation instructions.