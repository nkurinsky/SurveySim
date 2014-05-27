//Anna Sajina
//9/27/2012
//This header file holds the values of constants

#ifndef CONSTANTS_H
#define CONSTANTS_H

//this includes some basic math definitions including
//M_PI for pi
#define _USE_MATH_DEFINES

//general astronomy/physics constants
#define MPC_TO_METER (3.08568025e22) // m/MPc
#define Wm2Hz_TO_mJy (1.e29) // W/m2/Hz to mJy
#define LSUN (3.839e26) // solar luminosity
#define speed_of_light (2.998e8) // c [m/s]

//cosmology
#define Ho (73.8) //km/s/Mpc, based on Riess et al. 2011, ApJ
#define OM (0.272) // from Komatsu et al. 2011
#define OL (0.728) // from Komatsu et al. 2011
#define d_hub (4062.33) // Hubble distance, i.e. c/Ho [Mpc]

#endif
