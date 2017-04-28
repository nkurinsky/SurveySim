#include "SurveySim.h"
#include "functions.h"
#include "mc_util.h"
#include <stdio.h>

using namespace std; 

int main(int argc, char** argv){
	SurveySim sim(argc, argv);
	sim.run();
	sim.save();
}