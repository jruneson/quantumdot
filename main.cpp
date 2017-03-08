#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>

//include "GLE-OPT.h"
//include "spline.hpp"
#include "simulation.hpp"
#include "parameters.hpp"
//#include "interaction.hpp"
//#include "polymer.hpp"
//#include "bias.hpp"
//#include "point.hpp"
//#include "observable.hpp"



int main()
{
	Parameters params;
	params.read_file("parameters.cfg");
	Simulation sim(params);
	sim.setup();
	sim.run();
	return 0;
}

/*
 x Point 
 x Polymer
 x Simulation
 x Interaction, basic
 x Observable class
 x Test
 x print certain observables to separate files
 x probability distribution (1D)
 x read parameters from file 
 x Thermostat
 x Timer
 x Test
 x write project description
 x ask Marco about thermostat
 x find error, then put back the new observable class
 * Potential + kinetic energy (virial)
 * fix reading in vector of obs to print every turn
 * Bias
 * Test
 * Spline
 * Test
 * read functions from file
 * Interaction, several particles
 * histogram at 3d
 */
 
 /*
  * Note that mass must be the same of all particles. Otherwise things have to be 
  * updated in GLE and force calculation (curvature is different)
  * 
  * 
  */
