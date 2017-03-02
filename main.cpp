#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>

//include "GLE-OPT.h"
//include "spline.hpp"
#include "simulation.hpp"
#include "interaction.hpp"
#include "polymer.hpp"
#include "bias.hpp"
#include "point.hpp"
#include "observable.hpp"



int main()
{
	Parameters params;
	Simulation sim(params);
	sim.setup();
	sim.run();
	return 0;
}

/*
 x Point 
 x Polymer
 * Simulation
 * Interaction, basic
 * Observable class
 * Test
 * Interaction, several particles
 * read parameters from file (also functions)
 * Thermostat
 * Test
 * Bias
 * Test
 * Spline
 * Test
 */
