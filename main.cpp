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





/*Force Interaction::force(Point last, Point p, Point next)
{
	return ext_force(p) + spring_force(last,p,next) + bias_force(p) + int_force(p);
}*/


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
 * Test
 * Interaction, several particles
 * read parameters from file
 * Thermostat
 * Test
 * Bias
 * Test
 * Spline
 * Test
 */
