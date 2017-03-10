#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>

#include "simulation.hpp"
#include "parameters.hpp"



int main()
{
	std::vector<double> times = {100000};
	std::vector<double> Ps = {10,30,50,70};
	Parameters params;
	params.read_file("parameters.cfg");
	std::ofstream results_file("results.dat");
	results_file.precision(8);
	results_file << "%P\tTTot";
	for(int id : params.to_measure)
		results_file << "\tObsId " << id << "\t\tError\t";
	results_file << std::endl;
	for(double ttot : times)
	{
		params.total_time = ttot;
		for(int P : Ps)
		{
			params.num_beads = P;
			params.calculate_dependencies();
			Simulation sim(params, results_file);
			sim.setup();
			sim.run();
		}
	}
	results_file.close();
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
 x Potential + kinetic energy (virial)
 x fix reading in vector of obs to print every turn
 x turn missing Amatrix into an error
 * centroid virial
 * fermion and boson
 * Bias
 * Test
 * Spline
 * Test
 * read functions from file - or choose function with an integer input
 * Interaction, several particles
 * histogram at 3d
 */
 
 /*
  * Note that mass must be the same of all particles. Otherwise things have to be 
  * updated in GLE and force calculation (curvature is different)
  * 
  * 
  */
