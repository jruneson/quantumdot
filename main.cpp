#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>

#include "simulation.hpp"
#include "parameters.hpp"



int main()
{
	//std::vector<double> betas = {10000};
	std::vector<double> betas = {6500,10000,15000,20000,50000,100000,200000};
	std::vector<double> taus = {5000};
	//std::vector<double> taus = {500,1000,2000,3000,5000,10000};
	//std::vector<double> taus = {2000, 5000, 10000, 20000, 30000, 50000};
	Parameters params;
	params.read_file("configuration.cfg");
	std::ofstream results_file("results.dat");
	results_file.precision(8);
	results_file << "%tau\tbeta";
	for(int id : params.to_measure)
		results_file << "\tObsId " << id << "\t\tError\t";
	results_file << std::endl;
	for(double beta : betas)
	{
		params.beta = beta;
		for(double tau : taus)
		{
			params.tau = tau;
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
 x centroid virial
 x fermion and boson
 x make P depend on beta
 x units
 * maybe separate out error estimation from exchange factor
 * read config.xyz in/out
 x check conv with tau 
 * clean up logfiles
 * use python instead of matlab, to also run program from there
 * Bias
 * Test
 * Spline
 * Test
 * read from XML-file or similar
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
