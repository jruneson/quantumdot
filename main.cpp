#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>

#include "simulation.hpp"
#include "parameters.hpp"



int main()
{
	std::vector<double> betas = {5.0};
	std::vector<double> taus = {1.0,0.5,0.2,0.15,0.1,0.075};
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
 * why are fluctuations of e_pot larger for smaller beta?
 * fermion and boson
 * make P depend on beta
 * units
 * maybe separate out error estimation from exchange factor
 * read config.xyz in/out
 * check conv with tau 
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
