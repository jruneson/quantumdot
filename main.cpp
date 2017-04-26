#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>

#include "simulation.hpp"
#include "parameters.hpp"



int main(int argc, char* argv[])
{
	bool continue_sim = false;
	if(argc > 1)
	{
		if((std::string(argv[1])=="-c")||(std::string(argv[1])=="--continue"))
		{
			continue_sim = true;
		}
	}
	std::vector<double> betas = {0.05,0.1,0.2,0.3,0.5,1,2,3,4};
	//std::vector<double> betas = {2.0};
	//std::vector<double> taus = {0.1};
	//std::vector<double> betas = {1.0};
	//std::vector<double> taus = {2.0,1.0,0.5,0.4,0.3,0.2,0.15,0.1,0.067,0.05,0.04};
	//std::vector<double> taus = {1.0};
	Parameters params;
	params.read_file("configuration.cfg");
	params.calculate_dependencies();
	std::ofstream results_file;
	if(!continue_sim)
	{
		results_file.open("results.dat");
		results_file << "%tau";
		for(int id : params.to_measure)
			results_file << "\tObsId " << id << "\t\tError\t";
		results_file << std::endl;
	}
	else
		results_file.open("results.dat", std::ios_base::app);
	results_file.precision(10);
	/*for(auto beta : betas)
	{
		params.beta = beta;
		params.calculate_dependencies();*/
	Simulation sim(params, results_file,continue_sim);
	sim.setup();
	sim.run();
	//}
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
 x maybe separate out error estimation from exchange factor
 x read config.xyz in/out
 x check conv with tau 
 x clean up logfiles
 x use python instead of matlab, to also run program from there
 x Bias
 x Test
 x Spline
 x Test
 * replace gle ptr with static obj
 x Interaction, several particles
 * Check energy estimators are correct
 * Check LJ is reasonable
 * Run LJ with+wo metad B, then F
 * histogram at 3d
 */
 
 /*
  * Note that mass must be the same of all particles. Otherwise things have to be 
  * updated in GLE and force calculation (curvature is different)
  * 
  * 
  */
