
#include <vector>

#ifndef PARAMS_HPP
#define PARAMS_HPP

class Parameters{
public:
	//read from file
	
	
	//independent parameters
	const int num_parts=2;
	const int num_beads=30;
	const int dim=3;
	const double dt=0.1;
	const double beta = 5;
	const int max_blocks=20;
	const int num_steps = 1000; //per block
	const int steps_per_sample = 5;
	
	const double hbar = 1;
	const double mass = 1;
	const double curvature = 1;
	const double charge = 1;
	const double diel_const = 1;
	const double length_scale = 1;
	
	//dependent parameters
	const int num_samples = (int) num_steps / steps_per_sample;
	const double spring_const = num_beads*mass/(hbar*hbar*beta*beta);
	
	
	//potential
	
	//observables to measure
	const std::vector<int> to_measure = {22};
	const std::vector<int> to_print_every_sample = {22};
};



#endif
