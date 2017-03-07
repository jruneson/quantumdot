
#include <vector>
#include <string>
#include <fstream>

#ifndef PARAMS_HPP
#define PARAMS_HPP

class Parameters{
public:
	//read from file
	void read_file(std::string);
	
	//independent parameters
	int num_parts=2;
	int num_beads=30;
	int dim=1;
	double dt=0.01;
	double beta = 5.0;
	int max_blocks=20;
	int num_steps = 10000; //per block
	int steps_per_sample = 5;
	bool with_thermostat = true;
	int thermalization_steps = 5000;
	int num_bins = 2000;
	
	double hbar = 1;
	double mass = 1;
	double curvature = 1;
	double charge = 1;
	double diel_const = 1;
	double length_scale = 1;
	
	//dependent parameters
	double temperature;
	int num_samples;
	double spring_const;
	double hist_size;
	
	//potential
	
	//observables to measure
	const std::vector<int> to_measure = {10, 22};
	const std::vector<int> to_print_every_sample = {10, 22};
};



#endif
