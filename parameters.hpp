
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
	int num_parts;
	int num_beads;
	int dim;
	double dt;
	double beta;
	int max_blocks;
	int num_steps; //per block
	int steps_per_sample;
	bool with_thermostat;
	int thermalization_steps;
	int num_bins;
	
	double hbar;
	double mass;
	double curvature;
	double charge;
	double diel_const;
	double length_scale;
	
	//dependent parameters
	double temperature;
	int num_samples;
	double spring_const;
	double hist_size;
	
	//potential
	
	//observables to measure
	std::vector<int> to_measure;
	std::vector<int> to_print_every_sample;
};



#endif
