
#include <vector>
#include <string>
#include <fstream>

#ifndef PARAMS_HPP
#define PARAMS_HPP

class Parameters{
public:
	//read from file
	void read_file(std::string);
	void calculate_dependencies();
	
	//independent parameters
	int num_parts;
	int dim;
	double tau; //in Ha^{-1}
	double beta; //in Ha^{-1}
	
	int max_blocks;
	double total_time;
	int steps_per_sample;
	int thermalization_steps;
	int num_bins;
	int sign;
	double tolerance;
	int steps_in_highest_mode;
	bool with_thermostat;
	bool using_input_file;
	
	double hw; //in Ha
	double m_hbar2 = 3.674932248e-5; //m/hbar^2 in  a_0^{-2} meV^{-1}
	double hbar = 0.65821195; //in meV ps
	double mass; //in kg
	double charge;
	double diel_const;
	double length_scale; //in a_0
	
	//dependent parameters
	int num_beads;
	double curvature; //in Ha a_0^{-2}
	double dt_md;
	double dt_2m;
	double temperature;
	int num_samples;
	int num_steps; //per block
	double spring_const; 
	double exc_const;
	double exc_der_const;
	double kin_offset;
	double virial_offset;
	double hist_size;
	
	//potential
	
	//observables to measure
	std::vector<int> to_measure;
	std::vector<int> to_print_every_sample;
};



#endif
