
#include <vector>
#include <iostream>
#include <fstream>
#include <chrono>
#include <ctime>

#include "parameters.hpp"
#include "polymer.hpp"
#include "interaction.hpp"
#include "bias.hpp"
#include "observable.hpp"

#ifndef SIM_HPP
#define SIM_HPP

class Simulation{
public:
	Simulation(const Parameters&);
	void setup();
	void run();
	void run_block();
	void verlet_step();
	void zero_avgs();
	void measure();
	void update_avgs();
	void update_screen();
	void print_to_file();
	void stop();
	
	double potential_energy();
	
	double coll_var() const;
	Force grad_coll_var() const;
	


private:
	const int num_parts;
	std::vector<Polymer> polymers;
	Bias bias;
	//Spline v_spline;
	//Spline f_spline;
	const double dt;
	Interaction interac;
	const double length_scale;
	bool finished;
	int block;
	const int max_blocks;
	const int num_samples;
	const int steps_per_sample;
	//int taken_samples;
	
	//double e_pot;
	//double e_pot_avg;
	//double e_pot_sq;
	
	Observable e_pot;
	
	int bar_width;
	int progress;
	
	std::ofstream logfile;

};

#endif
