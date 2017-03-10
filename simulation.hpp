
#include <vector>
#include <iostream>
#include <fstream>
#include <chrono>
#include <ctime>
#include <map>
#include <algorithm>

#include "parameters.hpp"
#include "polymer.hpp"
#include "interaction.hpp"
#include "bias.hpp"
#include "observable.hpp"
#include "GLE-OPT.hpp"
#include "timer.hpp"

#ifndef SIM_HPP
#define SIM_HPP



class Simulation{
public:
	Simulation(const Parameters&, std::ofstream&);
	void setup();
	void thermalize();
	void run();
	void run_block();
	void verlet_step();
	void reset_obs();
	void measure();
	void update_avgs();
	void update_histogram();
	int calc_bin(const double&);
	void update_screen();
	void print_to_file();
	void stop();
		
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
	const double temperature;
	bool finished;
	double time;
	int block;
	const double total_time;
	const int max_blocks;
	const int num_samples;
	const int thermalization_steps;
	const int steps_per_sample;	
	const int num_bins;
	const double hist_size;
	GLE* gle;
	const bool thermostat_on;
	
	std::map<int,Observable> obs;	
	std::vector<double> histogram; //1D probability density. If dim>1, only the first coordinate is considered.
	std::vector<double> histogram_avg;
	std::vector<double> histogram_sq_avg;
	
	int bar_width;
	int progress;
	Timer timer;
	
	std::ofstream logfile;
	std::ofstream& res_file;

};



#endif
