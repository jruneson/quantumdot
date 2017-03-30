
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
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
	Simulation(const Parameters&, std::ofstream&, bool);
	void setup();
	void run();
	
private:
	const int num_parts;
	std::vector<Polymer> polymers;
	const double dt;
	Interaction interac;
	const double length_scale;
	const double temperature;
	bool finished;
	double time;
	int block;
	double total_time;
	int num_blocks;
	const int num_samples;
	const int thermalization_steps;
	const int steps_per_sample;	
	const int num_bins;
	const double hist_size;
	GLE* gle;
	//const bool using_input_file;
	const bool thermostat_on;
	Bias bias;
	const double bias_update_time;
	double bias_update_counter;
	int iteration_nbr;
	const bool cont_sim; //continue a previous simulation
	//const double tolerance;
	
	const double beta;
	const double tau;
	const int sign;
	const double exc_const;
	
	double exchange_factor;
	double exc_sum;
	double exc_avg;
	double exc_avg_sq;
	
	std::map<int,Observable> obs;	
	std::vector<double> histogram; //1D probability density. If dim>1, only the first coordinate is considered.
	std::vector<double> histogram_avg;
	std::vector<double> histogram_sq_avg;
	
	int bar_width;
	int progress;
	Timer timer;
	double movie_start_time=1000;
	double movie_end_time=1010;
	
	std::ofstream logfile;
	std::ofstream& res_file;
	std::ofstream exc_file;
	std::ofstream cv_file;
	std::ofstream rew_factor_file;
	std::ofstream vmd_file;
	std::ofstream vmd_file2;

	void read_input_coords();
	void initialize_coords_simple();
	void read_old_measurements();
	void thermalize();
	void run_block();
	void verlet_step();
	void reset_obs();
	void measure();
	void update_avgs();
	void update_histogram();
	int calc_bin(double);
	void update_screen();
	void print_to_file();
	void print_vmd();
	void stop();
	void print_config();
	void update_exc();
	double simple_uncertainty(double,double) const;
	double weighted_uncertainty(double,double) const;
		
	//double coll_var() const;
	//Force grad_coll_var() const;

};



#endif
