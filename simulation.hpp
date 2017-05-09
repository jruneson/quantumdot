
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

	void run();
	
private:
	const int num_parts;
	std::vector<Polymer> polymers;
	const double dt_md;
	const double dt_sample;
	Interaction interac;
	const double length_scale;
	const double temperature;
	bool finished;
	double overall_time; //since thermalization is finished
	double time_sampled;
	int block; //number of measured blocks
	int samples; //number of samples taken in a block
	double sampling_time; //time to sample
	double sampling_time_per_block;
	double non_sampling_time; //only used if not cont_sim
	int num_blocks;
	const int thermalization_steps;
	//const int steps_per_sample;	
	const int num_bins;
	const double hist_size;
	GLE gle;
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
	double exc_sq; // <Gamma^2> within a block
	double exc_sq_avg; // <Gamma^2> over all blocks, for standard deviation calculation
	double e_s_sum;
	double e_s_avg; //<exp(-s)>
	double e_s_avg_sq;
	
	std::map<int,Observable> obs;	
	std::vector<double> histogram; //1D probability density. If dim>1, only the first coordinate is considered.
	std::vector<double> histogram_avg;
	std::vector<double> histogram_sq_avg;
	std::vector<std::vector<double>> histogram_1p;
	std::vector<std::vector<double>> histogram_1p_avg;
	std::vector<std::vector<double>> histogram_1p_sq_avg;
	std::vector<double> histogram_delta_e;
	double hist_de_resolution=0.5;
	double hist_de_max=50;
	double hist_de_min=-50;
	double hist_de_width;
	double hist_de_num_bins;
	
	
	int bar_width;
	int progress;
	Timer timer;
	double movie_start_time;
	double movie_end_time;
	
	std::ofstream logfile;
	std::ofstream& res_file;
	std::ofstream exc_file;
	std::ofstream cv_file;
	std::ofstream rew_factor_file;
	std::ofstream vmd_file;
	std::ofstream vmd_file2;

	void setup();
	void read_input_coords();
	void initialize_coords_simple();
	void read_old_measurements();
	void thermalize();
	void run_wo_sampling();
	void run_block();
	void verlet_step();
	//void reset_obs();
	void measure();
	void update_avgs();
	void update_histogram();
	int calc_bin(double,int,double);
	//int calc_bin_1p(double);
	void update_screen();
	void print_to_logfile();
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
