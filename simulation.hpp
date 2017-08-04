
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <ctime>
#include <map>
#include <algorithm>
#include <random>

#include "parameters.hpp"
#include "polymer.hpp"
#include "interaction.hpp"
#include "bias.hpp"
#include "observable.hpp"
#include "GLE-OPT.hpp"
#include "timer.hpp"
#include "graph.hpp"

#ifndef SIM_HPP
#define SIM_HPP



class Simulation{
public:
	Simulation(const Parameters&, std::ofstream&, bool, const std::vector<Graph>&);

	void run();
	
private:
	const int num_parts;
	const int num_beads;
	const int dim;
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
	const int thermalization_time;
	const int num_bins;
	const double hist_size;
	GLE* gle;
	const bool thermostat_on;
	Bias bias;
	const double bias_update_time;
	double bias_update_counter;
	double permutation_switch_counter;
	const double permutation_trial_time;
	int iteration_nbr;
	const bool cont_sim; //continue a previous simulation
	const bool more_output;
	bool printed_warning;
	const bool allow_perm_switch;
	const double wigner_parameter;
	
	const double beta;
	const double tau;
	const int sign;
	const double exc_const;
	const int spin;
	
	double exchange_factor;
	double pos_weight;
	double neg_weight;
	double exc_sum;
	double exc_last_block;
	double exc_avg;
	double exc_avg_sq;
	double exc_sq; // <W^2> within a block
	double exc_sq_avg; // <W^2> over all blocks, for standard deviation calculation
	double e_s;
	double e_s_sum;
	double e_s_avg; //<exp(-s)>
	double e_s_avg_sq;
	double sgn;
	double sgn_sum;
	double sgn_last_block;
	double sgn_avg;
	double sgn_avg_sq;
	
	std::vector<Graph> graphs;
	int current_graph_id;
	std::map<int,Observable> obs;	
	std::vector<double> histogram; //pair correlation function. 
	std::vector<double> histogram_avg;
	std::vector<double> histogram_sq_avg;
	std::vector<std::vector<double>> histogram_1d; //1d probability density. One vector per dimenstion.
	std::vector<std::vector<double>> histogram_1d_avg;
	std::vector<std::vector<double>> histogram_1d_sq_avg;
	std::vector<std::vector<double>> histogram_2d;
	std::vector<std::vector<double>> histogram_2d_avg;
	std::vector<std::vector<double>> histogram_2d_sq_avg;
	std::vector<std::vector<double>> pair_distr_2d;
	std::vector<std::vector<double>> pair_distr_2d_avg;
	std::vector<std::vector<double>> pair_distr_2d_sq_avg;
	std::vector<std::vector<double>> pair_distr_1d_proj;
	std::vector<std::vector<double>> pair_distr_1d_proj_avg;
	std::vector<std::vector<double>> pair_distr_1d_proj_sq_avg;
	
	const int num_bins_2d;
	std::vector<double> cv_hist;
	std::vector<double> exc_fac_hist;
	std::vector<double> weight_en_hist;
	const double cv_hist_min;
	const double cv_hist_max;
	double cv_hist_width;
	const double cv_hist_res; //resolution, i.e. bin size
	int cv_hist_num_bins;
	std::vector<double> histogram_delta_e;
	std::vector<double> hist_c;
	const double hist_1d_min;
	const double hist_size_1d;
	double hist_c_min=-3;
	double hist_c_max=0;
	double hist_c_resolution=0.1;
	double hist_c_num_bins;
	
	
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
	std::ofstream file_fsum;
	
	std::random_device rd;
	std::mt19937 mt;
	std::uniform_real_distribution<double> uni_distr;
	std::uniform_int_distribution<int> int_distr;

	void setup();
	void read_input_coords();
	void initialize_coords_simple();
	void read_old_measurements();
	void thermalize();
	void run_wo_sampling();
	void run_block();
	void verlet_step();
	void measure();
	void update_avgs();
	void update_histogram();
	int calc_bin(double,int,double);
	void update_screen();
	void print_to_logfile();
	double fermi_dirac(double);
	void print_vmd();
	void stop();
	void print_config();
	void update_exc(bool);
	double exc_exponent(int) const;
	double simple_uncertainty(double,double) const;
	double weighted_uncertainty(double,double) const;
	void try_permutation_change();
		
};



#endif
