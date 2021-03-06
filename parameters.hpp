
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
	bool connected;
	
	int num_blocks;
	double sampling_time;
	double dt_md;
	double dt_md_slow;
	double dt_sample;
	double non_sampling_time; //is set to zero in simulation.cpp if simulation is continuing a previous one
	int thermalization_time;
	int num_bins;
	int num_bins_2d;
	double hist_size_in_r_star;
	int sign;
	int spin; //In halves 
	bool spin_proj; //Tells if "spin" means z projection or total spin
	bool with_thermostat;
	bool metad_on;
	bool more_output;
	
	int cv_id;
	double gauss_width;
	double bias_factor;
	double first_height_in_kBT;
	double first_height;
	double bias_update_time;
	int biased_graph; //connected graph if N=2
	int reference_graph; //disconnected graph if N=2
	double permutation_trial_time;
	bool allow_perm_switch;
	
	double wigner_parameter; //electrostatic energy scale / hw
	double hw; //in meV
	double hwx;
	double hwy;
	double m_hbar2 = 0.013234213;//3.674932248e-5; //m/hbar^2 in  nm^{-2} meV^{-1}
	double hbar = 0.65821195; //in meV ps
	double mass; //in kg
	double charge; //in e
	double diel_const; // 4\pi\eps0\eps_r in meV^{-1} e^2 a0^{-1}
	double screening_factor;
	double electrost_factor; //e^2/diel_const * screening_factor
	
	int interaction_id;
	double lj_length; //Lennard-Jones sigma, in a_0
	double lj_energy; //Lennard-Jones epsilon, in meV
	
	int wall_id; //0=no wall, 1=upper, 2=lower
	double wall_pos;
	double wall_energy;
	
	double cv_hist_max;
	double cv_hist_min;
	double cv_hist_res;
	
	//dependent parameters
	int num_beads;
	double curvature; //in Ha a_0^{-2}
	double curvature_x;
	double curvature_y; //not used it d<2
	double curvature_z; //not used if d<3
	int steps_per_sample;
	double dt_2m;
	double temperature;
	int num_samples;
	int num_steps; //per block
	double spring_const; 
	double exc_const;
	double mass_factor;
	double exc_der_const;
	double kin_offset;
	double virial_offset;
	double length_scale; //r_star in nm
	double hist_size;
	double spline_step;

	
	
	//observables to measure
	std::vector<int> to_measure;
	std::vector<int> to_print_every_sample;
};



#endif
