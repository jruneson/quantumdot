
#include <vector>

#include "parameters.hpp"
#include "polymer.hpp"
#include "interaction.hpp"
#include "bias.hpp"

#ifndef SIM_HPP
#define SIM_HPP

class Simulation{
public:
	Simulation(const Parameters&);
	void setup();
	void run();
	void run_block();
	void verlet_step();
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
	const double dt_;
	Interaction i;
};

#endif
