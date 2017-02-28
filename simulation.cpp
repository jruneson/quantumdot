#include "simulation.hpp"


//Simulation::Simulation() : dt_(0.1), dt_2m_(0.05){}

Simulation::Simulation(const Parameters& params)
	: dt_(params.dt), num_parts(params.num_parts)
{
	for(int n=0; n<num_parts; ++n)
		polymers.push_back(Polymer(params));
}






void Simulation::setup() 
{}

void Simulation::run() 
{}

void Simulation::run_block()
{}

void Simulation::verlet_step()
{
	for(Polymer& pol : polymers)
	{
		pol.update_vels();
		pol.move();
		i.update_forces(pol)
		pol.update_vels();
	}
}

void Simulation::print_to_file()
{}


void Simulation::stop()
{}
	
