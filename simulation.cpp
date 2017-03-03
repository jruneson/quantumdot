

#include "simulation.hpp"


//Simulation::Simulation() : dt_(0.1), dt_2m_(0.05){}

Simulation::Simulation(const Parameters& params)
	: dt(params.dt), num_parts(params.num_parts), interac(Interaction(params)),
	  length_scale(params.length_scale), max_blocks(params.max_blocks),
	  num_samples(params.num_samples), steps_per_sample(params.steps_per_sample),
	  to_print(params.to_print_every_sample)
{
	for(int n=0; n<num_parts; ++n)
		polymers.push_back(Polymer(params));
	finished = false;
	block = 0;
	time = 0;
	bar_width=70;
	progress=0;
	for(int id : params.to_measure)
		obs.insert(std::pair<int,Observable>(id,Observable(id)));
}




void Simulation::setup() 
{
	for(int n=0; n<polymers.size(); ++n)
	{
		Polymer& pol = polymers[n];
		for(int bead=0; bead<pol.num_beads; ++bead)
		{
			for(int d=0; d<pol[0].size(); ++d)
			{
				pol[bead][d] = length_scale * ((double) bead/pol.num_beads - 0.5) * ((double) n/num_parts + 0.1);
				pol.vels[bead][d] = 0.0;
				pol.forces[bead][d] = 0.0;
			}
		}
	}
	std::cout << polymers.size() << "\t" << polymers[0].num_beads << "\t" << polymers[0][0].size() 
			<< "\t" << max_blocks << std::endl;
	std::time_t t = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
	logfile.open("logfile"); //add time
	logfile << std::ctime(&t) << std::endl;
	logfile.precision(8);
}

void Simulation::run() 
{
	while(!finished)
		run_block();
	stop();
}

void Simulation::run_block()
{
	reset_obs();
	for(int s=0; s<num_samples; ++s)
	{
		for(int step=0; step<steps_per_sample; ++step)
		{
			//thermostat
			verlet_step();
			//thermostat
		}
		time += dt*steps_per_sample;
		measure();
	}
	update_avgs();
	++block;
	update_screen();
	print_to_file();
	if(block>=max_blocks)
		finished = true;
}

void Simulation::verlet_step()
{
	for(Polymer& pol : polymers)
	{
		pol.update_vels();
		pol.move();
	}
	interac.update_forces(polymers);
	for(Polymer& pol: polymers)
		pol.update_vels();

}

void Simulation::reset_obs()
{
	for(auto& ob : obs)
		ob.second.set_zero();
}

void Simulation::measure()
{
	for(auto& pair : obs)
		pair.second.measure(polymers,interac);
		//add print to separate files for id:s in to_print
}

void Simulation::update_avgs()
{
	for(auto& ob : obs)
		ob.second.update_avg(num_samples);
}

void Simulation::update_screen()
{
	std::cout << "[";
	progress = (int) bar_width*block/max_blocks;
	for(int i=0; i<bar_width; ++i)
	{
		if(i<progress) std::cout << "=";
		else if(i==progress) std::cout << ">";
		else std::cout << " ";
	}
	std::cout << "]" << "\r";
	std::cout.flush();
}


void Simulation::print_to_file()
{
	logfile << block;
	for(auto& ob : obs)
		logfile <<	"\t" << ob.second.get_avg()/block;
	logfile << std::endl;
}


void Simulation::stop()
{
	logfile << "Name\t\t\tValue\t\tError" << std::endl;
	for(auto& pair : obs)
	{
		Observable& ob = pair.second;
		ob.normalize_avg(block);
		logfile << ob.get_name() <<"\t" << ob.get_avg() << "\t" << ob.std_dev(block) << std::endl;
	}
	std::cout << std::endl;
}
	


