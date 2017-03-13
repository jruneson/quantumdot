
#include <stdexcept>
#include "simulation.hpp"


//Simulation::Simulation() : dt_(0.1), dt_2m_(0.05){}

Simulation::Simulation(const Parameters& params, std::ofstream& _res_file)
	: dt(params.dt_md), num_parts(params.num_parts), interac(Interaction(params)),
	  length_scale(params.length_scale), max_blocks(params.max_blocks),
	  num_samples(params.num_samples), steps_per_sample(params.steps_per_sample),
	  num_bins(params.num_bins), hist_size(params.hist_size),
	  temperature(params.temperature), thermalization_steps(params.thermalization_steps),
	  thermostat_on(params.with_thermostat), res_file(_res_file),
	  total_time(params.total_time),tolerance(params.tolerance),
	  beta(params.beta), sign(params.sign), exc_const(params.exc_const)
{
	for(int n=0; n<num_parts; ++n)
		polymers.push_back(Polymer(params));
	finished = false;
	block = 0;
	time = 0;
	bar_width=70;
	progress=0;
	std::vector<int> to_print = params.to_print_every_sample;
	for(int id : params.to_measure)
	{
		auto pair = obs.insert(std::pair<int,Observable>(id,Observable(id,params)));
		if(std::find(to_print.begin(), to_print.end(), id) != to_print.end())
			pair.first->second.set_print_on();		
	}
	histogram.assign(params.num_bins,0);
	histogram_avg.assign(params.num_bins,0);
	histogram_sq_avg.assign(params.num_bins,0);
}




void Simulation::setup() 
{
	timer.start();
	for(int n=0; n<polymers.size(); ++n)
	{
		Polymer& pol = polymers[n];
		for(int bead=0; bead<pol.num_beads; ++bead)
		{
			for(int d=0; d<pol[0].size(); ++d)
			{
				pol[bead][d] = length_scale * ((double) bead/pol.num_beads - 0.5) * ((double) n/num_parts + 0.1);
				pol.vels[bead][d] = 0.5;
				pol.forces[bead][d] = 0.0;
			}
		}
	}
	std::cout << "time = " << total_time << "\tP = " << polymers[0].num_beads 
				<< "\t dt = " << dt << "\t beta = " << beta << std::endl;
	std::time_t t = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
	logfile.open("logfile_P"+std::to_string(polymers[0].num_beads));
	logfile << std::ctime(&t);
	logfile << "P=" << polymers[0].num_beads << " dim=" << polymers[0][0].size() << " dt=" << dt 
			<< " maximal time=" << total_time
			//<< " actual time=" << dt*max_blocks*num_samples*steps_per_sample
			<< " samples=" << max_blocks*num_samples << " T=" << temperature << std::endl << std::endl;
	logfile << "Block";
	for(auto& pair : obs)
		logfile << "\t" << pair.second.get_name();
	logfile << std::endl;
	logfile.precision(8);
	try 
	{
		gle = new GLE(polymers, dt, temperature, polymers[0].mass, polymers[0].num_beads, 
					num_parts, polymers[0][0].size(),thermostat_on);
	}
	catch(const std::exception& e)
	{
		std::cout << "exception: " << e.what() << std::endl;
		return;
	}
		
	thermalize();
	update_exc();
}

void Simulation::thermalize()
{
	for(int step=0; step<thermalization_steps; ++step)
	{
		gle->run();
		verlet_step();
		gle->run();
	}
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
			gle->run();
			verlet_step();
			gle->run();
		}
		time += dt*steps_per_sample;
		update_exc();
		measure();
	}
	update_avgs();
	++block;
	print_to_file();
	update_screen();
	if((block>=max_blocks)||(converged()))
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
		pair.second.measure(polymers,interac,time,exchange_factor);
	update_histogram();
}

void Simulation::update_avgs()
{
	double tmp = exc_sum/num_samples;
	exc_avg = (exc_avg*block + tmp) / (block+1);
	exc_avg_sq = (exc_avg_sq*block + tmp*tmp)/(block+1);
	exc_sum = 0;
	for(auto& ob : obs)
		ob.second.update_avg(num_samples);
	double hist_norm=0;
	for(int bin=0; bin<num_bins; ++bin)
		hist_norm += histogram[bin];
	hist_norm *= hist_size/num_bins;
	for(int bin=0; bin<num_bins; ++bin)
	{
		histogram_avg[bin] += histogram[bin]/hist_norm;
		histogram_sq_avg[bin] += std::pow(histogram[bin]/hist_norm,2);
		histogram[bin]=0;
	}
}

void Simulation::update_histogram()
{
	for(const auto& pol : polymers)
	{
		for(int bead=0; bead<pol.num_beads; ++bead)
		{
			int bin = calc_bin(pol[bead][0]);
			if(bin>=0 && bin<num_bins)
				histogram[bin]++;
		}
	}
}

int Simulation::calc_bin(double coord)
{
	return round((0.5*hist_size + coord)*num_bins/hist_size);
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
		logfile <<	"\t" << ob.second.get_value()/num_samples;
	logfile << std::endl;
}

bool Simulation::converged()
{
	/*
	if(block>=5)
	{
		bool tmp = true;
		for(const auto& pair : obs)
		{
			const auto& ob = pair.second;
			tmp*=(abs(ob.std_dev()/ob.get_avg()) <= tolerance);
		}
		return tmp;
	}*/
	return false;
}

void Simulation::stop()
{
	logfile << "Name\t\tValue\t\tError with " << block << " blocks" << std::endl;
	res_file << polymers[0].num_beads << "\t" << beta;
	for(const auto& pair : obs)
	{
		const auto& ob = pair.second;
		logfile << ob.get_name() << "\t" << ob.get_avg(exc_avg) << "\t" 
				<< ob.std_dev(exc_avg,exc_avg_sq) << std::endl;
		res_file << "\t" << ob.get_avg(exc_avg) << "\t" << ob.std_dev(exc_avg,exc_avg_sq);
	}
	res_file << std::endl;
	std::ofstream hist_file("Prob_distribution.dat");
	for(int bin = 0; bin<num_bins; ++bin)
	{
		histogram_avg[bin] /= block;
		histogram_sq_avg[bin] /= block;
		hist_file << hist_size*((double) bin /num_bins - 0.5) << "\t" << histogram_avg[bin]
					<< "\t" << std::sqrt((histogram_sq_avg[bin]-std::pow(histogram_avg[bin],2))/(block-1)) 
					<< std::endl;
	}
	logfile.close();
	delete gle;
	timer.stop();
	std::cout << std::endl << timer.duration() << " s" << std::endl;
}

void Simulation::update_exc()
{
	double tmp = 0;
	for(int bead=0; bead<polymers[0].num_beads; ++bead)
		tmp += std::exp( - exc_const * (polymers[0][bead]-polymers[1][bead])*(polymers[0][bead+1]-polymers[1][bead+1]));
	exchange_factor = 1.0 + sign*tmp/polymers[0].num_beads;
	exc_sum += exchange_factor;
}
