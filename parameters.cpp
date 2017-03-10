#include <sstream>
#include <iostream>
#include <cmath>
#include "parameters.hpp"


void Parameters::read_file(std::string filename)
{
	std::ifstream file(filename);
	std::string name;
	file >> name >> num_parts;
	file >> name >> num_beads;
	file >> name >> dim;
	file >> name >> beta;
	file >> name >> max_blocks;
	file >> name >> total_time;
	file >> name >> steps_per_sample;
	file >> name >> thermalization_steps;
	file >> name >> num_bins;
	int to_bool;
	file >> name >> to_bool;
	with_thermostat = (to_bool != 0);

	
	file >> name >> hbar;
	file >> name >> mass;
	file >> name >> curvature;
	file >> name >> charge;
	file >> name >> diel_const;
	file >> name >> length_scale;
	
	std::string line;
	int id;
	getline(file, line);
	getline(file, line);
	std::istringstream iss(line);
	iss >> name;
	while(iss >> id)
		to_measure.push_back(id);
	getline(file, line);
	std::istringstream iss2(line);
	iss2 >> name;
	while(iss2 >> id)
		to_print_every_sample.push_back(id);
}

void Parameters::calculate_dependencies()
{
	dt_md = 2*M_PI * std::pow(curvature/mass + 3.0*num_beads/(hbar*hbar*beta*beta),-0.5) * 0.05; //at least 10dt within the highest freq mode
	//dt_sample = dt_md * steps_per_sample;
	temperature = 1.0/beta;
	num_steps = (int) total_time / (dt_md * max_blocks); //per block
	num_samples = (int) num_steps / steps_per_sample; //per block
	spring_const = num_beads*mass/(hbar*hbar*beta*beta);
	hist_size = length_scale * 10;
	
}


