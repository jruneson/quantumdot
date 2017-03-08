#include <sstream>
#include <iostream>
#include "parameters.hpp"


void Parameters::read_file(std::string filename)
{
	std::ifstream file(filename);
	std::string name;
	file >> name >> num_parts;
	file >> name >> num_beads;
	file >> name >> dim;
	file >> name >> dt;
	file >> name >> beta;
	file >> name >> max_blocks;
	file >> name >> num_steps;
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

	temperature = 1.0/beta;
	num_samples = (int) num_steps / steps_per_sample;
	spring_const = num_beads*mass/(hbar*hbar*beta*beta);
	hist_size = length_scale * 10;
	
}


