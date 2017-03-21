#include <sstream>
#include <iostream>
#include <cmath>
#include "parameters.hpp"


void Parameters::read_file(std::string filename)
{
	std::ifstream file(filename);
	std::string name;
	file >> name >> num_parts;
	file >> name >> dim;
	file >> name >> tau;
	file >> name >> beta;
	file >> name >> max_blocks;
	file >> name >> total_time;
	file >> name >> steps_per_sample;
	file >> name >> thermalization_steps;
	file >> name >> num_bins;
	file >> name >> sign;
	file >> name >> tolerance;
	file >> name >> steps_in_highest_mode;
	int to_bool;
	file >> name >> to_bool;
	with_thermostat = (to_bool != 0);
	file >> name >> to_bool;
	using_input_file = (to_bool != 0);
	file >> name >> to_bool;
	metad_on = (to_bool != 0);

	file >> name >> cv_id;
	file >> name >> gauss_width;
	file >> name >> bias_factor;
	file >> name >> bias_update_time;

	file >> name >> hw;
	file >> name >> mass;
	m_hbar2 *= mass;
	mass *= 9.10938356e-31;
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
	num_beads = round(beta/tau);
	if(num_beads<1)
		num_beads=1;
	curvature = m_hbar2 * hw*hw;
	dt_md = 2*M_PI * std::pow(1 + 4.0*num_beads/(hw*hw*beta*beta),-0.5) * hbar/hw
			* 1.0/steps_in_highest_mode; //at least 10dt within the highest freq mode
	//change dt_md depending on what observables are measured : 0.02 for ekin, 0.1 or 0.05 otherwise
	dt_2m = dt_md / (2.0*mass) * 5.7214765779e-26; //containing conversion factor from meV/a_0 to kg a_0 ps^{-2}
	temperature = 1.0/beta * 11.60452205;
	first_height = 1.0/beta;
	num_steps = (int) total_time / (dt_md * max_blocks); //per block
	num_samples = (int) num_steps / steps_per_sample; //per block
	spring_const = num_beads*m_hbar2/(beta*beta);
	exc_const = num_beads*m_hbar2/beta;
	exc_der_const = -sign * m_hbar2/(beta*beta);
	kin_offset = num_beads*dim/(2.0*beta);
	virial_offset = dim/(2.0*beta); //not sure about the dim factor
	hist_size = length_scale * 500;	
	
	std::cout << "gamma: " << bias_factor << std::endl;
}


