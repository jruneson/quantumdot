#include <sstream>
#include <iostream>
#include <cmath>
#include "parameters.hpp"


void Parameters::read_file(std::string filename)
{
	std::ifstream file(filename);
	std::string name;
	std::string line;
	int to_bool; int tmp;
	double tmp2;
	while(std::getline(file, line))
	{
		std::istringstream iss(line);
		if(iss >> name)
		{
			if(name=="num_parts")
				iss >> num_parts;
			else if(name=="dim")
				iss >> dim;
			else if(name=="tau")
				iss >> tau;
			else if(name=="beta")
				iss >> beta;
			else if(name=="num_blocks")
				iss >> num_blocks;
			else if(name=="sampling_time")
				iss >> sampling_time;
			else if(name=="dt_md")
				iss >> dt_md;
			else if(name=="dt_md_slow")
				iss >> dt_md_slow;
			else if(name=="dt_sample")
				iss >> dt_sample;
			else if(name=="non_sampling_time")
				iss >> non_sampling_time; 
			else if(name=="thermalization_time")
			{
				iss >> tmp2;
				thermalization_steps = tmp2/dt_md;
			}
			else if(name=="num_bins")
				iss >> num_bins;
			else if(name=="sign")
				iss >> sign;
			else if(name=="with_thermostat")
			{
				iss >> to_bool;
				with_thermostat = (to_bool != 0);
			}
			else if(name=="metad_on")
			{
				iss >> to_bool;
				metad_on = (to_bool != 0);
			}
			else if(name=="cv_id")
				iss >> cv_id;
			else if(name=="gauss_width")
				iss >> gauss_width;
			else if(name=="bias_factor")
				iss >> bias_factor;
			else if(name=="first_height")
				iss >> first_height;
			else if(name=="bias_update_time")
				iss >> bias_update_time;
			else if(name=="hw")
				iss >> hw;
			else if(name=="mass_in_m_e")
			{	
				iss >> mass;	
				m_hbar2 *= mass;
				mass *= 9.10938356e-31;
			}
			else if(name=="charge_in_e")
				iss >> charge;
			else if(name=="diel_const_rel")
				iss >> diel_const;
			else if(name=="to_measure")
				while(iss >> tmp)
					to_measure.push_back(tmp);
			else if(name=="to_print_every_sample")
				while(iss >> tmp)
					to_print_every_sample.push_back(tmp);
			else if(name=="lj_length")
				iss >> lj_length;
			else if(name=="lj_energy")
				iss >> lj_energy;
			else if(name=="interaction_id")
				iss >> interaction_id;
			else if(name=="wall_id")
				iss >> wall_id;
			else if(name=="wall_pos")
				iss >> wall_pos;
			else if(name=="wall_energy")
				iss >> wall_energy;
			else
				std::cout << name << " is not an allowed parameter" << std::endl;
		}
	}
	if(dt_md_slow==0)
		dt_md_slow = dt_md;
					
				
	/*		
	file >> name >> num_parts;
	file >> name >> dim;
	file >> name >> tau;
	file >> name >> beta;
	file >> name >> num_blocks;
	file >> name >> total_time;
	file >> name >> dt_md;
	file >> name >> dt_sample;
	file >> name >> thermalization_steps;
	file >> name >> num_bins;
	file >> name >> sign;
	//file >> name >> tolerance;
	//file >> name >> steps_in_highest_mode;
	int to_bool;
	file >> name >> to_bool;
	with_thermostat = (to_bool != 0);
	//file >> name >> to_bool;
	//using_input_file = (to_bool != 0);
	file >> name >> to_bool;
	metad_on = (to_bool != 0);

	file >> name >> cv_id;
	file >> name >> gauss_width;
	file >> name >> bias_factor;
	file >> name >> first_height;
	file >> name >> bias_update_time;

	file >> name >> hw;
	file >> name >> mass;

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
		to_print_every_sample.push_back(id);*/
}

void Parameters::calculate_dependencies()
{
	num_beads = round(beta/tau);
	if(num_beads<1)
		num_beads=1;
	curvature = m_hbar2 * hw*hw;
	//dt_md = 2*M_PI * std::pow(1 + 4.0*num_beads/(hw*hw*beta*beta),-0.5) * hbar/hw
	//		* 1.0/steps_in_highest_mode; //at least 10dt within the highest freq mode
	//dt_md = 0.005;
	//steps_per_sample = round(dt_sample/dt_md);
	dt_2m = dt_md / (2.0*mass) * 5.7214765779e-26; //containing conversion factor from meV/a_0 to kg a_0 ps^{-2}
	temperature = 1.0/beta * 11.60452205;
	first_height = first_height/beta;
	//num_steps = (int) sampl_time / (dt_md * num_blocks); //per block
	//num_samples = (int) num_steps / steps_per_sample; //per block
	spring_const = num_beads*m_hbar2/(beta*beta);
	exc_const = num_beads*m_hbar2/beta;
	exc_der_const = -sign * m_hbar2/(beta*beta);
	kin_offset = num_parts*num_beads*dim/(2.0*beta);
	virial_offset = dim/(2.0*beta); //not sure about the dim factor
	length_scale = std::sqrt(1.0/(m_hbar2*hw));
	hist_size = length_scale * 4;	
}


