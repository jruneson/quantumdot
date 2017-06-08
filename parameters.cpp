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
	dt_md_slow = 0;
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
			else if(name=="connected")
			{
				iss >> to_bool;
				connected = (to_bool != 0);
				//if(connected)
				//	throw std::runtime_error("It's not safe to set connected=1 at the moment");
			}
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
				iss >> thermalization_time;
			}
			else if(name=="num_bins")
				iss >> num_bins;
			else if(name=="num_bins_2d")
				iss >> num_bins_2d;
			else if(name=="sign")
				iss >> sign;
			else if(name=="spin_times_two")
			{
				iss >> spin;
				if( (num_parts % 2) != (spin % 2)) //If num_parts and spin are not both even or both odd
				//if(((num_parts % 2) && (!(total_spin % 2))) || ((!(num_parts % 2)) && (total_spin%2))) //If number of particles is odd and spin is even
					throw std::runtime_error("Invalid total spin!");
				if(spin > num_parts)
					throw std::runtime_error("Total spin is too high!");
			}
			else if(name=="use_projection")
			{
				iss >> to_bool;
				spin_proj = (to_bool != 0);
			}
			else if(name=="biased_graph")
				iss >> biased_graph;
			else if(name=="with_thermostat")
			{
				iss >> to_bool;
				with_thermostat = (to_bool != 0);
			}
			else if(name=="metad_on")
			{
				iss >> to_bool;
				metad_on = (to_bool != 0);
				if(sign==0)
					metad_on=false;
			}
			else if(name=="cv_id")
				iss >> cv_id;
			else if(name=="gauss_width")
				iss >> gauss_width;
			else if(name=="bias_factor")
				iss >> bias_factor;
			else if(name=="first_height_in_kBT")
				iss >> first_height_in_kBT;
			else if(name=="bias_update_time")
				iss >> bias_update_time;
			else if(name=="wigner_parameter")
				iss >> wigner_parameter;
			else if(name=="hwx")
				iss >> hwx;
			else if(name=="hwy")
				iss >> hwy;
			else if(name=="mass_in_m_e")
			{	
				iss >> mass;	
				m_hbar2 *= mass;
				mass *= 9.10938356e-31;
			}
			else if(name=="charge_in_e")
				iss >> charge;
			else if(name=="diel_const_rel")
			{
				iss >> diel_const;
				diel_const *= 6.94461548e-4;//3.674932e-5;
			}
			else if(name=="screening_factor")
				iss >> screening_factor;
			else if(name=="to_measure")
				while(iss >> tmp)
					to_measure.push_back(tmp);
			else if(name=="to_print_every_sample")
				while(iss >> tmp)
					to_print_every_sample.push_back(tmp);
			else if(name=="more_output")
			{	
				iss >> to_bool;
				more_output = (to_bool != 0);
			}
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
}


void Parameters::calculate_dependencies()
{
	num_beads = round(beta/tau);
	if(num_beads<1)
		num_beads=1;
	electrost_factor = screening_factor*charge*charge/diel_const;
	//hw = electrost_factor*electrost_factor*m_hbar2/(wigner_parameter*wigner_parameter);
	hw = std::sqrt((hwx*hwx+hwy*hwy)/2);
	wigner_parameter = electrost_factor*std::sqrt(m_hbar2/hw);
	//electrost_factor = 0;
	//hwx = 5.1;//4.23;
	//hwy = 5.1;//5.84;
	//electrost_factor *= screening_factor;
	curvature_x = m_hbar2 * hwx*hwx;
	curvature_y = m_hbar2 * hwy*hwy;
	//dt_md = 2*M_PI * std::pow(1 + 4.0*num_beads/(hw*hw*beta*beta),-0.5) * hbar/hw
	//		* 1.0/steps_in_highest_mode; //at least 10dt within the highest freq mode
	//dt_md = 0.005;
	//steps_per_sample = round(dt_sample/dt_md);
	dt_2m = dt_md / (2.0*mass) * 1.60217662e-28;//5.7214765779e-26; //containing conversion factor from meV/nm to kg nm ps^{-2}
	temperature = 1.0/beta * 11.60452205;
	first_height = first_height_in_kBT/beta;
	//num_steps = (int) sampl_time / (dt_md * num_blocks); //per block
	//num_samples = (int) num_steps / steps_per_sample; //per block
	spring_const = num_beads*m_hbar2/(beta*beta);
	exc_const = num_beads*m_hbar2/(beta);
	exc_der_const = sign * num_beads*m_hbar2/(beta*beta);
	kin_offset = num_parts*num_beads*dim/(2.0*beta);
	//if(connected)
	virial_offset = dim/(2.0*beta);
	//else
	//	virial_offset = num_parts*dim/(2.0*beta);
	length_scale = std::sqrt(1.0/(m_hbar2*hw));
	hist_size = length_scale * 4;	
	
	//std::cout << hwx << "\t" << hwy << "\t" << wigner_parameter << "\t" << virial_offset << std::endl;
}


