
#include "observable.hpp"

/*
Observable::Observable()
{
	value,avg,avg_sq=0;
}*/

Observable::Observable(int _id, double _beta) : id(_id), beta(_beta)
{ 
	value, avg, avg_sq = 0;
	print = false;
}


void Observable::set_zero()
{
	value = 0;
}

void Observable::update_avg(const int& num_samples)
{
	double tmp = value/num_samples;
	avg += tmp;
	avg_sq += tmp*tmp;
}

double Observable::get_avg() const
{
	return avg;
}

double Observable::get_value() const
{
	return value;
}

void Observable::normalize_avg(const int& num_blocks)
{
	//std::cout << avg;
	avg /= num_blocks;
	avg_sq /= num_blocks;
	//std::cout << "\t" << avg << std::endl;
}

double Observable::std_dev(const double& n) const
{
	return std::sqrt((avg_sq - avg*avg)*n/(n-1.0)); //note that avg_sq is not normalized
}

/*void Observable::operator+=(double to_add)
{
	value += to_add;
}*/

std::string Observable::get_name()
{
	switch(id)
	{
		case 0:
			return "Pot_energy";
		case 1:
			return "Kinetic_energy";
		case 2:
			return "Total_energy";
		case 3:
			return "Kin_en_virial";
		case 10:
			return "X_coord_n1p1";
		case 20:
			return "Pot_energy_cl";
		case 21:
			return "Kin_energy_cl";
		case 22:
			return "Total_energy_cl";
		default:
			return "Unknown observable";
	}
}

void Observable::set_print_on()
{
	print=true;
	file.open(get_name()+".dat");
	file.precision(8);
}

void Observable::measure(const std::vector<Polymer>& polymers, Interaction& interac, const double& time)
{
	double tmp = 0;
	for(int n=0; n<polymers.size(); ++n)
	{
		const Polymer& pol = polymers[n];
		double tmp2 = 0;
		switch(id)
		{
			case 0:
				tmp2 = potential_energy(pol,interac);
				break;
			case 1:
				tmp2 = kinetic_energy(pol,interac);
				break;
			case 2:
				tmp2 = total_energy(pol,interac);
				break;
			case 3:
				tmp2 = kinetic_energy_virial(pol,interac);
				break;
			case 10:
				if(n==0)
					tmp2 = pol[0][0]; //x-coordinate of bead 0
				break;
			case 20:
				tmp2 = potential_energy_cl(pol,interac);
				break;
			case 21:
				tmp2 = kinetic_energy_cl(pol);	
				break;
			case 22:
				tmp2 = total_energy_cl(pol,interac);
				break;
		}
		tmp += tmp2;
	}
	value += tmp;
	if(print)
		print_measure(tmp,time);
}

void Observable::print_measure(const double& measured_value, const double& time)
{
	file << time << "\t" << measured_value << std::endl;
}





double Observable::potential_energy(const Polymer& pol, const Interaction& interac)
{
	double tmp = 0;
	for(int bead=0; bead<pol.num_beads; ++bead)
		tmp += interac.ext_potential(pol[bead]);
	tmp /= pol.num_beads;
	return tmp;
}

double Observable::kinetic_energy(const Polymer& pol, const Interaction& interac)
{
	double tmp = 0;
	for(int bead=0; bead<pol.num_beads; ++bead)
		tmp += pol[bead].sqdist(pol[bead+1]);
	double offset = pol.num_beads*pol[0].size()/(2*beta);
	return offset - 0.5 * interac.get_spring_const() * tmp;
}

double Observable::total_energy(const Polymer& pol, const Interaction& interac)
{
	return potential_energy(pol, interac) + kinetic_energy(pol, interac);
}

double Observable::kinetic_energy_virial(const Polymer& pol, const Interaction& interac)
{
	double tmp = 0;
	for(int bead=0; bead<pol.num_beads; ++bead)
		tmp -= pol[bead]*interac.ext_force(pol[bead]); //minus to cancel minus in force expression
	return tmp/(2*pol.num_beads);
}

double Observable::potential_energy_cl(const Polymer& pol, const Interaction& interac)
{
	double tmp = 0;
	double spr_c = 0.5*interac.get_spring_const();
	for(int bead=0; bead<pol.num_beads; ++bead)
		tmp += interac.ext_potential(pol[bead])/pol.num_beads + spr_c*pol[bead].sqdist(pol[bead+1]);
	return tmp;
}

double Observable::kinetic_energy_cl(const Polymer& pol)
{
	double tmp = 0;
	for(int bead=0; bead<pol.num_beads; ++bead)
		tmp += pol.vels[bead]*pol.vels[bead];
	tmp *= 0.5*pol.mass;
	return tmp;
}

double Observable::total_energy_cl(const Polymer& pol, const Interaction& interac)
{
	return potential_energy_cl(pol, interac) + kinetic_energy_cl(pol); //note, take 1 away!
}



