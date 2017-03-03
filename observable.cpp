
#include <cmath>
#include "observable.hpp"

/*
Observable::Observable()
{
	value,avg,avg_sq=0;
}*/

Observable::Observable(int id_to_set) : id(id_to_set)
{ value, avg, avg_sq = 0;}


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

void Observable::normalize_avg(const int& num_blocks)
{
	avg /= num_blocks;
	avg_sq /= num_blocks;	
}

double Observable::std_dev(const int& n) const
{
	return std::sqrt((avg_sq-avg*avg)/(n-1));
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
			return "Potential energy";
		case 20:
			return "Potential energy (cl)";
		case 21:
			return "Kinetic energy (cl)";
		case 22:
			return "Total energy (cl)";
		default:
			return "Unknown observable";
	}
}

double Observable::measure(const std::vector<Polymer>& polymers, Interaction& interac)
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
			case 20:
				tmp2 = potential_energy_cl(pol,interac);
			case 21:
				tmp2 = kinetic_energy_cl(pol);		
			case 22:
				tmp2 = total_energy_cl(pol,interac);
		}
		tmp += tmp2;
	}
	value += tmp;
	return tmp;
}

double Observable::potential_energy(const Polymer& pol, const Interaction& interac)
{
	return potential_energy_cl(pol,interac) / pol.num_beads;
}		

double Observable::potential_energy_cl(const Polymer& pol, const Interaction& interac)
{
	double tmp = 0;
	for(int bead=0; bead<pol.num_beads; ++bead)
		tmp += interac.ext_potential(pol[bead]);
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
	return potential_energy_cl(pol, interac) + kinetic_energy_cl(pol);
}
