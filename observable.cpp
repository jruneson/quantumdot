
#include "observable.hpp"


Observable::Observable(int _id, const Parameters& params) : id(_id),
				kin_offset(params.kin_offset),virial_offset(params.virial_offset),
				exc_const(params.exc_const), exc_der_const(params.exc_der_const)
{ 
	value = 0;
	avg = 0;
	avg_sq = 0;
	blocks = 0;
	print = false;
}

void Observable::set_zero()
{
	value = 0;
}

void Observable::update_avg(int num_samples)
{
	double tmp = value/num_samples;
	avg = (avg*blocks + tmp)/(blocks+1);
	avg_sq = (avg_sq*blocks + tmp*tmp)/(blocks+1);
	++blocks;
}

double Observable::get_avg(double exc_avg) const
{
	return avg/exc_avg;
}

double Observable::get_value() const
{
	return value;
}
/*
void Observable::normalize_avg(const int& num_blocks)
{
//	std::cout << avg << "\t" << num_blocks;
	avg /= num_blocks;
	avg_sq /= num_blocks;
//	std::cout << "\t" << avg << std::endl;
}*/

double Observable::std_dev(double exc_avg, double exc_avg_sq) const
{
	if(blocks>=2)
	{
		double num_error = std::sqrt((avg_sq - avg*avg)*blocks/(blocks-1.0)); //note that avg_sq is not normalized
		double den_error = std::sqrt((exc_avg_sq - exc_avg*exc_avg)*blocks/(blocks-1.0));
		return avg/exc_avg * (num_error/avg + den_error/exc_avg);		
		//error = num/den * (num_err/num + den_err/den)
	}
	return 0;
}


std::string Observable::get_name() const
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

void Observable::measure(const std::vector<Polymer>& polymers, Interaction& interac, 
						double time, double exc_factor)
{
	double tmp = 0;
	for(int n=0; n<polymers.size(); ++n)
	{
		const Polymer& pol = polymers[n];
		double tmp2 = 0;
		switch(id)
		{
			case 0:
				tmp2 = potential_energy(pol,interac) * exc_factor;
				break;
			case 1:
				tmp2 = kinetic_energy(pol,interac) * exc_factor;
				break;
			case 2:
				tmp2 = total_energy(pol,interac) * exc_factor;
				break;
			case 3:
				tmp2 = kinetic_energy_virial(pol,interac) * exc_factor;
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
	//tmp *= exc_factor; //consider moving this out of the switch-case
	switch(id)
	{
		case 1:
			tmp += exc_der(polymers);
			break;
		case 2:
			tmp += exc_der(polymers);
		default:
			break;
	}
	value += tmp;
	if(print)
		print_measure(tmp,time);
}

void Observable::print_measure(double measured_value, double time)
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
	return kin_offset - 0.5 * interac.get_spring_const() * tmp;
}

double Observable::total_energy(const Polymer& pol, const Interaction& interac)
{
	return potential_energy(pol, interac) + kinetic_energy(pol, interac);
}

double Observable::kinetic_energy_virial(const Polymer& pol, const Interaction& interac)
{
	double tmp = 0;
	Point mean_point(pol[0].size());
	for(int bead=0; bead<pol.num_beads; ++bead)
		mean_point += pol[bead];
	mean_point *= 1.0/pol.num_beads;
	for(int bead=0; bead<pol.num_beads; ++bead)
		tmp -= (pol[bead]-mean_point)*interac.ext_force(pol[bead]); //minus to cancel minus in force expression
	return virial_offset + tmp/(2*pol.num_beads);
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
	return potential_energy_cl(pol, interac) + kinetic_energy_cl(pol);
}

double Observable::exc_der(const std::vector<Polymer>& pols)
{
	double tmp=0;
	for(int bead=0; bead<pols[0].num_beads; ++bead)
	{
		double scalar_prod = (pols[0][bead]-pols[1][bead])*(pols[0][bead+1]-pols[1][bead+1]);
		tmp += scalar_prod * std::exp(-exc_const*scalar_prod);
	}
	return exc_der_const * tmp;
}


