
#include "observable.hpp"


Observable::Observable(int _id, const Parameters& params) : id(_id),
				kin_offset(params.kin_offset),virial_offset(params.virial_offset),
				exc_const(params.exc_const), exc_der_const(params.exc_der_const)
{ 
	value = 0;
	avg = 0;
	avg_sq = 0;
	weighted_avg = 0;
	weighted_avg_sq = 0;
	blocks = 0;
	print = false;
}

void Observable::set_zero()
{
	value = 0;
}

void Observable::update_avg(int num_samples, double exc_avg)
{
	double tmp = value/num_samples;
	avg = (avg*blocks + tmp)/(blocks+1.0);
	avg_sq = (avg_sq*blocks + tmp*tmp)/(blocks+1.0);
	weighted_avg = (weighted_avg*blocks + tmp/exc_avg)/(blocks+1.0);
	weighted_avg_sq = (weighted_avg_sq*blocks + std::pow(tmp/exc_avg,2))/(blocks+1.0);
	++blocks;
}

int Observable::get_id() const
{
	return id;
}

double Observable::get_avg() const
{
	return avg;
}

double Observable::get_avg_sq() const
{
	return avg_sq;
}

double Observable::get_weighted_avg() const
{
	return weighted_avg;
	//return avg/exc_avg;
}

double Observable::get_weighted_avg_sq() const
{
	return weighted_avg_sq;
	//return avg_sq/(exc_avg*exc_avg);
}

void Observable::set_avgs(double avg_, double avg_sq_, double w_avg_, double w_avg_sq_, double blocks_)
{
	avg = avg_;
	avg_sq = avg_sq_;
	weighted_avg = w_avg_;
	weighted_avg_sq = w_avg_sq_;
	blocks = blocks_;
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
		case 23:
			return "Twopart_energy";
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

void Observable::measure(const std::vector<Polymer>& polymers, const Interaction& interac, 
						double time, double exc_factor, double rew_factor)
{
	double tmp = 0;
	switch(id)
	{
		case 0:
			tmp = potential_energy(polymers, interac);
			/*for(int m=0; m<n; ++m)
			{
				const Polymer& pol2 = polymers[m];
				tmp2 += potential_energy(pol,pol2,interac);
			}*/
			break;
		case 1:
			tmp = kinetic_energy(polymers,interac);
			break;
		case 2:
			tmp = total_energy(polymers,interac);
			break;
		case 3:
			tmp = kinetic_energy_virial(polymers,interac);
			break;
		case 10:
			tmp = polymers[0][0][0]; //x-coordinate of bead 0
			break;
		case 20:
			tmp = potential_energy_cl(polymers,interac); //note that exc_factor needs to be omitted for the last three
			break;
		case 21:
			tmp = kinetic_energy_cl(polymers);	
			break;
		case 22:
			tmp = total_energy_cl(polymers,interac);
			break;
		case 23:
			tmp = interparticle_energy(polymers,interac);
	}
	tmp *= exc_factor;
	switch(id)
	{
		case 1:
			tmp += exc_der(polymers);
			break;
		case 3:
			tmp += exc_der_virial(polymers);
			break;
		default:
			break;
	}
	tmp *= rew_factor;
	value += tmp;
	if(print)
		print_measure(tmp,time);
}

void Observable::print_measure(double measured_value, double time)
{
	file << time << "\t" << measured_value << std::endl;
}




double Observable::potential_energy(const std::vector<Polymer>& pols, const Interaction& interac)
{
	double tmp = 0;
	for(int n=0; n<pols.size(); ++n)
	{
		for(int m=0; m<n; ++m)
		{
			for(int bead=0; bead<pols[0].num_beads; ++bead)
			{
				tmp += interac.ext_potential(pols[n][bead]) + interac.ext_potential(pols[m][bead]);
				tmp += interac.two_particle_pot(pols[n][bead],pols[m][bead]);
			}
		}
	}
	tmp /= pols[0].num_beads;
	return tmp;
}

double Observable::kinetic_energy(const std::vector<Polymer>& pols, const Interaction& interac)
{
	double tmp = 0;
	for(const auto& pol : pols)
	{
		for(int bead=0; bead<pol.num_beads; ++bead)
			tmp += pol[bead].sqdist(pol[bead+1]);
	}
	return kin_offset - 0.5 * interac.get_spring_const() * tmp;
}

double Observable::total_energy(const std::vector<Polymer>& pols, const Interaction& interac)
{
	return potential_energy(pols, interac) + kinetic_energy_virial(pols, interac)
			+ interparticle_energy(pols, interac);
}

double Observable::kinetic_energy_virial(const std::vector<Polymer>& pols, const Interaction& interac)
{
	double tmp = 0;
	for(int n=0; n<pols.size(); ++n)
	{
		Point mean_point(pols[n][0].size());
		for(int bead=0; bead<pols[n].num_beads; ++bead)
			mean_point += pols[n][bead];
		mean_point *= 1.0/pols[n].num_beads;
		for(int bead=0; bead<pols[n].num_beads; ++bead)
		{
			const Point& p = pols[n][bead];
			tmp += (p-mean_point)*interac.ext_force(p);
			for(int m=0; m<pols.size(); ++m)
				if(m!=n)
					tmp += (p-mean_point)*interac.two_particle_force(p,pols[m][bead]);
		}
	}
	return virial_offset + (-1)*tmp/(2*pols[0].num_beads); //minus sign since force functions above calculate minus grad V
}

double Observable::interparticle_energy(const std::vector<Polymer>& pols, const Interaction& interac)
{
	double tmp=0;
	for(int bead=0; bead<pols[0].num_beads; ++bead)
		tmp += interac.two_particle_pot(pols[0][bead],pols[1][bead]);
	return tmp/pols[0].num_beads;
}

double Observable::potential_energy_cl(const std::vector<Polymer>& pols, const Interaction& interac)
{
	double tmp = 0;
	double spr_c = 0.5*interac.get_spring_const();
	for(const Polymer& pol : pols)
	{
		for(int bead=0; bead<pol.num_beads; ++bead)
			tmp += interac.ext_potential(pol[bead])/pol.num_beads + spr_c*pol[bead].sqdist(pol[bead+1]);
	}
	return tmp;
}

double Observable::kinetic_energy_cl(const std::vector<Polymer>& pols)
{
	double tmp = 0;
	for(const Polymer& pol : pols)
	{
		for(int bead=0; bead<pol.num_beads; ++bead)
			tmp += pol.vels[bead]*pol.vels[bead];
	}
	tmp *= 0.5*pols[0].mass * 1.747801e25;
	return tmp;
}

double Observable::total_energy_cl(const std::vector<Polymer>& pols, const Interaction& interac)
{
	return potential_energy_cl(pols, interac) + kinetic_energy_cl(pols) + interparticle_energy(pols,interac);
}

double Observable::exc_der(const std::vector<Polymer>& pols) const
{
	if(pols.size()==1)
		return 0;
	double tmp=0;
	for(int bead=0; bead<pols[0].num_beads; ++bead)
	{
		double sc_prod = scalar_product(pols,bead);
		tmp += sc_prod * std::exp(-exc_const*sc_prod);
	}
	return exc_der_const * tmp;
}

double Observable::exc_der_virial(const std::vector<Polymer>& pols) const
{
	if(pols.size()==1)
		return 0;
	Point tmp(pols[0][0].size());
	Point mean0(pols[0][0].size());
	Point mean1(pols[1][0].size());
	for(int bead=0; bead<pols[0].num_beads; ++bead)
	{
		mean0 += pols[0][bead];
		mean1 += pols[1][bead];
		double sc_prod = scalar_product(pols,bead);
		tmp += (pols[0][bead]-pols[1][bead])*std::exp(-exc_const*sc_prod);
	}
	mean0 /= pols[0].num_beads;
	mean1 /= pols[1].num_beads;
	return  exc_der_const*(mean0-mean1)*tmp;
}

double Observable::scalar_product(const std::vector<Polymer>& pols, int bead) const
{
	return (pols[0][bead]-pols[1][bead])*(pols[0][bead+1]-pols[1][bead+1]);
}

