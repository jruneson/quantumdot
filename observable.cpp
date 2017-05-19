
#include "observable.hpp"


Observable::Observable(int _id, const Parameters& params) : id(_id),
				kin_offset(params.kin_offset),virial_offset(params.virial_offset),
				exc_const(params.exc_const), exc_der_const(params.exc_der_const),
				sign(params.sign)
{ 
	last_value = 0;
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

double Observable::get_last_value() const
{
	return last_value;
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
		case 24:
			return "Spring_pot_en";
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
			break;
		case 24:
			tmp = spring_energy_cl(polymers,interac);
			break;
	}
	double bare_value = tmp;
	tmp *= exc_factor;
	switch(id)
	{
		case 1:
			tmp += exc_der(polymers);
			break;
		case 2:
		case 3:
			tmp += exc_der_virial(polymers,interac);
			break;
		default:
			break;
	}
	tmp *= rew_factor;
	last_value = tmp;
	value += tmp;
	if(print)
		print_measure(time,tmp,bare_value);
}

void Observable::print_measure(double time, double weighted_value, double bare_value)
{
	file << time << "\t" << weighted_value << "\t" << bare_value << std::endl;
}




double Observable::potential_energy(const std::vector<Polymer>& pols, const Interaction& interac) const
{
	double tmp = 0;
	for(int n=0; n<pols.size(); ++n)
	{
		for(int bead=0; bead<pols[0].num_beads; ++bead)
		{
			tmp += interac.ext_potential(pols[n][bead]);
			for(int m=0; m<n; ++m)
			{
				tmp += interac.two_particle_pot(pols[n][bead],pols[m][bead]);
			}
		}
	}
	tmp /= pols[0].num_beads;
	return tmp;
}

double Observable::kinetic_energy(const std::vector<Polymer>& pols, const Interaction& interac) const
{
	double tmp = 0;
	for(int n=0; n<pols.size(); ++n)
	//for(const auto& pol : pols)
	{
		const auto& pol = pols[n];
		if(pol.connected)
		{
			for(int bead=0; bead<pol.num_beads-1; ++bead)
				tmp += pol[bead].sqdist(pol[bead+1]);
			const auto& other_pol = pols[pols.size()-1-n];
			tmp += pol[pol.num_beads-1].sqdist(other_pol[0]);
		}
		else
			for(int bead=0; bead<pol.num_beads; ++bead)
				tmp += pol[bead].sqdist(pol[bead+1]);
	}
	return kin_offset - 0.5 * interac.get_spring_const() * tmp;
}

double Observable::total_energy(const std::vector<Polymer>& pols, const Interaction& interac) const
{
	return potential_energy(pols, interac) + kinetic_energy_virial(pols, interac);
}

double Observable::kinetic_energy_virial(const std::vector<Polymer>& pols, const Interaction& interac) const
{
	double tmp = 0;
	/*if(pols[0].connected)
	{
		Point mean_point(pols[0][0]);
		for(int n=0; n<pols.size(); ++n)
			mean_point += pols[n].mean();
		mean_point /= pols.size();
		for(int n=0; n<pols.size(); ++n)
		{
			for(int bead=0; bead<pols[0].num_beads; bead++)
			{
				const Point& p = pols[n][bead];
				tmp += (p-mean_point)*interac.ext_force(p);
				for(int m=0; m<pols.size(); ++m)
					if(m!=n)
						tmp += (p-mean_point)*interac.two_particle_force(p,pols[m][bead]);
			}
		}
	}
	else
	{*/
		for(int n=0; n<pols.size(); ++n)
		{
			const Point& mean_point = pols[n].mean();
			for(int bead=0; bead<pols[n].num_beads; ++bead)
			{
				const Point& p = pols[n][bead];
				tmp += (p-mean_point)*interac.ext_force(p);
				for(int m=0; m<pols.size(); ++m)
					if(m!=n)
						tmp += (p-mean_point)*interac.two_particle_force(p,pols[m][bead]); //double-counting?
			}
		}
	//}
	//if(pols[0].connected)
	//	return virial_offset + (-1)*tmp/(2*pols[0].num_beads) + exc_der_const*scalar_product(pols,pols[0].num_beads-1);
	return virial_offset + (-1)*tmp/(2*pols[0].num_beads); //minus sign since force functions above calculate minus grad V
}

	
double Observable::interparticle_energy(const std::vector<Polymer>& pols, const Interaction& interac) const
{
	double tmp=0;
	for(int bead=0; bead<pols[0].num_beads; ++bead)
		tmp += interac.two_particle_pot(pols[0][bead],pols[1][bead]);
	return tmp/pols[0].num_beads;
}

double Observable::potential_energy_cl(const std::vector<Polymer>& pols, const Interaction& interac) const
{
	double tmp = 0;
	for(const Polymer& pol : pols)
	{
		for(int bead=0; bead<pol.num_beads; ++bead)
			tmp += interac.ext_potential(pol[bead])/pol.num_beads;
	}
	tmp += spring_energy_cl(pols, interac);
	return tmp;
}

double Observable::spring_energy_cl(const std::vector<Polymer>& pols, const Interaction& interac) const
{
	double tmp=0;
	double spr_c = 0.5*interac.get_spring_const();
	for(const Polymer& pol : pols)
	{
		for(int bead=0; bead<pol.num_beads; ++bead)
			tmp += spr_c*pol[bead].sqdist(pol[bead+1]);
	}
	return tmp;
}

double Observable::kinetic_energy_cl(const std::vector<Polymer>& pols) const
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

double Observable::total_energy_cl(const std::vector<Polymer>& pols, const Interaction& interac) const
{
	return potential_energy_cl(pols, interac) + kinetic_energy_cl(pols) + interparticle_energy(pols,interac);
}

double Observable::exc_der(const std::vector<Polymer>& pols) const
{
	if(exc_der_const==0)
		return 0;
	double tmp=0;
	double sc_prod = scalar_product(pols,pols[0].num_beads-1);
	if(pols[0].connected)
	{
		/*for(int bead=0; bead<pols[0].num_beads; ++bead)
		{
			double sc_prod = scalar_product(pols,bead);
			tmp += (-1)*sc_prod*std::exp(exc_const*sc_prod);
		}*/

		tmp = sc_prod*std::exp(exc_const*sc_prod);
	}
	else
	{
		/*for(int bead=0; bead<pols[0].num_beads; ++bead)
		{
			double sc_prod = scalar_product(pols,bead);
			tmp += sc_prod * std::exp(-exc_const*sc_prod);
		}*/
		tmp = (-1)*sc_prod*std::exp(-exc_const*sc_prod);
	}
	return exc_der_const * tmp;
}

double Observable::exc_der_virial(const std::vector<Polymer>& pols, const Interaction& interac) const
{
	if(exc_der_const==0)
		if(!pols[0].connected)
			return virial_offset;
	double tmp=0;
	Point total_mean(pols[0][0].size());
	for(int n=0; n<pols.size(); ++n)
		total_mean += pols[n].mean();
	total_mean /= pols.size();
	int num_beads = pols[0].num_beads;
	//double sc_prod = scalar_product(pols,num_beads-1);
	//const Point& tmp = (pols[0][num_beads-1]-pols[1][num_beads-1]+pols[0][0]-pols[1][0]);
	Force pot_der(pols[0][0].size());
	for(int n=0; n<pols.size(); ++n)
	{
		pot_der.set_zero();
		const Point& mean = pols[n].mean();
		for(int bead=0; bead<pols[n].num_beads; ++bead)
		{
			const Point& p = pols[n][bead];
			pot_der += interac.ext_force(p);
			for(int m=0; m<pols.size(); ++m)
				if(m!=n)
					pot_der += interac.two_particle_force(p,pols[m][bead]);
		}
		tmp -= pot_der*(pols[n].mean()-total_mean); //minus to compensate minus sign of force
	}
	tmp /= (2*num_beads);
	if(pols[0].connected)
		return virial_offset*std::exp(exc_const*scalar_product(pols,num_beads-1))+sign*tmp;
	else
		return virial_offset + sign*tmp*std::exp(-exc_const*scalar_product(pols,num_beads-1));
	//if(pols[0].connected)
	//	return 0.5*std::abs(exc_der_const)*(mean0-mean1)*tmp;
	//else
	//	tmp *= std::exp(-exc_const*sc_prod);
	/*{
		for(int bead=0; bead<pols[0].num_beads; ++bead)
		{
			mean0 += pols[0][bead];
			mean1 += pols[1][bead];
			double sc_prod = scalar_product(pols,bead);
			tmp += (pols[0][bead]-pols[1][bead])*std::exp(-exc_const*sc_prod);
		}
		mean0 /= pols[0].num_beads;
		mean1 /= pols[1].num_beads;
	}*/
	//return (-1)*0.5*exc_der_const*std::exp(-exc_const*sc_prod)*(mean0-mean1)*tmp;
}

double Observable::scalar_product(const std::vector<Polymer>& pols, int bead) const
{
	return (pols[0][bead]-pols[1][bead])*(pols[0][bead+1]-pols[1][bead+1]);
}

/*double Observable::scalar_product_conn(const std::vector<Polymer>& pols, int bead, int num_beads) const
{
	return (pols[0][bead]-pols[0][bead+num_beads])*(pols[0][bead+1]-pols[0][bead+num_beads+1]);
}*/

