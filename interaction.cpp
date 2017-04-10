#include "interaction.hpp"

Interaction::Interaction(Parameters params) 
	: charge(params.charge), diel_const(params.diel_const),
	  curvature(params.curvature), spring_const(params.spring_const),
	  lj_energy24(24*params.lj_energy), 
	  lj_length_sq(params.lj_length*params.lj_length),
	  rcut2(2.5*2.5*params.lj_length*params.lj_length),
	  interaction_id(params.interaction_id),
	  dt_fast(params.dt_md), dt_slow(params.dt_md_slow)
{
	time_since_slow_update = 0;
	
	double tmp = std::pow(lj_length_sq/rcut2,3);
	vcut = lj_energy24/6.0 * tmp*(tmp-1);
}

/*
double Interaction::potential(const Point& p)
{
	return 0;
}*/

double Interaction::ext_potential(const Point& p) const
{
	return 0.5 * curvature * p.sqdist0(); //note that this must be changed if masses are to be unequal
}

double Interaction::two_particle_pot(const Point& p1, const Point& p2) const
{
	switch(interaction_id)
	{
		case 1:
		{
			double sqdist = p1.sqdist(p2);
			if(sqdist>rcut2)
				return 0;
			double tmp = std::pow(lj_length_sq/sqdist,3);
			return lj_energy24/6.0*tmp*(tmp-1.0)-vcut;
		}
		default:
			return 0;
	}
}

/*
double Interaction::spring_potential(const Point& last, const Point& p, const Point& next) 
{
	return 0;
}
double Interaction::bias_potential(const Point& p) 
{
	return 0;
}*/

//double int_potential(Point p);

/*
void Interaction::update_one_pol_forces(Polymer& pol)
{
	for(int bead=0; bead<pol.num_beads; ++bead)
	{
		pol.forces[bead] = ext_force(pol[bead])/pol.num_beads
						+ spring_force(pol[bead-1],pol[bead],pol[bead+1])
						+ bias.calc_force(pol[bead]);
	}
}*/

/*Force Interaction::force(const Point& last, const Point& p, const Point& next)
{
	return ext_force(p) + spring_force(last,p,next) + bias_force(p);
}*/

Force Interaction::ext_force(const Point& p) const
{
	return (-curvature) * p;
}
Force Interaction::spring_force(const Point& last, const Point& p, const Point& next) 
{
	return spring_const*(last - 2*p + next);
}

Force Interaction::two_particle_force(const Point& p1, const Point& p2) const
{
	switch(interaction_id)
	{
		case 1: //Lennard-Jones
		{
			double sqdist = p1.sqdist(p2);
			if(sqdist>rcut2)
				break;
			double tmp = std::pow(lj_length_sq/sqdist,3); // (sigma/r)^6
			//std::cout << (p1-p2).dist0() << "\t" << std::sqrt(sqdist) << std::endl;
			//if(std::sqrt(sqdist)<10000)
			//	std::cout << std::sqrt(sqdist) << "\t" << spring_const << "\t" << (lj_energy24*(2*tmp*tmp-tmp)/sqdist) << std::endl;
			return (p1-p2)*(lj_energy24*(2*tmp*tmp-tmp)/sqdist);
		}
	}
	return Point(p1.size()); //No inter-particle interaction
}


//Force int_force(Point p);

void Interaction::update_forces(std::vector<Polymer>& polymers, const Bias& bias)
{
	update_fast_forces(polymers);
	time_since_slow_update += dt_fast;
	if(time_since_slow_update >= dt_slow)
	{
		for(int n=0; n<polymers.size(); ++n)
		{
			Polymer& pol = polymers[n];
			for(int bead=0; bead<pol.num_beads; ++bead)
			{
				pol.forces[bead] = ext_force(pol[bead])/pol.num_beads
					+ spring_force(pol[bead-1],pol[bead],pol[bead+1])
					+ bias.calc_force(polymers,bead,n);
				//if(pol.fast_forces[bead].sqdist0() > 0.001)
					//std::cout << pol.forces[bead].sqdist0() << "\t" << pol.fast_forces[bead].sqdist0() << std::endl;
			}
		}
		time_since_slow_update = 0;
		/*int bead=0;
		const Polymer& pol=polymers[0];
		double dist = polymers[0][bead].dist(polymers[1][bead]);
		if(dist < 27)
		{
			//double sqdist = dist*dist;
			//double tmp = std::pow(lj_length_sq/sqdist,3); // (sigma/r)^6
			std::cout << dist << "\t" //<< tmp << "\t" << ((polymers[0][bead]-polymers[1][bead])*(lj_energy24*(2*tmp*tmp-tmp)/sqdist)).dist0() << "\t"
					  << (ext_force(pol[bead])/pol.num_beads).dist0() << "\t"
					  << (spring_force(pol[bead-1],pol[bead],pol[bead+1])).dist0() << "\t"
					  << (two_particle_force(polymers[0][bead],polymers[1][bead])).dist0()
					  << std::endl;
			std::cout << dist << "\t" << ext_potential(pol[bead])/pol.num_beads << "\t"
					  << spring_const*(pol[bead]-pol[bead+1])*(pol[bead]-pol[bead+1]) << "\t"
					  << two_particle_pot(polymers[0][bead],polymers[1][bead]) << std::endl;
		}*/
	}
	//if(polymers[0].fast_forces[0].dist0() > 1)
	//	std::cout << polymers[0].forces[0].dist0() << "\t" << polymers[0].fast_forces[0].dist0() << std::endl;
}

void Interaction::update_fast_forces(std::vector<Polymer>& polymers)
{
	if(interaction_id!=0)
	{
		for(int n=0; n<polymers.size(); ++n)
		{
			for(int m=0; m<n; ++m)
			{
				for(int bead=0; bead<polymers[0].num_beads; ++bead)
				{
					Force f = two_particle_force(polymers[n][bead],polymers[m][bead])/(polymers[0].num_beads);
					polymers[n].fast_forces[bead] = f;
					polymers[m].fast_forces[bead] = f*(-1);
					/*if(f.dist0()>3)
						std::cout << f.dist0() << std::endl;*/
				}
			}
		}
	}
}

double Interaction::get_spring_const() const
{
	return spring_const;
}
