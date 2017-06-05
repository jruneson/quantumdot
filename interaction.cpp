#include "interaction.hpp"

Interaction::Interaction(Parameters params) 
	: electrost_factor(params.electrost_factor),
	  curvature(params.dim), curvature_x(params.curvature_x),
	  curvature_y(params.curvature_y),
	  curvature_z(params.curvature_z),
	  spring_const(params.spring_const),
	  lj_energy24(24*params.lj_energy), 
	  lj_length_sq(params.lj_length*params.lj_length),
	  rcut2(2.5*2.5*params.lj_length*params.lj_length),
	  interaction_id(params.interaction_id),
	  dt_fast(params.dt_md), dt_slow(params.dt_md_slow)
{
	time_since_slow_update = 0;
	curvature[0] = curvature_x;
	int d = curvature.size();
	if(d>=2)
		curvature[1] = curvature_y;
	if(d>=3)
		curvature[2] = curvature_z;
	
	double tmp = std::pow(lj_length_sq/rcut2,3);
	vcut = lj_energy24/6.0 * tmp*(tmp-1);
}


double Interaction::ext_potential(const Point& p) const
{
	double tmp = 0;
	for(int d=0; d<p.size(); ++d)
		tmp += curvature[d]*std::pow(p[d],2);
	return 0.5*tmp;
	//return 0.5 * curvature * p.sqdist0(); //note that this must be changed if masses are to be unequal
}

double Interaction::two_particle_pot(const Point& p1, const Point& p2) const
{
	switch(interaction_id)
	{
		case 1: //Lennard-Jones
		{
			double sqdist = p1.sqdist(p2);
			if(sqdist>rcut2)
				return 0;
			double tmp = std::pow(lj_length_sq/sqdist,3);
			return lj_energy24/6.0*tmp*(tmp-1.0)-vcut;
		}
		case 2: //Electrostatic interaction
		{
			return electrost_factor/p1.dist(p2);
		}
		default:
			return 0;
	}
}

Force Interaction::ext_force(const Point& p) const
{
	Force f = Force(p.size());
	for(int d=0; d<f.size(); ++d)
		f[d] = -curvature[d]*p[d];
	return f;
}
Force Interaction::spring_force(const Point& last, const Point& p, const Point& next) 
{
	return spring_const*(last - 2*p + next);
}

Force Interaction::two_particle_force(const Point& p1, const Point& p2) const
{
	//On bead p1 by bead p2
	switch(interaction_id)
	{
		case 1: //Lennard-Jones
		{
			double sqdist = p1.sqdist(p2);
			if(sqdist>rcut2)
				break;
			double tmp = std::pow(lj_length_sq/sqdist,3); // (sigma/r)^6
			return (p1-p2)*(lj_energy24*(2*tmp*tmp-tmp)/sqdist);
		}
		case 2: //Electrostatic interaction
		{
			return (p1-p2)*(electrost_factor/std::pow(p1.dist(p2),3));
		}
		default:
			return Point(p1.size()); //No inter-particle interaction
	}
}


//Force int_force(Point p);

void Interaction::update_forces(std::vector<Polymer>& polymers, const Bias& bias)
{
	/*update_fast_forces(polymers);
	time_since_slow_update += dt_fast;
	if(time_since_slow_update >= dt_slow)
	{*/
	for(int n=0; n<polymers.size(); ++n)
	{
		Polymer& pol = polymers[n];
		for(int bead=0; bead<pol.num_beads; ++bead)
		{
			for(int m=0; m<n; ++m)
			{
				Force f = two_particle_force(pol[bead],polymers[m][bead])/polymers[0].num_beads;
				pol.twopart_forces[bead] = f;
				polymers[m].twopart_forces[bead] = f*(-1);
			}
			pol.ext_forces[bead] = ext_force(pol[bead])/pol.num_beads;
			pol.spring_forces[bead] = spring_force(pol[bead-1],pol[bead],pol[bead+1]);
			pol.bias_forces[bead] = bias.calc_force(polymers,bead,n);
		}

			
		if(polymers[0].connected && polymers.size()==2)
		{
			int num_beads = pol.num_beads;
			const Polymer& other_pol = polymers[polymers.size()-1-n];
			pol.spring_forces[0] = spring_force(other_pol[num_beads-1],pol[0],pol[1]);
			pol.spring_forces[num_beads-1] = spring_force(pol[num_beads-2],pol[num_beads-1],other_pol[0]);
			
			//pol.forces[0] += spring_force(other_pol[num_beads-1],pol[0],pol[1]);
			//pol.forces[num_beads-1] += spring_force(pol[num_beads-2],pol[num_beads-1],other_pol[0]);
			//std::cout << spring_force(other_pol[num_beads-1],pol[0],pol[1]).dist0() << "\t"
			//		  << spring_force(pol[num_beads-1],pol[0],pol[1]).dist0() << std::endl;
		}
		/*else
		{
			for(int bead=0; bead<pol.num_beads; ++bead)
				pol.forces[bead] += spring_force(pol[bead-1],pol[bead],pol[bead+1]);
		}*/
	}
	//}
}

/*
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
				}
			}
		}
	}
}*/


double Interaction::get_spring_const() const
{
	return spring_const;
}
