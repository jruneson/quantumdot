#include "interaction.hpp"

Interaction::Interaction(Parameters params) 
	: charge(params.charge), diel_const(params.diel_const),
	  curvature(params.curvature), spring_const(params.spring_const) {}

/*
double Interaction::potential(const Point& p)
{
	return 0;
}*/

double Interaction::ext_potential(const Point& p) const
{
	return 0.5 * curvature * p.sqdist0(); //note that this must be changed if masses are to be unequal
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


void Interaction::update_one_pol_forces(Polymer& pol)
{
	for(int bead=0; bead<pol.num_beads; ++bead)
	{
		pol.forces[bead] = ext_force(pol[bead])/pol.num_beads
						+ spring_force(pol[bead-1],pol[bead],pol[bead+1])
						+ bias_force(pol[bead]);
	}
}

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

Force Interaction::bias_force(const Point& p) 
{
	return Force(p.size());
}

//Force int_force(Point p);

void Interaction::update_forces(std::vector<Polymer>& polymers)
{
	for(int n=0; n<polymers.size(); ++n)
	{
		Polymer& pol = polymers[n];
		update_one_pol_forces(pol);
		//two_pol_forces
	}
}

double Interaction::get_spring_const() const
{
	return spring_const;
}
