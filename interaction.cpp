#include "interaction.hpp"

Interaction::Interaction(Parameters params) 
	: charge(params.charge), diel_const(params.diel_const){}


double Interaction::potential(Polymer pol) const
{
}

double Interaction::ext_potential(const Point& p) const
{
	return 0.5 * curvature * p.sqdist0();
}

double spring_potential(Point last, Point p, Point next);
double bias_potential(Point p);
//double int_potential(Point p);



Force Interaction::force(Point last,Point p,Point next)
{
	return p;
}


Force Interaction::force(const Point& last, const Point& p, const Point& next)
{
	return ext_force(p) + spring_force(last,p,next) + bias_force(p);
}

Force Interaction::ext_force(Point p)
{
	return -curvature * p;
}
Force spring_force(Point last, Point p, Point next);
Force bias_force(Point p);
//Force int_force(Point p);

void update_forces(Polymer& pol)
