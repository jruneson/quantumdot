
#include "point.hpp"
#include "parameters.hpp"

#ifndef INTERACTION_HPP
#define INTERACTION_HPP

class Interaction{
public:
	Interaction(Parameters);

	double potential(Point) const;
	double ext_potential(Point);
	double spring_potential(Point, Point, Point);
	double bias_potential(Point);
//	double int_potential(Point);
	
	Force force(Point,Point,Point);
	Force ext_force(Point);
	Force spring_force(Point, Point, Point);
	Force bias_force(Point);
//	Force int_force(Point);
	
	void update_forces(Polymer& pol)
	
private:
	const double curvature; //=mw^2 for a harmonic oscillator
	const double charge;
	const double diel_const;
};

#endif
