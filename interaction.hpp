
#include "point.hpp"
#include "parameters.hpp"
#include "polymer.hpp"

#ifndef INTERACTION_HPP
#define INTERACTION_HPP

class Interaction{
public:
	Interaction(Parameters);

	double potential(const Point&);
	double ext_potential(const Point&);
	double spring_potential(const Point&, const Point&, const Point&);
	double bias_potential(const Point&);
//	double int_potential(Point);
	
//	Force force(const Point&,const Point&,const Point&);
	Force ext_force(const Point&);
	Force spring_force(const Point&, const Point&, const Point&) ;
	Force bias_force(const Point&);
//	Force int_force(Point);
	
	void update_forces(std::vector<Polymer>&);
	void update_one_pol_forces(Polymer&);
	
private:
	const double curvature; //=mw^2 for a harmonic oscillator
	const double spring_const;
	const double charge;
	const double diel_const;
};

#endif
