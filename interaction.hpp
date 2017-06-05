
#include "point.hpp"
#include "parameters.hpp"
#include "polymer.hpp"
#include "bias.hpp"

#ifndef INTERACTION_HPP
#define INTERACTION_HPP

class Interaction{
public:
	explicit Interaction(Parameters);

//	double potential(const Point&) const;
	double ext_potential(const Point&) const;
	double two_particle_pot(const Point&, const Point&) const;
//	double spring_potential(const Point&, const Point&, const Point&);
//	double bias_potential(const Point&);
//	double int_potential(Point);
	
//	Force force(const Point&,const Point&,const Point&);
	Force ext_force(const Point&) const;
	Force spring_force(const Point&, const Point&, const Point&) ;
//	Force bias_force(const Point&);
	Force two_particle_force(const Point&, const Point&) const;
	
	void update_forces(std::vector<Polymer>&, const Bias&);
	void update_fast_forces(std::vector<Polymer>&);
//	void update_one_pol_forces(Polymer&);
	
	double get_spring_const() const;
	
private:
	Point curvature; //=mw^2 for a harmonic oscillator
	const double curvature_x;
	const double curvature_y;
	const double curvature_z;
	const double spring_const;
	const double electrost_factor;
	//const double charge;
	//const double diel_const;
	
	const int interaction_id;
	const double lj_energy24;
	const double lj_length_sq;
	const double rcut2;
	double vcut;
	
	const double dt_fast;
	const double dt_slow;
	double time_since_slow_update;
};

#endif
