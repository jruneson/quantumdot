
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include "polymer.hpp"
#include "interaction.hpp"
#include "point.hpp"
#include "parameters.hpp"
#include "bias.hpp"

#ifndef OBSERVABLE_HPP
#define OBSERVABLE_HPP

class Observable {
public:
	Observable(int,const Parameters&);
	
	void measure(const std::vector<Polymer>&, const Interaction&, double, 
				double, double);
	void print_measure(double, double);
	
	void set_zero();
	void update_avg(int);
	int get_id() const;
	double get_value() const;
	double get_avg() const;
	double get_avg_sq() const;
	double get_weighted_avg(double) const;
	double get_weighted_avg_sq(double) const;
	void set_avgs(double,double,double);
	double std_dev(double,double) const;
	std::string get_name() const;
	
	void set_print_on();

	
private:
	double value;
	double avg;
	double avg_sq;
	std::ofstream file;
	const int id;
	double blocks;
	
	const double kin_offset;
	const double virial_offset;
	const double exc_const;
	const double exc_der_const;

	
	double potential_energy(const Polymer&, const Interaction&);
	double kinetic_energy(const Polymer&, const Interaction&);
	double total_energy(const Polymer&, const Interaction&);
	double kinetic_energy_virial(const Polymer&, const Interaction&);
	double potential_energy_cl(const Polymer&, const Interaction&);
	double kinetic_energy_cl(const Polymer&);
	double total_energy_cl(const Polymer&, const Interaction&);
	
	double exc_der(const std::vector<Polymer>&) const;
	double exc_der_virial(const std::vector<Polymer>&) const;
	double scalar_product(const std::vector<Polymer>&, int) const;
	
	bool print;
		
};

#endif
