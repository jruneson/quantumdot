
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
				double, double, double);
	void print_measure(double, double, double);
	
	void set_zero();
	void update_avg(int,double);
	int get_id() const;
	double get_last_value() const;
	double get_value() const;
	double get_avg() const;
	double get_avg_sq() const;
	double get_weighted_avg() const;
	double get_weighted_avg_sq() const;
	void set_avgs(double,double,double,double,double);
	double std_dev(double,double) const;
	std::string get_name() const;
	
	void set_print_on();

	
private:
	double last_value;
	double value;
	double avg;
	double avg_sq;
	double weighted_avg;
	double weighted_avg_sq;
	std::ofstream file;
	const int id;
	double blocks;
	
	const double kin_offset;
	const double virial_offset;
	const double exc_const;
	const double exc_der_const;
	const int sign;

	
	double potential_energy(const std::vector<Polymer>&, const Interaction&) const;
	double kinetic_energy(const std::vector<Polymer>&, const Interaction&) const;
	double total_energy(const std::vector<Polymer>&, const Interaction&) const;
	double kinetic_energy_virial(const std::vector<Polymer>&, const Interaction&) const;
	double interparticle_energy(const std::vector<Polymer>&, const Interaction&) const;
	double potential_energy_cl(const std::vector<Polymer>&, const Interaction&) const;
	double kinetic_energy_cl(const std::vector<Polymer>&) const;
	double total_energy_cl(const std::vector<Polymer>&, const Interaction&) const;
	double spring_energy_cl(const std::vector<Polymer>&, const Interaction&) const;
	
	double exc_der(const std::vector<Polymer>&) const;
	double exc_der_virial(const std::vector<Polymer>&, const Interaction&, double) const;
	double scalar_product(const std::vector<Polymer>&, int) const;
	double scalar_product_conn(const std::vector<Polymer>&, int, int) const;
	
	bool print;
		
};

#endif
