
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include "polymer.hpp"
#include "interaction.hpp"
#include "point.hpp"



#ifndef OBSERVABLE_HPP
#define OBSERVABLE_HPP

class Observable {
public:
	Observable(int,double);
	
	void measure(const std::vector<Polymer>&, Interaction&, const double&);
	void print_measure(const double&, const double&);
	
	void set_zero();
	void update_avg(const int&);
	double get_value() const;
	double get_avg() const;
	double std_dev() const;
	std::string get_name() const;
	
	void set_print_on();
	
private:
	double value;
	double avg;
	double avg_sq;
	std::ofstream file;
	const int id;
	const double beta;
	int blocks;
	
	double potential_energy(const Polymer&, const Interaction&);
	double kinetic_energy(const Polymer&, const Interaction&);
	double total_energy(const Polymer&, const Interaction&);
	double kinetic_energy_virial(const Polymer&, const Interaction&);
	double potential_energy_cl(const Polymer&, const Interaction&);
	double kinetic_energy_cl(const Polymer&);
	double total_energy_cl(const Polymer&, const Interaction&);
	
	bool print;
		
};

#endif
