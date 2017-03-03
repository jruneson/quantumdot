
#include <string>
#include <vector>
#include <iostream>
#include "polymer.hpp"
#include "interaction.hpp"


#ifndef OBSERVABLE_HPP
#define OBSERVABLE_HPP

class Observable {
public:
	//Observable();
	Observable(int);
	
	double measure(const std::vector<Polymer>&, Interaction&);
	
	void set_zero();
	void update_avg(const int&);
	void normalize_avg(const int&);
	double get_avg() const;
	double std_dev(const int&) const;
	std::string get_name();
	
	void operator+=(double);

private:
	double value;
	double avg;
	double avg_sq;
	const int id;
	
	double potential_energy(const Polymer&, const Interaction&);
	double potential_energy_cl(const Polymer&, const Interaction&);
	double kinetic_energy_cl(const Polymer&);
	double total_energy_cl(const Polymer&, const Interaction&);
		
};

#endif
