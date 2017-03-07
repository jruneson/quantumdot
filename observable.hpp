
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include "polymer.hpp"
#include "interaction.hpp"


#ifndef OBSERVABLE_HPP
#define OBSERVABLE_HPP

class Observable {
public:
	//Observable();
	Observable(int);
	
	void measure(const std::vector<Polymer>&, Interaction&, const double&);
	void print_measure(const double&, const double&);
	
	void set_zero();
	void update_avg(const int&);
	void normalize_avg(const int&);
	double get_avg() const;
	double std_dev(const int&) const;
	std::string get_name();
	
	void set_print_on();
	
	//void operator+=(double);

private:
	double value;
	double avg;
	double avg_sq;
	std::ofstream file;
	const int id;
	
	double potential_energy(const Polymer&, const Interaction&);
	double potential_energy_cl(const Polymer&, const Interaction&);
	double kinetic_energy_cl(const Polymer&);
	double total_energy_cl(const Polymer&, const Interaction&);
	
	bool print;
		
};

/*class DistributionObservable : public Observable {
public:
	DistributionObservable(int,int,double);
	
	void measure(const std::vector<Polymer>&);
private:
	std::vector<double> histogram;
	
	int calc_bin();
	
};*/

#endif
