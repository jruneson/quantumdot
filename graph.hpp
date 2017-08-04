#include <vector>
#include <tuple>
#include <cmath>
#include <iostream>
#include "parameters.hpp"
#include "polymer.hpp"
#include "point.hpp"

#ifndef GRAPH_HPP
#define GRAPH_HPP

class Graph{
public:
	Graph(const Parameters&, int graph_id);
	
	double get_weight(const std::vector<Polymer>&, const Graph&, bool) const;
	double get_weight_signed(const std::vector<Polymer>&, const Graph&) const;
	Force get_grad_weight(const std::vector<Polymer>&, const Graph&, bool, int, int) const;
	//double get_energy_diff(const std::vector<Polymer>&) const;
	//Force get_energy_diff_grad(const std::vector<Polymer>&, int,int) const;
	double energy_diff(const std::vector<Polymer>&, const Graph&) const;
	Force energy_diff_grad(const std::vector<Polymer>&, const Graph&, int, int) const;
	double energy_absolute(const std::vector<Polymer>&) const;
	
	int get_id() const;
	int get_mult() const;
	int get_sign() const;
	std::vector<std::vector<int>> get_chains() const;
	
	std::vector<std::pair<int,int>> get_exchange_pairs() const;
	
protected:
	std::vector<std::pair<int,int>> exchange_pairs; 	

private:
	std::vector<std::vector<int>> chains;
	const int id;
	const int num_parts;
	int exponent_sign;
	int mult; //multiplicity
	bool positive;
	const double exc_const;
	const double mass_factor;
	int exc_bead;
	int exc_bead2;
	const bool spin_proj; //True=>spin means "z projection", false>=spin means total spin
	
	int cyclic(int,int) const;
	double calc_exponent(const std::vector<Polymer>&, const Graph&) const;
	Force grad_exponent(const std::vector<Polymer>&, const Graph&, int, int) const;
	double weight(const std::vector<Polymer>&,const Graph&) const;
	double single_spring(const std::vector<Polymer>&,int,int) const;

};



#endif
