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
	
	double get_weight(const std::vector<Polymer>&, bool) const;
	double get_weight_signed(const std::vector<Polymer>&) const;
	Force get_grad_weight(const std::vector<Polymer>&, bool, int, int) const;
	double get_energy_diff(const std::vector<Polymer>&) const;
	Force get_energy_diff_grad(const std::vector<Polymer>&,int,int) const;
	
	int get_id() const;
	int get_mult() const;
	int get_sign() const;
	std::vector<std::vector<int>> get_chains() const;

private:
	std::vector<std::vector<std::pair<int,int>>> exchange_pairs; 
	std::vector<std::vector<int>> chains;
	const int id;
	int exponent_sign;
	int mult; //multiplicity
	bool positive;
	const double exc_const;
	int exc_bead;
	const bool spin_proj; //True=>spin means "z projection", false>=spin means total spin
	
	int cyclic(int,int) const;
	double calc_exponent(const std::vector<Polymer>&, const std::vector<std::pair<int,int>>&) const;
	Force grad_exponent(const std::vector<Polymer>&, const std::vector<std::pair<int,int>>&, int, int) const;


};



#endif
