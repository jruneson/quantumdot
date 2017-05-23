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
	Force get_grad_weight(const std::vector<Polymer>&, bool, int, int) const;

	int get_id() const;
	int get_mult() const;
	int get_sign() const;

private:
	std::vector<std::vector<std::pair<int,int>>> exchange_pairs; 
	const int id;
	int mult; //multiplicity
	bool positive;
	const double exc_const;
	int exc_bead;
	
	int cyclic(int,int) const;

};



#endif
