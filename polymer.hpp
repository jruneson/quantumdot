#include <vector>
#include <iostream>

#include "parameters.hpp"
#include "point.hpp"

#ifndef POLYMER_HPP
#define POLYMER_HPP

class Polymer{
public:
	//Polymer();
	Polymer(const Parameters&);
	
	const Point& operator[](int) const; //use this operator to reach coordinate Points
	Point& operator[](int);
	
	void move();
	void update_vels();
	
	int num_beads;
	std::vector<Point> coords;
	std::vector<Point> vels;
	std::vector<Force> forces;
	const double dt_2m;
	const double dt;
	const double mass;
};

#endif
