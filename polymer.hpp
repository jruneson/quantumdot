#include <vector>
#include <iostream>

#include "parameters.hpp"
#include "point.hpp"

#ifndef POLYMER_HPP
#define POLYMER_HPP

class Polymer{
public:
	//Polymer();
	explicit Polymer(const Parameters&);
	
	const Point& operator[](int) const; //use this operator to reach coordinate Points
	Point& operator[](int);
	
	void move();
	void update_vels();
	Point mean() const; //centroid position
	Force get_potential_force(int) const;
	
	int num_beads;
	std::vector<Point> coords;
	std::vector<Point> vels;
	std::vector<Force> ext_forces;
	std::vector<Force> twopart_forces;
	std::vector<Force> spring_forces;
	std::vector<Force> bias_forces;
	const double dt_2m;
	const double dt;
	const double mass;
	bool connected;
};

#endif
