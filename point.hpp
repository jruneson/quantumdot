
#include <vector>
#include <cmath>

#include "parameters.hpp"

#ifndef POINT_HPP
#define POINT_HPP


class Point{ //storing components of coordinates, velocities, etc
public:
	Point(const Parameters&);
	Point(const int);

	int dim;
	std::vector<double> point_coords;
	
	double dist(const Point&) const; //distance to point
	double sqdist(const Point&) const; //squared distance to point
	double sqdist0() const; //squared distance to origin
	
	Point operator+(const Point&) const;
	void operator+=(const Point&);
	Point operator*(const double) const;
};


typedef Point Force; //using the same name for Forces (the functionality is the same, but the name Point for forces might be misleading)


#endif
