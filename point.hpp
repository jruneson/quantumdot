
#include <vector>
#include <cmath>
#include <iostream>

#include "parameters.hpp"

#ifndef POINT_HPP
#define POINT_HPP


class Point{ //storing components of coordinates, velocities, etc
public:
	explicit Point(const Parameters&);
	explicit Point(int);

	//int dim;
	std::vector<double> point_coords;
	
	double dist(const Point&) const; //distance to point
	double sqdist(const Point&) const; //squared distance to point
	double dist0() const;
	double sqdist0() const; //squared distance to origin
	int size() const;
	void set_zero();

	double& operator[](int);
	const double& operator[](int) const;
	void operator+=(const Point&);
	void operator-=(const Point&);
	void operator*=(double);
	void operator/=(double);
};

Point operator+(const Point&, const Point&);
Point operator-(const Point&, const Point&);
Point operator*(const Point&, double);
Point operator*(double, const Point&);
Point operator/(const Point&, double);
double operator*(const Point&, const Point&); //scalar product


typedef Point Force; //using the same class for Forces (the functionality is the same, but the name Point for forces might be misleading)


#endif
