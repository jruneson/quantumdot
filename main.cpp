#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <array>

#include "GLE-OPT.h"
#include "spline.hpp"

struct Point{ //storing components of coordinates, velocities, etc
	std::array<double,int dim> point_coords;
	double dist(const Point &p) const; //distance to point
	operator+(const Point &p);
	operator+=(const Point &p);
	operator*(const double c);
};

typedef Point Force;

template <>
Point const std::array<Point, int>::operator[]() //make specialization for points
{
	//make cyclic indices
}


class Polymer{
public:
	Polymer(int,int);
	Point const operator[](); //use this operator to reach coordinate Points
	void move();

private:
	std::array<Point,num_beads> coords;
	std::array<Point,num_beads> vels;
	std::array<Point,num_beads> forces;
	const double mass; //if we want it to be different between polymers
};

class Bias{
public:
	double calcBias(const double) const; //calculate bias explicitly (only used when creating spline)
	double calcBiasDer(const double) const;
	
private:
	std::vector<double> heights;
	std::vector<double> cv_centers;
	double bias_factor;
	double gauss_width;
	int sign;
};

class Interaction{
public:
	double potential(Point p) const;
	double ext_potential(Point p);
	double spring_potential(Point last, Point p, Point next);
	double bias_potential(Point p);
	double int_potential(Point p);
	
	Force force(Point p);
	Force ext_force(Point p);
	Force spring_force(Point last, Point p, Point next);
	Force bias_force(Point p);
	Force int_force(Point p);
	
private:
	double charge;
};


class Simulation{
public:
	void setup();
	void run();
	void run_block();
	void verlet_step();
	void print_to_file();
	void stop();
	

	
	double coll_var() const;
	Force grad_coll_var() const;

private:
	std::array<Polymer,N> polymers;
	Bias bias;
	Spline v_spline;
	Spline f_spline;
	const double dt;
	const double dt_2m;
	const Interaction i;
}

Simulation::verlet_step()
{
	for(Polymer p : polymers)
		for(int bead=0; bead<num_beads; ++bead)
		{	
			p.vels[bead] += force(p[bead-1], p[bead], p[bead+1]) * dt_2m;
			p[bead] += mass*p.vels[bead] * dt;
			p.vels[bead] += force(p[bead-1], p[bead], p[bead+1]) * dt_2m;
		}
	}
}

Simulation::force(Point last, Point p, Point next)
{
	return ext_force(p) + spring_force(last,p,next) + bias_force(p) + int_force(p);
}


int main()
{
	Simulation sim;
	sim.setup();
	sim.run();
	return 0;
}
