

#include "point.hpp"

//Point::Point(const Parameters& params) : point_coords(std::vector<double>(params.dim,0)){}
Point::Point(const Parameters& params)
{
	point_coords.assign(params.dim,0);
}

Point::Point(const int& dim)
{
	point_coords.assign(dim,0);
}


double Point::dist(const Point& p) const
{
	return std::sqrt(sqdist(p));
}

double Point::sqdist(const Point& p) const
{
	double tmp = 0;
	for(int d=0; d<point_coords.size(); ++d)
		tmp += std::pow((p.point_coords[d]-point_coords[d]),2);
	return tmp;
}

double Point::sqdist0() const
{
	double tmp=0;
	for(int d=0; d<point_coords.size(); ++d)
		tmp += std::pow(point_coords[d],2);
	return tmp;
}

double& Point::operator[](const int& index) 
{
	return point_coords[index];
}

const double& Point::operator[](const int& index) const
{
	return point_coords[index];
}


void Point::operator+=(const Point& p)
{
	for(int d=0; d<point_coords.size(); ++d)
		point_coords[d] += p.point_coords[d];
}

int Point::size() const
{
	return point_coords.size();
}



//outside class
Point operator+(const Point& p1, const Point& p2)
{
	Point p_new(p1.size());
	for(int d=0; d<p1.size(); ++d)
		p_new.point_coords[d]=p1.point_coords[d]+p2.point_coords[d];
	return p_new;
}

Point operator-(const Point& p1, const Point& p2)
{
	return p1 + (-1)*p2;
}

Point operator*(const Point& p, const double& c)
{
	Point p_new(p.size());
	for(int d=0; d<p.point_coords.size(); ++d)
		p_new.point_coords[d] = p.point_coords[d] * c;
	return p_new;
}

Point operator*(const double& c, const Point& p)
{
	return p*c;
}

Point operator/(const Point& p, const double& c)
{
	return p*(1.0/c);
}
