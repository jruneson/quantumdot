

#include "point.hpp"

Point::Point(const Parameters& params) : dim(params.dim) {}

Point::Point(const int dims) : dim(dims) {}




double Point::dist(const Point& p) const
{
	return std::sqrt(sqdist(p));
}

double Point::sqdist(const Point& p) const
{
	double tmp = 0;
	for(int d=0; d<dim; ++d)
		tmp += std::pow((p.point_coords[d]-point_coords[d]),2);
	return tmp;
}

double Point::sqdist0() const
{
	double tmp=0;
	for(int d=0; d<dim; ++d)
		tmp += std::pow(point_coords[d],2);
	return tmp;
}

Point Point::operator+(const Point& p) const
{
	Point p_new(dim);
	for(int d=0; d<dim; ++d)
		p_new.point_coords[d]=point_coords[d]+p.point_coords[d];
	return p_new;
}

void Point::operator+=(const Point& p)
{
	for(int d=0; d<dim; ++d)
		point_coords[d] += p.point_coords[d];
}

Point Point::operator*(const double c) const
{
	Point p_new(dim);
	for(int d=0; d<dim; ++d)
		p_new.point_coords[d] = point_coords[d] * c;
	return p_new;
}

/*Point Point::operator=(const Point& p) const
{
	for
	
}*/
