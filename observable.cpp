
#include <cmath>
#include "observable.hpp"


Observable::Observable()
{
	value,avg,avg_sq=0;
}

void Observable::measure(const std::vector<Polymer>& pols)
{
	
}


void Observable::set_zero()
{
	value = 0;
}

void Observable::update_avg(const int& num_samples)
{
	double tmp = value/num_samples;
	avg += tmp;
	avg_sq += tmp*tmp;
}

double Observable::get_avg() const
{
	return avg;
}

void Observable::normalize_avg(const int& num_blocks)
{
	avg /= num_blocks;
	avg_sq /= num_blocks;	
}

double Observable::std_dev(const int& n) const
{
	return std::sqrt((avg_sq-avg*avg)/(n-1));
}

void Observable::operator+=(double to_add)
{
	value += to_add;
}

