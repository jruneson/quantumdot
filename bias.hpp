
#include <vector>

#ifndef BIAS_HPP
#define BIAS_HPP


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

#endif
