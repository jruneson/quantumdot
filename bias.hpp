
#include <vector>
#include <fstream>
#include "polymer.hpp"
#include "parameters.hpp"
#include "spline.hpp"


#ifndef BIAS_HPP
#define BIAS_HPP


class Bias{
public:
	Bias(const Parameters& params);
		
	void update_cv(const std::vector<Polymer>&, double, const double);
	double get_cv() const;

	void update_bias(const std::vector<Polymer>&,double,double);
	
	Force calc_force(const std::vector<Polymer>&, int, int) const;
	double get_rew_factor() const;
	const double get_gauss_width() const;
	const bool metad_on;
	std::vector<double> heights;
	std::vector<double> cv_centers;
	
private:
	double cv;
	double transient; //c(t)
	double rew_factor; // exp(beta(V(s)-c(t)))
	const int id;
	const double bias_factor;
	const double gauss_width;
	const double first_height;
	const int sign;
	const double exc_const;
	const double exponent_factor;
	
	Spline v_spline;
	Spline vder_spline;
	const double spline_step;
	
	double calc_bias(double) const; //calculate bias explicitly (only used when creating spline)
	double calc_bias_der(double) const;
	double coll_var(const std::vector<Polymer>&) const;
	double coll_var_der(const std::vector<Polymer>&) const;
	double scalar_product(const std::vector<Polymer>&, int) const;
	double gaussian(double, double, double) const;
	Force cv_grad(const std::vector<Polymer>&, int, int) const;
	double sum_exp(const std::vector<Polymer>&) const;
	Force two_terms(const std::vector<Polymer>&, int, int) const;
	void update_transient(double);
	
	std::ofstream cv_file;
	std::ofstream heights_file;
	std::ofstream cv_centers_file;
	std::ofstream rew_factor_file;

};

#endif
