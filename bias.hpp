
#include <vector>
#include <fstream>
#include "polymer.hpp"
#include "parameters.hpp"
#include "spline.hpp"
#include "graph.hpp"

#ifndef BIAS_HPP
#define BIAS_HPP


class Bias{
public:
	Bias(const Parameters&, bool, const std::vector<Graph>&);
		
	void update_cv(const std::vector<Polymer>&);
	void update_cv_rew(const std::vector<Polymer>&, double);
	double get_cv() const;

	void update_bias(const std::vector<Polymer>&,double,double);
	void restore_splines_transient(double);
	
	double energy_diff(const std::vector<Polymer>& pols) const;
	Force calc_force(const std::vector<Polymer>&, int, int) const;
	double get_rew_factor() const;
	double get_rew_factor_avg() const;
	void set_new_block();
	void update_rew_factor_avg(int);
	double get_rew_factor_block() const;
	void set_rew_factor_avg(double,double);
	const double get_gauss_width() const;
	double get_count() const;
	const bool metad_on;
	void add_gaussian(double,double);
	std::vector<double> get_heights() const;
	std::vector<double> get_centers() const;
	
private:
	double cv;
	double transient; //c(t)
	double rew_factor; // exp(beta(V(s)-c(t)))
	double rew_factor_avg; //total instantaneous average of all previous values
	double rew_factor_block; //sum in current block
	double count; //counts #measurements of rew_factor
	const int id;
	const double bias_factor;
	const double gauss_width;
	const double first_height;
	const int biased_graph;
	const int sign;
	const double exc_const;
	const double exponent_factor;
	std::vector<double> heights;
	std::vector<double> cv_centers;
	double regularization;
	
	const std::vector<Graph>& graphs;
	
	Spline v_spline;
	Spline vder_spline;
	const double spline_step;
	
	double calc_bias(double) const; //calculate bias explicitly (only used when creating spline)
	double calc_bias_der(double) const;
	double coll_var(const std::vector<Polymer>&) const;
	double coll_var_der(const std::vector<Polymer>&) const;
	double scalar_product(const std::vector<Polymer>&, int) const;
	double scalar_product_conn(const std::vector<Polymer>&, int, int) const;
	double gaussian(double, double, double) const;
	Force cv_grad(const std::vector<Polymer>&, int, int) const;
	double sum_exp(const std::vector<Polymer>&, int) const;
	double sum_exp_distcorr(const std::vector<Polymer>&) const;
	Force two_terms(const std::vector<Polymer>&, int, int) const;
	//Force two_terms_conn(const std::vector<Polymer>&, int, int, int) const;
	double sq_distAB(const std::vector<Polymer>&) const;
	void create_splines();
	void update_transient(double);
	double wall_force_magn(double) const;
	double wall_potential(double) const;
	
	//std::ofstream heights_file;
	std::ofstream cv_centers_file;
	
	const double wall_id;
	const double wall_energy;
	const double wall_pos;

};

#endif
