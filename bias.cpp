
#include "bias.hpp"
#include "point.hpp"
#include <algorithm>

Bias::Bias(const Parameters& params) : id(params.cv_id), sign(params.sign),
			gauss_width(params.gauss_width), bias_factor(params.bias_factor),
			exc_const(params.exc_const), first_height(params.first_height),
			exponent_factor(1.0/(2*params.gauss_width*params.gauss_width)),
			metad_on(params.metad_on),spline_step(gauss_width/10.0)
{
	cv_centers_file.open("cv_centers.dat");
	heights_file.open("heights.dat");
	v_spline = Spline(spline_step);
	vder_spline = Spline(spline_step);
	rew_factor = 1;
	rew_factor_avg = 1;
	count = 0;
}
	
double Bias::coll_var(const std::vector<Polymer>& pols) const
{
	double tmp = 0;
	switch(id)
	{
		case 1:
			tmp = sum_exp(pols);
			return 1.0+sign*tmp/pols[0].num_beads;
		case 2:
			tmp = sum_exp(pols);
			return std::log(tmp);
		default:
			return pols[0][0][0];
	}
}

double Bias::sum_exp(const std::vector<Polymer>& pols) const
{
	double tmp = 0;
	for(int bead=0; bead<pols[0].num_beads; ++bead)
		tmp += std::exp(-exc_const * scalar_product(pols,bead));
	return tmp;
}

Force Bias::cv_grad(const std::vector<Polymer>& pols, int bead, int part) const
{
	Force tmp(pols[0][0].size());
	switch(id)
	{
		case 1:
			tmp += two_terms(pols, bead, part);
			return tmp * exc_const/pols[part].num_beads * sign*(-1)*std::pow(-1,part);
		case 2:
			tmp += two_terms(pols, bead, part);
			tmp /= sum_exp(pols);
			return tmp * exc_const * sign*std::pow(-1,part);
		default:
			return tmp;
	}
}

Force Bias::two_terms(const std::vector<Polymer>& pols, int bead, int part) const
{
	return (pols[0][bead+1]-pols[1][bead+1])*std::exp(-exc_const*scalar_product(pols,bead)) +
		   (pols[0][bead-1]-pols[1][bead-1])*std::exp(-exc_const*scalar_product(pols,bead-1));			
}

Force Bias::calc_force(const std::vector<Polymer>& pols, int bead, int part) const
{	
	if(!metad_on)
		return Force(pols[0][0].size());
	//std::cout << "before eval" << std::endl;
	//double tmp = vder_spline.eval_spline(cv);
	//std::cout << "after eval" << std::endl;
	return vder_spline.eval_spline(cv)*cv_grad(pols,bead,part);
}

double Bias::scalar_product(const std::vector<Polymer>& pols, int bead) const
{
	return (pols[0][bead]-pols[1][bead])*(pols[0][bead+1]-pols[1][bead+1]);
}

double Bias::calc_bias(double cv) const
{
	double tmp = 0;
	for(int index=0; index<heights.size(); ++index)
		tmp += gaussian(cv,cv_centers[index],heights[index]);
	return tmp;
}

double Bias::calc_bias_der(double cv) const
{
	double tmp=0;
	for(int index=0; index<heights.size(); ++index)
		tmp += (cv-cv_centers[index])*gaussian(cv,cv_centers[index],heights[index]);
	return tmp/std::pow(gauss_width,2);
}

double Bias::gaussian(double cv, double center, double height) const
{
	return height*std::exp(- std::pow(cv-center,2)*exponent_factor);
}

void Bias::update_bias(const std::vector<Polymer>& pols, double beta, double t)
{
	if(metad_on)
	{
		double cv = coll_var(pols);
		cv_centers.push_back(cv);
		double h = first_height * std::exp(-beta/(bias_factor-1.0) *calc_bias(cv));
		heights.push_back(h);
		cv_centers_file << t << "\t" << cv << std::endl;
		heights_file << t << "\t" << h << std::endl;	
		update_transient(beta);
		create_splines();
	}
}

void Bias::create_splines()
{
	if(metad_on)
	{
		std::vector<double> vs;
		std::vector<double> vders;
		std::vector<double> cvs;
		auto it = std::min_element(std::begin(cv_centers),std::end(cv_centers));
		double min = *it - 5*gauss_width;
		it = std::max_element(std::begin(cv_centers),std::end(cv_centers));
		double max = *it + 5*gauss_width;
		for(double s=min; s<=max; s+=spline_step)
		{
			cvs.push_back(s);
			vs.push_back(calc_bias(s));
			vders.push_back(calc_bias_der(s));
		}
		v_spline.create_spline(cvs,vs);
		vder_spline.create_spline(cvs,vders);
	}
}


void Bias::update_transient(double beta)
{
	int num_samples = 2000;
	double smin = v_spline.get_min();
	double smax = v_spline.get_max();
	double step = (smax-smin)/num_samples;
	const double num_const = beta*bias_factor/(bias_factor-1.0);
	const double den_const = bias_factor/(bias_factor-1.0);
	double tmp_avg, tmp_int_num, tmp_int_den = 0;
	for(double s = smin; s<smax; s+=step)
		tmp_avg += v_spline.eval_spline(s);
	tmp_avg *= step/(smax-smin);
	for(double s = smin; s<smax; s+=step)
	{
		tmp_int_num += std::exp(num_const*(v_spline.eval_spline(s)-tmp_avg));
		tmp_int_den += std::exp(den_const*(v_spline.eval_spline(s)-tmp_avg));
	}
	transient = beta*tmp_avg + std::log(tmp_int_num/tmp_int_den);
}

void Bias::update_cv(const std::vector<Polymer>& pols, const double beta)
{
	cv = coll_var(pols);
	rew_factor = std::exp(beta*(v_spline.eval_spline(cv)-transient));
	rew_factor_avg = (rew_factor_avg*count + rew_factor)/(count+1.0);
	++count;
}

void Bias::set_rew_factor_avg(double rew_factor_avg_, double count_)
{
	rew_factor_avg = rew_factor_avg_;
	count = count_;
}

double Bias::get_cv() const
{
	return cv;
}

double Bias::get_rew_factor() const
{
	return rew_factor;
}

double Bias::get_rew_factor_avg() const
{
	return rew_factor_avg;
}

const double Bias::get_gauss_width() const
{
	return gauss_width;
}

double Bias::get_count() const
{
	return count;
}

void Bias::add_gaussian(double height, double center)
{
	heights.push_back(height);
	cv_centers.push_back(center);
}

std::vector<double> Bias::get_heights() const
{
	return heights;
}

std::vector<double> Bias::get_centers() const
{
	return cv_centers;
}