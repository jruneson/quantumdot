
#include "bias.hpp"
#include "point.hpp"
#include <algorithm>

Bias::Bias(const Parameters& params, bool cont_sim) : id(params.cv_id), sign(params.sign),
			gauss_width(params.gauss_width), bias_factor(params.bias_factor),
			exc_const(params.exc_const), first_height(params.first_height),
			exponent_factor(1.0/(2*params.gauss_width*params.gauss_width)),
			metad_on(params.metad_on),spline_step(gauss_width/10.0),
			wall_id(params.wall_id),wall_energy(params.wall_energy),
			wall_pos(params.wall_pos)
{
	if(cont_sim)
	{
		cv_centers_file.open("cv_centers.dat",std::ios_base::app);
		heights_file.open("heights.dat",std::ios_base::app);		
	}
	else
	{
		cv_centers_file.open("cv_centers.dat");
		heights_file.open("heights.dat");
	}
	v_spline = Spline(spline_step);
	vder_spline = Spline(spline_step);
	rew_factor = 1;
	rew_factor_avg = 1;
	rew_factor_block = 0;
	transient=0;
	count = 0;
	regularization = 1e-5;
}

double Bias::energy_diff(const std::vector<Polymer>& pols) const
{
	/*if(pols[0].connected)
	{
		int num_beads_1p = pols[0].num_beads/2;
		return -std::log(sum_exp(pols,num_beads_1p)/num_beads_1p);
	}
	else
	{*/
		//int num_beads = pols[0].num_beads;
		//return -std::log(sum_exp(pols,num_beads)/num_beads);
	//}
	/*if(pols[0].connected)
		return exc_const*scalar_product(pols,pols[0].num_beads);
	else*/
		return exc_const*scalar_product(pols,pols[0].num_beads-1);
}
	
double Bias::coll_var(const std::vector<Polymer>& pols, double e_s) const
{
	double tmp = 0;
	double dist;
	double s;
	switch(id)
	{
		case 1:
			if(pols[0].connected)
				return e_s+sign;
			else
				return 1+sign*e_s;
			//tmp = std::exp(-exc_const*scalar_product(pols,pols[0].num_beads-1));
			//if(pols[0].connected)
				//return 1.0/tmp + sign;
			//return 1.0+sign*tmp;
			//tmp = sum_exp(pols,pols[0].num_beads);
			//return 1.0+sign*tmp/pols[0].num_beads;
		case 2: //energy-difference cv
			/*if(pols[0].connected)
				return std::log(e_s);
			else
			{*/
				//std::cout << std::log(std::abs(exc_factor-1)) << std::endl; 
			return -std::log(e_s);
			//}
			//return energy_diff(pols);
		case 3: //distance-corrected cv
			tmp = sum_exp_distcorr(pols);
			return -std::log(tmp/pols[0].num_beads);
		case 4:
			if(pols[0].connected)
				return -std::log(std::abs(e_s+sign)+regularization);
			else
				return -std::log(std::abs(1+sign*e_s)+regularization);
			/*dist = sq_distAB(pols);
			for(int bead=0; bead<pols[0].num_beads; ++bead)
				tmp += std::exp(-exc_const*(scalar_product(pols,bead)-dist));
			return -std::log(tmp/pols[0].num_beads);*/
		case 5:
			s = energy_diff(pols);
			return 0.5*(1+s+(1-s)*std::tanh(s));
		case 6:
			s = energy_diff(pols);
			return 0.5*(1.1*s-0.9*s*std::tanh(s));
		case 7:
			s = energy_diff(pols);
			return 0.5*(1.1*s+20-(0.9*s-20)*std::tanh(0.2*(s-20)));
		default:
			std::cout << "Not a valid CV option" << std::endl;
			return pols[0][0][0];
	}
}

double Bias::sum_exp(const std::vector<Polymer>& pols, int num_beads) const
{
	double tmp = 0;
	/*if(pols[0].connected)
	{
		for(int bead=0; bead<num_beads; ++bead)
			tmp += std::exp(-exc_const * scalar_product_conn(pols,bead,num_beads));
	}
	else
	{*/
	if(pols[0].connected)
		for(int bead=0; bead<num_beads; ++bead)
			tmp += std::exp(exc_const * scalar_product(pols,bead));
	else
		for(int bead=0; bead<num_beads; ++bead)
			tmp += std::exp(-exc_const * scalar_product(pols,bead));
	//}
	return tmp;
}

double Bias::sum_exp_distcorr(const std::vector<Polymer>& pols) const
{
	double tmp=0;
	for(int bead=0; bead<pols[0].num_beads; ++bead)
	{
		Point tmp2 = pols[0][bead]-pols[1][bead]-(pols[0][bead+1]-pols[1][bead+1]);
		tmp += std::exp(0.5*exc_const*tmp2.sqdist0());
	}
	return tmp;
}

double Bias::sq_distAB(const std::vector<Polymer>& pols) const
{
	double tmp=0;
	for(int bead=0; bead<pols[0].num_beads; ++bead)
		tmp += (pols[0][bead]-pols[1][bead])*(pols[0][bead]-pols[1][bead]);
	tmp /= (pols[0].num_beads);
	return tmp;
}
/*
double Bias::sq_distAB(const std::vector<Polymer>& pols) const
{
	Point tmp(pols[0][0].size());
	for(int bead=0; bead<pols[0].num_beads; ++bead)
		tmp += pols[0][bead]-pols[1][bead];
	tmp /= pols[0].num_beads;
	return tmp*tmp;
}*/

Force Bias::cv_grad(const std::vector<Polymer>& pols, int bead, int part,double e_s) const
{
	Force tmp(pols[0][0].size());
	Force tmp2(pols[0][0].size());
	double tmp3;
	int num_beads = pols[0].num_beads;
	int tmp_bead;
	int sign2=1;
	switch(id)
	{
		case 1:
			if(pols[0].connected)
				return std::pow(-1,part) * exc_const/num_beads * two_terms(pols,bead,part);
			else
				return -sign * std::pow(-1,part) * exc_const/num_beads * two_terms(pols,bead,part);
			/*if(bead==0)
				tmp_bead = pols[part].num_beads-1;
			else if(bead==pols[part].num_beads-1)
				tmp_bead = 0;
			else
				return tmp;
			tmp = exc_const*(pols[part][tmp_bead]-pols[1-part][tmp_bead]);
			if(pols[part].connected)
				return (-1)*sign*std::exp(-exc_const*scalar_product(pols,pols[part].num_beads-1))*tmp;
			else
				return sign*std::exp(exc_const*scalar_product(pols,pols[part].num_beads-1))*tmp;*/
			//tmp += two_terms(pols, bead, part);
			//return tmp * exc_const/pols[part].num_beads * std::pow(-1,part);
		case 2:
			tmp += two_terms(pols,bead,part);
			/*if(pols[0].connected)
				tmp /= (num_beads*e_s);
			else*/
			tmp /= (num_beads*e_s);
			//tmp /= sum_exp(pols,num_beads);
			//std::cout << (num_beads*e_s) << "\t" << sum_exp(pols,num_beads) << std::endl;
			return exc_const*std::pow(-1,part)*tmp;
			/*if(bead==0)
				tmp_bead = pols[part].num_beads-1;
			else if(bead==pols[part].num_beads-1)
				tmp_bead = 0;
			else
				return tmp; 
			return exc_const*(pols[part][tmp_bead]-pols[1-part][tmp_bead]);*/
		case 5:
		case 6:
		case 7:
			/*if(pols[0].connected)
			{
				num_beads /= 2;
				tmp += two_terms_conn(pols, bead, part,num_beads);
			}
			else */
			tmp += two_terms(pols,bead,part);
			tmp /= sum_exp(pols,num_beads);
			if(id==5)
			{
				double t = std::tanh(cv);
				tmp *= 0.5*(1-t+(1-cv)*(1-t*t));
			}
			else if(id==6)
			{
				double t = std::tanh(cv);
				tmp *= 0.5*(1.1-0.9*(t + cv*(1-t*t)));
			}
			else if(id==7)
			{
				double t = std::tanh(0.2*(cv-20));
				tmp *= 0.5*(1.1-0.9*t-0.2*(0.9*cv-20)*(1-t*t));
			}
			return tmp * (exc_const * std::pow(-1,part));
		case 3:
			tmp2 = pols[0][bead]-pols[1][bead]-(pols[0][bead+1]-pols[1][bead+1]);
			tmp = tmp2*std::exp(0.5*exc_const*tmp2.sqdist0());
			tmp2 = pols[0][bead]-pols[1][bead]-(pols[0][bead-1]-pols[1][bead-1]);
			tmp += tmp2*std::exp(0.5*exc_const*tmp2.sqdist0());
			return tmp*(-std::pow(-1,part)*exc_const/sum_exp_distcorr(pols));
		case 4:
			tmp = two_terms(pols,bead,part);
			tmp3 = std::pow(-1,part)*exc_const/num_beads;
			if(pols[0].connected)
			{
				if((e_s+sign)<0)
					sign2 = -1;
				return (-tmp3 / (e_s + sign + sign2*regularization)) * tmp;
			}
			else
			{
				if((1+sign*e_s)<0)
					sign2 = -1;
				return (sign*tmp3 / (1+sign*e_s + sign2*regularization)) * tmp;
			}
			/*tmp2 = pols[0][bead+1]-pols[1][bead+1]-(pols[0][bead]-pols[1][bead])*2.0/pols[0].num_beads;
			tmp = tmp2*std::exp(-exc_const*scalar_product(pols,bead));
			tmp2 = pols[0][bead-1]-pols[1][bead-1]-(pols[0][bead]-pols[1][bead])*2.0/pols[0].num_beads;
			tmp += tmp2*std::exp(-exc_const*scalar_product(pols,bead-1));
			return tmp*(exc_const*std::pow(-1,part)/sum_exp(pols,pols[0].num_beads));*/
		default:
			return tmp;
	}
}

Force Bias::two_terms(const std::vector<Polymer>& pols, int bead, int part) const
{
	if(pols[0].connected)
	{
		return (pols[0][bead+1]-pols[1][bead+1])*std::exp(exc_const*scalar_product(pols,bead)) +
			   (pols[0][bead-1]-pols[1][bead-1])*std::exp(exc_const*scalar_product(pols,bead-1));
	}
	else
	{
		return (pols[0][bead+1]-pols[1][bead+1])*std::exp(-exc_const*scalar_product(pols,bead)) +
		       (pols[0][bead-1]-pols[1][bead-1])*std::exp(-exc_const*scalar_product(pols,bead-1));			
	}
}

/*Force Bias::two_terms_conn(const std::vector<Polymer>& pols, int bead, int part) const
{
	return (pols[0][bead+1]-pols[0][bead+1])*std::exp(-exc_const*scalar_product_conn(pols,bead,num_beads)) +
		   (pols[0][bead-1]-pols[0][bead-1])*std::exp(-exc_const*scalar_product_conn(pols,bead-1,num_beads));
}*/

Force Bias::calc_force(const std::vector<Polymer>& pols, int bead, int part,double e_s) const
{	
	if(!metad_on)
		return Force(pols[0][0].size());
	//std::cout << "before eval" << std::endl;
	//double tmp = vder_spline.eval_spline(cv);
	//std::cout << "after eval" << std::endl;
	return (vder_spline.eval_spline(cv)+wall_force_magn(cv))*cv_grad(pols,bead,part,e_s);
}

double Bias::scalar_product(const std::vector<Polymer>& pols, int bead) const
{
	return (pols[0][bead]-pols[1][bead])*(pols[0][bead+1]-pols[1][bead+1]);
}

/*double Bias::scalar_product_conn(const std::vector<Polymer>& pols, int bead, int num_beads) const
{
	return (pols[0][bead]-pols[0][bead+num_beads])*(pols[0][bead+1]-pols[0][bead+num_beads+1]);
}*/

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
	/*if(id==5)
	{
		double tan = std::tanh(cv);
		tmp *= 0.5*(1-tan+(1-cv)*(1-tan*tan));
	}*/
	return tmp/std::pow(gauss_width,2);
}

double Bias::gaussian(double cv, double center, double height) const
{
	return height*std::exp(- std::pow(cv-center,2)*exponent_factor);
}

void Bias::update_bias(const std::vector<Polymer>& pols, double beta, double t, double e_s)
{
	if(metad_on)
	{
		double cv_now = coll_var(pols,e_s);
		double h = first_height * std::exp(-beta/(bias_factor-1.0)*calc_bias(cv_now));
		cv_centers.push_back(cv_now);
		heights.push_back(h);
		cv_centers_file << t << "\t" << cv << std::endl;
		heights_file << t << "\t" << h << std::endl;
		create_splines();
		update_transient(beta);
	}
	update_cv(pols,e_s);
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
		//assign zeros
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
	double smin = v_spline.get_min();
	double smax = v_spline.get_max();
	double step = 0.01*gauss_width;
	int num_samples = round((smax-smin)/step);
	const double num_const = beta*bias_factor/(bias_factor-1.0);
	const double den_const = beta/(bias_factor-1.0);
	double tmp_avg=0;
	double tmp_int_num=0;
	double tmp_int_den = 0;
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

void Bias::update_cv(const std::vector<Polymer>& pols, double e_s)
{
	cv = coll_var(pols, e_s);
}

void Bias::update_cv_rew(const std::vector<Polymer>& pols, double beta, double e_s)
{
	cv = coll_var(pols,e_s);
	if(metad_on)
	{
		rew_factor = std::exp(beta*(v_spline.eval_spline(cv)+wall_potential(cv)-transient));
		/*rew_factor *= std::exp(beta*wall_potential(cv));
		if(std::exp(beta*wall_potential(cv))>1.1)
			std::cout << rew_factor << "\t" << std::exp(beta*wall_potential(cv)) << std::endl;*/
	}
	else
		rew_factor = 1;
	rew_factor_avg = (rew_factor_avg*count + rew_factor)/(count+1.0);
	rew_factor_block += rew_factor;
	++count;
}

void Bias::set_rew_factor_avg(double rew_factor_avg_, double count_)
{
	rew_factor_avg = rew_factor_avg_;
	count = count_;
}

double Bias::wall_force_magn(double cv) const
{
	if(wall_id==0)
		return 0;
	if(wall_id==1)
		if(cv>wall_pos)
			return -wall_energy*(cv-wall_pos);
	if(wall_id==2)
		if(cv<wall_pos)
			return -wall_energy*(wall_pos-cv);
	return 0;
}

double Bias::wall_potential(double cv) const
{
	if(wall_id==0)
		return 0;
	if(wall_id==1)
		if(cv>wall_pos)
			return 0.5*wall_energy*std::pow(cv-wall_pos,2);
	if(wall_id==2)
		if(cv<wall_pos)
			return 0.5*wall_energy*std::pow(wall_pos-cv,2);
	return 0;
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

void Bias::set_new_block()
{
	rew_factor_block = 0;
}

double Bias::get_rew_factor_block() const
{
	return rew_factor_block;
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

void Bias::restore_splines_transient(double beta) //Used when a simulation is continued. Must be called after the old gaussians are added
{
	create_splines();
	update_transient(beta);
}

std::vector<double> Bias::get_heights() const
{
	return heights;
}

std::vector<double> Bias::get_centers() const
{
	return cv_centers;
}
