
#include "bias.hpp"
#include "point.hpp"
#include <algorithm>

Bias::Bias(const Parameters& params, bool cont_sim, const std::vector<Graph>& graphs_in) : 
			id(params.cv_id), sign(params.sign), graphs(graphs_in),
			gauss_width(params.gauss_width), bias_factor(params.bias_factor),
			exc_const(params.exc_const), first_height(params.first_height),
			exponent_factor(1.0/(2*params.gauss_width*params.gauss_width)),
			gauss2(params.gauss_width*params.gauss_width),
			metad_on(params.metad_on),spline_step(params.spline_step),
			v_spline(Spline(params.spline_step,0)),vder_spline(Spline(params.spline_step,0)),
			wall_id(params.wall_id),wall_energy(params.wall_energy),
			wall_pos(params.wall_pos), biased_graph(params.biased_graph),
			reference_graph(params.reference_graph),
			border_region(8*params.gauss_width)
{
	if(cont_sim)
	{
		cv_centers_file.open("cv_centers.dat",std::ios_base::app);
		//heights_file.open("heights.dat",std::ios_base::app);		
	}
	else
	{
		cv_centers_file.open("cv_centers.dat");
		//heights_file.open("heights.dat");
	}
	rew_factor = 1;
	rew_factor_avg = 0;
	rew_factor_block = 0;
	transient=0;
	count = 0;
	regularization = 1e-4;
}

double Bias::energy_diff(const std::vector<Polymer>& pols) const
{
	return graphs[biased_graph].energy_diff(pols,graphs[reference_graph]);
}
	
double Bias::coll_var(const std::vector<Polymer>& pols) const
{
	if(sign==0)
		return 0;
	double tmp = 0;
	double dist;
	double s;
	//double pos_weight;
	//double neg_weight;
	switch(id)
	{
		/*case 1:
			if(pols[0].connected)
				return e_s+sign;
			else
				return 1+sign*e_s;
			//tmp = std::exp(-exc_const*scalar_product(pols,pols[0].num_beads-1));
			//if(pols[0].connected)
				//return 1.0/tmp + sign;
			//return 1.0+sign*tmp;
			//tmp = sum_exp(pols,pols[0].num_beads);
			//return 1.0+sign*tmp/pols[0].num_beads;*/
		case 2:
			return std::log(pos_weight + neg_weight);
		case 3: //distance-corrected cv
			tmp = sum_exp_distcorr(pols);
			return -std::log(tmp/pols[0].num_beads);
		case 4:
			return std::log(std::exp(-graphs[0].energy_absolute(pols))
							+std::exp(-graphs[0].energy_absolute(pols)));
		case 5:
			return -std::log(neg_weight/pos_weight);
		case 6:
			return graphs[1].energy_diff(pols,graphs[0]);
			/*s = energy_diff(pols);
			return 0.5*(1.1*s-0.9*s*std::tanh(s));
		case 7:
			/*pos_weight = 0;
			neg_weight = 0;
			for(const Graph& graph : graphs)
			{
				pos_weight += graph.get_weight(pols,true);
				neg_weight += graph.get_weight(pols,false);
			}*/
			//return (pos_weight-neg_weight)/(pos_weight+neg_weight);
			/*s = energy_diff(pols);
			return 0.5*(1.1*s+20-(0.9*s-20)*std::tanh(0.2*(s-20)));*/
		case 8:
			return graphs[biased_graph].energy_diff(pols,graphs[reference_graph]);
		case 9:
			return std::log(1+neg_weight/pos_weight);
		default:
			std::cout << "Not a valid CV option" << std::endl;
			return pols[0][0][0];
	}
}

void Bias::set_weights(double pos, double neg)
{
	pos_weight = pos;
	neg_weight = neg;
}

double Bias::sum_exp(const std::vector<Polymer>& pols, int num_beads) const
{
	double tmp = 0;
	if(pols[0].connected)
		for(int bead=0; bead<num_beads; ++bead)
			tmp += std::exp(exc_const * scalar_product(pols,bead));
	else
		for(int bead=0; bead<num_beads; ++bead)
			tmp += std::exp(-exc_const * scalar_product(pols,bead));
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


Force Bias::cv_grad(const std::vector<Polymer>& pols, int bead, int part) const
{
	Force tmp(pols[0][0].size());
	Force tmp2(pols[0][0].size());
	double tmp3;
	//double pos_weight=0;
	//double neg_weight=0;
	int num_beads = pols[0].num_beads;
	int tmp_bead;
	int sign2=1;
	switch(id)
	{
		/*case 1:
			if(pols[0].connected)
				return std::pow(-1,part) * exc_const/num_beads * two_terms(pols,bead,part);
			else
				return -sign * std::pow(-1,part) * exc_const/num_beads * two_terms(pols,bead,part);*/
		case 2:
			for(const Graph& graph : graphs)
			{
				tmp += graph.get_grad_weight(pols,graphs[0], true,bead,part); //grad W+
				tmp2 += graph.get_grad_weight(pols,graphs[0], false,bead,part); //grad W-
			}
			return (tmp + tmp2)/(pos_weight+neg_weight);
		case 3:
			tmp2 = pols[0][bead]-pols[1][bead]-(pols[0][bead+1]-pols[1][bead+1]);
			tmp = tmp2*std::exp(0.5*exc_const*tmp2.sqdist0());
			tmp2 = pols[0][bead]-pols[1][bead]-(pols[0][bead-1]-pols[1][bead-1]);
			tmp += tmp2*std::exp(0.5*exc_const*tmp2.sqdist0());
			return tmp*(-std::pow(-1,part)*exc_const/sum_exp_distcorr(pols));
		case 4:
			for(const Graph& graph : graphs)
			{
				tmp += graph.get_grad_weight(pols,graphs[current_graph_id], true,bead,part); //grad W+
				tmp2 += graph.get_grad_weight(pols,graphs[current_graph_id], false,bead,part); //grad W-
			}
			return (tmp + tmp2)/(pos_weight+neg_weight);			
		case 5:
			for(const Graph& graph : graphs)
			{
				tmp += graph.get_grad_weight(pols,graphs[current_graph_id],true,bead,part); //grad W+
				tmp2 += graph.get_grad_weight(pols,graphs[current_graph_id],false,bead,part); //grad W-
			}
			return tmp/pos_weight - tmp2/neg_weight;
		case 6:
			return graphs[1].energy_diff_grad(pols,graphs[0],bead,part);
		case 7:
			for(const Graph& graph : graphs)
			{
				tmp += graph.get_grad_weight(pols,graphs[current_graph_id],true,bead,part);
				tmp2 += graph.get_grad_weight(pols,graphs[current_graph_id],false,bead,part);
			}
			return (tmp*neg_weight - tmp2*pos_weight)/std::pow(pos_weight+neg_weight,2);
		case 8:
			return graphs[biased_graph].energy_diff_grad(pols,graphs[reference_graph],bead,part);
		case 9:
			for(const Graph& graph : graphs)
			{
				tmp += graph.get_grad_weight(pols,graphs[current_graph_id],true,bead,part);
				tmp2 += graph.get_grad_weight(pols,graphs[current_graph_id],false,bead,part);
			}
			return (tmp2 - neg_weight/pos_weight * tmp)/(pos_weight + neg_weight);
		default:
			std::cout << "No gradCV is implemented!" << std::endl;
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


Force Bias::calc_force(const std::vector<Polymer>& pols, int bead, int part) const
{	
	if((!metad_on))
		return Force(pols[0][0].size());
	return (vder_spline.eval_spline(cv)+wall_force_magn(cv))*cv_grad(pols,bead,part);
}

double Bias::scalar_product(const std::vector<Polymer>& pols, int bead) const
{
	return (pols[0][bead]-pols[1][bead])*(pols[0][bead+1]-pols[1][bead+1]);
}

double Bias::calc_bias(double s) const
{
	double tmp = 0;
	for(int index=0; index<heights.size(); ++index)
		tmp += gaussian(s,cv_centers[index],heights[index]);
	return tmp;
}

double Bias::calc_bias2(double s) const
{
	return gaussian(s,latest_cv_center,latest_height) + v_spline.eval_spline(s);
}

double Bias::calc_bias_der(double s) const
{
	double tmp=0;
	for(int index=0; index<heights.size(); ++index)
		tmp += (s-cv_centers[index])*gaussian(s,cv_centers[index],heights[index]);
	return tmp/gauss2;
}

double Bias::calc_bias_der2(double s) const
{
	double tmp =(s-latest_cv_center)*gaussian(s,latest_cv_center,latest_height)/gauss2
				 + vder_spline.eval_spline(s);
	return tmp;
}

double Bias::gaussian(double s, double center, double height) const
{
	if(std::abs(s-center)>border_region)
		return 0;
	return height*std::exp(- std::pow(s-center,2)*exponent_factor);
}

void Bias::update_bias(const std::vector<Polymer>& pols, double beta, double t)
{
	if(metad_on)
	{
		latest_cv_center = coll_var(pols);
		latest_height = first_height * std::exp(-beta/(bias_factor-1.0)*v_spline.eval_spline(latest_cv_center));
		cv_centers_file << t << "\t" << latest_cv_center << "\t" << latest_height << std::endl;
		create_splines();
		update_transient(beta);
	}
	update_cv(pols);
}

void Bias::create_splines()
{
	if(metad_on)
	{
		std::vector<double> vs;
		std::vector<double> vders;
		std::vector<double> cvs;
		//double min = 0;
		//double max = 0;
		double min = v_spline.get_min();
		double max = v_spline.get_max();
		if((latest_cv_center) < (min+border_region))
			min = latest_cv_center - border_region;
		if((latest_cv_center) > (max-border_region))
			max = latest_cv_center + border_region;
		int num_steps_for_spline = std::floor((max-min)/spline_step)+1;
		//if(num_steps_for_spline != num_steps_for_spline)
		//	std::cout << "Have a look at the Bias::create_splines method!" << std::endl;
		cvs.resize(num_steps_for_spline);
		vders.resize(num_steps_for_spline);
		vs.resize(num_steps_for_spline);
		/*auto it = std::min_element(std::begin(cv_centers),std::end(cv_centers));
		double min = *it - 5*gauss_width;
		it = std::max_element(std::begin(cv_centers),std::end(cv_centers));
		double max = *it + 5*gauss_width;*/
		min = std::round(min/spline_step)*spline_step;
		double s = min;
		for(int i=0; i<num_steps_for_spline; ++i)
		{
			s = min+i*spline_step;
			cvs[i]=s;
			vs[i]=calc_bias2(s);
			vders[i]=calc_bias_der2(s);
		}
		v_spline.create_spline(cvs,vs);
		vder_spline.create_spline(cvs,vders);
	}
}


void Bias::update_transient(double beta)
{
	double smin = v_spline.get_min()+2*gauss_width;
	double smax = v_spline.get_max()-2*gauss_width;
	double step = 0.02*gauss_width;
	//int num_steps = round((smax-smin)/step);
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
		double diff = v_spline.eval_spline(s)-tmp_avg;
		tmp_int_num += std::exp(num_const*diff);
		tmp_int_den += std::exp(den_const*diff);
	}
	transient = beta*tmp_avg + std::log(tmp_int_num/tmp_int_den);
}

void Bias::update_cv(const std::vector<Polymer>& pols)
{
	cv = coll_var(pols);
}

void Bias::update_cv_rew(const std::vector<Polymer>& pols, double beta)
{
	cv = coll_var(pols);
	if(metad_on)
	{
		rew_factor = std::exp(beta*(v_spline.eval_spline(cv)+wall_potential(cv))-transient);
	}
	else
		rew_factor = 1;
	//rew_factor_avg = (rew_factor_avg*count + rew_factor)/(count+1.0);
	rew_factor_block += rew_factor;
	++count;
}

void Bias::update_rew_factor_avg(int block)
{
	rew_factor_avg = (rew_factor_avg*block + rew_factor_block/count)/(block+1.0);
	count = 0;
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
			return wall_energy*(wall_pos-cv);
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

void Bias::set_current_graph(int g_id)
{
	current_graph_id = g_id;
}
