#include <algorithm>
#include "graph.hpp"


Graph::Graph(const Parameters& params, int graph_id) //Chain lengths given as number of chains with length 1, followed by 2, 3 etc.
 : exc_const(params.exc_const), id(graph_id), spin_proj(params.spin_proj)
{
	int num_parts = params.num_parts;
	int spin = params.spin;
	exc_bead = params.num_beads-1;
	if(params.connected)
		exponent_sign = 1;
	else
		exponent_sign = -1;
	mult = 0;
	positive = true;
	if(num_parts==1)
		mult = 1;
	if(num_parts==2)
	{
		//std::vector<int> graph1 = {2,0};
		//std::vector<int> graph2 = {0,1};
		mult = 1;
		if(graph_id==1) //{2,0}
		{
			for(int i=0; i<=1; ++i)
			{
				std::vector<int> tmp_ch = {i};
				chains.push_back(tmp_ch);
			} 
			//chains.push_back(std::vector<int>(0));
			//chains.push_back(std::vector<int>(1));
		}
		if(graph_id==2) // {0,1}
		{
			std::vector<int> tmp_ch = {0,1};
			chains.push_back(tmp_ch);
			
			//replace by pair calculation method
			std::vector<std::pair<int,int>> tmp;
			tmp.push_back(std::pair<int,int>(0,1));
			exchange_pairs.push_back(tmp);
			
			if(spin==2)
				positive = false;
		}
	}
	if(num_parts==3)
	{
		if(graph_id==1) //{3,0,0}
		{
			mult = 1;
			for(int i=0; i<3; ++i)
			{
				std::vector<int> tmp_ch = {i};
				chains.push_back(tmp_ch);
			}
		}
		if(graph_id==2) //{1,1,0}
		{
			std::vector<int> tmp_ch = {0};
			chains.push_back(tmp_ch);
			tmp_ch = {1,2};
			chains.push_back(tmp_ch);
			std::vector<std::pair<int,int>> tmp;
			tmp.push_back(std::pair<int,int>(1,2));//cyclic(n+1,num_parts)));
			exchange_pairs.push_back(tmp);
			if(spin==3)
			{
				//both for spin_proj true and false
				mult = 3;
				positive = false;
			}
		}
		if(graph_id==3) //{0,0,1}
		{
			std::vector<int> tmp_ch = {0,1,2};
			chains.push_back(tmp_ch);
			std::vector<std::pair<int,int>> tmp;
			for(int n=0; n<num_parts; ++n)
				tmp.push_back(std::pair<int,int>(n,cyclic(n+1,num_parts)));
			exchange_pairs.push_back(tmp);
			if(spin==1)
			{
				if(spin_proj)
					mult = 0;
				else
				{	
					mult = 1;
					positive = false;
				}
			}
			if(spin==3)
			{
				if(spin_proj)
					mult = 1;
				else
					mult = 2;
			}
		}		
	}
}

int Graph::cyclic(int i, int n) const
{
	return (i+n)%n;
}

double Graph::get_weight(const std::vector<Polymer>& pols, bool use_positive) const
{
	if(positive == use_positive)
		return std::abs(get_weight_signed(pols));
	return 0;
}

double Graph::get_weight_signed(const std::vector<Polymer>& pols) const
{
	if(exchange_pairs.size()==0)
		return 1;
	double tmp = 0;
	for(const auto& list : exchange_pairs)
		tmp += std::exp(exponent_sign * calc_exponent(pols,list));
	return get_sign()*mult*tmp/exchange_pairs.size();
}

double Graph::calc_exponent(const std::vector<Polymer>& pols, const std::vector<std::pair<int,int>>& list) const
{
	double tmp = 0;
	for(const auto& pair : list)
	{
		tmp += pols[std::get<0>(pair)][exc_bead].sqdist(pols[std::get<1>(pair)][exc_bead+1]);
		tmp += pols[std::get<1>(pair)][exc_bead].sqdist(pols[std::get<0>(pair)][exc_bead+1]);
		tmp -= pols[std::get<0>(pair)][exc_bead].sqdist(pols[std::get<0>(pair)][exc_bead+1]);
		tmp -= pols[std::get<1>(pair)][exc_bead].sqdist(pols[std::get<1>(pair)][exc_bead+1]);
	}
	return 0.5*exc_const*tmp;
}

Force Graph::get_grad_weight(const std::vector<Polymer>& pols, bool use_positive, int bead, int part) const
{
	Force tmp(pols[0][0].size());
	if(positive == use_positive)
	{
		if(exchange_pairs.size()==0)
			return tmp;
		for(const auto& list : exchange_pairs)
			tmp += std::exp(exponent_sign*calc_exponent(pols,list))*grad_exponent(pols,list,bead,part);
		return mult*exc_const/exchange_pairs.size() * tmp;	
	}		
	return tmp;
}

Force Graph::grad_exponent(const std::vector<Polymer>& pols, const std::vector<std::pair<int,int>>& list, 
						   int bead, int part) const
{
	Force tmp(pols[0][0].size());
	int other_bead;
	int other_part;
	if(bead==exc_bead)
		other_bead = cyclic(exc_bead+1,pols[0].num_beads);
	else if(bead==cyclic(exc_bead+1,pols[0].num_beads))
		other_bead = exc_bead;
	else
		return tmp;
	auto it = std::find_if(list.begin(), list.end(), [&](const std::pair<int,int>& p){return p.first == part;});
	if(it == list.end())
	{	
		it = std::find_if(list.begin(), list.end(), [&](const std::pair<int,int>& p){return p.second == part;});
		if(it == list.end())
			return tmp;
		other_part = it->first;
	}
	else
		other_part = it->second;
	return pols[part][other_bead]-pols[other_part][other_bead];
}

double Graph::get_energy_diff(const std::vector<Polymer>& pols) const
{
	double tmp = 0;
	for(const auto& list : exchange_pairs)
		tmp += calc_exponent(pols,list);
	return -tmp*exponent_sign;
}

Force Graph::get_energy_diff_grad(const std::vector<Polymer>& pols, int bead, int part) const
{
	Force tmp(pols[0][0].size());
	/*int other_bead;
	if(bead==pols[0].num_beads-1)
		other_bead = 0;
	else if(bead==0)
		other_bead = pols[0].num_beads-1;*/
	//return (pols[part][other_bead] - pols[pols.size()-1-part][other_bead])*exc_const;
	
	for(const auto& list : exchange_pairs)
		tmp += grad_exponent(pols,list,bead,part);
	return tmp*exc_const*exponent_sign*(-1);
}

int Graph::get_id() const
{
	return id;
}
	
int Graph::get_mult() const
{
	return mult;
}

int Graph::get_sign() const
{
	if(positive)
		return 1;
	return -1;
}

std::vector<std::vector<int>> Graph::get_chains() const
{
	return chains;
}
