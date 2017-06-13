#include <algorithm>
#include "graph.hpp"


Graph::Graph(const Parameters& params, int graph_id) //Chain lengths given as number of chains with length 1, followed by 2, 3 etc.
 : exc_const(params.exc_const), id(graph_id), spin_proj(params.spin_proj),
	num_parts(params.num_parts)
{
	int spin = params.spin;
	exc_bead = params.num_beads-1;
	exc_bead2 = 0;
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
		if(graph_id==0) //{2,0}
		{
			for(int i=0; i<=1; ++i)
			{
				std::vector<int> tmp_ch;
				tmp_ch.push_back(i);
				chains.push_back(tmp_ch);
				exchange_pairs.push_back(std::pair<int,int>(i,i));
			} 
			//chains.push_back(std::vector<int>(0));
			//chains.push_back(std::vector<int>(1));
		}
		if(graph_id==1) // {0,1}
		{
			std::vector<int> tmp_ch = {0,1};
			chains.push_back(tmp_ch);
			
			//replace by pair calculation method
			//std::vector<std::pair<int,int>> tmp;
			exchange_pairs.push_back(std::pair<int,int>(1,1));
			exchange_pairs.push_back(std::pair<int,int>(0,0));
			
			if(spin==2)
				positive = false;
		}
	}
	if(num_parts==3)
	{
		if(graph_id==0) //{3,0,0}
		{
			mult = 1;
			for(int i=0; i<3; ++i)
			{
				std::vector<int> tmp_ch = {i};
				chains.push_back(tmp_ch);
				exchange_pairs.push_back(std::pair<int,int>(i,i));
			}
		}
		if(graph_id==1) //{1,1,0}
		{
			std::vector<int> tmp_ch = {0};
			chains.push_back(tmp_ch);
			tmp_ch = {1,2};
			chains.push_back(tmp_ch);
			//std::vector<std::pair<int,int>> tmp;
			exchange_pairs.push_back(std::pair<int,int>(0,0));
			exchange_pairs.push_back(std::pair<int,int>(2,2));//cyclic(n+1,num_parts)));
			exchange_pairs.push_back(std::pair<int,int>(1,1));
			if(spin==1)
				if(spin_proj)
				{
					mult = 1;
					positive = false;
				}
			if(spin==3)
			{
				//both for spin_proj true and false
				mult = 3;
				positive = false;
			}
		}
		if(graph_id==2) //{0,0,1}
		{
			std::vector<int> tmp_ch = {0,1,2};
			chains.push_back(tmp_ch);
			//std::vector<std::pair<int,int>> tmp;
			//for(int n=0; n<num_parts; ++n)
				//tmp.push_back(std::pair<int,int>(n,cyclic(n+1,num_parts)));
			exchange_pairs.push_back(std::pair<int,int>(2,1));
			exchange_pairs.push_back(std::pair<int,int>(0,2));
			exchange_pairs.push_back(std::pair<int,int>(1,0));	
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
				mult = 2; //for both with and without spin_proj
			}
		}		
	}
}

int Graph::cyclic(int i, int n) const
{
	return i%n;
}

double Graph::get_weight(const std::vector<Polymer>& pols,const Graph& current_graph, bool use_positive) const
{
	if(positive == use_positive)
		return weight(pols,current_graph);
	return 0;
}

double Graph::get_weight_signed(const std::vector<Polymer>& pols,const Graph& current_graph) const
{
	return get_sign()*weight(pols,current_graph);
}

double Graph::weight(const std::vector<Polymer>& pols,const Graph& current_graph) const
{
	//if(exchange_pairs.size()==0)
	//	return 1;
	//double tmp = 0;
	//for(const auto& list : exchange_pairs)
	double tmp = std::exp(-calc_exponent(pols,current_graph));
	return mult*tmp;
}

double Graph::calc_exponent(const std::vector<Polymer>& pols,const Graph& current_graph) const
{
	double tmp = 0;
	//for(const auto& pair : exchange_pairs)
	for(int n=0; n<num_parts; ++n)
	{
		tmp += single_spring(pols,n,std::get<1>(exchange_pairs[n]));
		tmp += single_spring(pols,std::get<0>(exchange_pairs[n]),n);
		tmp -= single_spring(pols,n,std::get<1>(current_graph.exchange_pairs[n]));
		tmp -= single_spring(pols,std::get<0>(current_graph.exchange_pairs[n]),n);
	}
	return 0.25*exc_const*tmp; //extra 0.5 to avoid double counting
}

double Graph::single_spring(const std::vector<Polymer>& pols,int n, int m) const
{
	return pols[n][exc_bead].sqdist(pols[m][exc_bead2]);
}

Force Graph::get_grad_weight(const std::vector<Polymer>& pols, const Graph& current_graph, 
							 bool use_positive, int bead, int part) const
{
	Force tmp(pols[0][0].size());
	if(positive == use_positive)
	{
		//if(exchange_pairs.size()==0)
		//	return tmp;
		//for(const auto& list : exchange_pairs)
		tmp = std::exp(-calc_exponent(pols,current_graph))*grad_exponent(pols,current_graph,bead,part);
		return -mult * tmp;	
	}		
	return tmp;
}

Force Graph::grad_exponent(const std::vector<Polymer>& pols, const Graph& current_graph,
						   int bead, int part) const
{
	Force tmp(pols[0][0].size());
	int other_bead;
	if(bead==exc_bead)
	{
		other_bead = exc_bead2;
		tmp += pols[std::get<1>(current_graph.exchange_pairs[part])][other_bead];
		tmp -= pols[std::get<1>(exchange_pairs[part])][other_bead];
	}
	else if(bead==exc_bead2)
	{
		other_bead = exc_bead;
		tmp += pols[std::get<0>(current_graph.exchange_pairs[part])][other_bead];
		tmp -= pols[std::get<0>(exchange_pairs[part])][other_bead];
	}
	else
		return tmp;
	//for(int n=0; n<num_parts; ++n)
	//{
	/*	tmp += pols[part][exc_bead]-pols[std::get<1>(exchange_pairs[part])][exc_bead2];
		tmp += pols[std::get<0>(exchange_pairs[part])][exc_bead] - pols[part][exc_bead2];
		tmp -= pols[part][exc_bead]-pols[std::get<1>(current_graph.exchange_pairs[part])][exc_bead2];
		tmp -= pols[std::get<0>(current_graph.exchange_pairs[part])][exc_bead]-pols[part][exc_bead2];*/
	//}
	return tmp*exc_const;
}

double Graph::energy_diff(const std::vector<Polymer>& pols, const Graph& graph) const
{
	return calc_exponent(pols,graph);
}

Force Graph::energy_diff_grad(const std::vector<Polymer>& pols, const Graph& graph,
							   int bead, int part) const
{
	return grad_exponent(pols,graph,bead,part);
}


/*Force Graph::grad_exponent(const std::vector<Polymer>& pols, const Graph& current_graph, 
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
	
	auto it = std::find_if(exchange_pairs.begin(), exchange_pairs.end(), [&](const std::pair<int,int>& p){return p.first == part;});
	if(it == exchange_pairs.end())
	{	
		it = std::find_if(list.begin(), list.end(), [&](const std::pair<int,int>& p){return p.second == part;});
		if(it == list.end())
			return tmp;
		other_part = it->first;
	}
	else
		other_part = it->second;
	return pols[part][other_bead]-pols[other_part][other_bead];
}*/

/*
double Graph::get_energy_diff(const std::vector<Polymer>& pols,const Graph&) const
{
	double tmp = 0;
	//for(const auto& list : exchange_pairs)
	tmp += calc_exponent(pols,current_graph);
	return tmp;
}*/

/*Force Graph::get_energy_diff_grad(const std::vector<Polymer>& pols, int bead, int part) const
{
	Force tmp(pols[0][0].size());*/
	/*int other_bead;
	if(bead==pols[0].num_beads-1)
		other_bead = 0;
	else if(bead==0)
		other_bead = pols[0].num_beads-1;*/
	//return (pols[part][other_bead] - pols[pols.size()-1-part][other_bead])*exc_const;
	
	/*for(const auto& list : exchange_pairs)
		tmp += grad_exponent(pols,list,bead,part);
	return tmp*exc_const;
}*/

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
