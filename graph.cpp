
#include "graph.hpp"


Graph::Graph(const Parameters& params, int graph_id) //Chain lengths given as number of chains with length 1, followed by 2, 3 etc.
 : exc_const(params.exc_const), id(graph_id)
{
	int num_parts = params.num_parts;
	int spin = params.total_spin;
	exc_bead = params.num_beads-1;
	mult = 0;
	positive = true;
	if(num_parts==1)
		mult = 1;
	if(num_parts==2)
	{
		//std::vector<int> graph1 = {2,0};
		//std::vector<int> graph2 = {0,1};
		mult = 1;
		/*if(graph_id==1) //{2,0}
			//do nothing
		*/
		if(graph_id==2) // {0,1}
		{
			std::vector<std::pair<int,int>> tmp;
			tmp.push_back(std::pair<int,int>(0,1));
			exchange_pairs.push_back(tmp);
			if(spin==2)
				positive = false;
		}
	}
	if(num_parts==3)
	{
		/*std::vector<int> graph1 = {3,0,0};
		std::vector<int> graph2 = {1,1,0};				
		std::vector<int> graph3 = {0,0,1};*/
		if(graph_id==1) //{3,0,0}
			mult = 1;
		if(graph_id==2) //{1,1,0}
		{
			for(int n=0; n<num_parts; ++n) 
			{
				std::vector<std::pair<int,int>> tmp;
				tmp.push_back(std::pair<int,int>(n,cyclic(n+1,num_parts)));
				exchange_pairs.push_back(tmp);
			}
			if(spin==3)
			{
				mult = 3;
				positive = false;
			}
		}
		if(graph_id==3) //{0,0,1}
		{
			std::vector<std::pair<int,int>> tmp;
			for(int n=0; n<num_parts; ++n)
				tmp.push_back(std::pair<int,int>(n,cyclic(n+1,num_parts)));
			exchange_pairs.push_back(tmp);
			if(spin==1)
			{
				mult = 1;
				positive = false;
			}
			if(spin==3)
				mult = 2;
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
	{
		if(exchange_pairs.size()==0)
			return 1;
		double tmp = 0;
		for(const auto& list : exchange_pairs)
		{
			double tmp2 = 0;
			for(const auto& pair : list)
			{
				tmp2 += pols[std::get<0>(pair)][exc_bead].sqdist(pols[std::get<1>(pair)][exc_bead+1]);
				tmp2 -= pols[std::get<0>(pair)][exc_bead].sqdist(pols[std::get<0>(pair)][exc_bead+1]);
			}
			tmp += std::exp(-0.5*exc_const*tmp2);
		}
		//std::cout << mult*tmp << "\t" << exchange_pairs.size() << std::endl;

		return mult*tmp/exchange_pairs.size();
	}
	return 0;
}

Force Graph::get_grad_weight(const std::vector<Polymer>& pols, bool use_positive, int bead, int part) const
{
	Force tmp(pols[0][0].size());
	
	return mult*exc_const/exchange_pairs.size() * tmp;
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
