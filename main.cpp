#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>

#include "simulation.hpp"
#include "parameters.hpp"
#include "graph.hpp"

std::vector<Graph> generate_graphs(const Parameters& params)
{
	std::vector<Graph> graphs;
	if(params.sign==0)
		graphs.push_back(Graph(params,0));
	else
	{
		int num_graphs = params.num_parts; //Needs to be modified for N>3
		for(int id=0; id<num_graphs; ++id)
			graphs.push_back(Graph(params,id));
	}
	return graphs;
}


int main(int argc, char* argv[])
{
	bool continue_sim = false;
	if(argc > 1)
	{
		if((std::string(argv[1])=="-c")||(std::string(argv[1])=="--continue"))
		{
			continue_sim = true;
		}
	}
	//std::vector<double> betas = {0.067,0.134,0.2,0.335,0.67,1.005,1.34,1.675,2.01};
	//{0.5,0.7,1.0,1.5,2.0,3.0,4.0}
	//std::vector<double> taus = {0.067};
	//std::vector<double> betas = {1.0};
	//std::vector<double> gammas = {0.862, 0.9, 0.95, 1.0};
	//std::vector<double> taus = {2.0,1.0,0.5,0.3,0.2,0.15,0.1};
	//std::vector<double> taus = {1.0,0.5,0.3,0.2,0.1,0.067,0.05};//,0.067,0.05,0.04};
	//std::vector<double> taus = {0.067};
	Parameters params;
	try
	{
		params.read_file("configuration.cfg");
	}
	catch(const std::exception& e)
	{
		std::cout << "exception: " << e.what() << std::endl;
		return 0;
	}
	params.calculate_dependencies();
	std::ofstream results_file;
	if(!continue_sim)
	{
		results_file.open("results.dat");
		results_file << "%beta";
		for(int id : params.to_measure)
			results_file << "\tObsId " << id << "\t\tError\t";
		results_file << std::endl;
	}
	else
		results_file.open("results.dat", std::ios_base::app);
	results_file.precision(10);

	std::vector<Graph> graphs;
	
	if(0) //If you want to do several beta in the same simulation
	{
		std::vector<double> betas = {0.2,0.3,0.5,1.0,2.0,3.0,4.0};
		for(auto beta : betas)
		{
			params.beta = beta;
			params.calculate_dependencies();
			graphs = generate_graphs(params);
			std::cout << beta << ",,\t" << params.beta << std::endl;
			Simulation sim(params, results_file, continue_sim, graphs);
			sim.run();
		}
	}
	if(0) //If several tau in the same simulation
	{
		std::vector<double> taus = {1.0,0.3,0.2,0.15,0.1,0.067,0.05};
		for(auto tau : taus)
		{
			params.tau = tau;
			params.calculate_dependencies();
			graphs = generate_graphs(params);
			Simulation sim(params, results_file, continue_sim, graphs);
			sim.run();
		}
	}
	if(1) //Standard option: pick beta and tau from the configuration file
	{
		graphs = generate_graphs(params);
		Simulation sim(params, results_file, continue_sim, graphs);
		sim.run();
	}
	
	results_file.close();
	return 0;
}
 
 /*
  * Note that mass must be the same of all particles. Otherwise things have to be 
  * updated in the GLE and force calculation (curvature is different)
  * 
  */
