
#include <stdexcept>
#include <cmath>
#include "simulation.hpp"


//Simulation::Simulation() : dt_(0.1), dt_2m_(0.05){}

Simulation::Simulation(const Parameters& params, std::ofstream& _res_file, bool continue_sim,
					   const std::vector<Graph>& graphs_in)
	: dt_md(params.dt_md), dt_sample(params.dt_sample), num_parts(params.num_parts), interac(Interaction(params)),
	  length_scale(params.length_scale), num_blocks(params.num_blocks),
	  //steps_per_sample(params.steps_per_sample),
	  num_beads(params.num_beads), dim(params.dim),
	  num_bins(params.num_bins), num_bins_2d(params.num_bins_2d), hist_size(params.hist_size),
	  hist_size_1d(params.hist_size*2), hist_1d_min(-params.hist_size),
	  temperature(params.temperature), thermalization_time(params.thermalization_time),
	  thermostat_on(params.with_thermostat), res_file(_res_file),
	  non_sampling_time(params.non_sampling_time),
	  sampling_time(params.sampling_time), sampling_time_per_block(params.sampling_time/params.num_blocks),
	  beta(params.beta), sign(params.sign), exc_const(params.exc_const),
	  tau(params.tau), bias(Bias(params,continue_sim,graphs_in)),
	  graphs(graphs_in), spin(params.spin),
	  bias_update_time(params.bias_update_time), cont_sim(continue_sim),
	  allow_perm_switch(params.allow_perm_switch), permutation_trial_time(params.permutation_trial_time),
	  more_output(params.more_output), wigner_parameter(params.wigner_parameter)
{
	for(int n=0; n<num_parts; ++n)
		polymers.push_back(Polymer(params));

	//
	std::cout << "Graph\tMult\tSign\tChainlengths" << std::endl;
	for(const Graph& graph: graphs)
	{
		std::cout << graph.get_id() << "\t" << graph.get_mult() << "\t" << graph.get_sign() << "\t";
		const auto& chains = graph.get_chains();
		for(const std::vector<int>& chain : chains)
		{
			std::cout << chain.size() << " ";
			/*for(int pol : chain)
				std::cout << pol << ", ";
			std::cout << "; ";*/
		}
		std::cout << std::endl;
	}
	//
	current_graph_id = 0;
	if(polymers[0].connected)
		current_graph_id=1;
	
	finished = false;
	block = 0;
	time_sampled = 0;
	samples = 0;
	overall_time = 0;
	bar_width=70;
	progress=0;
	exchange_factor=1.0;
	exc_sum=exc_avg=exc_avg_sq=exc_sq=exc_sq_avg=e_s_sum=e_s_avg=e_s_avg_sq=0;
	sgn_sum=sgn_avg=sgn_avg_sq=0;
	std::vector<int> to_print = params.to_print_every_sample;
	for(int id : params.to_measure)
	{
		auto pair = obs.insert(std::pair<int,Observable>(id,Observable(id,params)));
		if(std::find(to_print.begin(), to_print.end(), id) != to_print.end())
			pair.first->second.set_print_on();		
	}
	histogram.assign(params.num_bins,0);
	histogram_avg.assign(params.num_bins,0);
	histogram_sq_avg.assign(params.num_bins,0);
	//int d = polymers[0][0].size();
	histogram_1d.assign(dim,histogram);
	histogram_1d_avg.assign(dim,histogram);
	histogram_1d_sq_avg.assign(dim,histogram);
	
	std::vector<double> tmp_2d;
	tmp_2d.assign(num_bins_2d,0);
	histogram_2d.assign(num_bins_2d, tmp_2d);
	histogram_2d_avg.assign(num_bins_2d, tmp_2d);
	histogram_2d_sq_avg.assign(num_bins_2d, tmp_2d);	
	pair_distr_2d.assign(num_bins_2d,tmp_2d);
	pair_distr_2d_avg.assign(num_bins_2d,tmp_2d);
	pair_distr_2d_sq_avg.assign(num_bins_2d,tmp_2d);
	pair_distr_1d_proj.assign(dim,tmp_2d);
	pair_distr_1d_proj_avg.assign(dim,tmp_2d);
	pair_distr_1d_proj_sq_avg.assign(dim,tmp_2d);
	
	cv_hist_num_bins = (int) std::round((cv_hist_max-cv_hist_min)/cv_hist_res)+1;
	cv_hist.assign(cv_hist_num_bins,0);
	weight_en_hist.assign(cv_hist_num_bins,0);
	exc_fac_hist.assign(cv_hist_num_bins,0);
	cv_hist_width = cv_hist_max-cv_hist_min;
	
	//hist_de_num_bins = round((hist_de_max-hist_de_min)/hist_de_resolution);
	//histogram_delta_e.assign(hist_de_num_bins,0);
	//hist_de_width=hist_de_max-hist_de_min;
	hist_c_num_bins=round((hist_c_max-hist_c_min)/hist_c_resolution) + 1;
	hist_c.assign(hist_c_num_bins,0);
	
	bias_update_counter = 0;
	permutation_switch_counter = 0;
	movie_start_time = 0;//non_sampling_time + sampling_time/2.0;
	movie_end_time = 50;// movie_start_time + dt_sample*10000;
	try 
	{
		gle = new GLE(polymers, dt_md, temperature, polymers[0].mass, polymers[0].num_beads, 
					num_parts, polymers[0][0].size(),thermostat_on);
	}
	catch(const std::exception& e)
	{
		std::cout << "exception: " << e.what() << std::endl;
		return;
	}
	
	mt = std::mt19937(rd());
	uni_distr = std::uniform_real_distribution<double>(0.0,1.0);
	int_distr = std::uniform_int_distribution<int>(0,graphs.size()-1);
	
	printed_warning = false;
}




void Simulation::setup() 
{
	timer.start();
	if(cont_sim)
	{
		read_input_coords();
		read_old_measurements();
		exc_file.open("exc_factor.dat", std::ios_base::app);
		cv_file.open("cv.dat", std::ios_base::app);
		rew_factor_file.open("rew_factor.dat", std::ios_base::app);
		vmd_file.open("vmd.xyz", std::ios_base::app);
		vmd_file2.open("vmd2.xyz", std::ios_base::app);
		file_fsum.open("fsum_N"+std::to_string(2-(int) polymers[0].connected)+"_P"+std::to_string(num_beads)+".dat", std::ios_base::app);
	}
	else
	{
		iteration_nbr = 0;
		initialize_coords_simple();
		exc_file.open("exc_factor.dat");
		cv_file.open("cv.dat");
		rew_factor_file.open("rew_factor.dat");
		vmd_file.open("vmd.xyz");
		vmd_file2.open("vmd2.xyz");
		file_fsum.open("fsum_N"+std::to_string(2-(int) polymers[0].connected)+"_P"+std::to_string(num_beads)+".dat");
		file_fsum << "0\t";
		for(int bin=0; bin<hist_c_num_bins; ++bin)
			file_fsum << hist_c_min + bin*hist_c_resolution << "\t\t";
		file_fsum << std::endl;
	}
	exc_file.precision(8);
	cv_file.precision(8);
	rew_factor_file.precision(8);
	vmd_file.precision(8);
	vmd_file2.precision(8);
	std::cout << iteration_nbr << std::endl;
	std::cout << "cumulated time = " << non_sampling_time + sampling_time << ", P = " << num_beads 
				<< ", spin = " << spin << ", dt = " << dt_md << ", conn = " << (int)polymers[0].connected
				<< ", perm switch: " << (int) allow_perm_switch << std::endl;
	std::time_t t = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
	logfile.open("logfile_P"+std::to_string(num_beads)+"_"+std::to_string(iteration_nbr));
	logfile << std::ctime(&t);
	logfile << "P=" << num_beads << " dim=" << dim << " dt=" << dt_md 
			<< " sampl_time=" << sampling_time
			<< " samples=" << std::round(sampling_time/dt_sample) 
			<< " T=" << temperature << " RW=" << wigner_parameter << std::endl << std::endl;
	logfile << "Block\tExc_const\tAvg_sign";
	for(auto& pair : obs)
		logfile << "\t" << pair.second.get_name();
	logfile << std::endl;
	logfile.precision(8);
	file_fsum.precision(8);
	if(!cont_sim)
	{	
		thermalize();
		run_wo_sampling();
	}
}

void Simulation::read_input_coords()
{
	std::ifstream coords_file("coords.xyz");
	std::ifstream vels_file("vels.xyz");
	std::string line,line2;
	getline(coords_file,line); getline(coords_file,line);
	getline(vels_file,line); getline(vels_file,line);
	for(int n=0; n<polymers.size(); ++n)
	{
		Polymer& pol = polymers[n];
		for(int bead=0; bead<pol.num_beads; ++bead)
		{
			getline(coords_file,line);
			std::istringstream iss(line);
			getline(vels_file,line2);
			std::istringstream iss2(line2);
			double number;
			for(int d=0; d<dim; ++d)
			{
				iss >> number;
				pol[bead][d] = number;
				iss2 >> number;
				pol.vels[bead][d] = number;
			}
		}
	}		
}

void Simulation::initialize_coords_simple()
{
	for(int n=0; n<polymers.size(); ++n)
	{
		Polymer& pol = polymers[n];
		for(int bead=0; bead<num_beads; ++bead)
		{
			for(int d=0; d<dim; ++d)
			{
				pol[bead][d] = length_scale * (1.5 + bead*0.01) * (1.0 + 0.5*n) * std::pow(-1,n);
				//pol[bead][d] = length_scale * ((double) bead/pol.num_beads + 1.5) * std::pow(-1,n);
				pol.vels[bead][d] = 0.0;
				//pol.forces[bead][d] = 0.0;
			}
		}
	}
}

void Simulation::read_old_measurements()
{
	std::ifstream infile("measurements.dat");
	std::string line;
	while(std::getline(infile, line))
	{
		std::istringstream iss(line);
		std::string name;
		int obs_id; double tmp1; double tmp2; double tmp3; double tmp4;
		int bin;
		double prev_blocks;
		if(iss >> name)
		{
			if(name=="iteration_nbr")
			{
				iss >> iteration_nbr;
				iteration_nbr++;
			}
			else if(name=="times")
			{
				iss >> non_sampling_time >> time_sampled;
				prev_blocks = time_sampled/sampling_time * num_blocks;
				block += prev_blocks;
				num_blocks += prev_blocks;
				sampling_time += time_sampled;
				overall_time += non_sampling_time+time_sampled;
				std::cout << "time in this run\t" << time_sampled << std::endl;
			}
			else if(name=="gamma")
			{
				iss >> exc_avg >> exc_avg_sq >> exc_sq_avg;
				std::cout << name << "\t" << exc_avg << "\t" << exc_avg_sq << "\t" << exc_sq_avg << std::endl;
			}
			else if(name=="bias_reweight")
			{
				iss >> tmp1 >> tmp2;
				bias.set_rew_factor_avg(tmp1,tmp2);
				std::cout << name << "\t" << tmp1 << "\t" << tmp2 << std::endl;
			}
			else if(name=="obs")
			{
				iss >> obs_id >> tmp1 >> tmp2 >> tmp3 >> tmp4; //avg, avg_sq, weighted_avg, weighted_avg_sq
				obs.at(obs_id).set_avgs(tmp1,tmp2,tmp3,tmp4,prev_blocks);
			}
			else if(name=="gaussian")
			{
				iss >> tmp1 >> tmp2; //height, center
				bias.add_gaussian(tmp1, tmp2);
			}
			else if(name=="prob_dist")
			{
				iss >> bin >> tmp1 >> tmp2;
				histogram_avg[bin] = tmp1;
				histogram_sq_avg[bin] = tmp2;
			}
			else
				std::cout << name << " is not a valid variable in measurement.dat" << std::endl;
		}
	}
	bias.restore_splines_transient(beta);
	//update_exc(false);
	bias.update_cv(polymers);
	infile.close();
}


void Simulation::thermalize()
{
	double time_thermalized = 0;
	while(time_thermalized<thermalization_time)
	{
		verlet_step();
		time_thermalized += dt_md;
	}
}

void Simulation::run_wo_sampling()
{
	double tmp2 = dt_md/(0.1*non_sampling_time);
	double time_before_sampl=0;
	while(time_before_sampl<non_sampling_time)
	{
		//update_exc(false);
		verlet_step();
		time_before_sampl += dt_md;
		bias_update_counter += dt_md;
		permutation_switch_counter += dt_md;
		overall_time += dt_md;
		update_exc(false);
		if(bias_update_counter >= bias_update_time)
		{
			bias.update_bias(polymers,beta,overall_time);
			//bias.update_cv(polymers,beta);
			bias_update_counter = 0;
		}
		if(permutation_switch_counter >= permutation_trial_time)
		{
			try_permutation_change();
			permutation_switch_counter = 0;
		}
		double tmp = time_before_sampl/(0.1*non_sampling_time);
		if(std::fmod(tmp,1.0)<tmp2)
			update_screen();

		/*if(more_output)
		{
			exc_file << overall_time << "\t" << exchange_factor << std::endl;
			cv_file << overall_time << "\t" << bias.get_cv() << "\t"
					<< bias.energy_diff(polymers) << "\t" << current_graph_id
					<< "\t" << polymers[0].bias_forces[0][0] << std::endl;
			rew_factor_file << overall_time << "\t" << bias.get_rew_factor() << std::endl;
		}*/
	}
}

void Simulation::run() 
{
	setup();
	while(!finished)
		run_block();
	stop();
}

void Simulation::run_block()
{
	//update_exc(false);
	double time_since_last_sample=0;
	while(time_sampled<sampling_time_per_block*(block+1))
	{
		while(time_since_last_sample<dt_sample) //(int step=0; step<steps_per_sample; ++step)
		{
			verlet_step();
			time_since_last_sample += dt_md;
		}
		time_since_last_sample = 0;
		++samples;
		time_sampled += dt_sample;
		overall_time += dt_sample;
		bias.update_cv_rew(polymers,beta);
		update_exc(true);
		measure();
		bias_update_counter += dt_sample;
		permutation_switch_counter += dt_sample;
		if(bias_update_counter >= bias_update_time) //check if remainder is 0
		{
			bias.update_bias(polymers,beta,overall_time);
			//bias.update_cv(polymers,exchange_factor);
			bias_update_counter = 0;
		}
		if(permutation_switch_counter >= permutation_trial_time)
		{
			try_permutation_change();
			permutation_switch_counter = 0;
		}
		if((polymers[0][0][0]!=polymers[0][0][0])&&(!printed_warning))
		{
			std::cout << overall_time << "\t" << "warning, coordinates are NaN!" << std::endl;
			printed_warning = true;
		}
	}
	update_avgs();
	++block;
	print_to_logfile();
	samples=0;
	bias.set_new_block();
	update_screen();
	if(block>=num_blocks)
		finished = true;
}

void Simulation::verlet_step()
{
	gle->run();
	for(Polymer& pol : polymers)
	{
		pol.update_vels();
		pol.move();
	}
	update_exc(false);
	bias.update_cv(polymers);
	interac.update_forces(polymers,bias);
	for(Polymer& pol: polymers)
		pol.update_vels();
	gle->run();
}


void Simulation::measure()
{
	for(auto& pair : obs)
		pair.second.measure(polymers,interac,overall_time,exchange_factor,bias.get_rew_factor(),graphs,current_graph_id);	
	update_histogram();
	if((!printed_warning)&&(bias.get_cv()!=bias.get_cv()))
	{
		std::cout << "time: " << overall_time << "\tWarning, CV is NaN!" << std::endl;
		printed_warning = true;
	}
	if(more_output)
	{
		exc_file << overall_time << "\t" << exchange_factor << std::endl;
		cv_file << overall_time << "\t" << bias.get_cv() << "\t"
				<< bias.energy_diff(polymers) << "\t" << current_graph_id
				<< "\t" << polymers[0].bias_forces[0][0] << std::endl;
					//bias.energy_diff(polymers) << //std::endl;
					// << std::endl;
		rew_factor_file << overall_time << "\t" << bias.get_rew_factor() << std::endl;
	}
	if(overall_time>movie_start_time && overall_time<movie_end_time)
		print_vmd();
}

void Simulation::update_avgs()
{
	double tmp = exc_sum/bias.get_rew_factor_block();
	//double tmp2 = e_s_sum/samples;
	exc_avg = (exc_avg*block + tmp) / (block+1.0);
	exc_avg_sq = (exc_avg_sq*block + tmp*tmp)/(block+1.0);
	exc_sq_avg = (exc_sq_avg*block + exc_sq)/(block+1.0);
	exc_sq = 0;
	exc_sum = 0;
	double tmp_sgn = sgn_sum/bias.get_rew_factor_block();
	sgn_avg = (sgn_avg*block + tmp_sgn)/(block+1.0);
	sgn_avg_sq = (sgn_avg_sq*block + tmp_sgn*tmp_sgn)/(block+1.0);
	sgn_sum = 0;
	//e_s_avg = (e_s_avg*block + tmp2)/(block+1.0);
	//e_s_avg_sq = (e_s_avg_sq*block + tmp2*tmp2)/(block+1.0);
	//e_s_sum = 0;
	for(auto& ob : obs)
	{
		ob.second.update_avg(samples,tmp);
		ob.second.set_zero();
	}
	/*double hist_norm=0;
	for(int bin=0; bin<num_bins; ++bin)
		hist_norm += histogram[bin];
	hist_norm *= hist_size/num_bins;*/
	for(int bin=0; bin<num_bins; ++bin)
	{
		histogram_avg[bin] = (histogram_avg[bin]*block + histogram[bin]/tmp)/(block+1.0);
		histogram_sq_avg[bin] = (histogram_sq_avg[bin]*block + std::pow(histogram[bin]/tmp,2))/(block+1.0);
		histogram[bin]=0;
		for(int d=0; d<polymers[0][0].size(); ++d)
		{
			histogram_1d_avg[d][bin] = (histogram_1d_avg[d][bin]*block + histogram_1d[d][bin]/tmp)/(block+1.0);
			histogram_1d_sq_avg[d][bin] = (histogram_1d_sq_avg[d][bin]*block + std::pow(histogram_1d[d][bin]/tmp,2))/(block+1.0);
			histogram_1d[d][bin]=0;
		}
	}
	for(int bin1=0; bin1<num_bins_2d; ++bin1)
	{
		for(int bin2=0; bin2<num_bins_2d; ++bin2)
		{
			histogram_2d_avg[bin1][bin2] = (histogram_2d_avg[bin1][bin2] + histogram_2d[bin1][bin2]/tmp)/(block+1.0);
			histogram_2d_sq_avg[bin1][bin2] = (histogram_2d_sq_avg[bin1][bin2] + std::pow(histogram_2d[bin1][bin2]/tmp,2))/(block+1.0);
			histogram_2d[bin1][bin2] = 0;
			pair_distr_2d_avg[bin1][bin2] = (pair_distr_2d_avg[bin1][bin2] + pair_distr_2d[bin1][bin2]/tmp)/(block+1.0);
			pair_distr_2d_sq_avg[bin1][bin2] = (pair_distr_2d_sq_avg[bin1][bin2] + std::pow(pair_distr_2d[bin1][bin2]/tmp,2))/(block+1.0);
			pair_distr_2d[bin1][bin2] = 0;
		}
		for(int d=0; d<dim; ++d)
		{
			pair_distr_1d_proj_avg[d][bin1] = (pair_distr_1d_proj_avg[d][bin1] + pair_distr_1d_proj[d][bin1]/tmp)/(block+1.0);
			pair_distr_1d_proj_sq_avg[d][bin1] = (pair_distr_1d_proj_sq_avg[d][bin1] + std::pow(pair_distr_1d_proj[d][bin1]/tmp,2))/(block+1.0);
			pair_distr_1d_proj[d][bin1] = 0;
		}
	}
	bias.update_rew_factor_avg(block);
}

void Simulation::update_histogram()
{
	double tmp=0;
	int bin;
	int bin1;
	int bin2;
	double rew_factor = bias.get_rew_factor();
	double cv = bias.get_cv();
	double weight = exchange_factor*rew_factor;
	for(int bead=0; bead<polymers[0].num_beads; ++bead)
	{
		//if(polymers[0].connected)
		//	tmp = polymers[0][bead].dist(polymers[0][bead+polymers[0].num_beads/2]);		
		//else
		tmp = polymers[0][bead].dist(polymers[1][bead]);
		bin = calc_bin(tmp,num_bins,hist_size);
		if(bin>=0 && bin<num_bins)
			histogram[bin] += weight;
		for(int n=0; n<polymers.size(); ++n)
		{
			for(int d=0; d<polymers[n][bead].size(); ++d)
			{
				tmp = polymers[n][bead][d];
				bin = calc_bin(tmp - hist_1d_min,num_bins,hist_size_1d);
				if(d==0)
					bin1 = calc_bin(tmp - hist_1d_min,num_bins_2d,hist_size_1d);
				if(d==1)
					bin2 = calc_bin(tmp - hist_1d_min,num_bins_2d,hist_size_1d);
				if(bin>=0 && bin<num_bins)
					histogram_1d[d][bin] += weight;
			}
			if(polymers[n][bead].size()==2)
			{
				if(bin1>=0 && bin2>= 0 && bin1<num_bins_2d && bin2<num_bins_2d)
					histogram_2d[bin1][bin2] += weight;
				int pc_bin1 = calc_bin(polymers[0][bead][0]-polymers[1][bead][0]-hist_1d_min,num_bins_2d, hist_size_1d);
				int pc_bin2 = calc_bin(polymers[0][bead][1]-polymers[1][bead][1]-hist_1d_min,num_bins_2d, hist_size_1d);
				if(pc_bin1>=0 && pc_bin2>=0 && pc_bin1<num_bins_2d && pc_bin2<num_bins_2d)
				{	
					pair_distr_2d[pc_bin1][pc_bin2] += weight;
					pair_distr_1d_proj[0][pc_bin1] += weight;
					pair_distr_1d_proj[1][pc_bin2] += weight;
				}
			}
		}
		/*bin = calc_bin(polymers[0][bead].dist0());
		if(bin>=0 && bin<num_bins)
			histogram_1p[0][bin] += 1;
		bin = calc_bin(polymers[1][bead].dist0());
		if(bin>=0 && bin<num_bins)
			histogram_1p[1][bin] += 1;*/
	}
	
	double cv_rel = cv-cv_hist_min;
	bin = calc_bin(cv_rel,cv_hist_num_bins,cv_hist_width);
	if((bin<cv_hist_num_bins)&&(bin>=0))
	{
		cv_hist[bin] += rew_factor;
		exc_fac_hist[bin] += weight;
		weight_en_hist[bin] += obs.at(2).get_last_value();
	}
	
	/*cv_rel = cv-hist_de_min;
	bin = calc_bin(cv_rel,hist_de_num_bins,hist_de_width);
	if((bin<hist_de_num_bins)&&(bin>=0))
		histogram_delta_e[bin] += rew_factor;*/
		
	double fd_argument;
	cv = bias.energy_diff(polymers);
	for(bin=0; bin<hist_c_num_bins; ++bin)
	{
		fd_argument = cv+(hist_c_min+bin*hist_c_resolution); //DeltaU + C
		if(polymers[0].connected)
			fd_argument *= (-1);
		hist_c[bin] += fermi_dirac(fd_argument) * bias.get_rew_factor();
	}
}

double Simulation::fermi_dirac(double x)
{
	return 1.0/(1.0+std::exp(x));
}

int Simulation::calc_bin(double dist, int nbins, double hsize)
{
	return round(dist*nbins/hsize);
}

void Simulation::update_screen()
{
	std::cout << "[";
	progress = (int) bar_width*overall_time/(non_sampling_time + sampling_time);
	for(int i=0; i<bar_width; ++i)
	{
		if(i<progress) std::cout << "=";
		else if(i==progress) std::cout << ">";
		else std::cout << " ";
	}
	std::cout << "]" << "\r";
	std::cout.flush();
}

void Simulation::print_to_logfile()
{
	logfile << block << "\t" << exc_avg;
	logfile << "\t" << sgn_avg;
	for(auto& ob : obs)
		logfile <<	"\t" << ob.second.get_weighted_avg(); 
	logfile << std::endl;
	double rew_factor_block = bias.get_rew_factor_block() / samples;
	file_fsum << block << "\t";
	//std::cout << rew_factor_block << "\t" << samples;
	for(int bin = 0; bin<hist_c_num_bins; ++bin)
	{
		file_fsum << hist_c[bin]/rew_factor_block << "\t";
		//std::cout << hist_c[bin] << "\t";
		hist_c[bin] = 0;
	}
	file_fsum << std::endl;
	//std::cout << std::endl; //
	//print C histogram, one row for each block
}

void Simulation::stop()
{
	print_config();
	logfile << "Name\t\tMeanValue\tErrorOfMean\tStandDev" << std::endl;
	//double rew_norm = bias.get_rew_factor_avg();
	logfile << "Exc_factor\t" << exc_avg << "\t" << simple_uncertainty(exc_avg,exc_avg_sq) 
			<< "\t" << std::sqrt(exc_sq_avg-exc_avg*exc_avg) << std::endl; // sqrt(n/(n-1)) unnecessary for large n
	logfile << "Avg_sign\t" << sgn_avg << "\t" << simple_uncertainty(sgn_avg,sgn_avg_sq) << std::endl;
	//logfile << "Exp_en_diff\t" << e_s_avg/rew_norm << "\t" << simple_uncertainty(e_s_avg,e_s_avg_sq)/rew_norm
	//		<< std::endl;
	res_file << beta;//polymers[0].num_beads;
	//res_file << 2-int(polymers[0].connected);//polymers[0].num_beads;//sampling_time;
	for(const auto& pair : obs)
	{
		const auto& ob = pair.second;
		//double avg = ob.get_avg();
		//double avg_sq = ob.get_avg_sq();
		double w_avg = ob.get_weighted_avg();
		double w_avg_sq = ob.get_weighted_avg_sq();
		logfile << ob.get_name() << "\t" << w_avg << "\t" << simple_uncertainty(w_avg,w_avg_sq) << "\t"
				<< std::endl;
				//<< weighted_uncertainty(w_avg,w_avg_sq) << std::endl;
		res_file << "\t" << w_avg << "\t" << simple_uncertainty(w_avg,w_avg_sq);
	}
	res_file << std::endl;
	
	std::ofstream pair_file("Pair_correlation.dat");
	std::ofstream hist_file_1d("Prob_dist1d.dat");

	for(int bin = 0; bin<num_bins; ++bin)
	{
		pair_file << hist_size*((double) bin/num_bins) << "\t" << histogram_avg[bin]
					<< "\t" << simple_uncertainty(histogram_avg[bin],histogram_sq_avg[bin]) << std::endl;
		hist_file_1d << hist_size_1d*((double) bin/num_bins) + hist_1d_min;
		for(int d=0; d<polymers[0][0].size(); ++d)
			hist_file_1d << "\t" << histogram_1d_avg[d][bin] << "\t"
						 << simple_uncertainty(histogram_1d_avg[d][bin],histogram_1d_sq_avg[d][bin]);
		hist_file_1d << std::endl;
	}
	pair_file.close();
	hist_file_1d.close();
	
	if(polymers[0][0].size()==2)
	{
		std::ofstream hist_file_2d("Prob_dist2d.dat");
		std::ofstream hist_file_2derr("Prob_dist2d_err.dat");
		hist_file_2d << "-1";
		hist_file_2derr << "-1";
		for(int bin1=0; bin1<num_bins_2d; ++bin1)
		{
			hist_file_2d << "\t" << hist_size_1d*((double) bin1/num_bins_2d) + hist_1d_min;
			hist_file_2derr << "\t" << hist_size_1d*((double) bin1/num_bins_2d) + hist_1d_min;		
		}
		hist_file_2d << std::endl;
		hist_file_2derr << std::endl;
		for(int bin2=0; bin2<num_bins_2d; ++bin2)
		{
			hist_file_2d << hist_size_1d*((double) bin2/num_bins_2d) + hist_1d_min;
			hist_file_2derr << hist_size_1d*((double) bin2/num_bins_2d) + hist_1d_min;
			for(int bin1=0; bin1<num_bins_2d; ++bin1)
			{
				hist_file_2d << "\t" << histogram_2d_avg[bin1][bin2];
				hist_file_2derr << "\t" << simple_uncertainty(histogram_2d_avg[bin1][bin2],histogram_2d_sq_avg[bin1][bin2]);
			}
			hist_file_2d << std::endl;
			hist_file_2derr << std::endl;
		} 
		hist_file_2d.close();
		hist_file_2derr.close();
		
		std::ofstream pair_corr_2d("Pair_corr2d.dat");
		std::ofstream hist_file_2d_other_format("Prob_dist2d_for_gnu.dat");
		std::ofstream pair_corr_1d("Pair_corr1d_proj.dat");
		for(int bin1=0; bin1<num_bins_2d; ++bin1)
		{
			double r1 = hist_size_1d*((double) bin1/num_bins_2d) + hist_1d_min;
			for(int bin2=0; bin2<num_bins_2d; ++bin2)
			{
				double r2 = hist_size_1d*((double) bin2/num_bins_2d) + hist_1d_min;
				pair_corr_2d << r1 << "\t" << r2 << "\t" << pair_distr_2d_avg[bin1][bin2] << "\t" 
							 << simple_uncertainty(pair_distr_2d_avg[bin1][bin2],pair_distr_2d_sq_avg[bin1][bin2]) << std::endl;
				hist_file_2d_other_format << r1 << "\t" << r2 << "\t" << histogram_2d_avg[bin1][bin2] << "\t" 
							 << simple_uncertainty(histogram_2d_avg[bin1][bin2],histogram_2d_sq_avg[bin1][bin2]) << std::endl;
			}
			pair_corr_1d << r1 << "\t" << pair_distr_1d_proj_avg[0][bin1] << "\t" << simple_uncertainty(pair_distr_1d_proj_avg[0][bin1],pair_distr_1d_proj_sq_avg[0][bin1])
						 << "\t" << pair_distr_1d_proj_avg[1][bin1] << "\t" << simple_uncertainty(pair_distr_1d_proj_avg[1][bin1],pair_distr_1d_proj_sq_avg[1][bin1]) << std::endl;
			
		}
		pair_corr_2d.close();
		hist_file_2d_other_format.close();
		pair_corr_1d.close();
	}
	
	std::ofstream cv_hist_file("CV_distributions.dat");
	//std::ofstream weighted_en_file("Energy_distr.dat");
	double rew_factor_avg = bias.get_rew_factor_avg();
	double norm=0;
	for(int bin=0; bin<cv_hist_num_bins; ++bin)
		norm += cv_hist[bin];
	norm *= cv_hist_res;
	for(int bin=0; bin<cv_hist_num_bins; ++bin)
	{
		double s = cv_hist_width*((double) bin/cv_hist_num_bins) + cv_hist_min;
		//double bde_entry = bde_hist[bin]/norm;
		//double exc_fac_entry = exc_fac_hist[bin]/norm;
		cv_hist_file << s << "\t" << cv_hist[bin]/norm << "\t" << exc_fac_hist[bin]/norm 
					  << "\t" << weight_en_hist[bin]/norm 
					  << "\t" << bias.calc_bias(s) << "\t" << bias.calc_bias_der(s) << std::endl;
		/*	bde_hist_file << std::exp(s) + sign << std::endl;
		else
			bde_hist_file << (1 + sign*std::exp(-s))*hist_entry << std::endl;*/
		//weighted_en_file << s << "\t" << weight_en_hist[bin]/(rew_factor_avg*bde_hist_res) << std::endl;
	}
	cv_hist_file.close();
	//weighted_en_file.close();
	
	//std::string particles;
	/*if(polymers[0].connected)
		particles = "1";
	else
		particles = "2";
	std::ofstream de_file("DeltaE_hist_N"+particles+".dat");
	de_file.precision(10);
	double normalization=0;
	for(int bin=0; bin<hist_de_num_bins; ++bin)
		normalization += histogram_delta_e[bin];
	normalization *= hist_cv_resolution;
	for(int bin = 0; bin<hist_de_num_bins; ++bin)
		de_file << hist_de_width*((double) bin/hist_de_num_bins) + hist_de_min 
				<< "\t" << histogram_delta_e[bin]/normalization << std::endl;*/
					
	delete gle;
	
	timer.stop();
	logfile << "Finished in " << timer.duration() << " s" << std::endl;
	logfile.close();

	std::cout << std::endl << timer.duration() << " s" << std::endl;
}

void Simulation::print_config()
{
	std::ofstream coord_file("coords.xyz");
	std::ofstream vel_file("vels.xyz");
	coord_file << num_parts*polymers[0].num_beads << std::endl << std::endl;
	vel_file << num_parts*polymers[0].num_beads << std::endl << std::endl;
	for(int n=0; n<polymers.size(); ++n)
	{
		for(int bead=0; bead<polymers[n].num_beads; ++bead)
		{
			Point& p = polymers[n][bead];
			Point& v = polymers[n].vels[bead];
			for(int d=0; d<p.size(); ++d)
			{
				coord_file << p[d] << "\t";
				vel_file << v[d] << "\t";
			}
			coord_file << std::endl;
			vel_file << std::endl;
		}
	}
	coord_file.close();
	vel_file.close();
	std::ofstream outfile("measurements.dat");
	outfile.precision(10);
	outfile << "iteration_nbr\t" << iteration_nbr << std::endl;
	outfile << "times\t" << non_sampling_time << "\t" << time_sampled << std::endl;
	outfile << "gamma\t" << exc_avg << "\t" << exc_avg_sq << "\t" << exc_sq_avg << std::endl;
	outfile << "bias_reweight\t" << bias.get_rew_factor_avg() << "\t" << bias.get_count() << std::endl;
	for(auto& pair : obs)
		outfile << "obs\t" << pair.second.get_id() << "\t" << pair.second.get_avg()
				<< "\t" << pair.second.get_avg_sq() << "\t" << pair.second.get_weighted_avg()
				<< "\t" << pair.second.get_weighted_avg_sq() << std::endl;
	std::vector<double> heights = bias.get_heights();
	std::vector<double> centers = bias.get_centers();
	for(int i=0; i<heights.size(); ++i)
		outfile << "gaussian\t" << heights[i] << "\t" << centers[i] << std::endl;
	for(int bin=0; bin<num_bins; ++bin)
		outfile << "prob_dist\t" << bin << "\t" << histogram_avg[bin] 
				<< "\t" << histogram_sq_avg[bin] << std::endl;
}

void Simulation::print_vmd()
{
	vmd_file << polymers[0].num_beads << std::endl << std::endl;
	if(!polymers[0].connected)
		vmd_file2 << polymers[1].num_beads << std::endl << std::endl;
	int dim = polymers[0][0].size();
	for(int bead=0; bead<polymers[0].num_beads; ++bead)
	{
		vmd_file << "O";
		vmd_file2 << "N";
		for(int d=0; d<dim; ++d)
		{
			vmd_file << "\t" << polymers[0][bead][d];
			if(!polymers[0].connected)
			vmd_file2 << "\t" << polymers[1][bead][d];
		}
		if(dim==2)
		{
			vmd_file << "\t" << 0;
			vmd_file2<< "\t" << 0;
		}
		if(dim==1)
		{
			vmd_file << "\t" << 0 << "\t" << 0;
			vmd_file2<< "\t" << 0 << "\t" << 0;
		}
		vmd_file << std::endl;
		vmd_file2<< std::endl;
	}
}

double Simulation::exc_exponent(int bead) const
{
	return exc_const * (polymers[0][bead]-polymers[1][bead])*(polymers[0][bead+1]-polymers[1][bead+1]);
}

void Simulation::update_exc(bool count_to_average)
{
	if(sign==0)
	{
		exchange_factor = 1.0;
	}
	//double tmp = 0;
	/*if(polymers[0].connected)
	{
		for(int bead=0; bead<num_beads; ++bead)
			tmp += std::exp(exc_exponent(bead));
		e_s = tmp/num_beads;
		exchange_factor = e_s + sign;
	}
	else
	{
		for(int bead=0; bead<num_beads; ++bead)
			tmp += std::exp(-exc_exponent(bead));
		//e_s = std::exp(-exc_exponent(polymers[0].num_beads-1));
		e_s = tmp/num_beads;
		exchange_factor = 1.0 + sign*e_s;
	}*/
	pos_weight = 0;
	neg_weight = 0;
	for(const Graph& graph : graphs)
	{
		pos_weight += graph.get_weight(polymers,graphs[current_graph_id],true);
		neg_weight += graph.get_weight(polymers,graphs[current_graph_id],false);
		//std::cout << graph.get_weight(polymers,graphs[current_graph_id],true)<< "\t";
	}
	//std::cout << std::endl;
	exchange_factor = pos_weight - neg_weight;
	sgn = (pos_weight - neg_weight)/(pos_weight+neg_weight);
	//e_s = std::abs(exchange_factor-1);
	//std::cout << pos_weight << "\t" << neg_weight << std::endl; // 
	if(count_to_average)
	{
		double exc_weighted = exchange_factor * bias.get_rew_factor();
		//double e_s_weighted = e_s * bias.get_rew_factor();
		exc_sum += exc_weighted;
		exc_sq = (exc_sq*samples + exc_weighted*exc_weighted)/(samples+1.0);
		//e_s_sum += e_s_weighted;
		sgn_sum += sgn*bias.get_rew_factor();
	}
	bias.set_weights(pos_weight,neg_weight);
}

double Simulation::simple_uncertainty(double avg, double avg_sq) const
{
	if(block>=2)
		return std::sqrt((avg_sq - avg*avg)/(block));
	return 0;
}

void Simulation::try_permutation_change()
{
	if(allow_perm_switch)
	{
		int graph_id_to_try=-1;
		do
		{
			graph_id_to_try = int_distr(mt);
			//std::cout << graph_id_to_try << " ";
		}	
		while(graph_id_to_try == current_graph_id);	
		//std::cout << std::endl;
		double en_diff = graphs[graph_id_to_try].energy_diff(polymers,graphs[current_graph_id]);
		if((en_diff<0)||(uni_distr(mt) < std::exp(-beta*en_diff)))
		{
			//std::cout << "switched from diagram " << current_graph_id<< " to diagram " << graph_id_to_try << std::endl;
			current_graph_id = graph_id_to_try;
			if(current_graph_id==1)
				polymers[0].connected=true;
			if(current_graph_id==0)
				polymers[0].connected=false;
			//for(int i=0; i<100; ++i)
			//	verlet_step();
		}
	}
}


/*
double Simulation::weighted_uncertainty(double avg, double avg_sq) const
{
	//standard error of the mean is std/sqrt(blocks)
	if(block>=2)
	{
		//avg*= exc_avg; //not weighted w.r.t. Gamma 
		//avg_sq *= exc_avg*exc_avg; //not weighted w.r.t. Gamma
		//exc_avg /= rew_avg;
		//exc_avg_sq /= (rew_avg*rew_avg);
		double num_error = simple_uncertainty(avg,avg_sq);
		double den_error = simple_uncertainty(exc_avg,exc_avg_sq);
		return std::abs(avg/exc_avg) * (std::abs(num_error/avg) + std::abs(den_error/exc_avg));
		//error = num/den * (num_err/num + den_err/den)
		//relative errors do not need to be normalized by rew_factor_avg of bias
	}
	return 0;
}*/
