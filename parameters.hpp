
#ifndef PARAMS_HPP
#define PARAMS_HPP

class Parameters{
public:
	//read from file
	const int num_parts=2;
	const int num_beads=30;
	const int dim=2;
	const double dt=0.1;
	
	const double mass = 1;
	const double curvature = 1;
	const double charge = 1;
	const double diel_const = 1;
};

#endif
