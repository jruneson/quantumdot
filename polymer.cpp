
#include "polymer.hpp"


Polymer::Polymer(const Parameters& params) : num_beads(params.num_beads) 
{
	for(int bead=0; bead<num_beads; ++bead)
	{
		coords.push_back(Point(params));
		vels.push_back(Point(params));
		forces.push_back(Force(params));
		dt_2m = params.dt / (2.0 * params.mass);
		mdt = params.mass * params.dt;
	}
}

const Point& Polymer::operator[](int bead) const
{
	bead = (bead + num_beads) % num_beads;
	return coords[bead];
}

Point& Polymer::operator[](int bead)
{
	bead = (bead + num_beads) % num_beads;
	return coords[bead];
}


void Polymer::move()
{
	for(int bead=0; bead<num_beads; ++bead)
		coords[bead] += vels[bead] * mdt;
}

void Polymer::update_vels()
{
	for(int bead=0; bead<num_beads; ++bead)
		vels[bead] += forces[bead] * dt_2m;
}
