
#include "polymer.hpp"


Polymer::Polymer(const Parameters& params) : num_beads(params.num_beads),
		dt(params.dt_md), dt_2m(params.dt_2m), mass(params.mass), connected(params.connected)
{
	for(int bead=0; bead<num_beads; ++bead)
	{
		coords.push_back(Point(params));
		vels.push_back(Point(params));
		ext_forces.push_back(Force(params));
		twopart_forces.push_back(Force(params));
		bias_forces.push_back(Force(params));
		spring_forces.push_back(Force(params));
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
		coords[bead] += vels[bead] * dt;
}

void Polymer::update_vels()
{
	for(int bead=0; bead<num_beads; ++bead)
		vels[bead] += (ext_forces[bead]+spring_forces[bead]+twopart_forces[bead]+bias_forces[bead]) * dt_2m;
}

Point Polymer::mean() const
{
	Point mean_point(coords[0].size());
	for(int bead=0; bead<num_beads; ++bead)
		mean_point += coords[bead];
	mean_point *= 1.0/num_beads;
	return mean_point;
}

Force Polymer::get_potential_force(int bead) const //forces without bias
{
	return ext_forces[bead]+twopart_forces[bead];
}
