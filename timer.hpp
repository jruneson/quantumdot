
#include <chrono>
#ifndef TIMER_HPP
#define TIMER_HPP

struct Timer
{
	using time_point_t = std::chrono::high_resolution_clock::time_point;
	time_point_t tstart, tend;
	void start() {tstart = std::chrono::high_resolution_clock::now();}
	void stop() {tend = std::chrono::high_resolution_clock::now();}
	double duration() const {return std::chrono::duration<double>(tend-tstart).count();}
};

#endif
