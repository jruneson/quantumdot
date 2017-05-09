
#ifndef SPLINE_HPP
#define SPLINE_HPP

#include <vector>
#include <array>

class Spline{
public:
	Spline();
	explicit Spline(double);

	
	void create_spline(const std::vector<double>&, const std::vector<double>&);
	//void delete_spline();
	void build_spline();
	double eval_spline(double) const;
	double get_min() const;
	double get_max() const;
	double get_step() const;
	
private:
	std::array<std::vector<double>,4> f_; //spline matrix with 4*npoints entries. The [0][i] entries are the values.
	int npoints_;
	double min_;
	double max_;
	bool created;
	double step_;
	
	
};

#endif
