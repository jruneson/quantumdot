
#ifndef OBSERVABLE_HPP
#define OBSERVABLE_HPP

class Observable {
public:
	Observable();
	//Observable(Parameters params);
	
	void set_zero();
	void update_avg(const int&);
	void normalize_avg(const int&);
	double get_avg() const;
	double std_dev(const int&) const;
	void operator+=(double);

private:
	double value;
	double avg;
	double avg_sq;
	
};

#endif
