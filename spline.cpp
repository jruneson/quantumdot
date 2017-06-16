
#include "spline.hpp"
#include <vector>
#include <algorithm>
#include <iostream>
#include <cmath>

/*Spline::Spline()
{
	created = false;
}*/

Spline::Spline(double step, double ref_point) : step_(step), ref_point_(ref_point_) 
{
	created = false;
}

void Spline::create_spline(const std::vector<double>& xs, const std::vector<double>& ys) //bool selects between force and pot
{
	//assume that xs is sorted from smallest to largest
	if(xs.size()!=ys.size())
		std::cout << "Spline: inputs must be of the same length" << std::endl;
	min_ = xs[0];
	max_ = xs.back();
	npoints_ = xs.size();
	f_[0] = ys;
	f_[1].assign(npoints_,0); //the last three lines are just to make size right, could be done differently
	f_[2].assign(npoints_,0);
	f_[3].assign(npoints_,0);
	build_spline();
}

/*void Spline::update_spline(const std::vector<double>& xs, const std::vector<double>& ys)
{
	if(xs.size()!=ys.size())
		std::cout << "Spline: inputs must be of the same length" << std::endl;
	if(xs[0] < min_)
		min_=xs[0];
	if(xs.back() > max_)
		max_ = xs.back();
	
}*/


/*void buildVSpline()
{
    double tw=1./12;
    double th=2./3;
    for(int i=2;i<Npoints-2;i++)
        f[i][1]=tw*f[i-2][0]-th*f[i-1][0]+th*f[i+1][0]-tw*f[i+2][0];
   
    f[1][1]=0.5*(f[2][0]-f[0][0]);
    f[Npoints-2][1]=0.5*(f[Npoints-1][0]-f[Npoints-3][0]);
    f[0][1]=2*f[1][1]-f[2][1];
    f[Npoints-1][1]=0;
    for(int i=0;i<Npoints-1;i++)
    {
        f[i][2]=3*(f[i+1][0]-f[i][0])-2*f[i][1]-f[i+1][1];
        f[i][3]=f[i+1][1]+f[i][1]+2*(f[i][0]-f[i+1][0]);
    }
    f[Npoints-1][2]=2*f[Npoints-2][2]-f[Npoints-3][2];
    f[Npoints-1][3]=2*f[Npoints-2][3]-f[Npoints-3][3];
}

void Spline::build_spline()
{
	double tw=1./12;
	double th=2./3;
	for(int i=2; i<n_points_-2; i++)
		f_[1][i]=tw*f_[i-2][0]-th*f_[0][i-1]+th*f_[0][i+1]-tw*f_[0][i+2];
	f_[1][1] = 0.5*(f_[0][2]-f_[0][0]);
	f_[1][npoints_-2] = 0.5*(f_[0][npoints_-1]-f_[0][npoints_-3]);
	f_[1][0] = 2*f_[1][1]-f_[1][2];
	f_[1][npoints_-1] = 0;
	for(int i=0; i<npoints_-1; i++)
	{
		f_[2][i] = 3*(f_[0][i+1]-f[0][i])-2*f_[1][i]-f_[1][i+1];
		f_[3][i] = f_[1][i+1]+f_[1][i]+2*(f_[0][i]-f_[0][i+1]);
	}
	f_[2][npoints_-1] = 2*f_[2][npoints_-2]-f_[2][npoints_-3];
	f_[3][npoints_-1] = 2*f_[3][npoints_-2]-f_[3][npoints_-3];
	created = true;
}*/

void Spline::build_spline()
{
	double tw = 1./12;
	double th = 2./3;
	for(int i=2; i<npoints_-2; i++)
		f_[1][i]=tw*f_[0][i-2]-th*f_[0][i-1]+th*f_[0][i+1]-tw*f_[0][i+2];
	f_[1][1] = 0.5*(f_[0][2]-f_[0][0]);
	f_[1][npoints_-2] = 0.5*(f_[0][npoints_-1]-f_[0][npoints_-3]);
	f_[1][0] = 2*f_[1][1]-f_[1][2];
	f_[1][npoints_-1] = 0;
	for(int i=0; i<npoints_-1; i++)
	{
		f_[2][i] = 3*(f_[0][i+1]-f_[0][i])-2*f_[1][i]-f_[1][i+1];
		f_[3][i] = f_[1][i+1]+f_[1][i]+2*(f_[0][i]-f_[0][i+1]);
	}
	f_[2][npoints_-1] = 2*f_[2][npoints_-2] - f_[2][npoints_-3];
	f_[3][npoints_-1] = 2*f_[3][npoints_-2] - f_[3][npoints_-3];
	created = true;
}

double Spline::eval_spline(double r) const
{
	if(!created)
		return 0;
	if((r>=max_)||(r<=min_)||(r!=r)) //last means r=nan
		return 0;
	const double dist = r-min_;
	const double nmb = dist/step_;
	double intpart = 0;
	const double fracpart = std::modf(nmb,(&intpart));
	const int i = (int) (intpart);
	//if((i > npoints_-1) || (i<0))
	//	std::cout << "Something is wrong in the eval_spline method" << std::endl;
	return f_[0][i] + fracpart*(f_[1][i]+fracpart*(f_[2][i] + fracpart*f_[3][i]));
}

double Spline::get_min() const
{
	if(!created)
		return -1000*step_;
	return min_;
}

double Spline::get_max() const
{
	if(!created)
		return 1000*step_;
	return max_;
}

double Spline::get_step() const
{
	return step_;
}

double Spline::get_ref_point() const
{
	return ref_point_;
}

bool Spline::is_created() const
{
	return created;
}
