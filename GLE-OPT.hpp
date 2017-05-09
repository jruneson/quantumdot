#include <vector>
#include "polymer.hpp"
#include "point.hpp"

#ifndef __GLEOPT__
#define __GLEOPT__

class RNumgen;
class GLE
{
	public:
		//send in full dt, it will divide by 2 inside
		GLE();
		GLE(int np, double** pvec, double _dt, double _T, double _mass, int _dim);
		GLE(int np,int nn, double*** pvec, double _dt, double _T, double _mass, int _dim);
		GLE(double* pvec, double _dt, double _T, double _mass, int _dim=1);
		GLE(std::vector<Polymer>&, const double&, const double&, const double&, 
			const int&, const int&, const int&, const bool);
		
		void run();
		
	private:
		void computeMatrices();
		double** psvec;
		int snum,ssize,psize;
		double dt, T,mass;
		int dim;
		double **Ap, **Cp, **Tmatrix, **Smatrix;
		RNumgen* gen;
		bool is_on;
		
		
		// linear algebra stuff
		void complex_matrix_inversion(double*** A, double*** Ainv, int n);
		double** pseudoInverse(double** M,int m,int n,double thres, double& cnumb);
		void real_matrix_squareroot(double** M, double** Minv, int n);
		void real_matrix_exponential(double** M, double** Mexp, int n, double fac);
		void real_matrix_product(double** A,double** B,double** P,int nr1, int nc1, int nr2, int nc2);
		void complex_matrix_product(double*** A,double** B,double*** P,int nr1, int nc1, int nr2, int nc2);
		void complex_matrix_product(double** A,double*** B,double*** P,int nr1, int nc1, int nr2, int nc2);
		void complex_matrix_product(double*** A,double*** B,double*** P,int nr1, int nc1, int nr2, int nc2);
		double complex_product(double* s1,double* s2,double* s3,int wh);
		double complex_division(double* s1,double* s2,int wh);
		double complex_product(double* s1,double* s2,int wh);
		
		
};

typedef class GLE GLE;


#endif
