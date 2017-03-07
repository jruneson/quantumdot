#ifndef __SIM__
#define __SIM__
#include<iostream>
#include<fstream>
#include<cmath>
#include<cfloat>
#include<complex>
#include<limits>

#define INCREASING 1
#define DECREASING 0
using namespace std;
class Symmeig{
	public:
		int n;
		double** z;
		double* d;   // d : eigenvalues,  z : columns are eigenvectors
		double* e;
		bool yesvecs;
		
		Symmeig(double** m, int size, bool yvec=true, int isort=0)
		{
			n = size;
			yesvecs=yvec;
			z = new double*[n];
			d = new double[n];
			e = new double[n];
			for(int i=0;i<n;i++)
			{
				d[i]=e[i]=0;
				z[i]= new double[n];
				for(int j=0;j<n;j++)
					z[i][j]=m[i][j];
			}
			
			tred2();
			tqli();
			if(isort==DECREASING)
				sort(DECREASING);
			else
				sort(INCREASING);
		};
		
		Symmeig(const char* fname, int size, bool yvec=true, int isort=0)
		{
			n = size;
			yesvecs=yvec;
			z = new double*[n];
			d = new double[n];
			e = new double[n];
			ifstream inp(fname);
			for(int i=0;i<n;i++)
			{
				d[i]=e[i]=0;
				z[i]= new double[n];
				for(int j=0;j<n;j++)
					inp>>z[i][j];
			}
			inp.close();
			
			tred2();
			tqli();
			if(isort==DECREASING)
				sort(DECREASING);
			else
				sort(INCREASING);
		};
		
		
		
		void sort(int);
		void tred2();
		void tqli();
		double pythag(const double a, const double b);
		
		~Symmeig()
		{
			delete d;
			delete e;
			for(int i=0;i<n;i++)
				delete z[i];
			delete z;
		}
		
		double SIGN(double a, double b)
		{
			int s = 1;
			if(b<0)
				s=-1;
			return s*abs(a);
		};

};


typedef class Symmeig Symmeig;

#endif
