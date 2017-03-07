#ifndef __SVD__
#define __SVD__

class SVD{
    public:
    	SVD(double** a, int _m, int _n);
	~SVD();
	int m,n;
    	double **u,**v;
    	double *w;
    	double eps, tsh;
    private:
    void decompose();
    void reorder();
    double pythag(const double a, const double b);
    double SQR(double);
    double SIGN(double,double);
    double MAX(double,double);
    double MIN(double,double);
    int MAX(int,int);
    int MIN(int,int);
};

typedef class SVD SVD;

#endif
