#ifndef __UNSIM__
#define __UNSIM__

class Unsymmeig {
	public:
    int n;
    double **a,**zz;
    double** wri;
    double* scale;
    int* perm;
    bool yesvecs,hessen;
    Unsymmeig(double ***aa,int _n,bool hessenb);
    ~Unsymmeig();
    void balance();
    void elmhes();
    void eltran();
    void hqr();
    void hqr2();
    void balbak();
    void SWAP(double&,double&);
    int MAX(int,int);
    double SIGN(double,double);
};

typedef class Unsymmeig Unsymmeig;

#endif
