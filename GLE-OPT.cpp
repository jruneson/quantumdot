#include<iostream>
#include<fstream>

#include "GLE-OPT.hpp"
#include "symmeig.hpp"
#include "unsymmeig.hpp"
#include "svd.hpp"
#include "numgen.hpp"

//When is snum initialized??
//Box-Muller sampling - bad choice - what can be improved?

using namespace std;

GLE::GLE(std::vector<Polymer>& pols, const double& timestep, const double& temperature,
		const double& part_mass, const int& num_beads, const int& num_parts, const int& dims,
		const bool thermostat_on)
{
	psize = num_beads*num_parts;
	T = temperature;
	dt = timestep;
	mass = part_mass; //note, must be the same for all particles
	dim = dims;
	is_on = thermostat_on;
	computeMatrices();
	ssize = psize*(snum+1)*dim; //size of total vector of p*dim momenta, and snum extra variables for each.
	psvec = new double*[ssize];
	int pknt1=0;
	int pknt2=0;
	int pknt3=0;	
	for(int i=0;i<ssize;i+=(snum+1))
	{
		Polymer& pol = pols[pknt2];
		// [bead][part][dim]
		psvec[i] = &pol.vels[pknt1][pknt3];
		pknt3++;
		if(pknt3==dim)
		{
			pknt3=0;
			pknt2++;
		}
		if(pknt2==num_parts)
		{
			pknt2=0;
			pknt1++;
		}
		for(int j=0;j<snum;j++)
		{
			psvec[i+j+1]= new double;
			*psvec[i+j+1]=0;
		}
	}
	gen = new RNumgen();
}

/*
GLE::GLE(int np, double** pvec, double _dt, double _T, double _mass, int _dim)
{
	psize=np; //num_beads
	T=_T;
	dt=_dt;
	mass=_mass;
	dim=_dim;
	computeMatrices();
	ssize = psize*(snum+1)*dim;
	psvec = new double*[ssize];
	int pknt1=0;
	int pknt2=0;
	for(int i=0;i<ssize;i+=(snum+1))
	{
		psvec[i]=&pvec[pknt1][pknt2];
		pknt2++;
		if(pknt2==dim)
		{
			pknt2=0;
			pknt1++;
		}
		for(int j=0;j<snum;j++)
		{
			psvec[i+j+1]= new double;
			*psvec[i+j+1]=0;
		}
	}
	gen = new RNumgen();
	
}*/



GLE::GLE(int np, int nn, double*** pvec, double _dt, double _T, double _mass, int _dim)
{
	//np = num_beads
	//nn = num_particles
	//T = temperature
	psize=np*nn; //total number of beads
	T=_T;
	dt=_dt;
	mass=_mass;
	dim=_dim;
	is_on = true;
	computeMatrices();
	ssize = psize*(snum+1)*dim; //size of total vector of p*dim momenta, and snum extra variables for each.
	psvec = new double*[ssize];
	int pknt1=0;
	int pknt2=0;
	int pknt3=0;
	
	for(int i=0;i<ssize;i+=(snum+1))
	{
		// [bead][part][dim]
		psvec[i]=&pvec[pknt1][pknt2][pknt3];
		pknt3++;
		if(pknt3==dim)
		{
			pknt3=0;
			pknt2++;
		}
		if(pknt2==nn)
		{
			pknt2=0;
			pknt1++;
		}
		for(int j=0;j<snum;j++)
		{
			psvec[i+j+1]= new double;
			*psvec[i+j+1]=0;
		}
	}
	gen = new RNumgen();
}


/*
GLE::GLE(double* pvec, double _dt, double _T, double _mass, int _dim)
{
	psize=1;
	T=_T;
	dt=_dt;
	mass=_mass;
	dim=_dim;
	computeMatrices();
	ssize = psize*(snum+1)*dim;
	psvec = new double*[ssize];
	int pknt2=0;
	for(int i=0;i<ssize;i+=(snum+1))
	{
		psvec[i]=&pvec[pknt2];
		pknt2++;
		for(int j=0;j<snum;j++)
		{
			psvec[i+j+1]= new double;
			*psvec[i+j+1]=0;
		}
	}
	gen = new RNumgen();
	
}*/

void GLE::run()
{
	if(is_on)
	{
		double xi[ssize];
		double v[ssize];
		for(int i=0;i<ssize;i++)
			xi[i]=gen->gaussGen(0,2);   // box-muller gaussian sampling, argument is 2sigma (bad choice)
		
		// psnew = T*psvec + S*xi
		
		for(int i=0;i<ssize;i++)
		{
			v[i]=0;
			int imod=i%(snum+1);
			int nmod=i/(snum+1);
			for(int j=0;j<snum+1;j++)
			{
				v[i] += Tmatrix[imod][j]*(*psvec[nmod*(snum+1)+j]);  // T*psvec part
				v[i] += Smatrix[imod][j]*xi[nmod*(snum+1)+j];	// S*xi part
			}		
		}
		for(int i=0;i<ssize;i++)
			(*psvec[i])=v[i];	
	}
}

void GLE::computeMatrices()
{
	int s;
	double temp=T;
	double boltzmann=1;
	
	ifstream inp("Amatrix.dat");
	if(!(inp.good()))
		cout << "Warning, Amatrix.dat is missing!" << endl;
	inp>>s;
	s--;
	snum=s;
	Ap = new double*[s+1];
	double** Dp = new double*[s+1];
	Tmatrix = new double*[s+1];
	
	for(int i=0;i<s+1;i++)
	{
		Ap[i]= new double[s+1];
		Dp[i]= new double[s+1];
		Tmatrix[i]= new double[s+1];
		for(int j=0;j<s+1;j++)
			Ap[i][j]=0;	
	}	
	
	for(int i=0;i<s+1;i++)
		for(int j=0;j<s+1;j++)
			inp>>Ap[i][j];
	
	
	for(int i=0;i<s+1;i++)
		for(int j=0;j<s+1;j++)
			Dp[i][j]=boltzmann*temp*(Ap[i][j]+Ap[j][i])/mass;

	
	
	// T matrix : T = exp(-0.5*dt*Ap)
	
	//cout<<"***** T matrix:"<<endl;
	real_matrix_exponential(Ap,Tmatrix,s+1,-0.5*dt);
	
	
	double **T2;  // T2 = exp(-0.5*dt*Ap^T)  it is a matrix required for S.
	
	T2 = new double*[s+1];
	for(int i=0;i<s+1;i++)
	{
		T2[i]= new double[s+1];
		for(int j=0;j<s+1;j++)
			T2[i][j]=Tmatrix[j][i];
	}
	
	// kronecker products I x Ap  and Ap x I and vector representation of Dp
	double Dpvec[(s+1)*(s+1)];
	double **kron,**kron_inverse;
	kron= new double*[(s+1)*(s+1)];
	kron_inverse= new double*[(s+1)*(s+1)];
	double Id[s+1][s+1];
	
	for(int i=0;i<(s+1)*(s+1);i++)
	{
		kron[i]= new double[(s+1)*(s+1)];
		kron_inverse[i]= new double[(s+1)*(s+1)];
	}
	int kn=0;
	for(int i=0;i<s+1;i++)
	{
		for(int j=0;j<s+1;j++)
		{
			Id[i][j]=0;
			Dpvec[kn++]=Dp[j][i];
		}
		Id[i][i]=1;
	}	
	
	
	

	for(int i=0;i<s+1;i++)
		for(int j=0;j<s+1;j++)
		{
			double a1 = Ap[i][j];
			double a2 = Id[i][j];
			for(int ii=0;ii<s+1;ii++)
				for(int jj=0;jj<s+1;jj++)
					kron[i*(s+1)+ii][j*(s+1)+jj]=a1*Id[ii][jj]+a2*Ap[ii][jj];
				
			
		}
		
	//cout<<"Inverting kronecher matrix..."<<endl;	
	
	//real_matrix_inversion(kron,kron_inverse,(s+1)*(s+1));
	double cnu;
	kron_inverse=pseudoInverse(kron,(s+1)*(s+1),(s+1)*(s+1),1E-12,cnu);
	
	double cpvec[(s+1)*(s+1)];
	for(int i=0;i<(s+1)*(s+1);i++)
	{
		cpvec[i]=0;
		for(int j=0;j<(s+1)*(s+1);j++)
			cpvec[i]+=kron_inverse[i][j]*Dpvec[j];
	}	
	
	
	double** Cp = new double*[s+1];
	kn=0;
	for(int i=0;i<s+1;i++)
		Cp[i]= new double[s+1];
	
	for(int i=0;i<s+1;i++)
	{
		for(int j=0;j<s+1;j++)
		{
			Cp[j][i]=cpvec[kn++];
			if(abs(Cp[j][i])<1E-12)
				Cp[j][i]=0;
		}		
	}
	
	// products for SS^T
	double p1[s+1][s+1], p2[s+1][s+1];
	double** F= new double*[s+1];
	for(int i=0;i<s+1;i++)
		for(int j=0;j<s+1;j++)
		{
			p1[i][j]=0;
			for(int k=0;k<s+1;k++)
				p1[i][j]+=Tmatrix[i][k]*Cp[k][j];
		}
	
	for(int i=0;i<s+1;i++)
	{
		F[i]= new double[s+1];
		for(int j=0;j<s+1;j++)
		{
			p2[i][j]=0;
			for(int k=0;k<s+1;k++)
				p2[i][j]+=p1[i][k]*T2[k][j];
			
			F[i][j]=Cp[i][j]-p2[i][j];
			
		}
	}
	
	
	
	
	
	Smatrix = new double*[s+1];
	for(int i=0;i<s+1;i++)
		Smatrix[i]= new double[s+1];
	
	for(int i=0;i<s+1;i++)
	{
		for(int j=0;j<s+1;j++)
		{
			if(abs(F[i][j])<1E-16)
				F[i][j]=0;   
		}
	}
		
	real_matrix_squareroot(F,Smatrix,s+1);

}






double GLE::complex_product(double* s1,double* s2,int wh)
{
	if(wh==0)  // Re
		return s1[0]*s2[0]-s1[1]*s2[1];
	else
		return s1[0]*s2[1]+s1[1]*s2[0];
}

double GLE::complex_division(double* s1,double* s2,int wh)
{
	double mod = s2[0]*s2[0]+s2[1]*s2[1];
	if(wh==0)
	{
		return (s1[0]*s2[0]+s1[1]*s2[1])/mod;
	}
	else
	{
		return (s1[1]*s2[0]-s1[0]*s2[1])/mod;
	}
}

double GLE::complex_product(double* s1,double* s2,double* s3,int wh)
{
	double re12= s1[0]*s2[0]-s1[1]*s2[1];
	double im12= s1[0]*s2[1]+s1[1]*s2[0];
	
	if(wh==0)  // Re
		return re12*s3[0]-im12*s3[1];
	else
		return re12*s3[1]+im12*s3[0];
}

void GLE::complex_matrix_product(double*** A,double*** B,double*** P,int nr1, int nc1, int nr2, int nc2)
{
	for(int i=0;i<nr1;i++)
	{
		for(int j=0;j<nc2;j++)
		{
			P[i][j][0]=0;
			P[i][j][1]=0;
			
			for(int k=0;k<nr2;k++){
				P[i][j][0]+= complex_product(A[i][k],B[k][j],0);
				P[i][j][1]+= complex_product(A[i][k],B[k][j],1);
			}
		}
	}
	
}

void GLE::complex_matrix_product(double** A,double*** B,double*** P,int nr1, int nc1, int nr2, int nc2)
{
	for(int i=0;i<nr1;i++)
	{
		for(int j=0;j<nc2;j++)
		{
			P[i][j][0]=0;
			P[i][j][1]=0;
			
			for(int k=0;k<nr2;k++){
				P[i][j][0]+= A[i][k]*B[k][j][0];//complex_product(A[i][k],B[k][j],0);
				P[i][j][1]+= A[i][k]*B[k][j][1];//complex_product(A[i][k],B[k][j],1);
			}
		}
	}
	
}

void GLE::complex_matrix_product(double*** A,double** B,double*** P,int nr1, int nc1, int nr2, int nc2)
{
	for(int i=0;i<nr1;i++)
	{
		for(int j=0;j<nc2;j++)
		{
			P[i][j][0]=0;
			P[i][j][1]=0;
			
			for(int k=0;k<nr2;k++){
				P[i][j][0]+= A[i][k][0]*B[k][j];//complex_product(A[i][k],B[k][j],0);
				P[i][j][1]+= A[i][k][1]*B[k][j];//complex_product(A[i][k],B[k][j],1);
			}
		}
	}
	
}

void GLE::real_matrix_product(double** A,double** B,double** P,int nr1, int nc1, int nr2, int nc2)
{
	for(int i=0;i<nr1;i++)
	{
		for(int j=0;j<nc2;j++)
		{
			P[i][j]=0;
			
			for(int k=0;k<nr2;k++)
				P[i][j]+= A[i][k]*B[k][j];
			
		}
	}
	
}

void GLE::real_matrix_exponential(double** M, double** Mexp, int n,double fac)
{
	Unsymmeig m(&M,n,false);
	
	// eigenvalues M.wri[i]
	// eigenvectors M.zz[][]

       double ***U, ***Uinv;
       U= new double**[n];
       Uinv= new double**[n];
       
        int cn=0;
        
	for(int i=0;i<n;i++)
        {
		U[i]= new double*[n];
		Uinv[i]= new double*[n];
		for(int j=0;j<n;j++)
		{
			U[i][j]=new double[2];
			Uinv[i][j]=new double[2];
		}
	}
	
	for(int i=0;i<n;i++)
        {
		if(m.wri[i][1]==0)
                {
                        for(int j=0;j<n;j++)
                        {
                                U[j][i][1]=0;
                                U[j][i][0]=m.zz[j][cn];
				
                        }
                        cn++;
                }
                else
                {
                        for(int j=0;j<n;j++)
                        {
                                U[j][i][0]=m.zz[j][cn];
                                U[j][i][1]=m.zz[j][cn+1];
                                U[j][i+1][0]=m.zz[j][cn];
                                U[j][i+1][1]=-m.zz[j][cn+1];
				
                        }
                        i++;
			cn++;
			cn++;
                }
        }
	


	
	// computing Uinv
	
	complex_matrix_inversion(U,Uinv,n);
		
	
	
	
	
	double expo_eig[n][2], one[2]={1,0};
	for(int i=0;i<n;i++)
	{
		expo_eig[i][0]=exp(fac*m.wri[i][0])*cos(fac*m.wri[i][1]);
		expo_eig[i][1]=exp(fac*m.wri[i][0])*sin(fac*m.wri[i][1]);
		
	}
	
		
	// U*D
	double ud[n][n][2];
	
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<n;j++)
		{
			ud[i][j][0]=complex_product(U[i][j],expo_eig[j],0);
			ud[i][j][1]=complex_product(U[i][j],expo_eig[j],1);
		}
	}
	
	// ud*Uinv
	double fs[n][n][2];
//	cout<<endl<<"**** exp(-0.5dtM)="<<endl;
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<n;j++)
		{
			fs[i][j][0]=fs[i][j][1]=0;
			for(int k=0;k<n;k++)
			{
				fs[i][j][0] += complex_product(ud[i][k],Uinv[k][j],0);
				fs[i][j][1] += complex_product(ud[i][k],Uinv[k][j],1);
			}
			Mexp[i][j]=fs[i][j][0];
			//cout<<fs[i][j][0]<<"+"<<fs[i][j][1]<<"i   ";
		}
	//	cout<<endl;
	
	}
	//cout<<endl;
}
void GLE::real_matrix_squareroot(double** M, double** Minv, int n)
{
	Unsymmeig m(&M,n,false);
	
	// eigenvalues M.wri[i]
	// eigenvectors M.zz[][]
	
	
        double ***U, ***Uinv;
       U= new double**[n];
       Uinv= new double**[n];
	for(int i=0;i<n;i++)
        {
		U[i]= new double*[n];
		Uinv[i]= new double*[n];
		for(int j=0;j<n;j++)
		{
			U[i][j]=new double[2];
			Uinv[i][j]=new double[2];
		}
	}       

        int cn=0;
        for(int i=0;i<n;i++)
        {
		if(m.wri[i][1]==0)
                {
                        for(int j=0;j<n;j++)
                        {
                                U[j][i][1]=0;
                                U[j][i][0]=m.zz[j][cn];
			}
                        cn++;
                }
                else
                {
                        for(int j=0;j<n;j++)
                        {
                                U[j][i][0]=m.zz[j][cn];
                                U[j][i][1]=m.zz[j][cn+1];
                                U[j][i+1][0]=m.zz[j][cn];
                                U[j][i+1][1]=-m.zz[j][cn+1];
			}
                        i++;
			cn++;
			cn++;
                }
        }
	
	
	// inversion of U
	complex_matrix_inversion(U,Uinv,n);
	
	double sqr_eig[n][2];
	for(int i=0;i<n;i++)
	{
		double sq = m.wri[i][0];
		sqr_eig[i][0]=sqr_eig[i][1]=0;
		if(sq>0)
			sqr_eig[i][0]=sqrt(sq);
		else
			sqr_eig[i][1]=sqrt(-sq);
	}
	
		
	// U*D
	double ud[n][n][2];
	
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<n;j++)
		{
			ud[i][j][0]=complex_product(U[i][j],sqr_eig[j],0);
			ud[i][j][1]=complex_product(U[i][j],sqr_eig[j],1);
		}
	}
	
	// ud*Uinv
	double fs[n][n][2];
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<n;j++)
		{
			fs[i][j][0]=fs[i][j][1]=0;
			for(int k=0;k<n;k++)
			{
				fs[i][j][0] += complex_product(ud[i][k],Uinv[k][j],0);
				fs[i][j][1] += complex_product(ud[i][k],Uinv[k][j],1);
			}
			Minv[i][j]=fs[i][j][0];
		}
	
	}
	
}


double** GLE::pseudoInverse(double** M,int m,int n,double thres, double& cnumb)
{
        double **out = new double*[n];
        for(int i=0;i<n;i++)
                out[i]= new double[m];

        SVD mdec(M,m,n);
        int msing = m;
        if(m>n)
                msing=n;
        int scut=msing;
        for(int i=0;i<msing;i++)
       	 	if(mdec.w[i]<thres)
                {
                        scut=i;
                        break;
                }
       
	double singmax=mdec.w[0];
        double singmin=mdec.w[scut-1];
        cnumb=singmax/singmin;
        double **t1 = new double*[n];
        for(int i=0;i<n;i++)
                t1[i]= new double[m];
        for(int i=0;i<n;i++)
                for(int j=0;j<m;j++)
                {
                        if(i>=scut)
                                t1[i][j]=0;
                        else
                                t1[i][j]=(1./mdec.w[i])*mdec.u[j][i];
                }

        for(int i=0;i<n;i++)
                for(int j=0;j<m;j++)
                {
                        out[i][j]=0;
                        for(int k=0;k<n;k++)
                                out[i][j]+=mdec.v[i][k]*t1[k][j];
                }
        for(int i=0;i<n;i++)
                delete t1[i];
        delete t1;
        return out;

}


void GLE::complex_matrix_inversion(double*** A, double*** Ainv, int n)
{
	// with a block form real matrix and SVD
	
	double **Ab = new double*[2*n];
	for(int i=0;i<2*n;i++)
		Ab[i]= new double[2*n];
	
	for(int i=0;i<n;i++)
		for(int j=0;j<n;j++)
		{
			Ab[i][j]=A[i][j][0];
			Ab[i][j+n]=-A[i][j][1];
			Ab[i+n][j]=A[i][j][1];
			Ab[i+n][j+n]=A[i][j][0];
		}
	
	
	
	SVD m(Ab,2*n,2*n);
	
	double W[n],U[n][n][2],V[n][n][2];
	int cn=0;
	for(int i=0;i<n;i++)
	{
		W[i]=m.w[cn];
		for(int j=0;j<n;j++)
		{
			U[j][i][0]=m.u[j][cn];
			U[j][i][1]=m.u[j+n][cn];
			V[j][i][0]=m.v[j][cn];
			V[j][i][1]=m.v[j+n][cn];
		}
		cn+=2;
		
	}
	
	
	double p1[n][n][2];
	for(int i=0;i<n;i++)
		for(int j=0;j<n;j++)
		{
			p1[i][j][0]=V[i][j][0]*(1./W[j]);
			p1[i][j][1]=V[i][j][1]*(1./W[j]);
		}
		
	double Uadj[n][n][2];
	for(int i=0;i<n;i++)
		for(int j=0;j<n;j++)
		{
			Uadj[i][j][0]=U[j][i][0];
			Uadj[i][j][1]=-U[j][i][1];
		}
		
	for(int i=0;i<n;i++)
		for(int j=0;j<n;j++)
		{
			Ainv[i][j][0]=Ainv[i][j][1]=0;
			for(int k=0;k<n;k++)
			{
				Ainv[i][j][0]+=complex_product(p1[i][k],Uadj[k][j],0);
				Ainv[i][j][1]+=complex_product(p1[i][k],Uadj[k][j],1);
			}
		}	
}
