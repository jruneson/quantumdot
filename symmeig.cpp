#include "symmeig.hpp"

void Symmeig::sort(int icr)
{
	double sgn = 1;
	if(icr==INCREASING)
	 sgn=-1;
	 
	int i;
	for(int j=1;j<n;j++)
	{
		double x = d[j];
		double eig[n];
		for(int ii=0;ii<n;ii++)
			eig[ii]=z[ii][j];
			
		for(i=j-1;i>=0;i--)
		{
			if(sgn*d[i]>=sgn*x) break;
			d[i+1]=d[i];
			for(int k=0;k<n;k++)
				z[k][i+1]=z[k][i];
		}
		d[i+1]=x;
		for(int k=0;k<n;k++)
			z[k][i+1]=eig[k];
	}
}
void Symmeig::tred2()
{
	int l,k,j,i;
	double scale,hh,h,g,f;
	for(i=n-1;i>0;i--){
		l=i-1;
		h=scale=0;
		if(l>0){
			for(k=0;k<i;k++)
				scale+=abs(z[i][k]);
			if(scale==0.0)
				e[i]=z[i][l];
			else {
				for(k=0;k<i;k++){
					z[i][k]/=scale;
					h+= z[i][k]*z[i][k];
				}
				f=z[i][l];
				g=(f>=0.0 ? -sqrt(h) : sqrt(h));
				e[i]=scale*g;
				h -= f*g;
				z[i][l]=f-g;
				f=0;
				for(j=0;j<i;j++){
					if(yesvecs)
						z[j][i]=z[i][j]/h;
					g=0;
					for(k=0;k<j+1;k++)
						g+= z[j][k]*z[i][k];
					for(k=j+1;k<i;k++)
						g+= z[k][j]*z[i][k];
					e[j]=g/h;
					f+= e[j]*z[i][j];
				}
				hh=f/(h+h);
				for(j=0;j<i;j++){
					f=z[i][j];
					e[j]=g=e[j]-hh*f;
					for(k=0;k<j+1;k++)
						z[j][k]-=(f*e[k]+g*z[i][k]);
				}
			}
		} else
			e[i]=z[i][l];
		d[i]=h;
	}
	if(yesvecs) d[0]=0;
	e[0]=0;
	for(i=0;i<n;i++){
		if(yesvecs){
			if(d[i]!=0.0){
				for(j=0;j<i;j++){
					g=0.0;
					for(k=0;k<i;k++)
						g+= z[i][k]*z[k][j];
					for(k=0;k<i;k++)
						z[k][j]-=g*z[k][i];
				}
			}
			d[i]=z[i][i];
			z[i][i]=1;
			for(j=0;j<i;j++) z[j][i]=z[i][j]=0;
		} else {
			d[i]=z[i][i];
		}
	}
}
void Symmeig::tqli()
{
	int m,l,iter,i,k;
	double s,r,p,g,f,dd,c,b;
	const double EPS = numeric_limits<double>::epsilon();
	for(i=1;i<n;i++) 
		e[i-1]=e[i];
	e[n-1]=0;
	for(l=0;l<n;l++){
		iter=0;
		do{
			for(m=l;m<n-1;m++){
				dd=abs(d[m])+abs(d[m+1]);
				if(abs(e[m])<=EPS*dd) break;
			}
			if(m!=l){
				if(iter++==30) cout<<"Too many iterations in tqli"<<endl;
				g=(d[l+1]-d[l])/(2*e[l]);
				r=pythag(g,1.0);
				g=d[m]-d[l]+e[l]/(g+SIGN(r,g));
				s=c=1;
				p=0;
				for(i=m-1;i>=l;i--){
					f=s*e[i];
					b=c*e[i];
					e[i+1]=(r=pythag(f,g));
					if(r==0.0){
						d[i+1]-=p;
						e[m]=0;
						break;
					}
					s=f/r;
					c=g/r;
					g=d[i+1]-p;
					r=(d[i]-g)*s+2*c*b;
					d[i+1]=g+(p=s*r);
					g=c*r-b;
					if (yesvecs){
						for(k=0;k<n;k++){
							f=z[k][i+1];
							z[k][i+1]=s*z[k][i]+c*f;
							z[k][i]=c*z[k][i]-s*f;
						}
					}
				}
				if(r==0.0&&i>=1) continue;
				d[l] -=p;
				e[l] = g;
				e[m]=0;
			}
		}while(m!=l);
	}
}

double Symmeig::pythag(const double a, const double b)
{
	double absa = abs(a);
	double absb = abs(b);
	return (absa > absb ? absa*sqrt(1+pow(absb/absa,2)) :
		( absb == 0.0 ? 0.0 : absb*sqrt(1+pow(absa/absb,2))));
}
