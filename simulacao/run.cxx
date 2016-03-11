#include"Fluid.h"

int main(){

  Fluid Air(1.4,2000,0,2);

  double *a=new double[2000];
  for(int i=0;i<2000;++i)a[i]=1.2754;

  double *b=new double[2000];
  for(int i=0;i<2000;++i)b[i]=0;

  double *c=new double[2000];
  for(int i=0;i<2000;++i)c[i]=2.5e5+1.5e5*exp(-((double)i-1000)*((double)i-1000)/(2*35*35));

  Air.SetRho(a);Air.SetRhou(b);Air.SetEnergy(c);
  Air.Evolution2("evolution1.run",5e-3,5);

}
