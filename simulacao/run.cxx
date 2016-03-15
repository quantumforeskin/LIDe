#include"Fluid.h"

int main(){

  Fluid Air(1.4,2501,0,2500);

  long double *a=new long double[2501];
  for(int i=0;i<2501;++i)a[i]=1.2754;

  long double *b=new long double[2501];
  for(int i=0;i<2501;++i)b[i]=0;

  long double *c=new long double[2501];
  for(int i=0;i<2501;++i)c[i]=2.5e5+2e7*exp(-((double)i-1250)*((double)i-1250)/(2*50*50));

  Air.SetRho(a);Air.SetRhou(b);Air.SetEnergy(c);
  Air.Evolution3("evolution3.run",1e-5,2);

}
