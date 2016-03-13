#include"Fluid.h"

int main(){

  Fluid Air(1.4,250,0,0.2);

  long double *a=new long double[250];
  for(int i=0;i<250;++i)a[i]=1.2754;

  long double *b=new long double[250];
  for(int i=0;i<250;++i)b[i]=12.754*exp(-((long double)i-250)*((long double)i-250)/(2*10*10));

  long double *c=new long double[250];
  for(int i=0;i<250;++i)c[i]=2.5e5;

  Air.SetRho(a);Air.SetRhou(b);Air.SetEnergy(c);
  Air.Evolution1("evolution1.run",5e-3,2);

}
