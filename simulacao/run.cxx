#include"Fluid.h"

int main(){

  int n=1.5e6;double stepx=60000./(n-1);
  Fluid Air(1.4,n,0,60000);

  ld *a=new ld[n];for(int i=0;i<n;++i)a[i]=1.2754;
  ld *c=new ld[n];for(int i=0;i<n;++i)c[i]=2.5e5+2.29e12/sqrt(200*M_PI/9)*exp(-9./200*pow(i*stepx-29700,2));

  Air.SetRho(a);Air.SetEnergy(c);
  Air.Evolution3("evolution_final.run",5e-3,2);
}
