#include"Fluid.h"

int main(){

// 20m diametro
// 12-13t massa
// 29.7km altitude
// K approx 1.8PJ
// v 19.16 km/sec
// 30km de simulacao

// definicoes da colisao
  ld xcol=0, ycol=0, zcol=29700, zbot=0, ztop=60000;
  ld xbot=-30000, xtop=30000, ybot=-30000, ytop=30000;
  ld disp=10./3, stepx=0.02;

// quero 20cm de grelha, donde 60km/(20cm)=3 k/c=3e5
  ld ap=1.4;int gridz=3000000;
  Fluid Air(ap,gridz,zbot,ztop);

// densidade do ar
  ld *a=new ld[gridz];
  for(int i=0;i<gridz;++i)a[i]=1.2754;

// velocidade do ar
  ld *b=new ld[gridz];
  for(int i=0;i<gridz;++i)b[i]=0;

// energia do ar
  ld *c=new ld[gridz];
  for(int i=0;i<gridz;++i)c[i]=2.5e5+2.29e12/sqrt(2*M_PI*disp*disp)*exp(-(stepx*i-zcol)*(stepx*i-zcol)/(2*disp*disp));

// evolucao
  Air.SetRho(a);Air.SetRhou(b);Air.SetEnergy(c);
  Air.Evolution3("evolution3.run",1e-5,2);

}
