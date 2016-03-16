#include"Fluid.h"

int main(){

// 20m diametro
// 12-13t massa
// 29.7km altitude
// K approx 1.8PJ
// v 19.16 km/sec
// 30km de simulacao

  ld ap=1.4, ad=0., aw=2500;
  Fluid Air(ap,2000,ad,aw);

  ld *a=new ld[2501];
  for(int i=0;i<2501;++i)a[i]=1.2754;

  ld *b=new ld[2501];
  for(int i=0;i<2501;++i)b[i]=0;

  ld *c=new ld[2501];
  for(int i=0;i<2501;++i)c[i]=2.5e5+2e7*exp(-((double)i-1250)*((double)i-1250)/(2*50*50));
/*
  int nx, ny, nz;// celulas
  ld ***E=new ld**[nx];
  for(int i=0;i<nx;++i)E[i]=new ld*[ny];
  for(int i=0;i<nx;++i)for(int j=0;j<ny;++j)E[i][j]=new ld[nz];
  double xcol, ycol, zcol;// colisao
  double dx, dy, dz;// passo 
  double disp;// raio
  for(int i=0;i<nx;++i)for(int j=0;j<ny;++j)for(int k=0;k<nz;++k)E[i][j][k]=Energy/sqrt(2*M_PI*disp*disp)*exp(-((dx*i-xcol)*(dx*i-col)+(dy*j-ycol)*(dy*j-ycol)+(dz*k-zcol)*(dz*k-zcol))/(2*disp*disp));
*/
  //Air.SetRho(a);Air.SetRhou(b);Air.SetEnergy(c);
  Air.Evolution1c("evolution1c.run",1e-5,2);

}
