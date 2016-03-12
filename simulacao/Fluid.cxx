#include"Fluid.h"

//___________________________________________________________
Fluid::Fluid(double adb, int nx, double bx, double ux){

 //  constructor, sets up a grid in x from bx to ux, with nx cells. It 
 // also sets up addiabatic to adb

  /*  Sets static variables  */
  numberx=nx;
  addiabatic=adb;
  upperx=ux;
  bottomx=bx;

  /*  Allocates memory for internal variables  */
  rho=new double[numberx];
  for(int i=0;i<numberx;++i)rho[i]=0;

  rhou=new double[numberx];
  for(int i=0;i<numberx;++i)rhou[i]=0;

  rhov=new double[numberx];
  for(int i=0;i<numberx;++i)rhov[i]=0;

  rhow=new double[numberx];
  for(int i=0;i<numberx;++i)rhow[i]=0;

  Energy=new double[numberx];
  for(int i=0;i<numberx;++i)Energy[i]=0;

}

//___________________________________________________________
Fluid::~Fluid(){

 // destructor

  /*  Deallocates memory  */
  delete []rho;
  delete []rhou;delete []rhov;delete []rhow;
  delete []Energy;

}

//___________________________________________________________
void Fluid::Evolution1(const char *output, double time, int nsteps){

 //  runs evolution of fluid through time with nsteps steps in time, 
 // writing the result into output, using simple numeric method

  ofstream aaa(output);
  double stept=time/(nsteps-1), stepx=(upperx-bottomx)/(numberx-1);

  for(int i=0;i<numberx;++i)if(i%20==0)// writes initial condition
    aaa<<"0\t"<<i*stepx<<"\t"<<rho[i]<<"\t"<<rhou[i]<<"\t"<<rhov[i]<<"\t"<<rhow[i]<<"\t"<<Energy[i]<<"\t"<<Pressure(i)<<endl;

  double **aux=new double*[5];
  for(int i=0;i<5;++i)aux[i]=new double[numberx];

  for(int i=1;i<nsteps;++i){// runs evolution in time

    for(int j=0;j<numberx;++j){// calculates evolution in the spacial grid

      double *f1=GetF(j-1), *f2=GetF(j+1);

      aux[0][j]=rho[j]   -stept/stepx*(f1[0]-f2[0]);
      aux[1][j]=rhou[j]  -stept/stepx*(f1[1]-f2[1]);
      aux[2][j]=rhov[j]  -stept/stepx*(f1[2]-f2[2]);
      aux[3][j]=rhow[j]  -stept/stepx*(f1[3]-f2[3]);
      aux[4][j]=Energy[j]-stept/stepx*(f1[4]-f2[4]);

    }

    for(int j=0;j<numberx;++j){
      rho[j]   =aux[0][j];
      rhou[j]  =aux[1][j];
      rhov[j]  =aux[2][j];
      rhow[j]  =aux[3][j];
      Energy[j]=aux[4][j];
    }

    for(int j=0;j<numberx;++j)if(j%20==0)// writes output into file
      aaa<<i<<"\t"<<j*stepx<<"\t"<<rho[j]<<"\t"<<rhou[j]<<"\t"<<rhov[j]<<"\t"<<rhow[j]<<"\t"<<Energy[j]<<"\t"<<Pressure(j)<<endl;

  }

}

//___________________________________________________________
void Fluid::Evolution2(const char *output, double time, int nsteps){

 //  runs evolution of fluid through time with nsteps steps in time, 
 // writing the result into output, using optimal numeric method

  ofstream aaa(output);
  double stept=time/(nsteps-1), stepx=(upperx-bottomx)/(numberx-1);

  for(int i=0;i<numberx;++i)if(i%10==0)// writes initial condition
    aaa<<"0\t"<<i*stepx<<"\t"<<rho[i]<<"\t"<<rhou[i]<<"\t"<<rhov[i]<<"\t"<<rhow[i]<<"\t"<<Energy[i]<<"\t"<<Pressure(i)<<endl;

  double **aux=new double*[5];
  for(int i=0;i<5;++i)aux[i]=new double[numberx];

  for(int i=1;i<nsteps;++i){// runs evolution in time

    for(int j=0;j<numberx;++j){// calculates evolution in the spacial grid

      double *f1=GetFb(j-1), *f2=GetFt(j);
      if(j==0){f1[0]=0;f1[1]=0;f1[2]=0;f1[3]=0;f1[4]=0;}
      if(j==numberx-1){f2[0]=0;f2[1]=0;f2[2]=0;f2[3]=0;f2[4]=0;}

      aux[0][j]=rho[j]   -stept/stepx*(f2[0]+f1[0]);
      aux[1][j]=rhou[j]  -stept/stepx*(f2[1]+f1[1]);
      aux[2][j]=rhov[j]  -stept/stepx*(f2[2]+f1[2]);
      aux[3][j]=rhow[j]  -stept/stepx*(f2[3]+f1[3]);
      aux[4][j]=Energy[j]-stept/stepx*(f2[4]+f1[4]);

    }

    for(int j=0;j<numberx;++j){
      rho[j]   =aux[0][j];
      rhou[j]  =aux[1][j];
      rhov[j]  =aux[2][j];
      rhow[j]  =aux[3][j];
      Energy[j]=aux[4][j];
    }

    for(int j=0;j<numberx;++j)if(j%10==0)// writes output into file
      aaa<<i*stept<<"\t"<<j*stepx<<"\t"<<rho[j]<<"\t"<<rhou[j]<<"\t"<<rhov[j]<<"\t"<<rhow[j]<<"\t"<<Energy[j]<<"\t"<<Pressure(j)<<endl;

  }

}

//___________________________________________________________
double* Fluid::GetF(int ix) const {

 // returns F(q) at cell ix

  if(ix<0)ix=0;if(ix>numberx-1)ix=numberx-1;

  double *f=new double[5];

  double u=rhou[ix]/rho[ix];// u speed
  double v=rhov[ix]/rho[ix];// v speed
  double w=rhow[ix]/rho[ix];// w speed

  /*  f vector  */
  f[0]=rho[ix]*u;
  f[1]=rho[ix]*u*u+Pressure(ix);
  f[2]=rho[ix]*u*v;
  f[3]=rho[ix]*u*w;
  f[4]=(Energy[ix]+Pressure(ix))*u;

  return f;

}

//___________________________________________________________
double* Fluid::GetG(int ix) const {

 // returns G(q) at cell ix

  double *g=new double[5];

  double u=rhou[ix]/rho[ix];// u speed
  double v=rhov[ix]/rho[ix];// v speed
  double w=rhow[ix]/rho[ix];// w speed

  /*  g vector  */
  g[0]=rho[ix]*v;
  g[1]=rho[ix]*u*v;
  g[2]=rho[ix]*v*v+Pressure(ix);
  g[3]=rho[ix]*v*w;
  g[4]=(Energy[ix]+Pressure(ix))*v;

  return g;

}

//___________________________________________________________
double* Fluid::GetH(int ix) const {

 // returns H(q) at cell ix

  double *h=new double[5];

  double u=rhou[ix]/rho[ix];// u speed
  double v=rhov[ix]/rho[ix];// v speed
  double w=rhow[ix]/rho[ix];// w speed

  /*  h vector  */
  h[0]=rho[ix]*w;
  h[1]=rho[ix]*u*w;
  h[2]=rho[ix]*v*w;
  h[3]=rho[ix]*w*w+Pressure(ix);
  h[4]=(Energy[ix]+Pressure(ix))*w;

  return h;

}

//___________________________________________________________
double Fluid::Pressure(int ix) const {

 // returns pressure at cell ix

  double u=rhou[ix]/rho[ix];// u speed
  double v=rhov[ix]/rho[ix];// v speed
  double w=rhow[ix]/rho[ix];// w speed

  double pressure=(addiabatic-1)*(Energy[ix]-rho[ix]/2*(u*u+v*v+w*w));
  return pressure;

}

//___________________________________________________________
double* Fluid::RoeAverage(int ix) const {

 // returns row average at cell ix

  double *q=new double[5];

  double pl=rho[ix], pr=rho[ix+1];// densities
  double El=Energy[ix], Er=Energy[ix+1];// energies
  double Pl=Pressure(ix), Pr=Pressure(ix+1);// pressures

  double ul=rhou[ix]/rho[ix], ur=rhou[ix+1]/rho[ix+1];// u speeds
  double vl=rhov[ix]/rho[ix], vr=rhov[ix+1]/rho[ix+1];// v speeds
  double wl=rhow[ix]/rho[ix], wr=rhow[ix+1]/rho[ix+1];// w speeds

  /*  roe average  */
  q[0]=(pl+pr)/2;
  q[1]=(sqrt(pl)*ul+sqrt(pr)*ur)/(sqrt(pl)+sqrt(pr));
  q[2]=(sqrt(pl)*vl+sqrt(pr)*vr)/(sqrt(pl)+sqrt(pr));
  q[3]=(sqrt(pl)*wl+sqrt(pr)*wr)/(sqrt(pl)+sqrt(pr));
  q[4]=((El+Pl)/sqrt(pl)+(Er+Pr)/sqrt(pr))/(sqrt(pl)+sqrt(pr));

  /*  checks physical significate  */
  double Emin=0.05/(addiabatic-1)+q[0]*(q[1]*q[1]+q[2]*q[2]+q[3]*q[3])/2;
  if(q[0]<0.05)q[0]=0.05;
  if(q[4]<Emin)q[4]=Emin;

  return q;

}

//___________________________________________________________
double* Fluid::Eigenvalues(double *roe) const {

 // returns wavespeeds of roe average

  double *lam=new double[5];

  double pressure=(addiabatic-1)/addiabatic*roe[0]*(roe[4]-roe[1]*roe[1]/2);// pressure
  double c=sqrt(addiabatic*pressure/roe[0]);// speed sound

  lam[0]=roe[1]-c;
  lam[1]=roe[1];
  lam[2]=roe[1];
  lam[3]=roe[1];
  lam[4]=roe[1]+c;

  return lam;

}

//___________________________________________________________
double** Fluid::Eigenvectors(double *roe) const {

 // returns eigenvectors of F(roe)

  double **r=new double*[5];
  for(int i=0;i<5;++i)r[i]=new double[5];

  double pressure=(addiabatic-1)/addiabatic*roe[0]*(roe[4]-roe[1]*roe[1]/2);// pressure
  double c=sqrt(addiabatic*pressure/roe[0]);// speed sound

  r[0][0]=1;
  r[1][0]=roe[1]-c;
  r[2][0]=roe[2];
  r[3][0]=roe[3];
  r[4][0]=roe[4]-roe[1]*c;

  r[0][1]=1;
  r[1][1]=roe[1];
  r[2][1]=roe[2];
  r[3][1]=roe[3];
  r[4][1]=(roe[1]*roe[1]+roe[2]*roe[2]+roe[3]*roe[3])/2;

  r[0][2]=0;
  r[1][2]=0;
  r[2][2]=1;
  r[3][2]=0;
  r[4][2]=roe[2];

  r[0][3]=0;
  r[1][3]=0;
  r[2][3]=0;
  r[3][3]=1;
  r[4][3]=roe[3];

  r[0][4]=1;
  r[1][4]=roe[1]+c;
  r[2][4]=roe[2];
  r[3][4]=roe[3];
  r[4][4]=roe[4]+roe[1]*c;

  return r;

}

//___________________________________________________________
double* Fluid::GetAlpha(int ix, double *roe) const {

 // return alpha coefficents of cell ix

  double *alpha=new double[5];
  double *deltaq=new double[5];

  /*  shift in q  */
  deltaq[0]=rho[ix+1]   -rho[ix];
  deltaq[1]=rhou[ix+1]  -rhou[ix];
  deltaq[2]=rhov[ix+1]  -rhov[ix];
  deltaq[3]=rhow[ix+1]  -rhow[ix];
  deltaq[4]=Energy[ix+1]-Energy[ix];

  double pressure=(addiabatic-1)/addiabatic*roe[0]*(roe[4]-roe[1]*roe[1]/2);// pressure
  double c=sqrt(addiabatic*pressure/roe[0]);// speed of sound

  double u=roe[1], v=roe[2], w=roe[3], H=roe[4];

  double ct=c*(2*H-u*u-v*v-w*w);
  double norm=u*u+v*v+w*w;

  /*  calculates inverse of eigenvector matrix  */
  double **m=new double*[5];
  for(int i=0;i<5;++i)m[i]=new double[5];

  m[0][0]=(2*H*u+(c-u)*norm)/(2*ct);
  m[0][1]=2*c*(H-norm)/ct;
  m[0][2]=c*v*(norm-2*H)/ct;
  m[0][3]=c*w*(norm-2*H)/ct;
  m[0][4]=(norm*(c+u)-2*H*u)/(2*ct);

  m[1][0]=(norm-2*c*u-2*H)/(2*ct);
  m[1][1]=2*c*u/ct;
  m[1][2]=0;
  m[1][3]=0;
  m[1][4]=(2*H-2*c*u-norm)/(2*ct);

  m[2][0]=-c*v/ct;
  m[2][1]=2*c*v/ct;
  m[2][2]=c*(2*H-norm)/ct;
  m[2][3]=0;
  m[2][4]=-c*v/ct;

  m[3][0]=-c*w/ct;
  m[3][1]=2*c*w/ct;
  m[3][2]=0;
  m[3][3]=c*(2*H-norm)/ct;
  m[3][4]=-c*w/ct;

  m[4][0]=c/ct;
  m[4][1]=-2*c/ct;
  m[4][2]=0;
  m[4][3]=0;
  m[4][4]=c/ct;

  /*  calculates alpha vector  */
  for(int i=0;i<5;++i)
    for(int j=0;j<5;++j)
      alpha[i]=m[i][j]*deltaq[j];

  return alpha;

}

//___________________________________________________________
double** Fluid::GetWaves(double *alpha, double **r) const {

 // returns waves at cell ix

  double **w=new double*[5];
  for(int i=0;i<5;++i)w[i]=new double[5];

  for(int i=0;i<5;++i)
    for(int j=0;j<5;++j)
      w[i][j]=alpha[j]*r[i][j];

  return w;

}

//___________________________________________________________
double* Fluid::GetFb(int ix) const {

 // returns F-(q) at cell ix

  double *roe=RoeAverage(ix);
  double *lambda=Eigenvalues(roe);
  double **r=Eigenvectors(roe); 
  double *alpha=GetAlpha(ix,roe);
  double **waves=GetWaves(alpha,r);

  double *result=new double[5];
  for(int i=0;i<5;++i)result[i]=0;

  for(int i=0;i<5;++i)if(lambda[i]<0)
    for(int j=0;j<5;++j)
      result[j]+=lambda[i]*waves[i][j];

  return result;

}

#include<iomanip>
//___________________________________________________________
double* Fluid::GetFt(int ix) const {

 // returns F+(q) at cell ix

  double *roe=RoeAverage(ix);
  double *lambda=Eigenvalues(roe);
  double **r=Eigenvectors(roe); 
  double *alpha=GetAlpha(ix,roe);
  double **waves=GetWaves(alpha,r);

  double *result=new double[5];
  for(int i=0;i<5;++i)result[i]=0;

  for(int i=0;i<5;++i)if(lambda[i]>0)
    for(int j=0;j<5;++j)
      result[j]+=lambda[i]*waves[i][j];

  return result;

}

//___________________________________________________________
double* Fluid::GetEnergy() const {

 // returns energy

  double *ret=new double[numberx];
  for(int i=0;i<numberx;++i)
    ret[i]=Energy[i];
  return ret;

}

//___________________________________________________________
double* Fluid::GetRho() const {

 // returns density

  double *ret=new double[numberx];
  for(int i=0;i<numberx;++i)
    ret[i]=rho[i];
  return ret;

}

//___________________________________________________________
double* Fluid::GetRhou() const {

 // returns momentum in x

  double *ret=new double[numberx];
  for(int i=0;i<numberx;++i)
    ret[i]=rhou[i];
  return ret;

}

//___________________________________________________________
double* Fluid::GetRhov() const {

 // returns momentum in y

  double *ret=new double[numberx];
  for(int i=0;i<numberx;++i)
    ret[i]=rhov[i];
  return ret;

}

//___________________________________________________________
double* Fluid::GetRhow() const {

 // returns momentum in z

  double *ret=new double[numberx];
  for(int i=0;i<numberx;++i)
    ret[i]=rhow[i];
  return ret;

}

//___________________________________________________________
void Fluid::SetEnergy(double *E){

 // sets energy to E

  for(int i=0;i<numberx;++i)
    Energy[i]=E[i];

}

//___________________________________________________________
void Fluid::SetRho(double *p){

 // sets rho to p

  for(int i=0;i<numberx;++i)
    rho[i]=p[i];

}

//___________________________________________________________
void Fluid::SetRhou(double *u){

 // sets momentum in x to u

  for(int i=0;i<numberx;++i)
    rhou[i]=u[i];

}

//___________________________________________________________
void Fluid::SetRhov(double *v){

 // sets momentum in y to v

  for(int i=0;i<numberx;++i)
    rhov[i]=v[i];

}

//___________________________________________________________
void Fluid::SetRhow(double *w){

 // sets momentum in z to w

  for(int i=0;i<numberx;++i)
    rhow[i]=w[i];

}









