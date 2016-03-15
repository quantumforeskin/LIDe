#include"Fluid.h"

//___________________________________________________________
Fluid::Fluid(long double adb, int nx, long double bx, long double ux){

 //  constructor, sets up a grid in x from bx to ux, with nx cells. It 
 // also sets up addiabatic to adb

  /*  Sets static variables  */
  numberx=nx;
  addiabatic=adb;
  upperx=ux;
  bottomx=bx;

  /*  Allocates memory for internal variables  */
  rho=new long double[numberx];
  for(int i=0;i<numberx;++i)rho[i]=0;

  rhou=new long double[numberx];
  for(int i=0;i<numberx;++i)rhou[i]=0;

  rhov=new long double[numberx];
  for(int i=0;i<numberx;++i)rhov[i]=0;

  rhow=new long double[numberx];
  for(int i=0;i<numberx;++i)rhow[i]=0;

  Energy=new long double[numberx];
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
void Fluid::Evolution1(const char *output, long double time, int nsteps){

 //  runs evolution of fluid through time with nsteps steps in time, 
 // writing the result into output, using simple numeric method

  ofstream aaa(output);
  long double stept=time/(nsteps-1), stepx=(upperx-bottomx)/(numberx-1);

  for(int i=0;i<numberx;++i)// writes initial condition
    aaa<<"0\t"<<i*stepx<<"\t"<<rho[i]<<"\t"<<rhou[i]<<"\t"<<rhov[i]<<"\t"<<rhow[i]<<"\t"<<Energy[i]<<"\t"<<Pressure(i)<<endl;

  long double **aux=new long double*[5];
  for(int i=0;i<5;++i)aux[i]=new long double[numberx];

  for(int i=1;i<nsteps;++i){// runs evolution in time

    for(int j=0;j<numberx;++j){// calculates evolution in the spacial grid

      long double *f1=GetF(j-1), *f2=GetF(j+1);

      cout<<f1[1]<<" "<<f2[1]<<"\t"<<f1[4]<<" "<<f2[4]<<endl;

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

    for(int j=0;j<numberx;++j)// writes output into file
      aaa<<i<<"\t"<<j*stepx<<"\t"<<rho[j]<<"\t"<<rhou[j]<<"\t"<<rhov[j]<<"\t"<<rhow[j]<<"\t"<<Energy[j]<<"\t"<<Pressure(j)<<endl;

  }

  cout<<endl<<stept<<" "<<stepx<<endl<<endl;

}

//___________________________________________________________
void Fluid::Evolution2(const char *output, long double time, int nsteps){

 //  runs evolution of fluid through time with nsteps steps in time, 
 // writing the result into output, using optimal numeric method

  ofstream aaa(output);
  long double stept=time/(nsteps-1), stepx=(upperx-bottomx)/(numberx-1);

  for(int i=0;i<numberx;++i)if(i%10==0)// writes initial condition
    aaa<<"0\t"<<i*stepx<<"\t"<<rho[i]<<"\t"<<rhou[i]<<"\t"<<rhov[i]<<"\t"<<rhow[i]<<"\t"<<Energy[i]<<"\t"<<Pressure(i)<<endl;

  long double **aux=new long double*[5];
  for(int i=0;i<5;++i)aux[i]=new long double[numberx];

  for(int i=1;i<nsteps;++i){// runs evolution in time

    for(int j=0;j<numberx;++j){// calculates evolution in the spacial grid

      long double *f1=GetFb(j-1), *f2=GetFt(j);
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
void Fluid::Evolution3(const char *output, long double time, int nsteps){

 //  runs evolution of fluid through time with nsteps steps in time, 
 // writing the result into output, using optimal numeric method

  ofstream aaa(output);
  aaa<<setprecision(10);
  long double stepx=(upperx-bottomx)/(numberx-1), stept=time/(nsteps-1);

  for(int i=0;i<numberx;++i)// writes initial condition
    aaa<<"0\t"<<i*stepx<<"\t"<<rho[i]<<"\t"<<rhou[i]<<"\t"<<rhov[i]<<"\t"<<rhow[i]<<"\t"<<Energy[i]<<"\t"<<Pressure(i)<<endl;

  double **aux=new double*[5];// auxiliary vector
  for(int i=0;i<5;++i)aux[i]=new double[numberx+1];

  for(int i=1;i<nsteps;++i){// runs evolution

    for(int j=0;j<numberx;++j){// calculates f(q) on grid
      long double *g=GetF(j);
      aux[0][j]=g[0];
      aux[1][j]=g[1];
      aux[2][j]=g[2];
      aux[3][j]=g[3];
      aux[4][j]=g[4];

    }

    for(int j=1;j<numberx-1;++j){
      aaa<<1<<"\t"<<j*stepx<<"\t"<<flush;
      aaa<<(aux[0][j+1]-aux[0][j-1])/(2*stepx)<<"\t"<<flush;
      aaa<<(aux[1][j+1]-aux[1][j-1])/(2*stepx)<<"\t"<<flush;
      aaa<<(aux[2][j+1]-aux[2][j-1])/(2*stepx)<<"\t"<<flush;
      aaa<<(aux[3][j+1]-aux[3][j-1])/(2*stepx)<<"\t"<<flush;
      aaa<<(aux[4][j+1]-aux[4][j-1])/(2*stepx)<<endl;
    }
    for(int j=1;j<numberx-1;++j){
      aaa<<2<<"\t"<<j*stepx<<"\t"<<flush;
      aaa<<(aux[0][j+1]-aux[0][j])/(stepx)<<"\t"<<flush;
      aaa<<(aux[1][j+1]-aux[1][j])/(stepx)<<"\t"<<flush;
      aaa<<(aux[2][j+1]-aux[2][j])/(stepx)<<"\t"<<flush;
      aaa<<(aux[3][j+1]-aux[3][j])/(stepx)<<"\t"<<flush;
      aaa<<(aux[4][j+1]-aux[4][j])/(stepx)<<endl;
    }
/*
    rho[0]   +=stept*(aux[0][1]-aux[0][0])/stepx;
    rhou[0]  +=stept*(aux[1][1]-aux[1][0])/stepx;
    rhov[0]  +=stept*(aux[2][1]-aux[2][0])/stepx;
    rhow[0]  +=stept*(aux[3][1]-aux[3][0])/stepx;
    Energy[0]+=stept*(aux[4][1]-aux[4][0])/stepx;
    for(int j=1;j<numberx-1;++j){
      rho[j]   +=stept*(aux[0][j+1]-aux[0][j-1])/(2*stepx);
      rhou[j]  +=stept*(aux[1][j+1]-aux[1][j-1])/(2*stepx);
      rhov[j]  +=stept*(aux[2][j+1]-aux[2][j-1])/(2*stepx);
      rhow[j]  +=stept*(aux[3][j+1]-aux[3][j-1])/(2*stepx);
      Energy[j]+=stept*(aux[4][j+1]-aux[4][j-1])/(2*stepx);
    }
    rho[numberx-1] +=stept*(aux[0][numberx-1]-aux[0][numberx-2])/stepx;
    rhou[numberx-1]+=stept*(aux[1][numberx-1]-aux[1][numberx-2])/stepx;
    rhov[numberx-1]+=stept*(aux[2][numberx-1]-aux[2][numberx-2])/stepx;
    rhow[numberx-1]+=stept*(aux[3][numberx-1]-aux[3][numberx-2])/stepx;
    Energy[numberx-1]+=stept*(aux[4][numberx-1]-aux[4][numberx-2])/stepx;
*/
    //for(int j=0;j<numberx;++j)
      //aaa<<i<<"\t"<<j*stepx<<"\t"<<rho[j]<<"\t"<<rhou[j]<<"\t"<<rhov[j]<<"\t"<<rhow[j]<<"\t"<<Energy[j]<<"\t"<<Pressure(j)<<endl;

  }

}

//___________________________________________________________
long double* Fluid::GetF(int ix) const {

 // returns F(q) at cell ix

  if(ix<0)ix=0;if(ix>numberx-1)ix=numberx-1;

  long double *f=new long double[5];

  long double u=rhou[ix]/rho[ix];// u speed
  long double v=rhov[ix]/rho[ix];// v speed
  long double w=rhow[ix]/rho[ix];// w speed

  /*  f vector  */
  f[0]=rho[ix]*u;
  f[1]=rho[ix]*u*u+Pressure(ix);
  f[2]=rho[ix]*u*v;
  f[3]=rho[ix]*u*w;
  f[4]=(Energy[ix]+Pressure(ix))*u;

  return f;

}

//___________________________________________________________
long double* Fluid::GetG(int ix) const {

 // returns G(q) at cell ix

  long double *g=new long double[5];

  long double u=rhou[ix]/rho[ix];// u speed
  long double v=rhov[ix]/rho[ix];// v speed
  long double w=rhow[ix]/rho[ix];// w speed

  /*  g vector  */
  g[0]=rho[ix]*v;
  g[1]=rho[ix]*u*v;
  g[2]=rho[ix]*v*v+Pressure(ix);
  g[3]=rho[ix]*v*w;
  g[4]=(Energy[ix]+Pressure(ix))*v;

  return g;

}

//___________________________________________________________
long double* Fluid::GetH(int ix) const {

 // returns H(q) at cell ix

  long double *h=new long double[5];

  long double u=rhou[ix]/rho[ix];// u speed
  long double v=rhov[ix]/rho[ix];// v speed
  long double w=rhow[ix]/rho[ix];// w speed

  /*  h vector  */
  h[0]=rho[ix]*w;
  h[1]=rho[ix]*u*w;
  h[2]=rho[ix]*v*w;
  h[3]=rho[ix]*w*w+Pressure(ix);
  h[4]=(Energy[ix]+Pressure(ix))*w;

  return h;

}

//___________________________________________________________
long double Fluid::Pressure(int ix) const {

 // returns pressure at cell ix

  long double u=rhou[ix]/rho[ix];// u speed
  long double v=rhov[ix]/rho[ix];// v speed
  long double w=rhow[ix]/rho[ix];// w speed

  long double pressure=(addiabatic-1)*(Energy[ix]-rho[ix]/2*(u*u+v*v+w*w));
  return pressure;

}

//___________________________________________________________
long double* Fluid::RoeAverage(int ix) const {

 // returns row average at cell ix

  long double *q=new long double[5];

  long double pl=rho[ix], pr=rho[ix+1];// densities
  long double El=Energy[ix], Er=Energy[ix+1];// energies
  long double Pl=Pressure(ix), Pr=Pressure(ix+1);// pressures

  long double ul=rhou[ix]/rho[ix], ur=rhou[ix+1]/rho[ix+1];// u speeds
  long double vl=rhov[ix]/rho[ix], vr=rhov[ix+1]/rho[ix+1];// v speeds
  long double wl=rhow[ix]/rho[ix], wr=rhow[ix+1]/rho[ix+1];// w speeds

  /*  roe average  */
  q[0]=(pl+pr)/2;
  q[1]=(sqrt(pl)*ul+sqrt(pr)*ur)/(sqrt(pl)+sqrt(pr));
  q[2]=(sqrt(pl)*vl+sqrt(pr)*vr)/(sqrt(pl)+sqrt(pr));
  q[3]=(sqrt(pl)*wl+sqrt(pr)*wr)/(sqrt(pl)+sqrt(pr));
  q[4]=((El+Pl)/sqrt(pl)+(Er+Pr)/sqrt(pr))/(sqrt(pl)+sqrt(pr));

  /*  checks physical significate  */
  long double Emin=0.05/(addiabatic-1)+q[0]*(q[1]*q[1]+q[2]*q[2]+q[3]*q[3])/2;
  if(q[0]<0.05)q[0]=0.05;
  if(q[4]<Emin)q[4]=Emin;

  return q;

}

//___________________________________________________________
long double* Fluid::Eigenvalues(long double *roe) const {

 // returns wavespeeds of roe average

  long double *lam=new long double[5];

  long double pressure=(addiabatic-1)/addiabatic*roe[0]*(roe[4]-roe[1]*roe[1]/2);// pressure
  long double c=sqrt(addiabatic*pressure/roe[0]);// speed sound

  lam[0]=roe[1]-c;
  lam[1]=roe[1];
  lam[2]=roe[1];
  lam[3]=roe[1];
  lam[4]=roe[1]+c;

  return lam;

}

//___________________________________________________________
long double** Fluid::Eigenvectors(long double *roe) const {

 // returns eigenvectors of F(roe)

  long double **r=new long double*[5];
  for(int i=0;i<5;++i)r[i]=new long double[5];

  long double pressure=(addiabatic-1)/addiabatic*roe[0]*(roe[4]-roe[1]*roe[1]/2);// pressure
  long double c=sqrt(addiabatic*pressure/roe[0]);// speed sound

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
long double* Fluid::GetAlpha(int ix, long double *roe) const {

 // return alpha coefficents of cell ix

  long double *alpha=new long double[5];
  long double *deltaq=new long double[5];

  /*  shift in q  */
  deltaq[0]=rho[ix+1]   -rho[ix];
  deltaq[1]=rhou[ix+1]  -rhou[ix];
  deltaq[2]=rhov[ix+1]  -rhov[ix];
  deltaq[3]=rhow[ix+1]  -rhow[ix];
  deltaq[4]=Energy[ix+1]-Energy[ix];

  long double pressure=(addiabatic-1)/addiabatic*roe[0]*(roe[4]-roe[1]*roe[1]/2);// pressure
  long double c=sqrt(addiabatic*pressure/roe[0]);// speed of sound

  long double u=roe[1], v=roe[2], w=roe[3], H=roe[4];

  long double ct=c*(2*H-u*u-v*v-w*w);
  long double norm=u*u+v*v+w*w;

  /*  calculates inverse of eigenvector matrix  */
  long double **m=new long double*[5];
  for(int i=0;i<5;++i)m[i]=new long double[5];

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
long double** Fluid::GetWaves(long double *alpha, long double **r) const {

 // returns waves at cell ix

  long double **w=new long double*[5];
  for(int i=0;i<5;++i)w[i]=new long double[5];

  for(int i=0;i<5;++i)
    for(int j=0;j<5;++j)
      w[i][j]=alpha[j]*r[i][j];

  return w;

}

//___________________________________________________________
long double* Fluid::GetFb(int ix) const {

 // returns F-(q) at cell ix

  long double *roe=RoeAverage(ix);
  long double *lambda=Eigenvalues(roe);
  long double **r=Eigenvectors(roe); 
  long double *alpha=GetAlpha(ix,roe);
  long double **waves=GetWaves(alpha,r);

  long double *result=new long double[5];
  for(int i=0;i<5;++i)result[i]=0;

  for(int i=0;i<5;++i)if(lambda[i]<0)
    for(int j=0;j<5;++j)
      result[j]+=lambda[i]*waves[i][j];

  return result;

}

#include<iomanip>
//___________________________________________________________
long double* Fluid::GetFt(int ix) const {

 // returns F+(q) at cell ix

  long double *roe=RoeAverage(ix);
  long double *lambda=Eigenvalues(roe);
  long double **r=Eigenvectors(roe); 
  long double *alpha=GetAlpha(ix,roe);
  long double **waves=GetWaves(alpha,r);

  long double *result=new long double[5];
  for(int i=0;i<5;++i)result[i]=0;

  for(int i=0;i<5;++i)if(lambda[i]>0)
    for(int j=0;j<5;++j)
      result[j]+=lambda[i]*waves[i][j];

  return result;

}

//___________________________________________________________
long double* Fluid::GetEnergy() const {

 // returns energy

  long double *ret=new long double[numberx];
  for(int i=0;i<numberx;++i)
    ret[i]=Energy[i];
  return ret;

}

//___________________________________________________________
long double* Fluid::GetRho() const {

 // returns density

  long double *ret=new long double[numberx];
  for(int i=0;i<numberx;++i)
    ret[i]=rho[i];
  return ret;

}

//___________________________________________________________
long double* Fluid::GetRhou() const {

 // returns momentum in x

  long double *ret=new long double[numberx];
  for(int i=0;i<numberx;++i)
    ret[i]=rhou[i];
  return ret;

}

//___________________________________________________________
long double* Fluid::GetRhov() const {

 // returns momentum in y

  long double *ret=new long double[numberx];
  for(int i=0;i<numberx;++i)
    ret[i]=rhov[i];
  return ret;

}

//___________________________________________________________
long double* Fluid::GetRhow() const {

 // returns momentum in z

  long double *ret=new long double[numberx];
  for(int i=0;i<numberx;++i)
    ret[i]=rhow[i];
  return ret;

}

//___________________________________________________________
void Fluid::SetEnergy(long double *E){

 // sets energy to E

  for(int i=0;i<numberx;++i)
    Energy[i]=E[i];

}

//___________________________________________________________
void Fluid::SetRho(long double *p){

 // sets rho to p

  for(int i=0;i<numberx;++i)
    rho[i]=p[i];

}

//___________________________________________________________
void Fluid::SetRhou(long double *u){

 // sets momentum in x to u

  for(int i=0;i<numberx;++i)
    rhou[i]=u[i];

}

//___________________________________________________________
void Fluid::SetRhov(long double *v){

 // sets momentum in y to v

  for(int i=0;i<numberx;++i)
    rhov[i]=v[i];

}

//___________________________________________________________
void Fluid::SetRhow(long double *w){

 // sets momentum in z to w

  for(int i=0;i<numberx;++i)
    rhow[i]=w[i];

}









