#include"Fluid.h"

//___________________________________________________________
Fluid::Fluid(ld adb, int nx, ld bx, ld ux){

 //  constructor, sets up a grid in x from bx to ux, with nx cells. It 
 // also sets up addiabatic to adb

  /*  Sets static variables  */
  numberx=nx;
  addiabatic=adb;
  upperx=ux;
  bottomx=bx;

  /*  Allocates memory for internal variables  */
  rho=new ld[numberx];
  for(int i=0;i<numberx;++i)rho[i]=0;

  rhou=new ld[numberx];
  for(int i=0;i<numberx;++i)rhou[i]=0;

  rhov=new ld[numberx];
  for(int i=0;i<numberx;++i)rhov[i]=0;

  rhow=new ld[numberx];
  for(int i=0;i<numberx;++i)rhow[i]=0;

  Energy=new ld[numberx];
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
void Fluid::Evolution1(const char *output, ld time, int nsteps){

 //  runs evolution of fluid through time with nsteps steps in time, 
 // writing the result into output, using simple numeric method

  ofstream aaa(output);
  ld stept=time/(nsteps-1), stepx=(upperx-bottomx)/(numberx-1);
/*
  for(int i=0;i<numberx;++i)// writes initial condition
    aaa<<"0\t"<<i*stepx<<"\t"<<rho[i]<<"\t"<<rhou[i]<<"\t"<<rhov[i]<<"\t"<<rhow[i]<<"\t"<<Energy[i]<<"\t"<<Pressure(i)<<endl;
*/
  ld **aux=new ld*[5];
  for(int i=0;i<5;++i)aux[i]=new ld[numberx];

  for(int i=1;i<nsteps;++i){// runs evolution in time
    for(int j=0;j<numberx;++j){// calculates evolution in the spacial grid

      ld *f1=GetF(j), *f2=GetF(j+1);

      aux[0][j]=rho[j]   -stept/stepx*(f1[0]-f2[0]);
      aux[1][j]=rhou[j]  -stept/stepx*(f1[1]-f2[1]);
      aux[2][j]=rhov[j]  -stept/stepx*(f1[2]-f2[2]);
      aux[3][j]=rhow[j]  -stept/stepx*(f1[3]-f2[3]);
      aux[4][j]=Energy[j]-stept/stepx*(f1[4]-f2[4]);

      aaa<<j*stepx<<"\t"<<(f2[1]-f1[1])/stepx<<endl;

    }

    for(int j=0;j<numberx;++j){
      rho[j]   =aux[0][j];
      rhou[j]  =aux[1][j];
      rhov[j]  =aux[2][j];
      rhow[j]  =aux[3][j];
      Energy[j]=aux[4][j];
    }
/*
    for(int j=0;j<numberx;++j)// writes output into file
      aaa<<i<<"\t"<<j*stepx<<"\t"<<rho[j]<<"\t"<<rhou[j]<<"\t"<<rhov[j]<<"\t"<<rhow[j]<<"\t"<<Energy[j]<<"\t"<<Pressure(j)<<endl;
*/
  }

}

//___________________________________________________________
void Fluid::Evolution2(const char *output, ld time, int nsteps){

 //  runs evolution of fluid through time with nsteps steps in time, 
 // writing the result into output, using optimal numeric method

  ofstream aaa(output);
  ld stepx=(upperx-bottomx)/(numberx-1), stept=time/(nsteps-1);
/*
  // writes initial condition
  aaa<<"0\t0\t0\t"<<rho[0]<<"\t"<<rhou[0]<<"\t"<<rhov[0]<<"\t"<<rhow[0]<<"\t"<<Energy[0]<<"\t"<<Pressure(0)<<endl;
  for(int i=1;i<numberx-1;++i)if(fabs(Energy[i]-Energy[i-1])>1||fabs(rho[i]-rho[i-1])>1e-3||fabs(rhou[i]-rhou[i-1])>0.1)
    aaa<<"0\t0\t"<<i*stepx<<"\t"<<rho[i]<<"\t"<<rhou[i]<<"\t"<<rhov[i]<<"\t"<<rhow[i]<<"\t"<<Energy[i]<<"\t"<<Pressure(i)<<endl;
  aaa<<"0\t0\t"<<(numberx-1)*stepx<<"\t"<<rho[numberx-1]<<"\t"<<rhou[numberx-1]<<"\t"<<rhov[numberx-1]<<"\t"<<rhow[numberx-1]<<"\t"<<Energy[numberx-1]<<"\t"<<Pressure(numberx-1)<<endl;
*/
  double **aux=new double*[5];// auxiliary vector
  for(int i=0;i<5;++i)aux[i]=new double[numberx+1];

  int up=clock();
  for(int i=1;i<nsteps;++i){// runs evolution

    for(int j=0;j<numberx;++j){// calculates f(q) on grid
      ld *g=GetF(j);
      aux[0][j]=g[0];
      aux[1][j]=g[1];
      aux[2][j]=g[2];
      aux[3][j]=g[3];
      aux[4][j]=g[4];
    }

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
      if(rho[j]<0.05)rho[j]=0.05;
      if(Pressure(j)<0.05)Energy[j]=(rhou[j]*rhou[j]+rhov[j]*rhov[j]+rhow[j]*rhow[j])/(2*rho[j])+0.05/(addiabatic-1);
    }
    rho[numberx-1]   +=stept*(aux[0][numberx-1]-aux[0][numberx-2])/stepx;
    rhou[numberx-1]  +=stept*(aux[1][numberx-1]-aux[1][numberx-2])/stepx;
    rhov[numberx-1]  +=stept*(aux[2][numberx-1]-aux[2][numberx-2])/stepx;
    rhow[numberx-1]  +=stept*(aux[3][numberx-1]-aux[3][numberx-2])/stepx;
    Energy[numberx-1]+=stept*(aux[4][numberx-1]-aux[4][numberx-2])/stepx;
/*
    if(nsteps>50){
      int avb=nsteps/50;
      if(i%avb==0){
        aaa<<i<<"\t"<<i*stept<<"\t0\t"<<rho[0]<<"\t"<<rhou[0]<<"\t"<<rhov[0]<<"\t"<<rhow[0]<<"\t"<<Energy[0]<<"\t"<<Pressure(0)<<endl;
        for(int j=0;j<numberx;++j)if(fabs(Energy[j]-Energy[j-1])>1||fabs(rho[j]-rho[j-1])>1e-3||fabs(rhou[j]-rhou[j-1])>0.1)
          aaa<<i<<"\t"<<i*stept<<"\t"<<j*stepx<<"\t"<<rho[j]<<"\t"<<rhou[j]<<"\t"<<rhov[j]<<"\t"<<rhow[j]<<"\t"<<Energy[j]<<"\t"<<Pressure(j)<<endl;
        aaa<<i<<"\t"<<i*stept<<"\t"<<numberx*stepx-stepx<<"\t"<<rho[numberx-1]<<"\t"<<rhou[numberx-1]<<"\t"<<rhov[numberx-1]<<"\t"<<rhow[numberx-1]<<"\t"<<Energy[numberx-1]<<"\t"<<Pressure(numberx-1)<<endl;
      }
    }
    else {
      aaa<<i<<"\t"<<i*stept<<"\t0\t"<<rho[0]<<"\t"<<rhou[0]<<"\t"<<rhov[0]<<"\t"<<rhow[0]<<"\t"<<Energy[0]<<"\t"<<Pressure(0)<<endl;
      for(int j=0;j<numberx;++j)if(fabs(Energy[j]-Energy[j-1])>1||fabs(rho[j]-rho[j-1])>1e-3||fabs(rhou[j]-rhou[j-1])>0.1)
        aaa<<i<<"\t"<<i*stept<<"\t"<<j*stepx<<"\t"<<rho[j]<<"\t"<<rhou[j]<<"\t"<<rhov[j]<<"\t"<<rhow[j]<<"\t"<<Energy[j]<<"\t"<<Pressure(j)<<endl;
      aaa<<i<<"\t"<<i*stept<<"\t"<<numberx*stepx-stepx<<"\t"<<rho[numberx-1]<<"\t"<<rhou[numberx-1]<<"\t"<<rhov[numberx-1]<<"\t"<<rhow[numberx-1]<<"\t"<<Energy[numberx-1]<<"\t"<<Pressure(numberx-1)<<endl;
    }
*/
  }
  int down=clock();
  //cout<<((double)(down-up))/CLOCKS_PER_SEC<<endl;

}

//___________________________________________________________
ld* Fluid::GetF(int ix) const {

 // returns F(q) at cell ix

  if(ix<0)ix=0;if(ix>numberx-1)ix=numberx-1;

  ld *f=new ld[5];

  ld u=rhou[ix]/rho[ix];// u speed
  ld v=rhov[ix]/rho[ix];// v speed
  ld w=rhow[ix]/rho[ix];// w speed

  /*  f vector  */
  f[0]=rho[ix]*u;
  f[1]=rho[ix]*u*u+Pressure(ix);
  f[2]=rho[ix]*u*v;
  f[3]=rho[ix]*u*w;
  f[4]=(Energy[ix]+Pressure(ix))*u;

  return f;

}

//___________________________________________________________
ld* Fluid::GetG(int ix) const {

 // returns G(q) at cell ix

  ld *g=new ld[5];

  ld u=rhou[ix]/rho[ix];// u speed
  ld v=rhov[ix]/rho[ix];// v speed
  ld w=rhow[ix]/rho[ix];// w speed

  /*  g vector  */
  g[0]=rho[ix]*v;
  g[1]=rho[ix]*u*v;
  g[2]=rho[ix]*v*v+Pressure(ix);
  g[3]=rho[ix]*v*w;
  g[4]=(Energy[ix]+Pressure(ix))*v;

  return g;

}

//___________________________________________________________
ld* Fluid::GetH(int ix) const {

 // returns H(q) at cell ix

  ld *h=new ld[5];

  ld u=rhou[ix]/rho[ix];// u speed
  ld v=rhov[ix]/rho[ix];// v speed
  ld w=rhow[ix]/rho[ix];// w speed

  /*  h vector  */
  h[0]=rho[ix]*w;
  h[1]=rho[ix]*u*w;
  h[2]=rho[ix]*v*w;
  h[3]=rho[ix]*w*w+Pressure(ix);
  h[4]=(Energy[ix]+Pressure(ix))*w;

  return h;

}

//___________________________________________________________
ld Fluid::Pressure(int ix) const {

 // returns pressure at cell ix

  ld u=rhou[ix]/rho[ix];// u speed
  ld v=rhov[ix]/rho[ix];// v speed
  ld w=rhow[ix]/rho[ix];// w speed

  ld pressure=(addiabatic-1)*(Energy[ix]-rho[ix]/2*(u*u+v*v+w*w));
  return pressure;

}

//___________________________________________________________
ld* Fluid::GetRoeAverage(int ix) const {

 // returns roe average at cell ix

  double pl=rho[ix],      pr=rho[ix+1];
  double ul=rhou[ix]/pl,  ur=rhou[ix+1]/pr;
  double vl=rhov[ix]/pl,  vr=rhov[ix+1]/pr;
  double wl=rhow[ix]/pl,  wr=rhow[ix+1]/pr;
  double el=Energy[ix],   er=Energy[ix+1];
  double Pl=Pressure(ix), Pr=Pressure(ix+1);
  double spl=sqrt(pl),    spr=sqrt(pr);

  ld *Q=new ld[5];
  Q[0]=(pl+pr)/2;
  Q[1]=(spl*ul+spr*ur)/(spl+spr);
  Q[2]=(spl*vl+spr*vr)/(spl+spr);
  Q[3]=(spl*wl+spr*wr)/(spl+spr);
  Q[4]=((el+Pl)/spl+(er+Pr)/spr)/(spl+spr);

  return Q;

}

//___________________________________________________________
ld* Fluid::GetLambda(ld *Q) const {

 // returns eigenvalues with roe average Q

  double u=Q[1], p=Q[0], 
}

//___________________________________________________________
ld* Fluid::GetEnergy() const {

 // returns energy

  ld *ret=new ld[numberx];
  for(int i=0;i<numberx;++i)
    ret[i]=Energy[i];
  return ret;

}

//___________________________________________________________
ld* Fluid::GetRho() const {

 // returns density

  ld *ret=new ld[numberx];
  for(int i=0;i<numberx;++i)
    ret[i]=rho[i];
  return ret;

}

//___________________________________________________________
ld* Fluid::GetRhou() const {

 // returns momentum in x

  ld *ret=new ld[numberx];
  for(int i=0;i<numberx;++i)
    ret[i]=rhou[i];
  return ret;

}

//___________________________________________________________
ld* Fluid::GetRhov() const {

 // returns momentum in y

  ld *ret=new ld[numberx];
  for(int i=0;i<numberx;++i)
    ret[i]=rhov[i];
  return ret;

}

//___________________________________________________________
ld* Fluid::GetRhow() const {

 // returns momentum in z

  ld *ret=new ld[numberx];
  for(int i=0;i<numberx;++i)
    ret[i]=rhow[i];
  return ret;

}

//___________________________________________________________
void Fluid::SetEnergy(ld *E){

 // sets energy to E

  for(int i=0;i<numberx;++i)
    Energy[i]=E[i];

}

//___________________________________________________________
void Fluid::SetRho(ld *p){

 // sets rho to p

  for(int i=0;i<numberx;++i)
    rho[i]=p[i];

}

//___________________________________________________________
void Fluid::SetRhou(ld *u){

 // sets momentum in x to u

  for(int i=0;i<numberx;++i)
    rhou[i]=u[i];

}

//___________________________________________________________
void Fluid::SetRhov(ld *v){

 // sets momentum in y to v

  for(int i=0;i<numberx;++i)
    rhov[i]=v[i];

}

//___________________________________________________________
void Fluid::SetRhow(ld *w){

 // sets momentum in z to w

  for(int i=0;i<numberx;++i)
    rhow[i]=w[i];

}









