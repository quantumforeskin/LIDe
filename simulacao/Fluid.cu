#include"Fluid.h"

//___________________________________________________________
__global__ void callcuda(float *cad, float *ct, float *cx, int *cnumber, float *caux1, float *caux2, float *caux3, float *caux4, float *caux5, float *crho, float *crhou, float *crhov, float *crhow, float *cenergy){

  int index=blockIdx.x*blockDim.x+threadIdx.x;
  int undindex=index-1, topindex=index+1;
  if(topindex>cnumber[0]-1)topindex=cnumber[0]-1;
  if(undindex<0)undindex=0;

  float uu=crhou[undindex]/crho[undindex];
  float ut=crhou[topindex]/crho[topindex];
  float undpressure=(cad[0]-1)*(cenergy[undindex]-crho[undindex]*uu*uu/2);
  float toppressure=(cad[0]-1)*(cenergy[topindex]-crho[topindex]*ut*ut/2);

  /*  calculates f1_r and f1_l  */
  float f10=crhou[undindex];
  float f11=crhou[undindex]*uu+undpressure;
  float f12=crhov[undindex]*uu;
  float f13=crhow[undindex]*uu;
  float f14=(cenergy[undindex]+undpressure)*uu;

  float f20=crhou[topindex];
  float f21=crhou[topindex]*ut+toppressure;
  float f22=crhov[topindex]*ut;
  float f23=crhow[topindex]*ut;
  float f24=(cenergy[topindex]+toppressure)*ut;

  /*  calculates evolution  */
  caux1[index]=crho[index]   -ct[0]/(cx[0])*(f10-f20);
  caux2[index]=crhou[index]  -ct[0]/(cx[0])*(f11-f21);
  caux3[index]=crhov[index]  -ct[0]/(cx[0])*(f12-f22);
  caux4[index]=crhow[index]  -ct[0]/(cx[0])*(f13-f23);
  caux5[index]=cenergy[index]-ct[0]/(cx[0])*(f14-f24);

}

//___________________________________________________________
void Fluid::Evolution1c(const char *output, long double time, int nsteps){

 //  runs evolution of fluid through time with nsteps steps in time, 
 // writing the result into output, using simple numeric method

  ofstream aaa(output);
  float stept=time/(nsteps-1), stepx=(upperx-bottomx)/(numberx-1);

  for(int i=0;i<numberx;++i)if(i%20==0)// writes initial condition
    aaa<<"0\t"<<i*stepx<<"\t"<<rho[i]<<"\t"<<rhou[i]<<"\t"<<rhov[i]<<"\t"<<rhow[i]<<"\t"<<Energy[i]<<"\t"<<Pressure(i)<<endl;

  /*  Allocates memory in CUDA and copies static variables  */

  int *cnumber;
  float *cad, *ct, *cx;
  float *crho, *crhou,*crhov,*crhow,*cenergy;
  float *caux1,*caux2,*caux3,*caux4,*caux5;

  if(cudaSuccess!=cudaMalloc((void**)&cnumber,sizeof(int)))cout<<"memoria"<<endl;
  if(cudaSuccess!=cudaMemcpy(cnumber,&numberx,sizeof(int),cudaMemcpyHostToDevice))cout<<"copia"<<endl;
  if(cudaSuccess!=cudaMalloc((void**)&cad,       sizeof(float)))cout<<"memoria"<<endl;
  if(cudaSuccess!=cudaMemcpy(cnumber,&addiabatic,sizeof(float),cudaMemcpyHostToDevice))cout<<"copia"<<endl;
  if(cudaSuccess!=cudaMalloc((void**)&ct,sizeof(float)))cout<<"memoria"<<endl;
  if(cudaSuccess!=cudaMemcpy(ct,&stept,  sizeof(float),cudaMemcpyHostToDevice))cout<<"copia"<<endl;
  if(cudaSuccess!=cudaMalloc((void**)&cx,sizeof(float)))cout<<"memoria"<<endl;
  if(cudaSuccess!=cudaMemcpy(cx,&stepx,  sizeof(float),cudaMemcpyHostToDevice))cout<<"copia"<<endl;

  if(cudaSuccess!=cudaMalloc((void**)&caux1  ,numberx*sizeof(float)))cout<<"memoria"<<endl;
  if(cudaSuccess!=cudaMalloc((void**)&caux2  ,numberx*sizeof(float)))cout<<"memoria"<<endl;
  if(cudaSuccess!=cudaMalloc((void**)&caux3  ,numberx*sizeof(float)))cout<<"memoria"<<endl;
  if(cudaSuccess!=cudaMalloc((void**)&caux4  ,numberx*sizeof(float)))cout<<"memoria"<<endl;
  if(cudaSuccess!=cudaMalloc((void**)&caux5  ,numberx*sizeof(float)))cout<<"memoria"<<endl;
  if(cudaSuccess!=cudaMalloc((void**)&crho   ,numberx*sizeof(float)))cout<<"memoria"<<endl;
  if(cudaSuccess!=cudaMalloc((void**)&crhou  ,numberx*sizeof(float)))cout<<"memoria"<<endl;
  if(cudaSuccess!=cudaMalloc((void**)&crhov  ,numberx*sizeof(float)))cout<<"memoria"<<endl;
  if(cudaSuccess!=cudaMalloc((void**)&crhow  ,numberx*sizeof(float)))cout<<"memoria"<<endl;
  if(cudaSuccess!=cudaMalloc((void**)&cenergy,numberx*sizeof(float)))cout<<"memoria"<<endl;

  /*  Runs evolution  */
  for(int i=1;i<nsteps;++i){

    /*  copies current status to cuda  */
    if(cudaSuccess!=cudaMemcpy(crho   ,rho   ,numberx*sizeof(float),cudaMemcpyHostToDevice))cout<<"copytodevice"<<endl;
    if(cudaSuccess!=cudaMemcpy(crhou  ,rhou  ,numberx*sizeof(float),cudaMemcpyHostToDevice))cout<<"copytodevice"<<endl;
    if(cudaSuccess!=cudaMemcpy(crhov  ,rhov  ,numberx*sizeof(float),cudaMemcpyHostToDevice))cout<<"copytodevice"<<endl;
    if(cudaSuccess!=cudaMemcpy(crhow  ,rhow  ,numberx*sizeof(float),cudaMemcpyHostToDevice))cout<<"copytodevice"<<endl;
    if(cudaSuccess!=cudaMemcpy(cenergy,Energy,numberx*sizeof(float),cudaMemcpyHostToDevice))cout<<"copytodevice"<<endl;

    //int threadsize=64;// calculates evolution
    callcuda<<<1,numberx>>>(cad,ct,cx,cnumber,caux1,caux2,caux3,caux4,caux5,crho,crhou,crhov,crhow,cenergy);

    /*  copies result into host  */
    if(cudaSuccess!=cudaMemcpy(rho,&caux1,numberx*sizeof(float),cudaMemcpyDeviceToHost))cout<<"copyfromdevice"<<endl;
    if(cudaSuccess!=cudaMemcpy(rhou,&caux2,numberx*sizeof(float),cudaMemcpyDeviceToHost))cout<<"copyfromdevice"<<endl;
    if(cudaSuccess!=cudaMemcpy(rhov,&caux3,numberx*sizeof(float),cudaMemcpyDeviceToHost))cout<<"copyfromdevice"<<endl;
    if(cudaSuccess!=cudaMemcpy(rhow,&caux4,numberx*sizeof(float),cudaMemcpyDeviceToHost))cout<<"copyfromdevice"<<endl;
    if(cudaSuccess!=cudaMemcpy(Energy,&caux5,numberx*sizeof(float),cudaMemcpyDeviceToHost))cout<<"copyfromdevice"<<endl;

    for(int j=0;j<numberx;++j)if(j%20==0)// writes output
      aaa<<i*stept<<"\t"<<j*stepx<<"\t"<<rho[j]<<"\t"<<rhou[j]<<"\t"<<rhov[j]<<"\t"<<rhow[j]<<"\t"<<Energy[j]<<"\t"<<Pressure(j)<<endl;

  }

}

