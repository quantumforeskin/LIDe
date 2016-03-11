#include"Fluid.h"

//___________________________________________________________
__global__ void callcuda(float *cad, float *ct, float *cx, int *cnumber, float *caux1, float *caux2, float *caux3, float *caux4, float *caux5, float *crho, float *crhou, float *crhov, float *crhow, float *cenergy){

  int index=blockIdx.x*blockDim.x+threadIdx.x;
  int undindex=index-1, topindex=index+1;
  if(topindex>*cnumber-1)topindex=*cnumber-1;
  if(undindex<0)undindex=0;

  float uu=crhou[undindex]/crho[undindex];
  float ut=crhou[topindex]/crho[topindex];
  float undpressure=(*cad-1)*(cenergy[undindex]-crho[undindex]*uu*uu/2);
  float toppressure=(*cad-1)*(cenergy[topindex]-crho[topindex]*ut*ut/2);

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

  caux1[index]=crho[index]   -*ct/(*cx)*(f10-f20);
  caux2[index]=crhou[index]  -*ct/(*cx)*(f11-f21);
  caux3[index]=crhov[index]  -*ct/(*cx)*(f12-f22);
  caux4[index]=crhow[index]  -*ct/(*cx)*(f13-f23);
  caux5[index]=cenergy[index]-*ct/(*cx)*(f14-f24);

}

//___________________________________________________________
void Fluid::Evolution1c(const char *output, double time, int nsteps){

 //  runs evolution of fluid through time with nsteps steps in time, 
 // writing the result into output, using simple numeric method

  ofstream aaa(output);
  float stept=time/(nsteps-1), stepx=(upperx-bottomx)/(numberx-1);

  for(int i=0;i<numberx;++i)
    aaa<<"0\t"<<i*stepx<<"\t"<<rho[i]<<"\t"<<rhou[i]<<"\t"<<rhov[i]<<"\t"<<rhow[i]<<"\t"<<Energy[i]<<"\t"<<Pressure(i)<<endl;

  int *cnumber;
  float *cad, *ct, *cx;
  float *crho, *crhou,*crhov,*crhow,*cenergy;
  float *caux1,*caux2,*caux3,*caux4,*caux5;

  cudaMalloc((void**)&cnumber,sizeof(int));
  cudaMemcpy(cnumber,&numberx,sizeof(int),cudaMemcpyHostToDevice);
  cudaMalloc((void**)&cad,       sizeof(float));
  cudaMemcpy(cnumber,&addiabatic,sizeof(float),cudaMemcpyHostToDevice);
  cudaMalloc((void**)&ct,sizeof(float));
  cudaMemcpy(ct,&stept,  sizeof(float),cudaMemcpyHostToDevice);
  cudaMalloc((void**)&cx,sizeof(float));
  cudaMemcpy(cx,&stepx,  sizeof(float),cudaMemcpyHostToDevice);

  cudaMalloc((void**)&caux1  ,numberx*sizeof(float));
  cudaMalloc((void**)&caux2  ,numberx*sizeof(float));
  cudaMalloc((void**)&caux3  ,numberx*sizeof(float));
  cudaMalloc((void**)&caux4  ,numberx*sizeof(float));
  cudaMalloc((void**)&caux5  ,numberx*sizeof(float));
  cudaMalloc((void**)&crho   ,numberx*sizeof(float));
  cudaMalloc((void**)&crhou  ,numberx*sizeof(float));
  cudaMalloc((void**)&crhov  ,numberx*sizeof(float));
  cudaMalloc((void**)&crhow  ,numberx*sizeof(float));
  cudaMalloc((void**)&cenergy,numberx*sizeof(float));

  for(int i=1;i<nsteps;++i){

    cudaMemcpy(crho   ,&rho   ,numberx*sizeof(float),cudaMemcpyHostToDevice);
    cudaMemcpy(crhou  ,&rhou  ,numberx*sizeof(float),cudaMemcpyHostToDevice);
    cudaMemcpy(crhov  ,&rhov  ,numberx*sizeof(float),cudaMemcpyHostToDevice);
    cudaMemcpy(crhow  ,&rhow  ,numberx*sizeof(float),cudaMemcpyHostToDevice);
    cudaMemcpy(cenergy,&Energy,numberx*sizeof(float),cudaMemcpyHostToDevice);

    int threadsize=64;
    callcuda<<<numberx/threadsize+1,threadsize>>>(cad,ct,cx,cnumber,caux1,caux2,caux3,caux4,caux5,crho,crhou,crhov,crhow,cenergy);

    cudaMemcpy(rho,&caux1,numberx*sizeof(float),cudaMemcpyDeviceToHost);
    cudaMemcpy(rhou,&caux2,numberx*sizeof(float),cudaMemcpyDeviceToHost);
    cudaMemcpy(rhov,&caux3,numberx*sizeof(float),cudaMemcpyDeviceToHost);
    cudaMemcpy(rhow,&caux4,numberx*sizeof(float),cudaMemcpyDeviceToHost);
    cudaMemcpy(Energy,&caux5,numberx*sizeof(float),cudaMemcpyDeviceToHost);

    for(int j=0;j<numberx;++j)
      aaa<<i*stept<<"\t"<<j*stepx<<"\t"<<rho[j]<<"\t"<<rhou[j]<<"\t"<<rhov[j]<<"\t"<<rhow[j]<<"\t"<<Energy[j]<<"\t"<<Pressure(j)<<endl;

  }

}

