#include"Fluid.h"

__global__ void cudacop(float *c1, float *c2, float *c3, float *c4, float *c5){

  int index=threadIdx.x;

  c1[index]=0;
  c2[index]=0;
  c3[index]=0;
  c4[index]=0;
  c5[index]=0;

}

void Fluid::Evolution1c(const char *output, ld time, int steps){

  float *caux1, *caux2, *caux3, *caux4, *caux5;
  float  *aux1,  *aux2,  *aux3,  *aux4,  *aux5;

  aux1=new float[numberx];aux2=new float[numberx];
  aux3=new float[numberx];aux4=new float[numberx];
  aux5=new float[numberx];

  if(cudaSuccess!=cudaMalloc((void**)&caux1,numberx*sizeof(float)))cout<<"oi"<<endl;
  if(cudaSuccess!=cudaMalloc((void**)&caux2,numberx*sizeof(float)))cout<<"oi"<<endl;
  if(cudaSuccess!=cudaMalloc((void**)&caux3,numberx*sizeof(float)))cout<<"oi"<<endl;
  if(cudaSuccess!=cudaMalloc((void**)&caux4,numberx*sizeof(float)))cout<<"oi"<<endl;
  if(cudaSuccess!=cudaMalloc((void**)&caux5,numberx*sizeof(float)))cout<<"oi"<<endl;

  cudacop<<<1,numberx>>>(caux1,caux2,caux3,caux4,caux5);

  if(cudaSuccess!=cudaMemcpy(aux1,caux1,numberx*sizeof(float),cudaMemcpyDeviceToHost))cout<<"oi"<<endl;
  if(cudaSuccess!=cudaMemcpy(aux2,caux2,numberx*sizeof(float),cudaMemcpyDeviceToHost))cout<<"oi"<<endl;
  if(cudaSuccess!=cudaMemcpy(aux3,caux3,numberx*sizeof(float),cudaMemcpyDeviceToHost))cout<<"oi"<<endl;
  if(cudaSuccess!=cudaMemcpy(aux4,caux4,numberx*sizeof(float),cudaMemcpyDeviceToHost))cout<<"oi"<<endl;
  if(cudaSuccess!=cudaMemcpy(aux5,caux5,numberx*sizeof(float),cudaMemcpyDeviceToHost))cout<<"oi"<<endl;

  cout<<aux1[0]<<" "<<aux1[1]<<endl;
  cout<<aux2[0]<<" "<<aux2[1]<<endl;
  cout<<aux3[0]<<" "<<aux3[1]<<endl;
  cout<<aux4[0]<<" "<<aux4[1]<<endl;
  cout<<aux5[0]<<" "<<aux5[1]<<endl;

}
