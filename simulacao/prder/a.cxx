#include<iostream>
#include<cstdlib>
#include<cstdio>
#include<cmath>
using namespace std;

double fff(double x){

  double a=-18./200*(1.4-1)*2.29e12/sqrt(200*M_PI/9);
  return a*(x-50)*exp(-9./200*pow(x-50,2));
  
}

int main(){

int *N=new int[30];
N[0]=10;      N[1]=20;      N[2]=30;      N[3]=50;      N[4]=75;
N[5]=100;     N[6]=200;     N[7]=300;     N[8]=500;     N[9]=750;
N[10]=1000;   N[11]=2000;   N[12]=3000;   N[13]=5000;   N[14]=7500;
N[15]=10000;  N[16]=20000;  N[17]=30000;  N[18]=50000;  N[19]=75000;
N[20]=100000; N[21]=200000; N[22]=300000; N[23]=500000; N[24]=750000;
N[25]=1000000;N[26]=2000000;N[27]=3000000;N[28]=5000000;N[29]=7500000;

const char **C=new const char*[30];
C[0]="10.run";      C[1]="20.run";      C[2]="30.run";      C[3]="50.run";      C[4]="75.run";
C[5]="100.run";     C[6]="200.run";     C[7]="300.run";     C[8]="500.run";     C[9]="750.run";
C[10]="1000.run";   C[11]="2000.run";   C[12]="3000.run";   C[13]="5000.run";   C[14]="7500.run";
C[15]="10000.run";  C[16]="20000.run";  C[17]="30000.run";  C[18]="50000.run";  C[19]="75000.run";
C[20]="100000.run"; C[21]="200000.run"; C[22]="300000.run"; C[23]="500000.run"; C[24]="750000.run";
C[25]="1000000.run";C[26]="2000000.run";C[27]="3000000.run";C[28]="5000000.run";C[29]="7500000.run";

for(int j=0;j<30;++j){
  FILE *aaa=fopen(C[j],"r+");
  double *x=new double[N[j]];
  double *y=new double[N[j]];
  for(int i=0;i<N[j];++i){
    fscanf(aaa,"%lf",&x[i]);
    fscanf(aaa,"%lf",&y[i]);
  }
  double chi=0;int u=0;
  for(int i=0;i<N[j];++i)if(y[i]>1){
    double teor=fff(x[i]);++u;
    chi+=pow((teor-y[i])/y[i],2);
  }
  cout<<N[j]<<"\t"<<chi/u<<endl;
}
}

