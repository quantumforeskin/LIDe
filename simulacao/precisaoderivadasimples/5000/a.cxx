#include<iostream>
#include<cstdio>
#include<cmath>
using namespace std;

double f(double x, double v){

  double fff=-9.86664*(x-25)*exp(-(x-25)*(x-25)*9/200)*pow(10,9);
  return v/fff;

}

int main(){

  FILE *aaa=fopen("evolution3_5000.run","r+");
  double x, y;
  for(int i=0;i<5000;++i){
    fscanf(aaa,"%lf",&x);fscanf(aaa,"%lf",&x);
    fscanf(aaa,"%lf",&x);fscanf(aaa,"%lf",&x);
    fscanf(aaa,"%lf",&x);fscanf(aaa,"%lf",&x);
    fscanf(aaa,"%lf",&x);fscanf(aaa,"%lf",&x);
  }
  fscanf(aaa,"%lf",&x);fscanf(aaa,"%lf",&x);fscanf(aaa,"%lf",&x);
  fscanf(aaa,"%lf",&x);fscanf(aaa,"%lf",&x);fscanf(aaa,"%lf",&x);
  fscanf(aaa,"%lf",&x);fscanf(aaa,"%lf",&x);fscanf(aaa,"%lf",&y);
  fscanf(aaa,"%lf",&x);fscanf(aaa,"%lf",&x);
  cout<<f(y,x)<<endl;

}
