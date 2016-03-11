#ifndef __FLUID_
#define __FLUID_

#include<iostream>
#include<fstream>
#include<cmath>
using namespace std;

class Fluid {

 public:

  /*  CONSTRUCTOR AND DESTRUCTOR  */
  Fluid(double,int,double,double);
  ~Fluid();

  /*  TIME EVOLUTION  */
  void Evolution1(const char*,double,int);
  void Evolution2(const char*,double,int);
  void Evolution1c(const char*,double,int);

  /*  PHYSICAL VARIABLES  */
  double* GetF(int) const;
  double* GetG(int) const;
  double* GetH(int) const;
  double Pressure(int) const;
  double* RoeAverage(int) const;
  double* Eigenvalues(double*) const;
  double** Eigenvectors(double*) const;
  double* GetAlpha(int,double*) const;
  double** GetWaves(double*,double**) const;
  double* GetFb(int) const;
  double* GetFt(int) const;

  /*  INTERNAL VARIABLE FUNCTIONS  */
  double* GetEnergy() const;
  double* GetRho() const;
  double* GetRhou() const;
  double* GetRhov() const;
  double* GetRhow() const;
  void SetEnergy(double*);
  void SetRho(double*);
  void SetRhou(double*);
  void SetRhov(double*);
  void SetRhow(double*);

 private:

  double addiabatic;		// addiabatic constant
  int numberx;			// number of cells in x
  double bottomx, upperx;	// limits in x
  double *rho;			// density in cell
  double *rhou, *rhov, *rhow;	// momentum in cell
  double *Energy;		// energy in cell

};

#endif
