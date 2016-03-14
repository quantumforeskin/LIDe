#ifndef __FLUID_
#define __FLUID_

#include<iostream>
#include<fstream>
#include<iomanip>
#include<cmath>
using namespace std;

class Fluid {

 public:

  /*  CONSTRUCTOR AND DESTRUCTOR  */
  Fluid(long double,int,long double,long double);
  ~Fluid();

  /*  TIME EVOLUTION  */
  void Evolution1(const char*,long double,int);
  void Evolution2(const char*,long double,int);
  void Evolution3(const char*,long double,int);
  void Evolution1c(const char*,long double,int);

  /*  PHYSICAL VARIABLES  */
  long double* GetF(int) const;
  long double* GetG(int) const;
  long double* GetH(int) const;
  long double Pressure(int) const;
  long double* RoeAverage(int) const;
  long double* Eigenvalues(long double*) const;
  long double** Eigenvectors(long double*) const;
  long double* GetAlpha(int,long double*) const;
  long double** GetWaves(long double*,long double**) const;
  long double* GetFb(int) const;
  long double* GetFt(int) const;

  /*  INTERNAL VARIABLE FUNCTIONS  */
  long double* GetEnergy() const;
  long double* GetRho() const;
  long double* GetRhou() const;
  long double* GetRhov() const;
  long double* GetRhow() const;
  void SetEnergy(long double*);
  void SetRho(long double*);
  void SetRhou(long double*);
  void SetRhov(long double*);
  void SetRhow(long double*);

 private:

  long double addiabatic;		// addiabatic constant
  int numberx;			// number of cells in x
  long double bottomx, upperx;	// limits in x
  long double *rho;			// density in cell
  long double *rhou, *rhov, *rhow;	// momentum in cell
  long double *Energy;		// energy in cell

};

#endif
