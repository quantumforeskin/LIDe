#ifndef __FLUID_
#define __FLUID_

#include<iostream>
#include<fstream>
#include<iomanip>
#include<cmath>
using namespace std;

typedef long double ld;

class Fluid {

 public:

  /*  CONSTRUCTOR AND DESTRUCTOR  */
  Fluid(ld,int,ld,ld);
  ~Fluid();

  /*  TIME EVOLUTION  */
  void Evolution1(const char*,ld,int);
  void Evolution2(const char*,ld,int);
  void Evolution3(const char*,ld,int);
  void Evolution1c(const char*,ld,int);

  /*  PHYSICAL VARIABLES  */
  ld* GetF(int) const;
  ld* GetG(int) const;
  ld* GetH(int) const;
  ld Pressure(int) const;
  ld* RoeAverage(int) const;
  ld* Eigenvalues(ld*) const;
  ld** Eigenvectors(ld*) const;
  ld* GetAlpha(int,ld*) const;
  ld** GetWaves(ld*,ld**) const;
  ld* GetFb(int) const;
  ld* GetFt(int) const;

  /*  INTERNAL VARIABLE FUNCTIONS  */
  ld* GetEnergy() const;
  ld* GetRho() const;
  ld* GetRhou() const;
  ld* GetRhov() const;
  ld* GetRhow() const;
  void SetEnergy(ld*);
  void SetRho(ld*);
  void SetRhou(ld*);
  void SetRhov(ld*);
  void SetRhow(ld*);

 private:

  ld addiabatic;		// addiabatic constant
  int numberx;			// number of cells in x
  ld bottomx, upperx;	// limits in x
  ld *rho;			// density in cell
  ld *rhou, *rhov, *rhow;	// momentum in cell
  ld *Energy;		// energy in cell

};

#endif
