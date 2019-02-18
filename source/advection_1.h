/* This is the header file associated to advection_1.cpp in which the prototypes for the functions
 * contained in that file are declared and any variables defined there to be used throughout the
 * other files are declared as external.  Any other header files which must be linked to for the
 * functions here are also included.
 *
 */

//************************//
//     INCLUDE GUARDS     //
//************************//

#ifndef ADVECTION_1_H_
#define ADVECTION_1_H_

//************************//
//        INCLUDES        //
//************************//

#include "LP_ompi.h"																					// allows the libraries included, macros defined and external variables declared in LP_ompi.h to be used in the advection_1 functions

//************************//
//   EXTERNAL VARIABLES   //
//************************//

extern double wt[5];																					// weights for Gaussian quadrature
extern double vt[5];																					// node values for Gaussian quadrature over the interval [-1,1]

//************************//
//   FUNCTION PROTOTYPES  //
//************************//

double Gridv(double m);

double Gridx(double m);

//double rho_x(double x, vector<double>& U, int i);

//double rho(vector<double>& U, int i);

//double computePhi_x_0(vector<double>& U);

//double computePhi(double *U, double x, int ix);	STILL NEEDS TO BE ADDED FOR 2D

//void PrintPhiVals(double *U, FILE *phifile);		STILL NEEDS TO BE ADDED FOR 2D

//double computeC_rho(vector<double>& U, int i);

//double Int_Int_rho(vector<double>& U, int i);

//double Int_Int_rho1st(vector<double>& U, int i);

//double Int_E(vector<double>& U, int i);

//double Int_E1st(vector<double>& U, int i);

//double Int_fE(vector<double>& U, int i_mod, int j);

double Int_fE1(vector<double>& U, int i, int j);

double Int_fE2(vector<double>& U, int i, int j);

//double Int_E2nd(vector<double>& U, int i);

double I1(vector<double>& U_vals0, int k, int l);

double I2(vector<double>& U_vals0, int k, int l);

double I3_x1(vector<double>& U_vals0, int k, int l);

double I3_x2(vector<double>& U_vals0, int k, int l);

double I5_v1(vector<double>& U_vals0, int k, int l);

double I5_v2(vector<double>& U_vals0, int k, int l);

//void computeH(double *U);

void RK3(vector<double>& U_vals, vector<double>& POTC, vector<double>& phix, vector<double>& phiy);

void RK3_NoField(vector<double>& U_vals);

#endif /* ADVECTION_1_H_ */
