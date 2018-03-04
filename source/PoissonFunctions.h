/* This is the header file associated to PoissonFunctions.cpp in which the prototypes for the
 * functions contained in that file are declared.  Any other header files which must be linked to
 * for the functions here are also included.
 *
 *  Created on: Feb 14, 2018
 */

//************************//
//     INCLUDE GUARDS     //
//************************//

#ifndef POISSONFUNCTIONS_H_
#define POISSONFUNCTIONS_H_

//************************//
//        INCLUDES        //
//************************//

#include "LP_ompi.h"																					// allows the libraries included, macros defined and external variables declared in LP_ompi.h to be used in the advection_1 functions

//************************//
//   FUNCTION PROTOTYPES  //
//************************//

void adim(void);

void coef_pois(void);

void config(void);

double fle(int k, double x);

double fled(int k, double x);

void pois2d(vector<double>& Ut, vector<double>& POTC, vector<double>& phix, vector<double>& phiy);

void setup(void);

void setup_matrix(void);

void setup_pois(void);

void transmat(void);

#endif /* POISSONFUNCTIONS_H_ */
