/* This is the header file associated to MomentCalculations.cpp in which the prototypes for the
 * functions contained in that file are declared.  Any other header files which must be linked to
 * for the functions here are also included.
 *
 */

//************************//
//     INCLUDE GUARDS     //
//************************//

#ifndef MOMENTCALCULATIONS_H_
#define MOMENTCALCULATIONS_H_

//************************//
//        INCLUDES        //
//************************//

#include "LP_ompi.h"																					// allows the libraries included, macros defined and external variables declared in LP_ompi.h to be used in the MomentCalculations functions
#include "advection_1.h"																				// allows the external variables and function prototypes declared in advection_1.h to be used in the MomentCalculations functions

//************************//
//   FUNCTION PROTOTYPES  //
//************************//

double computeMass(vector<double>& U_vals);

void computeMomentum(vector<double>& U_vals, double *a);

double computeKiE(vector<double>& U_vals);

//double computeKiEratio(vector<double>& U_vals, int *NegVals);

double computeEleE(vector<double>& phix_vals, vector<double>& phiy_vals);

#endif /* MOMENTCALCULATIONS_H_ */
