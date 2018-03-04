/* This is the header file associated to MarginalCreaetion.cpp in which the prototypes for the
 * functions contained in that file are declared.  Any other header files which must be linked to
 * for the functions here are also included.
 *
 */

//************************//
//     INCLUDE GUARDS     //
//************************//

#ifndef MARGINALCREATION_H_
#define MARGINALCREATION_H_

//************************//
//        INCLUDES        //
//************************//

#include "LP_ompi.h"																					// allows the libraries included, macros defined and external variables declared in LP_ompi.h to be used in the MarginalCreation functions
#include "advection_1.h"																				// allows the external variables and function prototypes declared in advection_1.h to be used in the MarginalCreation functions

//************************//
//   FUNCTION PROTOTYPES  //
//************************//

double f_marg(vector<double>& U_vals0, int i, int j1, double x, double v1);

double f_marg_2D(vector<double>& U_vals0, int i1, int j1, double x1, double v1);

void PrintMarginalLoc(FILE *margfile_x1v1, FILE *margfile_x1x2);

void PrintMarginal(vector<double>& U_vals, FILE *margfile_x1v1, FILE *margfile_x1x2);

#endif /* MARGINALCREATION_H_ */
