/* This is the header file associated to PrintFunctions.cpp in which the prototypes for the
 * functions contained in that file are declared.  Any other header files which must be linked to
 * for the functions here are also included.
 *
 *  Created on: Feb 14, 2018
 */

//************************//
//     INCLUDE GUARDS     //
//************************//

#ifndef POISSONPRINT_H_
#define POISSONPRINT_H_

//************************//
//        INCLUDES        //
//************************//

#include "LP_ompi.h"																					// allows the libraries included, macros defined and external variables declared in TestSolver.h to be used in the PrintFunctions functions

//************************//
//   FUNCTION PROTOTYPES  //
//************************//

void PrintFieldData(FILE *phifile, FILE *Ex1file, FILE *Ex2file, vector<double>& POTC, vector<double>& phix, vector<double>& phiy);				// function to print the values of the potential and the field in the x1 & x2 directions in the file tagged as phifile, Ex1file & Ex2file, respectively, at the given timestep

void PrintFieldLoc(FILE *phifile, FILE *Ex1file, FILE *Ex2file);																				// function to print the values of x & v1 which the marginal will be evaluated at in the first two rows of the file with tag margfile (subsequent rows will be the values of the marginal at given timesteps)

#endif /* PRINTFUNCTIONS_H_ */
