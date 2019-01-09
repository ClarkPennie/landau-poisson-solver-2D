/* This is the header file associated to PoissonVariables.cpp in which all variables to be used
 * throughout the various files are declared as external.  Any other header files which must be
 * linked to for the functions here are also included.
 *
 *  Created on: Feb 14, 2018
 */

//************************//
//     INCLUDE GUARDS     //
//************************//

#ifndef POISSONVARIABLES_H_
#define POISSONVARIABLES_H_

//************************//
//        INCLUDES        //
//************************//

#include "LP_ompi.h"																					// allows the libraries included, macros defined and external variables declared in LP_ompi.h to be used in the advection_1 functions

//************************//
//   EXTERNAL VARIABLES   //
//************************//

extern double tps;

// constants (defined in adim)
extern double PD, PEP, PVE, CDS, CDOS;

// parameters (defined in config)
extern double DX[NXD], DY[NYD], CX[NXD], CY[NYD], DOP[NXD][NYD];

// The potential (in the refined grid!)
extern double POT[N2XD][N2YD];

/*
  where, set R = [x_{i-1/2}, x_{i+1/2}] * [y_{j-1/2}, y_{j+1/2}]
                           (the cell in (x,y) plane),

                   1         ( (
  Iaaa[i][j] = ----------    | |  bbb(x,y) dx dy
               measure(R)    ) )
                               R
and
     aaa         =>        bbb(x,y)
     ===================================

     E_X               E_x(x,y)

     XE_X              E_x(x,y) 2(x - x_i)/dx_i

     YE_X              E_x(x,y) 2(y - y_j)/dy_j

     XX_EX             E_x(x,y) 2(x - x_i)/dx_i 2(x - x_i)/dx_i

     XYE_X             E_x(x,y) 2(x - x_i)/dx_i 2(y - y_j)/dy_j

     YYE_X             E_x(x,y) 2(y - y_j)/dy_j 2(y - y_j)/dy_j

*/

// parameters (defined in setup)
//int NX,NY,mppois,mpto,ldab;
extern int mppois,mpto,ldab;
extern double xa,xb,xc,xd,ya,yb,pleft,pbot,pright,ptop;
extern double pybc; // DEBUG Fix y boundary conditions

//sub and super band for the big linear system, defined in setup
extern int nsub,nsuper;

//defined in setup
extern int k01;

extern int *ipiv,*info,ndim;

//compuational constants (defined in setup)
extern double c11;

//gauss points, defined in setup
extern double gauss[10][3], gaut[10][10];

//coef defined in setup
extern double epsil[NXD][NYD2];

//right hand side of the linear system, defined in coef
extern double ff[MPD][NXD][NYD2];

//right hand side of the big linear system, defined in coef
extern double poisb[MPD*NXD*NYD2];

//left hand side of the big linear system, defined in coef
//size of ab is (2*nsub+nsuper+1), (mp+1)nx
//be very careful with this, I overdefine the size, be careful when using c
//transfer it to 1d matrix, cab
//double ab[7*MPD+5][MPD*NXD],cab[(7*MPD+5)*MPD*NXD];
extern double poisab[4*MPD*NYD2][MPD*NXD*NYD2],*poiscab;

//compuational array for the basis (defined in setup)
extern double ai[MPD], aii[MPD][MPD];

//array to multiply in each element, defined in coef
extern double coc[MPD][MPD],col[MPD][MPD],cor[MPD][MPD],cott[MPD][MPD],cob[MPD][MPD];

//defined in coef, related to the legendre polynomial
extern double bl[MPD][MGD], br[MPD][MGD], bb[MPD][MGD], bt[MPD][MGD];

extern double  cx[MPD][MGD][MGD],cy[MPD][MGD][MGD],cxy[MPD][MGD][MGD];

extern double cvx[MPD][MPD],cvy[MPD][MPD],cyll[MPD][MPD],cylr[MPD][MPD],cyrl[MPD][MPD],cyrr[MPD][MPD],
       cxll[MPD][MPD],cxlr[MPD][MPD],cxrl[MPD][MPD],cxrr[MPD][MPD];

extern double acvx[MPD][MPD],acvy[MPD][MPD],acyll[MPD][MPD],acylr[MPD][MPD],acyrl[MPD][MPD],acyrr[MPD][MPD],
       acxll[MPD][MPD],acxlr[MPD][MPD],acxrl[MPD][MPD],acxrr[MPD][MPD];

#endif /* POISSONVARIABLES_H_ */
