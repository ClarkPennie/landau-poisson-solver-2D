/* This is the header file associated to the main file LP_ompi.cpp in which all libraries required
 * by the program are included, all macros (to decide the behaviour of a given run) are defined and
 * all variables to be used throughout the various files are declared as external.
 *
 */

//************************//
//     INCLUDE GUARDS     //
//************************//

#ifndef LP_OMPI_H_
#define LP_OMPI_H_

//************************//
//        LIBRARIES       //
//************************//

#include <mpi.h>																					// allows all MPI routines to be used
#include <stdio.h>																					// allows the object type FILE (an object suitable for storing information for a file stream) to be used, as well as the functions printf, sprintf, fopen, fread, fclose
#include <malloc.h>																					// allows malloc to be used
#include <math.h>																					// allows sqrt to be used as well as the value of M_PI
#include <stdlib.h>																					// allows malloc & free to be used
#include <omp.h>																					// allows all OpenMP routines to be used
#include <fftw3.h>																					// allows the Fast Fourier Transform to be used
#ifdef HAVE_MKL
#include <mkl.h>
#elif HAVE_OPENBLAS
#include <lapacke.h>
#endif
#include <iostream>
#include <vector>
using std::vector;
using std::cout;
using std::endl;


//************************//
//         MACROS         //
//************************//

// CHOOSE WHETHER OR NOT TO USE MPI:
#define UseMPI 																						// define the macro MPI (UNCOMMENT IF THE CODE SHOULD UTILISE MPI)

// CHOOSE WHICH VARIATION OF THE CODE TO RUN:
#define Damping																						// define the macro Damping (UNCOMMENT IF BEING RUN FOR THE LANDAU DAMPING PROBLEM)
//#define FullandLinear 																			// define the macro FullandLinear (UNCOMMENT IF THE LINEAR ELE-ION COLLISION OPERATOR IS NEEDED; OTHERWISE, ONLY ELE-ELE COLLISIONS)
//#define TwoStream																					// define the macro TwoStream (UNCOMMENT IF BEING RUN FOR THE TWO STREAM PROBLEM)
//#define FourHump																					// define the macro FourHump (UNCOMMENT IF BEING RUN FOR THE FOUR HUMP IC PROBLEM)
//#define TwoHump																					// define the macro TwoHump (UNCOMMENT IF BEING RUN FOR THE TWO HUMP IC PROBLEM)

// CHOOSE IF THIS IS THE INITIAL RUN OR A SUBSEQUENT RUN:
#define First																						// define the macro First (UNCOMMENT IF RUNNING THE CODE FOR THE FIRST TIME)
//#define Second																					// define the macro Second (UNCOMMENT IF PICKING UP DATA FROM A PREVIOUS RUN)

/*POISSON SOLVER MACROS*/

#define MPD 7																						// define the macro MPD (the max order of LDG scheme; e.g. p0=1,p1=1+2=3, p2=6, etc.) and set its value
#define MGD 7																						// define the macro MGD (the max gauss point used i) and set its value
#define NXD 51																						// define the macro NXD (the dimension of arrays in the x1 direction, where NXD > NX+1) and set its value (Jose's value: 16)
#define NYD 51																						// define the macro NYD (the dimension of some arrays in the x2 direction, where NYD > NY+1) and set its value (Jose's value: 16)
#define NYD2 105																					// define the macro NYD2 (the dimension of some arrays in the x2 direction) and set its value (Jose's value: 20)

// dimensions for potential and related matrices
// We must define N2XD > 2*(NX + 2) and N2YD > 2*(1 + 7*NY/6)
#define N2XD 120																					// define the macro N2XD (the dimension of refined arrays in the x1 direction, where N2XD > 2(NX+2) and set its value (Jose's value: 44)
#define N2YD 120																					// define the macro N2YD (the dimension of refined arrays in the x2 direction, where N2YD > 2(NY+2) and set its value (Jose's value: 44


//************************//
//   EXTERNAL VARIABLES   //
//************************//

extern double PI;																					// declare PI and set it to M_PI (the value stored in the library math.h)
extern int M;																						// declare M (the number of collision invarients) and set it equal to 5
extern int Nx, Nv, nT, N; 											 								// declare Nx (no. of x discretised points), Nv (no. of v discretised point), nT (no. of time discretised points) & N (no. of nodes in the spectral method) and setting all their values
extern int size_x, size_v, size, size_ft; 															// declare size_v (no. of total v discretised points in 3D) and set it to Nv^3, size (the total no. of discretised points) and set it to size_v*Nx & size_ft (total no. of spectral discretised points in 3D) and set it to N*N*N

extern int NX, NY, NYREAL;																			// declare NX (no. of x1 discretised points for the Poisson solver), NY (no. of x2 discretised points for the Poisson solver) & NYREAL (no. of x2 discretised points for the Poisson solver if there is an oxide-silicon region on top)

extern double A_amp, k_wave;																		// declare A_amp (the amplitude of the perturbing wave) & k_wave (the wave number of the perturbing wave) and set their values
extern double Lx, Lv;																				// declare Lx (for 0 < x < Lx) and set it to & Lv (for -Lv < v < Lv in the advection problem) and set their values
extern double dv, dx; 																				// declare dv (the velocity stepsize) and set it to 2Lv/Nv & dx (the space stepsize) and set it to Lx/Nx
extern double L_v, R_v, L_eta;																		// declare L_v (for -Lv < v < Lv in the collision problem) and set it to Lv, R_v (for v in B_(R_v) in the collision problem) and set it to Lv & L_eta (for Fourier space, -L_eta < eta < L_eta)
extern double h_eta, h_v;																			// declare h_eta (the Fourier stepsize) & h_v (also the velocity stepsize but for the collision problem)
extern double nu, dt, nthread; 																		// declare nu (1/knudson#) and set it to 0.1, dt (the timestep) and set it to 0.004 & nthread (the number of OpenMP threads) and set it to 16

extern double *v, *eta;																				// declare v (the velocity variable) & eta (the Fourier space variable)
extern double *wtN;																					// declare wtN (the trapezoidal rule weights to be used)
extern double scale, scale3, scaleL, scalex, scalev;												// declare scale (the 1/sqrt(2pi) factor appearing in Gaussians), scale (the 1/(sqrt(2pi))^3 factor appearing in the Maxwellian), scaleL (the volume of the velocity domain) and set it to 8Lv^3 & scalev (the volume of a discretised velocity element) and set it to dv^3
extern double **C1, **C2;																			// declare pointers to matrices C1 (the real part of the conservation matrix C) & C2 (the imaginary part of the conservation matrix C), CCt (of dimension 5x5) & CCt_linear (of dimension 2x2)
extern double CCt[5*5], CCt_linear[2*2];															// declare matrices CCt (C*C^T, for the conservation matrix C) & CCt_linear (C*C^T, for the conservation matrix C, in the two species collision operator)
extern double lamb[5], lamb_linear[2];																// declare the arrays lamb (to hold 5 values) & lamb_linear (to hold 2 values)

extern vector<double> U1;
extern double *Utmp, *output_buffer_vp;//, **H;														// declare pointers to Utmp (both used to help store the values in U, declared later) & output_buffer_vp (a buffer used when sending the data between MPI processes during the VP method)
extern double *Q, *f1, *Q1, *Utmp_coll;//*f2, *f3;													// declare pointers to Q (the discretised collision operator), f1 (used to help store the solution during the collisional problem), Q1 (used in calculation of the collision operator) & Utmp_coll (used to store calculations from the RK4 method used in the collisional problem)
extern fftw_complex *Q1_fft, *Q2_fft, *Q3_fft;														// declare pointers to the complex numbers Q1_fft, Q2_fft & Q3_fft (involved in storing the FFT of Q)

#ifdef FullandLinear																				// only do this if FullandLinear was defined
extern fftw_complex *Q1_fft_linear, *Q2_fft_linear, *Q3_fft_linear;									// declare pointers to the complex numbers Q1_fft_linear, Q2_fft_linear & Q3_fft_linear (involved in storing the FFT of the two species collison operator Q)
#endif

extern fftw_complex *fftIn, *fftOut;																// declare pointers to the FFT variables fftIn (a vector to be to have the FFT applied to it) & fftOut (the output of an FFT)
//extern double IntM[10];																				// declare an array IntM to hold 10 double variables
//#pragma omp threadprivate(IntM)																	// start the OpenMP parallel construct to start the threads which will run in parallel, passing IntM to each thread as private variables which will have their contents deleted when the threads finish (doesn't seem to be doing anything since no {} afterwards???)

extern double ce, *cp, *intE, *intE1, *intE2;														// declare ce and pointers to cp, intE, intE1 & intE2 (precomputed quantities for advections)

extern vector<double> IE_X, IXE_X, IYE_X, IXXE_X, IXYE_X, IYYE_X;									// declare vectors to store the integrals of the field in the x direction in each space cell (see PoissonVariables.cpp for an explanation of each vector)
extern vector<double> IE_Y, IXE_Y, IYE_Y, IXXE_Y, IXYE_Y, IYYE_Y;									// declare vectors to store the integrals of the field in the y direction in each space cell(see PoissonVariables.cpp for an explanation of each vector)
extern vector<double> SE_X, SE_Y;																	// declare vectors to store the sign of the components of the electric field in the each space cell

extern fftw_plan p_forward; 																		// declare the fftw_plan p_forward (an object which contains all the data which allows fftw3 to compute the FFT)
extern fftw_plan p_backward; 																		// declare the fftw_plan p_backward (an object which contains all the data which allows fftw3 to compute the inverse FFT)
extern fftw_complex *temp;																			// declare a pointer to complex number temp (a temporary array used when calculating the FFT)

extern int myrank_mpi, nprocs_mpi, nprocs_Nx;														// declare myrank_mpi (the rank of the current MPI process running), nprocs_mpi (the total number of MPI processes) & nprocs_Nx (the amount of MPI processes used for the collisionless VP problem)
extern int chunksize_dg, chunksize_ft, chunk_Nx;													// declare chunksize_dg (the amount of data each processes works on during the DG method in the VP problem), chunksize_ft (the amount of data each process works on during the collisional problem) & chunk_Nx (the number of space-steps sent to each process in the collisional problem)

extern int *fNegVals;																				// declare fNegVals (to store where DG solution goes negative - a 1 if negative and a 0 if positive)
extern double *fAvgVals;																			// declare fAvgVals (to store the average values of f on each cell)
extern double *fEquiVals;																			// declare f_equivals (to store the equilibrium solution)

//extern double a[3];

//************************//
//        INCLUDES        //
//************************//

#include "advection_1.h"																			// allows computeMass, computeMomentum, computeKiE, computeEleE & RK3 to be used
#include "SetInit_1.h"																				// allows TrapezoidalRule, SetInit_LD, SetInit_4H, SetInit_2H & setInit_spectral to be used
#include "conservationRoutines.h"         															// allows createCCtAndPivot & conserveAllMoments to be used
#include "collisionRoutines_1.h"            														// allows generate_conv_weights, generate_conv_weights_linear, computeQ & RK4 to be used
#include "MomentCalculations.h"																		// allows computeMass, computeMomentum, computeKiE, computeKiERatio, computeEleE to be used
#include "EntropyCalculations.h"																	// allows computeEntropy, computeEntropy_wAvg & computeRelEntropy to be used
#include "MarginalCreation.h"																		// allows PrintMarginalLoc & PrintMarginal to be used
//#include "EquilibriumSolution.h"																	// allows ExportRhoQuadVals, ComputeEquiVals & PrintEquiVals to be used
#include "NegativityChecks.h"																		// allows computeCellAvg, FindNegVals & CheckNegVals to be used
#include "PoissonVariables.h"																		// where the global variables to be used throughout the various files are defined as external
#include "PoissonFunctions.h"																		// allows adim, config, setup, setup_matrix, setup_pois & pois2d to be used
#include "PoissonPrint.h"																			// allows PrintFieldLoc & PrintFieldData to be used


#endif /* LP_OMPI_H_ */
