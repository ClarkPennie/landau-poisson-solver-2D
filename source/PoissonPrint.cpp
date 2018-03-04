/* This is the source file which contains the subroutines necessary for printing values of the
 * potential and the field, suitable for plotting.
 *
 * Functions included: f_marg, PrintMarginalLoc, PrintMarginal
 *
 *  Created on: Feb 5, 2018
 */

#include "PoissonPrint.h"												// PrintFunctions.h is where the prototypes for the functions contained in this file are declared and any variables defined here to be used throughout the other files are declared as external

void PrintFieldLoc(FILE *phifile, FILE *Ex1file, FILE *Ex2file)			// function to print the values of x & v1 which the marginal will be evaluated at in the first two rows of the file with tag margfile (subsequent rows will be the values of the marginal at given timesteps)
{
	double x_val;														// declare x_val (the x value at the left edge of a given cell, in either direction)

	for(int i1=0; i1<Nx+1; i1++)
	{
		x_val = i1*dx;													// set x_val to the value of x at the left edge of the i1-th space cell
		for(int i2=0; i2<Nx+1; i2++)
		{
			fprintf(phifile, "%11.8g  ", x_val);						// in the file tagged as phifile, print the x1 coordinate
			fprintf(Ex1file, "%11.8g  ", x_val);						// in the file tagged as Ex1file, print the x1 coordinate
			fprintf(Ex2file, "%11.8g  ", x_val);						// in the file tagged as Ex2file, print the x1 coordinate
		}
	}
	fprintf(phifile, "\n");												// print a new line in the file tagged as phifile
	fprintf(Ex1file, "\n");												// print a new line in the file tagged as Ex1file
	fprintf(Ex2file, "\n");												// print a new line in the file tagged as Ex2file
	for(int i1=0; i1<Nx+1; i1++)
	{
		x_val = i1*dx;													// set x_val to the value of x at the left edge of the i1-th space cell
		for(int i2=0; i2<Nx+1; i2++)
		{
			x_val = i2*dx;												// set x_val to the value of x at the left edge of the i2-th space cell
			fprintf(phifile, "%11.8g  ", x_val);						// in the file tagged as phifile, print the x2 coordinate
			fprintf(Ex1file, "%11.8g  ", x_val);						// in the file tagged as Ex1file, print the x2 coordinate
			fprintf(Ex2file, "%11.8g  ", x_val);						// in the file tagged as Ex2file, print the x2 coordinate
		}
	}
	fprintf(phifile, "\n");												// print a new line in the file tagged as phifile
	fprintf(Ex1file, "\n");												// print a new line in the file tagged as Ex1file
	fprintf(Ex2file, "\n");												// print a new line in the file tagged as Ex2file
}

void PrintFieldData(FILE *phifile, FILE *Ex1file, FILE *Ex2file, vector<double>& POTC, vector<double>& phix, vector<double>& phiy)										// function to print the values of the potential and the field in the x1 & x2 directions in the file tagged as phifile, Ex1file & Ex2file, respectively, at the given timestep
{
	double phi_val, Ex1_val, Ex2_val;																			// declare x_0 (the x value at the left edge of a given cell), x_val (the x value to be evaluated at), phi_val (the value of phi evaluated at x_val), rho_val (the value of rho evaluated at x_val) & ddx (the space between x values)
	int ii, ii_b, ii_l, ii_bl;

	for(int i1=0; i1<Nx+1; i1++)
	{
		for(int i2=0; i2<Nx+1; i2++)
		{
			ii = i1*Nx + i2;
			ii_b = i1*Nx + i2-1;
			ii_l = (i1-1)*Nx + i2;
			ii_bl = (i1-1)*Nx + i2-1;
			if(i1<Nx && i2<Nx)
			{
				phi_val = POTC[3*ii] - 0.5*POTC[3*ii+1] - 0.5*POTC[3*ii+2];
				Ex1_val = phix[3*ii] - 0.5*phix[3*ii+1] - 0.5*phix[3*ii+2];
				Ex2_val = phiy[3*ii] - 0.5*phiy[3*ii+1] - 0.5*phiy[3*ii+2];
			}
			else if(i2==Nx && i1<Nx)
			{
				phi_val = POTC[3*ii_b] - 0.5*POTC[3*ii_b+1] + 0.5*POTC[3*ii_b+2];
				Ex1_val = phix[3*ii_b] - 0.5*phix[3*ii_b+1] + 0.5*phix[3*ii_b+2];
				Ex2_val = phiy[3*ii_b] - 0.5*phiy[3*ii_b+1] + 0.5*phiy[3*ii_b+2];
			}
			else if(i1==Nx && i2<Nx)
			{
				phi_val = POTC[3*ii_l] + 0.5*POTC[3*ii_l+1] - 0.5*POTC[3*ii_l+2];
				Ex1_val = phix[3*ii_l] + 0.5*phix[3*ii_l+1] - 0.5*phix[3*ii_l+2];
				Ex2_val = phiy[3*ii_l] + 0.5*phiy[3*ii_l+1] - 0.5*phiy[3*ii_l+2];
			}
			else
			{
				phi_val = POTC[3*ii_bl] + 0.5*POTC[3*ii_bl+1] + 0.5*POTC[3*ii_bl+2];
				Ex1_val = phix[3*ii_bl] + 0.5*phix[3*ii_bl+1] + 0.5*phix[3*ii_bl+2];
				Ex2_val = phiy[3*ii_bl] + 0.5*phiy[3*ii_bl+1] + 0.5*phiy[3*ii_bl+2];
			}

			fprintf(phifile, "%11.8g  ", phi_val);																// in the file tagged as phifile, print the value of the potential phi
			fprintf(Ex1file, "%11.8g  ", Ex1_val);																// in the file tagged as phifile, print the value of the potential phi
			fprintf(Ex2file, "%11.8g  ", Ex2_val);																// in the file tagged as phifile, print the value of the potential phi
		}
	}
	fprintf(phifile, "\n");																						// in the file tagged as phifile, print a new line
	fprintf(Ex1file, "\n");																						// in the file tagged as Ex1file, print a new line
	fprintf(Ex2file, "\n");																						// in the file tagged as Ex2file, print a new line
}
