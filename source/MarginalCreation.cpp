/* This is the source file which contains the subroutines necessary for producing marginals of the
 * solutions, suitable for plotting.
 *
 * Functions included: f_marg, PrintMarginalLoc, PrintMarginal
 *
 *  Created on: Nov 15, 2017
 */

#include "MarginalCreation.h"																	// MarginalCreation.h is where the prototypes for the functions contained in this file are declared

double f_marg(vector<double>& U_vals0, int i, int j1, double x, double v1)
{
	int j2, j3, k0, j2N, k;																		// declare j2, j3 (the indices of the velocity cell in the v2 & v3 directions), k0 (to store the value of i*Nv^3 + j1_Nv^2), j2N (to store the value of j2*N) & k (the index of the cell in U)
	double x_dif, v1_dif, retn;																	// declare x_dif (to store x - x_i), v1_dif (to store v1 - v_j1) & retn (the value of the marginal evaluated at the given x & v1 to be returned at the end

	x_dif = x -  Gridx((double)i);																// set x_dif to x - x_i
	v1_dif = v1 - Gridv((double)j1);															// set v1_dif to v1 - v_j1
	k0 = i*Nv*Nv*Nv + j1*Nv*Nv;																	// set k0 to i*Nv^3 + j1*Nv^2
	retn = 0;																					// initialise fM_val at 0, since the value is calculated through a sum

	//#pragma omp parallel for schedule(dynamic) private(j2,j2N,j3,k)  shared(U,Nv,k0,dx,dv,x_dif,v1_dif) reduction(+:retn)
	for(j2=0; j2<Nv; j2++)
	{
		j2N = j2*Nv;																			// set j2N to j2*Nv
		for(j3=0; j3<Nv; j3++)
		{
			k = k0 + j2N + j3;																	// set k to i*Nv^3 + j1*Nv^2 + j2*Nv + j3
			retn += dv*dv*U_vals0[7*k+0] + dv*dv*U_vals0[7*k+1]*x_dif/dx + dv*U_vals0[7*k+3]*v1_dif
							+ U_vals0[7*k+6]*(v1_dif*v1_dif +dv*dv/6);								// add dv*dv*U[6*k+0] + dv*dv*U[6*k+1]*x_dif/dx + dv*U[6*k+2]*v1_dif + U[6*k+5]*(v1_dif*v1_dif +dv*dv/6) for the given j2 & j3 in the sum for retn
		}
	}
	return retn;																				// return the value of the marginal evaluated at x & v1
}

double f_marg2D_x1v1(vector<double>& U_vals0, int i1, int j1, double x1, double v1)
{
	int i2, j2, j3, k0, i2NNN, j2N, k;															// declare i2 (the index of the space cell in the x2 direction), j2, j3 (the indices of the velocity cell in the v2 & v3 directions), k0 (to store the value of i*Nv^3 + j1_Nv^2), i2NNN (to store the value of i2*N^3), j2N (to store the value of j2*N) & k (the index of the cell in U)
	double x1_dif, v1_dif, retn;																// declare x1_dif (to store x1 - x_i1), v1_dif (to store v1 - v_j1) & retn (the value of the marginal evaluated at the given x & v1 to be returned at the end

	x1_dif = x1 -  Gridx((double)i1);															// set x1_dif to x1 - x_i1
	v1_dif = v1 - Gridv((double)j1);															// set v1_dif to v1 - v_j1
	k0 = i1*Nx*size_v + j1*Nv*Nv;																// set k0 to i1*NxNv^3 + j1*Nv^2
	retn = 0;																					// initialise fM_val at 0, since the value is calculated through a sum

	//#pragma omp parallel for schedule(dynamic) private(j2,j2N,j3,k)  shared(U,Nv,k0,dx,dv,x_dif,v1_dif) reduction(+:retn)
	for(i2=0; i2<Nx; i2++)
	{
		i2NNN = i2*size_v;																		// set i2NNN to i2*Nv^3
		for(j2=0; j2<Nv; j2++)
		{
			j2N = j2*Nv;																		// set j2N to j2*Nv
			for(j3=0; j3<Nv; j3++)
			{
				k = k0 + i2NNN + j2N + j3;																// set k to i*Nv^3 + j1*Nv^2 + j2*Nv + j3
				retn += dx*dv*dv*U_vals0[7*k+0] + dv*dv*U_vals0[7*k+1]*x1_dif + dx*dv*U_vals0[7*k+3]*v1_dif
						+ dx*U_vals0[7*k+6]*(v1_dif*v1_dif +dv*dv/6.);								// add dx*dv*dv*U_vals0[7*k+0] + dv*dv*U_vals0[7*k+1]*x1_dif + dx*dv*U_vals0[7*k+3]*v1_dif + dx*U_vals0[7*k+6]*(v1_dif*v1_dif +dv*dv/6.) for the given j2 & j3 in the sum for retn
			}
		}
	}
	return retn;																				// return the value of the marginal evaluated at x & v1
}

double f_marg2D_x1x2(vector<double>& U_vals0, int i1, int i2, double x1, double x2)
{
	int j1, j2, j3, k0, j1NN, j2N, k;															// declare j1, j2, j3 (the indices of the velocity cell in the v1, v2 & v3 directions), k0 (to store the value of i1*NxNv^3 + i2*Nv^3), j1NN (to store the value of j1*Nv^2), j2N (to store the value of j2*N) & k (the index of the cell in U)
	double x1_dif, x2_dif, retn;																// declare x1_dif (to store x1 - x_i1), x2_dif (to store x2 - x_i2) & retn (the value of the marginal evaluated at the given x1 & x2 to be returned at the end

	x1_dif = x1 -  Gridx((double)i1);															// set x1_dif to x1 - x_i1
	x2_dif = x2 - Gridv((double)i2);															// set x2_dif to x2 - x_i2
	k0 = (i1*Nx + i2)*size_v;																	// set k0 to i1*NxNv^3 + i2*Nv^3
	retn = 0;																					// initialise fM_val at 0, since the value is calculated through a sum

	//#pragma omp parallel for schedule(dynamic) private(j2,j2N,j3,k)  shared(U,Nv,k0,dx,dv,x_dif,v1_dif) reduction(+:retn)
	for(j1=0; j1<Nx; j1++)
	{
		j1NN = j1*Nv*Nv;																		// set j1NN to j1*Nv^2
		for(j2=0; j2<Nv; j2++)
		{
			j2N = j2*Nv;																		// set j2N to j2*Nv
			for(j3=0; j3<Nv; j3++)
			{
				k = k0 + j1NN + j2N + j3;														// set k to i1*NxNv^3 + i2*Nv^3 + j1*Nv^2 + j2*Nv + j3
				retn += U_vals0[7*k+0] + U_vals0[7*k+1]*x1_dif/dx + U_vals0[7*k+2]*x2_dif/dx
						+ U_vals0[7*k+6]/4.;													// add U_vals0[7*k+0] + U_vals0[7*k+1]*x1_dif/dx + U_vals0[7*k+2]*x2_dif/dx + U_vals0[7*k+6]/4. for the given j1, j2 & j3 in the sum for retn
			}
		}
	}
	retn = retn*scalev;																			// multiply the sum bydv^3
	return retn;																				// return the value of the marginal evaluated at x & v1
}

double f_marg_2D(vector<double>& U_vals0, int i1, int j1, double x1, double v1)
{
	int i2, j2, j3, k0, i2NNN, j2N, k;															// declare i2 (the index of the space cell in the x1 direction) j2, j3 (the indices of the velocity cell in the v2 & v3 directions), k0 (to store the value of i*Nv^3 + j1_Nv^2), i2NNN (to store the value of i2*Nv^3) j2N (to store the value of j2*Nv) & k (the index of the cell in U)
	double x1_dif, v1_dif, retn;																	// declare x_dif (to store x - x_i), v1_dif (to store v1 - v_j1) & retn (the value of the marginal evaluated at the given x & v1 to be returned at the end

	x1_dif = x1 - Gridx((double)i1);															// set x1_dif to x1 - x1_i
	v1_dif = v1 - Gridv((double)j1);															// set v1_dif to v1 - v_j1
	k0 = i1*Nx*size_v + j1*Nv*Nv;																// set k0 to i1*Nx*Nv^3 + j1*Nv^2
	retn = 0;																					// initialise fM_val at 0, since the value is calculated through a sum

	//#pragma omp parallel for schedule(dynamic) private(j2,j2N,j3,k)  shared(U,Nv,k0,dx,dv,x_dif,v1_dif) reduction(+:retn)
	for(i2=0; i2<Nx; i2++)
	{
		i2NNN = i2*size_v;																		// set i2NNN to i2*Nv^3
		for(j2=0; j2<Nv; j2++)
		{
			j2N = j2*Nv;																		// set j2N to j2*Nv
			for(j3=0; j3<Nv; j3++)
			{
				k = k0 + i2NNN + j2N + j3;														// set k to i1*Nx*Nv^3 + i2*Nv^3 + j1*Nv^2 + j2*Nv + j3
				retn += dx*dv*dv*U_vals0[7*k+0] + dv*dv*U_vals0[7*k+1]*x1_dif + dx*dv*U_vals0[7*k+3]*v1_dif
						+ dx*U_vals0[7*k+6]*(v1_dif*v1_dif +dv*dv/6);								// add dv*dv*U[6*k+0] + dv*dv*U[6*k+1]*x_dif/dx + dv*U[6*k+2]*v1_dif + U[6*k+5]*(v1_dif*v1_dif +dv*dv/6) for the given j2 & j3 in the sum for retn
			}
		}
		return retn;																				// return the value of the marginal evaluated at x & v1
	}
}

void PrintMarginalLoc(FILE *margfile_x1v1, FILE *margfile_x1x2)									// function to print the values of x & v1 which the marginal will be evaluated at in the first two rows of the file with tag margfile (subsequent rows will be the values of the marginal at given timesteps)
{
	int i1, i2, j1, np, nx, nv;																	// declare i1, i2 (the indices of the space cells in the x1 & x2 directions),  j1 (the index of the velocity cell in the v1 direction), np (the number of points to evaluate in a given space/velocity cell), nx (a counter for the points in the space cell) & nv (a counter for the points in the velocity cell)
	double x1_0, x1_val, x2_0, x2_val, v1_0, v1_val, ddx, ddv;									// declare x1_0 (the x1 value at the left edge of a given cell), x1_val (the x1 value to be evaluated at), x2_0 (the x2 value at the left edge of a given cell), x2_val (the x2 value to be evaluated at), declare v1_0 (the v1 value at the left edge of a given cell), v1_val (the v1 value to be evaluated at), ddx (the space between x values) & ddv (the space between v1 values)

	np = 4;																						// set np to 4
	ddx = dx/np;																				// set ddx to the space cell width divided by np
	ddv = dv/np;																				// set ddv to the velocity cell width divided by np
	for(i1=0; i1<Nx; i1++)
	{
		x1_0 = Gridx((double)i1 - 0.5);															// set x1_0 to the value of x1 at the left edge of the i1-th space cell
		for (nx=0; nx<np; nx++)
		{
			x1_val = x1_0 + nx*ddx;																// set x1_val to x1_0 plus nx increments of width ddx
			for(j1=0; j1<Nv; j1++)
			{
				for (nv=0; nv<np; nv++)
				{
					fprintf(margfile_x1v1, "%11.8g  ", x1_val);									// in the file tagged as margfile_x1v1, print the x1 coordinate
				}
			}
			for(i2=0; i2<Nx; i2++)
			{
				for (nv=0; nv<np; nv++)
				{
					fprintf(margfile_x1x2, "%11.8g  ", x1_val);									// in the file tagged as margfile_x1x2, print the x1 coordinate
				}
			}
		}
	}
	fprintf(margfile_x1v1, "\n");																// print a new line in the file tagged as margfile_x1v1
	fprintf(margfile_x1x2, "\n");																// print a new line in the file tagged as margfile_x1x2
	for(i1=0; i1<Nx; i1++)
	{
		for (nx=0; nx<np; nx++)
		{
			for(j1=0; j1<Nv; j1++)
			{
				v1_0 = Gridv((double)j1 - 0.5);													// set v1_0 to the value of v1 at the left edge of the j1-th velocity cell in the v1 direction
				for (nv=0; nv<np; nv++)
				{
					v1_val = v1_0 + nv*ddv;														// set v1_val to v1_0 plus nv increments of width ddv
					fprintf(margfile_x1v1, "%11.8g  ", v1_val);									// in the file tagged as margfile_x1v1, print the v1 coordinate
				}
			}
			for(i2=0; i2<Nx; i2++)
			{
				x2_0 = Gridx((double)i2 - 0.5);													// set x2_0 to the value of x2 at the left edge of the i2-th velocity cell in the i2 direction
				for (nv=0; nv<np; nv++)
				{
					x2_val = x2_0 + nv*ddx;														// set x2_val to x2_0 plus nv increments of width ddv
					fprintf(margfile_x1x2, "%11.8g  ", x2_val);									// in the file tagged as margfile_x1x2, print the x2 coordinate
				}
			}
		}
	}
	fprintf(margfile_x1v1, "\n");																// print a new line in the file tagged as margfile_x1v1
	fprintf(margfile_x1x2, "\n");																// print a new line in the file tagged as margfile_x1x2
}

void PrintMarginal(vector<double>& U_vals, FILE *margfile_x1v1, FILE *margfile_x1x2)													// function to print the values of the marginal in the file tagged as margfile at the given timestep
{
	int i1, i2, j1, np, nx, nv;																	// declare i1, i2 (the indices of the space cell in the x1 & x2 direction),  j1 (the index of the velocity cell in the v1 direction), np (the number of points to evaluate in a given space/velocity cell), nx (a counter for the points in the space cell) & nv (a counter for the points in the velocity cell)
	double x1_0, x1_val, x2_0, x2_val, v1_0, v1_val, fM_val, ddx, ddv;							// declare x1_0 (the x1 value at the left edge of a given cell), x1_val (the x1 value to be evaluated at), x2_0 (the x2 value at the left edge of a given cell), x2_val (the x2 value to be evaluated at), declare v1_0 (the v1 value at the left edge of a given cell), v1_val (the v1 value to be evaluated at), fM_val (the value of the marginal evaluated at (x_val, v1_val), ddx (the space between x values) & ddv (the space between v1 values)

	np = 4;																						// set np to 4
	ddx = dx/np;																				// set ddx to the space cell width divided by np
	ddv = dv/np;																				// set ddv to the velocity cell width divided by np
	for(i1=0; i1<Nx; i1++)
	{
		x1_0 = Gridx((double)i1 - 0.5);															// set x1_0 to the value of x1 at the left edge of the i-th space cell
		for (nx=0; nx<np; nx++)
		{
			x1_val = x1_0 + nx*ddx;																// set x1_val to x1_0 plus nx increments of width ddx
			for(j1=0; j1<Nv; j1++)
			{
				v1_0 = Gridv((double)j1 - 0.5);													// set v1_0 to the value of v1 at the left edge of the j1-th velocity cell in the v1 direction
				for (nv=0; nv<np; nv++)
				{
					v1_val = v1_0 + nv*ddv;														// set v1_val to v1_0 plus nv increments of width ddv
					fM_val = f_marg2D_x1v1(U_vals, i1, j1, x1_val, v1_val);						// calculate the value of the marginal, evaluated at x_val & v1_val by using the function in the space cell i1 in the x1 direction and velocity cell j1 in the v1 direction
					fprintf(margfile_x1v1, "%11.8g  ", fM_val);									// in the file tagged as fmarg, print the value of the marginal f_M(t, x, v1)
				}
			}
			for(i2=0; i2<Nx; i2++)
			{
				x2_0 = Gridx((double)i2 - 0.5);													// set x2_0 to the value of x2 at the left edge of the i2-th velocity cell in the i2 direction
				for (nv=0; nv<np; nv++)
				{
					x2_val = x2_0 + nv*ddx;														// set x2_val to x2_0 plus nv increments of width ddv
					fM_val = f_marg2D_x1x2(U_vals, i1, i2, x1_val, x2_val);						// calculate the value of the marginal, evaluated at x1_val & x2_val by using the function in the space cell (i1, i2)
					fprintf(margfile_x1x2, "%11.8g  ", fM_val);										// in the file tagged as fmarg, print the value of the marginal f_M(t, x, v1)
				}
			}
		}
	}
	fprintf(margfile_x1v1, "\n");																// print a new line in the file tagged as margfile_x1v1
	fprintf(margfile_x1x2, "\n");																// print a new line in the file tagged as margfile_x1x2
}
