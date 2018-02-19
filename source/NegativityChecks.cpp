/* This is the source file which contains the subroutines necessary for checking where the solution
 * loses positivity.
 *
 * Functions included: computeCellAvg, FindNegVals, FindNegVals_old, CheckNegVals
 *
 *  Created on: Nov 15, 2017
 */

#include "NegativityChecks.h"																										// NegativityChecks.h is where the prototypes for the functions contained in this file are declared

double computeCellAvg(vector<double>& U_vals0, int i1, int i2, int j1, int j2, int j3)															// function to calculate the average value of the approximate function f (with DG coefficients in U) on the cell I_(i1,i2) x K_(j1,j2,j3), namely (1/cell_volume)*int_(I_(i1,i2) x K_(j1,j2,j3)) f dxdv = (1/(dx*dv^3))*int_(I_i x K_(j1,j2,j3)) f dxdv
{
	int k, nx1, nx2, nv1, nv2, nv3;																									// declare k (the location of the given cell in U), nx1, nx2 (counters for the space cell) & nv1, nv2, nv3 (counters for the velocity cell)
	double x1_0, x2_0, v1_0, v2_0, v3_0, x1_val, x2_val, v1_val, v2_val, v3_val, f_val, avg;										// declare x1_0, x2_0, v1_0, v2_0, v3_0 (to store the values in the middle of the current cell), x1_val, x2_val, v1_val, v2_val, v3_val (to store the x & v values to be evaluated at), f_val (to store the value of the function evaluated at the current x & v values) & avg (to store the average to be returned)
	k = Nv*(Nv*(Nv*(Nx*i1 + i2) + j1) + j2) + j3;																					// set k to i1*Nx*Nv^3 + i2*Nv^3 + j1*Nv^2 + j2*Nv + j3
	x1_0 = Gridx((double)i1);																										// set x1_0 to the value of x at the center of the i1th space cell
	x2_0 = Gridx((double)i2);																										// set x2_0 to the value of x at the center of the i2th space cell
	v1_0 = Gridv((double)j1);																										// set v1_0 to the value of v1 at the center of the j1th velocity cell in that direction
	v2_0 = Gridv((double)j2);																										// set v2_0 to the value of v2 at the center of the j2th velocity cell in that direction
	v3_0 = Gridv((double)j3);																										// set v3_0 to the value of v3 at the center of the j3th velocity cell in that direction
	avg = 0;																														// set avg to 0 originally

	//#pragma omp parallel for private(nx,nv1,nv2,nv3,x_val,v1_val,v2_val,v3_val,f_val) shared(U,dv,dx,vt,wt,x_0,v1_0,v2_0,v3_0,k) reduction(+:avg)
	for(nx1=0;nx1<5;nx1++)																											// loop through the five quadrature points in the x1 direction of the cell
	{
		x1_val = x1_0 + 0.5*vt[nx1]*dx;																								// set x1_val to the nx1-th quadrature point in the cell
		for(nx2=0;nx2<5;nx2++)																										// loop through the five quadrature points in the x2 direction of the cell
		{
			x2_val = x2_0 + 0.5*vt[nx2]*dx;																							// set x2_val to the nx2-th quadrature point in the cell
			for(nv1=0;nv1<5;nv1++)																									// loop through the five quadrature points in the v1 direction of the cell
			{
				v1_val = v1_0 + 0.5*vt[nv1]*dv;																						// set v1_val to the nv1-th quadrature point in the cell
				for(nv2=0;nv2<5;nv2++)																								// loop through the five quadrature points in the v2 direction of the cell
				{
					v2_val = v2_0 + 0.5*vt[nv2]*dv;																					// set v2_val to the nv2-th quadrature point in the cell
					for(nv3=0;nv3<5;nv3++)																							// loop through the five quadrature points in the v3 direction of the cell
					{
						v3_val = v3_0 + 0.5*vt[nv3]*dv;																				// set v3_val to the nv3-th quadrature point in the cell
						f_val = U_vals0[k*7+0] + U_vals0[k*7+1]*(x1_val-x1_0)/dx + U_vals0[k*7+2]*(x2_val-x2_0)/dx + U_vals0[k*7+3]*(v1_val-v1_0)/dv +
								U_vals0[k*7+4]*(v2_val-v2_0)/dv + U_vals0[k*7+5]*(v3_val-v3_0)/dv +
								U_vals0[k*7+6]*(((v1_val-v1_0)/dv)*((v1_val-v1_0)/dv)+((v2_val-v2_0)/dv)*((v2_val-v2_0)/dv)
										+((v3_val-v3_0)/dv)*((v3_val-v3_0)/dv));													// set f_val to the evaluation of the approximation at (x1_val,x2_val,v1_val,v2_val,v3_val)
						avg += wt[nx1]*wt[nx2]*wt[nv1]*wt[nv2]*wt[nv3]*f_val;														// add f(x1_val,x2_val,v1_val,v2_val,v3_val) times the quadrature weights corresponding to the current x & v values to the current quadrature value
					}
				}
			}
		}
		avg = avg*0.5*0.5*0.5*0.5*0.5;																								// multiply the quadrature result by 0.5*0.5*0.5*0.5*0.5, since each quadrature is done over intervals of width dv (three times) or dx (twice) instead of the standard interval [-1,1] (should be multiplied by dx^2*dv^3 also but then, to find the average, need to divide by dx^2*dv^3)
		return avg;																													// return the average value of f on the cell
	}
}

void FindNegVals(vector<double>& U_vals, int *NegVals, double *AvgVals)																			// function to find out the cells in which the approximation from U is negative on average and stores the cell locations in NegVals
{
	int i1, i2, j1, j2, j3, i1NNNN, i2NNN, j1NN, j2N, k, nx1, nx2, nv1, nv2, nv3;													// declare i1, i2 (the indeces of the space cell), j1, j2, j3 (the indices of the velocity cell), iNNN (to store i*Nv^3), j1NN (to store j1*Nv^2), j2N (to store j2*Nv), k (the location of the given cell in U), nx (a counter for the space cell) & nv1, nv2, nv3 (counters for the velocity cell)
	double x_val, v1_val, v2_val, v3_val, f_avg;																					// declare x1_val, x2_val, v1_val, v2_val, v3_val (to store the x & v values to be evaluated at) & f_val (to store the value of the function evaluated at the current x & v values)
	#pragma omp parallel for private(i1,i2,j1,j2,j3,k,i1NNNN,i2NNN,j1NN,j2N,f_avg) shared(U_vals,NegVals,AvgVals,size_v,Nv)
	for(i1=0;i1<Nx;i1++)																											// loop through the space cells in the x1 direction
	{
		i1NNNN = i1*Nx*size_v;																										// set i1NNNN to i1*Nx*Nv^3
		for(i2=0;i2<Nx;i2++)																										// loop through the space cells in the x2 direction
		{
			i2NNN = i2*size_v;																										// set i2NNN to i2*Nv^3
			for(j1=0;j1<Nv;j1++)																									// loop through the velocity cells in the v1 direction
			{
				j1NN = j1*Nv*Nv;																									// set j1NN to j1*Nv^2
				for(j2=0;j2<Nv;j2++)																								// loop through the velocity cells in the v2 direction
				{
					j2N = j2*Nv;																									// set j2N to j2*Nv
					for(j3=0;j3<Nv;j3++)																							// loop through the velocity cells in the v3 direction
					{
						k = i1NNNN + i2NNN + j1NN + j2N + j3;																		// set k to i1*Nx*Nv^3 + i2*Nv^3 + j1*Nv^2 + j2*Nv + j3
						NegVals[k] = 0;																								// set NegVals[k] to 0, which assumes at first that there will be no negative values in the kth cell
						f_avg = computeCellAvg(U_vals,i1,i2,j1,j2,j3);																	// calculate the average value of the approximate solution in the current cell and set it to f_avg
						AvgVals[k] = f_avg;																							// store f_avg in AvgVals[k]
						if(f_avg < 0)																								// check if this value was negative
						{
							NegVals[k] = 1;																							// if so, set NegVals[k] to 1 to indicate that there was a negative value in this cell
						}
					}
				}
			}
		}
	}
}

void FindNegVals_old(vector<double>& U_vals, int *NegVals)																						// function to find out the cells in which the approximation from U turns negative and stores the cell locations in NegVals
{
	int i1, i2, j1, j2, j3, i1NNNN, i2NNN, j1NN, j2N, k, nx1, nx2, nv1, nv2, nv3;													// declare i1, i2 (the indeces of the space cell), j1, j2, j3 (the indices of the velocity cell), iNNN (to store i*Nv^3), j1NN (to store j1*Nv^2), j2N (to store j2*Nv), k (the location of the given cell in U), nx (a counter for the space cell) & nv1, nv2, nv3 (counters for the velocity cell)
	double x1_0, x2_0, v1_0, v2_0, v3_0, x1_val, x2_val, v1_val, v2_val, v3_val, f_val;												// declare x1_0, x2_0, v1_0, v2_0, v3_0 (to store the coordinates in the middle of the current cell), x1_val, x2_val, v1_val, v2_val, v3_val (to store the x & v values to be evaluated at) & f_val (to store the value of the function evaluated at the current x & v values)
	for(i1=0;i1<Nx;i1++)																											// loop through the space cells
	{
		i1NNNN = i1*Nx*size_v;																										// set i1NNNN to i1*Nx*Nv^3
		x1_0 = Gridx((double)i1);																									// set x1_0 to the value of x1 at the center of the i1th space cell
		for(i2=0;i2<Nx;i2++)																										// loop through the space cells in the x1 direction
		{
			i2NNN = i2*size_v;																										// set i2NNN to i2*Nv^3
			x2_0 = Gridx((double)i2);																								// set x2_0 to the value of x2 at the center of the i2th space cell
			for(j1=0;j1<Nv;j1++)																									// loop through the velocity cells in the v1 direction
			{
				j1NN = j1*Nv*Nv;																									// set j1NN to j1*Nv^2
				v1_0 = Gridv((double)j1);																							// set v1_0 to the value of v1 at the center of the j1th velocity cell in that direction
				for(j2=0;j2<Nv;j2++)																								// loop through the velocity cells in the v2 direction
				{
					j2N = j2*Nv;																									// set j2N to j2*Nv
					v2_0 = Gridv((double)j2);																						// set v2_0 to the value of v2 at the center of the j2th velocity cell in that direction
					for(j3=0;j3<Nv;j3++)																							// loop through the velocity cells in the v3 direction
					{
						k = i1NNNN + i2NNN + j1NN + j2N + j3;																		// set k to i1*Nx*Nv^3 + i2*Nv^3 + j1*Nv^2 + j2*Nv + j3
						v3_0 = Gridv((double)j3);																					// set v3_0 to the value of v3 at the center of the j3th velocity cell in that direction
						NegVals[k] = 0;																								// set NegVals[k] to 0, which assumes at first that there will be no negative values in the kth cell
						for(nx1=0;nx1<5;nx1++)																						// loop through the five quadrature points in the x1 direction of the cell
						{
							x1_val = x1_0 + 0.5*vt[nx1]*dx;																			// set x1_val to the nx1-th quadrature point in the cell
							for(nx2=0;nx2<5;nx2++)																					// loop through the five quadrature points in the x2 direction of the cell
							{
								x2_val = x2_0 + 0.5*vt[nx2]*dx;																		// set x2_val to the nx2-th quadrature point in the cell
								for(nv1=0;nv1<5;nv1++)																				// loop through the five quadrature points in the v1 direction of the cell
								{
									v1_val = v1_0 + 0.5*vt[nv1]*dv;																	// set v1_val to the nv1-th quadrature point in the cell
									for(nv2=0;nv2<5;nv2++)																			// loop through the five quadrature points in the v2 direction of the cell
									{
										v2_val = v2_0 + 0.5*vt[nv2]*dv;																// set v2_val to the nv2-th quadrature point in the cell
										for(nv3=0;nv3<5;nv3++)																		// loop through the five quadrature points in the v3 direction of the cell
										{
											v3_val = v3_0 + 0.5*vt[nv3]*dv;															// set v3_val to the nv3-th quadrature point in the cell
											f_val = U_vals[k*7+0] + U_vals[k*7+1]*(x1_val-x1_0)/dx + U_vals[k*7+2]*(x2_val-x2_0)/dx + U_vals[k*7+3]*(v1_val-v1_0)/dv +
													U_vals[k*7+4]*(v2_val-v2_0)/dv + U_vals[k*7+5]*(v3_val-v3_0)/dv +
													U_vals[k*7+6]*(((v1_val-v1_0)/dv)*((v1_val-v1_0)/dv)+((v2_val-v2_0)/dv)*((v2_val-v2_0)/dv)
															+((v3_val-v3_0)/dv)*((v3_val-v3_0)/dv));								// set f_val to the evaluation of the approximation at (x1_val,x2_val,v1_val,v2_val,v3_val)
											if(f_val < 0)																			// check if this value was negative
											{
												NegVals[k] = 1;																		// if so, set NegVals[k] to 1 to indicate that there was a negative value in this cell
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}

void CheckNegVals(int *NegVals, double *AvgVals)																			// function to find out the cells in which the approximation from U turns negative and stores the cell locations in NegVals
{
	int i1, i2, j1, j2, j3, i1NNNN, i2NNN, j1NN, j2N, k;																			// declare i1, i2 (the index of the space cell), j1, j2, j3 (the indices of the velocity cell), iNNN (to store i*Nv^3), j1NN (to store j1*Nv^2), j2N (to store j2*Nv), k (the location of the given cell in U), nx (a counter for the space cell) & nv1, nv2, nv3 (counters for the velocity cell)
	for(i1=0;i1<Nx;i1++)																											// loop through the space cells in the x1 direction
	{
		i1NNNN = i1*Nx*size_v;																										// set i1NNNN to i1*Nx*Nv^3
		for(i2=0;i2<Nx;i2++)																										// loop through the space cells in the x2 direction
		{
			i2NNN = i2*size_v;																										// set i2NNN to i2*Nv^3
			for(j1=0;j1<Nv;j1++)																									// loop through the velocity cells in the v1 direction
			{
				j1NN = j1*Nv*Nv;																									// set j1NN to j1*Nv^2
				for(j2=0;j2<Nv;j2++)																								// loop through the velocity cells in the v2 direction
				{
					j2N = j2*Nv;																									// set j2N to j2*Nv
					for(j3=0;j3<Nv;j3++)																							// loop through the velocity cells in the v3 direction
					{
						k = i1NNNN + i2NNN + j1NN + j2N + j3;																		// set k to i1*Nx*Nv^3 + i2*Nv3 + j1*Nv^2 + j2*Nv + j3
						printf("%d ~ (%d, %d, %d, %d, %d) f_avg = %g \n", NegVals[k], i1, i2, j1, j2, j3, AvgVals[k]);				// print the value of NegVals[k], the current cell indices and the value of the average value on the cell in the output file
						printf(" \n");
					}
				}
			}
		}
	}
}
