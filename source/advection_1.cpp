/* This is the source file which contains the subroutines necessary for solving the collisionless advection
 * problem resulting from time-splitting, as well as some which are necessary for subroutines for calculating
 * moments or entropy, etc.
 *
 * Functions included: Gridv, Gridx, rho_x, rho, computePhi_x_0, computePhi, PrintPhiVals, computeC_rho, Int_Int_rho
 * Int_Int_rho1st, Int_E, Int_E1st, Int_E2nd, Int_fE, I1, I2, I3, I5, computeH, RK3
 *
 */

#include "advection_1.h"																					// advection_1.h is where the prototypes for the functions contained in this file are declared and any variables defined here to be used throughout the other files are declared as external

double wt[5]={0.5688888888888889, 0.4786286704993665, 0.4786286704993665,0.2369268850561891, 0.2369268850561891};				// weights for Gaussian quadrature
double vt[5]={0., -0.5384693101056831,0.5384693101056831,-0.9061798459386640,0.9061798459386640};								// node values for Gaussian quadrature over the interval [-1,1]

double Gridv(double m){ //v in [-Lv,Lv]
	return (-Lv+(m+0.5)*dv);
}

double Gridx(double m){ // x in [0,Lx]  (returns the x value at the mth discrete space-step
	return (m+0.5)*dx;
}

/* 1D
double rho_x(double x, vector<double>& U, int i) // for x in I_i			// MAY NEVER BE CALLED
{
  int j, k;
  double tmp=0.;
  //#pragma omp parallel for shared(U) reduction(+:tmp)
  for(j=0;j<size_v;j++){
	k=i*size_v + j;
	tmp += U[k*7+0] + U[k*7+1]*(x-Gridx((double)i))/dx + U[k*7+6]/4.;
  }
  
  return dv*dv*dv*tmp;
}

double rho(vector<double>& U, int i) //  \int f(x,v) dxdv in I_i * K_j		// MAY NEVER BE CALLED
{
  int j, k;
  double tmp=0.;
 // #pragma omp parallel for shared(U) reduction(+:tmp)
  for(j=0;j<size_v;j++){
	k=i*size_v + j;
	tmp += U[k*7+0] + U[k*7+6]/4.;
  }
  
  return dx*dv*dv*dv*tmp;
}

double computePhi_x_0(vector<double>& U) // compute the constant coefficient of x in phi, which is actually phi_x(0) (Calculate C_E in the paper -between eq. 52 & 53?)
{
	int i, j, k, m, q;
	double tmp=0.;

	//#pragma omp parallel for private(j,q,m,k) shared(U) reduction(+:tmp) //reduction may change the final result a little bit
	for(j=0;j<size_v;j++){
		//j = j1*Nv*Nv + j2*Nv + j3;
		for(q=0;q<Nx;q++){
			for(m=0;m<q;m++){			
				k=m*Nx*size_v + j;	// assuming k = i1*Nx*size_v + i2*size_v + j (choose i2 = 0 as should be same for all i2 in 1D)
				tmp += U[k*7+0] + U[k*7+6]/4.;
			}
			k=q*Nx*size_v + j;	// assuming k = i1*Nx*size_v + i2*size_v + j (choose i2 = 0 as should be same for all i2 in 1D)
			tmp += 0.5*(U[k*7+0] + U[k*7+6]/4.) - U[k*7+1]/12.;
		}
	}
	tmp = tmp*scalev*dx*dx;

	return 0.5*Lx - tmp/Lx;
}

double computeC_rho(vector<double>& U, int i) // sum_m=0..i-1 int_{I_m} rho(z)dz   (Calculate integral of rho_h(z,t) from 0 to x as in eq. 53)
{
	double retn=0.;
	int j, k, m;
	for (m=0;m<i;m++){ //BUG: was "m < i-1"
		for(j=0;j<size_v;j++){
			k = m*Nx*size_v + j;	// assuming k = i1*Nx*size_v + i2*size_v + j (choose i2 = 0 as should be same for all i2 in 1D)
			retn += U[k*7+0] + U[k*7+6]/4.;
		}
	}
	
	retn *= dx*scalev;
	return retn;
}

double Int_Int_rho(vector<double>& U, int i) // \int_{I_i} [ \int^{x}_{x_i-0.5} rho(z)dz ] dx
{
  int j, k;
  double retn=0.;
  for(j=0;j<size_v;j++){
    k=i*Nx*size_v + j;	// assuming k = i1*Nx*size_v + i2*size_v + j (choose i2 = 0 as should be same for all i2 in 1D)
    retn += 0.5*(U[k*7+0] + U[k*7+6]/4.) - U[k*7+1]/12.;
  }
  
  return retn*dx*dx*scalev;  
}

double Int_Int_rho1st(double *U, int i)// \int_{I_i} [(x-x_i)/delta_x * \int^{x}_{x_i-0.5} rho(z)dz ] dx			// MAY NEVER BE CALLED
{
 int j, k;
 double retn=0.;
 for(j=0;j<size_v;j++){
    k=i*Nx*size_v + j;	// assuming k = i1*Nx*size_v + i2*size_v + j (choose i2 = 0 as should be same for all i2 in 1D)
    retn += (U[k*7+0] + U[k*7+6]/4.)/12.;
  }
  return retn*dx*dx*scalev; 
}
*/
/*double Int_Cumulativerho(double **U, int i)// \int_{I_i} [ \int^{x}_{0} rho(z)dz ] dx
{
  double retn=0., cp, tmp;
  cp = computeC_rho(U,i);
  tmp = Int_Int_rho(U,i);
  
  retn = dx*cp + tmp;
  return retn;  
}

double Int_Cumulativerho_sqr(double **U, int i)// \int_{I_i} [ \int^{x}_{0} rho(z)dz ]^2 dx
{
  double retn=0., cp, tmp1, tmp2, tmp3, c1=0., c2=0.;
  int j, k;
  cp = computeC_rho(U,i);
  tmp1 = cp*cp*dx;
  tmp2 = 2*cp*Int_Int_rho(U,i);
  for(j=0;j<size_v;j++){
    k=i*size_v + j;
    c1 += U[k][0] + U[k][5]/4.;
    c2 += U[k][1];
  }
  c2 *= dx/2.;
  tmp3 = pow(dv, 6)* ( c1*c1*dx*dx*dx/3. + c2*c2*dx/30. + c1*c2*(dx*Gridx((double)i)/6. - dx*dx/4.) );
  retn = tmp1 + tmp2 + tmp3;
  return retn;  
}*/
/* 1D
double Int_E(vector<double>& U, int i) // \int_i E dx      // Function to calculate the integral of E_h w.r.t. x over the interval I_i = [x_(i-1/2), x_(i+1/2))
{
	int m, j, k;
	double tmp=0., result;
	//#pragma omp parallel for shared(U) reduction(+:tmp)
	for(j=0;j<size_v;j++){
		for(m=0;m<i;m++){	
			k=m*Nx*size_v + j;	// assuming k = i1*Nx*size_v + i2*size_v + j (choose i2 = 0 as should be same for all i2 in 1D)
			tmp += U[k*7+0] + U[k*7+6]/4.;
		}
		k=i*Nx*size_v + j;	// assuming k = i1*Nx*size_v + i2*size_v + j (choose i2 = 0 as should be same for all i2 in 1D)
		tmp += 0.5*(U[k*7+0] + U[k*7+6]/4.) - U[k*7+1]/12.;		
	}

	//ce = computePhi_x_0(U);
	result = -ce*dx - tmp*dx*dx*scalev + Gridx((double)i)*dx;

	return result;
}

double Int_E1st(vector<double>& U, int i) // \int_i E*(x-x_i)/delta_x dx
{
	int j, k;
	double tmp=0., result;
	//#pragma omp parallel for reduction(+:tmp)
	for(j=0;j<size_v;j++){
		k=i*Nx*size_v + j;	// assuming k = i1*Nx*size_v + i2*size_v + j (choose i2 = 0 as should be same for all i2 in 1D)
		tmp += U[k*7+0] + U[k*7+6]/4.;
	}
	tmp = tmp*scalev;
	
	result = (1-tmp)*dx*dx/12.;

	return result;
}

double Int_fE_1DE(vector<double>& U, int i, int j) // \int f * E(f) dxdv on element I_i * K_j
{
	double retn=0.;
	int k;
	k = i*Nx*size_v + j;	// assuming k = i1*Nx*size_v + i2*size_v + j (choose i2 = 0 as should be same for all i2 in 1D)

	//retn = (U[k][0] + U[k][5]/4.)*Int_E(U,i) + U[k][1]*Int_E1st(U,i);
	retn = (U[k*7+0] + U[k*7+6]/4.)*intE[i] + U[k*7+1]*intE1[i];
	return retn*scalev;
}
*/

double Int_fE1(vector<double>& U, int i, int j) // \int f * E_1(f) dxdv on element I_i * K_j
{
	double retn=0.;
	int k;
	k = i*size_v + j;	// assuming k = i1*Nx*size_v + i2*size_v + j = i*size_v + j, for i = i_mod = i1*Nx + i2 (choose i2 = 0 as should be same for all i2 in 1D)

	//retn = (U[k*7+0] + U[k*7+6]/4.)*IE_X[i] + U[k*7+1]*intE1[i];	//1D E
	retn = (U[k*7+0] + U[k*7+6]/4.)*IE_X[i] + U[k*7+1]*IXE_X[i] + U[k*7+2]*IYE_X[i];	//2D E
	return retn*scalev;
}

double Int_fE2(vector<double>& U, int i, int j) // \int f * E_2(f) dxdv on element I_i * K_j
{
	double retn=0.;
	int k;
	k = i*size_v + j;	// assuming k = i1*Nx*size_v + i2*size_v + j = i*size_v + j, for i = i_mod = i1*Nx + i2 (choose i2 = 0 as should be same for all i2 in 1D)

	retn = (U[k*7+0] + U[k*7+6]/4.)*IE_Y[i] + U[k*7+1]*IXE_Y[i] + U[k*7+2]*IYE_Y[i];	//2D E
	return retn*scalev;
}

/* 1D
double Int_E2nd(vector<double>& U, int i) // \int_i E* [(x-x_i)/delta_x]^2 dx
{
    int m, j, j1, j2, j3, k;
    double c1=0., c2=0., result;
  
    //cp = computeC_rho(U,i); ce = computePhi_x_0(U);

    for(j=0;j<size_v;j++){
	    k=i*Nx*size_v + j;	// assuming k = i1*Nx*size_v + i2*size_v + j (choose i2 = 0 as should be same for all i2 in 1D)
	    c1 += U[k*7+0] + U[k*7+6]/4.;
	    c2 += U[k*7+1];
    }
    c2 *= dx/2.;				
    
    result = (-cp[i] - ce+scalev*(c1*Gridx(i-0.5) + 0.25*c2))*dx/12. + (1-scalev*c1)*dx*Gridx((double)i)/12. - scalev*c2*dx/80.; //BUG: missed -cp

    return result;
}
*/

double I1_1D(vector<double>& U_vals0, int k, int l) // Calculate the first inegtral in H_(i,j), namely \int v1*f*phi_x dxdv
{
  double result;
  int i1, i2, j1, j2, j3; // k=i1*Nx*Nv^3 + i2*Nv^3 + (j1*Nv*Nv + j2*Nv + j3)								// declare i1 & i2 (the space cell coordinatex), j1, j2, j3 (the coordinates of the velocity cell)
  int j_mod = k%size_v;																						// declare and calculate j_mod (the remainder when k is divided by size_v = Nv^3 - used to help determine the values of i1, i2, j1, j2 & j3 from the value of k)
  int i_mod = (k-j_mod)/size_v;																				// declare and calculate i_mod (the remainder value when k has j_mod subtracted and divided through by size_v - used to help determine the values of i1 & i2)
  j3 = j_mod%Nv;																							// calculate j3 for the given k
  j2 = ((j_mod-j3)%(Nv*Nv))/Nv;																				// calculate j2 for the given k
  j1 = (j_mod-j3-j2*Nv)/(Nv*Nv);																			// calculate j1 for the given k
  i2 = i_mod%Nx;																							// calculate i2 for the given k
  i1 = (i_mod - i2)/Nx;																						// calculate i1 for the given k

  if(l==1) result = dv*dv*dv*( Gridv((double)j1)*U_vals0[k*7+0] + dv*U_vals0[k*7+3]/12. + U_vals0[k*7+6]*Gridv((double)j1)/4.);
  else result=0.;

  return result;
}

double I1(vector<double>& U_vals0, int k, int l) // Calculate the first inegtral in H_(i,j), namely \int v1*f*phi_x dxdv
{
	double result;
	int i1, i2, j1, j2, j3; // k=i1*Nx*Nv^3 + i2*Nv^3 + (j1*Nv*Nv + j2*Nv + j3)								// declare i1 & i2 (the space cell coordinatex), j1, j2, j3 (the coordinates of the velocity cell)
	int j_mod = k%size_v;																						// declare and calculate j_mod (the remainder when k is divided by size_v = Nv^3 - used to help determine the values of i1, i2, j1, j2 & j3 from the value of k)
	int i_mod = (k-j_mod)/size_v;																				// declare and calculate i_mod (the remainder value when k has j_mod subtracted and divided through by size_v - used to help determine the values of i1 & i2)
	j3 = j_mod%Nv;																							// calculate j3 for the given k
	j2 = ((j_mod-j3)%(Nv*Nv))/Nv;																				// calculate j2 for the given k
	j1 = (j_mod-j3-j2*Nv)/(Nv*Nv);																			// calculate j1 for the given k
	i2 = i_mod%Nx;																							// calculate i2 for the given k
	i1 = (i_mod - i2)/Nx;																						// calculate i1 for the given k
	double v_j1, v_j2;																						// declare v_j1 & v_j2 (the value at the center of the velocity grid in the v1 & v2 directions)

	/* NEED TO MULTIPLY ALL THESE BY dx AFTER I FUNCTIONS IMPLEMENTED! */
	if(l==1)
	{
	  v_j1 = Gridv((double)j1);																				// set v_j1 to value of v_1 at the center of cell j_1
	  result = scalev*(U_vals0[k*7+0]*v_j1 + U_vals0[k*7+3]*dv/12. + U_vals0[k*7+6]*v_j1/4.);
	}
	else if(l==2)
	{
	  v_j2 = Gridv((double)j2);
	  result = scalev*(U_vals0[k*7+0]*v_j2 + U_vals0[k*7+4]*dv/12. + U_vals0[k*7+6]*v_j2/4.);
	}
	else result=0.;

	return result;
}

double I2_1D(vector<double>& U_vals0, int k, int l) // Calculate the fourth integral in H_(i,j), namely \int E*f*phi_v1 dxdv
{
  double result;
  int j_mod = k%size_v;																						// declare and calculate j_mod (the remainder when k is divided by size_v = Nv^3 - used to help determine the values of i1, i2, j1, j2 & j3 from the value of k)
  int i_mod = (k-j_mod)/size_v;		// k = i_mod*Nv^3 + j_mod												// declare and calculate i_mod (the remainder value when k has j_mod subtracted and divided through by size_v - used to help determine the values of i1 & i2)

  //if(l==3) result = Int_fE_1DE(U_vals0,i1,j_mod)/dv;	//1D E
  //else if(l==6) result = U_vals0[k*7+3]*dv*dv*intE[i1]/6.;	//1D E
  if(l==3) result = Int_fE1(U_vals0,i_mod,j_mod)/dv;	//2D E
  else if(l==6) result = U_vals0[k*7+3]*dv*dv*IE_X[i_mod]/6.;	//2D E
  else result = 0.;
  
  return result;
}

double I2(vector<double>& U_vals0, int k, int l) // Calculate the fourth integral in H_(i,j), namely \int E*f*phi_v1 dxdv
{
	double result;
	int j_mod = k%size_v;																						// declare and calculate j_mod (the remainder when k is divided by size_v = Nv^3 - used to help determine the values of i1, i2, j1, j2 & j3 from the value of k)
	int i_mod = (k-j_mod)/size_v;		// k = i_mod*Nv^3 + j_mod												// declare and calculate i_mod (the remainder value when k has j_mod subtracted and divided through by size_v - used to help determine the values of i1 & i2)

	if(l==3)
	{
	  result = Int_fE1(U_vals0,i_mod,j_mod)/dv;	//2D E
	}
	else if(l==4)
	{
	  result = Int_fE2(U_vals0,i_mod,j_mod)/dv;	//2D E
	}
	else if(l==6)
	{
	  result = dv*dv*(U_vals0[k*7+3]*IE_X[i_mod] + U_vals0[k*7+4]*IE_Y[i_mod])/6.;	//2D E
	}
	else result = 0.;

	return result;
}

double I3_1D(vector<double>& U_vals0, int k, int l) 																			// Calculate the difference of the second and third integrals in H_(i,j), namely \int_j v1*gh*phi dv at interface x=x_i+1/2 - \int_j v1*gh*phi dv at interface x=x_i-1/2
{
	double result, ur, ul;																					// declare result (the result of the integral to be returned), ur (used in the evaluation of gh^+/- on the right space cell edge) & ul (used in the evaluation of gh^+/- on the left space cell edge)
	int i1, i2, j1, j2, j3, iil, iir, kkl, kkr; 															// declare i1 & i2 (the space cell coordinatex), j1, j2, j3 (the coordinates of the velocity cell), iil (the cell from which the flux is flowing on the left of cell i in space), iir (the cell from which the flux is flowing on the right of cell i in space), kkl (the global index of the cell with coordinate (iil, j1, j2, j3)) & kkr (the global index of the cell with coordinate (iir, j1, j2, j3))
	int j_mod = k%size_v;																					// declare and calculate j_mod (the remainder when k is divided by size_v = Nv^3 - used to help determine the values of i1, i2, j1, j2 & j3 from the value of k)
	int i_mod = (k-j_mod)/size_v;																			// declare and calculate i_mod (the remainder value when k has j_mod subtracted and divided through by size_v - used to help determine the values of i1 & i2)
	j3 = j_mod%Nv;																							// calculate j3 for the given k
	j2 = ((j_mod-j3)%(Nv*Nv))/Nv;																			// calculate j2 for the given k
	j1 = (j_mod-j3-j2*Nv)/(Nv*Nv);																			// calculate j1 for the given k
	i2 = i_mod%Nx;																							// calculate i2 for the given k
	i1 = (i_mod - i2)/Nx;																					// calculate i1 for the given k

	if(j1<Nv/2)																								// do this if j1 < Nv/2 (so that the velocity in the v1 direction is negative)
	{
//		iir=i2+1; iil=i2; 																					// set iir to the value of i2+1 and iil to the value of i2 (as here the flow of information is from right to left so that gh^+ must be used at the cell edges)
		iir=i1+1; iil=i1; 																					// set iir to the value of i2+1 and iil to the value of i2 (as here the flow of information is from right to left so that gh^+ must be used at the cell edges)
		if(iir==Nx)iir=0; //periodic bc																		// if iir = Nx (the maximum value that can be obtained, since i = 0,1,...,Nx-1) and so this cell is at the right boundary, requiring information from the non-existent cell with space index Nx, since there are periodic boundary conditions, set iir = 0 and use the cell with space index 0 (i.e. the cell at the left boundary)
//		kkr=i1*Nx*size_v + iir*size_v + j_mod; 																// calculate the value of kkr for this value of iir
		kkr=iir*Nx*size_v + i2*size_v + j_mod; 																// calculate the value of kkr for this value of iir
		kkl=k;																								// set kkl to k (since iil = i2)
		ur = -U_vals0[kkr*7+1]; 																					// set ur to the negative of the coefficient of the basis function with shape l which is non-zero in the cell with global index kkr (which corresponds to the evaluation of gh^+ at the right boundary and -ve since v_1 < 0 in here)
		ul = -U_vals0[kkl*7+1];																					// set ul to the negative of the coefficient of the basis function with shape l which is non-zero in the cell with global index kkl	(which corresponds to the evaluation of gh^+ at the left boundary and -ve since phi < 0 here)
	}
	else																									// do this if j1 >= Nv/2 (so that the velocity in the v1 direction is non-negative)
	{
//		iir=i2; iil=i2-1;																					// set iir to the value of i2 and iil to the value of i2-1 (as here the flow of information is from left to right so that gh^- must be used at the cell edges)
		iir=i1; iil=i1-1;																					// set iir to the value of i2 and iil to the value of i2-1 (as here the flow of information is from left to right so that gh^- must be used at the cell edges)
		if(iil==-1)iil=Nx-1; // periodic bc																	// if iil = -1 (the minimum value that can be obtained, since i = 0,1,...,Nx-1) and so this cell is at the left boundary, requiring information from the non-existent cell with space index -1, since there are periodic boundary conditions, set iil = Nx-1 and use the cell with space index Nx-1 (i.e. the cell at the right boundary)
		kkr=k; 																								// set kkr to k (since iir = i1)
//		kkl=i1*Nx*size_v + iil*size_v + j_mod; 																// calculate the value of kkl for this value of iil
		kkl=iil*Nx*size_v + i2*size_v + j_mod; 																// calculate the value of kkl for this value of iil
		ur = U_vals0[kkr*7+1];																					// set ur to the value of the coefficient of the basis function with shape l which is non-zero in the cell with global index kkr (which corresponds to the evaluation of gh^- at the right boundary and +ve since v_1 >= 0 in here)
		ul = U_vals0[kkl*7+1];																					// set ul to the value of the coefficient of the basis function with shape l which is non-zero in the cell with global index kkl (which corresponds to the evalutaion of gh^- at the left boundary and +ve since v_r >= 0 in here)
	}

	if(l==0)result = dv*dv*dv*( (U_vals0[kkr*7+0]+0.5*ur - U_vals0[kkl*7+0]-0.5*ul)*Gridv((double)j1) + (U_vals0[kkr*7+3]-U_vals0[kkl*7+3])*dv/12. + (U_vals0[kkr*7+6]-U_vals0[kkl*7+6])*Gridv((double)j1)/4.);					// calculate \int_j v1*gh*phi dv at interface x=x_i+1/2 - \int_j v1*gh*phi dv at interface x=x_i-1/2 for the basis function with shape 0 (i.e. constant) which is non-zero in the cell with global index k
	if(l==1)result = 0.5*dv*dv*dv*( (U_vals0[kkr*7+0]+0.5*ur + U_vals0[kkl*7+0]+0.5*ul)*Gridv((double)j1) + (U_vals0[kkr*7+3]+U_vals0[kkl*7+3])*dv/12. + (U_vals0[kkr*7+6]+U_vals0[kkl*7+6])*Gridv((double)j1)/4.);				// calculate \int_j v1*gh*phi dv at interface x=x_i+1/2 - \int_j v1*gh*phi dv at interface x=x_i-1/2 for the basis function with shape 1 (i.e. linear in x) which is non-zero in the cell with global index k
	if(l==3)result = dv*dv*(( (U_vals0[kkr*7+0]-U_vals0[kkl*7+0])*dv*dv + (ur-ul)*0.5*dv*dv + (U_vals0[kkr*7+3]-U_vals0[kkl*7+3])*dv*Gridv((double)j1))/12. + (U_vals0[kkr*7+6]-U_vals0[kkl*7+6])*dv*dv*19./720.);				// calculate \int_j v1*gh*phi dv at interface x=x_i+1/2 - \int_j v1*gh*phi dv at interface x=x_i-1/2 for the basis function with shape 2 (i.e. linear in v_1) which is non-zero in the cell with global index k
	if(l==4)result = (U_vals0[kkr*7+4]-U_vals0[kkl*7+4])*Gridv((double)j1)*dv*dv*dv/12.;																												// calculate \int_j v1*gh*phi dv at interface x=x_i+1/2 - \int_j v1*gh*phi dv at interface x=x_i-1/2 for the basis function with shape 3 (i.e. linear in v_2) which is non-zero in the cell with global index k
	if(l==5)result = (U_vals0[kkr*7+5]-U_vals0[kkl*7+5])*Gridv((double)j1)*dv*dv*dv/12.;																												// calculate \int_j v1*gh*phi dv at interface x=x_i+1/2 - \int_j v1*gh*phi dv at interface x=x_i-1/2 for the basis function with shape 4 (i.e. linear in v_3) which is non-zero in the cell with global index k
	if(l==6)result = dv*dv*dv*((U_vals0[kkr*7+0] + 0.5*ur - U_vals0[kkl*7+0]-0.5*ul)*Gridv((double)j1)/4. + (U_vals0[kkr*7+3]-U_vals0[kkl*7+3])*dv*19./720. + (U_vals0[kkr*7+6]-U_vals0[kkl*7+6])*Gridv((double)j1)*19./240.);	// calculate \int_j v1*gh*phi dv at interface x=x_i+1/2 - \int_j v1*gh*phi dv at interface x=x_i-1/2 for the basis function with shape 5 (i.e. modulus of v) which is non-zero in the cell with global index k

	return result;
}

double I3_x1(vector<double>& U_vals0, int k, int l) 																			// Calculate the difference of the second and third integrals in H_(i,j), namely \int_j v1*gh*phi dv at interface x=x_i+1/2 - \int_j v1*gh*phi dv at interface x=x_i-1/2
{
	double result, ur, ul;																					// declare result (the result of the integral to be returned), ur (used in the evaluation of gh^+/- on the right space cell edge) & ul (used in the evaluation of gh^+/- on the left space cell edge)
	int i1, i2, j1, j2, j3, iil, iir, kkl, kkr; 															// declare i1 & i2 (the space cell coordinatex), j1, j2, j3 (the coordinates of the velocity cell), iil (the cell from which the flux is flowing on the left of cell i in space), iir (the cell from which the flux is flowing on the right of cell i in space), kkl (the global index of the cell with coordinate (iil, j1, j2, j3)) & kkr (the global index of the cell with coordinate (iir, j1, j2, j3))
	int j_mod = k%size_v;																					// declare and calculate j_mod (the remainder when k is divided by size_v = Nv^3 - used to help determine the values of i1, i2, j1, j2 & j3 from the value of k)
	int i_mod = (k-j_mod)/size_v;																			// declare and calculate i_mod (the remainder value when k has j_mod subtracted and divided through by size_v - used to help determine the values of i1 & i2)
	j3 = j_mod%Nv;																							// calculate j3 for the given k
	j2 = ((j_mod-j3)%(Nv*Nv))/Nv;																			// calculate j2 for the given k
	j1 = (j_mod-j3-j2*Nv)/(Nv*Nv);																			// calculate j1 for the given k
	i2 = i_mod%Nx;																							// calculate i2 for the given k
	i1 = (i_mod - i2)/Nx;																					// calculate i1 for the given k
	double v_j1 = Gridv((double)j1);																		// declare v_j1 and set it to value of v_1 at the center of cell j_1

	if(j1<Nv/2)																								// do this if j1 < Nv/2 (so that the velocity in the v1 direction is negative)
	{
//		iir=i2+1; iil=i2; 																					// set iir to the value of i2+1 and iil to the value of i2 (as here the flow of information is from right to left so that gh^+ must be used at the cell edges)
		iir=i1+1; iil=i1; 																					// set iir to the value of i2+1 and iil to the value of i2 (as here the flow of information is from right to left so that gh^+ must be used at the cell edges)
		if(iir==Nx)iir=0; //periodic bc																		// if iir = Nx (the maximum value that can be obtained, since i = 0,1,...,Nx-1) and so this cell is at the right boundary, requiring information from the non-existent cell with space index Nx, since there are periodic boundary conditions, set iir = 0 and use the cell with space index 0 (i.e. the cell at the left boundary)
//		kkr=i1*Nx*size_v + iir*size_v + j_mod; 																// calculate the value of kkr for this value of iir
		kkr=iir*Nx*size_v + i2*size_v + j_mod; 																// calculate the value of kkr for this value of iir
		kkl=k;																								// set kkl to k (since iil = i2)
		ur = -U_vals0[kkr*7+1]; 																			// set ur to the negative of the coefficient of the basis function with shape l which is non-zero in the cell with global index kkr (which corresponds to the evaluation of gh^+ at the right boundary and -ve since v_1 < 0 in here)
		ul = -U_vals0[kkl*7+1];																				// set ul to the negative of the coefficient of the basis function with shape l which is non-zero in the cell with global index kkl	(which corresponds to the evaluation of gh^+ at the left boundary and -ve since phi < 0 here)
	}
	else																									// do this if j1 >= Nv/2 (so that the velocity in the v1 direction is non-negative)
	{
//		iir=i2; iil=i2-1;																					// set iir to the value of i2 and iil to the value of i2-1 (as here the flow of information is from left to right so that gh^- must be used at the cell edges)
		iir=i1; iil=i1-1;																					// set iir to the value of i2 and iil to the value of i2-1 (as here the flow of information is from left to right so that gh^- must be used at the cell edges)
		if(iil==-1)iil=Nx-1; // periodic bc																	// if iil = -1 (the minimum value that can be obtained, since i = 0,1,...,Nx-1) and so this cell is at the left boundary, requiring information from the non-existent cell with space index -1, since there are periodic boundary conditions, set iil = Nx-1 and use the cell with space index Nx-1 (i.e. the cell at the right boundary)
		kkr=k; 																								// set kkr to k (since iir = i1)
//		kkl=i1*Nx*size_v + iil*size_v + j_mod; 																// calculate the value of kkl for this value of iil
		kkl=iil*Nx*size_v + i2*size_v + j_mod; 																// calculate the value of kkl for this value of iil
		ur = U_vals0[kkr*7+1];																				// set ur to the value of the coefficient of the basis function with shape l which is non-zero in the cell with global index kkr (which corresponds to the evaluation of gh^- at the right boundary and +ve since v_1 >= 0 in here)
		ul = U_vals0[kkl*7+1];																				// set ul to the value of the coefficient of the basis function with shape l which is non-zero in the cell with global index kkl (which corresponds to the evalutaion of gh^- at the left boundary and +ve since v_r >= 0 in here)
	}
  
	/* NEED TO MULTIPLY ALL THESE BY dx AFTER I FUNCTIONS IMPLEMENTED! */
	if(l==0)
	{
		result = scalev*(((U_vals0[kkr*7+0] - U_vals0[kkl*7+0]) + 0.5*(ur - ul) + 0.25*(U_vals0[kkr*7+6] - U_vals0[kkl*7+6]))*v_j1 + (U_vals0[kkr*7+3] - U_vals0[kkl*7+3])*dv/12.);					// calculate \int_j v1*gh*phi dv at interface x=x_i+1/2 - \int_j v1*gh*phi dv at interface x=x_i-1/2 for the basis function with shape 0 (i.e. constant) which is non-zero in the cell with global index k
	}
	else if(l==1)
	{
		result = 0.5*scalev*(((U_vals0[kkr*7+0] + U_vals0[kkl*7+0]) + 0.5*(ur + ul) + 0.25*(U_vals0[kkr*7+6] + U_vals0[kkl*7+6]))*v_j1 + (U_vals0[kkr*7+3] + U_vals0[kkl*7+3])*dv/12.);					// calculate \int_j v1*gh*phi dv at interface x=x_i+1/2 - \int_j v1*gh*phi dv at interface x=x_i-1/2 for the basis function with shape 0 (i.e. constant) which is non-zero in the cell with global index k
	}
	else if(l==2)
	{
		result = scalev*(U_vals0[kkr*7+2] - U_vals0[kkl*7+2])*v_j1/12.;																												// calculate \int_j v1*gh*phi dv at interface x=x_i+1/2 - \int_j v1*gh*phi dv at interface x=x_i-1/2 for the basis function with shape 3 (i.e. linear in v_2) which is non-zero in the cell with global index k
	}
	else if(l==3)
	{
		result = scalev*(((U_vals0[kkr*7+0] - U_vals0[kkl*7+0]) + 0.5*(ur - ul) + (U_vals0[kkr*7+6] - U_vals0[kkl*7+6])*19./60.)*dv + (U_vals0[kkr*7+3] - U_vals0[kkl*7+3])*v_j1)/12.;				// calculate \int_j v1*gh*phi dv at interface x=x_i+1/2 - \int_j v1*gh*phi dv at interface x=x_i-1/2 for the basis function with shape 2 (i.e. linear in v_1) which is non-zero in the cell with global index k
	}
	else if(l==4)
	{
		result = scalev*(U_vals0[kkr*7+4] - U_vals0[kkl*7+4])*v_j1/12.;																												// calculate \int_j v1*gh*phi dv at interface x=x_i+1/2 - \int_j v1*gh*phi dv at interface x=x_i-1/2 for the basis function with shape 3 (i.e. linear in v_2) which is non-zero in the cell with global index k
	}
	else if(l==5)
	{
		result = scalev*(U_vals0[kkr*7+5] - U_vals0[kkl*7+5])*v_j1/12.;																												// calculate \int_j v1*gh*phi dv at interface x=x_i+1/2 - \int_j v1*gh*phi dv at interface x=x_i-1/2 for the basis function with shape 4 (i.e. linear in v_3) which is non-zero in the cell with global index k
	}
	else if(l==6)
	{
		result = scalev*(((U_vals0[kkr*7+0] - U_vals0[kkl*7+0]) + 0.5*(ur - ul) + (U_vals0[kkr*7+6] - U_vals0[kkl*7+6])*19./60.)*v_j1 + (U_vals0[kkr*7+3] - U_vals0[kkl*7+3])*dv*19./180.)/4.;	// calculate \int_j v1*gh*phi dv at interface x=x_i+1/2 - \int_j v1*gh*phi dv at interface x=x_i-1/2 for the basis function with shape 5 (i.e. modulus of v) which is non-zero in the cell with global index k
	}

	return result;
}

double I3_x2(vector<double>& U_vals0, int k, int l) 																			// Calculate the difference of the second and third integrals in H_(i,j), namely \int_j v1*gh*phi dv at interface x=x_i+1/2 - \int_j v1*gh*phi dv at interface x=x_i-1/2
{
	double result, ur, ul;																					// declare result (the result of the integral to be returned), ur (used in the evaluation of gh^+/- on the right space cell edge) & ul (used in the evaluation of gh^+/- on the left space cell edge)
	int i1, i2, j1, j2, j3, iil, iir, kkl, kkr; 															// declare i1 & i2 (the space cell coordinatex), j1, j2, j3 (the coordinates of the velocity cell), iil (the cell from which the flux is flowing on the left of cell i in space), iir (the cell from which the flux is flowing on the right of cell i in space), kkl (the global index of the cell with coordinate (iil, j1, j2, j3)) & kkr (the global index of the cell with coordinate (iir, j1, j2, j3))
	int j_mod = k%size_v;																					// declare and calculate j_mod (the remainder when k is divided by size_v = Nv^3 - used to help determine the values of i1, i2, j1, j2 & j3 from the value of k)
	int i_mod = (k-j_mod)/size_v;																			// declare and calculate i_mod (the remainder value when k has j_mod subtracted and divided through by size_v - used to help determine the values of i1 & i2)
	j3 = j_mod%Nv;																							// calculate j3 for the given k
	j2 = ((j_mod-j3)%(Nv*Nv))/Nv;																			// calculate j2 for the given k
//	j1 = (j_mod-j3-j2*Nv)/(Nv*Nv);																			// calculate j1 for the given k
	i2 = i_mod%Nx;																							// calculate i2 for the given k
	i1 = (i_mod - i2)/Nx;																					// calculate i1 for the given k
	double v_j2 = Gridv((double)j2);																		// declare v_j2 and set it to value of v_2 at the center of cell j_2

	if(j2<Nv/2)																								// do this if j1 < Nv/2 (so that the velocity in the v1 direction is negative)
	{
		iir=i2+1; iil=i2; 																					// set iir to the value of i2+1 and iil to the value of i2 (as here the flow of information is from right to left so that gh^+ must be used at the cell edges)
//		iir=i1+1; iil=i1; 																					// set iir to the value of i2+1 and iil to the value of i2 (as here the flow of information is from right to left so that gh^+ must be used at the cell edges)
		if(iir==Nx)iir=0; //periodic bc																		// if iir = Nx (the maximum value that can be obtained, since i = 0,1,...,Nx-1) and so this cell is at the right boundary, requiring information from the non-existent cell with space index Nx, since there are periodic boundary conditions, set iir = 0 and use the cell with space index 0 (i.e. the cell at the left boundary)
		kkr=i1*Nx*size_v + iir*size_v + j_mod; 																// calculate the value of kkr for this value of iir
//		kkr=iir*Nx*size_v + i2*size_v + j_mod; 																// calculate the value of kkr for this value of iir
		kkl=k;																								// set kkl to k (since iil = i2)
		ur = -U_vals0[kkr*7+2]; 																			// set ur to the negative of the coefficient of the basis function with shape l which is non-zero in the cell with global index kkr (which corresponds to the evaluation of gh^+ at the right boundary and -ve since v_1 < 0 in here)
		ul = -U_vals0[kkl*7+2];																				// set ul to the negative of the coefficient of the basis function with shape l which is non-zero in the cell with global index kkl	(which corresponds to the evaluation of gh^+ at the left boundary and -ve since phi < 0 here)
	}
	else																									// do this if j1 >= Nv/2 (so that the velocity in the v1 direction is non-negative)
	{
		iir=i2; iil=i2-1;																					// set iir to the value of i2 and iil to the value of i2-1 (as here the flow of information is from left to right so that gh^- must be used at the cell edges)
//		iir=i1; iil=i1-1;																					// set iir to the value of i2 and iil to the value of i2-1 (as here the flow of information is from left to right so that gh^- must be used at the cell edges)
		if(iil==-1)iil=Nx-1; // periodic bc																	// if iil = -1 (the minimum value that can be obtained, since i = 0,1,...,Nx-1) and so this cell is at the left boundary, requiring information from the non-existent cell with space index -1, since there are periodic boundary conditions, set iil = Nx-1 and use the cell with space index Nx-1 (i.e. the cell at the right boundary)
		kkr=k; 																								// set kkr to k (since iir = i1)
		kkl=i1*Nx*size_v + iil*size_v + j_mod; 																// calculate the value of kkl for this value of iil
//		kkl=iil*Nx*size_v + i2*size_v + j_mod; 																// calculate the value of kkl for this value of iil
		ur = U_vals0[kkr*7+2];																				// set ur to the value of the coefficient of the basis function with shape l which is non-zero in the cell with global index kkr (which corresponds to the evaluation of gh^- at the right boundary and +ve since v_1 >= 0 in here)
		ul = U_vals0[kkl*7+2];																				// set ul to the value of the coefficient of the basis function with shape l which is non-zero in the cell with global index kkl (which corresponds to the evalutaion of gh^- at the left boundary and +ve since v_r >= 0 in here)
	}

	/* NEED TO MULTIPLY ALL THESE BY dx AFTER I FUNCTIONS IMPLEMENTED! */
	if(l==0)
	{
		result = scalev*(((U_vals0[kkr*7+0] - U_vals0[kkl*7+0]) + 0.5*(ur - ul) + 0.25*(U_vals0[kkr*7+6] - U_vals0[kkl*7+6]))*v_j2 + (U_vals0[kkr*7+4] - U_vals0[kkl*7+4])*dv/12.);					// calculate \int_j v1*gh*phi dv at interface x=x_i+1/2 - \int_j v1*gh*phi dv at interface x=x_i-1/2 for the basis function with shape 0 (i.e. constant) which is non-zero in the cell with global index k
	}
	else if(l==1)
	{
		result = scalev*(U_vals0[kkr*7+1] - U_vals0[kkl*7+1])*v_j2/12.;																												// calculate \int_j v1*gh*phi dv at interface x=x_i+1/2 - \int_j v1*gh*phi dv at interface x=x_i-1/2 for the basis function with shape 3 (i.e. linear in v_2) which is non-zero in the cell with global index k
	}
	else if(l==2)
	{
		result = 0.5*scalev*(((U_vals0[kkr*7+0] + U_vals0[kkl*7+0]) + 0.5*(ur + ul) + 0.25*(U_vals0[kkr*7+6] + U_vals0[kkl*7+6]))*v_j2 + (U_vals0[kkr*7+4] + U_vals0[kkl*7+4])*dv/12.);					// calculate \int_j v1*gh*phi dv at interface x=x_i+1/2 - \int_j v1*gh*phi dv at interface x=x_i-1/2 for the basis function with shape 0 (i.e. constant) which is non-zero in the cell with global index k
	}
	else if(l==3)
	{
		result = scalev*(U_vals0[kkr*7+3] - U_vals0[kkl*7+3])*v_j2/12.;																												// calculate \int_j v1*gh*phi dv at interface x=x_i+1/2 - \int_j v1*gh*phi dv at interface x=x_i-1/2 for the basis function with shape 3 (i.e. linear in v_2) which is non-zero in the cell with global index k
	}
	else if(l==4)
	{
		result = scalev*(((U_vals0[kkr*7+0] - U_vals0[kkl*7+0]) + 0.5*(ur - ul) + (U_vals0[kkr*7+6] - U_vals0[kkl*7+6])*19./60.)*dv + (U_vals0[kkr*7+4] - U_vals0[kkl*7+4])*v_j2)/12.;				// calculate \int_j v1*gh*phi dv at interface x=x_i+1/2 - \int_j v1*gh*phi dv at interface x=x_i-1/2 for the basis function with shape 2 (i.e. linear in v_1) which is non-zero in the cell with global index k
	}
	else if(l==5)
	{
		result = scalev*(U_vals0[kkr*7+5] - U_vals0[kkl*7+5])*v_j2/12.;																												// calculate \int_j v1*gh*phi dv at interface x=x_i+1/2 - \int_j v1*gh*phi dv at interface x=x_i-1/2 for the basis function with shape 4 (i.e. linear in v_3) which is non-zero in the cell with global index k
	}
	else if(l==6)
	{
		result = scalev*(((U_vals0[kkr*7+0] - U_vals0[kkl*7+0]) + 0.5*(ur - ul) + (U_vals0[kkr*7+6] - U_vals0[kkl*7+6])*19./60.)*v_j2 + (U_vals0[kkr*7+4] - U_vals0[kkl*7+4])*dv*19./180.)/4.;	// calculate \int_j v1*gh*phi dv at interface x=x_i+1/2 - \int_j v1*gh*phi dv at interface x=x_i-1/2 for the basis function with shape 5 (i.e. modulus of v) which is non-zero in the cell with global index k
	}

	return result;
}

double I5_1D(vector<double>& U_vals0, int k, int l) 	// Calculate the difference of the fifth and sixth integrals in H_(i,j), namely \int_i E*f*phi dx at interface v1==v_j+1/2 - \int_i E*f*phi dx at interface v1==v_j-1/2
{
	double result, ur, ul;																			// declare result (the result of the integral to be returned), ur (used in the evaluation of gh^+/- on the right velocity cell edge) & ul (used in the evaluation of gh^+/- on the left velocity cell edge)
	int i1, i2, j1, j2, j3, j1r, j1l, kkr, kkl; 													// declare i1, i2 (the space cell coordinates), j1, j2, j3 (the coordinates of the velocity cell), j1l (the cell from which the flux is flowing on the left of the cell with coordinate j1 in the v1 direction), j1r (the cell from which the flux is flowing on the right of the cell with coordinate j1 in the v1 direction), kkl (the global index of the cell with coordinate (i, j1l, j2, j3)) & kkr (the global index of the cell with coordinate (i, j1r, j2, j3))
	int j_mod = k%size_v;																			// declare and calculate j_mod (the remainder when k is divided by size_v = Nv^3 - used to help determine the values of i1, i2, j1, j2 & j3 from the value of k)
	int i_mod = (k-j_mod)/size_v;																	// declare and calculate i_mod (the remainder value when k has j_mod subtracted and divided through by size_v - used to help determine the values of i1 & i2)
	j3 = j_mod%Nv;																					// calculate j3 for the given k
	j2 = ((j_mod-j3)%(Nv*Nv))/Nv;																	// calculate j2 for the given k
	j1 = (j_mod-j3-j2*Nv)/(Nv*Nv);																	// calculate j1 for the given k
	i2 = i_mod%Nx;																					// calculate i2 for the given k
	i1 = (i_mod - i2)/Nx;																			// calculate i1 for the given k
  
	if(IE_X[i_mod]>0)	//2D E																					// do this if the average direction of the field E over the space cell i is positive
	//if(intE[i1] > 0)	//1D E
	{
		j1r=j1+1;  j1l=j1;																			// set j1r to the value of j1+1 and j1l to the value of j1 (as here the the average flow of the field is from left to right so that gh^- must be used at the cell edges, as information flows against the field)
		kkr=i_mod*size_v + (j1r*Nv*Nv + j2*Nv + j3);													// calculate the value of kkr for this value of j1r
		kkl=k; 																						// set kkl to k (since j1l = j1)
		if(j1r<Nv)ur = -U_vals0[kkr*7+3];																	// if j1r is not Nv (so that this cell is not receiving information from the right boundary), set ur to the negative of the coefficient of the basis function with shape 2 which is non-zero in the cell with global index kkr (which corresponds to the evaluation of gh^- at the right boundary and -ve since the field is negative here?) - note that if the cell was receiving information from the right boundary then gh^- = 0 here so ur is not needed
		ul = -U_vals0[kkl*7+3];																			// set ul to the negative of the coefficient of the basis function with shape 2 which is non-zero in the cell with global index kkl (which corresponds to the evaluation of gh^- at the right boundary and -ve since phi < 0 here)
	}
	else																							// do this if the average direction of the field E over the space cell i is non-positive
	{
		j1r=j1; j1l=j1-1;																			// set j1r to the value of j1 and j1l to the value of j1-1 (as here the the average flow of the field is from right to left so that gh^+ must be used at the cell edges, as information flows against the field)
		kkr=k;																						// set kkr to k (since j1r = j1)
		kkl=i_mod*size_v + (j1l*Nv*Nv + j2*Nv + j3);													// calculate the value of kkl for this value of j1l
		ur = U_vals0[kkr*7+3];																			// set ur to the the coefficient of the basis function with shape 2 which is non-zero in the cell with global index kkr (which corresponds to the evaluation of gh^+ at the left boundary and +ve since phi > 0 here)
		if(j1l>-1)ul = U_vals0[kkl*7+3];																	// if j1l is not -1 (so that this cell is not receiving information from the left boundary), set ul to the coefficient of the basis function with shape 2 which is non-zero in the cell with global index kkl (which corresponds to the evaluation of gh^+ at the left boundary and +ve since phi < 0 here and being subtracted?) - note that if the cell was receiving information from the left boundary then gh^+ = 0 here so ul is not needed
	}

	/* 1D E
	if(l==0)																						// calculate \int_i E*f*phi dx at interface v1==v_j+1/2 - \int_i E*f*phi dx at interface v1==v_j-1/2 for the basis function with shape 0 (i.e. constant) which is non-zero in the cell with global index k
	{
		if(j1r<Nv && j1l>-1) result = dv*dv*(U_vals0[kkr*7+0] + 0.5*ur + U_vals0[kkr*7+6]*5./12.- U_vals0[kkl*7+0] - 0.5*ul - U_vals0[kkl*7+6]*5./12.)*intE[i1] + dv*dv*(U_vals0[kkr*7+1]-U_vals0[kkl*7+1])*intE1[i1];		// this is the value at an interior cell
		else if(j1r<Nv)result =   dv*dv*(U_vals0[kkr*7+0] + 0.5*ur + U_vals0[kkr*7+6]*5./12.)*intE[i1] + dv*dv*U_vals0[kkr*7+1]*intE1[i1];																	// this is the value at the cell at the left boundary in v1 (so that the integral over the left edge is zero)
		else if(j1l>-1)result = dv*dv*(- U_vals0[kkl*7+0] - 0.5*ul - U_vals0[kkl*7+6]*5./12.)*intE[i1] - dv*dv*U_vals0[kkl*7+1]*intE1[i1];																	// this is the value at the cell at the right boundary in v1 (so that the integral over the right edge is zero)
	}
	if(l==1)																						// calculate \int_i E*f*phi dx at interface v1==v_j+1/2 - \int_i E*f*phi dx at interface v1==v_j-1/2 for the basis function with shape 1 (i.e. linear in x) which is non-zero in the cell with global index k
	{
		if(j1r<Nv && j1l>-1) result = dv*dv*( (U_vals0[kkr*7+0] + 0.5*ur + U_vals0[kkr*7+6]*5./12. - U_vals0[kkl*7+0] - 0.5*ul - U_vals0[kkl*7+6]*5./12.)*intE1[i1] + (U_vals0[kkr*7+1] - U_vals0[kkl*7+1])*intE2[i1] );		// this is the value at an interior cell
		else if(j1r<Nv)result=dv*dv*( (U_vals0[kkr*7+0] + 0.5*ur + U_vals0[kkr*7+6]*5./12.)*intE1[i1] + U_vals0[kkr*7+1]*intE2[i1] );																		// this is the value at the cell at the left boundary in v1 (so that the integral over the left edge is zero)
		else if(j1l>-1)result=dv*dv*( (- U_vals0[kkl*7+0] - 0.5*ul - U_vals0[kkl*7+6]*5./12.)*intE1[i1] - U_vals0[kkl*7+1]*intE2[i1] );																		// this is the value at the cell at the right boundary in v1 (so that the integral over the right edge is zero)
	}
	if(l==3)																						// calculate \int_i E*f*phi dx at interface v1==v_j+1/2 - \int_i E*f*phi dx at interface v1==v_j-1/2 for the basis function with shape 2 (i.e. linear in v1) which is non-zero in the cell with global index k
	{
		if(j1r<Nv && j1l>-1)result = 0.5*(dv*dv*(U_vals0[kkr*7+0] + 0.5*ur + U_vals0[kkr*7+6]*5./12.+ U_vals0[kkl*7+0] + 0.5*ul + U_vals0[kkl*7+6]*5./12.)*intE[i1] + dv*dv*(U_vals0[kkr*7+1]+U_vals0[kkl*7+1])*intE1[i1]);	// this is the value at an interior cell
		else if(j1r<Nv)result = 0.5*(dv*dv*(U_vals0[kkr*7+0] + 0.5*ur + U_vals0[kkr*7+6]*5./12.)*intE[i1] + dv*dv*U_vals0[kkr*7+1]*intE1[i1]);																// this is the value at the cell at the left boundary in v1 (so that the integral over the left edge is zero)
		else if(j1l>-1)result = 0.5*(dv*dv*(U_vals0[kkl*7+0] + 0.5*ul + U_vals0[kkl*7+6]*5./12.)*intE[i1] + dv*dv*U_vals0[kkl*7+1]*intE1[i1]);    															// this is the value at the cell at the right boundary in v1 (so that the integral over the right edge is zero)
	}
	if(l==4)																						// calculate \int_i E*f*phi dx at interface v1==v_j+1/2 - \int_i E*f*phi dx at interface v1==v_j-1/2 for the basis function with shape 3 (i.e. linear in v2) which is non-zero in the cell with global index k
	{
		if(j1r<Nv && j1l>-1)result = (U_vals0[kkr*7+4]-U_vals0[kkl*7+4])*intE[i1]*dv*dv/12.;						// this is the value at an interior cell
		else if(j1r<Nv)result = U_vals0[kkr*7+4]*intE[i1]*dv*dv/12.;										// this is the value at the cell at the left boundary in v1 (so that the integral over the left edge is zero)
		else if(j1l>-1)result = -U_vals0[kkl*7+4]*intE[i1]*dv*dv/12.;										// this is the value at the cell at the right boundary in v1 (so that the integral over the right edge is zero)
	}
	if(l==5)																						// calculate \int_i E*f*phi dx at interface v1==v_j+1/2 - \int_i E*f*phi dx at interface v1==v_j-1/2 for the basis function with shape 4 (i.e. linear in v3) which is non-zero in the cell with global index k
	{
		if(j1r<Nv && j1l>-1)result = (U_vals0[kkr*7+5]-U_vals0[kkl*7+5])*intE[i1]*dv*dv/12.;						// this is the value at an interior cell
		else if(j1r<Nv)result=U_vals0[kkr*7+5]*intE[i1]*dv*dv/12.;											// this is the value at the cell at the left boundary in v1 (so that the integral over the left edge is zero)
		else if(j1l>-1)result=-U_vals0[kkl*7+5]*intE[i1]*dv*dv/12.;										// this is the value at the cell at the right boundary in v1 (so that the integral over the right edge is zero)
	}
	if(l==6)																						// calculate \int_i E*f*phi dx at interface v1==v_j+1/2 - \int_i E*f*phi dx at interface v1==v_j-1/2 for the basis function with shape 5 (i.e. the modulus of v) which is non-zero in the cell with global index k
	{
		if(j1r<Nv && j1l>-1)result = dv*dv*( ((U_vals0[kkr*7+0] + 0.5*ur - U_vals0[kkl*7+0] - 0.5*ul)*5./12. + (U_vals0[kkr*7+6]- U_vals0[kkl*7+6])*133./720.)*intE[i1] + (U_vals0[kkr*7+1] - U_vals0[kkl*7+1])*intE1[i1]*5./12. ); //BUG: coefficient of U[k][5] was 11/48 insteadof 133/720		// this is the value at an interior cell
		else if(j1r<Nv)result= dv*dv*( ((U_vals0[kkr*7+0] + 0.5*ur)*5./12. + U_vals0[kkr*7+6]*133./720.)*intE[i1] + U_vals0[kkr*7+1]*intE1[i1]*5./12. );													// this is the value at the cell at the left boundary in v1 (so that the integral over the left edge is zero)
		else if(j1l>-1)result=-dv*dv*( ((U_vals0[kkl*7+0] + 0.5*ul)*5./12. + U_vals0[kkl*7+6]*133./720.)*intE[i1] + U_vals0[kkl*7+1]*intE1[i1]*5./12. );													// this is the value at the cell at the right boundary in v1 (so that the integral over the right edge is zero)
	}
	*/

	/* 2D E */
	if(l==0)																						// calculate \int_i E*f*phi dx at interface v1==v_j+1/2 - \int_i E*f*phi dx at interface v1==v_j-1/2 for the basis function with shape 0 (i.e. constant) which is non-zero in the cell with global index k
	{
		if(j1r<Nv && j1l>-1) result = dv*dv*(U_vals0[kkr*7+0] + 0.5*ur + U_vals0[kkr*7+6]*5./12.- U_vals0[kkl*7+0] - 0.5*ul - U_vals0[kkl*7+6]*5./12.)*IE_X[i_mod] + dv*dv*(U_vals0[kkr*7+1]-U_vals0[kkl*7+1])*IXE_X[i_mod];		// this is the value at an interior cell
		else if(j1r<Nv)result =   dv*dv*(U_vals0[kkr*7+0] + 0.5*ur + U_vals0[kkr*7+6]*5./12.)*IE_X[i_mod] + dv*dv*U_vals0[kkr*7+1]*IXE_X[i_mod];																	// this is the value at the cell at the left boundary in v1 (so that the integral over the left edge is zero)
		else if(j1l>-1)result = dv*dv*(- U_vals0[kkl*7+0] - 0.5*ul - U_vals0[kkl*7+6]*5./12.)*IE_X[i_mod] - dv*dv*U_vals0[kkl*7+1]*IXE_X[i_mod];																	// this is the value at the cell at the right boundary in v1 (so that the integral over the right edge is zero)
	}
	if(l==1)																						// calculate \int_i E*f*phi dx at interface v1==v_j+1/2 - \int_i E*f*phi dx at interface v1==v_j-1/2 for the basis function with shape 1 (i.e. linear in x) which is non-zero in the cell with global index k
	{
		if(j1r<Nv && j1l>-1) result = dv*dv*( (U_vals0[kkr*7+0] + 0.5*ur + U_vals0[kkr*7+6]*5./12. - U_vals0[kkl*7+0] - 0.5*ul - U_vals0[kkl*7+6]*5./12.)*IXE_X[i_mod] + (U_vals0[kkr*7+1] - U_vals0[kkl*7+1])*IXXE_X[i_mod] );		// this is the value at an interior cell
		else if(j1r<Nv)result=dv*dv*( (U_vals0[kkr*7+0] + 0.5*ur + U_vals0[kkr*7+6]*5./12.)*IXE_X[i_mod] + U_vals0[kkr*7+1]*IXXE_X[i_mod] );																		// this is the value at the cell at the left boundary in v1 (so that the integral over the left edge is zero)
		else if(j1l>-1)result=dv*dv*( (- U_vals0[kkl*7+0] - 0.5*ul - U_vals0[kkl*7+6]*5./12.)*IXE_X[i_mod] - U_vals0[kkl*7+1]*IXXE_X[i_mod] );																		// this is the value at the cell at the right boundary in v1 (so that the integral over the right edge is zero)
	}
	if(l==3)																						// calculate \int_i E*f*phi dx at interface v1==v_j+1/2 - \int_i E*f*phi dx at interface v1==v_j-1/2 for the basis function with shape 2 (i.e. linear in v1) which is non-zero in the cell with global index k
	{
		if(j1r<Nv && j1l>-1)result = 0.5*(dv*dv*(U_vals0[kkr*7+0] + 0.5*ur + U_vals0[kkr*7+6]*5./12.+ U_vals0[kkl*7+0] + 0.5*ul + U_vals0[kkl*7+6]*5./12.)*IE_X[i_mod] + dv*dv*(U_vals0[kkr*7+1]+U_vals0[kkl*7+1])*IXE_X[i_mod]);	// this is the value at an interior cell
		else if(j1r<Nv)result = 0.5*(dv*dv*(U_vals0[kkr*7+0] + 0.5*ur + U_vals0[kkr*7+6]*5./12.)*IE_X[i_mod] + dv*dv*U_vals0[kkr*7+1]*IXE_X[i_mod]);																// this is the value at the cell at the left boundary in v1 (so that the integral over the left edge is zero)
		else if(j1l>-1)result = 0.5*(dv*dv*(U_vals0[kkl*7+0] + 0.5*ul + U_vals0[kkl*7+6]*5./12.)*IE_X[i_mod] + dv*dv*U_vals0[kkl*7+1]*IXE_X[i_mod]);    															// this is the value at the cell at the right boundary in v1 (so that the integral over the right edge is zero)
	}
	if(l==4)																						// calculate \int_i E*f*phi dx at interface v1==v_j+1/2 - \int_i E*f*phi dx at interface v1==v_j-1/2 for the basis function with shape 3 (i.e. linear in v2) which is non-zero in the cell with global index k
	{
		if(j1r<Nv && j1l>-1)result = (U_vals0[kkr*7+4]-U_vals0[kkl*7+4])*IE_X[i_mod]*dv*dv/12.;						// this is the value at an interior cell
		else if(j1r<Nv)result = U_vals0[kkr*7+4]*IE_X[i_mod]*dv*dv/12.;										// this is the value at the cell at the left boundary in v1 (so that the integral over the left edge is zero)
		else if(j1l>-1)result = -U_vals0[kkl*7+4]*IE_X[i_mod]*dv*dv/12.;										// this is the value at the cell at the right boundary in v1 (so that the integral over the right edge is zero)
	}
	if(l==5)																						// calculate \int_i E*f*phi dx at interface v1==v_j+1/2 - \int_i E*f*phi dx at interface v1==v_j-1/2 for the basis function with shape 4 (i.e. linear in v3) which is non-zero in the cell with global index k
	{
		if(j1r<Nv && j1l>-1)result = (U_vals0[kkr*7+5]-U_vals0[kkl*7+5])*IE_X[i_mod]*dv*dv/12.;						// this is the value at an interior cell
		else if(j1r<Nv)result=U_vals0[kkr*7+5]*IE_X[i_mod]*dv*dv/12.;											// this is the value at the cell at the left boundary in v1 (so that the integral over the left edge is zero)
		else if(j1l>-1)result=-U_vals0[kkl*7+5]*IE_X[i_mod]*dv*dv/12.;										// this is the value at the cell at the right boundary in v1 (so that the integral over the right edge is zero)
	}
	if(l==6)																						// calculate \int_i E*f*phi dx at interface v1==v_j+1/2 - \int_i E*f*phi dx at interface v1==v_j-1/2 for the basis function with shape 5 (i.e. the modulus of v) which is non-zero in the cell with global index k
	{
		if(j1r<Nv && j1l>-1)result = dv*dv*( ((U_vals0[kkr*7+0] + 0.5*ur - U_vals0[kkl*7+0] - 0.5*ul)*5./12. + (U_vals0[kkr*7+6]- U_vals0[kkl*7+6])*133./720.)*IE_X[i_mod] + (U_vals0[kkr*7+1] - U_vals0[kkl*7+1])*IXE_X[i_mod]*5./12. ); //BUG: coefficient of U[k][5] was 11/48 insteadof 133/720		// this is the value at an interior cell
		else if(j1r<Nv)result= dv*dv*( ((U_vals0[kkr*7+0] + 0.5*ur)*5./12. + U_vals0[kkr*7+6]*133./720.)*IE_X[i_mod] + U_vals0[kkr*7+1]*IXE_X[i_mod]*5./12. );													// this is the value at the cell at the left boundary in v1 (so that the integral over the left edge is zero)
		else if(j1l>-1)result=-dv*dv*( ((U_vals0[kkl*7+0] + 0.5*ul)*5./12. + U_vals0[kkl*7+6]*133./720.)*IE_X[i_mod] + U_vals0[kkl*7+1]*IXE_X[i_mod]*5./12. );													// this is the value at the cell at the right boundary in v1 (so that the integral over the right edge is zero)
	}
	/**/
  	return result;
}

double I5_v1(vector<double>& U_vals0, int k, int l) 	// Calculate the difference of the fifth and sixth integrals in H_(i,j), namely \int_i E*f*phi dx at interface v1==v_j+1/2 - \int_i E*f*phi dx at interface v1==v_j-1/2
{
	double result, ur, ul;																			// declare result (the result of the integral to be returned), ur (used in the evaluation of gh^+/- on the right velocity cell edge) & ul (used in the evaluation of gh^+/- on the left velocity cell edge)
	int i1, i2, j1, j2, j3, j1r, j1l, kkr, kkl; 													// declare i1, i2 (the space cell coordinates), j1, j2, j3 (the coordinates of the velocity cell), j1l (the cell from which the flux is flowing on the left of the cell with coordinate j1 in the v1 direction), j1r (the cell from which the flux is flowing on the right of the cell with coordinate j1 in the v1 direction), kkl (the global index of the cell with coordinate (i, j1l, j2, j3)) & kkr (the global index of the cell with coordinate (i, j1r, j2, j3))
	int j_mod = k%size_v;																			// declare and calculate j_mod (the remainder when k is divided by size_v = Nv^3 - used to help determine the values of i1, i2, j1, j2 & j3 from the value of k)
	int i_mod = (k-j_mod)/size_v;																	// declare and calculate i_mod (the remainder value when k has j_mod subtracted and divided through by size_v - used to help determine the values of i1 & i2)
	j3 = j_mod%Nv;																					// calculate j3 for the given k
	j2 = ((j_mod-j3)%(Nv*Nv))/Nv;																	// calculate j2 for the given k
	j1 = (j_mod-j3-j2*Nv)/(Nv*Nv);																	// calculate j1 for the given k

	if(IE_X[i_mod]>0)	//2D E																					// do this if the average direction of the field E over the space cell i is positive
	{
		j1r=j1+1;  j1l=j1;																			// set j1r to the value of j1+1 and j1l to the value of j1 (as here the the average flow of the field is from left to right so that gh^- must be used at the cell edges, as information flows against the field)
		kkr=i_mod*size_v + (j1r*Nv*Nv + j2*Nv + j3);													// calculate the value of kkr for this value of j1r
		kkl=k; 																						// set kkl to k (since j1l = j1)
		if(j1r<Nv)ur = -U_vals0[kkr*7+3];																	// if j1r is not Nv (so that this cell is not receiving information from the right boundary), set ur to the negative of the coefficient of the basis function with shape 2 which is non-zero in the cell with global index kkr (which corresponds to the evaluation of gh^- at the right boundary and -ve since the field is negative here?) - note that if the cell was receiving information from the right boundary then gh^- = 0 here so ur is not needed
		ul = -U_vals0[kkl*7+3];																			// set ul to the negative of the coefficient of the basis function with shape 2 which is non-zero in the cell with global index kkl (which corresponds to the evaluation of gh^- at the right boundary and -ve since phi < 0 here)
	}
	else																							// do this if the average direction of the field E over the space cell i is non-positive
	{
		j1r=j1; j1l=j1-1;																			// set j1r to the value of j1 and j1l to the value of j1-1 (as here the the average flow of the field is from right to left so that gh^+ must be used at the cell edges, as information flows against the field)
		kkr=k;																						// set kkr to k (since j1r = j1)
		kkl=i_mod*size_v + (j1l*Nv*Nv + j2*Nv + j3);													// calculate the value of kkl for this value of j1l
		ur = U_vals0[kkr*7+3];																			// set ur to the the coefficient of the basis function with shape 2 which is non-zero in the cell with global index kkr (which corresponds to the evaluation of gh^+ at the left boundary and +ve since phi > 0 here)
		if(j1l>-1)ul = U_vals0[kkl*7+3];																	// if j1l is not -1 (so that this cell is not receiving information from the left boundary), set ul to the coefficient of the basis function with shape 2 which is non-zero in the cell with global index kkl (which corresponds to the evaluation of gh^+ at the left boundary and +ve since phi < 0 here and being subtracted?) - note that if the cell was receiving information from the left boundary then gh^+ = 0 here so ul is not needed
	}

	if(l==0)																						// calculate \int_i E*f*phi dx at interface v1==v_j+1/2 - \int_i E*f*phi dx at interface v1==v_j-1/2 for the basis function with shape 0 (i.e. constant) which is non-zero in the cell with global index k
	{
		if(j1r<Nv && j1l>-1)
		{
			result = dv*dv*(((U_vals0[kkr*7+0] - U_vals0[kkl*7+0]) + 0.5*(ur - ul) + (U_vals0[kkr*7+6] - U_vals0[kkl*7+6])*5./12.)*IE_X[i_mod] + (U_vals0[kkr*7+1] - U_vals0[kkl*7+1])*IXE_X[i_mod] + (U_vals0[kkr*7+2] - U_vals0[kkl*7+2])*IYE_X[i_mod]);		// this is the value at an interior cell
		}
		else if(j1r<Nv)
		{
			result = dv*dv*((U_vals0[kkr*7+0] + 0.5*ur + U_vals0[kkr*7+6]*5./12.)*IE_X[i_mod] + U_vals0[kkr*7+1]*IXE_X[i_mod] + U_vals0[kkr*7+2]*IYE_X[i_mod]);		// this is the value at the cell at the left boundary in v1 (so that the integral over the left edge is zero)
		}
		else if(j1l>-1)
		{
			result = dv*dv*((-U_vals0[kkl*7+0] - 0.5*ul - U_vals0[kkl*7+6]*5./12.)*IE_X[i_mod] - U_vals0[kkl*7+1]*IXE_X[i_mod] - U_vals0[kkl*7+2]*IYE_X[i_mod]);	// this is the value at the cell at the right boundary in v1 (so that the integral over the right edge is zero)
		}
	}
	if(l==1)																						// calculate \int_i E*f*phi dx at interface v1==v_j+1/2 - \int_i E*f*phi dx at interface v1==v_j-1/2 for the basis function with shape 1 (i.e. linear in x1) which is non-zero in the cell with global index k
	{
		if(j1r<Nv && j1l>-1)
		{
			result = dv*dv*(((U_vals0[kkr*7+0] - U_vals0[kkl*7+0]) + 0.5*(ur - ul) + (U_vals0[kkr*7+6] - U_vals0[kkl*7+6])*5./12.)*IXE_X[i_mod] + (U_vals0[kkr*7+1] - U_vals0[kkl*7+1])*IXXE_X[i_mod] + (U_vals0[kkr*7+2] - U_vals0[kkl*7+2])*IXYE_X[i_mod]);		// this is the value at an interior cell
		}
		else if(j1r<Nv)
		{
			result = dv*dv*((U_vals0[kkr*7+0] + 0.5*ur + U_vals0[kkr*7+6]*5./12.)*IXE_X[i_mod] + U_vals0[kkr*7+1]*IXXE_X[i_mod] + U_vals0[kkr*7+2]*IXYE_X[i_mod]);		// this is the value at the cell at the left boundary in v1 (so that the integral over the left edge is zero)
		}
		else if(j1l>-1)
		{
			result = dv*dv*((-U_vals0[kkl*7+0] - 0.5*ul - U_vals0[kkl*7+6]*5./12.)*IXE_X[i_mod] - U_vals0[kkl*7+1]*IXXE_X[i_mod] - U_vals0[kkl*7+2]*IXYE_X[i_mod]);	// this is the value at the cell at the right boundary in v1 (so that the integral over the right edge is zero)
		}
	}
	if(l==2)																						// calculate \int_i E*f*phi dx at interface v1==v_j+1/2 - \int_i E*f*phi dx at interface v1==v_j-1/2 for the basis function with shape 2 (i.e. linear in x2) which is non-zero in the cell with global index k
	{
		if(j1r<Nv && j1l>-1)
		{
			result = dv*dv*(((U_vals0[kkr*7+0] - U_vals0[kkl*7+0]) + 0.5*(ur - ul) + (U_vals0[kkr*7+6] - U_vals0[kkl*7+6])*5./12.)*IYE_X[i_mod] + (U_vals0[kkr*7+1] - U_vals0[kkl*7+1])*IXYE_X[i_mod] + (U_vals0[kkr*7+2] - U_vals0[kkl*7+2])*IYYE_X[i_mod]);		// this is the value at an interior cell
		}
		else if(j1r<Nv)
		{
			result = dv*dv*((U_vals0[kkr*7+0] + 0.5*ur + U_vals0[kkr*7+6]*5./12.)*IYE_X[i_mod] + U_vals0[kkr*7+1]*IXYE_X[i_mod] + U_vals0[kkr*7+2]*IYYE_X[i_mod]);		// this is the value at the cell at the left boundary in v1 (so that the integral over the left edge is zero)
		}
		else if(j1l>-1)
		{
			result = dv*dv*((-U_vals0[kkl*7+0] - 0.5*ul - U_vals0[kkl*7+6]*5./12.)*IYE_X[i_mod] - U_vals0[kkl*7+1]*IXYE_X[i_mod] - U_vals0[kkl*7+2]*IYYE_X[i_mod]);	// this is the value at the cell at the right boundary in v1 (so that the integral over the right edge is zero)
		}
	}
	if(l==3)																						// calculate \int_i E*f*phi dx at interface v1==v_j+1/2 - \int_i E*f*phi dx at interface v1==v_j-1/2 for the basis function with shape 2 (i.e. linear in v1) which is non-zero in the cell with global index k
	{
		if(j1r<Nv && j1l>-1)
		{
			result = 0.5*dv*dv*(((U_vals0[kkr*7+0] + U_vals0[kkl*7+0]) + 0.5*(ur + ul) + (U_vals0[kkr*7+6] + U_vals0[kkl*7+6])*5./12.)*IE_X[i_mod] + (U_vals0[kkr*7+1] + U_vals0[kkl*7+1])*IXE_X[i_mod] + (U_vals0[kkr*7+2] + U_vals0[kkl*7+2])*IYE_X[i_mod]);		// this is the value at an interior cell
		}
		else if(j1r<Nv)
		{
			result = 0.5*dv*dv*((U_vals0[kkr*7+0] + 0.5*ur + U_vals0[kkr*7+6]*5./12.)*IE_X[i_mod] + U_vals0[kkr*7+1]*IXE_X[i_mod] + U_vals0[kkr*7+2]*IYE_X[i_mod]);		// this is the value at the cell at the left boundary in v1 (so that the integral over the left edge is zero)
		}
		else if(j1l>-1)
		{
			result = 0.5*dv*dv*((U_vals0[kkl*7+0] + 0.5*ul + U_vals0[kkl*7+6]*5./12.)*IE_X[i_mod] + U_vals0[kkl*7+1]*IXE_X[i_mod] + U_vals0[kkl*7+2]*IYE_X[i_mod]);	// this is the value at the cell at the right boundary in v1 (so that the integral over the right edge is zero)
		}
	}
	if(l==4)																						// calculate \int_i E*f*phi dx at interface v1==v_j+1/2 - \int_i E*f*phi dx at interface v1==v_j-1/2 for the basis function with shape 3 (i.e. linear in v2) which is non-zero in the cell with global index k
	{
		if(j1r<Nv && j1l>-1)
		{
			result = dv*dv*(U_vals0[kkr*7+4] - U_vals0[kkl*7+4])*IE_X[i_mod]/12.;						// this is the value at an interior cell
		}
		else if(j1r<Nv)
		{
			result = dv*dv*U_vals0[kkr*7+4]*IE_X[i_mod]/12.;										// this is the value at the cell at the left boundary in v1 (so that the integral over the left edge is zero)
		}
		else if(j1l>-1)
		{
			result = -dv*dv*U_vals0[kkl*7+4]*IE_X[i_mod]/12.;										// this is the value at the cell at the right boundary in v1 (so that the integral over the right edge is zero)
		}
	}
	if(l==5)																						// calculate \int_i E*f*phi dx at interface v1==v_j+1/2 - \int_i E*f*phi dx at interface v1==v_j-1/2 for the basis function with shape 4 (i.e. linear in v3) which is non-zero in the cell with global index k
	{
		if(j1r<Nv && j1l>-1)
		{
			result = dv*dv*(U_vals0[kkr*7+5] - U_vals0[kkl*7+5])*IE_X[i_mod]/12.;						// this is the value at an interior cell
		}
		else if(j1r<Nv)
		{
			result = dv*dv*U_vals0[kkr*7+5]*IE_X[i_mod]/12.;										// this is the value at the cell at the left boundary in v1 (so that the integral over the left edge is zero)
		}
		else if(j1l>-1)
		{
			result = -dv*dv*U_vals0[kkl*7+5]*IE_X[i_mod]/12.;										// this is the value at the cell at the right boundary in v1 (so that the integral over the right edge is zero)
		}
	}
	if(l==6)																						// calculate \int_i E*f*phi dx at interface v1==v_j+1/2 - \int_i E*f*phi dx at interface v1==v_j-1/2 for the basis function with shape 5 (i.e. the modulus of v) which is non-zero in the cell with global index k
	{
		if(j1r<Nv && j1l>-1)
		{
			result = dv*dv*(((U_vals0[kkr*7+0] - U_vals0[kkl*7+0])*5./12.+ (ur - ul)*5./24. + (U_vals0[kkr*7+6] - U_vals0[kkl*7+6])*133./720.)*IE_X[i_mod] + (U_vals0[kkr*7+1] - U_vals0[kkl*7+1])*IXE_X[i_mod]*5./12. + (U_vals0[kkr*7+2] - U_vals0[kkl*7+2])*IYE_X[i_mod]*5./12.);		// this is the value at an interior cell
		}
		else if(j1r<Nv)
		{
			result = dv*dv*((U_vals0[kkr*7+0]*5./12. + ur*5./24. + U_vals0[kkr*7+6]*133./720.)*IE_X[i_mod] + U_vals0[kkr*7+1]*IXE_X[i_mod]*5./12. + U_vals0[kkr*7+2]*IYE_X[i_mod]*5./12.);		// this is the value at the cell at the left boundary in v1 (so that the integral over the left edge is zero)
		}
		else if(j1l>-1)
		{
			result = dv*dv*((-U_vals0[kkl*7+0]*5./12. - ul*5./24. - U_vals0[kkl*7+6]*133./720.)*IE_X[i_mod] - U_vals0[kkl*7+1]*IXE_X[i_mod]*5./12. - U_vals0[kkl*7+2]*IYE_X[i_mod]*5./12.);	// this is the value at the cell at the right boundary in v1 (so that the integral over the right edge is zero)
		}
	}
  	return result;
}

double I5_v2(vector<double>& U_vals0, int k, int l) 	// Calculate the difference of the fifth and sixth integrals in H_(i,j), namely \int_i E*f*phi dx at interface v1==v_j+1/2 - \int_i E*f*phi dx at interface v1==v_j-1/2
{
	double result, ur, ul;																			// declare result (the result of the integral to be returned), ur (used in the evaluation of gh^+/- on the right velocity cell edge) & ul (used in the evaluation of gh^+/- on the left velocity cell edge)
	int i1, i2, j1, j2, j3, j2r, j2l, kkr, kkl; 													// declare i1, i2 (the space cell coordinates), j1, j2, j3 (the coordinates of the velocity cell), j2l (the cell from which the flux is flowing on the left of the cell with coordinate j1 in the v1 direction), j2r (the cell from which the flux is flowing on the right of the cell with coordinate j1 in the v1 direction), kkl (the global index of the cell with coordinate (i, j2l, j2, j3)) & kkr (the global index of the cell with coordinate (i, j2r, j2, j3))
	int j_mod = k%size_v;																			// declare and calculate j_mod (the remainder when k is divided by size_v = Nv^3 - used to help determine the values of i1, i2, j1, j2 & j3 from the value of k)
	int i_mod = (k-j_mod)/size_v;																	// declare and calculate i_mod (the remainder value when k has j_mod subtracted and divided through by size_v - used to help determine the values of i1 & i2)
	j3 = j_mod%Nv;																					// calculate j3 for the given k
	j2 = ((j_mod-j3)%(Nv*Nv))/Nv;																	// calculate j2 for the given k
	j1 = (j_mod-j3-j2*Nv)/(Nv*Nv);																	// calculate j1 for the given k

	if(IE_Y[i_mod]>0)	//2D E																					// do this if the average direction of the field E over the space cell i is positive
	{
		j2r=j2+1;  j2l=j2;																			// set j2r to the value of j1+1 and j2l to the value of j1 (as here the the average flow of the field is from left to right so that gh^- must be used at the cell edges, as information flows against the field)
		kkr=i_mod*size_v + (j1*Nv*Nv + j2r*Nv + j3);													// calculate the value of kkr for this value of j2r
		kkl=k; 																						// set kkl to k (since j2l = j1)
		if(j2r<Nv)ur = -U_vals0[kkr*7+4];																	// if j2r is not Nv (so that this cell is not receiving information from the right boundary), set ur to the negative of the coefficient of the basis function with shape 2 which is non-zero in the cell with global index kkr (which corresponds to the evaluation of gh^- at the right boundary and -ve since the field is negative here?) - note that if the cell was receiving information from the right boundary then gh^- = 0 here so ur is not needed
		ul = -U_vals0[kkl*7+4];																			// set ul to the negative of the coefficient of the basis function with shape 2 which is non-zero in the cell with global index kkl (which corresponds to the evaluation of gh^- at the right boundary and -ve since phi < 0 here)
	}
	else																							// do this if the average direction of the field E over the space cell i is non-positive
	{
		j2r=j2; j2l=j2-1;																			// set j2r to the value of j1 and j2l to the value of j1-1 (as here the the average flow of the field is from right to left so that gh^+ must be used at the cell edges, as information flows against the field)
		kkr=k;																						// set kkr to k (since j2r = j1)
		kkl=i_mod*size_v + (j1*Nv*Nv + j2l*Nv + j3);													// calculate the value of kkl for this value of j2l
		ur = U_vals0[kkr*7+4];																			// set ur to the the coefficient of the basis function with shape 2 which is non-zero in the cell with global index kkr (which corresponds to the evaluation of gh^+ at the left boundary and +ve since phi > 0 here)
		if(j2l>-1)ul = U_vals0[kkl*7+4];																	// if j2l is not -1 (so that this cell is not receiving information from the left boundary), set ul to the coefficient of the basis function with shape 2 which is non-zero in the cell with global index kkl (which corresponds to the evaluation of gh^+ at the left boundary and +ve since phi < 0 here and being subtracted?) - note that if the cell was receiving information from the left boundary then gh^+ = 0 here so ul is not needed
	}

	if(l==0)																						// calculate \int_i E*f*phi dx at interface v1==v_j+1/2 - \int_i E*f*phi dx at interface v1==v_j-1/2 for the basis function with shape 0 (i.e. constant) which is non-zero in the cell with global index k
	{
		if(j2r<Nv && j2l>-1)
		{
			result = dv*dv*(((U_vals0[kkr*7+0] - U_vals0[kkl*7+0]) + 0.5*(ur - ul) + (U_vals0[kkr*7+6] - U_vals0[kkl*7+6])*5./12.)*IE_Y[i_mod] + (U_vals0[kkr*7+1] - U_vals0[kkl*7+1])*IXE_Y[i_mod] + (U_vals0[kkr*7+2] - U_vals0[kkl*7+2])*IYE_Y[i_mod]);		// this is the value at an interior cell
		}
		else if(j2r<Nv)
		{
			result = dv*dv*((U_vals0[kkr*7+0] + 0.5*ur + U_vals0[kkr*7+6]*5./12.)*IE_Y[i_mod] + U_vals0[kkr*7+1]*IXE_Y[i_mod] + U_vals0[kkr*7+2]*IYE_Y[i_mod]);		// this is the value at the cell at the left boundary in v1 (so that the integral over the left edge is zero)
		}
		else if(j2l>-1)
		{
			result = dv*dv*((-U_vals0[kkl*7+0] - 0.5*ul - U_vals0[kkl*7+6]*5./12.)*IE_Y[i_mod] - U_vals0[kkl*7+1]*IXE_Y[i_mod] - U_vals0[kkl*7+2]*IYE_Y[i_mod]);	// this is the value at the cell at the right boundary in v1 (so that the integral over the right edge is zero)
		}
	}
	if(l==1)																						// calculate \int_i E*f*phi dx at interface v1==v_j+1/2 - \int_i E*f*phi dx at interface v1==v_j-1/2 for the basis function with shape 1 (i.e. linear in x1) which is non-zero in the cell with global index k
	{
		if(j2r<Nv && j2l>-1)
		{
			result = dv*dv*(((U_vals0[kkr*7+0] - U_vals0[kkl*7+0]) + 0.5*(ur - ul) + (U_vals0[kkr*7+6] - U_vals0[kkl*7+6])*5./12.)*IXE_Y[i_mod] + (U_vals0[kkr*7+1] - U_vals0[kkl*7+1])*IXXE_Y[i_mod] + (U_vals0[kkr*7+2] - U_vals0[kkl*7+2])*IXYE_Y[i_mod]);		// this is the value at an interior cell
		}
		else if(j2r<Nv)
		{
			result = dv*dv*((U_vals0[kkr*7+0] + 0.5*ur + U_vals0[kkr*7+6]*5./12.)*IXE_Y[i_mod] + U_vals0[kkr*7+1]*IXXE_Y[i_mod] + U_vals0[kkr*7+2]*IXYE_Y[i_mod]);		// this is the value at the cell at the left boundary in v1 (so that the integral over the left edge is zero)
		}
		else if(j2l>-1)
		{
			result = dv*dv*((-U_vals0[kkl*7+0] - 0.5*ul - U_vals0[kkl*7+6]*5./12.)*IXE_Y[i_mod] - U_vals0[kkl*7+1]*IXXE_Y[i_mod] - U_vals0[kkl*7+2]*IXYE_Y[i_mod]);	// this is the value at the cell at the right boundary in v1 (so that the integral over the right edge is zero)
		}
	}
	if(l==2)																						// calculate \int_i E*f*phi dx at interface v1==v_j+1/2 - \int_i E*f*phi dx at interface v1==v_j-1/2 for the basis function with shape 2 (i.e. linear in x2) which is non-zero in the cell with global index k
	{
		if(j2r<Nv && j2l>-1)
		{
			result = dv*dv*(((U_vals0[kkr*7+0] - U_vals0[kkl*7+0]) + 0.5*(ur - ul) + (U_vals0[kkr*7+6] - U_vals0[kkl*7+6])*5./12.)*IYE_Y[i_mod] + (U_vals0[kkr*7+1] - U_vals0[kkl*7+1])*IXYE_Y[i_mod] + (U_vals0[kkr*7+2] - U_vals0[kkl*7+2])*IYYE_Y[i_mod]);		// this is the value at an interior cell
		}
		else if(j2r<Nv)
		{
			result = dv*dv*((U_vals0[kkr*7+0] + 0.5*ur + U_vals0[kkr*7+6]*5./12.)*IYE_Y[i_mod] + U_vals0[kkr*7+1]*IXYE_Y[i_mod] + U_vals0[kkr*7+2]*IYYE_Y[i_mod]);		// this is the value at the cell at the left boundary in v1 (so that the integral over the left edge is zero)
		}
		else if(j2l>-1)
		{
			result = dv*dv*((-U_vals0[kkl*7+0] - 0.5*ul - U_vals0[kkl*7+6]*5./12.)*IYE_Y[i_mod] - U_vals0[kkl*7+1]*IXYE_Y[i_mod] - U_vals0[kkl*7+2]*IYYE_Y[i_mod]);	// this is the value at the cell at the right boundary in v1 (so that the integral over the right edge is zero)
		}
	}
	if(l==3)																						// calculate \int_i E*f*phi dx at interface v1==v_j+1/2 - \int_i E*f*phi dx at interface v1==v_j-1/2 for the basis function with shape 2 (i.e. linear in v1) which is non-zero in the cell with global index k
	{
		if(j2r<Nv && j2l>-1)
		{
			result = dv*dv*(U_vals0[kkr*7+3] - U_vals0[kkl*7+3])*IE_Y[i_mod]/12.;						// this is the value at an interior cell
		}
		else if(j2r<Nv)
		{
			result = dv*dv*U_vals0[kkr*7+3]*IE_Y[i_mod]/12.;										// this is the value at the cell at the left boundary in v1 (so that the integral over the left edge is zero)
		}
		else if(j2l>-1)
		{
			result = -dv*dv*U_vals0[kkl*7+3]*IE_Y[i_mod]/12.;										// this is the value at the cell at the right boundary in v1 (so that the integral over the right edge is zero)
		}
	}
	if(l==4)																						// calculate \int_i E*f*phi dx at interface v1==v_j+1/2 - \int_i E*f*phi dx at interface v1==v_j-1/2 for the basis function with shape 3 (i.e. linear in v2) which is non-zero in the cell with global index k
	{
		if(j2r<Nv && j2l>-1)
		{
			result = 0.5*dv*dv*(((U_vals0[kkr*7+0] + U_vals0[kkl*7+0]) + 0.5*(ur + ul) + (U_vals0[kkr*7+6] + U_vals0[kkl*7+6])*5./12.)*IE_Y[i_mod] + (U_vals0[kkr*7+1] + U_vals0[kkl*7+1])*IXE_Y[i_mod] + (U_vals0[kkr*7+2] + U_vals0[kkl*7+2])*IYE_Y[i_mod]);		// this is the value at an interior cell
		}
		else if(j2r<Nv)
		{
			result = 0.5*dv*dv*((U_vals0[kkr*7+0] + 0.5*ur + U_vals0[kkr*7+6]*5./12.)*IE_Y[i_mod] + U_vals0[kkr*7+1]*IXE_Y[i_mod] + U_vals0[kkr*7+2]*IYE_Y[i_mod]);		// this is the value at the cell at the left boundary in v1 (so that the integral over the left edge is zero)
		}
		else if(j2l>-1)
		{
			result = 0.5*dv*dv*((U_vals0[kkl*7+0] + 0.5*ul + U_vals0[kkl*7+6]*5./12.)*IE_Y[i_mod] + U_vals0[kkl*7+1]*IXE_Y[i_mod] + U_vals0[kkl*7+2]*IYE_Y[i_mod]);	// this is the value at the cell at the right boundary in v1 (so that the integral over the right edge is zero)
		}
	}
	if(l==5)																						// calculate \int_i E*f*phi dx at interface v1==v_j+1/2 - \int_i E*f*phi dx at interface v1==v_j-1/2 for the basis function with shape 4 (i.e. linear in v3) which is non-zero in the cell with global index k
	{
		if(j2r<Nv && j2l>-1)
		{
			result = dv*dv*(U_vals0[kkr*7+5] - U_vals0[kkl*7+5])*IE_Y[i_mod]/12.;						// this is the value at an interior cell
		}
		else if(j2r<Nv)
		{
			result = dv*dv*U_vals0[kkr*7+5]*IE_Y[i_mod]/12.;										// this is the value at the cell at the left boundary in v1 (so that the integral over the left edge is zero)
		}
		else if(j2l>-1)
		{
			result = -dv*dv*U_vals0[kkl*7+5]*IE_Y[i_mod]/12.;										// this is the value at the cell at the right boundary in v1 (so that the integral over the right edge is zero)
		}
	}
	if(l==6)																						// calculate \int_i E*f*phi dx at interface v1==v_j+1/2 - \int_i E*f*phi dx at interface v1==v_j-1/2 for the basis function with shape 5 (i.e. the modulus of v) which is non-zero in the cell with global index k
	{
		if(j2r<Nv && j2l>-1)
		{
			result = dv*dv*(((U_vals0[kkr*7+0] - U_vals0[kkl*7+0])*5./12.+ (ur - ul)*5./24. + (U_vals0[kkr*7+6] - U_vals0[kkl*7+6])*133./720.)*IE_Y[i_mod] + (U_vals0[kkr*7+1] - U_vals0[kkl*7+1])*IXE_Y[i_mod]*5./12. + (U_vals0[kkr*7+2] - U_vals0[kkl*7+2])*IYE_Y[i_mod]*5./12.);		// this is the value at an interior cell
		}
		else if(j2r<Nv)
		{
			result = dv*dv*((U_vals0[kkr*7+0]*5./12. + ur*5./24. + U_vals0[kkr*7+6]*133./720.)*IE_Y[i_mod] + U_vals0[kkr*7+1]*IXE_Y[i_mod]*5./12. + U_vals0[kkr*7+2]*IYE_Y[i_mod]*5./12.);		// this is the value at the cell at the left boundary in v1 (so that the integral over the left edge is zero)
		}
		else if(j2l>-1)
		{
			result = dv*dv*((-U_vals0[kkl*7+0]*5./12. - ul*5./24. - U_vals0[kkl*7+6]*133./720.)*IE_Y[i_mod] - U_vals0[kkl*7+1]*IXE_Y[i_mod]*5./12. - U_vals0[kkl*7+2]*IYE_Y[i_mod]*5./12.);	// this is the value at the cell at the right boundary in v1 (so that the integral over the right edge is zero)
		}
	}
  	return result;
}

#ifdef UseMPI
/*
void computeH(double *H, double *U)// H_k(i,j)(f, E, phi_l)  					// MAY NEVER BE CALLED
{
  int k, l; // k=i*Nv^3 + (j1*Nv*Nv + j2*Nv + j3)
  double tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tp0, tp5;
 
  #pragma omp parallel for schedule(dynamic) private(tp0, tp5, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6) shared(H, U)
  for(k=chunksize_dg*myrank_mpi;k<chunksize_dg*(myrank_mpi+1);k++){        
	tp0 = I1(U,k,0) - I2(U,k,0) - I3(U,k,0) + I5(U,k,0);
	tp5 = I1(U,k,5) - I2(U,k,5) - I3(U,k,5) + I5(U,k,5);
	H[(k%chunksize_dg)*7+0] = (19*tp0/4. - 15*tp5)/(dx*dv*dv*dv);
	H[(k%chunksize_dg)*7+6] = (60*tp5 - 15*tp0)/(dx*dv*dv*dv);
    for(l=1;l<6;l++){
	  tmp1=I1(U,k,l); tmp2=I2(U,k,l); tmp3=I3(U,k,l);  tmp5=I5(U,k,l); 
	  H[(k%chunksize_dg)*7+l] = (tmp1-tmp2-tmp3+tmp5)*12./(dx*dv*dv*dv);
	}	
  }  
}
*/

void RK3(vector<double>& U_vals, vector<double>& POTC, vector<double>& phix, vector<double>& phiy) // RK3 for f_t = H(f)
{
  int i, k, l, k_local;
  double tp0, tp1, tp2, tp3, tp4, tp5, H[7];//, tp0, tp5, tmp1, tmp2, tmp3, tmp5;
 
  MPI_Status status;
  
  /* 1D
  ce = computePhi_x_0(U_vals);

  #pragma omp parallel for private(i) shared(U_vals,cp, intE, intE1)
  for(i=0;i<Nx;i++){
    cp[i] = computeC_rho(U_vals,i); intE[i] = Int_E(U_vals,i); intE1[i] = Int_E1st(U_vals,i);
  }
  
  #pragma omp parallel for private(i) shared(intE2)
  for(i=0;i<Nx;i++){
    intE2[i] = Int_E2nd(U_vals,i); // BUG: Int_E2nd() require knowldege of cp
  }
  */
  /*
  if(myrank_mpi == 0)
  {
	  for(k=0; k<size; k++)
	  {
		  for(l=0; l<7; l++)
		  {
			  printf("U_vals[%d*7+%d] = %g \n", k, l, U_vals[k*7+l]);
		  }
	  }
  }
  */
  pois2d(U_vals, POTC, phix, phiy);																// solve Poisson's equation, using the DG coefficients stored in U, storing the coefficients of the potential, field in the x1 direction & field in the x2 direction in POTC, phix & phiy, respectively

  /* 1D
  for(i=0;i<Nx;i++)
  {
	  printf("intE[%d] = %g, IE_X[%d] = %g, dx*IE_X[%d] = %g \n", i, intE[i], i, IE_X[i*Nx], i, dx*IE_X[i*Nx]);
	  printf("intE1[%d] = %g, IXE_X[%d] = %g, dx*IXE_X[%d] = %g \n", i, intE1[i], i, IXE_X[i*Nx], i, dx*IXE_X[i*Nx]);
	  printf("intE2[%d] = %g, IXXE_X[%d] = %g, dx*IXXE_X[%d] = %g \n\n", i, intE2[i], i, IXXE_X[i*Nx], i, dx*IXXE_X[i*Nx]);
  }
  */

  #pragma omp parallel for schedule(dynamic)  private(H,k, k_local, l, tp0, tp1, tp2, tp3, tp4, tp5) shared(U_vals, Utmp)
  for(k=chunksize_dg*myrank_mpi;k<chunksize_dg*(myrank_mpi+1);k++){
    k_local = k%chunksize_dg;
    
    tp0=I1(U_vals,k,0)-I2(U_vals,k,0)-I3_x1(U_vals,k,0)-I3_x2(U_vals,k,0)+I5_v1(U_vals,k,0)+I5_v2(U_vals,k,0);
    tp1=I1(U_vals,k,1)-I2(U_vals,k,1)-I3_x1(U_vals,k,1)-I3_x2(U_vals,k,1)+I5_v1(U_vals,k,1)+I5_v2(U_vals,k,1);
    tp2=I1(U_vals,k,3)-I2(U_vals,k,3)-I3_x1(U_vals,k,3)-I3_x2(U_vals,k,3)+I5_v1(U_vals,k,3)+I5_v2(U_vals,k,3);
    tp3=I1(U_vals,k,4)-I2(U_vals,k,4)-I3_x1(U_vals,k,4)-I3_x2(U_vals,k,4)+I5_v1(U_vals,k,4)+I5_v2(U_vals,k,4);
    tp4=I1(U_vals,k,5)-I2(U_vals,k,5)-I3_x1(U_vals,k,5)-I3_x2(U_vals,k,5)+I5_v1(U_vals,k,5)+I5_v2(U_vals,k,5);
    tp5=I1(U_vals,k,6)-I2(U_vals,k,6)-I3_x1(U_vals,k,6)-I3_x2(U_vals,k,6)+I5_v1(U_vals,k,6)+I5_v2(U_vals,k,6);

    //H[k_local][0] = (19*tp[0]/4. - 15*tp[5])/dx/scalev;
    //H[k_local][5] = (60*tp[5] - 15*tp[0])/dx/scalev;	
    H[0] = (19*tp0/4. - 15*tp5)/dx/scalev;
    H[6] = (60*tp5 - 15*tp0)/dx/scalev;
    //for(l=1;l<5;l++)H[l] = tp[l]*12./dx/scalev;;//H[k_local][l] = tp[l]*12./dx/scalev;
    H[1] = tp1*12./dx/scalev; H[3] = tp2*12./dx/scalev; H[4] = tp3*12./dx/scalev; H[5] = tp4*12./dx/scalev;
    H[2] = 0; // termporarily necessary as stray uninitiated values may have been messing with U results

    for(l=0;l<7;l++) Utmp[k_local*7+l] = U_vals[k*7+l] + dt*H[l];
  } 

  if(myrank_mpi == 0) {
    //dump the weights we've computed into U1
    for(k=0;k<chunksize_dg;k++) {
	for(l=0;l<7;l++) U1[k*7+l] = Utmp[k*7+l];
    } 
    //receive from all other processes
    for(i=1;i<nprocs_mpi;i++) {
	    MPI_Recv(output_buffer_vp, chunksize_dg*7, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &status); // receive weight from other processes consecutively, rank i=1..numNodes-1, ensuring the weights are stored in the file consecutively !
	    for(k=0;k<chunksize_dg;k++) {
			for(l=0;l<7;l++)U1[(k + i*chunksize_dg)*7+l] = output_buffer_vp[k*7+l];
		}			
    }
  }
  else  MPI_Send(Utmp, chunksize_dg*7, MPI_DOUBLE, 0, myrank_mpi, MPI_COMM_WORLD);
  
  
  MPI_Bcast(&U1[0], size*7, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  /*
  if(myrank_mpi == 0)
  {
	  for(k=0; k<size; k++)
	  {
		  for(l=0; l<7; l++)
		  {
			  printf("U1[%d*7+%d] = %g \n", k, l, U1[k*7+l]);
		  }
	  }
  }
  */
  /////////////////// 1st step of RK3 done//////////////////////////////////////////////////////// 
    
  /* 1D
  ce = computePhi_x_0(U1); 

  #pragma omp parallel for private(i) shared(U1,cp, intE, intE1)
  for(i=0;i<Nx;i++){
    cp[i] = computeC_rho(U1,i); intE[i] = Int_E(U1,i); intE1[i] = Int_E1st(U1,i); 
  }
  
  #pragma omp parallel for private(i) shared(U1,intE2)
  for(i=0;i<Nx;i++){
    intE2[i] = Int_E2nd(U1,i); // BUG: Int_E2nd() require knowldege of cp 
  }
  */

  pois2d(U1, POTC, phix, phiy);																// solve Poisson's equation, using the DG coefficients stored in U1, storing the coefficients of the potential, field in the x1 direction & field in the x2 direction in POTC, phix & phiy, respectively

  #pragma omp parallel for schedule(dynamic) private(H, k, k_local, l, tp0, tp1, tp2, tp3, tp4, tp5)  shared(U1,Utmp)
  for(k=chunksize_dg*myrank_mpi;k<chunksize_dg*(myrank_mpi+1);k++){      
    k_local = k%chunksize_dg;
    
    tp0=I1(U1,k,0)-I2(U1,k,0)-I3_x1(U1,k,0)-I3_x2(U1,k,0)+I5_v1(U1,k,0)+I5_v2(U1,k,0);
    tp1=I1(U1,k,1)-I2(U1,k,1)-I3_x1(U1,k,1)-I3_x2(U1,k,1)+I5_v1(U1,k,1)+I5_v2(U1,k,1);
    tp2=I1(U1,k,3)-I2(U1,k,3)-I3_x1(U1,k,3)-I3_x2(U1,k,3)+I5_v1(U1,k,3)+I5_v2(U1,k,3);
    tp3=I1(U1,k,4)-I2(U1,k,4)-I3_x1(U1,k,4)-I3_x2(U1,k,4)+I5_v1(U1,k,4)+I5_v2(U1,k,4);
    tp4=I1(U1,k,5)-I2(U1,k,5)-I3_x1(U1,k,5)-I3_x2(U1,k,5)+I5_v1(U1,k,5)+I5_v2(U1,k,5);
    tp5=I1(U1,k,6)-I2(U1,k,6)-I3_x1(U1,k,6)-I3_x2(U1,k,6)+I5_v1(U1,k,6)+I5_v2(U1,k,6);

    //H[k_local][0] = (19*tp[0]/4. - 15*tp[5])/dx/scalev;
    //H[k_local][5] = (60*tp[5] - 15*tp[0])/dx/scalev;	
    H[0] = (19*tp0/4. - 15*tp5)/dx/scalev;
    H[6] = (60*tp5 - 15*tp0)/dx/scalev;
    //for(l=1;l<5;l++)H[l] = tp[l]*12./dx/scalev;;//H[k_local][l] = tp[l]*12./dx/scalev;
    H[1] = tp1*12./dx/scalev; H[3] = tp2*12./dx/scalev; H[4] = tp3*12./dx/scalev; H[5] = tp4*12./dx/scalev;
    H[2] = 0; // termporarily necessary as stray uninitiated values may have been messing with U results

    for(l=0;l<7;l++) Utmp[k_local*7+l] = 0.75*U_vals[k*7+l] + 0.25*U1[k*7+l] + 0.25*dt*H[l];
  }    
  if(myrank_mpi == 0) {
    //dump the weights we've computed into U1
    for(k=0;k<chunksize_dg;k++) {
	      for(l=0;l<7;l++) U1[k*7+l] = Utmp[k*7+l];
    } 
    //receive from all other processes
    for(i=1;i<nprocs_mpi;i++) {      
		MPI_Recv(output_buffer_vp, chunksize_dg*7, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &status); // receive weight from other processes consecutively, rank i=1..numNodes-1, ensuring the weights are stored in the file consecutively !
		for(k=0;k<chunksize_dg;k++) {
			for(l=0;l<7;l++)U1[(k + i*chunksize_dg)*7+l] = output_buffer_vp[k*7+l];
        }
    }
  }
  else MPI_Send(Utmp, chunksize_dg*7, MPI_DOUBLE, 0, myrank_mpi, MPI_COMM_WORLD);
  
  
  MPI_Bcast(&U1[0], size*7, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  /////////////////// 2nd step of RK3 done//////////////////////////////////////////////////////// 
   
  /* 1D
  ce = computePhi_x_0(U1);

  #pragma omp parallel for private(i) shared(U1,cp, intE, intE1)
  for(i=0;i<Nx;i++){
    cp[i] = computeC_rho(U1,i); intE[i] = Int_E(U1,i); intE1[i] = Int_E1st(U1,i); 
  }
  
  #pragma omp parallel for private(i) shared(U1,intE2)
  for(i=0;i<Nx;i++){
    intE2[i] = Int_E2nd(U1,i); // BUG: Int_E2nd() require knowldege of cp 
  }
  */

  pois2d(U1, POTC, phix, phiy);																// solve Poisson's equation, using the DG coefficients stored in U, storing the coefficients of the potential, field in the x1 direction & field in the x2 direction in POTC, phix & phiy, respectively

  #pragma omp parallel for schedule(dynamic) private(H, k, k_local, l, tp0, tp1, tp2, tp3, tp4, tp5)  shared(U1,Utmp)
  for(k=chunksize_dg*myrank_mpi;k<chunksize_dg*(myrank_mpi+1);k++){      
    k_local = k%chunksize_dg;
  
    tp0=I1(U1,k,0)-I2(U1,k,0)-I3_x1(U1,k,0)-I3_x2(U1,k,0)+I5_v1(U1,k,0)+I5_v2(U1,k,0);
    tp1=I1(U1,k,1)-I2(U1,k,1)-I3_x1(U1,k,1)-I3_x2(U1,k,1)+I5_v1(U1,k,1)+I5_v2(U1,k,1);
    tp2=I1(U1,k,3)-I2(U1,k,3)-I3_x1(U1,k,3)-I3_x2(U1,k,3)+I5_v1(U1,k,3)+I5_v2(U1,k,3);
    tp3=I1(U1,k,4)-I2(U1,k,4)-I3_x1(U1,k,4)-I3_x2(U1,k,4)+I5_v1(U1,k,4)+I5_v2(U1,k,4);
    tp4=I1(U1,k,5)-I2(U1,k,5)-I3_x1(U1,k,5)-I3_x2(U1,k,5)+I5_v1(U1,k,5)+I5_v2(U1,k,5);
    tp5=I1(U1,k,6)-I2(U1,k,6)-I3_x1(U1,k,6)-I3_x2(U1,k,6)+I5_v1(U1,k,6)+I5_v2(U1,k,6);

    H[0] = (19*tp0/4. - 15*tp5)/dx/scalev;
    H[6] = (60*tp5 - 15*tp0)/dx/scalev;
    //for(l=1;l<5;l++)H[l] = tp[l]*12./dx/scalev;;//H[k_local][l] = tp[l]*12./dx/scalev;
    H[1] = tp1*12./dx/scalev; H[3] = tp2*12./dx/scalev; H[4] = tp3*12./dx/scalev; H[5] = tp4*12./dx/scalev;
    H[2] = 0; // termporarily necessary as stray uninitiated values may have been messing with U results

    for(l=0;l<7;l++) Utmp[k_local*7+l] = U_vals[k*7+l]/3. + U1[k*7+l]*2./3. + dt*H[l]*2./3.;
  }    
  if(myrank_mpi == 0) {
    //dump the weights we've computed into U1
    for(k=0;k<chunksize_dg;k++) {
	      for(l=0;l<7;l++) U_vals[k*7+l] = Utmp[k*7+l];
    } 
    //receive from all other processes
    for(i=1;i<nprocs_mpi;i++) {      
	    MPI_Recv(output_buffer_vp, chunksize_dg*7, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &status); // receive weight from other processes consecutively, rank i=1..numNodes-1, ensuring the weights are stored in the file consecutively !
	    for(k=0;k<chunksize_dg;k++) {
		  for(l=0;l<7;l++)U_vals[(k + i*chunksize_dg)*7+l] = output_buffer_vp[k*7+l];
        } 
    }
  }
  else MPI_Send(Utmp, chunksize_dg*7, MPI_DOUBLE, 0, myrank_mpi, MPI_COMM_WORLD);
  
  MPI_Bcast(&U_vals[0], size*7, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD); 
  /////////////////// 3rd step of RK3 done//////////////////////////////////////////////////////// 
}


#else
/*
void computeH(double *U)// H_k(i,j)(f, E, phi_l)  
{
  int i, k, l; // k=i*Nv^3 + (j1*Nv*Nv + j2*Nv + j3)
  double tp0, tp1, tp2, tp3, tp4, tp5;
 

  ce = computePhi_x_0(U); 
  // #pragma omp barrier
  #pragma omp parallel for private(i) shared(U,cp, intE, intE1)
  for(i=0;i<Nx;i++){
    cp[i] = computeC_rho(U,i); intE[i] = Int_E(U,i); intE1[i] = Int_E1st(U,i); 
  }
  
  #pragma omp parallel for private(i) shared(U, intE2)
  for(i=0;i<Nx;i++){
    intE2[i] = Int_E2nd(U,i); // BUG: Int_E2nd() require knowldege of cp 
  }
  
  #pragma omp parallel for schedule(dynamic) private(tp0, tp1, tp2, tp3, tp4, tp5, k,l) shared(U, H)
  for(k=0;k<size;k++){
   // double tp[6];

  // #pragma omp critical
   // {
      /*for(l=0;l<6;l++){
	tp[l]=I1(U,k,l)-I2(U,k,l)-I3(U,k,l)+I5(U,k,l);		
      } *//*
   // }
    tp0=I1(U,k,0)-I2(U,k,0)-I3(U,k,0)+I5(U,k,0);
    tp1=I1(U,k,1)-I2(U,k,1)-I3(U,k,1)+I5(U,k,1);
    tp2=I1(U,k,3)-I2(U,k,3)-I3(U,k,3)+I5(U,k,3);
    tp3=I1(U,k,4)-I2(U,k,4)-I3(U,k,4)+I5(U,k,4);
    tp4=I1(U,k,5)-I2(U,k,5)-I3(U,k,5)+I5(U,k,5);
    tp5=I1(U,k,6)-I2(U,k,6)-I3(U,k,6)+I5(U,k,6);
     
    H[k*7+0] = (19*tp0/4. - 15*tp5)/dx/scalev;
    H[k*7+6] = (60*tp5 - 15*tp0)/dx/scalev;
    //for(l=1;l<5;l++)H[k][l] = tp[l]*12./dx/scalev;

    H[k*7+1] = tp1*12./dx/scalev; H[k*7+3] = tp2*12./dx/scalev; H[k*7+4] = tp3*12./dx/scalev; H[k*7+5] = tp4*12./dx/scalev;
  }
  
}
*/

void RK3(vector<double>& U_vals) // RK3 for f_t = H(f)
{
  int k, l;
  //double **H = (double **)malloc(size*sizeof(double *));
  //for (l=0;l<size;l++) H[l] = (double*)malloc(6*sizeof(double));
  //double **U1 = (double **)malloc(size*sizeof(double *));
  //for (l=0;l<size;l++) U1[l] = (double*)malloc(6*sizeof(double));
  
  computeH(U_vals);
   #pragma omp parallel for private(k,l) shared(U_vals)
  for(k=0;k<size;k++){				  
	  for(l=0;l<7;l++) U1[k*7+l] = U_vals[k*7+l] + dt*H[k*7+l];
  }
 
  computeH(U1);
   #pragma omp parallel for private(k,l) shared(U_vals)
  for(k=0;k<size;k++){						  
	  for(l=0;l<7;l++) U1[k*7+l] = 0.75*U_vals[k*7+l] + 0.25*U1[k*7+l] + 0.25*dt*H[k*7+l];
  }

  computeH(U1);
   #pragma omp parallel for private(k,l) shared(U_vals)
  for(k=0;k<size;k++){						  
	for(l=0;l<7;l++) U_vals[k*7+l] = U_vals[k*7+l]/3. + U1[k*7+l]*2./3. + dt*H[k*7+l]*2./3.;
  }
  //free(H); free(U1);
}
#endif
