#include "MomentCalculations.h"																					// MomentCalculations.h is where the prototypes for the functions contained in this file are declared and any variables defined here to be used throughout the other files are declared as external

double computeMass(vector<double>& U_vals)
{
	int k;
	double tmp=0.;
	#pragma omp parallel for private(k) shared(U_vals) reduction(+:tmp)
	for(k=0;k<size;k++)
	{
		tmp += U_vals[k*7+0] + U_vals[k*7+6]/4.;
	}

	return tmp*dx*dx*scalev;
}

void computeMomentum(vector<double>& U_vals, double *a)
{
	int k, i, j,j1,j2,j3; //k=i*Nv*Nv*Nv + (j1*Nv*Nv + j2*Nv + j3);
	double tmp1=0., tmp2=0., tmp3=0.;
	a[0]=0.; a[1]=0.; a[2]=0.; // the three momentum
	#pragma omp parallel for private(k,i,j,j1,j2,j3) shared(U_vals) reduction(+:tmp1, tmp2, tmp3)  //reduction directive may change the result a little bit
	for(k=0;k<size;k++)
	{
		j=k%size_v; i=(k-j)/size_v;
		j3=j%Nv; j2=((j-j3)%(Nv*Nv))/Nv; j1=(j-j3-j2*Nv)/(Nv*Nv);
		tmp1 += Gridv((double)j1)*dv*U_vals[k*7+0] + U_vals[k*7+3]*dv*dv/12. + U_vals[k*7+6]*Gridv((double)j1)*dv/4.;
		tmp2 += Gridv((double)j2)*dv*U_vals[k*7+0] + U_vals[k*7+4]*dv*dv/12. + U_vals[k*7+6]*Gridv((double)j2)*dv/4.;
		tmp3 += Gridv((double)j3)*dv*U_vals[k*7+0] + U_vals[k*7+5]*dv*dv/12. + U_vals[k*7+6]*Gridv((double)j3)*dv/4.;
	}
	a[0]=tmp1*dx*dx*dv*dv; a[1]=tmp2*dx*dx*dv*dv; a[2]=tmp3*dx*dx*dv*dv;
}

double computeKiE(vector<double>& U_vals)
{
	int k, i, j,j1,j2,j3;
	double tmp=0., tp=0., tp1=0.;
	//#pragma omp parallel for private(k,i,j,j1,j2,j3,tp, tp1) shared(U_vals) reduction(+:tmp)
	for(k=0;k<size;k++)
	{
		j=k%size_v; i=(k-j)/size_v;
		j3=j%Nv; j2=((j-j3)%(Nv*Nv))/Nv; j1=(j-j3-j2*Nv)/(Nv*Nv);
		tp1 = Gridv((double)j1)*Gridv((double)j1) + Gridv((double)j2)*Gridv((double)j2) + Gridv((double)j3)*Gridv((double)j3);
		tmp += U_vals[k*7+0]*(tp1 + dv*dv/4.)*dv + (Gridv(j1)*U_vals[k*7+3]+Gridv(j2)*U_vals[k*7+4]+Gridv(j3)*U_vals[k*7+5])*dv*dv/6. + U_vals[k*7+6]*( dv*dv*dv*19./240. + tp1*dv/4.);
	}
	tmp *= dx*dx*dv*dv;
	return 0.5*tmp;
}

/* 1D
double computeEleE_1D(vector<double>& U_vals)
{
  int k, i, j;
  double retn, tmp1=0., tmp2=0., tmp3=0., tmp4=0., tmp5=0., tmp6=0., tmp7=0., tp1, tp2, c;
  double ce1, cp1;
  ce1 = computePhi_x_0(U_vals);

  tmp1 = ce1*ce1*Lx;
  tmp2 = Lx*Lx*Lx/3.; tmp3 = -ce1*Lx*Lx;

  //#pragma omp parallel for private(j,k, i, tp1, tp2, cp1, c) shared(U) reduction(+:tmp4, tmp5, tmp6)
  for(i=0;i<Nx;i++){
    c = Int_Int_rho(U_vals,i);
    cp1 = computeC_rho(U_vals,i);
    tmp4 += dx*cp1 + c;
    tmp5 += dx*Gridx((double)i)*cp1;
    tp1=0.; tp2=0.;
    for(j=0;j<size_v;j++){
      k = i*size_v + j;
      tp1 += (U_vals[k*7+0] + U_vals[k*7+6]/4.);
      tp2 += U_vals[k*7+1];
    }
    tmp5 += scalev* (tp1*( (pow(Gridx(i+0.5), 3) - pow(Gridx(i-0.5), 3))/3. - Gridx(i-0.5)*Gridx((double)i)*dx ) - tp2 * dx*dx*Gridx((double)i)/12.);

    tp2 *= dx/2.;
    tmp6 +=  cp1*cp1*dx + 2*cp1*c + pow(dv, 6)* ( tp1*tp1*dx*dx*dx/3. + tp2*tp2*dx/30. - tp1*tp2*dx*dx/6.);//+ tp1*tp2*(dx*Gridx((double)i)/6. - dx*dx/4.) ); //Int_Cumulativerho_sqr(i);
  }
  retn = tmp1 + tmp2 + tmp3 + 2*ce1*tmp4 - 2*tmp5 + tmp6;
  return 0.5*retn;
}
*/

double computeEleE(vector<double>& phix_vals, vector<double>& phiy_vals)							// function to compute the electric energy (EleE = 0.5*\int_{Omega_x} |E(x,t)|^2 dx) using the coefficients of the field stored in phix_vals & phiy_vals
{
	double energy = 0;																				// initialise the energy at 0
	for(int k=0; k<size_x; k++)																		// loop through each space cell
	{
		energy += phix_vals[3*k]*phix_vals[3*k] + phiy_vals[3*k]*phiy_vals[3*k];					// add the contribution from integrating the constant basis function on the current cell
		energy += (phix_vals[3*k+1]*phix_vals[3*k+1] + phiy_vals[3*k+1]*phiy_vals[3*k+1])/12;		// add the contribution from integrating the basis function linear in x1 on the current cell
		energy += (phix_vals[3*k+2]*phix_vals[3*k+2] + phiy_vals[3*k+2]*phiy_vals[3*k+2])/12;		// add the contribution from integrating the basis function linear in x2 on the current cell
	}
	energy = 0.5*scalex*energy;																		// multiply the sum by 0.5 as well as dx^2 for the volume of each cell
	return energy;																					// return the value of the energy calculated
}
