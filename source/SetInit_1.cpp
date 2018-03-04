#include "SetInit_1.h"																					// SetInit_1.h is where the prototypes for the functions contained in this file are declared and any variables defined here to be used throughout the other files are declared as external

//static double vt[4] = {-0.3399810435848562648026658,0.3399810435848562648026658,-0.8611363115940525752239465,0.8611363115940525752239465};
//static double wt[4] = {0.6521451548625461426269361,0.6521451548625461426269361,0.3478548451374538573730639,0.3478548451374538573730639};
//double wt[5]={0.5688888888888889, 0.4786286704993665, 0.4786286704993665,0.2369268850561891, 0.2369268850561891};
//double vt[5]={0., -0.5384693101056831,0.5384693101056831,-0.9061798459386640,0.9061798459386640};
	
double sinc(double x)
{
  double result;
  
  if(x==0.0)
    {
      result = 1.0;
    }
  else
    {
      result = sin(x)/x;
    }
  
  return result;
}

/*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$*/

void trapezoidalRule(int nPoints, double *weight)
{
	int i;
	
	weight[0] = 0.5;
	weight[nPoints-1] = weight[0];
	
	for(i=1;i<nPoints-1;i++) weight[i] = 1.0;
}

double f_TS(double v1, double v2, double v3) //two-stream instability initial
{
  double r2=v1*v1+v2*v2+v3*v3;
  return r2*exp(-r2/2)/(2*PI*sqrt(2*PI))/3.;
}

double f_2Gauss(double v1, double v2, double v3)
{	
  double retn, sig=M_PI/10;
  retn = 0.5*(exp(-((v1-2*sig)*(v1-2*sig)+v2*v2+v3*v3)/(2*sig*sig))+exp(-((v1+2*sig)*(v1+2*sig)+v2*v2+v3*v3)/(2*sig*sig)))/(2*M_PI*sig*sig*sqrt(2*M_PI*sig*sig));
  return retn;
}

double Mw(double v1, double v2, double v3)
{
  double r2, T, retn;
  T=0.9; // when testing the nonlinear damping, T was chosen too small that the "effective grid" is  not fine enough 
  r2=v1*v1+v2*v2+v3*v3;
  retn = exp(-r2/(2*T))/(2*PI*T*sqrt(2*T*PI));
 // if(v1<=dv && v1>=-dv)retn=0.5/(dv*Lv*Lv);
  //else retn=0.;
  return retn;
}

void SetInit_LD(vector<double>& U_vals)																			// function to calculate the DG coefficients for the initial condition for Landau Damping
{
    int i1, i2, j1, j2, j3, k, m1,m2,m3,nt=5;																	// declare i1, i2 (to represent cell (i1,i2) in x-space), j1, j2, j3 (to represent cell (j1,j2,j3) in v-space), k (the index of cell (i,j1,j2,j3) in U), m1, m2, m3 (counters for the Gaussian quadrature in 3D) & nt (the number of points in the quadrature)
    double a=A_amp, c=k_wave;																					// declare a (the amplitude of cosine wave) and set it to A_amp & c (the frequency of the cosine wave) and set it to k_wave
    double tp, tp0, tp5, tmp0, tmp1, tmp2, tmp3, tmp4;															// declare tp, tp0, tmp0, tmp1, tmp2, tmp3, tmp4 (temporary values while calculating the quadrature for the integral w.r.t. v)
    //#pragma omp parallel for private(k,j1,j2,j3,i,tmp0, tmp1, tmp2, tmp3, tmp4, tp0, tp5, tp) shared(U)
    for(j1=0;j1<Nv;j1++)																						// loop through all the velocity cells
    {
    	for(j2=0;j2<Nv;j2++)
    	{
    		for(j3=0;j3<Nv;j3++)
    		{
    			tmp0=0.; tmp1=0.; tmp2=0.; tmp3=0.; tmp4=0.;													// initialise tmp0, tmp1, tmp2, tmp3 & tmp4 at 0 for a new quadrature integral to calculate int_Kj Mw(v)*phi_(6k+l) dv, for l = 0, 2, 4, 5, 6
    			for(m1=0;m1<nt;m1++)																			// loop through the quadrature sum
    			{
    				for(m2=0;m2<nt;m2++)
    				{
    					for(m3=0;m3<nt;m3++)
    					{
							#ifdef Damping
							tp = wt[m1]*wt[m2]*wt[m3]*Mw(Gridv((double)j1)+0.5*dv*vt[m1], Gridv((double)j2)+0.5*dv*vt[m2], Gridv((double)j3)+0.5*dv*vt[m3]);		// calculate w_m1*w_m2*w_m3*Mw(v_m1,v_m2,v_m3), which appears in all quadrature integral approximations
							#endif

							#ifdef TwoStream
							tp = wt[m1]*wt[m2]*wt[m3]*f_2Gauss(Gridv((double)j1)+0.5*dv*vt[m1], Gridv((double)j2)+0.5*dv*vt[m2], Gridv((double)j3)+0.5*dv*vt[m3]);
							#endif

							tmp0 += tp;																																// add tp to tmp0 (for the integral int_Kj Mw(v)*phi_6k dv)
							tmp1 += tp*0.5*vt[m1];																													// add tp*v_m1/2 to tmp1 (for the integral int_Kj Mw(v)*phi_(6k+2) dv)
							tmp2 += tp*0.5*vt[m2];																													// add tp*v_m2/2 to tmp2 (for the integral int_Kj Mw(v)*phi_(6k+3) dv)
							tmp3 += tp*0.5*vt[m3];																													// add tp*v_m3/2 to tmp3 (for the integral int_Kj Mw(v)*phi_(6k+4) dv)
							tmp4 += tp*0.25*(vt[m1]*vt[m1] + vt[m2]*vt[m2]+ vt[m3]*vt[m3]);																			// add tp*((v_m1/2)^2 + (v_m2/2)^2 + (v_m3/2)^2) to tmp4 (for the integral int_Kj Mw(v)*phi_(6k+5) dv)
    					}
    				}
    			}
    			tmp0 = tmp0*0.5*0.5*0.5; tmp1 = tmp1*0.5*0.5*0.5; tmp2 = tmp2*0.5*0.5*0.5; tmp3 = tmp3*0.5*0.5*0.5; tmp4 = tmp4*0.5*0.5*0.5;						// multiply tmp0, tmp1, tmp2, tmp3 & tmp4 by (1/2)^3 to represent the fact that quadrature isn't done over [-1, 1] (should also multiply by dv^3 but this cancels with 1/dv^3 later)
    			for(i1=0;i1<Nx;i1++)																																// loop through the space cells in the x1 direction|
    			{
        			for(i2=0;i2<Nx;i2++)		// NOTE THAT EVERY i2 HERE SHOULD BE i1 BUT NOT WHILE STILL 1D														// loop through the space cells in the x2 direction
        			{
        				k=i1*Nx*size_v +i2*size_v + (j1*Nv*Nv + j2*Nv + j3);																										// calculate the index of cell (i1,i2,j1,j2,j3) in U
        				tp0 = (dx + (sin(c*Gridx(i1+0.5)) - sin(c*Gridx(i1-0.5)))*a/c)*tmp0/dx;																		// calculate b_7k = (int_Ii (1 + Acos(kx)) dx)*(int_Kj Mw(v)*phi_7k(v) dv) (NOTE: int_(Omega_i) (1 + a*cos(c*x)) dx = dx + (sin(c*x_(i+0.5)) - sin(c*x_(i+0.5)))*(a/c))
        				tp5 = (dx + (sin(c*Gridx(i1+0.5)) - sin(c*Gridx(i1-0.5)))*a/c)*tmp4/dx;																		// calculate b_(7k+6) = (int_Ii (1 + Acos(kx)) dx)*(int_Kj Mw(v)*phi_(7k+6)(v) dv) (NOTE: int_(Omega_i) (1 + a*cos(c*x)) dx = dx + (sin(c*x_(i+0.5)) - sin(c*x_(i+0.5)))*(a/c))
        				U_vals[7*k+0] = 19*tp0/4. - 15*tp5;																												// calculate the coefficient U[7k]
        				U_vals[7*k+6] = 60*tp5 - 15*tp0;																													// calculate the coefficient U[7k+6]

        				U_vals[7*k+1] = (0.5*(sin(c*Gridx(i1+0.5)) + sin(c*Gridx(i1-0.5))) + (cos(c*Gridx(i1+0.5)) - cos(c*Gridx(i1-0.5)))/(c*dx))*(a/c)*tmp0*12./dx;	// calculate the coefficient U[7k+1] = 12*b_(7k+1) = 12*(int_Ii (1 + Acos(kx))*phi_(7k+1)(x) dx)*(int_Kj Mw(v) dv) (NOTE: int_(Omega_i) (1 + a*cos(c*x))*phi_(6k+1)(x) dx = (0.5*(sin(c*x_(i+0.5)) + sin(c*x_(i+0.5))) + (cos(c*x_(i+0.5)) - cos(c*x_(i+0.5)))/(c*dx))*(a/c))
        				U_vals[7*k+3] = (dx + (sin(c*Gridx(i1+0.5)) - sin(c*Gridx(i1-0.5)))*a/c)*tmp1*12/dx;																// calculate the coefficient U[7k+3] = 12*b_(7k+3) = 12*(int_Ii (1 + Acos(kx)) dx)*(int_Kj Mw(v)*phi_(6k+2)(v) dv) (NOTE: int_(Omega_i) (1 + a*cos(c*x)) dx = dx + (sin(c*x_(i+0.5)) - sin(c*x_(i+0.5)))*(a/c))
        				U_vals[7*k+4] = (dx + (sin(c*Gridx(i1+0.5)) - sin(c*Gridx(i1-0.5)))*a/c)*tmp2*12/dx;																// calculate the coefficient U[7k+4] = 12*b_(7k+4) = 12*(int_Ii (1 + Acos(kx)) dx)*(int_Kj Mw(v)*phi_(6k+3)(v) dv) (NOTE: int_(Omega_i) (1 + a*cos(c*x)) dx = dx + (sin(c*x_(i+0.5)) - sin(c*x_(i+0.5)))*(a/c))
        				U_vals[7*k+5] = (dx + (sin(c*Gridx(i1+0.5)) - sin(c*Gridx(i1-0.5)))*a/c)*tmp3*12/dx;																// calculate the coefficient U[7k+5] = 12*b_(7k+5) = 12*(int_Ii (1 + Acos(kx)) dx)*(int_Kj Mw(v)*phi_(6k+4)(v) dv) (NOTE: int_(Omega_i) (1 + a*cos(c*x)) dx = dx + (sin(c*x_(i+0.5)) - sin(c*x_(i+0.5)))*(a/c))
        			}
    			}
    		}
    	}
	}
}




#ifdef MPI
void setInit_spectral(double *U, double **f)
{
  int i, j1, j2, j3, k, l, m ,n;
  for(i=chunk_Nx*myrank_mpi;i<chunk_Nx*(myrank_mpi+1) && i<Nx;i++){  
    for(l=0;l<N;l++){
      j1 = (l*h_v)/dv; // integer part = floor() for non-negative integers.
      if(j1==Nv)j1=Nv-1; // let the right end point lie in the last element
      for(m=0;m<N;m++){
		j2 = (m*h_v)/dv;
		if(j2==Nv)j2=Nv-1;
		for(n=0;n<N;n++){
			j3 = (n*h_v)/dv;
			if(j3==Nv)j3=Nv-1;
			k=i*size_v + (j1*Nv*Nv + j2*Nv + j3); // determine in which element the Fourier nodes lie	  
			f[i%chunk_Nx][l*N*N+m*N+n] = U[k*7+0] + U[k*7+3]*(v[l]-Gridv((double)j1))/dv + U[k*7+4]*(v[m]-Gridv((double)j2))/dv + U[k*7+5]*(v[n]-Gridv((double)j3))/dv + U[k*7+6]*( ((v[l]-Gridv((double)j1))/dv)*((v[l]-Gridv((double)j1))/dv) + ((v[m]-Gridv((double)j2))/dv)*((v[m]-Gridv((double)j2))/dv) + ((v[n]-Gridv((double)j3))/dv)*((v[n]-Gridv((double)j3))/dv) );
		  //BUG: index was "l*N*N+m*N+n*N" !!!!!!
		}
      }
    }
  }   
}
#else
void setInit_spectral(vector<double>& U_vals, double **f)
{
  int i, j1, j2, j3, k, l, m ,n;  
  for(l=0;l<N;l++){
      j1 = (l*h_v)/dv; // integer part = floor() for non-negative integers.
      if(j1==Nv)j1=Nv-1; // let the right end point lie in the last element
      for(m=0;m<N;m++){
	  j2 = (m*h_v)/dv;
	  if(j2==Nv)j2=Nv-1;
		for(n=0;n<N;n++){
		  j3 = (n*h_v)/dv;
		  if(j3==Nv)j3=Nv-1;
		  for(i=0;i<Nx;i++){
		  k=i*size_v + (j1*Nv*Nv + j2*Nv + j3); // determine in which element the Fourier nodes lie	  
		  f[i][l*N*N+m*N+n] = U_vals[k*7+0] + U_vals[k*7+2]*(v[l]-Gridv((double)j1))/dv + U_vals[k*7+4]*(v[m]-Gridv((double)j2))/dv + U_vals[k*7+5]*(v[n]-Gridv((double)j3))/dv + U_vals[k*7+6]*( ((v[l]-Gridv((double)j1))/dv)*((v[l]-Gridv((double)j1))/dv) + ((v[m]-Gridv((double)j2))/dv)*((v[m]-Gridv((double)j2))/dv) + ((v[n]-Gridv((double)j3))/dv)*((v[n]-Gridv((double)j3))/dv) );
		  //BUG: index was "l*N*N+m*N+n*N" !!!!!!
		  }
        }
      }
    }   
}
#endif
