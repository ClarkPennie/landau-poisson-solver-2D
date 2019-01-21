/* This is the source file which contains all functions required for the LDG method
 * to solve Poisson's equation
 *
 * 				div(epsilon_r grad(Phi(x,y))) = R(x,y)
 *
 * for 0 < x < xlength & 0 < y < ylength, where xlegnth & ylength are currently set in conf2D.ci.
 *
 * Also, 	R(x,y) = PEP(\int_(Omega_v) f(x,v,t) dv - N_D(x)),
 *
 * where f is the solution to the Boltzmann/Landau equation and N_D is some doping profile, which
 * seems to be constant at the moment.
 *
 * All code here written by Yingda Cheng/Jose Morales, but collected
 * in to one file for convenience.
 *
 * Functions included: adim, coef_pois, config, fle, fled, pois2d, setup, setup_matrix, setup_pois, transmat
 *
 *  Created on: Feb 12, 2018
 */

#include "PoissonFunctions.h"																					// PoissonFunctions.h is where the prototypes for the functions contained in this file are declared and any variables defined here to be used throughout the other files are declared as external

void adim(void)
{
	// Dirac constant - J*s (Jose's value: 1.05459e-34)
	const double h = 1;

	// electron rest mass - kg (Jose's value: 9.1091e-31)
	const double mr = 1;

	// effective mass of electron (Jose's value: mr*.32)
	const double massa = 1;

	// crystal temperature - K (Jose's value: 300.)
	const double temp = 1;

	// Boltzmann constant - J/K (Jose's value: 1.3805e-23)
	const double kb = 1;

	// Proton charge - C (Jose's value: 1.60219e-19)
	const double q = 1;

	// dielectric constant - F/m (Jose's value: 8.85419e-12)
	const double eps0 = 1;

	// relative dielectric constant in silicon (Jose's value: 11.7)
	const double ces = 1;

	// relative dielectric constant in oxide silicon (Jose's value: 3.9)
	const double ceos = 1;

	// characteristic length -m (Jose's value: 1.e-6)
	const double lcar = 1;

	// characteristic electric field - V/m = N/C 1.e5)
	const double Ecar = 1;

	// parameter for the charge density (Jose's value: pow((sqrt(2.*massa*kb*temp)/h),3.))
	PD = 1;

	// the parameter of Poisson equation (c_p)
	PEP = lcar*lcar*q*PD/eps0;

	// the parameter cv
	PVE = 1./(lcar*Ecar);

	// relative dielectric constant in silicon
	CDS = ces;

	// relative dielectric constant in oxide silicon
	CDOS = ceos;

	// print data
	if (myrank_mpi == 0)
	{
		cout << "characteristic length = " << lcar << endl;
		cout << "parameter for the charge density = " << PD << endl;
		cout << "the parameter of Poisson equation = " << PEP << endl;
	}

	// .. also in a file
	/*
	FILE *dfi;
	dfi = fopen("param_adim.dc","w");
	fprintf(dfi,"characteristic length = %g \n",lcar);
	fprintf(dfi,"parameter for the charge density = %16.12g \n",PD);
	fprintf(dfi,"the parameter of Poisson equation = %16.12g \n",PEP);
	fclose(dfi);
	*/
}

void setup(void)
{
	int er=0;

	NYREAL=NY;		// *7/6;

	if (myrank_mpi == 0)
	{
		cout << "number of cells in x direction = " << NX << endl;
		if((NX+2) > NXD)
		{
			cout << "number of cells in x direction >= " << NXD << endl;
			er = 1;
		}
		cout << "number of cells in y direction = " << NY << endl;
		if((NY+2) > NYD)
		{
			cout << "number of cells in y direction >= " << NYD << endl;
			er = 2;
		}

		if(er>0)
		{
			cout << "error = " << er << " in setup" << endl;
			exit(1);
		}
	}
}

void config(void)
{
	int i, j;
	double xi;
	const double xlength = Lx;   // in micron;
	const double ylength = Lx;   // in micron;
	const double ax = 0.;
	const double ay = 0.;

	// grid in the plane (x,y): lengths of x and y cells

	for(i=1; i<=NX; ++i)
	{
		DX[i] = xlength/((double)NX);
	}
	for(j=1; j<=NYREAL; ++j)
	{
		DY[j] = ylength/((double)NY);
	}

	// the coordinates of the center of (x,y)-cell
	CX[1] = ax + DX[1]/2.;
	for(i=2; i<=NX; ++i)
	{
		CX[i] = CX[i-1] + DX[i-1]/2. + DX[i]/2.;
	}
	CY[1] = ay + DY[1]/2.;
	for(j=2; j<=NY; ++j)
	{
		CY[j] = CY[j-1] + DY[j-1]/2. + DY[j]/2.;
	}

	//define the y grid in the oxide
	for(j=NY+1; j<=NYREAL; ++j)
	{
		CY[j] = CY[j-1] + DY[j-1]/2. + DY[j]/2.;
	}

	for(i=1; i<=NX; ++i)
	{
		xi = CX[i];
		for(j=1; j<=NY; ++j)
		{
			DOP[i][j] = 1;		// constant background density of ions
		}
	}

	/*
	FILE *fd;
	fd = fopen("dop.dc","w");

	for(j=1; j<=NY; ++j)
	{
		for(i=1; i<=NX; ++i)
		{
			fprintf(fd,"%g %g %g\n",CX[i], CY[j], PD*DOP[i][j]);
		}
		fprintf(fd," \n");
	}
	fclose(fd);
	*/

	return;
}

void setup_pois(void){


double aic[21][21];
double s3;
int i0,i01,i02,i,j,iii,l,ll;



//left and right boundary of the domain
//  xa=0.;
//  xb=1.;

  xa=CX[1]-DX[1]/2.;
  xb=CX[NX]+DX[NX]/2.;

//gate position
  xc=xa+(xb-xa)/2.;
  xd=xa+(xb-xa)/2.;

//bottom and top boundary of the domain (larger than the silicon domain)
//  ya=0.;
//  yb=1.;
  ya=CY[1]-DY[1]/2.;
  yb=CY[NYREAL]+DY[NYREAL]/2.;

//dirichlet at left (source)
//  pleft=0.52354;
  pleft = 0;	// should be the same as Chenglong's periodic BC

//dirichlet at bottom
  pbot=0;

//dirichlet at right (drain)
//  pright=1.5235;
  pright = 0;	// should be the same as Chenglong's periodic BC

//dirichlet on part of the top, from xc to xd (gate)
//  ptop=1.06;

// NO GATE (Reduced to a point / lenght zero), so:

// Potential at gate point = Midpoint for Applied Potential Bias...? CHECK


  ptop = 0.5*(pleft + pright);


//   ptop=2.0;





//number of cells
//  NX=80;
//  NY=80;

//epsilon, dielectric coef. in the silicon region (bottom part)
for(i=1;i<=NX;i++){
for(j=1;j<=NY;j++){
epsil[i][j]=CDS;
}
}

//in silicon-oxide region (top part)
for(i=1;i<=NX;i++){
for(j=NY+1;j<=NYREAL;j++){
epsil[i][j]=CDOS;
}
}


// order of the LDG scheme
  mppois=1;
 mpto=(mppois+1)*(mppois+2)/2-1;


//DEFINE NSUB, NSUPER
//FROM YUNXIAN
nsub=(NYREAL+1)*(mpto+1)-1;
nsuper=(NYREAL+1)*(mpto+1)-1;
ldab=2*nsub+nsuper+1;





//constants needed in computation
       s3=sqrt(3.0)/6.0;

//constants in the numerical flux
	c11=1.;

//array for the integration
// read in the inverse of the mass matrix when the polynomial in the cell is
// orthogonal:
      aic[0][0]=1.;
      aic[0][1]=1.;
      aic[1][1]=12.0;
      aic[2][1]=12.0;
      aic[0][2]=1.0;
      aic[1][2]=12.0;
      aic[2][2]=12.0;
      aic[4][2]=144.0;
      aic[3][2]=180.0;
      aic[5][2]=180.0;
      aic[0][3]=1.0;
      aic[1][3]=12.0;
      aic[2][3]=12.0;
      aic[3][3]=180.0;
      aic[4][3]=144.0;
      aic[5][3]=180.0;
      aic[6][3]=2800.0;
      aic[7][3]=2160.0;
      aic[8][3]=2160.0;
      aic[9][3]=2800.0;
      aic[0][4]=1.;
      aic[1][4]=12.0;
      aic[2][4]=12.0;
      aic[3][4]=180.0;
      aic[4][4]=144.0;
      aic[5][4]=180.0;
      aic[6][4]=2800.0;
      aic[7][4]=2160.0;
      aic[8][4]=2160.0;
      aic[9][4]=2800.0;
      aic[10][4]=44100.0;
      aic[11][4]=33600.0;
      aic[12][4]=32400.0;
      aic[13][4]=33600.0;
      aic[14][4]=44100.0;


//ai is the matrix for orthogonal basis
      iii=(mppois+1)*(mppois+2)/2;
      for(i=0;i<=iii-1;i++)     ai[i]=aic[i][mppois];

// mppois=1 or 3 use Gauss points, otherwise, use Gauss Lobatto points. (why)
       if(mppois==1||mppois==3)  k01=mppois+1;
       else  k01=mppois+2;


//setup the points and weights for Gauss-Lobatto quadratures
//(.,1) are points, (.,2) are weights
	if(mppois==0) {
        gauss[1][1]=-0.5;
	gauss[2][1]=0.5;
        gauss[1][2]=0.5;
	gauss[2][2]=0.5;
	}

	else if(mppois==1){
//  the points of 4h order Gauss-Lobatto quadrature
//	gauss[1][1]=-0.5;
//	gauss[2][1]=0.5;
//	gauss[3][1]=0.0;
       gauss[1][1]=-s3;
       gauss[2][1]=s3;
//   coefficients of 4h order Gauss-Lobatto quadrature
//	gauss[1][2]=1.0/6.0;
//	gauss[2][2]=gauss[1][2];
//	gauss[3][2]=2.0/3.0;
       gauss[1][2]=0.5;
       gauss[2][2]=0.5;
        }

	else if(mppois==2){
//  the points of 6th order Gauss-Lobatto quadrature
	gauss[1][1]=-0.5;
	gauss[2][1]=0.5;
	gauss[3][1]=-sqrt(5.)/10.0;
	gauss[4][1]= sqrt(5.0)/10.0;
//   coefficients of 6th order Gauss-Lobatto quadrature
	gauss[1][2]=1.0/12.0;
	gauss[2][2]=gauss[1][2];
	gauss[3][2]=5.0/12.0;
	gauss[4][2]=gauss[3][2];
        }

        else{
        	cout << "order too high in gauss points setup" << endl;
        }

for(l=1;l<=k01;l++){
for(ll=1;ll<=k01;ll++){
     gaut[l][ll]=gauss[l][2]*gauss[ll][2];
}
}
}

// the function of Legendre polynomial
double fle(int k, double x){

  double z;
  if(k==0)
    z=1.0;
  else if(k==1)
    z=x;
  else if(k==2)
    z=x*x-1./12.;
  else if(k==3)
    z=x*x*x-0.15*x;
  else if(k==4)
    z=(x*x-3./14.0)*x*x+3.0/560.0;
  else{
	  cout << "order too high in fle" << endl;
    z=0.;
  }
     return z;
}



//the function of derivative of Legendre polynomial
double fled(int k, double x){
  double z;
  if(k==0)
    z=0.;
  else if(k==1)
    z=1.;
  else if(k==2)
    z=2.*x;
  else if(k==3)
    z=3.*x*x-0.15;
  else if(k==4)
    z=(4.0*x*x-3.0/7.0)*x;
  else{
    cout << "order too high in fled" << endl;
    z=0.;
  }
     return z;
}


//To get the LHS of the big linear system
void coef_pois(void){

  int i,j,k,ii,l,kk;
  int irow, icol, ibrow,m1;
  int mmp, ll, ix,jy, kx, ky, kxy;
  double aic[MPD][MPD];

      m1=nsub+nsuper+1;


// set up value of base function at [-+1/2][gauss[l][1]] and [gauss[l][1]][-+1/2]
// to be used for computation of flux at boundary of cell
   for(l=1;l<=k01;l++){

     for(k=0;k<=mppois;k++){
     for(kk=0;kk<=mppois-k;kk++){
      aic[k][kk]=fle(k,-0.5)*fle(kk,gauss[l][1]);
        }
      }
     for(i=0;i<=mppois;i++){
      mmp=i*(i+1)/2;
     for(j=0;j<=i;j++){
      bl[mmp+j][l]=aic[i-j][j];
     }
     }


    for(k=0;k<=mppois;k++){
      for(kk=0;kk<=mppois-k;kk++){
	aic[k][kk]=fle(k,0.5)*fle(kk,gauss[l][1]);
      }
    }
    for(i=0;i<=mppois;i++){
      mmp=i*(i+1)/2;
      for(j=0;j<=i;j++){
      br[mmp+j][l]=aic[i-j][j];
      }
    }

    for(k=0;k<=mppois;k++){
      for(kk=0;kk<=mppois-k;kk++){
      aic[k][kk]=fle(kk,-0.5)*fle(k,gauss[l][1]);
      }
    }
    for(i=0;i<=mppois;i++){
      mmp=i*(i+1)/2;
     for(j=0;j<=i;j++){
      bb[mmp+j][l]=aic[i-j][j];
     }
     }


    for(k=0;k<=mppois;k++){
      for(kk=0;kk<=mppois-k;kk++){
      aic[k][kk]=fle(kk,0.5)*fle(k,gauss[l][1]);
      }
    }
    for(i=0;i<=mppois;i++){
      mmp=i*(i+1)/2;
     for(j=0;j<=i;j++){
      bt[mmp+j][l]=aic[i-j][j];
     }
     }
}


// set up the value of the partial derivative of x of base function
// at gauss points in the cell to be used for computation of integral on the cell
    for(l=1;l<=k01;l++){
    for(ll=1;ll<=k01;ll++){

     for(k=0;k<=mppois;k++){
     for(kk=0;kk<=mppois-k;kk++){
      aic[k][kk]=fled(k,gauss[l][1])*fle(kk,gauss[ll][1]);
      }
      }
     for(i=0;i<=mppois;i++){
      mmp=i*(i+1)/2;
     for(j=0;j<=i;j++){
      cx[mmp+j][l][ll]=aic[i-j][j];
      }
      }


// set up the value of the partial derivative of y of base function
// at gauss points in the cell to be used for computation of integral on the cell
     for(k=0;k<=mppois;k++){
     for(kk=0;kk<=mppois-k;kk++){
      aic[k][kk]=fle(k,gauss[l][1])*fled(kk,gauss[ll][1]);
      }
      }
     for(i=0;i<=mppois;i++){
      mmp=i*(i+1)/2;
     for(j=0;j<=i;j++){
      cy[mmp+j][l][ll]=aic[i-j][j];
      }
      }

// set up the value of the base function at gauss points in the cell
// to be used for computation of point value on the cell
     for(k=0;k<=mppois;k++){
     for(kk=0;kk<=mppois-k;kk++){
      aic[k][kk]=fle(k,gauss[l][1])*fle(kk,gauss[ll][1]);
      }
      }
     for(i=0;i<=mppois;i++){
      mmp=i*(i+1)/2;
     for(j=0;j<=i;j++){
      cxy[mmp+j][l][ll]=aic[i-j][j];
      }
      }

      }
      }
/*
     for(k=0;k<=mp;k++){
     for(kk=0;kk<=mp-k;kk++){
      aic[k][kk]=fle(k,0.)*fle(kk,0.);
      }
      }
     for(i=0;i<=mp;i++){
      mmp=i*(i+1)/2;
     for(j=0;j<=i;j++){
      co[mmp+j]=aic[i-j][j];
      }
      }
*/






// prepare the coefficents for the preconditioner matrix.
     for(j=0;j<=mpto;j++){
     for(i=0;i<=mpto;i++){
	cvx[i][j]=0.;
	cvy[i][j]=0.;
        for(kx=1;kx<=k01;kx++){
        for(ky=1;ky<=k01;ky++){
	cvx[i][j]=cvx[i][j]-cxy[j][kx][ky]*cx[i][kx][ky]*gaut[kx][ky];
	cvy[i][j]=cvy[i][j]-cxy[j][kx][ky]*cy[i][kx][ky]*gaut[kx][ky];
	}
	}

	}
	}

     for(j=0;j<=mpto;j++){
     for(i=0;i<=mpto;i++){
	cyll[i][j]=0.;
	cylr[i][j]=0.;
	cyrl[i][j]=0.;
	cyrr[i][j]=0.;
	cxll[i][j]=0.;
	cxlr[i][j]=0.;
	cxrl[i][j]=0.;
	cxrr[i][j]=0.;
        for(kxy=1;kxy<=k01;kxy++){
	cyll[i][j]=cyll[i][j]+bl[i][kxy]*bl[j][kxy]*gauss[kxy][2];
	cylr[i][j]=cylr[i][j]+br[i][kxy]*bl[j][kxy]*gauss[kxy][2];
	cyrl[i][j]=cyrl[i][j]+bl[i][kxy]*br[j][kxy]*gauss[kxy][2];
	cyrr[i][j]=cyrr[i][j]+br[i][kxy]*br[j][kxy]*gauss[kxy][2];
	cxll[i][j]=cxll[i][j]+bb[i][kxy]*bb[j][kxy]*gauss[kxy][2];
	cxlr[i][j]=cxlr[i][j]+bt[i][kxy]*bb[j][kxy]*gauss[kxy][2];
	cxrl[i][j]=cxrl[i][j]+bb[i][kxy]*bt[j][kxy]*gauss[kxy][2];
	cxrr[i][j]=cxrr[i][j]+bt[i][kxy]*bt[j][kxy]*gauss[kxy][2];
	}
//	if(ix==1&&jy==1) printf("00, cyll %8.6g, cyrl %8.6g\n ", cyll[0][0], cyrl[0][0]);

	}
	}

//	if(ix==1&&jy==1) printf("star, br %8.6g %8.6g, bl %8.6g %8.6g\n ", br[0][1], br[0][2],  bl[0][1], bl[0][2] );
//	if(ix==1&&jy==1) printf("0, cyll %8.6g, cyrl %8.6g %8.6g\n ", cyll[0][0], cyrl[0][0],DY[jy]);
/*
     for(j=0;j<=mpto;j++){
	vxr[j][ix]=0.;
        for(kx=1;kx<=k01;kx++){
	vxr[j][ix]=vxr[j][ix]+bt[j][kx]*gauss[kx][2];
	}
	vxr[j][ix]=vxr[j][ix]*DX[ix];
	}
*/

     for(j=0;j<=mpto;j++){
     for(i=0;i<=mpto;i++){
	acvx[i][j]=ai[i]*cvx[i][j];
	acvy[i][j]=ai[i]*cvy[i][j];
	acyll[i][j]=ai[i]*cyll[i][j];
	acylr[i][j]=ai[i]*cylr[i][j];
	acyrl[i][j]=ai[i]*cyrl[i][j];
	acyrr[i][j]=ai[i]*cyrr[i][j];
	acxll[i][j]=ai[i]*cxll[i][j];
	acxlr[i][j]=ai[i]*cxlr[i][j];
	acxrr[i][j]=ai[i]*cxrr[i][j];
	acxrl[i][j]=ai[i]*cxrl[i][j];
	}
	}


     for(ix=1;ix<=NX;ix++){
       for(jy=1;jy<=NYREAL;jy++){
         ii=(ix-1)*NYREAL+jy;


/* Construct the preconditioner matrix.
* c col is the matrix multiplying uc(k,i-1,j),cob is the matrix multiplying  uc(k,i,j-1),
* cott is the matrix multiplying uc(k,i,j+1),cor is the matrix multiplying uc(k,i+1,j),
* coc is the matrix multiplying uc(k,i,j).
*/

     for(j=0;j<=mpto;j++){
     for(i=0;i<=mpto;i++){
	col[i][j]=0.;
	cob[i][j]=0.;
	coc[i][j]=0.;
	cott[i][j]=0.;
	cor[i][j]=0.;
     for(k=0;k<=mpto;k++){
//check index for epsil, why is Yunxian's code reversed(flux might be different)
	col[i][j]=col[i][j]-(cvx[i][k]-cyll[i][k])*epsil[ix][jy]*acyrl[k][j]*DY[jy]/DX[ix];
//	if(ix==1&&jy==1) printf("1, coc %8.6g, col %8.6g, cor %8.6g, cob %8.6g, cot %8.6g \n ",cvx[i][k], cyll[i][k], epsil[ix][jy], acyrl[k][j], DX[ix]);
	cob[i][j]=cob[i][j]-(cvy[i][k]-cxll[i][k])*epsil[ix][jy]*acxrl[k][j]*DX[ix]/DY[jy];
	cott[i][j]=cott[i][j]+cxlr[i][k]*epsil[ix][jy+1]*(acvy[k][j]+acxrr[k][j])*DX[ix]/DY[jy+1];
	cor[i][j]=cor[i][j]+cylr[i][k]*(acvx[k][j]+acyrr[k][j])*epsil[ix+1][jy]*DY[jy]/DX[ix+1];
	coc[i][j]=coc[i][j]
     +(cvx[i][k]-cyll[i][k])*epsil[ix][jy]*(acvx[k][j]+acyrr[k][j])*DY[jy]/DX[ix]
     +(cvy[i][k]-cxll[i][k])*epsil[ix][jy]*(acvy[k][j]+acxrr[k][j])*DX[ix]/DY[jy]
     -cylr[i][k]*acyrl[k][j]*epsil[ix+1][jy]*DY[jy]/DX[ix+1]-cxlr[i][k]*acxrl[k][j]*epsil[ix][jy+1]
     *DX[ix]/DY[jy+1];
	}
//	if(ix==1&&jy==1&&i==0&&j==1) printf("1, %8.6g, %8.6g \n ",col[0][1], cob[0][1]);
	col[i][j]=col[i][j]-c11*cyrl[i][j]*DY[jy];
	cob[i][j]=cob[i][j]-c11*cxrl[i][j]*DX[ix];
	coc[i][j]=coc[i][j]+c11*(cyrr[i][j]+cyll[i][j])*DY[jy]+c11*(cxrr[i][j]+cxll[i][j])*DX[ix];
	cott[i][j]=cott[i][j]-c11*cxlr[i][j]*DX[ix];
	cor[i][j]=cor[i][j]-c11*cylr[i][j]*DY[jy];





//Neumann boundary condition at jy=NYREAL, x \in [xa, xc]&[xd, xb]
//the top right corner need special treatment as below
     if(jy==NYREAL&&ix!=NX&&(CX[ix]<xc||CX[ix]>xd)){
	coc[i][j]=0.;
     for(k=0;k<=mpto;k++){
	coc[i][j]=coc[i][j]
     +(cvx[i][k]-cyll[i][k])*epsil[ix][jy]*(acvx[k][j]+acyrr[k][j])*DY[jy]/DX[ix]
     +(cvy[i][k]-cxll[i][k])*epsil[ix][jy]*(acvy[k][j]+acxrr[k][j])*DX[ix]/DY[jy]
     -cylr[i][k]*acyrl[k][j]*epsil[ix+1][jy]*DY[jy]/DX[ix+1];
      }
	coc[i][j]=coc[i][j]+c11*(cyrr[i][j]+cyll[i][j])*DY[jy]+c11*cxll[i][j]*DX[ix];
     }

//Dirichelet condition at jy=NYREAL, x \in [xc,xd]
//flip the flux in y direction at top end point (similar to 1d case)
     if(jy==NYREAL&&CX[ix]>=xc&&CX[ix]<=xd){
	coc[i][j]=0.;
	cob[i][j]=0.;
     for(k=0;k<=mpto;k++){
	cob[i][j]=cob[i][j]-(cvy[i][k]-cxll[i][k]+cxrr[i][k])*epsil[ix][jy]*acxrl[k][j]*DX[ix]/DY[jy];
	coc[i][j]=coc[i][j]
     +(cvx[i][k]-cyll[i][k])*epsil[ix][jy]*(acvx[k][j]+acyrr[k][j])*DY[jy]/DX[ix]
     +(cvy[i][k]-cxll[i][k]+cxrr[i][k])*epsil[ix][jy]*acvy[k][j]*DX[ix]/DY[jy]
     -cylr[i][k]*acyrl[k][j]*epsil[ix+1][jy]*DY[jy]/DX[ix+1];
      }
	cob[i][j]=cob[i][j]-c11*cxrl[i][j]*DX[ix];
	coc[i][j]=coc[i][j]+c11*(cyrr[i][j]+cyll[i][j])*DY[jy]+c11*(cxrr[i][j]+cxll[i][j])*DX[ix];
     }

     if(jy==NYREAL-1&&CX[ix]>=xc&&CX[ix]<=xd){
	cott[i][j]=0.;
     for(k=0;k<=mpto;k++){
	cott[i][j]=cott[i][j]+cxlr[i][k]*epsil[ix][jy+1]*acvy[k][j]*DX[ix]/DY[jy+1];
      }
	cott[i][j]=cott[i][j]-c11*cxlr[i][j]*DX[ix];
     }

//Dirichelet condition at ix=NX, needs to modify NX, NX-1
//flip the flux in x direction at the right end point (similar to 1d case)

     if(ix==NX&&jy!=NYREAL){
	coc[i][j]=0.;
	col[i][j]=0.;
     for(k=0;k<=mpto;k++){
	col[i][j]=col[i][j]-(cvx[i][k]-cyll[i][k]+cyrr[i][k])*epsil[ix][jy]*acyrl[k][j]*DY[jy]/DX[ix];
	coc[i][j]=coc[i][j]
     +(cvy[i][k]-cxll[i][k])*epsil[ix][jy]*(acvy[k][j]+acxrr[k][j])*DX[ix]/DY[jy]
     +(cvx[i][k]-cyll[i][k]+cyrr[i][k])*epsil[ix][jy]*acvx[k][j]*DY[jy]/DX[ix]-cxlr[i][k]*acxrl[k][j]*epsil[ix][jy+1]
     *DX[ix]/DY[jy+1];
      }
	col[i][j]=col[i][j]-c11*cyrl[i][j]*DY[jy];
	coc[i][j]=coc[i][j]+c11*(cyrr[i][j]+cyll[i][j])*DY[jy]+c11*(cxrr[i][j]+cxll[i][j])*DX[ix];
     }

     if(ix==NX&&jy==NYREAL){
	coc[i][j]=0.;
	col[i][j]=0.;
     for(k=0;k<=mpto;k++){
	col[i][j]=col[i][j]-(cvx[i][k]-cyll[i][k]+cyrr[i][k])*epsil[ix][jy]*acyrl[k][j]*DY[jy]/DX[ix];
	coc[i][j]=coc[i][j]
     +(cvy[i][k]-cxll[i][k])*epsil[ix][jy]*(acvy[k][j]+acxrr[k][j])*DX[ix]/DY[jy]
     +(cvx[i][k]-cyll[i][k]+cyrr[i][k])*epsil[ix][jy]*acvx[k][j]*DY[jy]/DX[ix];
      }
	col[i][j]=col[i][j]-c11*cyrl[i][j]*DY[jy];
	coc[i][j]=coc[i][j]+c11*(cyrr[i][j]+cyll[i][j])*DY[jy]+c11*cxll[i][j]*DX[ix];
     }


     if(ix==NX-1){
	cor[i][j]=0.;
     for(k=0;k<=mpto;k++){
	cor[i][j]=cor[i][j]+cylr[i][k]*acvx[k][j]*epsil[ix+1][jy]*DY[jy]/DX[ix+1];
      }
	cor[i][j]=cor[i][j]-c11*cylr[i][j]*DY[jy];
     }


//neumann boundary condition at jy=1
//ONLY WORKS for p^1

     if(jy==1){
//this has been changed recently to adapt the neumann boundary conditions
/*
	coc[i][j]=0.;
     for(k=0;k<=mpto;k++){
	coc[i][j]=coc[i][j]
     +(cvx[i][k]-cyll[i][k])*epsil[ix][jy]*(acvx[k][j]+acyrr[k][j])*DY[jy]/DX[ix]
     +cvy[i][k]*epsil[ix][jy]*(acvy[k][j]+acxrr[k][j]-acxll[k][j])*DX[ix]/DY[jy]
     -cylr[i][k]*acyrl[k][j]*epsil[ix+1][jy]*DY[jy]/DX[ix+1]-cxlr[i][k]*acxrl[k][j]*epsil[ix][jy+1]
     *DX[ix]/DY[jy+1];
	}
	coc[i][j]=coc[i][j]+c11*(cyrr[i][j]+cyll[i][j])*DY[jy]+c11*cxrr[i][j]*DX[ix];
*/
     for(k=0;k<=mpto;k++){
	coc[i][j]=coc[i][j]-(cvy[i][k]-cxll[i][k])*epsil[ix][jy]*(acvy[k][j]+acxrr[k][j])*DX[ix]/DY[jy]+cvy[i][k]*epsil[ix][jy]*(acvy[k][j]+acxrr[k][j]-acxll[k][j])*DX[ix]/DY[jy];
      }
     coc[i][j]=coc[i][j]-c11*cxll[i][j]*DX[ix];

     }


	}
	}

/*
	if(ix==2&&jy==2){
        for(j=0;j<=mpto;j++){
         for(i=0;i<=mpto;i++){
           printf(" coc i=%i, j=%i, %8.6g\n ",i,j,coc[i][j]);

}
}
        for(j=0;j<=mpto;j++){
         for(i=0;i<=mpto;i++){
           printf(" col i=%i, j=%i, %8.6g\n ",i,j,col[i][j]);

}
}
        for(j=0;j<=mpto;j++){
         for(i=0;i<=mpto;i++){
           printf(" cor i=%i, j=%i, %8.6g\n ",i,j,cor[i][j]);

}
}
        for(j=0;j<=mpto;j++){
         for(i=0;i<=mpto;i++){
           printf(" cob i=%i, j=%i, %8.6g\n ",i,j,cob[i][j]);

}
}
        for(j=0;j<=mpto;j++){
         for(i=0;i<=mpto;i++){
           printf(" cot i=%i, j=%i, %8.6g\n ",i,j,cott[i][j]);

}
}
}
*/

//	if(ix==1&&jy==1) printf("coc %8.6g, col %8.6g, cor %8.6g, cob %8.6g, cot %8.6g \n ",coc[0][0], col[0][0], cor[0][0], cob[0][0], cott[0][0]);
//	if(ix==1&&jy==1) printf("2, coc %8.6g, col %8.6g, cor %8.6g, cob %8.6g, cot %8.6g \n ",coc[0][0], col[0][0], cor[0][0], cob[0][0], cott[0][0]);
//	if(ix==1&&jy==1) printf("2, coc %8.6g, col %8.6g, cor %8.6g, cob %8.6g, cot %8.6g \n ",coc[0][1], col[0][1], cor[0][1], cob[0][1], cott[0][1]);
//	if(ix==1&&jy==1) printf("2, coc %8.6g, col %8.6g, cor %8.6g, cob %8.6g, cot %8.6g \n ",coc[1][0], col[1][0], cor[1][0], cob[1][0], cott[1][0]);
//	if(ix==NX&&jy==2) printf("2, coc %8.6g, col %8.6g, cor %8.6g, cob %8.6g, cot %8.6g \n ",coc[0][0], col[0][0], cor[0][0], cob[0][0], cott[0][0]);
//	if(ix==2&&jy==2) printf("2, coc %8.6g, col %8.6g, cor %8.6g, cob %8.6g, cot %8.6g \n ",coc[1][1], col[1][1], cor[1][1], cob[1][1], cott[1][1]);
//	if(ix==2&&jy==2) printf("2, coc %8.6g, col %8.6g, cor %8.6g, cob %8.6g, cot %8.6g \n ",coc[0][2], col[0][2], cor[0][2], cob[0][2], cott[0][2]);
//	if(ix==2&&jy==2) printf("2, coc %8.6g, col %8.6g, cor %8.6g, cob %8.6g, cot %8.6g \n ",coc[1][2], col[1][2], cor[1][2], cob[1][2], cott[1][2]);

// Change a matrix of block tridiagonal to a banded storage through general storage
        if(ii==1){
	 for(j=0;j<=mpto;j++){
         for(k=0;k<=mpto;k++){
         for(ll=0;ll<=1;ll++){
	   irow=j+1;
	   icol=ll*(mpto+1)+k+1;
	   ibrow=irow-icol+m1;
	   if (ll==0)    poisab[ibrow][icol]=coc[j][k];
	   else poisab[ibrow][icol]=cott[j][k];
         }
         irow=j+1;
	 icol=NYREAL*(mpto+1)+k+1;
	 ibrow=irow-icol+m1;
	 poisab[ibrow][icol]=cor[j][k];
       }
       }
       }//if ii=1

        else if(ii<=NYREAL-1){
	 for(j=0;j<=mpto;j++){
         for(k=0;k<=mpto;k++){
         for(ll=0;ll<=2;ll++){
	   irow=(ii-1)*(mpto+1)+j+1;
	   icol=(ii-2+ll)*(mpto+1)+k+1;
	   ibrow=irow-icol+m1;
	   if (ll==0)    poisab[ibrow][icol]=cob[j][k];
	   else if(ll==1) poisab[ibrow][icol]=coc[j][k];
	   else poisab[ibrow][icol]=cott[j][k];
         }
	   irow=(ii-1)*(mpto+1)+j+1;
	   icol=(ii-1+NYREAL)*(mpto+1)+k+1;
	 ibrow=irow-icol+m1;
	 poisab[ibrow][icol]=cor[j][k];
       }
       }

       }//elseif

        else if(ii==NYREAL){
	 for(j=0;j<=mpto;j++){
         for(k=0;k<=mpto;k++){
         for(ll=0;ll<=1;ll++){
	   irow=(ii-1)*(mpto+1)+j+1;
	   icol=(ii-2+ll)*(mpto+1)+k+1;
	   ibrow=irow-icol+m1;
	   if (ll==0)    poisab[ibrow][icol]=cob[j][k];
	   else  poisab[ibrow][icol]=coc[j][k];
         }
	   irow=(ii-1)*(mpto+1)+j+1;
	   icol=(ii-1+NYREAL)*(mpto+1)+k+1;
	 ibrow=irow-icol+m1;
	 poisab[ibrow][icol]=cor[j][k];
       }
       }
       }//elseif

        else if(ii<=NYREAL*(NX-1)){

         if(jy==1){
//no bottom cells
	 for(j=0;j<=mpto;j++){
         for(k=0;k<=mpto;k++){
         for(ll=1;ll<=2;ll++){
	   irow=(ii-1)*(mpto+1)+j+1;
	   icol=(ii-2+ll)*(mpto+1)+k+1;
	   ibrow=irow-icol+m1;
	   if (ll==1)    poisab[ibrow][icol]=coc[j][k];
	   else poisab[ibrow][icol]=cott[j][k];
         }
	   irow=(ii-1)*(mpto+1)+j+1;
	   icol=(ii-1+NYREAL)*(mpto+1)+k+1;
	 ibrow=irow-icol+m1;
	 poisab[ibrow][icol]=cor[j][k];

	   irow=(ii-1)*(mpto+1)+j+1;
	   icol=(ii-1-NYREAL)*(mpto+1)+k+1;
	 ibrow=irow-icol+m1;
	 poisab[ibrow][icol]=col[j][k];
       }
       }

         }



         else if(jy==NYREAL){
//no top cells
	 for(j=0;j<=mpto;j++){
         for(k=0;k<=mpto;k++){
         for(ll=0;ll<=1;ll++){
	   irow=(ii-1)*(mpto+1)+j+1;
	   icol=(ii-2+ll)*(mpto+1)+k+1;
	   ibrow=irow-icol+m1;
	   if (ll==0)    poisab[ibrow][icol]=cob[j][k];
	   else poisab[ibrow][icol]=coc[j][k];
         }
	   irow=(ii-1)*(mpto+1)+j+1;
	   icol=(ii-1+NYREAL)*(mpto+1)+k+1;
	 ibrow=irow-icol+m1;
	 poisab[ibrow][icol]=cor[j][k];

	   irow=(ii-1)*(mpto+1)+j+1;
	   icol=(ii-1-NYREAL)*(mpto+1)+k+1;
	 ibrow=irow-icol+m1;
	 poisab[ibrow][icol]=col[j][k];
       }
       }

         }

         else{
	 for(j=0;j<=mpto;j++){
         for(k=0;k<=mpto;k++){
         for(ll=0;ll<=2;ll++){
	   irow=(ii-1)*(mpto+1)+j+1;
	   icol=(ii-2+ll)*(mpto+1)+k+1;
	   ibrow=irow-icol+m1;
	   if (ll==0)    poisab[ibrow][icol]=cob[j][k];
	   else if(ll==1) poisab[ibrow][icol]=coc[j][k];
	   else poisab[ibrow][icol]=cott[j][k];
         }
	   irow=(ii-1)*(mpto+1)+j+1;
	   icol=(ii-1+NYREAL)*(mpto+1)+k+1;
	 ibrow=irow-icol+m1;
	 poisab[ibrow][icol]=cor[j][k];

	   irow=(ii-1)*(mpto+1)+j+1;
	   icol=(ii-1-NYREAL)*(mpto+1)+k+1;
	 ibrow=irow-icol+m1;
	 poisab[ibrow][icol]=col[j][k];
       }
       }
       }//else

       }//elseif


        else if(ii==NYREAL*(NX-1)+1){
//no bot, no right
	 for(j=0;j<=mpto;j++){
         for(k=0;k<=mpto;k++){
         for(ll=1;ll<=2;ll++){
	   irow=(ii-1)*(mpto+1)+j+1;
	   icol=(ii-2+ll)*(mpto+1)+k+1;
	   ibrow=irow-icol+m1;
	   if (ll==1)    poisab[ibrow][icol]=coc[j][k];
	   else poisab[ibrow][icol]=cott[j][k];
         }
         irow=(ii-1)*(mpto+1)+j+1;
	   icol=(ii-1-NYREAL)*(mpto+1)+k+1;
	 ibrow=irow-icol+m1;
	 poisab[ibrow][icol]=col[j][k];
       }
       }


       }//elseif

        else if(ii<NYREAL*NX){
	 for(j=0;j<=mpto;j++){
         for(k=0;k<=mpto;k++){
         for(ll=0;ll<=2;ll++){
	   irow=(ii-1)*(mpto+1)+j+1;
	   icol=(ii-2+ll)*(mpto+1)+k+1;
	   ibrow=irow-icol+m1;
	   if (ll==0)    poisab[ibrow][icol]=cob[j][k];
	   else if(ll==1) poisab[ibrow][icol]=coc[j][k];
	   else poisab[ibrow][icol]=cott[j][k];
         }
         irow=(ii-1)*(mpto+1)+j+1;
	   icol=(ii-1-NYREAL)*(mpto+1)+k+1;
	 ibrow=irow-icol+m1;
	 poisab[ibrow][icol]=col[j][k];
       }
       }

       }//elseif

      else{
//no top, no right
	 for(j=0;j<=mpto;j++){
         for(k=0;k<=mpto;k++){
         for(ll=0;ll<=1;ll++){
	   irow=(ii-1)*(mpto+1)+j+1;
	   icol=(ii-2+ll)*(mpto+1)+k+1;
	   ibrow=irow-icol+m1;
	   if (ll==0)    poisab[ibrow][icol]=cob[j][k];
	   else poisab[ibrow][icol]=coc[j][k];
         }
         irow=(ii-1)*(mpto+1)+j+1;
	   icol=(ii-1-NYREAL)*(mpto+1)+k+1;
	 ibrow=irow-icol+m1;
	 poisab[ibrow][icol]=col[j][k];
       }
       }





    }

   }
   }


}




void transmat(void){
// to transfer the LHS ab to a 1d matrix that works for C
//double check the size of ab
int i,j;
for(j=1;j<=(mpto+1)*NX*NYREAL;j++)
for(i=1;i<=2*nsub+nsuper+1;i++)
    poiscab[(j-1)*(2*nsub+nsuper+1)+i-1]=poisab[i][j];

}



void setup_matrix(void){
//setup the LHS of big linear system


ndim=(mpto+1)*NX*NYREAL;
ldab=2*nsub+nsuper+1;



ipiv=(int *)malloc(ndim*sizeof(int));
info=(int *)malloc(sizeof(int));
poiscab=(double *)malloc(ldab*ndim*sizeof(double));

coef_pois();
transmat();

#ifdef HAVE_MKL
dgbtrf(&ndim, &ndim, &nsub, &nsuper, poiscab, &ldab, ipiv, info);
#elif HAVE_OPENBLAS
dgbtrf_(&ndim, &ndim, &nsub, &nsuper, poiscab, &ldab, ipiv, info);
#else
printf("Error: unsupported BLAS configuration\n");
exit(1);
#endif
if(*info!=0) cout << "error in factorizing " << *info << endl;

}

void pois2d(vector<double>& Ut, vector<double>& POTC, vector<double>& phix, vector<double>& phiy)
{
	int i,j,k,kk,ix,jy,m,n,ii;
	int i1NxNvNvNv, i2NvNvNv, j1NvNv, j2Nv;
	int ii_l, j_l;
    double w;
    int nrhs;
    char trans='N';
    double rt[NXD][NYD2], rx[NXD][NYD2], ry[NXD][NYD2];
    double tempcol, tempcob;

    //RHS of poisson
    // The coefficients for the charge density. See Eq. (36).
    // Note that when integrating the DG approximation here, since it only uses centered linear basis
    // functions in velocity, they integrate to zero, so the solution is the space part times dv^3.
    for(int i1=0; i1<NX; i1++)
    {
    	i1NxNvNvNv = i1*Nx*Nv*Nv*Nv;
    	for(int i2=0; i2<NY; i2++)
    	{
        	i2NvNvNv = i2*Nv*Nv*Nv;
    		rt[i1+1][i2+1] = rx[i1+1][i2+1] = ry[i1+1][i2+1] = 0.;
    		for(int j1=0; j1<Nv; j1++)
    		{
    			j1NvNv = j1*Nv*Nv;
    			for(int j2=0; j2<Nv; j2++)
    			{
    				j2Nv = j2*Nv;
    				for(int j3=0; j3<Nv; j3++)
    				{
    					k = i1NxNvNvNv + i2NvNvNv + j1NvNv + j2Nv + j3;
    					rt[i1+1][i2+1] += scalev*(Ut[7*k+0] + 0.25*Ut[7*k+6]);
    					rx[i1+1][i2+1] += scalev*Ut[7*k+1];		// Note that Poisson basis functions are the same as Chenglong's even though Jose's are 2 times Chenglong's
    					ry[i1+1][i2+1] += scalev*Ut[7*k+2];
    				}
    			}
    		}
    	}
    }



    //for every cell, the integral of rhs times test function
    for(ix=1;ix<=NX;ix++)
    {
    	pybc = POT[ix][1];	 // DEBUG Fix y boundary conditions
    	/*
    	if(myrank_mpi==0)
    	{
    		printf("ix = %d: pybc = %g; ", ix, pybc);
    	}
    	*/
    	for(jy=1;jy<=NYREAL;jy++)
    	{
    		/*
			 for(k=0;k<=mpto;k++)
			 {
				 ff[k][ix][jy]=0.;
				 for(kx=1;kx<=k01;kx++)
				 {
					 for(ky=1;ky<=k01;ky++)
					 {
						 ff[k][ix][jy]=ff[k][ix][jy]
							 +rhs(CX[ix]+gauss[kx][1]*DX[ix],CY[jy]+gauss[ky][1]*DY[jy])
							 *cxy[k][kx][ky]*gaut[kx][ky];
					 }
				 }
    		 */

    		for(k=0;k<=mpto;k++)
			{
    			ff[k][ix][jy]=0.;
			}
    		//only works for p^1
    		//silicon region
    		//ff[l][ix][jy]=\int_{I_(ix,jy)} (rt[ix][jy]+rx[ix][jy](nu_1(x))+ry[ix][jy](nu_2(y)) - DOP[ix][jy]) nu_l(x,y)) dxdy
    		if(jy<=NY)
    		{
    			ff[0][ix][jy]=PEP*(rt[ix][jy]-DOP[ix][jy])*DX[ix]*DY[jy];
    			ff[1][ix][jy]=PEP*rx[ix][jy]*DX[ix]*DY[jy]/12.;
    			ff[2][ix][jy]=PEP*ry[ix][jy]*DX[ix]*DY[jy]/12.;
    		}

    		for(k=0;k<=mpto;k++)
    		{
    			//for left boundaries, dirichlet condition
    			if(ix==1)
    			{
    				tempcol=0.;
    				for(kk=0;kk<=mpto;kk++)
    				{
    					tempcol=tempcol-(cvx[k][kk]-cyll[k][kk])*epsil[ix][jy]*acyrl[kk][0]*DY[jy]/DX[ix];
    				}
    				tempcol=tempcol-c11*cyrl[k][0]*DY[jy];
    				ff[k][ix][jy]=ff[k][ix][jy]-tempcol*pleft;
    			}


    			//for bottom boundaries, dirichlet condition
    			/* try no Dirichlet on bottom...
    			if(jy==1)
    			{
    				tempcob=0.;
    				for(kk=0;kk<=mpto;kk++)
    				{
    					tempcob=tempcob-(cvy[k][kk]-cxll[k][kk])*epsil[ix][jy]*acxrl[kk][0]*DX[ix]/DY[jy];
    				}
    				tempcob=tempcob-c11*cxrl[k][0]*DX[ix];

//    				ff[k][ix][jy]=ff[k][ix][jy]-tempcob*pbot;	// DEBUG Fix y boundary conditions
    				ff[k][ix][jy]=ff[k][ix][jy]-tempcob*pybc;	// DEBUG Fix y boundary conditions
    			}
    			*/


    			//for right boundaries, dirichlet condition
    			if(ix==NX)
    			{
    				tempcob=0.;
    				for(kk=0;kk<=mpto;kk++)
    				{
    					tempcob=tempcob-(cvx[k][kk]+cyrr[k][kk]-cyll[k][kk])*epsil[ix][jy]*acyrr[kk][0]*DY[jy]/DX[ix];
    				}
    				tempcob=tempcob+c11*cyrr[k][0]*DY[jy];

    				ff[k][ix][jy]=ff[k][ix][jy]+tempcob*pright;
    			}

    			if(ix==NX-1)
    			{
    				tempcob=0.;
    				for(kk=0;kk<=mpto;kk++)
    				{
    					tempcob=tempcob-cylr[k][kk]*epsil[ix+1][jy]*acyrr[kk][0]*DY[jy]/DX[ix+1];
    				}


    				ff[k][ix][jy]=ff[k][ix][jy]+tempcob*pright;
    			}


    			//for top boundaries, dirichlet condition in [xc,xd]
    			if(jy==NYREAL&&CX[ix]>=xc&&CX[ix]<=xd)
    			{
    				tempcob=0.;
    				for(kk=0;kk<=mpto;kk++)
    				{
    					tempcob=tempcob-(cvy[k][kk]+cxrr[k][kk]-cxll[k][kk])*epsil[ix][jy]*acxrr[kk][0]*DX[ix]/DY[jy];
    				}
    				tempcob=tempcob+c11*cxrr[k][0]*DX[ix];

    				//ff[k][ix][jy]=ff[k][ix][jy]+tempcob*ptop;	// DEBUG Fix y boundary conditions
    				ff[k][ix][jy]=ff[k][ix][jy]+tempcob*ptop;	// DEBUG Fix y boundary conditions
    			}

    			if(jy==NYREAL-1&&CX[ix]>=xc&&CX[ix]<=xd)
    			{
    				tempcob=0.;
    				for(kk=0;kk<=mpto;kk++)
    				{
    					tempcob=tempcob-cxlr[k][kk]*epsil[ix][jy+1]*acxrr[kk][0]*DX[ix]/DY[jy+1];
    				}
    				//ff[k][ix][jy]=ff[k][ix][jy]+tempcob*ptop;	// DEBUG Fix y boundary conditions
    				ff[k][ix][jy]=ff[k][ix][jy]+tempcob*pybc;	// DEBUG Fix y boundary conditions
    			}
    		}
    	}

    	//  	if(ix==2&&jy==2) printf("2, coefr %8.6g, %8.6g,%8.6g\n ",ff[0][2][2],ff[1][2][2],ff[2][2][2]);

    }
    //s to change ff(k,i) to poisb(j)

    for(ix=1;ix<=NX;ix++)
    {
    	for(jy=1;jy<=NYREAL;jy++)
    	{
    		for(k=0;k<=mpto;k++)
    		{
    			ii=(ix-1)*NYREAL+jy;
    			j=(ii-1)*(mpto+1)+k;
    			poisb[j]=ff[k][ix][jy];
    			//    printf("%i, %i %i %i, %8.6g \n ",j, k,ix,jy,ff[k][ix][jy]);
    		}
    	}
    }

    //for(i=0;i<ldab*ndim;i++)
    //printf(" %i, %g \n", i, poiscab[i]);

    /* dgbtrs is a library function from mkl_lapack.h, which solves a system of linear equations with an
     * LU-factored band coefficient matrix A, with multiple right-hand sides contained as the columns of B.
     *
     * The general syntax is:
     *
     * dgbtrs (char trans, lapack_int n, lapack_int kl, lapack_int ku, lapack_int nrhs, const double* ab,
     * 			 lapack_int ldab, const lapack_int* ipiv, double* b, lapack_int ldb);
     *
     * where: 	- trans determines if the matrix is tranpose (N for no, T for yes, or C for yes, complex and adjoint)
     * 			- n is the number of rows in the matrices A & B
     * 			- kl is the number of sub-diagonals in the band of A
     * 			- ku is the number of super-diagonals in the band of A
     * 			- nrhs is the number of right-hand sides (i.e. columns of B)
     * 			- ab contains the elements of the LU factorisation of A, calculated by the mkl_lapack function dgbtrf
     * 			- ldab is the leading dimension of ab
     * 			- ipiv is the array containing the pivot indices of the LU factorisation, output by dgbtrf
     * 			- b is the matrix B of right-hand sides and is overwritten by the solution to the system
     * 			- ldb is the leading dimension of B
     * 			- info determines if the solve was successful (0 if yes, -i if there was an "illegal value")
     *
     * Here, the system to be solved is essentially A*POTC = ff, so there is only one RHS.
     *
     * - trans=N is set at the top, here
     * - n=ndim is set in coef_pois.c
     * - kl=nsub, ku=super & ldab are set in coef_pois.c
     * - nrhs=1 is set below
     * - poiscab is assembled and calculated by dgbtrf, along with ipiv in coef_pois.c
     * - poisb is assembled above using ff
     *
     */
    nrhs=1;
    #ifdef HAVE_MKL
    dgbtrs(&trans, &ndim, &nsub, &nsuper, &nrhs, poiscab, &ldab, ipiv, poisb, &ndim, info);
    #elif HAVE_OPENBLAS
    dgbtrs_(&trans, &ndim, &nsub, &nsuper, &nrhs, poiscab, &ldab, ipiv, poisb, &ndim, info);
    #else
    printf("Error: unsupported BLAS configuration\n");
    exit(1);
    #endif
    if(*info!=0)
    {
    	cout << "error in solving system " << *info << endl;
    }

	j=1;

	//transform the solution back
	for(ix=0;ix<NX;ix++)
	{
		for(jy=0;jy<NYREAL;jy++)
		{
			ii=ix*NYREAL+jy;
			for(k=0;k<mpto+1;k++)
			{
				POTC[3*ii+k]=poisb[3*ii+k];
			}
			POT[ix+1][jy+1]=POTC[3*ii];
		}
	}

	//boundaries

	//get E_X
	for(ix=2;ix<=NX-1;ix++)
	{
		for(jy=1;jy<=NYREAL;jy++)
		{
			ii = (ix-1)*NX + jy-1;
			ii_l = (ix-2)*NX + jy-1;
			for(k=0;k<=mpto;k++)
			{
				phix[3*ii+k]=0.;
				for(kk=0;kk<=mpto;kk++)
				{
					phix[3*ii+k]=phix[3*ii+k]+(acvx[k][kk]+acyrr[k][kk])*POTC[3*ii+kk]-acyrl[k][kk]*POTC[3*ii_l+kk];
				}
				phix[3*ii+k]=-phix[3*ii+k]/DX[ix];
			}
		}
	}

	for(jy=1;jy<=NYREAL;jy++)
	{
		ix=NX;
		ii = (ix-1)*NX + jy-1;
		ii_l = (ix-2)*NX + jy-1;
		for(k=0;k<=mpto;k++)
		{
			phix[3*ii+k]=0.;
			for(kk=0;kk<=mpto;kk++)
			{
				phix[3*ii+k]=phix[3*ii+k]+acvx[k][kk]*POTC[3*ii+kk]-acyrl[k][kk]*POTC[3*ii_l+kk];
			}

			phix[3*ii+k]=phix[3*ii+k]+acyrr[k][0]*pright;
			phix[3*ii+k]=-phix[3*ii+k]/DX[ix];
		}

		ix=1;
		ii = (ix-1)*NX + jy-1;
		for(k=0;k<=mpto;k++)
		{
			phix[3*ii+k]=0.;
			for(kk=0;kk<=mpto;kk++)
			{
				phix[3*ii+k]=phix[3*ii+k]+(acvx[k][kk]+acyrr[k][kk])*POTC[3*ii+kk];
			}
			//printf("k %8.6g\n ",acyrl[k][0]);
			phix[3*ii+k]=phix[3*ii+k]-acyrl[k][0]*pleft;
			phix[3*ii+k]=-phix[3*ii+k]/DX[ix];
		}
	}

	//E_y
	for(ix=1;ix<=NX;ix++)
	{
		for(jy=2;jy<=NYREAL;jy++)
		{
			ii = (ix-1)*NX + jy-1;
			ii_l = (ix-1)*NX + jy-2;
			/* jy=1 case taken care of next with Neumann BC!
			if(jy > 1) // May be necessary when jy = 1 as then ii_l < 0
			{
				ii_l = (ix-1)*NX + jy-2;
			}
			else
			{
				ii_l = ii;
			}
			*/
			for(k=0;k<=mpto;k++)
			{
				phiy[3*ii+k]=0.;
				for(kk=0;kk<=mpto;kk++)
				{
					phiy[3*ii+k]=phiy[3*ii+k]+(acvy[k][kk]+acxrr[k][kk])*POTC[3*ii+kk]-acxrl[k][kk]*POTC[3*ii_l+kk];
				}
				phiy[3*ii+k]=-phiy[3*ii+k]/DY[jy];
			}
		}
	}

	for(ix=1;ix<=NX;ix++)
	{
		//at bottom, neumann
		jy=1;
		ii = (ix-1)*NX + jy-1;
		for(k=0;k<=mpto;k++)
		{
			phiy[3*ii+k]=0.;
		}
		phiy[3*ii]=-POTC[3*ii+2]/DY[jy];
	}

	//top
	for(ix=1;ix<=NX;ix++)
	{
		//at top, neumann
		jy=NYREAL;
		ii = (ix-1)*NX + jy-1;
		for(k=0;k<=mpto;k++)
		{
			phiy[3*ii+k]=0.;
		}
		phiy[3*ii]=-POTC[3*ii+2]/DY[jy];	// There's a chance this not be negative...?

		/*
		jy=NYREAL;
		ii = (ix-1)*NX + jy-1;
		ii_l = (ix-1)*NX + jy-2;,,m,
		if(CX[ix]>=xc&&CX[ix]<=xd)
		{
			for(k=0;k<=mpto;k++)
			{
				phiy[3*ii+k]=0.;
				for(kk=0;kk<=mpto;kk++)
				{
					j = (mpto+1)*ii + kk;
					j_l = (mpto+1)*ii_l + kk;
					phiy[3*ii+k]=phiy[3*ii+k]+acvy[k][kk]*POTC[3*ii+kk]-acxrl[k][kk]*POTC[3*ii_l+kk];
				}
				phiy[3*ii+k]=phiy[3*ii+k]+acxrr[k][0]*pybc;
				phiy[3*ii+k]=-phiy[3*ii+k]/DY[jy];
			}
		}
		*/
	}

	for(ix=1;ix<=NX;ix++)
	{
		for(jy=1;jy<=NYREAL;jy++)
		{
			ii = (ix-1)*NX + jy-1;
			for(k=0;k<=mpto;k++)
			{
				phix[3*ii+k]=PVE*phix[3*ii+k];
				phiy[3*ii+k]=PVE*phiy[3*ii+k];
			}
		}
	}

	// evaluate the integrals related to the electric field
	// CASE 1D

	for(i=0; i<NX; ++i)
	{
		for(j=0; j<NY; ++j)
		{
			ii = i*NX + j;
			xb = CX[i+1] + DX[NX]/2.;

			/*
			IE_X[i][j] = phix[0][i][j];
			IXE_X[i][j] = phix[1][i][j]/6.;
			IYE_X[i][j] = phix[2][i][j]/6.;
			IXXE_X[i][j] =phix[0][i][j]/3.;
			IXYE_X[i][j] = 0.;
			IYYE_X[i][j] = phix[0][i][j]/3.;

			IE_Y[i][j] = phiy[0][i][j];
			IXE_Y[i][j] = phiy[1][i][j]/6.;
			IYE_Y[i][j] = phiy[2][i][j]/6.;
			IXXE_Y[i][j] = phiy[0][i][j]/3.;
			IXYE_Y[i][j] = 0.;
			IYYE_Y[i][j] = phiy[0][i][j]/3.;
			*/

			// WILL NEED TO CHANGE THIS ALL TO dx*dx* FOR 2D (or dx*dy*)
			IE_X[ii] = dx*dx*phix[3*ii];
			IXE_X[ii] = dx*dx*phix[3*ii+1]/12.;	// This one is giving a substantial error
			IYE_X[ii] = dx*dx*phix[3*ii+2]/12.;
			IXXE_X[ii] =dx*dx*phix[3*ii]/12.;
			IXYE_X[ii] = 0.;
			IYYE_X[ii] = dx*dx*phix[3*ii]/12.;

			IE_Y[ii] = dx*dx*phiy[3*ii];
			IXE_Y[ii] = dx*dx*phiy[3*ii+1]/12.;
			IYE_Y[ii] = dx*dx*phiy[3*ii+2]/12.;
			IXXE_Y[ii] = dx*dx*phiy[3*ii]/12.;
			IXYE_Y[ii] = 0.;
			IYYE_Y[ii] = dx*dx*phiy[3*ii]/12.;

			// the sign of the components of the electric field in the cells
			if(IE_X[ii]>0.)
			{
				SE_X[ii] = 1.;
			}
			else
			{
				SE_X[ii] = -1.;
			}
			if(IE_Y[ii]>0.)
			{
				SE_Y[ii] = 1.;
			}
			else
			{
				SE_Y[ii] = -1.;
			}
		}
	}
}

void InitPOT(void)
{
	for(int ix=1;ix<=NX;ix++)
	{
		for(int iy=1;iy<=NY;iy++)
		{
			POT[ix][iy] = 0;
		}
	}
}
