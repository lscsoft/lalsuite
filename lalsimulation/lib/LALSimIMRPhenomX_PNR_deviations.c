
/*
* Copyright (C) 2021 The University of Amsterdam
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
*  MA  02110-1301  USA
*/

#ifdef __cplusplus
extern "C" {
#endif

/**
* \author Lionel London
*/

//
#include "LALSimIMRPhenomX_PNR_deviations.h"

#ifndef _OPENMP
#define omp ignore
#endif

#ifndef PHENOMXHMDEBUG
#define DEBUG 0
#else
#define DEBUG 1
#endif


// Import usefuls
#include <math.h>



// MU1 fit implementation
// Header formatting for MU1 of (l,m)=(2,2) multipole
double XLALSimIMRPhenomXCP_MU1_l2m2( double theta, double eta, double a1 ){

	/*
	Hola, soy un codigo escribido por "4b_document_fits.py". Dat script is geschreven door een mens die niet kan worden genoemd.
	*/

	// // Flatten in eta
	// eta = ( eta >= 0.09876 ) ? eta : 0.09876;
	// a1  = ( a1  <= 0.8     ) ? a1  : 0.8;
	// a1  = ( a1  >= 0.2     ) ? a1  : 0.2;

	// Preliminaries
	double u = cos(theta);
	double u2 = u*u;
	double u3 = u2*u;
	double a12 = a1*a1;
	double eta2 = eta*eta;
	double eta3 = eta2*eta;
	double mu1;

	// Evaluate fit for this parameter
	mu1 = 3.10731588e+00*a1 + -2.98890174e-01 + 2.83580591e+00*eta + -2.69235962e-01*u + 4.58765298e+00*(u*eta) + -4.58961732e+01*(eta*a1) + -5.27239994e+00*a12 + 2.35767086e+00*(u*a1) + -1.66981400e+00*(u*u) + -6.65783904e+00*eta2 + 2.28592256e+02*(eta2*a1) + -2.85191772e+00*(u*a12) + 8.53735433e+01*(eta*a12) + -1.69564326e+01*(u*eta2) + 3.80240537e+01*(u2*eta) + -5.90965234e+01*(u*eta*a1) + -3.25953234e-01*(u2*u) + -3.71962074e+02*(eta3*a1) + 3.65316048e+02*(u*eta2*a1) + 3.69910553e+00*(u2*a12) + 2.43498635e+00*(u3*a1) + -4.48771815e+02*(eta2*a12) + -3.06134124e+01*(u2*eta*a1) + -2.49773688e+02*(u2*eta2) + 7.77424585e+01*(u*eta*a12) + -3.86558965e+01*(u2*eta*a12) + 7.62503928e+02*(eta3*a12) + 5.19092710e+02*(u2*eta2*eta) + -5.04214664e+02*(u*eta2*a12) + -6.11685145e+02*(u*eta3*a1) + 2.91303571e+02*(u2*eta2*a1) + -6.71291161e+00*(u3*eta*a1) + -3.34353185e+00*(u3*a12) + 4.49179407e+01*(u3*eta2*eta) + 1.10210016e+02*(u2*eta2*a12) + 9.08677395e+02*(u*eta3*a12) + -3.78457492e+01*(u3*eta2*a1) + -7.50921355e+02*(u2*eta3*a1) + 1.89942289e+01*(u3*eta*a12);

	// Return answer
	return mu1;

} // END of MU1 (2,2) fit implementation



// MU2 fit implementation
// Header formatting for MU2 of (l,m)=(2,2) multipole
double XLALSimIMRPhenomXCP_MU2_l2m2( double theta, double eta, double a1 ){

	/*
	Hola, soy un codigo escribido por "4b_document_fits.py". Dat script is geschreven door een mens die niet kan worden genoemd.
	*/

        // // Flatten in eta
	// eta = ( eta >= 0.09876 ) ? eta : 0.09876;
	// a1  = ( a1  <= 0.8     ) ? a1  : 0.8;
	// a1  = ( a1  >= 0.2     ) ? a1  : 0.2;

	// Preliminaries
	double u = cos(theta);
	double u2 = u*u;
	double u3 = u2*u;
	double a12 = a1*a1;
	double eta2 = eta*eta;
	double eta3 = eta2*eta;
	double mu2;

	// Evaluate fit for this parameter
	mu2 = -2.02502540e-01 + 8.72364229e+00*(eta*a1) + -1.58131236e+00*a12 + -3.87863737e+00*(u*u) + -5.05936693e+01*(eta2*a1) + -8.00671915e+00*(u*a12) + 9.12534447e+00*(u2*a1) + 5.33069685e+00*(eta*a12) + 8.51436158e+00*(eta2*eta) + 5.78604124e+01*(u2*eta) + 6.55310845e+00*(u*eta*a1) + -1.94647765e+00*(u2*u) + 4.08691562e+01*(u3*eta) + 7.52332356e+01*(eta3*a1) + -1.06817659e+01*(u2*a12) + -2.87950209e+00*(u3*a1) + -6.43240440e+01*(u2*eta*a1) + -2.49087968e+01*(u*eta2*eta) + -3.16185296e+02*(u2*eta2) + 1.06442533e+02*(u*eta*a12) + 1.03449590e+02*(u2*eta*a12) + 6.51115537e+02*(u2*eta2*eta) + -5.79811830e+02*(u*eta2*a12) + -2.76924732e+02*(u3*eta2) + 2.10132916e+01*(u3*eta*a1) + 2.70715432e+00*(u3*a12) + 6.47099541e+02*(u3*eta2*eta) + -2.05802266e+02*(u2*eta2*a12) + 1.03222489e+03*(u*eta3*a12) + -9.28493309e+01*(u3*eta2*a1) + 2.54661923e+02*(u2*eta3*a1);

	// Return answer
	return mu2;

} // END of MU2 (2,2) fit implementation



// MU3 fit implementation
// Header formatting for MU3 of (l,m)=(2,2) multipole
double XLALSimIMRPhenomXCP_MU3_l2m2( double theta, double eta, double a1 ){

	/*
	Hola, soy un codigo escribido por "4b_document_fits.py". Dat script is geschreven door een mens die niet kan worden genoemd.
	*/

       	// // Flatten in eta
	// eta = ( eta >= 0.09876 ) ? eta : 0.09876;
	// a1  = ( a1  <= 0.8     ) ? a1  : 0.8;
	// a1  = ( a1  >= 0.2     ) ? a1  : 0.2;

	// Preliminaries
	double u = cos(theta);
	double u2 = u*u;
	double u3 = u2*u;
	double a12 = a1*a1;
	double eta2 = eta*eta;
	double eta3 = eta2*eta;
	double mu3;

	// Evaluate fit for this parameter
	mu3 = 6.43824401e-02*a1 + -3.17858924e-03 + 6.54067059e-03*eta + 2.48066282e-02*u + -1.50482854e-01*(u*eta) + -8.51847106e-01*(eta*a1) + -1.14919973e-01*(u*a1) + -1.23872940e-02*(u*u) + 4.35501780e+00*(eta2*a1) + 1.51718799e-01*(u*a12) + -4.01788734e-02*(u2*a1) + 1.94815013e-01*(u*eta2) + 1.61029465e-01*(u2*eta) + 5.63051639e-01*(u*eta*a1) + -9.08508144e-03*(u2*u) + -7.78742176e+00*(eta3*a1) + 7.84548415e-02*(u3*a1) + 7.87294341e-01*(u2*eta*a1) + -9.60114345e-01*(u2*eta2) + -9.10606320e-01*(u*eta*a12) + -1.42496492e-01*(u2*eta*a12) + 1.67682276e+00*(u2*eta2*eta) + 1.14108488e+00*(u*eta2*a12) + -1.45121681e+00*(u*eta3*a1) + -3.90246382e+00*(u2*eta2*a1) + -1.88632810e-01*(u3*eta*a1) + -1.20970655e-01*(u3*a12) + 5.44663689e-01*(u3*eta2*eta) + 3.60304651e-01*(u2*eta2*a12) + -5.29498004e-01*(u3*eta2*a1) + 6.67574591e+00*(u2*eta3*a1) + 5.00625025e-01*(u3*eta*a12);

	// Return answer
	return mu3;

} // END of MU3 (2,2) fit implementation



// NU4 fit implementation
// Header formatting for NU4 of (l,m)=(2,2) multipole
double XLALSimIMRPhenomXCP_NU4_l2m2( double theta, double eta, double a1 ){

	/*
	Hola, soy un codigo escribido por "4b_document_fits.py". Dat script is geschreven door een mens die niet kan worden genoemd.
	*/

        //// Flatten in eta
	// eta = ( eta >= 0.09876 ) ? eta : 0.09876;
	// a1  = ( a1  <= 0.8     ) ? a1  : 0.8;
	// a1  = ( a1  >= 0.2     ) ? a1  : 0.2;

	// Preliminaries
	double u = cos(theta);
	double u2 = u*u;
	double u3 = u2*u;
	double a12 = a1*a1;
	double eta2 = eta*eta;
	double eta3 = eta2*eta;
	double nu4;

	// Evaluate fit for this parameter
	nu4 = -3.16028942e-02*a1 + 4.49481835e-03 + -4.67316682e-02*eta + 1.59179379e-02*u + -3.27118716e-01*(u*eta) + 3.77496559e-01*(eta*a1) + 3.25729551e-02*a12 + -7.76768528e-02*(u*a1) + -7.33656551e-03*(u*u) + 1.08635205e-01*eta2 + -1.32662456e+00*(eta2*a1) + 4.84795811e-02*(u*a12) + 7.64861112e-02*(u2*a1) + -3.45847253e-01*(eta*a12) + 1.90217409e+00*(u*eta2) + 1.74016014e+00*(u*eta*a1) + 1.38785202e+00*(eta3*a1) + -1.06852751e+01*(u*eta2*a1) + -8.33709907e-02*(u2*a12) + -7.93716551e-03*(u3*a1) + 8.34552086e-01*(eta2*a12) + -7.44879716e-01*(u2*eta*a1) + -3.35542027e+00*(u*eta2*eta) + 5.13747474e-01*(u2*eta2) + -1.42151090e+00*(u*eta*a12) + 8.28154534e-01*(u2*eta*a12) + -1.67075105e+00*(u2*eta2*eta) + 9.46433414e+00*(u*eta2*a12) + 1.97711659e+01*(u*eta3*a1) + 1.80163775e+00*(u2*eta2*a1) + 2.22117601e-01*(u3*eta2) + -9.13086579e-02*(u3*eta*a1) + 3.41917884e-02*(u3*a12) + -1.01537065e+00*(u3*eta2*eta) + -2.00129470e+00*(u2*eta2*a12) + -1.81767526e+01*(u*eta3*a12) + 5.87741796e-01*(u3*eta2*a1) + -1.52767109e-01*(u3*eta*a12);

	// Return answer
	return nu4;

} // END of NU4 (2,2) fit implementation



// NU5 fit implementation
// Header formatting for NU5 of (l,m)=(2,2) multipole
double XLALSimIMRPhenomXCP_NU5_l2m2( double theta, double eta, double a1 ){

	/*
	Hola, soy un codigo escribido por "4b_document_fits.py". Dat script is geschreven door een mens die niet kan worden genoemd.
	*/

   	// // Flatten in eta
	// eta = ( eta >= 0.09876 ) ? eta : 0.09876;
	// a1  = ( a1  <= 0.8     ) ? a1  : 0.8;
	// a1  = ( a1  >= 0.2     ) ? a1  : 0.2;

	// Preliminaries
	double u = cos(theta);
	double u2 = u*u;
	double u3 = u2*u;
	double a12 = a1*a1;
	double eta2 = eta*eta;
	double eta3 = eta2*eta;
	double nu5;

	// Evaluate fit for this parameter
	nu5 = -8.24852329e-03 + 1.29541487e-01*eta + -2.68106259e-02*u + 3.96908249e-01*(u*eta) + -4.57402308e-01*(eta*a1) + -6.96270582e-02*a12 + 8.87404766e-02*(u*a1) + -3.95509060e-02*(u*u) + -7.66209616e-01*eta2 + 3.52527019e+00*(eta2*a1) + -8.52495198e-02*(u*a12) + 1.28777123e-01*(u2*a1) + 1.34613777e+00*(eta*a12) + -1.97688956e+00*(u*eta2) + 1.33051669e+00*(eta2*eta) + 6.66123113e-01*(u2*eta) + -8.83292400e-01*(u*eta*a1) + 4.13508242e-03*(u2*u) + -6.34717777e+00*(eta3*a1) + 3.39853108e+00*(u*eta2*a1) + -5.33540972e-02*(u2*a12) + -3.96980662e-02*(u3*a1) + -8.95748509e+00*(eta2*a12) + -1.67674467e+00*(u2*eta*a1) + 3.19730879e+00*(u*eta2*eta) + -4.21512938e+00*(u2*eta2) + 3.45438117e-01*(u*eta*a12) + 2.93065248e-01*(u2*eta*a12) + 1.82938901e+01*(eta3*a12) + 9.11903434e+00*(u2*eta2*eta) + -5.03020035e+00*(u*eta3*a1) + 9.57710907e+00*(u2*eta2*a1) + -6.00764557e-01*(u3*eta2) + 4.25544610e-01*(u3*eta*a1) + 4.02874363e-02*(u3*a12) + 2.55062055e+00*(u3*eta2*eta) + -1.31746981e+00*(u3*eta2*a1) + -2.11785354e+01*(u2*eta3*a1) + -1.29990564e-01*(u3*eta*a12);

	// Return answer
	return nu5;

} // END of NU5 (2,2) fit implementation



// NU6 fit implementation
// Header formatting for NU6 of (l,m)=(2,2) multipole
double XLALSimIMRPhenomXCP_NU6_l2m2( double theta, double eta, double a1 ){

	/*
	Hola, soy un codigo escribido por "4b_document_fits.py". Dat script is geschreven door een mens die niet kan worden genoemd.
	*/

	// // Flatten in eta
	// eta = ( eta >= 0.09876 ) ? eta : 0.09876;
	// a1  = ( a1  <= 0.8     ) ? a1  : 0.8;
	// a1  = ( a1  >= 0.2     ) ? a1  : 0.2;

	// Preliminaries
	double u = cos(theta);
	double u2 = u*u;
	double u3 = u2*u;
	double a12 = a1*a1;
	double eta2 = eta*eta;
	double eta3 = eta2*eta;
	double nu6;

	// Evaluate fit for this parameter
	nu6 = -8.14163374e-02*a1 + 1.60753917e-02 + -3.14932520e-01*eta + 3.69589365e-02*u + -5.43366127e-01*(u*eta) + 1.65256072e+00*(eta*a1) + 1.57642062e-01*a12 + -3.53871924e-01*(u*a1) + 1.91773389e+00*eta2 + -1.01369676e+01*(eta2*a1) + 4.49846052e-01*(u*a12) + -3.23664263e-02*(u2*a1) + -2.77552262e+00*(eta*a12) + 3.18104872e+00*(u*eta2) + -3.81823248e+00*(eta2*eta) + -4.05623855e-02*(u2*eta) + 5.78266831e+00*(u*eta*a1) + 4.67460809e-02*(u2*u) + -9.03957684e-01*(u3*eta) + 1.97689861e+01*(eta3*a1) + -3.29265427e+01*(u*eta2*a1) + 4.68080513e-02*(u2*a12) + 1.57596950e+01*(eta2*a12) + 4.06667657e-01*(u2*eta*a1) + -6.37440535e+00*(u*eta2*eta) + -7.23117535e+00*(u*eta*a12) + -6.79662526e-01*(u2*eta*a12) + -2.89973575e+01*(eta3*a12) + 2.92509266e-01*(u2*eta2*eta) + 3.96497438e+01*(u*eta2*a12) + 6.16173322e+01*(u*eta3*a1) + 4.62060734e+00*(u3*eta2) + 5.63722799e-01*(u3*eta*a1) + -7.10443312e-02*(u3*a12) + -7.10237285e+00*(u3*eta2*eta) + 1.81706494e+00*(u2*eta2*a12) + -7.13125243e+01*(u*eta3*a12) + -2.12229389e+00*(u3*eta2*a1) + -3.62928100e+00*(u2*eta3*a1) + 2.39385816e-01*(u3*eta*a12);

	// Return answer
	return nu6;

} // END of NU6 (2,2) fit implementation



// ZETA1 fit implementation
// Header formatting for ZETA1 of (l,m)=(2,2) multipole
double XLALSimIMRPhenomXCP_ZETA1_l2m2( double theta, double eta, double a1 ){

	/*
	Hola, soy un codigo escribido por "4b_document_fits.py". Dat script is geschreven door een mens die niet kan worden genoemd.
	*/

	// // Flatten in eta
	// eta = ( eta >= 0.09876 ) ? eta : 0.09876;
	// a1  = ( a1  <= 0.8     ) ? a1  : 0.8;
	// a1  = ( a1  >= 0.2     ) ? a1  : 0.2;

	// Preliminaries
	double u = cos(theta);
	double u2 = u*u;
	double u3 = u2*u;
	double a12 = a1*a1;
	double eta2 = eta*eta;
	double eta3 = eta2*eta;
	double zeta1;

	// Evaluate fit for this parameter
	zeta1 = 9.14075603e-05*a1 + -6.64596933e-06 + 7.26570370e-05*eta + -2.65487354e-05*u + 6.13445304e-04*(u*eta) + -1.13031091e-03*(eta*a1) + -7.04782574e-05*a12 + 2.78657414e-04*(u*a1) + 2.93046524e-05*(u*u) + -1.62592535e-04*eta2 + 4.77274706e-03*(eta2*a1) + -3.28088748e-04*(u*a12) + -2.40358466e-04*(u2*a1) + 6.89104675e-04*(eta*a12) + -4.33992656e-03*(u*eta2) + -1.71570550e-04*(u2*eta) + -5.75897400e-03*(u*eta*a1) + -4.06577752e-05*(u2*u) + 7.19295439e-04*(u3*eta) + -7.10577068e-03*(eta3*a1) + 3.68151682e-02*(u*eta2*a1) + 1.70628764e-04*(u2*a12) + 2.97287666e-05*(u3*a1) + -1.58004432e-03*(eta2*a12) + 2.69942014e-03*(u2*eta*a1) + 9.27082591e-03*(u*eta2*eta) + 6.89186183e-03*(u*eta*a12) + -1.33029865e-03*(u2*eta*a12) + 1.14545659e-03*(u2*eta2*eta) + -4.28423981e-02*(u*eta2*a12) + -7.28752610e-02*(u*eta3*a1) + -1.14024958e-02*(u2*eta2*a1) + -3.59470348e-03*(u3*eta2) + -6.18248914e-04*(u3*eta*a1) + 5.76871905e-03*(u3*eta2*eta) + 2.67905006e-03*(u2*eta2*a12) + 8.18766554e-02*(u*eta3*a12) + 1.68471071e-03*(u3*eta2*a1) + 1.73091684e-02*(u2*eta3*a1) + 8.47582112e-05*(u3*eta*a12);

	// Return answer
	return zeta1;

} // END of ZETA1 (2,2) fit implementation



// ZETA2 fit implementation
// Header formatting for ZETA2 of (l,m)=(2,2) multipole
double XLALSimIMRPhenomXCP_ZETA2_l2m2( double theta, double eta, double a1 ){

	/*
	Hola, soy un codigo escribido por "4b_document_fits.py". Dat script is geschreven door een mens die niet kan worden genoemd.
	*/

	// // Flatten in eta
	// eta = ( eta >= 0.09876 ) ? eta : 0.09876;
	// a1  = ( a1  <= 0.8     ) ? a1  : 0.8;
	// a1  = ( a1  >= 0.2     ) ? a1  : 0.2;

	// Preliminaries
	double u = cos(theta);
	double u2 = u*u;
	double u3 = u2*u;
	double a12 = a1*a1;
	double eta2 = eta*eta;
	double eta3 = eta2*eta;
	double zeta2;

	// Evaluate fit for this parameter
	zeta2 = -1.93349374e+01*a1 + 1.21917535e+00 + 5.64084304e+00*u + -1.35818573e+02*(u*eta) + 2.03450770e+02*(eta*a1) + 1.91050318e+01*a12 + -3.45695270e+01*(u*a1) + -4.07600123e+00*(u*u) + -7.13247650e+01*eta2 + -6.59883962e+02*(eta2*a1) + 3.69257896e+01*(u*a12) + 4.19166290e+01*(u2*a1) + -1.83835664e+02*(eta*a12) + 8.83374450e+02*(u*eta2) + 2.33819330e+02*(eta2*eta) + 8.18248325e+02*(u*eta*a1) + 3.34325806e+00*(u2*u) + -1.82452234e+01*(u3*eta) + 5.76375582e+02*(eta3*a1) + -5.26052382e+03*(u*eta2*a1) + -3.57570544e+01*(u2*a12) + -1.54037927e+01*(u3*a1) + 4.36853250e+02*(eta2*a12) + -4.35345536e+02*(u2*eta*a1) + -1.69758341e+03*(u*eta2*eta) + 2.99182150e+02*(u2*eta2) + -8.86211351e+02*(u*eta*a12) + 3.47168913e+02*(u2*eta*a12) + -9.23263240e+02*(u2*eta2*eta) + 5.67490957e+03*(u*eta2*a12) + 1.00342283e+04*(u*eta3*a1) + 1.27778712e+03*(u2*eta2*a1) + 9.20983561e+01*(u3*eta*a1) + 1.53555516e+01*(u3*a12) + -8.06783834e+02*(u2*eta2*a12) + -1.07373527e+04*(u*eta3*a12) + -3.64839255e+01*(u3*eta2*a1) + -8.79147687e+02*(u2*eta3*a1) + -8.26009443e+01*(u3*eta*a12);

	// Return answer
	return zeta2;

} // END of ZETA2 (2,2) fit implementation



// NU0 fit implementation
// Header formatting for NU0 of (l,m)=(2,2) multipole
double XLALSimIMRPhenomXCP_NU0_l2m2( double theta, double eta, double a1 ){

	/*
	Hola, soy un codigo escribido por "4b_document_fits.py". Dat script is geschreven door een mens die niet kan worden genoemd.
	*/

	// // Flatten in eta
	// eta = ( eta >= 0.09876 ) ? eta : 0.09876;
	// a1  = ( a1  <= 0.8     ) ? a1  : 0.8;
	// a1  = ( a1  >= 0.2     ) ? a1  : 0.2;

	// Preliminaries
	double u = cos(theta);
	double u2 = u*u;
	double u3 = u2*u;
	double a12 = a1*a1;
	double eta2 = eta*eta;
	double eta3 = eta2*eta;
	double nu0;

	// Evaluate fit for this parameter
	nu0 = 9.69656479e+02*a1 + -1.40050910e+02 + 1.28604150e+03*eta + 1.21045131e+02*u + -8.58972371e+03*(eta*a1) + -3.80191737e+02*a12 + -2.62095295e+03*(u*a1) + 3.75573898e+02*(u*u) + -3.02176869e+03*eta2 + 1.91665811e+04*(eta2*a1) + 3.76489144e+03*(u*a12) + -2.80301545e+03*(u2*a1) + -3.28538542e+03*(u2*eta) + 3.46510614e+04*(u*eta*a1) + 3.56880662e+02*(u2*u) + -8.93916723e+03*(u3*eta) + -1.76931761e+05*(u*eta2*a1) + 2.34502372e+03*(u2*a12) + 1.13064791e+03*(u3*a1) + 2.54661712e+04*(eta2*a12) + 2.98001687e+04*(u2*eta*a1) + -1.27661553e+04*(u*eta2*eta) + -5.40526831e+04*(u*eta*a12) + -2.41699850e+04*(u2*eta*a12) + -7.73072478e+04*(eta3*a12) + 2.65990668e+04*(u2*eta2*eta) + 2.78587059e+05*(u*eta2*a12) + 3.35360241e+05*(u*eta3*a1) + -7.36180264e+04*(u2*eta2*a1) + 4.96525314e+04*(u3*eta2) + -3.75871238e+03*(u3*eta*a1) + -1.09191647e+03*(u3*a12) + -7.11798894e+04*(u3*eta2*eta) + 5.84802389e+04*(u2*eta2*a12) + -4.99111316e+05*(u*eta3*a12) + -7.99273344e+03*(u3*eta2*a1) + 5.11072226e+03*(u3*eta*a12);

	// Return answer
	return nu0;

} // END of NU0 (2,2) fit implementation

// MU1 fit implementation
// Header formatting for MU1 of (l,m)=(3,3) multipole
double XLALSimIMRPhenomXCP_MU1_l3m3( double theta, double eta, double a1 ){

	/*
	Hola, soy un codigo escribido por "4b_document_fits.py". Dat script is geschreven door een mens die niet kan worden genoemd.
	*/

	// Preliminaries
	double u = cos(theta);
	double u2 = u*u;
	double u3 = u2*u;
	UNUSED double u4 = u3*u;
	double a12 = a1*a1;
	double a13 = a12*a1;
	double eta2 = eta*eta;
	double mu1;

	// Evaluate fit for this parameter
	mu1 = -3.84796066e+00*a1 + -7.05817506e-02 + 1.21568676e-01*u + -3.32955057e+00*(u*a1) + 1.45279803e+01*a12 + 6.02824893e+00*(u*u) + 6.34966196e+01*(eta*a1) + 1.09295567e+01*(u*a12) + 3.79205216e+01*(u*eta*a1) + -8.59816637e+01*(u2*eta) + -2.95015094e+00*(u2*u) + -1.79406237e+01*(u2*a1) + 8.10569148e+00*(u*eta2) + -1.28362313e+01*(a12*a1) + -2.24119230e+02*(eta*a12) + -2.31493382e+02*(eta2*a1) + 8.01722640e+02*(eta2*a12) + 2.67574200e+01*(u3*a1) + -9.19576147e+00*(u3*u) + 1.94511999e+02*(eta*a13) + -8.92178321e+00*(u*a13) + 4.08796320e+01*(u3*eta) + 2.77108126e+02*(u2*eta2) + 2.33691206e+02*(u2*eta*a1) + -1.90906175e+02*(u*eta2*a1) + -1.42896519e+02*(u*eta*a12) + 6.59719120e+01*(u2*eta*a12) + 6.31279593e+02*(u*eta2*a12) + -6.03516802e+02*(u2*eta2*a1) + -6.35584153e+01*(u3*a12) + -3.76169532e+02*(u3*eta*a1) + -6.88708416e+02*(eta2*a13) + -1.48924737e+02*(u3*eta2) + 2.01218446e+01*(u2*a13) + 1.23850850e+02*(u*eta*a13) + 3.97214761e+01*(u4*a1) + 1.29916578e+02*(u4*eta) + 1.35427866e+03*(u3*eta2*a1) + -4.95972583e+01*(u4*a12) + -3.42665931e+02*(u2*eta*a13) + -4.08179558e+02*(u4*eta2) + -6.89874345e+02*(u2*eta2*a12) + -5.38995076e+02*(u4*eta*a1) + 9.11616099e+02*(u3*eta*a12) + -5.28188545e+02*(u*eta2*a13) + 4.45398586e+01*(u3*a13) + 1.43408019e+01*(u4*a13) + 1.50904278e+03*(u4*eta2*a1) + -6.51036477e+02*(u3*eta*a13) + 6.36604875e+02*(u4*eta*a12) + 1.49605186e+03*(u2*eta2*a13) + -3.29736033e+03*(u3*eta2*a12) + -1.43387906e+03*(u4*eta2*a12) + 2.37658315e+03*(u3*eta2*a13) + -1.47497570e+02*(u4*eta*a13);

	// Return answer
	return mu1;

} // END of MU1 (3,3) fit implementation



// MU2 fit implementation
// Header formatting for MU2 of (l,m)=(3,3) multipole
double XLALSimIMRPhenomXCP_MU2_l3m3( double theta, double eta, double a1 ){

	/*
	Hola, soy un codigo escribido por "4b_document_fits.py". Dat script is geschreven door een mens die niet kan worden genoemd.
	*/

	// Preliminaries
	double u = cos(theta);
	double u2 = u*u;
	double u3 = u2*u;
	UNUSED double u4 = u3*u;
	double a12 = a1*a1;
	double a13 = a12*a1;
	double eta2 = eta*eta;
	double mu2;

	// Evaluate fit for this parameter
	mu2 = 7.47818099e+00*a1 + -7.03752800e-01 + 1.18890118e+01*eta + -3.39315179e+00*u + 2.99494916e+01*(u*a1) + -1.70517129e+01*a12 + -1.15606124e+00*(u*u) + 3.70333685e+01*(u*eta) + -1.21612496e+02*(eta*a1) + -4.46727469e+01*eta2 + -6.95692459e+01*(u*a12) + -3.30667559e+02*(u*eta*a1) + 2.35627361e+00*(u2*eta) + 4.24956067e+00*(u2*u) + -9.80880323e+00*(u2*a1) + -8.27439034e+01*(u*eta2) + 1.16949376e+01*(a12*a1) + 2.81765277e+02*(eta*a12) + 4.37554921e+02*(eta2*a1) + -1.01344194e+03*(eta2*a12) + -3.98776713e+01*(u3*a1) + 2.88708815e+00*(u3*u) + -1.95334810e+02*(eta*a13) + 4.67635563e+01*(u*a13) + -4.11686620e+01*(u3*eta) + 2.26219351e+02*(u2*eta*a1) + 2.86543174e+01*(u2*a12) + 7.77340880e+02*(u*eta2*a1) + 7.61635064e+02*(u*eta*a12) + -5.63600753e+02*(u2*eta*a12) + -1.78703069e+03*(u*eta2*a12) + -7.12923184e+02*(u2*eta2*a1) + 9.72543354e+01*(u3*a12) + 4.14287030e+02*(u3*eta*a1) + -2.86733254e+01*(u4*eta) + 7.04420773e+02*(eta2*a13) + 6.52575573e+01*(u3*eta2) + -1.69940614e+01*(u2*a13) + -5.02734853e+02*(u*eta*a13) + -8.32035149e+02*(u3*eta2*a1) + -9.15864076e+00*(u4*a12) + 3.26091951e+02*(u2*eta*a13) + 8.60925038e+01*(u4*eta2) + 1.63112868e+03*(u2*eta2*a12) + -7.38065591e+01*(u4*eta*a1) + -1.03595569e+03*(u3*eta*a12) + 1.15307817e+03*(u*eta2*a13) + -6.80286416e+01*(u3*a13) + 3.88453028e+00*(u4*a13) + 1.84299386e+02*(u4*eta2*a1) + 7.28708130e+02*(u3*eta*a13) + 2.48344476e+02*(u4*eta*a12) + -8.26802223e+02*(u2*eta2*a13) + 2.22021515e+03*(u3*eta2*a12) + -4.77957901e+02*(u4*eta2*a12) + -1.58963329e+03*(u3*eta2*a13) + -1.08030307e+02*(u4*eta*a13);

	// Return answer
	return mu2;

} // END of MU2 (3,3) fit implementation



// MU3 fit implementation
// Header formatting for MU3 of (l,m)=(3,3) multipole
double XLALSimIMRPhenomXCP_MU3_l3m3( double theta, double eta, double a1 ){

	/*
	Hola, soy un codigo escribido por "4b_document_fits.py". Dat script is geschreven door een mens die niet kan worden genoemd.
	*/

	// Preliminaries
	double u = cos(theta);
	double u2 = u*u;
	double u3 = u2*u;
	UNUSED double u4 = u3*u;
	double a12 = a1*a1;
	double a13 = a12*a1;
	double eta2 = eta*eta;
	double mu3;

	// Evaluate fit for this parameter
	mu3 = -2.92826085e-01*a1 + 4.00831355e-02 + -6.78453540e-01*eta + 1.69841734e-02*u + -2.57312944e-01*(u*a1) + 6.39784806e-01*a12 + -2.31536877e-01*(u*u) + -2.41888246e-01*(u*eta) + 4.95548971e+00*(eta*a1) + 2.37429504e+00*eta2 + 9.89479593e-01*(u*a12) + 3.80025000e+00*(u*eta*a1) + 3.42384733e+00*(u2*eta) + 1.32362070e+00*(u2*a1) + 2.03296485e+00*(u*eta2) + -4.51848038e-01*(a12*a1) + -1.05977941e+01*(eta*a12) + -1.79880007e+01*(eta2*a1) + 3.88585826e+01*(eta2*a12) + 1.08419742e-01*(u3*a1) + 2.95605100e-01*(u3*u) + 7.24579590e+00*(eta*a13) + -1.00632694e+00*(u*a13) + -9.89315032e+00*(u2*eta2) + -2.02592474e+01*(u2*eta*a1) + -1.89650586e+00*(u2*a12) + -2.27713382e+01*(u*eta2*a1) + -1.49608491e+01*(u*eta*a12) + 2.88394384e+01*(u2*eta*a12) + 7.50486750e+01*(u*eta2*a12) + 5.76958608e+01*(u2*eta2*a1) + -7.14406060e-01*(u3*a12) + -1.50879383e+00*(u3*eta*a1) + -4.36702529e+00*(u4*eta) + -2.65052614e+01*(eta2*a13) + -1.92813117e+00*(u3*eta2) + 4.59142180e-01*(u2*a13) + 1.55287853e+01*(u*eta*a13) + -1.77690774e+00*(u4*a1) + 2.08333492e+01*(u3*eta2*a1) + 2.65579384e+00*(u4*a12) + -6.22168581e+00*(u2*eta*a13) + 1.20197689e+01*(u4*eta2) + -7.14793444e+01*(u2*eta2*a12) + 2.70295932e+01*(u4*eta*a1) + 1.06042548e+01*(u3*eta*a12) + -7.13323842e+01*(u*eta2*a13) + 8.70257518e-01*(u3*a13) + -7.10610270e-01*(u4*a13) + -7.26302148e+01*(u4*eta2*a1) + -1.33863275e+01*(u3*eta*a13) + -4.04045931e+01*(u4*eta*a12) + -7.40725231e+01*(u3*eta2*a12) + 9.21784618e+01*(u4*eta2*a12) + 7.41760537e+01*(u3*eta2*a13) + 1.01692139e+01*(u4*eta*a13);

	// Return answer
	return mu3;

} // END of MU3 (3,3) fit implementation



// MU4 fit implementation
// Header formatting for MU4 of (l,m)=(3,3) multipole
double XLALSimIMRPhenomXCP_MU4_l3m3( double theta, double eta, double a1 ){

	/*
	Hola, soy un codigo escribido por "4b_document_fits.py". Dat script is geschreven door een mens die niet kan worden genoemd.
	*/

	// Preliminaries
	double u = cos(theta);
	double u2 = u*u;
	double u3 = u2*u;
	UNUSED double u4 = u3*u;
	double a12 = a1*a1;
	double a13 = a12*a1;
	double eta2 = eta*eta;
	double mu4;

	// Evaluate fit for this parameter
	mu4 = 5.45551261e+01*a1 + -5.83827248e+00 + 5.81966619e+01*eta + -9.48923334e-01*u + 3.72899345e+00*(u*a1) + -1.38994806e+02*a12 + 1.00036537e+01*(u*u) + 3.04562293e+01*(u*eta) + -5.90879466e+02*(eta*a1) + -1.37615553e+02*eta2 + -1.99481881e+02*(u*eta*a1) + -6.91242427e+01*(u2*eta) + 1.11149781e+01*(u2*u) + -2.23772698e+02*(u2*a1) + 9.66405418e+01*(a12*a1) + 1.56222297e+03*(eta*a12) + 1.33032215e+03*(eta2*a1) + -3.62697305e+03*(eta2*a12) + -7.71277331e+01*(u3*a1) + 6.85036555e+00*(u3*u) + -1.10803299e+03*(eta*a13) + 3.94842342e+00*(u*a13) + -1.93469074e+02*(u3*eta) + 2.34388005e+03*(u2*eta*a1) + 6.56809120e+02*(u2*a12) + 3.84613679e+02*(u*eta*a12) + -7.22522187e+03*(u2*eta*a12) + -4.88515693e+03*(u2*eta2*a1) + 1.58527884e+02*(u3*a12) + 1.44620260e+03*(u3*eta*a1) + -1.78863857e+02*(u4*eta) + 2.59728329e+03*(eta2*a13) + 4.82980813e+02*(u3*eta2) + -4.93259096e+02*(u2*a13) + -3.36531843e+02*(u*eta*a13) + 7.26098680e+01*(u4*a1) + -3.72112274e+03*(u3*eta2*a1) + -3.02945607e+02*(u4*a12) + 5.52454100e+03*(u2*eta*a13) + 6.39855808e+02*(u4*eta2) + 1.62955355e+04*(u2*eta2*a12) + -3.18908657e+03*(u3*eta*a12) + 3.51554728e+02*(u*eta2*a13) + -1.13423452e+02*(u3*a13) + 2.49558837e+02*(u4*a13) + -2.04608317e+03*(u4*eta2*a1) + 2.32483877e+03*(u3*eta*a13) + 1.57948011e+03*(u4*eta*a12) + -1.26913975e+04*(u2*eta2*a13) + 8.52885695e+03*(u3*eta2*a12) + 1.23469390e+03*(u4*eta2*a12) + -6.52298706e+03*(u3*eta2*a13) + -1.55267929e+03*(u4*eta*a13);

	// Return answer
	return mu4;

} // END of MU4 (3,3) fit implementation



// NU4 fit implementation
// Header formatting for NU4 of (l,m)=(3,3) multipole
double XLALSimIMRPhenomXCP_NU4_l3m3( double theta, double eta, double a1 ){

	/*
	Hola, soy un codigo escribido por "4b_document_fits.py". Dat script is geschreven door een mens die niet kan worden genoemd.
	*/

	// Preliminaries
	double u = cos(theta);
	double u2 = u*u;
	double u3 = u2*u;
	UNUSED double u4 = u3*u;
	double a12 = a1*a1;
	double a13 = a12*a1;
	double eta2 = eta*eta;
	double nu4;

	// Evaluate fit for this parameter
	nu4 = 3.40013838e+02*a1 + -4.17016718e+01 + 5.08579504e+02*eta + -1.24767473e+01*u + 1.15267257e+02*(u*a1) + -8.66645766e+02*a12 + 3.28617079e+01*(u*u) + 2.43586995e+01*(u*eta) + -4.26289804e+03*(eta*a1) + -1.49491758e+03*eta2 + -2.02461336e+02*(u*a12) + -5.76425487e+02*(u*eta*a1) + 8.84169301e+01*(u2*u) + -5.65307237e+02*(u2*a1) + 1.98452373e+02*(u*eta2) + 6.44330519e+02*(a12*a1) + 1.07376871e+04*(eta*a12) + 1.23778355e+04*(eta2*a1) + -3.07813731e+04*(eta2*a12) + -6.03242124e+02*(u3*a1) + 3.04356715e+01*(u3*u) + -7.91116730e+03*(eta*a13) + 8.87519901e+01*(u*a13) + -9.64974292e+02*(u3*eta) + -1.26378281e+03*(u2*eta2) + 4.30324935e+03*(u2*eta*a1) + 1.89884930e+03*(u2*a12) + 9.63096376e+02*(u*eta*a12) + -1.77061941e+04*(u2*eta*a12) + -4.26844170e+03*(u2*eta2*a1) + 1.13882913e+03*(u3*a12) + 6.56751365e+03*(u3*eta*a1) + -8.70688504e+02*(u4*eta) + 2.25489703e+04*(eta2*a13) + 2.63477737e+03*(u3*eta2) + -1.61074203e+03*(u2*a13) + -2.59765577e+02*(u*eta*a13) + -1.79777245e+04*(u3*eta2*a1) + -5.38787951e+02*(u4*a12) + 1.60688664e+04*(u2*eta*a13) + 3.48393256e+03*(u4*eta2) + 3.38429907e+04*(u2*eta2*a12) + 3.19433450e+03*(u4*eta*a1) + -1.23148926e+04*(u3*eta*a12) + -4.31784001e+02*(u*eta2*a13) + -6.34889573e+02*(u3*a13) + 6.17705749e+02*(u4*a13) + -1.58904907e+04*(u4*eta2*a1) + 6.84581302e+03*(u3*eta*a13) + -3.51385116e+04*(u2*eta2*a13) + 3.39978704e+04*(u3*eta2*a12) + 1.45951339e+04*(u4*eta2*a12) + -1.94203565e+04*(u3*eta2*a13) + -3.31109664e+03*(u4*eta*a13);

	// Return answer
	return nu4;

} // END of NU4 (3,3) fit implementation



// NU5 fit implementation
// Header formatting for NU5 of (l,m)=(3,3) multipole
double XLALSimIMRPhenomXCP_NU5_l3m3( double theta, double eta, double a1 ){

	/*
	Hola, soy un codigo escribido por "4b_document_fits.py". Dat script is geschreven door een mens die niet kan worden genoemd.
	*/

	// Preliminaries
	double u = cos(theta);
	double u2 = u*u;
	double u3 = u2*u;
	UNUSED double u4 = u3*u;
	double a12 = a1*a1;
	double a13 = a12*a1;
	double eta2 = eta*eta;
	double nu5;

	// Evaluate fit for this parameter
	nu5 = -2.34201019e-01*a1 + 1.57947479e-02 + -4.63424723e-01*eta + -3.62638047e-02*u + 1.46394047e-01*(u*a1) + 3.98741298e-01*a12 + 1.23478549e-01*(u*u) + 8.42584471e-01*(u*eta) + 4.40840352e+00*(eta*a1) + 1.49164663e+00*eta2 + -1.32092407e-01*(u*a12) + -4.59948687e+00*(u*eta*a1) + -7.24090808e-01*(u2*eta) + 6.61979843e-02*(u2*u) + -6.16027245e-01*(u2*a1) + -2.89737525e+00*(u*eta2) + -2.93739518e-01*(a12*a1) + -9.13440506e+00*(eta*a12) + -1.50854982e+01*(eta2*a1) + 3.25867708e+01*(eta2*a12) + -2.54339746e-01*(u3*a1) + -1.73628971e-01*(u3*u) + 6.11601500e+00*(eta*a13) + -1.36160854e+00*(u3*eta) + 2.31093344e+00*(u2*eta2) + 7.44899921e-01*(u2*eta*a1) + 1.41007673e+00*(u2*a12) + 1.66754135e+01*(u*eta2*a1) + 7.92815922e+00*(u*eta*a12) + -1.07275000e+00*(u2*eta*a12) + -3.13084759e+01*(u*eta2*a12) + 2.05531046e-01*(u3*a12) + 7.04140887e+00*(u3*eta*a1) + 1.36171424e+00*(u4*eta) + -2.18718577e+01*(eta2*a13) + 4.41435255e+00*(u3*eta2) + -8.92233167e-01*(u2*a13) + -4.86562861e+00*(u*eta*a13) + 8.35427916e-01*(u4*a1) + -2.29938199e+01*(u3*eta2*a1) + -1.53357656e+00*(u4*a12) + -5.16778504e+00*(u4*eta2) + -6.33643518e+00*(u2*eta2*a12) + -3.53706296e+00*(u4*eta*a1) + -1.15068359e+01*(u3*eta*a12) + 2.07410794e+01*(u*eta2*a13) + 8.10019340e-01*(u4*a13) + 1.42674917e+01*(u4*eta2*a1) + 6.71217152e+00*(u3*eta*a13) + 3.86284661e+00*(u4*eta*a12) + 8.57638494e+00*(u2*eta2*a13) + 4.00126182e+01*(u3*eta2*a12) + -1.49674798e+01*(u4*eta2*a12) + -2.50754463e+01*(u3*eta2*a13);

	// Return answer
	return nu5;

} // END of NU5 (3,3) fit implementation



// NU6 fit implementation
// Header formatting for NU6 of (l,m)=(3,3) multipole
double XLALSimIMRPhenomXCP_NU6_l3m3( double theta, double eta, double a1 ){

	/*
	Hola, soy un codigo escribido por "4b_document_fits.py". Dat script is geschreven door een mens die niet kan worden genoemd.
	*/

	// Preliminaries
	double u = cos(theta);
	double u2 = u*u;
	double u3 = u2*u;
	UNUSED double u4 = u3*u;
	double a12 = a1*a1;
	double a13 = a12*a1;
	double eta2 = eta*eta;
	double nu6;

	// Evaluate fit for this parameter
	nu6 = -6.13299936e-01*a1 + 7.74945166e-02 + -9.92293859e-01*eta + -5.78117599e-02*u + 4.92166340e-01*(u*a1) + 1.46165576e+00*a12 + -6.73176061e-02*(u*u) + 9.09129706e-01*(u*eta) + 7.95284252e+00*(eta*a1) + 2.91356725e+00*eta2 + -1.18937142e+00*(u*a12) + -7.25685345e+00*(u*eta*a1) + 5.35830697e-01*(u2*eta) + 1.00792856e+00*(u2*a1) + -3.26954277e+00*(u*eta2) + -1.02001672e+00*(a12*a1) + -1.84831156e+01*(eta*a12) + -2.25535556e+01*(eta2*a1) + 5.15546641e+01*(eta2*a12) + -1.98044512e-01*(u3*a1) + -3.47283165e-02*(u3*u) + 1.26559281e+01*(eta*a13) + 7.77236846e-01*(u*a13) + -1.19418679e+01*(u2*eta*a1) + -2.53178466e+00*(u2*a12) + 2.49494198e+01*(u*eta2*a1) + 1.66810275e+01*(u*eta*a12) + 3.16948169e+01*(u2*eta*a12) + -5.53395022e+01*(u*eta2*a12) + 2.75031189e+01*(u2*eta2*a1) + 8.02012853e-01*(u3*a12) + 2.79773137e+00*(u3*eta*a1) + 6.55005753e-01*(u4*eta) + -3.50059319e+01*(eta2*a13) + 1.80306799e+00*(u2*a13) + -1.08164305e+01*(u*eta*a13) + -9.05094719e+00*(u3*eta2*a1) + -2.29554462e+01*(u2*eta*a13) + -2.35859284e+00*(u4*eta2) + -8.15190127e+01*(u2*eta2*a12) + -9.42415197e-01*(u4*eta*a1) + -1.05656892e+01*(u3*eta*a12) + 3.50558271e+01*(u*eta2*a13) + -6.77435810e-01*(u3*a13) + 5.53577694e+00*(u4*eta2*a1) + 8.75921812e+00*(u3*eta*a13) + 8.03289179e-01*(u4*eta*a12) + 6.20721994e+01*(u2*eta2*a13) + 3.24316273e+01*(u3*eta2*a12) + -5.66444592e+00*(u4*eta2*a12) + -2.59354599e+01*(u3*eta2*a13);

	// Return answer
	return nu6;

} // END of NU6 (3,3) fit implementation



// ZETA1 fit implementation
// Header formatting for ZETA1 of (l,m)=(3,3) multipole
double XLALSimIMRPhenomXCP_ZETA1_l3m3( double theta, double eta, double a1 ){

	/*
	Hola, soy un codigo escribido por "4b_document_fits.py". Dat script is geschreven door een mens die niet kan worden genoemd.
	*/

	// Preliminaries
	double u = cos(theta);
	double u2 = u*u;
	double u3 = u2*u;
	UNUSED double u4 = u3*u;
	double a12 = a1*a1;
	double a13 = a12*a1;
	double eta2 = eta*eta;
	double zeta1;

	// Evaluate fit for this parameter
	zeta1 = -2.04229254e-02*a1 + 2.32463376e-03 + -3.26611196e-02*eta + 5.54613059e-02*a12 + 5.97615154e-03*(u*u) + 3.11165343e-01*(eta*a1) + 1.14861293e-01*eta2 + -8.42724865e-02*(u2*eta) + 1.53672126e-03*(u*a12) + 2.13562451e-02*(u*eta*a1) + -7.82250711e-03*(u2*u) + -4.26208422e-02*(a12*a1) + -8.24789278e-01*(eta*a12) + -1.09735572e+00*(eta2*a1) + 2.84361222e+00*(eta2*a12) + 5.67868535e-02*(u3*a1) + -1.37057681e-02*(u3*u) + 6.28736844e-01*(eta*a13) + 1.02072525e-01*(u3*eta) + 2.74232998e-01*(u2*eta2) + -7.72407346e-02*(u2*a12) + -1.21999759e-01*(u*eta2*a1) + -9.38376631e-02*(u*eta*a12) + 1.03646838e+00*(u2*eta*a12) + 5.36631834e-01*(u*eta2*a12) + -1.20974593e-01*(u3*a12) + -7.59902876e-01*(u3*eta*a1) + 1.90410879e-01*(u4*eta) + -2.14502225e+00*(eta2*a13) + -3.01242096e-01*(u3*eta2) + 9.50104118e-02*(u2*a13) + 7.72491674e-02*(u*eta*a13) + 5.89485924e-02*(u4*a1) + 2.32663587e+00*(u3*eta2*a1) + -5.46860247e-02*(u4*a12) + -1.26587982e+00*(u2*eta*a13) + -5.68884508e-01*(u4*eta2) + -3.12594554e+00*(u2*eta2*a12) + -7.87092191e-01*(u4*eta*a1) + 1.66421687e+00*(u3*eta*a12) + -5.08025947e-01*(u*eta2*a13) + 8.30425618e-02*(u3*a13) + 2.24312749e+00*(u4*eta2*a1) + -1.17127699e+00*(u3*eta*a13) + 7.22404563e-01*(u4*eta*a12) + 3.78338598e+00*(u2*eta2*a13) + -5.29115134e+00*(u3*eta2*a12) + -2.04192729e+00*(u4*eta2*a12) + 3.83627132e+00*(u3*eta2*a13);

	// Return answer
	return zeta1;

} // END of ZETA1 (3,3) fit implementation



// ZETA2 fit implementation
// Header formatting for ZETA2 of (l,m)=(3,3) multipole
double XLALSimIMRPhenomXCP_ZETA2_l3m3( double theta, double eta, double a1 ){

	/*
	Hola, soy un codigo escribido por "4b_document_fits.py". Dat script is geschreven door een mens die niet kan worden genoemd.
	*/

	// Preliminaries
	double u = cos(theta);
	double u2 = u*u;
	double u3 = u2*u;
	UNUSED double u4 = u3*u;
	double a12 = a1*a1;
	double a13 = a12*a1;
	double eta2 = eta*eta;
	double zeta2;

	// Evaluate fit for this parameter
	zeta2 = 1.09732309e+03*a1 + -1.29956302e+02 + 1.71171206e+03*eta + -2.40434643e+01*u + 1.90189840e+02*(u*a1) + -2.83466689e+03*a12 + 6.89141243e+01*(u*eta) + -1.51839054e+04*(eta*a1) + -5.20266434e+03*eta2 + -2.55247517e+02*(u*a12) + -1.51769049e+03*(u*eta*a1) + 9.62136977e+02*(u2*eta) + 3.40492973e+02*(u2*u) + -1.44812934e+03*(u2*a1) + 2.11476872e+03*(a12*a1) + 3.87129738e+04*(eta*a12) + 4.78393934e+04*(eta2*a1) + -1.21366735e+05*(eta2*a12) + -2.43542648e+03*(u3*a1) + 3.24006670e+02*(u3*u) + -2.87255552e+04*(eta*a13) + -4.20052272e+03*(u3*eta) + -5.85607014e+03*(u2*eta2) + 1.30667775e+04*(u2*eta*a1) + 5.67744960e+03*(u2*a12) + 4.96146668e+03*(u*eta2*a1) + 2.88780982e+03*(u*eta*a12) + -6.21991851e+04*(u2*eta*a12) + -1.41087557e+04*(u*eta2*a12) + -2.33926292e+04*(u2*eta2*a1) + 4.67155539e+03*(u3*a12) + 3.11353523e+04*(u3*eta*a1) + -5.54586429e+03*(u4*eta) + 8.99725512e+04*(eta2*a13) + 1.27104408e+04*(u3*eta2) + -5.08830988e+03*(u2*a13) + -6.20582797e+02*(u*eta*a13) + -1.08127769e+03*(u4*a1) + -9.75198579e+04*(u3*eta2*a1) + 5.92129181e+04*(u2*eta*a13) + 1.90829556e+04*(u4*eta2) + 1.55585964e+05*(u2*eta2*a12) + 2.14966350e+04*(u4*eta*a1) + -6.18699333e+04*(u3*eta*a12) + 8.74120677e+03*(u*eta2*a13) + -2.71690257e+03*(u3*a13) + 1.11932540e+03*(u4*a13) + -7.63133173e+04*(u4*eta2*a1) + 3.75485915e+04*(u3*eta*a13) + -1.47864788e+04*(u4*eta*a12) + -1.58997119e+05*(u2*eta2*a13) + 2.01573722e+05*(u3*eta2*a12) + 6.94684577e+04*(u4*eta2*a12) + -1.28211338e+05*(u3*eta2*a13) + -5.41058820e+03*(u4*eta*a13);

	// Return answer
	return zeta2;

} // END of ZETA2 (3,3) fit implementation



// NU0 fit implementation
// Header formatting for NU0 of (l,m)=(3,3) multipole
double XLALSimIMRPhenomXCP_NU0_l3m3( double theta, double eta, double a1 ){

	/*
	Hola, soy un codigo escribido por "4b_document_fits.py". Dat script is geschreven door een mens die niet kan worden genoemd.
	*/

	// Preliminaries
	double u = cos(theta);
	double u2 = u*u;
	double u3 = u2*u;
	UNUSED double u4 = u3*u;
	double a12 = a1*a1;
	double a13 = a12*a1;
	double eta2 = eta*eta;
	double nu0;

	// Evaluate fit for this parameter
	nu0 = -1.37315476e+04*a1 + 1.85417933e+03 + -2.29264256e+04*eta + 8.49715181e+02*u + -2.54555598e+03*(u*a1) + 3.50368649e+04*a12 + 2.52832580e+03*(u*u) + -4.47310560e+03*(u*eta) + 1.83765758e+05*(eta*a1) + 6.49025176e+04*eta2 + -6.61950758e+03*(u*a12) + -5.29860928e+04*(u2*eta) + 3.13972938e+03*(u2*u) + -1.01208675e+04*(u2*a1) + 9.11725618e+03*(u*eta2) + -2.63888862e+04*(a12*a1) + -4.56265594e+05*(eta*a12) + -5.38777740e+05*(eta2*a1) + 1.32897091e+06*(eta2*a12) + -5.14813058e+04*(u3*a1) + -6.25916711e+03*(u3*u) + 3.38084277e+05*(eta*a13) + 1.07256261e+04*(u*a13) + -5.34056649e+04*(u3*eta) + 2.63888356e+05*(u2*eta2) + 2.66420643e+05*(u2*eta*a1) + -1.12065646e+03*(u2*a12) + 1.24613876e+05*(u*eta*a12) + -2.83733264e+05*(u2*eta*a12) + -2.91072717e+05*(u*eta2*a12) + -1.44041457e+06*(u2*eta2*a1) + 1.70362237e+05*(u3*a12) + 7.31758593e+05*(u3*eta*a1) + 8.94546863e+04*(u4*eta) + -9.82188226e+05*(eta2*a13) + 1.68885321e+05*(u3*eta2) + 1.55237057e+04*(u2*a13) + -1.47996898e+05*(u*eta*a13) + 3.07686912e+04*(u4*a1) + -2.18662877e+06*(u3*eta2*a1) + -2.33463167e+04*(u4*a12) + -3.40425969e+05*(u4*eta2) + 2.15682318e+06*(u2*eta2*a12) + -4.12811039e+05*(u4*eta*a1) + -2.29829471e+06*(u3*eta*a12) + 3.38544074e+05*(u*eta2*a13) + -1.37106456e+05*(u3*a13) + -7.62864480e+03*(u4*a13) + 1.60963550e+06*(u4*eta2*a1) + 1.81604885e+06*(u3*eta*a13) + 2.83804938e+05*(u4*eta*a12) + -8.30158123e+05*(u2*eta2*a13) + 6.73203864e+06*(u3*eta2*a12) + -1.52802555e+06*(u4*eta2*a12) + -5.26694866e+06*(u3*eta2*a13) + 1.29757183e+05*(u4*eta*a13);

	// Return answer
	return nu0;

} // END of NU0 (3,3) fit implementation

#ifdef __cplusplus
}
#endif
