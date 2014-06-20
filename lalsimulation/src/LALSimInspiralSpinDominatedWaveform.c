/*
 * Copyright (C) 2014 M. Tapai
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
 *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 *  MA  02111-1307  USA
 */

#include <math.h>
#include <complex.h>

#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/Date.h>
#include <lal/FrequencySeries.h>
#include <lal/TimeSeries.h>
#include <lal/TimeFreqFFT.h>
#include <lal/Units.h>
#include <lal/LALSimInspiral.h>
#include <lal/LALAdaptiveRungeKutta4.h>

#define LAL_SDW_ABSOLUTE_TOLERANCE 1.e-12
#define LAL_SDW_RELATIVE_TOLERANCE 1.e-12
#define LAL_SDW_NUM_VARIABLES 3
#define LAL_SDW_MAX_PN_PARAM 0.8


#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif


/**
 * Structure containing the prefered variabloes for Spin-Dominated waveforms.
 */
typedef struct tagLALSDWaveformParams {
	REAL8 totalmass; //total mass of the binary
	REAL8 nu;  // mass ratio, <1
	REAL8 chi1; // chi1 dimensionless spin parameter
	REAL8 dist; // distance to the source
	REAL8 kappa1; // angle between L and S1
	REAL8 beta1;  // angle between J and S1
	REAL8 theta; // angle between J and N
	REAL8 eps; // PN paramter
	REAL8 xi; // second small parameter
	REAL8 omega;
	int pnamp;
	int pnphase;
	REAL8 ccoeff0pn[4];
	REAL8 ccoeff0_5pn[22];
	REAL8 ccoeff1pn[30];
	REAL8 ccoeff1_5pn[34];
	REAL8 prevdomega;
} LALSDWaveformParams;

static INT4 XLALSpinDominatedWaveformStoppingTest(UNUSED REAL8 t, const REAL8 values[], REAL8 dvalues[], UNUSED void *mparams);
static INT4 XLALSpinDominatedWaveformDerivatives(UNUSED REAL8 t, const REAL8 values[], REAL8 dvalues[], void *mparams);
int XLALSpinDominatedWaveformConstantCoefficients (LALSDWaveformParams * params);

/**
 * Function for allocating memory for a matrix
 */
REAL8 * XLALDmatrix(INT8 nrh, INT8 nch);

REAL8 * XLALDmatrix(INT8 nrh, INT8 nch){
INT8 size = (nrh)*(nch)*sizeof(REAL8);
REAL8 *ptr = (REAL8 *)LALMalloc(size);
if (ptr != NULL) return ptr;
printf("malloc error");
return NULL;
}

/**
 * Function for freeing memory for a matrix
 */
void XLALFreeDmatrix(REAL8 *m);

void XLALFreeDmatrix(REAL8 *m){
LALFree(m);
}

/**
 * Function for calculating the constant coefficients of Spin-Dominated waveforms
 * See tables 1 to 5 in the appendix of Arxiv:1209.1722
 */
int XLALSpinDominatedWaveformConstantCoefficients (LALSDWaveformParams * params){

int i,j;
REAL8 *acoeff0pn, *b0coeff0pn, *d0coeff0pn, *acoeff0_5pn, *b0coeff0_5pn, *d0coeff0_5pn, *acoeff1pn, *b0coeff1pn, *d0coeff1pn, *b1coeff1pn, *d1coeff1pn,	*acoeff1_5pn, *b0coeff1_5pn, *d0coeff1_5pn, *b1coeff1_5pn, *d1coeff1_5pn;

REAL8 sk1=sin(params->kappa1);
REAL8 sk12 = sk1 * sk1;
REAL8 st = sin (params->theta) ;
REAL8 st2 = st*st;
REAL8 st4 = st2*st2;
REAL8 ck1=cos(params->kappa1);
REAL8 ct = cos (params->theta) ;
REAL8 ct2 = ct*ct;
REAL8 c2t = cos (2. * params->theta);
REAL8 kplus = ck1-1.;
REAL8 kmin = ck1+1.;

acoeff0pn=XLALDmatrix(2,2);
b0coeff0pn=XLALDmatrix(2,2);
d0coeff0pn=XLALDmatrix(2,2);
acoeff0_5pn=XLALDmatrix(11,2);
b0coeff0_5pn=XLALDmatrix(11,2);
d0coeff0_5pn=XLALDmatrix(11,2);
acoeff1pn=XLALDmatrix(15,2);
b0coeff1pn=XLALDmatrix(15,2);
d0coeff1pn=XLALDmatrix(15,2);
b1coeff1pn=XLALDmatrix(15,2);
d1coeff1pn=XLALDmatrix(15,2);
acoeff1_5pn=XLALDmatrix(17,2);
b0coeff1_5pn=XLALDmatrix(17,2);
d0coeff1_5pn=XLALDmatrix(17,2);
b1coeff1_5pn=XLALDmatrix(17,2);
d1coeff1_5pn=XLALDmatrix(17,2);
/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%   0 PN   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/
acoeff0pn[0+2*0]=-2.*kplus;
acoeff0pn[0+2*1]=2.*kmin;
acoeff0pn[1+2*0]=-kplus;
acoeff0pn[1+2*1]=kmin;
b0coeff0pn[0+2*0]=-1.;
b0coeff0pn[0+2*1]=-1.;
b0coeff0pn[1+2*0]=-2.;
b0coeff0pn[1+2*1]=-2.;
d0coeff0pn[0+2*0]=0.;
d0coeff0pn[0+2*1]=0.;
d0coeff0pn[1+2*0]=0.;
d0coeff0pn[1+2*1]=0.;
for (i=0; i<=1; i++){
	for (j=0; j<=1; j++){
		params->ccoeff0pn[i+2*j] = acoeff0pn[i+2*j] + sk12*(b0coeff0pn[i+2*j]+ck1*d0coeff0pn[i+2*j]);
	}
}
/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%   0.5 PN   %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/
acoeff0_5pn[0+11*0]=4.*kplus*(6.-st2);
acoeff0_5pn[0+11*1]=4.*kmin*(6.-st2);
acoeff0_5pn[1+11*0]=4.*kplus;
acoeff0_5pn[1+11*1]=4.*kmin;
acoeff0_5pn[2+11*0]=2.*kplus*(6.-st2);
acoeff0_5pn[2+11*1]=-2.*kmin*(6.-st2);
acoeff0_5pn[3+11*0]=12.*kplus;
acoeff0_5pn[3+11*1]=12.*kmin;
acoeff0_5pn[4+11*0]=2.*kplus*(2.*st2-3.);
acoeff0_5pn[4+11*1]=-2.*kmin*(2.*st2-3.);
acoeff0_5pn[5+11*0]=-2.*kplus;
acoeff0_5pn[5+11*1]=-2.*kmin;
acoeff0_5pn[6+11*0]=-2.*params->ccoeff0pn[1+2*0]*(6.-st2);
acoeff0_5pn[6+11*1]=2.*params->ccoeff0pn[1+2*1]*(6.-st2);
acoeff0_5pn[7+11*0]=2.*kplus;
acoeff0_5pn[7+11*1]=-2.*kmin;
acoeff0_5pn[8+11*0]=(44.-34.*st2+2.*(5.*st2-46.)*ck1);
acoeff0_5pn[8+11*1]=(44.-34.*st2-2.*(5.*st2-46.)*ck1);
acoeff0_5pn[9+11*0]=22.+46.*ck1;
acoeff0_5pn[9+11*1]=22.-46.*ck1;
acoeff0_5pn[10+11*0]=-2.*kplus*(3.-2.*st2);
acoeff0_5pn[10+11*1]=-2.*kmin*(3.-2.*st2);
b0coeff0_5pn[0+11*0]=(46.-5.*st2);
b0coeff0_5pn[0+11*1]=-(46.-5.*st2);
b0coeff0_5pn[1+11*0]=3.;
b0coeff0_5pn[1+11*1]=-3.;
b0coeff0_5pn[2+11*0]=(2.-3.*st2);
b0coeff0_5pn[2+11*1]=(2.-3.*st2);
b0coeff0_5pn[3+11*0]=23.;
b0coeff0_5pn[3+11*1]=-23.;
b0coeff0_5pn[4+11*0]=-c2t;
b0coeff0_5pn[4+11*1]=-c2t;
b0coeff0_5pn[5+11*0]=-4.;
b0coeff0_5pn[5+11*1]=4.;
b0coeff0_5pn[6+11*0]=0.;
b0coeff0_5pn[6+11*1]=0.;
b0coeff0_5pn[7+11*0]=3.;
b0coeff0_5pn[7+11*1]=3.;
b0coeff0_5pn[8+11*0]=-15.*(2.-3.*st2);
b0coeff0_5pn[8+11*1]=-15.*(2.-3.*st2);
b0coeff0_5pn[9+11*0]=-15.;
b0coeff0_5pn[9+11*1]=-15.;
b0coeff0_5pn[10+11*0]=-4.*(3.-2.*st2);
b0coeff0_5pn[10+11*1]=4.*(3.-2.*st2);
d0coeff0_5pn[0+11*0]=5.*(3.*st2-2.);
d0coeff0_5pn[0+11*1]=5.*(3.*st2-2.);
d0coeff0_5pn[1+11*0]=-1.;
d0coeff0_5pn[1+11*1]=-1.;
d0coeff0_5pn[2+11*0]=0.;
d0coeff0_5pn[2+11*1]=0.;
d0coeff0_5pn[3+11*0]=-5.;
d0coeff0_5pn[3+11*1]=-5.;
d0coeff0_5pn[4+11*0]=0.;
d0coeff0_5pn[4+11*1]=0.;
d0coeff0_5pn[5+11*0]=3.;
d0coeff0_5pn[5+11*1]=3.;
d0coeff0_5pn[6+11*0]=-3.*(2.-3.*st2);
d0coeff0_5pn[6+11*1]=-3.*(2.-3.*st2);
d0coeff0_5pn[7+11*0]=0.;
d0coeff0_5pn[7+11*1]=0.;
d0coeff0_5pn[8+11*0]=0.;
d0coeff0_5pn[8+11*1]=0.;
d0coeff0_5pn[9+11*0]=0.;
d0coeff0_5pn[9+11*1]=0.;
d0coeff0_5pn[10+11*0]=3.*c2t;
d0coeff0_5pn[10+11*1]=3.*c2t;
for (i=0; i<=10; i++){
	for (j=0; j<=1; j++){
		params->ccoeff0_5pn[i+11*j] = acoeff0_5pn[i+11*j] + sk12*(b0coeff0_5pn[i+11*j]+ck1*d0coeff0_5pn[i+11*j]);
	}
}
/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%   1 PN   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/
acoeff1pn[0+15*0]=8.*kplus;
acoeff1pn[0+15*1]=-8.*kmin;
acoeff1pn[1+15*0]=6.*kplus*(st2+5.);
acoeff1pn[1+15*1]=6.*kmin*(st2+5.);
acoeff1pn[2+15*0]=2.*kplus*(4.-st2);
acoeff1pn[2+15*1]=2.*kmin*(4.-st2);
acoeff1pn[3+15*0]=2.*kplus*(2.*st4+11.*st2-38.);
acoeff1pn[3+15*1]=-2.*kmin*(2.*st4+11.*st2-38.);
acoeff1pn[4+15*0]=6.*kplus*(3.*st2+5.);
acoeff1pn[4+15*1]=6.*kmin*(3.*st2+5.);
acoeff1pn[5+15*0]=2.*kplus*(4.*st2+19.);
acoeff1pn[5+15*1]=-2.*kmin*(4.*st2+19.);
acoeff1pn[6+15*0]=-2.*kplus*(3.*st2-4.);
acoeff1pn[6+15*1]=-2.*kmin*(3.*st2-4.);
acoeff1pn[7+15*0]=2.*kplus*(4.-st2);
acoeff1pn[7+15*1]=-2.*kmin*(4.-st2);
acoeff1pn[8+15*0]=6.*kplus*(5.+st2);
acoeff1pn[8+15*1]=-6.*kmin*(5.+st2);
acoeff1pn[9+15*0]=-4.*kplus;
acoeff1pn[9+15*1]=4.*kmin;
acoeff1pn[10+15*0]=kplus*(22.+29.*st2-16.*st4);
acoeff1pn[10+15*1]=kmin*(22.+29.*st2-16.*st4);
acoeff1pn[11+15*0]=2.*kplus;
acoeff1pn[11+15*1]=2.*kmin;
acoeff1pn[12+15*0]=6.*kplus*(3.*st2+5.);
acoeff1pn[12+15*1]=-6.*kmin*(3.*st2+5.);
acoeff1pn[13+15*0]=-kplus*(20.*st2+11.);
acoeff1pn[13+15*1]=-kmin*(20.*st2+11.);
acoeff1pn[14+15*0]=-2.*kplus*(3.*st2-4.);
acoeff1pn[14+15*1]=2.*kmin*(3.*st2-4.);
b0coeff1pn[0+15*0]=8.;
b0coeff1pn[0+15*1]=8.;
b0coeff1pn[1+15*0]=-18.+7.*st2;
b0coeff1pn[1+15*1]=18.-7.*st2;
b0coeff1pn[2+15*0]=-3.*st2+6.;
b0coeff1pn[2+15*1]=3.*st2-6.;
b0coeff1pn[3+15*0]=(-22.-29.*st2+16.*st4);
b0coeff1pn[3+15*1]=(-22.-29.*st2+16.*st4);
b0coeff1pn[4+15*0]=26.*st2-18.;
b0coeff1pn[4+15*1]=-26.*st2+18.;
b0coeff1pn[5+15*0]=11.+20.*st2;
b0coeff1pn[5+15*1]=11.+20.*st2;
b0coeff1pn[6+15*0]=-6.*st2+6.;
b0coeff1pn[6+15*1]=6.*st2-6.;
b0coeff1pn[7+15*0]=2.*(11.-5.*st2);
b0coeff1pn[7+15*1]=2.*(11.-5.*st2);
b0coeff1pn[8+15*0]=6.*(7.+9.*st2);
b0coeff1pn[8+15*1]=6.*(7.+9.*st2);
b0coeff1pn[9+15*0]=-11.;
b0coeff1pn[9+15*1]=-11.;
b0coeff1pn[10+15*0]=-3.*(8.-20.*st2+7.*st4);
b0coeff1pn[10+15*1]=3.*(8.-20.*st2+7.*st4);
b0coeff1pn[11+15*0]=3.;
b0coeff1pn[11+15*1]=-3.;
b0coeff1pn[12+15*0]=3.*(19.*st2+14.);
b0coeff1pn[12+15*1]=3.*(19.*st2+14.);
b0coeff1pn[13+15*0]=12.*c2t;
b0coeff1pn[13+15*1]=-12.*c2t;
b0coeff1pn[14+15*0]=(22.-21.*st2);
b0coeff1pn[14+15*1]=(22.-21.*st2);
d0coeff1pn[0+15*0]=-4.;
d0coeff1pn[0+15*1]=4.;
d0coeff1pn[1+15*0]=(6.-14.*st2);
d0coeff1pn[1+15*1]=(6.-14.*st2);
d0coeff1pn[2+15*0]=2.*(st2-1.);
d0coeff1pn[2+15*1]=2.*(st2-1.);
d0coeff1pn[3+15*0]=-2.*(8.-20.*st2+7.*st4);
d0coeff1pn[3+15*1]=2.*(8.-20.*st2+7.*st4);
d0coeff1pn[4+15*0]=(6.-7.*st2);
d0coeff1pn[4+15*1]=(6.-7.*st2);
d0coeff1pn[5+15*0]=(-16.*st2+8.);
d0coeff1pn[5+15*1]=(16.*st2-8.);
d0coeff1pn[6+15*0]=(3.*st2-2.);
d0coeff1pn[6+15*1]=(3.*st2-2.);
d0coeff1pn[7+15*0]=9.*(st2-2.);
d0coeff1pn[7+15*1]=-9.*(st2-2.);
d0coeff1pn[8+15*0]=3.*(18.-7.*st2);
d0coeff1pn[8+15*1]=-3.*(18.-7.*st2);
d0coeff1pn[9+15*0]=9.;
d0coeff1pn[9+15*1]=-9.;
d0coeff1pn[10+15*0]=4.*(2.-8.*st2+7.*st4);
d0coeff1pn[10+15*1]=4.*(2.-8.*st2+7.*st4);
d0coeff1pn[11+15*0]=-2.;
d0coeff1pn[11+15*1]=-2.;
d0coeff1pn[12+15*0]=6.*(9.-13.*st2);
d0coeff1pn[12+15*1]=-6.*(9.-13.*st2);
d0coeff1pn[13+15*0]=2.*(7.*st2-2.);
d0coeff1pn[13+15*1]=2.*(7.*st2-2.);
d0coeff1pn[14+15*0]=-18.*ct2;
d0coeff1pn[14+15*1]=18.*ct2;
b1coeff1pn[0+15*0]=-1.;
b1coeff1pn[0+15*1]=-1.;
b1coeff1pn[1+15*0]=0.;
b1coeff1pn[1+15*1]=0.;
b1coeff1pn[2+15*0]=0.;
b1coeff1pn[2+15*1]=0.;
b1coeff1pn[3+15*0]=-2.*(2.-8.*st2+7.*st4);
b1coeff1pn[3+15*1]=-2.*(2.-8.*st2+7.*st4);
b1coeff1pn[4+15*0]=0.;
b1coeff1pn[4+15*1]=0.;
b1coeff1pn[5+15*0]=(2.-7.*st2);
b1coeff1pn[5+15*1]=(2.-7.*st2);
b1coeff1pn[6+15*0]=0.;
b1coeff1pn[6+15*1]=0.;
b1coeff1pn[7+15*0]=-8.*ct2;
b1coeff1pn[7+15*1]=-8.*ct2;
b1coeff1pn[8+15*0]=8.*(3.-7.*st2);
b1coeff1pn[8+15*1]=8.*(3.-7.*st2);
b1coeff1pn[9+15*0]=4.;
b1coeff1pn[9+15*1]=4.;
b1coeff1pn[10+15*0]=0.;
b1coeff1pn[10+15*1]=0.;
b1coeff1pn[11+15*0]=0.;
b1coeff1pn[11+15*1]=0.;
b1coeff1pn[12+15*0]=-4.*(7.*st2-6.);
b1coeff1pn[12+15*1]=-4.*(7.*st2-6.);
b1coeff1pn[13+15*0]=0.;
b1coeff1pn[13+15*1]=0.;
b1coeff1pn[14+15*0]=-4.*(2.-3.*st2);
b1coeff1pn[14+15*1]=-4.*(2.-3.*st2);
for (i=0;i<15;i++){
	d1coeff1pn[i+15*0]=0.;
	d1coeff1pn[i+15*1]=0.;
}
for (i=0; i<=14; i++){
	for (j=0; j<=1; j++){
		params->ccoeff1pn[i+15*j] = acoeff1pn[i+15*j] + sk12*(b0coeff1pn[i+15*j]+ck1*d0coeff1pn[i+15*j]) + sk12*sk12*(b1coeff1pn[i+15*j]+ck1*d1coeff1pn[i+15*j]);
	}
}
/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%   1.5 PN   %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/
acoeff1_5pn[0+17*0]=4.*kplus*(st2-6.);
acoeff1_5pn[0+17*1]=-4.*kmin*(st2-6.);
acoeff1_5pn[1+17*0]=4.*kplus*(st4+42.*st2-166.);
acoeff1_5pn[1+17*1]=-4.*kmin*(st4+42.*st2-166.);
acoeff1_5pn[2+17*0]=16.*kplus;
acoeff1_5pn[2+17*1]=16.*kmin;
acoeff1_5pn[3+17*0]=8.*kplus*(st4+8.*st2-28.);
acoeff1_5pn[3+17*1]=8.*kmin*(st4+8.*st2-28.);
acoeff1_5pn[4+17*0]=8.*kplus*(-332.+94.*st2+st4);
acoeff1_5pn[4+17*1]=8.*kmin*(-332.+94.*st2+st4);
acoeff1_5pn[5+17*0]=8.*kplus*(38.-42.*st2-9.*st4);
acoeff1_5pn[5+17*1]=-8.*kmin*(38.-42.*st2-9.*st4);
acoeff1_5pn[6+17*0]=-16.*kplus*(152.-46.*st2-9.*st4);
acoeff1_5pn[6+17*1]=-16.*kmin*(152.-46.*st2-9.*st4);
acoeff1_5pn[7+17*0]=24.*kplus*(3.*st2-10.);
acoeff1_5pn[7+17*1]=-24.*kmin*(3.*st2-10.);
acoeff1_5pn[8+17*0]=-8.*kplus*(160.-204.*st2-63.*st4);
acoeff1_5pn[8+17*1]=-8.*kmin*(160.-204.*st2-63.*st4);
acoeff1_5pn[9+17*0]=4.*kplus*(3.-2.*st2);
acoeff1_5pn[9+17*1]=-4.*kmin*(3.-2.*st2);
acoeff1_5pn[10+17*0]=-8.*kplus*(14.+3.*st2);
acoeff1_5pn[10+17*1]=-8.*kmin*(14.+3.*st2);
acoeff1_5pn[11+17*0]=-16.*kplus*(15.*st2+76.);
acoeff1_5pn[11+17*1]=-16.*kmin*(15.*st2+76.);
acoeff1_5pn[12+17*0]=-8.*kplus*(5.*st2+166.);
acoeff1_5pn[12+17*1]=-8.*kmin*(5.*st2+166.);
acoeff1_5pn[13+17*0]=-8.*kplus*(80.+63.*st2);
acoeff1_5pn[13+17*1]=-8.*kmin*(80.+63.*st2);
acoeff1_5pn[14+17*0]=4.*kplus*(166.-125.*st2-8.*st4);
acoeff1_5pn[14+17*1]=-4.*kmin*(166.-125.*st2-8.*st4);
acoeff1_5pn[15+17*0]=-8.*kplus*(38.-61.*st2-24.*st4);
acoeff1_5pn[15+17*1]=8.*kmin*(38.-61.*st2-24.*st4);
acoeff1_5pn[16+17*0]=8.*kplus*(5.-4.*st2);
acoeff1_5pn[16+17*1]=-8.*kmin*(5.-4.*st2);
b0coeff1_5pn[0+17*0]=(5.*st2-6.);
b0coeff1_5pn[0+17*1]=(5.*st2-6.);
b0coeff1_5pn[1+17*0]=(18.*st4+252.*st2-188.);
b0coeff1_5pn[1+17*1]=(18.*st4+252.*st2-188.);
b0coeff1_5pn[2+17*0]=20.;
b0coeff1_5pn[2+17*1]=-20.;
b0coeff1_5pn[3+17*0]=(9.*st4-90.*st2+56.);
b0coeff1_5pn[3+17*1]=(-9.*st4+90.*st2-56.);
b0coeff1_5pn[4+17*0]=-4.*(1184.-172.*st2-7.*st4);
b0coeff1_5pn[4+17*1]=4.*(1184.-172.*st2-7.*st4);
b0coeff1_5pn[5+17*0]=2.*(46.+48.*st2-99.*st4);
b0coeff1_5pn[5+17*1]=2.*(46.+48.*st2-99.*st4);
b0coeff1_5pn[6+17*0]=-12.*(72.+110.*st2-63.*st4);
b0coeff1_5pn[6+17*1]=12.*(72.+110.*st2-63.*st4);
b0coeff1_5pn[7+17*0]=144.*(st2-2.);
b0coeff1_5pn[7+17*1]=144.*(st2-2.);
b0coeff1_5pn[8+17*0]=-3.*(-204.+406.*st2-189.*st4);
b0coeff1_5pn[8+17*1]=3.*(-204.+406.*st2-189.*st4);
b0coeff1_5pn[9+17*0]=3.-4.*st2;
b0coeff1_5pn[9+17*1]=3.-4.*st2;
b0coeff1_5pn[10+17*0]=28.-31.*st2;
b0coeff1_5pn[10+17*1]=-28.+31.*st2;
b0coeff1_5pn[11+17*0]=-432.-876.*st2;
b0coeff1_5pn[11+17*1]=432.+876.*st2;
b0coeff1_5pn[12+17*0]=-4.*(71.*st2+592.);
b0coeff1_5pn[12+17*1]=4.*(71.*st2+592.);
b0coeff1_5pn[13+17*0]=306.-651.*st2;
b0coeff1_5pn[13+17*1]=-306.+651.*st2;
b0coeff1_5pn[14+17*0]=2.*(94.-173.*st2-24.*st4);
b0coeff1_5pn[14+17*1]=2.*(94.-173.*st2-24.*st4);
b0coeff1_5pn[15+17*0]=-2.*(46.+25.*st2-180.*st4);
b0coeff1_5pn[15+17*1]=-2.*(46.+25.*st2-180.*st4);
b0coeff1_5pn[16+17*0]=48.*ct2;
b0coeff1_5pn[16+17*1]=48.*ct2;
d0coeff1_5pn[0+17*0]=0.;
d0coeff1_5pn[0+17*1]=0.;
d0coeff1_5pn[1+17*0]=(-6.*st4+72.*st2-20.);
d0coeff1_5pn[1+17*1]=(+6.*st4-72.*st2+20.);
d0coeff1_5pn[2+17*0]=-12.;
d0coeff1_5pn[2+17*1]=-12.;
d0coeff1_5pn[3+17*0]=(-15.*st4+22.*st2-8.);
d0coeff1_5pn[3+17*1]=(-15.*st4+22.*st2-8.);
d0coeff1_5pn[4+17*0]=(1920.-2832.*st2-84.*st4);
d0coeff1_5pn[4+17*1]=(1920.-2832.*st2-84.*st4);
d0coeff1_5pn[5+17*0]=6.*(10.-44.*st2+27.*st4);
d0coeff1_5pn[5+17*1]=-6.*(10.-44.*st2+27.*st4);
d0coeff1_5pn[6+17*0]=-4.*(88.-422.*st2+171.*st4);
d0coeff1_5pn[6+17*1]=-4.*(88.-422.*st2+171.*st4);
d0coeff1_5pn[7+17*0]=+12.*(14.-9.*st2);
d0coeff1_5pn[7+17*1]=-12.*(14.-9.*st2);
d0coeff1_5pn[8+17*0]=-9.*(28.-126.*st2+105.*st4);
d0coeff1_5pn[8+17*1]=-9.*(28.-126.*st2+105.*st4);
d0coeff1_5pn[9+17*0]=0.;
d0coeff1_5pn[9+17*1]=0.;
d0coeff1_5pn[10+17*0]=(9.*st2-4.);
d0coeff1_5pn[10+17*1]=(9.*st2-4.);
d0coeff1_5pn[11+17*0]=(-176.+756.*st2);
d0coeff1_5pn[11+17*1]=(-176.+756.*st2);
d0coeff1_5pn[12+17*0]=12.*(7.*st2+80.);
d0coeff1_5pn[12+17*1]=12.*(7.*st2+80.);
d0coeff1_5pn[13+17*0]=(-126.+189.*st2);
d0coeff1_5pn[13+17*1]=(-126.+189.*st2);
d0coeff1_5pn[14+17*0]=2.*(10.-41.*st2+36.*st4);
d0coeff1_5pn[14+17*1]=-2.*(10.-41.*st2+36.*st4);
d0coeff1_5pn[15+17*0]=-6.*(10.-49.*st2+44.*st4);
d0coeff1_5pn[15+17*1]=6.*(10.-49.*st2+44.*st4);
d0coeff1_5pn[16+17*0]=-4.*(7.-8.*st2);
d0coeff1_5pn[16+17*1]=4.*(7.-8.*st2);
b1coeff1_5pn[0+17*0]=0.;
b1coeff1_5pn[0+17*1]=0.;
b1coeff1_5pn[1+17*0]=(-15.*st4+12.*st2-2.);
b1coeff1_5pn[1+17*1]=(-15.*st4+12.*st2-2.);
b1coeff1_5pn[2+17*0]=-5.;
b1coeff1_5pn[2+17*1]=5.;
b1coeff1_5pn[3+17*0]=0.;
b1coeff1_5pn[3+17*1]=0.;
b1coeff1_5pn[4+17*0]=-(236.-294.*st2+21.*st4);
b1coeff1_5pn[4+17*1]=(236.-294.*st2+21.*st4);
b1coeff1_5pn[5+17*0]=3.*(6.-36.*st2+45.*st4);
b1coeff1_5pn[5+17*1]=3.*(6.-36.*st2+45.*st4);
b1coeff1_5pn[6+17*0]=-3.*(232.-510.*st2+243.*st4);
b1coeff1_5pn[6+17*1]=3.*(232.-510.*st2+243.*st4);
b1coeff1_5pn[7+17*0]=9.*(6.-5.*st2);
b1coeff1_5pn[7+17*1]=9.*(6.-5.*st2);
b1coeff1_5pn[8+17*0]=0.;
b1coeff1_5pn[8+17*1]=0.;
b1coeff1_5pn[9+17*0]=0.;
b1coeff1_5pn[9+17*1]=0.;
b1coeff1_5pn[10+17*0]=0.;
b1coeff1_5pn[10+17*1]=0.;
b1coeff1_5pn[11+17*0]=(-348.+591.*st2);
b1coeff1_5pn[11+17*1]=(+348.-591.*st2);
b1coeff1_5pn[12+17*0]=(273.*st2-118.);
b1coeff1_5pn[12+17*1]=(-273.*st2+118.);
b1coeff1_5pn[13+17*0]=0.;
b1coeff1_5pn[13+17*1]=0.;
b1coeff1_5pn[14+17*0]=(2.-13.*st2+12.*st4);
b1coeff1_5pn[14+17*1]=(2.-13.*st2+12.*st4);
b1coeff1_5pn[15+17*0]=-9.*(2.-13.*st2+12.*st4);
b1coeff1_5pn[15+17*1]=-9.*(2.-13.*st2+12.*st4);
b1coeff1_5pn[16+17*0]=-3.*(3.-4.*st2);
b1coeff1_5pn[16+17*1]=-3.*(3.-4.*st2);
d1coeff1_5pn[0+17*0]=0.;
d1coeff1_5pn[0+17*1]=0.;
d1coeff1_5pn[1+17*0]=0.;
d1coeff1_5pn[1+17*1]=0.;
d1coeff1_5pn[2+17*0]=1.;
d1coeff1_5pn[2+17*1]=1.;
d1coeff1_5pn[3+17*0]=0.;
d1coeff1_5pn[3+17*1]=0.;
d1coeff1_5pn[4+17*0]=(28.-126.*st2+105.*st4);
d1coeff1_5pn[4+17*1]=(28.-126.*st2+105.*st4);
d1coeff1_5pn[5+17*0]=0.;
d1coeff1_5pn[5+17*1]=0.;
d1coeff1_5pn[6+17*0]=27.*(8.-22.*st2+15.*st4);
d1coeff1_5pn[6+17*1]=27.*(8.-22.*st2+15.*st4);
d1coeff1_5pn[7+17*0]=0.;
d1coeff1_5pn[7+17*1]=0.;
d1coeff1_5pn[8+17*0]=0.;
d1coeff1_5pn[8+17*1]=0.;
d1coeff1_5pn[9+17*0]=0.;
d1coeff1_5pn[9+17*1]=0.;
d1coeff1_5pn[10+17*0]=0.;
d1coeff1_5pn[10+17*1]=0.;
d1coeff1_5pn[11+17*0]=(-243.*st2+108.);
d1coeff1_5pn[11+17*1]=(-243.*st2+108.);
d1coeff1_5pn[12+17*0]=7.*(2.-3.*st2);
d1coeff1_5pn[12+17*1]=7.*(2.-3.*st2);
d1coeff1_5pn[13+17*0]=0.;
d1coeff1_5pn[13+17*1]=0.;
d1coeff1_5pn[14+17*0]=0.;
d1coeff1_5pn[14+17*1]=0.;
d1coeff1_5pn[15+17*0]=0.;
d1coeff1_5pn[15+17*1]=0.;
d1coeff1_5pn[16+17*0]=0.;
d1coeff1_5pn[16+17*1]=0.;
for (i=0; i<=16; i++){
	for (j=0; j<=1; j++){
		params->ccoeff1_5pn[i+17*j] = acoeff1_5pn[i+17*j] + sk12*(b0coeff1_5pn[i+17*j]+ck1*d0coeff1_5pn[i+17*j]) + sk12*sk12*(b1coeff1_5pn[i+17*j]+ck1*d1coeff1_5pn[i+17*j]);
	}
}
XLALFreeDmatrix(acoeff0pn);
XLALFreeDmatrix(b0coeff0pn);
XLALFreeDmatrix(d0coeff0pn);
XLALFreeDmatrix(acoeff0_5pn);
XLALFreeDmatrix(b0coeff0_5pn);
XLALFreeDmatrix(d0coeff0_5pn);
XLALFreeDmatrix(acoeff1pn);
XLALFreeDmatrix(b0coeff1pn);
XLALFreeDmatrix(d0coeff1pn);
XLALFreeDmatrix(b1coeff1pn);
XLALFreeDmatrix(d1coeff1pn);
XLALFreeDmatrix(acoeff1_5pn);
XLALFreeDmatrix(b0coeff1_5pn);
XLALFreeDmatrix(d0coeff1_5pn);
XLALFreeDmatrix(b1coeff1_5pn);
XLALFreeDmatrix(d1coeff1_5pn);
return XLAL_SUCCESS;
}

/**
* Function building the wavefrom from the calculated parameters at a given time
* For the formulae see the appendix of Arxiv:1209.1722
*/
int XLALSpinDominatedWaveformBuild (LALSDWaveformParams *params, REAL8 expr[], REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, int idx);

int XLALSpinDominatedWaveformBuild (
	LALSDWaveformParams *params,  		/**< The SDW parameters */
	REAL8 expr[], 				/**< The 3 time dependent variables of the waveform at the time indexed by idx */
	REAL8TimeSeries **hplus,        	/**< +-polarization waveform */
	REAL8TimeSeries **hcross,      		/**< x-polarization waveform */
	int idx)
{
REAL8 *waveampcoeffs;
REAL8 ampcoeff;
REAL8 chi1_1=params->chi1;
REAL8 sk1=sin(params->kappa1);
REAL8 sk12 = sk1 * sk1;
REAL8 sk13 = sk12 * sk1;
REAL8 sk14 = sk12 * sk12;
REAL8 s2k1=sin( 2. * params->kappa1);
REAL8 st = sin (params->theta) ;
REAL8 st2 = st*st;
REAL8 st3 = st2*st;
REAL8 s2t = sin (2.* params->theta);
REAL8 st4 = st2*st2;
REAL8 ck1=cos(params->kappa1);
REAL8 c2k1=cos(2. * params->kappa1);
REAL8 ck12 = ck1 * ck1;
REAL8 ct = cos (params->theta) ;
REAL8 ct2 = ct*ct;
REAL8 c2t = cos (2. * params->theta);
REAL8 kplus = ck1-1.;
REAL8 kmin = ck1+1.;
REAL8 omega0=1;
REAL8 v=cbrt(LAL_G_SI*params->totalmass*expr[1]/LAL_C_SI/LAL_C_SI/LAL_C_SI);
params->eps=v*v;
REAL8 eps0pn1pncorr, eps0_5pn1pncorr, eps0pn1_5pncorr;
eps0_5pn1pncorr=params->eps*(1.+v*v*(1.-params->nu/(1.+ params->nu)/(1.+ params->nu)/3.)+v*v*v*(params->nu/(1.+ params->nu)/(1.+ params->nu)+2./3./(1.+ params->nu)/(1.+ params->nu)));
eps0pn1pncorr=params->eps*(1.+v*v*(1.-params->nu/(1.+ params->nu)/(1.+ params->nu)/3.)+v*v*v*(params->nu/(1.+ params->nu)/(1.+ params->nu)+2./3./(1.+ params->nu)/(1.+ params->nu)));
eps0pn1_5pncorr=params->eps*(1.+v*v*(1.-params->nu/(1.+ params->nu)/(1.+ params->nu)/3.)+v*v*v*(params->nu/(1.+ params->nu)/(1.+ params->nu)+2./3./(1.+ params->nu)/(1.+ params->nu))+v*v*v*v*(1.-65./12.*params->nu/(1.+ params->nu)/(1.+ params->nu)-params->chi1*params->chi1*params->nu*params->nu/(1.+ params->nu)/(1.+ params->nu)*(3.*ck12-1.)/4.));
params->xi=params->nu/sqrt(params->eps); // second small parameter
waveampcoeffs=XLALDmatrix(11,2); // amplitude coefficients
REAL8 epssqrt=sqrt(params->eps);
REAL8 phin_1=expr[0];
REAL8 phin_2=2.*expr[0];
REAL8 phin_3=3.*expr[0];
REAL8 phin_4=4.*expr[0];
REAL8 phin_5=5.*expr[0];
REAL8 phi1_1=LAL_PI/2.;
REAL8 psi_1 = 1. * expr[2];
REAL8 psi_2 = 2. * expr[2];
REAL8 psi_3 = 3. * expr[2];
REAL8 psi_4 = 4. * expr[2];
REAL8 psi_5 = 5. * expr[2];
/*
 computing the coefficients
*/
(*hplus)->data->data[idx]=0.;
(*hcross)->data->data[idx]=0.;
ampcoeff=LAL_G_SI*params->totalmass/LAL_C_SI/LAL_C_SI*params->eps*epssqrt*params->xi*2./params->dist;
switch (params->pnamp) {
	case (3):
	waveampcoeffs[8+11*0]=1./12288.*(12.*ct*sk1*st2*(cos(psi_3)*(1701.*(2.-3.*st2)*sk14+72.*sk12*(63.*st2+178.))+cos(psi_1)*(-14.*(2.-3.*st2)*sk14-8.*sk12*(7.*st2+162.)+16.*(st2+66.))-4375.*(2.-3.*st2)*sk14*cos(psi_5))+(2.*(st2-2.)*sk14*st3*(sin(phin_5+psi_1)*kplus+sin(phin_5-psi_1)*kmin)+4.*ct*sk13*st2*(cos(phin_4+psi_1)*params->ccoeff1_5pn[0+17*0]+cos(phin_4-psi_1)*params->ccoeff1_5pn[0+17*1])+16.*ct*sk1*(cos(phin_2+psi_1)*params->ccoeff1_5pn[1+17*0]+cos(phin_2-psi_1)*params->ccoeff1_5pn[1+17*1])+1250.*sk14*st*(105.*st4-126.*st2+28.)*(sin(phin_1+psi_5)*kplus+sin(phin_1-psi_5)*kmin)+625.*(st2-2.)*st3*(sin(phin_5+psi_5)*params->ccoeff1_5pn[2+17*0]+sin(phin_5-psi_5)*params->ccoeff1_5pn[2+17*1])+6.*sk12*st*(sin(phin_3+psi_1)*params->ccoeff1_5pn[3+17*0]+sin(phin_3-psi_1)*params->ccoeff1_5pn[3+17*1])+243.*(st2-2.)*sk12*st3*(sin(phin_5+psi_3)*params->ccoeff0_5pn[1+11*0]+sin(phin_5-psi_3)*params->ccoeff0_5pn[1+11*1])+4.*st*(sin(phin_1+psi_1)*params->ccoeff1_5pn[4+17*0]+sin(phin_1-psi_1)*params->ccoeff1_5pn[4+17*1])+5000.*ct*sk13*(15.*st4-12.*st2+2.)*(cos(phin_2+psi_5)*params->ccoeff0pn[0+2*0]+cos(phin_2-psi_5)*params->ccoeff0pn[0+2*1])-1250.*ct*sk1*st2*(5.*st2-6.)*(cos(phin_4+psi_5)*params->ccoeff1pn[0+15*0]+cos(phin_4-psi_5)*params->ccoeff1pn[0+15*1])+1875.*sk12*st*(8.-22.*st2+15.*st4)*(sin(phin_3+psi_5)*params->ccoeff0_5pn[1+11*0]+sin(phin_3-psi_5)*params->ccoeff0_5pn[1+11*1])+216.*ct*sk1*(cos(phin_2+psi_3)*params->ccoeff1_5pn[5+17*0]+cos(phin_2-psi_3)*params->ccoeff1_5pn[5+17*1])+27.*st*(sin(phin_3+psi_3)*params->ccoeff1_5pn[6+17*0]+sin(phin_3-psi_3)*params->ccoeff1_5pn[6+17*1])+54.*ct*sk1*st2*(cos(phin_4+psi_3)*params->ccoeff1_5pn[7+17*0]+cos(phin_4-psi_3)*params->ccoeff1_5pn[7+17*1])+54.*sk12*st*(sin(phin_1+psi_3)*params->ccoeff1_5pn[8+17*0]+sin(phin_1-psi_3)*params->ccoeff1_5pn[8+17*1])));
	waveampcoeffs[8+11*1]=1./6144.*(192.*ck1*sk1*st2*(sin(psi_1)*(64.-sk12*(7.*st2-6.)+4.*st2)+27.*sin(psi_3)*sk12*(7.*st2-6.))+(4.*sk13*st2*(sin(phin_4+psi_1)*params->ccoeff1_5pn[9+17*0]+sin(phin_4-psi_1)*params->ccoeff1_5pn[9+17*1])-2.*ct*sk14*st3*(cos(phin_5+psi_1)*kplus+cos(phin_5-psi_1)*kmin)-243.*ct*sk12*st3*(cos(phin_5+psi_3)*params->ccoeff0_5pn[1+11*0]+cos(phin_5-psi_3)*params->ccoeff0_5pn[1+11*1])-625.*ct*st3*(cos(phin_5+psi_5)*params->ccoeff1_5pn[2+17*0]+cos(phin_5-psi_5)*params->ccoeff1_5pn[2+17*1])+3.*sk12*s2t*(cos(phin_3+psi_1)*params->ccoeff1_5pn[10+17*0]+cos(phin_3-psi_1)*params->ccoeff1_5pn[10+17*1])+27.*ct*st*(cos(phin_3+psi_3)*params->ccoeff1_5pn[11+17*0]+cos(phin_3-psi_3)*params->ccoeff1_5pn[11+17*1])+1875.*ct*sk12*st*(4.-9.*st2)*(cos(phin_3+psi_5)*params->ccoeff0_5pn[1+11*0]+cos(phin_3-psi_5)*params->ccoeff0_5pn[1+11*1])+2.*s2t*(cos(phin_1+psi_1)*params->ccoeff1_5pn[12+17*0]+cos(phin_1-psi_1)*params->ccoeff1_5pn[12+17*1])+27.*sk12*s2t*(cos(phin_1+psi_3)*params->ccoeff1_5pn[13+17*0]+cos(phin_1-psi_3)*params->ccoeff1_5pn[13+17*1])+8.*sk1*(sin(phin_2+psi_1)*params->ccoeff1_5pn[14+17*0]+sin(phin_2-psi_1)*params->ccoeff1_5pn[14+17*1])+4375.*(2.-3.*st2)*sk14*s2t*(cos(phin_1+psi_5)*kplus+cos(phin_1-psi_5)*kmin)+108.*sk1*(sin(phin_2+psi_3)*params->ccoeff1_5pn[15+17*0]+sin(phin_2-psi_3)*params->ccoeff1_5pn[15+17*1])+162.*sk1*st2*(sin(phin_4+psi_3)*params->ccoeff1_5pn[16+17*0]+sin(phin_4-psi_3)*params->ccoeff1_5pn[16+17*1])-2500.*sk13*(2.-13.*st2+12.*st4)*(sin(phin_2+psi_5)*params->ccoeff0pn[0+2*0]+sin(phin_2-psi_5)*params->ccoeff0pn[0+2*1])-1250.*sk1*st2*(3.-4.*st2)*(sin(phin_4+psi_5)*params->ccoeff1pn[0+15*0]+sin(phin_4-psi_5)*params->ccoeff1pn[0+15*1])));
	waveampcoeffs[9+11*0]=chi1_1/2.*(4.*sk1*(ck1*sk1*cos(phin_2)-c2k1*ct*sin(phin_1)*st+ck1*sk1*st2*(6.*sin(psi_1)*sin(psi_1)-2.+sin(phin_1)*sin(phin_1)))+(2.*ct*sk1*st*((-3.*kplus-4.*sk12)*sin(phin_1+psi_2)+(3.*kmin-4.*sk12)*sin(phin_1-psi_2))+(st2-2.)*((-2.*kplus+(2.*ck1-3.)*sk12)*cos(phin_2+psi_2)+(-2.*kmin+(2.*ck1+3.)*sk12)*cos(phin_2-psi_2))));
	waveampcoeffs[9+11*1]=chi1_1*(-2.*cos(phin_1)*sk1*(st*c2k1+ct*s2k1*sin(phin_1))+(ct*((-2.*kplus+(2.*ck1-3.)*sk12)*sin(phin_2+psi_2)+(-2.*kmin+(2.*ck1+3.)*sk12)*sin(phin_2-psi_2))+sk1*st*((-3.*kplus-4.*sk12)*cos(phin_1+psi_2)+(3.*kmin-4.*sk12)*cos(phin_1-psi_2))));
	waveampcoeffs[10+11*0]=LAL_PI/2.*(6.*sk12*st2*cos(psi_2)+((st2-2.)*(cos(phin_2+psi_2)*params->ccoeff0pn[0+2*0]+cos(phin_2-psi_2)*params->ccoeff0pn[0+2*1])-2.*sk1*s2t*(sin(phin_1+psi_2)*kplus+sin(phin_1-psi_2)*kmin)))	+3.*log(expr[1]/omega0)*(4.*ct*sk1*st*(-kmin*cos(phin_1-psi_2)+kplus*cos(phin_1+psi_2))+(-st2+2.)*(2.*ck1-sk12+2.)*sin(phin_2-psi_2)+(2.*ck1+sk12-2.)*sin(phin_2+psi_2)+6.*sk12*sin(psi_2)*st2);
	waveampcoeffs[10+11*1]=LAL_PI*(ct*(sin(phin_2+psi_2)*params->ccoeff0pn[0+2*0]+sin(phin_2-psi_2)*params->ccoeff0pn[0+2*1])-2.*sk1*st*(cos(phin_1+psi_2)*kplus+cos(phin_1-psi_2)*kmin))+6.*log(expr[1]/omega0)*((2.*ck1-sk12+2.)*ct*cos(phin_2-psi_2)+(2.*ck1+sk12-2.)*ct*cos(phin_2+psi_2)+kmin*2.*sk1*st*sin(phin_1-psi_2)-kplus*2.*sk1*st*sin(phin_1+psi_2));
	// Highest order, only leading order of eps(omega) is needed.
	(*hplus)->data->data[idx]+=ampcoeff*(params->eps*epssqrt*(waveampcoeffs[8+11*0]+waveampcoeffs[9+11*0]+waveampcoeffs[10+11*0]));
	(*hcross)->data->data[idx]+=ampcoeff*(params->eps*epssqrt*(waveampcoeffs[8+11*1]+waveampcoeffs[9+11*1]+waveampcoeffs[10+11*1]));
	case (2):
	waveampcoeffs[4+11*0]=1./48.*(2.*sk12*st2*(5.*sk12*(7.*st2-6.)*(cos(psi_2)-4.*cos(psi_4))-2.*(15.*st2+51.)*cos(psi_2))+(16.*sk13*s2t*(7.*st2-3.)*(sin(phin_1+psi_4)*kplus+sin(phin_1-psi_4)*kmin)-(st2-2.)*sk12*st2*(cos(phin_4+psi_2)*params->ccoeff0pn[0+2*0]+cos(phin_4-psi_2)*params->ccoeff0pn[0+2*1])+4.*(st2-2.)*st2*(cos(phin_4+psi_4)*params->ccoeff1pn[0+15*0]+cos(phin_4-psi_4)*params->ccoeff1pn[0+15*1])+2.*sk1*s2t*((sin(phin_1+psi_2)*params->ccoeff1pn[1+15*0]+sin(phin_1-psi_2)*params->ccoeff1pn[1+15*1])+(sin(phin_3+psi_2)*params->ccoeff1pn[2+15*0]+sin(phin_3-psi_2)*params->ccoeff1pn[2+15*1])-8.*ct2*(sin(phin_3+psi_4)*params->ccoeff0_5pn[1+11*0]+sin(phin_3-psi_4)*params->ccoeff0_5pn[1+11*1]))+2.*(cos(phin_2+psi_2)*params->ccoeff1pn[3+15*0]+cos(phin_2-psi_2)*params->ccoeff1pn[3+15*1])-16.*sk12*(7.*st4-2.*(4.*st2-1.))*(cos(phin_2+psi_4)*params->ccoeff0pn[0+2*0]+cos(phin_2-psi_4)*params->ccoeff0pn[0+2*1])));
	waveampcoeffs[4+11*1]=1./24.*(60.*ck1*ct*sk12*sin(psi_2)*st2+(8.*sk13*st*(7.*st2-6.)*(cos(phin_1+psi_4)*kplus+cos(phin_1-psi_4)*kmin)-ct*sk12*st2*(sin(phin_4+psi_2)*params->ccoeff0pn[0+2*0]+sin(phin_4-psi_2)*params->ccoeff0pn[0+2*1])+2.*sk1*st*(cos(phin_1+psi_2)*params->ccoeff1pn[4+15*0]+cos(phin_1-psi_2)*params->ccoeff1pn[4+15*1])+4.*ct*st2*(sin(phin_4+psi_4)*params->ccoeff1pn[0+15*0]+sin(phin_4-psi_4)*params->ccoeff1pn[0+15*0])+2.*ct*(sin(phin_2+psi_2)*params->ccoeff1pn[5+15*0]+sin(phin_2-psi_2)*params->ccoeff1pn[5+15*1])-8.*ct*sk12*(7.*st2-2.)*(sin(phin_2+psi_4)*params->ccoeff0pn[0+2*0]+sin(phin_2-psi_4)*params->ccoeff0pn[0+2*1])-8.*(2.-3.*st2)*sk1*st*(cos(phin_3+psi_4)*params->ccoeff0_5pn[1+11*0]+cos(phin_3-psi_4)*params->ccoeff0_5pn[1+11*1])+2.*sk1*st*(cos(phin_3+psi_2)*params->ccoeff1pn[6+15*0]+cos(phin_3-psi_2)*params->ccoeff1pn[6+15*1])));
	waveampcoeffs[5+11*0]=1./24.*(ck1*sk1*st2*(4.*(15.*st2+51.)*cos(psi_2)+20.*sk12*(7.*st2-6.)*(4.*cos(psi_4)-cos(psi_2)))+(s2t*(sin(phin_3+psi_2)*params->ccoeff1pn[7+15*0]+sin(phin_3-psi_2)*params->ccoeff1pn[7+15*1])+s2t*(sin(phin_1+psi_2)*params->ccoeff1pn[8+15*0]+sin(phin_1-psi_2)*params->ccoeff1pn[8+15*1])+ct2*s2t*8.*(sin(phin_3+psi_4)*params->ccoeff1pn[9+15*0]+sin(phin_3-psi_4)*params->ccoeff1pn[9+15*1])+2.*sk1*(cos(phin_2+psi_2)*params->ccoeff1pn[10+15*0]+cos(phin_2-psi_2)*params->ccoeff1pn[10+15*1])+8.*sk12*s2t*(7.*st2-3.)*(sin(phin_1+psi_4)*(3.*kplus+4.*sk12)+sin(phin_1-psi_4)*(-3.*kmin+4.*sk12))+16.*sk1*(2.-8.*st2+7.*st4)*(cos(phin_2+psi_4)*params->ccoeff1pn[11+15*0]+cos(phin_2-psi_4)*params->ccoeff1pn[11+15*1])+sk1*st2*(st2-2.)*(cos(phin_4+psi_2)*params->ccoeff1pn[11+15*0]+cos(phin_4-psi_2)*params->ccoeff1pn[11+15*1]-8.*(cos(phin_4+psi_4)*params->ccoeff0_5pn[1+11*0]+cos(phin_4-psi_4)*params->ccoeff0_5pn[1+11*1]))));
	waveampcoeffs[5+11*1]=1./12.*(30.*ct*sk1*sin(psi_2)*st2*(3.*sk12-2.)+(st*(cos(phin_1+psi_2)*params->ccoeff1pn[12+15*0]+cos(phin_1-psi_2)*params->ccoeff1pn[12+15*1])+2.*ct*sk1*((sin(phin_2+psi_2)*params->ccoeff1pn[13+15*0]+sin(phin_2-psi_2)*params->ccoeff1pn[13+15*1])+4.*(7.*st2-2.)*(sin(phin_2+psi_4)*params->ccoeff1pn[11+15*0]+sin(phin_2-psi_4)*params->ccoeff1pn[11+15*1]))+ct*sk1*st2*((sin(phin_4+psi_2)*params->ccoeff1pn[11+15*0]+sin(phin_4-psi_2)*params->ccoeff1pn[11+15*1])-8.*(sin(phin_4+psi_4)*params->ccoeff0_5pn[1+11*0]+sin(phin_4-psi_4)*params->ccoeff0_5pn[1+11*1]))+st*((cos(phin_3+psi_2)*params->ccoeff1pn[14+15*0]+cos(phin_3-psi_2)*params->ccoeff1pn[14+15*1])+4.*(2.-3.*st2)*(cos(phin_3+psi_4)*params->ccoeff1pn[9+15*0]+cos(phin_3-psi_4)*params->ccoeff1pn[9+15*1]))+4.*sk12*st*(6.-7.*st2)*((-3.*kplus-4.*sk12)*cos(phin_1+psi_4)+(+3.*kmin-4.*sk12)*cos(phin_1-psi_4))));
	waveampcoeffs[6+11*0]=chi1_1/2.*st*(kplus*sin(phin_1+psi_1)-kmin*sin(phin_1-psi_1));
	waveampcoeffs[6+11*1]=chi1_1/2.*st*(2.*sk1*sin(psi_1)*st+ct*(kplus*cos(phin_1+psi_1)-kmin*cos(phin_1-psi_1)));
	waveampcoeffs[7+11*0]=chi1_1/2.*((-kplus*ct*sin(phin_2+psi_1-phi1_1)+kmin*ct*sin(phin_2-psi_1-phi1_1))+sk1*st*(sin(phin_1+psi_1)-sin(phin_1-psi_1)
+cos(phin_1+psi_1-phi1_1)-cos(phin_1-psi_1-phi1_1)));
	waveampcoeffs[7+11*1]=chi1_1/4.*(-st2*(2.*(sin(phi1_1)+2.)*ck1*sin(psi_1)+2.*cos(phi1_1)*cos(psi_1))+(2.*ct*sk1*st*(cos(phin_1+psi_1)-cos(phin_1-psi_1)-sin(phin_1+psi_1-phi1_1)+sin(phin_1-psi_1-phi1_1))+(st2-2.)*(kplus*cos(phin_2+psi_1-phi1_1)-kmin*cos(phin_2-psi_1-phi1_1))));
	// Since highest order is 1.5 PN and there is no 0.5 PN order correction to eps(omega), leading order eps is enough.
	(*hplus)->data->data[idx]+=ampcoeff*(params->eps*(waveampcoeffs[4+11*0]+waveampcoeffs[2+11*0]*4.*params->xi+params->beta1*waveampcoeffs[5+11*0]+waveampcoeffs[6+11*0]+waveampcoeffs[7+11*0]*params->beta1));
	(*hcross)->data->data[idx]+=ampcoeff*(params->eps*(waveampcoeffs[4+11*1]+waveampcoeffs[2+11*1]*4.*params->xi+params->beta1*waveampcoeffs[5+11*1]+waveampcoeffs[6+11*1]+waveampcoeffs[7+11*1]*params->beta1));
	case (1):
	waveampcoeffs[2+11*0]=1./64.*(4.*ct*sk1*st2*(-45.*sk12*cos(psi_3)+cos(psi_1)*(5.*sk12-4.))-sk12*st*((st2-2.)*(sin(phin_3+psi_1)*kplus+sin(phin_3-psi_1)*kmin)-45.*(2.-3.*st2)*(sin(phin_1+psi_3)*kplus+sin(phin_1-psi_3)*kmin))+st*(sin(phin_1+psi_1)*params->ccoeff0_5pn[0+11*0]+sin(phin_1-psi_1)*params->ccoeff0_5pn[0+11*1]-9.*(st2-2.)*(sin(phin_3+psi_3)*params->ccoeff0_5pn[1+11*0]+sin(phin_3-psi_3)*params->ccoeff0_5pn[1+11*1]))+2.*ct*sk1*(cos(phin_2+psi_1)*params->ccoeff0_5pn[2+11*0]+cos(phin_2-psi_1)*params->ccoeff0_5pn[2+11*1]+9.*(2-3.*st2)*(cos(phin_2+psi_3)*params->ccoeff0pn[0+2*0]+cos(phin_2-psi_3)*params->ccoeff0pn[0+2*1])));
	waveampcoeffs[2+11*1]=1./32.*(-16.*s2k1*st2*sin(psi_1)+ct*st*sk12*(cos(phin_3+psi_1)*kplus+cos(phin_3-psi_1)*kmin+45.*(cos(phin_1+psi_3)*kplus+cos(phin_1-psi_3)*kmin))+0.5*s2t*(cos(phin_1+psi_1)*params->ccoeff0_5pn[3+11*0]+cos(phin_1-psi_1)*params->ccoeff0_5pn[3+11*1]+9.*(cos(phin_3+psi_3)*params->ccoeff0_5pn[1+11*0]+cos(phin_3-psi_3)*params->ccoeff0_5pn[1+11*1]))+2.*sk1*(sin(phin_2+psi_1)*params->ccoeff0_5pn[4+11*0]+sin(phin_2-psi_1)*params->ccoeff0_5pn[4+11*1]-9.*c2t*(sin(phin_2+psi_3)*params->ccoeff0pn[0+2*0]+sin(phin_2-psi_3)*params->ccoeff0pn[0+2*1])));
	waveampcoeffs[3+11*0]=1./64.*(4.*ck1*ct*st2*(135.*cos(psi_3)*sk12+cos(psi_1)*(4.-15.*sk12))+(2.*ct*(9.*(2.-3.*st2)*(cos(phin_2+psi_3)*params->ccoeff0_5pn[5+11*0]+cos(phin_2-psi_3)*params->ccoeff0_5pn[5+11*1])+cos(phin_2+psi_1)*params->ccoeff0_5pn[6+11*0]+cos(phin_2-psi_1)*params->ccoeff0_5pn[6+11*1])-sk1*st*(st2-2.)*(27.*(sin(phin_3+psi_3)*params->ccoeff0pn[0+2*0]+sin(phin_3-psi_3)*params->ccoeff0pn[0+2*0])+sin(phin_3+psi_1)*params->ccoeff0_5pn[7+11*0]+sin(phin_3-psi_1)*params->ccoeff0_5pn[7+11*1])+sk1*st*(45.*(2.-3.*st2)*(sin(phin_1+psi_3)*params->ccoeff0_5pn[7+11*0]+sin(phin_1-psi_3)*params->ccoeff0_5pn[7+11*1])+sin(phin_1+psi_1)*params->ccoeff0_5pn[8+11*0]+sin(phin_1-psi_1)*params->ccoeff0_5pn[8+11*1])));
	waveampcoeffs[3+11*1]=1./32.*(32.*sin(psi_1)*st2*c2k1+(ct*sk1*st*(27.*(cos(phin_3+psi_3)*params->ccoeff0pn[0+2*0]+cos(phin_3-psi_3)*params->ccoeff0pn[0+2*1])+(cos(phin_3+psi_1)+45.*cos(phin_1+psi_3))*params->ccoeff0_5pn[7+11*0]+(cos(phin_3-psi_1)+45*cos(phin_1-psi_3))*params->ccoeff0_5pn[7+11*1]+(params->ccoeff0_5pn[9+11*0]*cos(phin_1+psi_1)+params->ccoeff0_5pn[9+11*1]*cos(phin_1-psi_1)))-(18.*c2t)*(sin(phin_2+psi_3)*params->ccoeff0_5pn[5+11*0]+sin(phin_2-psi_3)*params->ccoeff0_5pn[5+11*1])-2.*(sin(phin_2+psi_1)*params->ccoeff0_5pn[10+11*0]+sin(phin_2-psi_1)*params->ccoeff0_5pn[10+11*1])));
	//The 0.5 PN correction needs to include 1 PN correction of eps(omega) the amplitude is taken to 1.5 PN order
	if (params->pnamp == 3){
	(*hplus)->data->data[idx]+=ampcoeff/(params->eps*epssqrt)*(eps0_5pn1pncorr*eps0_5pn1pncorr*(waveampcoeffs[2+11*0]+params->beta1*waveampcoeffs[3+11*0]-2.*params->xi*waveampcoeffs[0+11*0]));
	(*hcross)->data->data[idx]+=ampcoeff/(params->eps*epssqrt)*(eps0_5pn1pncorr*eps0_5pn1pncorr*(waveampcoeffs[2+11*1]+params->beta1*waveampcoeffs[3+11*1]-2.*params->xi*waveampcoeffs[0+11*1]));
	}
	else{
	(*hplus)->data->data[idx]+=ampcoeff*(epssqrt*(waveampcoeffs[2+11*0]+params->beta1*waveampcoeffs[3+11*0]-2.*params->xi*waveampcoeffs[0+11*0]));
	(*hcross)->data->data[idx]+=ampcoeff*(epssqrt*(waveampcoeffs[2+11*1]+params->beta1*waveampcoeffs[3+11*1]-2.*params->xi*waveampcoeffs[0+11*1]));
	}
	case (0):
	waveampcoeffs[0+11*0]=0.25*((st2-2.)*(cos(phin_2+psi_2)*params->ccoeff0pn[0+2*0]+cos(phin_2-psi_2)*params->ccoeff0pn[0+2*1])-2.*sk1*s2t*(sin(phin_1+psi_2)*kplus+sin(phin_1-psi_2)*kmin)+6.*sk12*st2*cos(psi_2));
	waveampcoeffs[0+11*1]=0.5*(ct*(sin(phin_2+psi_2)*params->ccoeff0pn[0+2*0]+sin(phin_2-psi_2)*params->ccoeff0pn[0+2*1])-2.*st*sk1*(cos(phin_1+psi_2)*kplus+cos(phin_1-psi_2)*kmin));
	waveampcoeffs[1+11*0]=0.5*(s2t*(sin(phin_1+psi_2)*params->ccoeff0pn[1+2*0]+sin(phin_1-psi_2)*params->ccoeff0pn[1+2*1])+sk1*(st2-2.)*(cos(phin_2+psi_2)*kplus+cos(phin_2-psi_2)*kmin)-3.*s2k1*st2*cos(psi_2));
	waveampcoeffs[1+11*1]=ct*sk1*(sin(phin_2+psi_2)*kplus+sin(phin_2-psi_2)*kmin)+st*(cos(phin_1+psi_2)*params->ccoeff0pn[1+2*0]+cos(phin_1-psi_2)*params->ccoeff0pn[1+2*1]);
	// If the amplitude is taken to 1 PN order, the eps(omega) needs to include the 1 PN correction, if amplitude is taken to 1.5 PN order, that eps(omega) needs to include corrections up to 1.5 PN order
	if (params->pnamp == 3){
	(*hplus)->data->data[idx]+=ampcoeff/(params->eps*epssqrt)*eps0pn1_5pncorr*sqrt(eps0pn1_5pncorr)*(waveampcoeffs[0+11*0]+params->beta1*waveampcoeffs[1+11*0]);
	(*hcross)->data->data[idx]+=ampcoeff/(params->eps*epssqrt)*eps0pn1_5pncorr*sqrt(eps0pn1_5pncorr)*(waveampcoeffs[0+11*1]+params->beta1*waveampcoeffs[1+11*1]);
	}
	else if (params->pnamp == 2){
	(*hplus)->data->data[idx]+=ampcoeff/(params->eps*epssqrt)*eps0pn1_5pncorr*sqrt(eps0pn1pncorr)*(waveampcoeffs[0+11*0]+params->beta1*waveampcoeffs[1+11*0]);
	(*hcross)->data->data[idx]+=ampcoeff/(params->eps*epssqrt)*eps0pn1_5pncorr*sqrt(eps0pn1pncorr)*(waveampcoeffs[0+11*1]+params->beta1*waveampcoeffs[1+11*1]);
	}
	else{
	(*hplus)->data->data[idx]+=ampcoeff*(waveampcoeffs[0+11*0]+params->beta1*waveampcoeffs[1+11*0]);
	(*hcross)->data->data[idx]+=ampcoeff*(waveampcoeffs[0+11*1]+params->beta1*waveampcoeffs[1+11*1]);
	}
	break;
}
XLALFreeDmatrix(waveampcoeffs);
return XLAL_SUCCESS;
}

/**
 * Interface routine, calculating the prefered variables for the Spin-dominated waveforms
 */
int XLALSimInspiralSpinDominatedWaveformInterfaceTD(
	REAL8TimeSeries **hplus,        /**< +-polarization waveform */
	REAL8TimeSeries **hcross,       /**< x-polarization waveform */
	REAL8 deltaT,                   /**< sampling interval (s) */
	REAL8 m1,                       /**< mass of companion 1 (kg) */
	REAL8 m2,                       /**< mass of companion 2 (kg) */
	REAL8 fStart,                   /**< start GW frequency (Hz) */
	REAL8 fRef,			/**< end GW frequency (Hz) */
	REAL8 D,                        /**< distance of source (m) */
	REAL8 s1x,                      /**< initial value of S1x */
	REAL8 s1y,                      /**< initial value of S1y */
	REAL8 s1z,                      /**< initial value of S1z */
	REAL8 lnhatx,                   /**< initial value of LNhatx */
	REAL8 lnhaty,                   /**< initial value of LNhaty */
	REAL8 lnhatz,                   /**< initial value of LNhatz */
	int phaseO,                     /**< twice PN phase order */
	int amplitudeO,                 /**< twice PN amplitude order */
	REAL8 phiRef			/**< Reference phase at the Reference Frequency */
)
{
REAL8 totalmass, nu, chi1, beta1, kappa1, totalJ, S1, Jx, Jy, Jz, omega, eta, romega, v, LN, theta;
REAL8 alpha0, phin0, JxNx, JxNy, JxNz, inertiaxx, inertiaxy, inertiaxz, JxLx, JxLy, JxLz, LXYx, LXYy, LXYz;
int n;
totalmass=m1+m2;
if (m1 > m2) nu=m2/m1; else nu=m1/m2;
if (LAL_SDW_MAX_PN_PARAM < 100.*nu*nu) {
XLALPrintError("XLAL Error: Spin-dominated waveforms error: Please make sure that the total mass is higher than 45 solar mass, and mass ratio is lower than 0.03125. Also above 130 solar mass be aware that high starting frequency may result in termination right after start, due to high value of the pn parameter. \n");
XLAL_ERROR(XLAL_EDOM);
} //too high mass ratio for the waveform, abort
omega = fStart * LAL_PI;
eta = nu/(1.+nu)/(1.+nu);
chi1 = sqrt(s1x*s1x+s1y*s1y+s1z*s1z);
if (chi1 < 0.5) {
XLALPrintError("XLAL Error: Spin-dominated waveforms error: Please make sure that the dimensionless spin parameter is higher than 0.5 \n");
XLAL_ERROR(XLAL_EDOM);
}
kappa1 = acos((lnhatx*s1x+lnhaty*s1y+lnhatz*s1z)/chi1);
// Calculate the orbital angular momentum, up to 1.5 PN, with SO corrections
v=cbrt(LAL_G_SI*totalmass*omega/LAL_C_SI/LAL_C_SI/LAL_C_SI);
romega=LAL_G_SI*totalmass/LAL_C_SI/LAL_C_SI/v/v*(1.); //-(3.-eta)*v*v/3.-v*v*v/3.*(2./(1.+nu)/(1.+nu)+3.*eta)-v*v*v*v*(-eta*(19./4.+eta/9.)+1./4.*(-chi1*chi1*eta*nu)*(3.*cos(kappa1)*cos(kappa1)-1.)));
LN = eta * totalmass * romega * romega * omega;
// Calculate Spin magnitude, and the total angular momentum J
S1 = chi1 * LAL_G_SI / LAL_C_SI * totalmass * totalmass * eta / nu;
Jx = LN*lnhatx + S1 * s1x/chi1;
Jy = LN*lnhaty + S1 * s1y/chi1;
Jz = LN*lnhatz + S1 * s1z/chi1;
totalJ = sqrt(Jx*Jx + Jy*Jy + Jz*Jz);
// calculate the remaining angles
theta = acos(Jz/totalJ);
if (kappa1 == 0) {
phin0=3.*LAL_PI_2;
beta1 = 0;
} else {
beta1 = acos((Jx*s1x+Jy*s1y+Jz*s1z)/totalJ/chi1);
// calculating the initial value of the \phi_n variable
JxNx = Jx/totalJ*1.0;
JxNy = - Jx/totalJ*1.0;
JxNz = 0.;
inertiaxx = JxNy*Jz/totalJ-JxNz*Jy/totalJ;
inertiaxy = JxNz*Jx/totalJ-JxNx*Jz/totalJ;
inertiaxz = JxNx*Jy/totalJ-JxNy*Jx/totalJ;
JxLx = Jy/totalJ*lnhatz - Jz/totalJ*lnhaty;
JxLy = Jz/totalJ*lnhatx - Jx/totalJ*lnhatz;
JxLz = Jx/totalJ*lnhaty - Jy/totalJ*lnhatx;
LXYx = JxLy*Jz/totalJ-JxLz*Jy/totalJ;
LXYy = JxLz*Jx/totalJ-JxLy*Jz/totalJ;
LXYz = JxLx*Jy/totalJ-JxLy*Jx/totalJ;
alpha0 = acos(inertiaxx*LXYx + inertiaxy*LXYy + inertiaxz*LXYz);
phin0=3.*LAL_PI_2-alpha0;
}
// calling the SDW driver with the prefered variables
n = XLALSimInspiralSpinDominatedWaveformDriver(hplus, hcross, totalmass, nu, chi1, D, kappa1, beta1, theta, fStart, fRef, phaseO, amplitudeO, deltaT, phiRef, phin0);
return n;
}

/**
 * Function calculating the Spin-Dominated waveforms
 * This waveform is an inspiral only, 1 spin, precessing waveform.
 * For the formulae see the appendix of Arxiv:1209.1722
 */
int XLALSimInspiralSpinDominatedWaveformDriver(
	REAL8TimeSeries **hplus,        /**< +-polarization waveform */
	REAL8TimeSeries **hcross,       /**< x-polarization waveform */
	REAL8 totalmass,		/**< total mass of the binary */
	REAL8 nu,			/**< mass ratio */
	REAL8 chi1,			/**< dimensionless spin paramter */
	REAL8 D,			/**< Distance to the source */
	REAL8 kappa1,			/**< Angle span by S_1 and L */
	REAL8 beta1,			/**< Angle span by J and S_1 */
	REAL8 theta,			/**< Angle span by the line of sight and J */
	REAL8 fStart,			/**< Starting gravitational wave frequency*/
	REAL8 fRef,			/**< Ending gravitational wave frequency*/
	int phaseO,                     /**< twice PN phase order */
	int amplitudeO,                 /**< twice PN amplitude order */
	REAL8 deltaT,			/**< Sampling time interval */
	REAL8 phiRef,			/**< Reference phase at the Reference Frequency */
	REAL8 phin0			/**< Starting value of the \phi_n parameter */
)
{
int idx;
int n;
unsigned int i;
REAL8 phiShift;
LIGOTimeGPS tStart=LIGOTIMEGPSZERO;
/* check inputs for sanity */
if (*hplus) XLAL_ERROR(XLAL_EFAULT);
if (*hcross) XLAL_ERROR(XLAL_EFAULT);
if (deltaT <= 0) XLAL_ERROR(XLAL_EDOM);
if (totalmass < 0) XLAL_ERROR(XLAL_EDOM);
if (fStart <= 0) XLAL_ERROR(XLAL_EDOM);
/* set up the integrator*/
ark4GSLIntegrator *integrator= XLALAdaptiveRungeKutta4Init(LAL_SDW_NUM_VARIABLES,XLALSpinDominatedWaveformDerivatives,XLALSpinDominatedWaveformStoppingTest,LAL_SDW_ABSOLUTE_TOLERANCE,LAL_SDW_RELATIVE_TOLERANCE);
if (!integrator) {
	XLALPrintError("XLAL Error - %s: Cannot allocate integrator\n", __func__);
	XLAL_ERROR(XLAL_EFUNC);
}
/* stop the integration only when the test is true */
integrator->stopontestonly = 1;
LALSDWaveformParams params;
params.totalmass=totalmass;
params.nu=nu;
params.chi1=chi1;
params.dist=D;
params.kappa1=kappa1;
params.beta1=beta1;
params.theta=theta;
params.eps=0.;
params.xi=0.;
params.pnamp=amplitudeO;
params.pnphase=phaseO;
params.prevdomega=0.;
n = XLALSpinDominatedWaveformConstantCoefficients(&params);
if( n < 0 ) XLAL_ERROR(XLAL_EFUNC);
REAL8 yin[LAL_SDW_NUM_VARIABLES];
yin[0]=phin0;
yin[1]=fStart*LAL_PI;
yin[2]=0.;
REAL8Array *yout;
// estimating the length of the waveform
REAL8 length=5./256.*pow(fStart*LAL_PI,-8./3.)*(1+params.nu)*(1+params.nu)/params.nu*pow(LAL_G_SI*params.totalmass/LAL_C_SI/LAL_C_SI/LAL_C_SI,-5./3.);
INT4 intLen    = XLALAdaptiveRungeKutta4Hermite(integrator,(void *)&params,yin,0.0,length,deltaT,&yout);
UNUSED INT4 intReturn = integrator->returncode;
XLALAdaptiveRungeKutta4Free(integrator);
REAL8TimeSeries *phin  = XLALCreateREAL8TimeSeries( "PHI_N", &tStart, 0., deltaT, &lalDimensionlessUnit, intLen);
REAL8TimeSeries *omega = XLALCreateREAL8TimeSeries( "OMEGA", &tStart, 0., deltaT, &lalDimensionlessUnit, intLen);
REAL8TimeSeries *psi   = XLALCreateREAL8TimeSeries( "ORBITAL_PHASE", &tStart, 0., deltaT, &lalDimensionlessUnit, intLen);
for (idx=0;idx<intLen;idx++) {
	phin->data->data[idx]    = yout->data[intLen+idx];
	omega->data->data[idx]  = yout->data[2*intLen+idx];
	psi->data->data[idx]  = yout->data[3*intLen+idx];
}

if (fRef == 0){
	phiShift = phiRef - psi->data->data[0];
	for( i=0; i < psi->data->length; i++){
		psi->data->data[i] += phiShift;
        }
} else if (fRef == fStart) {
	phiShift = phiRef - psi->data->data[0];
	for( i=0; i < psi->data->length; i++){
		psi->data->data[i] += phiShift;
	}
} else {
	XLALPrintError("XLAL Error: Spin-dominated waveforms error: Please set the reference frequency as the starting frequency, Setting 0 will default to the starting frequency. \n");
	XLAL_ERROR(XLAL_EDOM);
}
if ((*hplus) && (*hcross)) {
	if ((*hplus)->data->length!=(*hcross)->data->length) {
	XLALPrintError("***  h+ and hx differ in length\n");
	XLAL_ERROR(XLAL_EFAILED);
} else {
	if ((int)(*hplus)->data->length<intLen) {
	XLALPrintError("*** ERROR: h+ and hx too short\n");
	XLAL_ERROR(XLAL_EFAILED);
} else {
	XLALGPSAdd(&((*hplus)->epoch),-intLen*deltaT);
	XLALGPSAdd(&((*hcross)->epoch),-intLen*deltaT);
}
	}
} else {
	XLALGPSAdd(&tStart,-intLen*deltaT);
	*hplus  = XLALCreateREAL8TimeSeries("H+", &tStart, 0.0, deltaT, &lalStrainUnit, intLen);
	*hcross = XLALCreateREAL8TimeSeries("Hx", &tStart, 0.0, deltaT, &lalStrainUnit, intLen);
	if(*hplus == NULL || *hcross == NULL)
	XLAL_ERROR(XLAL_ENOMEM);
}

REAL8 expr[3];
for (idx=0;idx<intLen;idx++) {
	expr[0]=phin->data->data[idx];
	expr[1]=omega->data->data[idx];
	expr[2]=psi->data->data[idx];
	n = XLALSpinDominatedWaveformBuild(&params, expr,hplus,hcross,idx);
	if( n < 0 ) XLAL_ERROR(XLAL_EFUNC);
}
XLALDestroyREAL8Array(yout);
XLALDestroyREAL8TimeSeries(phin);
XLALDestroyREAL8TimeSeries(omega);
XLALDestroyREAL8TimeSeries(psi);
return intLen;
}

/**
 * Function calculating the derivatives of the three time dependent variables of the Spin-Dominated waveforms (SDW)
 * The first paramter is \phi_n, Eq 27 of Arxiv:1005.5330, taken for 1 spin case, and integrated over an orbital period.
 * The second parameter is \omega, the derivative is taken from Arxiv: astro-ph/0504538, up to 2 PN orders with 1 spin. (In order to stay consistent with SDW)
 * The thirs parameter is the phase.
 */
static INT4 XLALSpinDominatedWaveformDerivatives(UNUSED REAL8 t,
					const REAL8 values[],
					REAL8 dvalues[],
					void *mparams)
{
  LALSDWaveformParams *params = (LALSDWaveformParams *) mparams;
// parameters required for the time derivatives
REAL8 v=cbrt(LAL_G_SI*params->totalmass*values[1]/LAL_C_SI/LAL_C_SI/LAL_C_SI);
params->eps=v*v;
params->xi=params->nu/sqrt(params->eps); // second small parameter
REAL8 phasecoeff=96./5.*params->nu/(1.+ params->nu)/(1.+ params->nu)*v*v*v*v*v*values[1]*values[1];
// Calculating the derivatives

dvalues[0]=0;
if (params->kappa1 == 0) {
dvalues[0]=0;
} else {
switch (params->pnphase) {
	case 4:
	dvalues[0]+=+3./2./LAL_G_SI/params->totalmass/sin(params->kappa1)*params->chi1*params->chi1*LAL_C_SI*LAL_C_SI*LAL_C_SI*params->eps*params->eps*params->eps*sqrt(params->eps)*((1.-2.*params->xi*sqrt(params->eps))*sin(params->kappa1)*cos(params->kappa1)+params->beta1*cos(params->kappa1)*cos(params->kappa1));
	case 3:
	dvalues[0]+=params->chi1*LAL_C_SI*LAL_C_SI*LAL_C_SI*params->eps*params->eps*params->eps/2./LAL_G_SI/params->totalmass*(-4.*cos(params->kappa1)*params->beta1+5*sqrt(params->eps)*params->xi-4.);
	case 2:
	case 1:
	case 0:
	break;
}
}

dvalues[1]=0;
switch (params->pnphase) {
	case 4:
	dvalues[1]+=phasecoeff*v*v*v*v*(34103./18144.+13661./2016.*params->nu/(1.+ params->nu)/(1.+ params->nu)+59./18.*params->nu/(1.+ params->nu)/(1.+ params->nu)*params->nu/(1.+ params->nu)/(1.+ params->nu));
	// QM and Self-Spin components taken together
	dvalues[1]+=phasecoeff*v*v*v*v*(5./2.*params->chi1*params->chi1/(1.+ params->nu)/(1.+ params->nu)*(3.*cos(params->kappa1)*cos(params->kappa1)-1.)+1./96.*params->chi1*params->chi1/(1.+ params->nu)/(1.+ params->nu)*(6.+sin(params->kappa1)*sin(params->kappa1)));
	case 3:
	dvalues[1]+=phasecoeff*v*v*v*(4.*LAL_PI);
	// SO component
	dvalues[1]+=phasecoeff*v*v*v*(-1/12.*cos(params->kappa1)*params->chi1*(113./(1.+ params->nu)/(1.+ params->nu)+75.*params->nu/(1.+ params->nu)/(1.+ params->nu)));
	case 2:
	dvalues[1]+=-phasecoeff*v*v*(743./336.+11./4.*params->nu/(1.+ params->nu)/(1.+ params->nu));
	case 1:
	case 0:
	dvalues[1]+=phasecoeff*(1.);
	break;
}
dvalues[2]=values[1];

return GSL_SUCCESS;
} /* end of XLALSpinDominatedWaveformDerivatives */

/**
* Stopping test for the Spin-Dominated waveforms. Using MECO, or the desired ending frequency. The maximum value of the PN parameter is set to 0.8.
*/
static INT4 XLALSpinDominatedWaveformStoppingTest(UNUSED REAL8 t, const REAL8 values[], REAL8 dvalues[], UNUSED void *mparams) {
LALSDWaveformParams *params = (LALSDWaveformParams *) mparams;
REAL8 v=cbrt(LAL_G_SI*params->totalmass*values[1]/LAL_C_SI/LAL_C_SI/LAL_C_SI);
REAL8 eta = params->nu/(1.+ params->nu)/(1.+ params->nu);
REAL8 omega  =  values[1];
REAL8 domega =  dvalues[1];
REAL8 d2omega = dvalues[1]-params->prevdomega;
params->prevdomega = dvalues[1];
REAL8 mecotest;
mecotest=0.;

switch (params->pnphase) {
	case 4:
	mecotest+=6.*v*v*v*v*(1./8.*(-27.+19.*eta-eta*eta/3.)-(3.*cos(params->kappa1)*cos(params->kappa1)-1)/2.*params->chi1*params->chi1*eta/params->nu);
	case 3:
	mecotest+=5.*v*v*v*(8./3.*eta/params->nu+2.*eta)*cos(params->kappa1)*params->chi1;
	case 2:
	mecotest+=-4.* v * v * (3. + eta/3.)/4.;
	case 1:
	case 0:
	mecotest+=2.;
	break;
}

if (mecotest<0) {
	XLALPrintWarning("** LALSimInspiralSDW WARNING **: MECO reached\n");
	return -1;
} else if (isnan(omega)) {
	XLALPrintWarning("** LALSimInspiralSDW WARNING **: omega is NAN\n");
	return -1;
} else if (v >= 1.) {
	XLALPrintWarning("** LALSimInspiralSDW WARNING **: PN parameter is too large\n");
	return -1;
} else if (domega < 0.0) {
	XLALPrintWarning("** LALSimInspiralSDW WARNING **: domega < 0\n");
	return -1;
} else if (d2omega <= 0.){
	XLALPrintWarning("** LALSimInspiralSDW WARNING **: d2omega < 0\n");
	return -1;
} else
	return GSL_SUCCESS;
} /* End of XLALSpinDominatedWaveformStoppingTest */
