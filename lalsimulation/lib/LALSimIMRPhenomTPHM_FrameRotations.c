/*
 * Copyright (C) 2020 Hector Estelles
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


/**
 * \author Hector Estelles
 */

#include <math.h>

/* LAL Header Files */
#include <lal/LALSimIMR.h>
#include <lal/LALConstants.h>
#include <lal/Date.h>
#include <lal/Units.h>

/* GSL Header Files */
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_roots.h>

#include <lal/XLALGSL.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>


#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_gamma.h>

#include <lal/LALSimInspiralPrecess.h>


/* Routines to perform frame rotations for a set of Spin-Weighted Spherical Harmonic modes. */

/****************************************/
/********* FUNCTION DEFINITION **********/
/****************************************/

/* Struct for storing recurrent squared roots in the Wigner coefficients */
typedef struct tagPhenomT_precomputed_sqrt
{
  REAL8 sqrt2, sqrt2half, sqrt3, sqrt5, sqrt6, sqrt7, sqrt10, sqrt14, sqrt15, sqrt21, sqrt30, sqrt35, sqrt70, sqrt210;
} PhenomT_precomputed_sqrt;

/* Struct to store Wignerd coefficients and powers of the exponentials of the precessing angles */
typedef struct tagPhenomTPWignerStruct
{
  REAL8 wignerdL2[5][5];
  REAL8 wignerdL3[7][7];
  REAL8 wignerdL4[9][9];
  REAL8 wignerdL5[11][11];

  COMPLEX16 expAlphaPowers[11];
  COMPLEX16 expGammaPowers[11];
} PhenomTPWignerStruct;

void IMRPhenomTPHM_SetPrecomputedSqrt(PhenomT_precomputed_sqrt *SQRT);

/* Function to set the PhenomTPWignerStruct given some Euler angles.
Sign of beta has to be passed independently since this information is lost in cosBeta.
globalRot flag specifies how many coefficients have to be computed, since for the time dependent 
rotation between the co-precessing and the J-frame not all are needed. */
void IMRPhenomTPHM_SetWignerDStruct(  PhenomTPWignerStruct *wS, /**< Wigner struct */
                                      PhenomT_precomputed_sqrt *SQRT, /**< Precomputed squared roots */
									                    REAL8 cosBeta,			  /**< cosinus of opening angle beta */
									                    REAL8 alpha,			  /**< precessing angle alpha */
									                    REAL8 gamma,			  /**< precessing angle gamma */
                                      INT4 LMAX,          /**< Maximum l mode number */
									                    INT4 sign,				  /**< Sign of beta */
									                    UINT4 globalRot			  /**< Compute coefficients for global or local rotation */
									                  );

/* Function to compute the WignerDMatrix. Analogous to XLALWignerDMatrix (which employs XLALWignerdMatrix), but Wignerd coefficients
here are read from the precomputed struct, in order to improve computational speed. */
COMPLEX16 PhenomTWignerDMatrix(
                                   int l,        /**< mode number l */
                                   int mp,        /**< mode number m' */
                                   int m,        /**< mode number m */
                                   PhenomTPWignerStruct *wS  /**< Wigner Struct, which contains Wignerd elements and powers of alpha and gamma */
    );

/* Function to rotate a given set of Spin-Weighted Spherical Harmonic mode time series with an time-dependent Euler rotation specified by three Euler angles.*/
int PhenomTPHM_RotateModes(
                SphHarmTimeSeries* h_lm, /**< spherical harmonic decomposed modes, modified in place */
                REAL8TimeSeries* alpha, /**< alpha Euler angle time series */
                REAL8TimeSeries* cosbeta, /**< beta Euler angle time series */
                REAL8TimeSeries* gam, /**< gamma Euler angle time series */
                PhenomT_precomputed_sqrt *SQRT /**< precomputed squared root factors */
);

/* Function to rotate a given set of Spin-Weighted Spherical Harmonic mode time series with an global Euler rotation specified by three Euler angles.*/
int PhenomTPHM_RotateModes_Global(
                SphHarmTimeSeries* h_lm, /**< spherical harmonic decomposed modes, modified in place */
                REAL8 alpha, /**< alpha Euler angle time series */
                REAL8 cosbeta, /**< beta Euler angle time series */
                REAL8 gam, /**< gamma Euler angle time series */
				        size_t length, /**< Length of waveform modes */
                PhenomT_precomputed_sqrt *SQRT /**< precomputed squared root factors */
);

/*********************************************/
/********* WIGNER-D MATRIX ROUTINES **********/
/*********************************************/

/* Function to precompute and store squared root factors for the Wigner coefficients */
void IMRPhenomTPHM_SetPrecomputedSqrt(PhenomT_precomputed_sqrt *SQRT)
{
  SQRT->sqrt2 = sqrt(2.0);
  SQRT->sqrt2half = sqrt(2.5);
  SQRT->sqrt3 = sqrt(3.0);
  SQRT->sqrt5 = sqrt(5.0);
  SQRT->sqrt6 = sqrt(6.0);
  SQRT->sqrt7 = sqrt(7.0);
  SQRT->sqrt10 = sqrt(10.0);
  SQRT->sqrt14 = sqrt(14.0);
  SQRT->sqrt15 = sqrt(15.0);
  SQRT->sqrt21 = sqrt(21.0);
  SQRT->sqrt30 = sqrt(30.0);
  SQRT->sqrt35 = sqrt(35.0);
  SQRT->sqrt70 = sqrt(70.0);
  SQRT->sqrt210 = sqrt(210.0);
}

/* Function to compute the WignerDMatrix. Analogous to XLALWignerDMatrix (which employs XLALWignerdMatrix), but Wignerd coefficients
here are read from the precomputed struct, in order to improve computational speed.

Format of the Winerd elements corresponds to dlm[idx] with idx from l-m' to l+m' and the corresponding element is given
the function WignerD[{l,m',m},0,-beta,0] in Mathematica 12. See eq. 5 of https://dcc.ligo.org/DocDB/0159/T1900080/008/STmodes.pdf for
the relation between Mathematica's WignerD and XLALWignerDMatrix.  */
void IMRPhenomTPHM_SetWignerDStruct(PhenomTPWignerStruct *wS, PhenomT_precomputed_sqrt *SQRT, REAL8 cosBeta, REAL8 alpha, REAL8 gamma, INT4 LMAX, INT4 sign, UINT4 globalRot)
{
	/*cos(beta/2) and sin(beta/2) powers*/

  REAL8 cBetah = sqrt(0.5*fabs(1 + cosBeta));
  REAL8 sBetah = sign*sqrt(0.5*fabs(1 - cosBeta));

  REAL8 cBetah2 = cBetah * cBetah;
  REAL8 cBetah3 = cBetah * cBetah2;
  REAL8 cBetah4 = cBetah * cBetah3;

  REAL8 sBetah2 = sBetah * sBetah;
  REAL8 sBetah3 = sBetah * sBetah2;
  REAL8 sBetah4 = sBetah * sBetah3;

  REAL8 C2mS2 = cBetah2 - sBetah2;

  REAL8 cBetah5, cBetah6, cBetah7, cBetah8, cBetah9, cBetah10;
  REAL8 sBetah5, sBetah6, sBetah7, sBetah8, sBetah9, sBetah10;

  cBetah5 = 0.0; cBetah6 = 0.0; cBetah7 = 0.0; cBetah8 = 0.0; cBetah9 = 0.0; cBetah10 = 0.0;
  sBetah5 = 0.0; sBetah6 = 0.0; sBetah7 = 0.0; sBetah8 = 0.0; sBetah9 = 0.0; sBetah10 = 0.0;


  /* L=2 */

  REAL8 d22[5]   = {sBetah4, 2.0*cBetah*sBetah3, SQRT->sqrt6*sBetah2*cBetah2, 2.0*cBetah3*sBetah, cBetah4};
  REAL8 d2m2[5]  = {d22[4],    -d22[3],      d22[2],     -d22[1],     d22[0]};

  REAL8 d21[5]   = {2.0*cBetah*sBetah3, 3.0*cBetah2*sBetah2 - sBetah4, SQRT->sqrt6*(cBetah3*sBetah - cBetah*sBetah3), cBetah2*(cBetah2 - 3.0*sBetah2), -2.0*cBetah3*sBetah};
  REAL8 d2m1[5]  = {-d21[4],   d21[3],     -d21[2],    d21[1],     -d21[0]};

  for(UINT4 i=0; i<5; i++){
      wS->wignerdL2[0][i] = d2m2[i];
      wS->wignerdL2[1][i] = d2m1[i];
      wS->wignerdL2[2][i] = 0.0;
      wS->wignerdL2[3][i] = d21[i];
      wS->wignerdL2[4][i] = d22[i];
    }

  switch(LMAX)
  {

    case 3:
    {
      cBetah5 = cBetah * cBetah4;
      cBetah6 = cBetah * cBetah5;

      sBetah5 = sBetah * sBetah4;
      sBetah6 = sBetah * sBetah5;

      REAL8 d33[7]   = {sBetah6, SQRT->sqrt6*cBetah*sBetah5, SQRT->sqrt15*cBetah2*sBetah4, 2.0*SQRT->sqrt5*cBetah3*sBetah3, SQRT->sqrt15*cBetah4*sBetah2, SQRT->sqrt6*cBetah5*sBetah, cBetah6};
      REAL8 d3m3[7]  = {d33[6],    -d33[5],     d33[4],      -d33[3],    d33[2],     -d33[1],    d33[0]};

      for(UINT4 i=0; i<7; i++){
        wS->wignerdL3[0][i] = d3m3[i];
        wS->wignerdL3[1][i] = 0.0;
        wS->wignerdL3[2][i] = 0.0;
        wS->wignerdL3[3][i] = 0.0;
        wS->wignerdL3[4][i] = 0.0;
        wS->wignerdL3[5][i] = 0.0;
        wS->wignerdL3[6][i] = d33[i];
        }

      break;
    }
    

    case 4:
    {
      cBetah5 = cBetah * cBetah4;
      cBetah6 = cBetah * cBetah5;
      cBetah7 = cBetah * cBetah6;
      cBetah8 = cBetah * cBetah7;

      sBetah5 = sBetah * sBetah4;
      sBetah6 = sBetah * sBetah5;
      sBetah7 = sBetah * sBetah6;
      sBetah8 = sBetah * sBetah7;

      REAL8 d33[7]   = {sBetah6, SQRT->sqrt6*cBetah*sBetah5, SQRT->sqrt15*cBetah2*sBetah4, 2.0*SQRT->sqrt5*cBetah3*sBetah3, SQRT->sqrt15*cBetah4*sBetah2, SQRT->sqrt6*cBetah5*sBetah, cBetah6};
      REAL8 d3m3[7]  = {d33[6],    -d33[5],     d33[4],      -d33[3],    d33[2],     -d33[1],    d33[0]};

      REAL8 d44[9]   = {sBetah8, 2.0*SQRT->sqrt2*cBetah*sBetah7, 2.0*SQRT->sqrt7*cBetah2*sBetah6, 2.0*SQRT->sqrt14*cBetah3*sBetah5, SQRT->sqrt70*cBetah4*sBetah4, 2.0*SQRT->sqrt14*cBetah5*sBetah3, 2.0*SQRT->sqrt7*cBetah6*sBetah2, 2.0*sqrt(2)*cBetah7*sBetah, cBetah8};
      REAL8 d4m4[9]  = {d44[8],     -d44[7],     d44[6],      -d44[5],     d44[4],    -d44[3],     d44[2],     -d44[1],    d44[0]};

      for(UINT4 i=0; i<7; i++){
        wS->wignerdL3[0][i] = d3m3[i];
        wS->wignerdL3[1][i] = 0.0;
        wS->wignerdL3[2][i] = 0.0;
        wS->wignerdL3[3][i] = 0.0;
        wS->wignerdL3[4][i] = 0.0;
        wS->wignerdL3[5][i] = 0.0;
        wS->wignerdL3[6][i] = d33[i];
        }

      for(UINT4 i=0; i<9; i++){
        wS->wignerdL4[0][i] = d4m4[i];
        wS->wignerdL4[1][i] = 0.0;
        wS->wignerdL4[2][i] = 0.0;
        wS->wignerdL4[3][i] = 0.0;
        wS->wignerdL4[4][i] = 0.0;
        wS->wignerdL4[5][i] = 0.0;
        wS->wignerdL4[6][i] = 0.0;
        wS->wignerdL4[7][i] = 0.0;
        wS->wignerdL4[8][i] = d44[i];
        } 

      break;   
    }

    case 5:
    {
      cBetah5 = cBetah * cBetah4;
      cBetah6 = cBetah * cBetah5;
      cBetah7 = cBetah * cBetah6;
      cBetah8 = cBetah * cBetah7;
      cBetah9 = cBetah * cBetah8;
      cBetah10 = cBetah * cBetah9;

      sBetah5 = sBetah * sBetah4;
      sBetah6 = sBetah * sBetah5;
      sBetah7 = sBetah * sBetah6;
      sBetah8 = sBetah * sBetah7;
      sBetah9 = sBetah * sBetah8;
      sBetah10 = sBetah * sBetah9;

      REAL8 d33[7]   = {sBetah6, SQRT->sqrt6*cBetah*sBetah5, SQRT->sqrt15*cBetah2*sBetah4, 2.0*SQRT->sqrt5*cBetah3*sBetah3, SQRT->sqrt15*cBetah4*sBetah2, SQRT->sqrt6*cBetah5*sBetah, cBetah6};
      REAL8 d3m3[7]  = {d33[6],    -d33[5],     d33[4],      -d33[3],    d33[2],     -d33[1],    d33[0]};

      REAL8 d44[9]   = {sBetah8, 2.0*SQRT->sqrt2*cBetah*sBetah7, 2.0*SQRT->sqrt7*cBetah2*sBetah6, 2.0*SQRT->sqrt14*cBetah3*sBetah5, SQRT->sqrt70*cBetah4*sBetah4, 2.0*SQRT->sqrt14*cBetah5*sBetah3, 2.0*SQRT->sqrt7*cBetah6*sBetah2, 2.0*sqrt(2)*cBetah7*sBetah, cBetah8};
      REAL8 d4m4[9]  = {d44[8],     -d44[7],     d44[6],      -d44[5],     d44[4],    -d44[3],     d44[2],     -d44[1],    d44[0]};

      REAL8 d55[11] = {sBetah10, SQRT->sqrt10*cBetah*sBetah9, 3*SQRT->sqrt5*cBetah2*sBetah8, 2*SQRT->sqrt30*cBetah3*sBetah7, SQRT->sqrt210*cBetah4*sBetah6, 6.0*SQRT->sqrt7*cBetah5*sBetah5,\
            SQRT->sqrt210*cBetah6*sBetah4, 2*SQRT->sqrt30*cBetah7*sBetah3, 3*SQRT->sqrt5*cBetah8*sBetah2, SQRT->sqrt10*cBetah9*sBetah, cBetah10};
      REAL8 d5m5[11] = {d55[10], -d55[9], d55[8], -d55[7], d55[6], -d55[5], d55[4], -d55[3], d55[2], -d55[1], d55[0]};

      for(UINT4 i=0; i<7; i++){
        wS->wignerdL3[0][i] = d3m3[i];
        wS->wignerdL3[1][i] = 0.0;
        wS->wignerdL3[2][i] = 0.0;
        wS->wignerdL3[3][i] = 0.0;
        wS->wignerdL3[4][i] = 0.0;
        wS->wignerdL3[5][i] = 0.0;
        wS->wignerdL3[6][i] = d33[i];
        }

      for(UINT4 i=0; i<9; i++){
        wS->wignerdL4[0][i] = d4m4[i];
        wS->wignerdL4[1][i] = 0.0;
        wS->wignerdL4[2][i] = 0.0;
        wS->wignerdL4[3][i] = 0.0;
        wS->wignerdL4[4][i] = 0.0;
        wS->wignerdL4[5][i] = 0.0;
        wS->wignerdL4[6][i] = 0.0;
        wS->wignerdL4[7][i] = 0.0;
        wS->wignerdL4[8][i] = d44[i];
        } 

      for(UINT4 i=0; i<11; i++){
        wS->wignerdL5[0][i] = d5m5[i];
        wS->wignerdL5[1][i] = 0.0;
        wS->wignerdL5[2][i] = 0.0;
        wS->wignerdL5[3][i] = 0.0;
        wS->wignerdL5[4][i] = 0.0;
        wS->wignerdL5[5][i] = 0.0;
        wS->wignerdL5[6][i] = 0.0;
        wS->wignerdL5[7][i] = 0.0;
        wS->wignerdL5[8][i] = 0.0;
        wS->wignerdL5[9][i] = 0.0;
        wS->wignerdL5[10][i] = d55[i];
        }

      break;
    }
  }


    /* For performing the global rotation between the J-frame and L0-frame, all coefficients for a given L have to be computed,
    since the list of modes in the J-frame contains all the modes for a given L */
    if(globalRot==1)
    {
    UNUSED REAL8 sinBeta = sign*sqrt(fabs(1.0 - cosBeta*cosBeta));

		UNUSED REAL8 cos2Beta = cosBeta*cosBeta - sinBeta*sinBeta;
		UNUSED REAL8 cos3Beta = cosBeta*(2.0*cos2Beta - 1.0);
		UNUSED REAL8 cos4Beta = pow(sinBeta,4) + pow(cosBeta,4) - 6.0*sinBeta*sinBeta*cosBeta*cosBeta;

    switch(LMAX)
    {
      case 2:
      {
        REAL8 d20[5]   = {SQRT->sqrt6*cBetah2*sBetah2 , SQRT->sqrt6*cBetah*sBetah*C2mS2 , 0.25*(1 + 3*(-4*cBetah2*sBetah2 + pow(C2mS2,2))) , -SQRT->sqrt6*cBetah*sBetah*C2mS2 , SQRT->sqrt6*cBetah2*sBetah2};

        for(UINT4 i=0; i<5; i++){
          wS->wignerdL2[2][i] = d20[i];
        }
        break;
      }

      case 3:
      {
        REAL8 d20[5]   = {SQRT->sqrt6*cBetah2*sBetah2 , SQRT->sqrt6*cBetah*sBetah*C2mS2 , 0.25*(1 + 3*(-4*cBetah2*sBetah2 + pow(C2mS2,2))) , -SQRT->sqrt6*cBetah*sBetah*C2mS2 , SQRT->sqrt6*cBetah2*sBetah2};

        REAL8 d32[7]   = {SQRT->sqrt6*cBetah*sBetah5, sBetah4*(5.0*cBetah2 - sBetah2), SQRT->sqrt10*sBetah3*(2.0*cBetah3 - cBetah*sBetah2), SQRT->sqrt30*cBetah2*(cBetah2 - sBetah2)*sBetah2, SQRT->sqrt10*cBetah3*(cBetah2*sBetah - 2.0*sBetah3), cBetah4*(cBetah2 - 5.0*sBetah2), -SQRT->sqrt6*cBetah5*sBetah};
        REAL8 d3m2[7]   = {-d32[6],d32[5],-d32[4],d32[3],-d32[2],d32[1],-d32[0]};

        REAL8 d31[7]   = {SQRT->sqrt15*cBetah2*sBetah4, SQRT->sqrt2half*cBetah*sBetah3*(1 + 3.0*C2mS2), 0.125*sBetah2*(13.0 + 20.0*C2mS2 + 15.0*(C2mS2*C2mS2 - 4.0*cBetah2*sBetah2)), 0.25*SQRT->sqrt3*cBetah*sBetah*(3.0 + 5.0*(C2mS2*C2mS2 - 4.0*cBetah2*sBetah2)),\
                          0.125*cBetah2*(13.0 - 20.0*C2mS2 + 15.0*(C2mS2*C2mS2 - 4.0*cBetah2*sBetah2)), -SQRT->sqrt2half*cBetah3*sBetah*(-1.0 + 3.0*C2mS2), SQRT->sqrt15*cBetah4*sBetah2};
        REAL8 d3m1[7]   = {d31[6], -d31[5], d31[4], -d31[3], d31[2], -d31[1], d31[0]};

        REAL8 d30[7]   = {2.0*SQRT->sqrt5*cBetah3*sBetah3, SQRT->sqrt30*cBetah2*sBetah2*C2mS2, 0.25*SQRT->sqrt3*cBetah*sBetah*(3.0 + 5.0*(C2mS2*C2mS2 - 4.0*cBetah2*sBetah2)), 0.125*(5.0*cos3Beta + 3*C2mS2),\
                          -0.25*SQRT->sqrt3*cBetah*sBetah*(3.0 + 5.0*(C2mS2*C2mS2 - 4.0*cBetah2*sBetah2)), SQRT->sqrt30*cBetah2*sBetah2*C2mS2, -2.0*SQRT->sqrt5*cBetah3*sBetah3};

        for(UINT4 i=0; i<5; i++){
          wS->wignerdL2[2][i] = d20[i];
        }
        for(UINT4 i=0; i<7; i++){
          wS->wignerdL3[1][i] = d3m2[i];
          wS->wignerdL3[2][i] = d3m1[i];
          wS->wignerdL3[3][i] = d30[i];
          wS->wignerdL3[4][i] = d31[i];
          wS->wignerdL3[5][i] = d32[i];
        }
        break;
      }

      case 4:
      {
        REAL8 d20[5]   = {SQRT->sqrt6*cBetah2*sBetah2 , SQRT->sqrt6*cBetah*sBetah*C2mS2 , 0.25*(1 + 3*(-4*cBetah2*sBetah2 + pow(C2mS2,2))) , -SQRT->sqrt6*cBetah*sBetah*C2mS2 , SQRT->sqrt6*cBetah2*sBetah2};

        REAL8 d32[7]   = {SQRT->sqrt6*cBetah*sBetah5, sBetah4*(5.0*cBetah2 - sBetah2), SQRT->sqrt10*sBetah3*(2.0*cBetah3 - cBetah*sBetah2), SQRT->sqrt30*cBetah2*(cBetah2 - sBetah2)*sBetah2, SQRT->sqrt10*cBetah3*(cBetah2*sBetah - 2.0*sBetah3), cBetah4*(cBetah2 - 5.0*sBetah2), -SQRT->sqrt6*cBetah5*sBetah};
        REAL8 d3m2[7]   = {-d32[6],d32[5],-d32[4],d32[3],-d32[2],d32[1],-d32[0]};

        REAL8 d31[7]   = {SQRT->sqrt15*cBetah2*sBetah4, SQRT->sqrt2half*cBetah*sBetah3*(1 + 3.0*C2mS2), 0.125*sBetah2*(13.0 + 20.0*C2mS2 + 15.0*(C2mS2*C2mS2 - 4.0*cBetah2*sBetah2)), 0.25*SQRT->sqrt3*cBetah*sBetah*(3.0 + 5.0*(C2mS2*C2mS2 - 4.0*cBetah2*sBetah2)),\
                          0.125*cBetah2*(13.0 - 20.0*C2mS2 + 15.0*(C2mS2*C2mS2 - 4.0*cBetah2*sBetah2)), -SQRT->sqrt2half*cBetah3*sBetah*(-1.0 + 3.0*C2mS2), SQRT->sqrt15*cBetah4*sBetah2};
        REAL8 d3m1[7]   = {d31[6], -d31[5], d31[4], -d31[3], d31[2], -d31[1], d31[0]};

        REAL8 d30[7]   = {2.0*SQRT->sqrt5*cBetah3*sBetah3, SQRT->sqrt30*cBetah2*sBetah2*C2mS2, 0.25*SQRT->sqrt3*cBetah*sBetah*(3.0 + 5.0*(C2mS2*C2mS2 - 4.0*cBetah2*sBetah2)), 0.125*(5.0*cos3Beta + 3*C2mS2),\
                          -0.25*SQRT->sqrt3*cBetah*sBetah*(3.0 + 5.0*(C2mS2*C2mS2 - 4.0*cBetah2*sBetah2)), SQRT->sqrt30*cBetah2*sBetah2*C2mS2, -2.0*SQRT->sqrt5*cBetah3*sBetah3};

        REAL8 d43[9]   = {2*SQRT->sqrt2*cBetah*sBetah7, 7*cBetah2*sBetah6-sBetah8, SQRT->sqrt14*(3*cBetah3*sBetah5-cBetah*sBetah7), SQRT->sqrt7*(5*cBetah4*sBetah4-3*cBetah2*sBetah6),\
                          2*5.916079783099616*(cBetah5*sBetah3-cBetah3*sBetah5), SQRT->sqrt7*(3*cBetah6*sBetah2-5*cBetah4*sBetah4), SQRT->sqrt14*(cBetah7*sBetah-3*cBetah5*sBetah3), cBetah8-7*cBetah6*sBetah2, -2.*SQRT->sqrt2*cBetah7*sBetah};
        REAL8 d4m3[9]   = {-d43[8], d43[7], -d43[6], d43[5],-d43[4],d43[3], -d43[2], d43[1], -d43[0]};

        REAL8 d42[9]   = {2*SQRT->sqrt7*cBetah2*sBetah6, SQRT->sqrt14*cBetah*sBetah5*(1.0 +2.0*C2mS2), sBetah4*(1.0 + 7.0*C2mS2 + 7.0*C2mS2*C2mS2), 0.5*SQRT->sqrt2*cBetah*sBetah3*(6.0 + 7.0*cos2Beta + 7.0*C2mS2),\
                          0.5*SQRT->sqrt2half*cBetah2*(5.0 + 7.0*cos2Beta)*sBetah2, 0.5*SQRT->sqrt2*cBetah3*sBetah*(6.0 + 7.0*cos2Beta - 7.0*C2mS2), cBetah4*(1.0 - 7.0*C2mS2 + 7.0*C2mS2*C2mS2), -SQRT->sqrt14*cBetah5*sBetah*(-1.0 +2.0*C2mS2), 2*SQRT->sqrt7*cBetah6*sBetah2};
        REAL8 d4m2[9]   = {d42[8], -d42[7], d42[6], -d42[5], d42[4], -d42[3], d42[2], -d42[1], d42[0]};

        REAL8 d41[9]   = {2*SQRT->sqrt14*cBetah3*sBetah5, SQRT->sqrt7*cBetah2*sBetah4*(1.0 + 4.0*C2mS2), 0.5*SQRT->sqrt2*cBetah*sBetah3*(6.0 + 7.0*cos2Beta + 7.0*C2mS2), 0.125*sBetah2*(15.0 +21.0*cos2Beta + 14.0*cos3Beta + 30.0*C2mS2),\
                          0.125*SQRT->sqrt5*cBetah*sBetah*(7.0*cos3Beta + 9.0*C2mS2), 0.125*cBetah2*(-15.0 + 30*cosBeta - 21.0*cos2Beta + 14.0*cos3Beta), 0.5*SQRT->sqrt2*cBetah3*sBetah*(-6.0 - 7.0*cos2Beta + 7.0*C2mS2),\
                          SQRT->sqrt7*cBetah4*sBetah2*(-1.0 + 4.0*C2mS2), -2*SQRT->sqrt14*cBetah5*sBetah3};
        REAL8 d4m1[9]   = {-d41[8], d41[7], -d41[6], d41[5], -d41[4], d41[3], -d41[2], d41[1], -d41[0]};

        REAL8 d40[9]   = {SQRT->sqrt70*cBetah4*sBetah4, 2*SQRT->sqrt35*cBetah3*sBetah3*C2mS2, 0.5*SQRT->sqrt2half*cBetah2*(5. + 7.*cos2Beta)*sBetah2, 0.125*SQRT->sqrt5*cBetah*sBetah*(7.*cos3Beta + 9.*C2mS2),\
                           0.015625*(9 + 20.*cos2Beta + 35.*cos4Beta), -0.125*SQRT->sqrt5*cBetah*sBetah*(7.*cos3Beta + 9.*C2mS2), 0.5*SQRT->sqrt2half*cBetah2*(5. + 7.*cos2Beta)*sBetah2, -2.*SQRT->sqrt35*cBetah3*sBetah3*C2mS2, SQRT->sqrt70*cBetah4*sBetah4};

        for(UINT4 i=0; i<5; i++){
          wS->wignerdL2[2][i] = d20[i];
        }
        for(UINT4 i=0; i<7; i++){
          wS->wignerdL3[1][i] = d3m2[i];
          wS->wignerdL3[2][i] = d3m1[i];
          wS->wignerdL3[3][i] = d30[i];
          wS->wignerdL3[4][i] = d31[i];
          wS->wignerdL3[5][i] = d32[i];
        }
        for(UINT4 i=0; i<9; i++){
          wS->wignerdL4[1][i] = d4m3[i];
          wS->wignerdL4[2][i] = d4m2[i];
          wS->wignerdL4[3][i] = d4m1[i];
          wS->wignerdL4[4][i] = d40[i];
          wS->wignerdL4[5][i] = d41[i];
          wS->wignerdL4[6][i] = d42[i];
          wS->wignerdL4[7][i] = d43[i];
        }
        break;
      }

      case 5:
      {
        REAL8 d20[5]   = {SQRT->sqrt6*cBetah2*sBetah2 , SQRT->sqrt6*cBetah*sBetah*C2mS2 , 0.25*(1 + 3*(-4*cBetah2*sBetah2 + pow(C2mS2,2))) , -SQRT->sqrt6*cBetah*sBetah*C2mS2 , SQRT->sqrt6*cBetah2*sBetah2};

        REAL8 d32[7]   = {SQRT->sqrt6*cBetah*sBetah5, sBetah4*(5.0*cBetah2 - sBetah2), SQRT->sqrt10*sBetah3*(2.0*cBetah3 - cBetah*sBetah2), SQRT->sqrt30*cBetah2*(cBetah2 - sBetah2)*sBetah2, SQRT->sqrt10*cBetah3*(cBetah2*sBetah - 2.0*sBetah3), cBetah4*(cBetah2 - 5.0*sBetah2), -SQRT->sqrt6*cBetah5*sBetah};
        REAL8 d3m2[7]   = {-d32[6],d32[5],-d32[4],d32[3],-d32[2],d32[1],-d32[0]};

        REAL8 d31[7]   = {SQRT->sqrt15*cBetah2*sBetah4, SQRT->sqrt2half*cBetah*sBetah3*(1 + 3.0*C2mS2), 0.125*sBetah2*(13.0 + 20.0*C2mS2 + 15.0*(C2mS2*C2mS2 - 4.0*cBetah2*sBetah2)), 0.25*SQRT->sqrt3*cBetah*sBetah*(3.0 + 5.0*(C2mS2*C2mS2 - 4.0*cBetah2*sBetah2)),\
                          0.125*cBetah2*(13.0 - 20.0*C2mS2 + 15.0*(C2mS2*C2mS2 - 4.0*cBetah2*sBetah2)), -SQRT->sqrt2half*cBetah3*sBetah*(-1.0 + 3.0*C2mS2), SQRT->sqrt15*cBetah4*sBetah2};
        REAL8 d3m1[7]   = {d31[6], -d31[5], d31[4], -d31[3], d31[2], -d31[1], d31[0]};

        REAL8 d30[7]   = {2.0*SQRT->sqrt5*cBetah3*sBetah3, SQRT->sqrt30*cBetah2*sBetah2*C2mS2, 0.25*SQRT->sqrt3*cBetah*sBetah*(3.0 + 5.0*(C2mS2*C2mS2 - 4.0*cBetah2*sBetah2)), 0.125*(5.0*cos3Beta + 3*C2mS2),\
                          -0.25*SQRT->sqrt3*cBetah*sBetah*(3.0 + 5.0*(C2mS2*C2mS2 - 4.0*cBetah2*sBetah2)), SQRT->sqrt30*cBetah2*sBetah2*C2mS2, -2.0*SQRT->sqrt5*cBetah3*sBetah3};

        REAL8 d43[9]   = {2*SQRT->sqrt2*cBetah*sBetah7, 7*cBetah2*sBetah6-sBetah8, SQRT->sqrt14*(3*cBetah3*sBetah5-cBetah*sBetah7), SQRT->sqrt7*(5*cBetah4*sBetah4-3*cBetah2*sBetah6),\
                          2*5.916079783099616*(cBetah5*sBetah3-cBetah3*sBetah5), SQRT->sqrt7*(3*cBetah6*sBetah2-5*cBetah4*sBetah4), SQRT->sqrt14*(cBetah7*sBetah-3*cBetah5*sBetah3), cBetah8-7*cBetah6*sBetah2, -2.*SQRT->sqrt2*cBetah7*sBetah};
        REAL8 d4m3[9]   = {-d43[8], d43[7], -d43[6], d43[5],-d43[4],d43[3], -d43[2], d43[1], -d43[0]};

        REAL8 d42[9]   = {2*SQRT->sqrt7*cBetah2*sBetah6, SQRT->sqrt14*cBetah*sBetah5*(1.0 +2.0*C2mS2), sBetah4*(1.0 + 7.0*C2mS2 + 7.0*C2mS2*C2mS2), 0.5*SQRT->sqrt2*cBetah*sBetah3*(6.0 + 7.0*cos2Beta + 7.0*C2mS2),\
                          0.5*SQRT->sqrt2half*cBetah2*(5.0 + 7.0*cos2Beta)*sBetah2, 0.5*SQRT->sqrt2*cBetah3*sBetah*(6.0 + 7.0*cos2Beta - 7.0*C2mS2), cBetah4*(1.0 - 7.0*C2mS2 + 7.0*C2mS2*C2mS2), -SQRT->sqrt14*cBetah5*sBetah*(-1.0 +2.0*C2mS2), 2*SQRT->sqrt7*cBetah6*sBetah2};
        REAL8 d4m2[9]   = {d42[8], -d42[7], d42[6], -d42[5], d42[4], -d42[3], d42[2], -d42[1], d42[0]};

        REAL8 d41[9]   = {2*SQRT->sqrt14*cBetah3*sBetah5, SQRT->sqrt7*cBetah2*sBetah4*(1.0 + 4.0*C2mS2), 0.5*SQRT->sqrt2*cBetah*sBetah3*(6.0 + 7.0*cos2Beta + 7.0*C2mS2), 0.125*sBetah2*(15.0 +21.0*cos2Beta + 14.0*cos3Beta + 30.0*C2mS2),\
                          0.125*SQRT->sqrt5*cBetah*sBetah*(7.0*cos3Beta + 9.0*C2mS2), 0.125*cBetah2*(-15.0 + 30*cosBeta - 21.0*cos2Beta + 14.0*cos3Beta), 0.5*SQRT->sqrt2*cBetah3*sBetah*(-6.0 - 7.0*cos2Beta + 7.0*C2mS2),\
                          SQRT->sqrt7*cBetah4*sBetah2*(-1.0 + 4.0*C2mS2), -2*SQRT->sqrt14*cBetah5*sBetah3};
        REAL8 d4m1[9]   = {-d41[8], d41[7], -d41[6], d41[5], -d41[4], d41[3], -d41[2], d41[1], -d41[0]};

        REAL8 d40[9]   = {SQRT->sqrt70*cBetah4*sBetah4, 2*SQRT->sqrt35*cBetah3*sBetah3*C2mS2, 0.5*SQRT->sqrt2half*cBetah2*(5. + 7.*cos2Beta)*sBetah2, 0.125*SQRT->sqrt5*cBetah*sBetah*(7.*cos3Beta + 9.*C2mS2),\
                           0.015625*(9 + 20.*cos2Beta + 35.*cos4Beta), -0.125*SQRT->sqrt5*cBetah*sBetah*(7.*cos3Beta + 9.*C2mS2), 0.5*SQRT->sqrt2half*cBetah2*(5. + 7.*cos2Beta)*sBetah2, -2.*SQRT->sqrt35*cBetah3*sBetah3*C2mS2, SQRT->sqrt70*cBetah4*sBetah4};

        REAL8 d54[11] = {SQRT->sqrt10*cBetah*sBetah9, sBetah8*(4.0 + 5.0*C2mS2), (3./sqrt(2))*cBetah*sBetah7*(3.0 +5.0*C2mS2), 2.0*SQRT->sqrt3*cBetah2*sBetah6*(2.0 + 5.0*C2mS2), SQRT->sqrt21*cBetah3*sBetah5*(1.0 + 5.0*C2mS2),\
                          3.0*SQRT->sqrt70*cBetah4*sBetah4*C2mS2, SQRT->sqrt21*cBetah5*sBetah3*(-1.0 + 5.0*C2mS2), 2.0*SQRT->sqrt3*cBetah6*sBetah2*(-2.0 + 5.0*C2mS2), (3./sqrt(2))*cBetah7*sBetah*(-3.0 +5.0*C2mS2), cBetah8*(-4.0 + 5.0*C2mS2), -SQRT->sqrt10*cBetah9*sBetah};
        REAL8 d5m4[11] = {-d54[10], d54[9], -d54[8], d54[7], -d54[6], d54[5], -d54[4], d54[3], -d54[2], d54[1], -d54[0]};

        REAL8 d53[11] = {3.0*SQRT->sqrt5*cBetah2*sBetah8, (3.0/SQRT->sqrt2)*cBetah*sBetah7*(3.0+5.0*C2mS2), 0.25*(13.0 + 54.0*C2mS2 + 45.0*C2mS2*C2mS2)*sBetah6, sqrt(1.5)*(1.0+12.0*C2mS2+15.0*C2mS2*C2mS2)*cBetah*sBetah5, 0.5*sqrt(10.5)*(-1.0+6.0*C2mS2+15.0*C2mS2*C2mS2)*cBetah2*sBetah4,\
                        0.25*SQRT->sqrt35*(7.0 + 9.0*cos2Beta)*cBetah3*sBetah3, 0.5*sqrt(10.5)*(-1.0-6.0*C2mS2+15.0*C2mS2*C2mS2)*cBetah4*sBetah2, sqrt(1.5)*(1.0-12.0*C2mS2+15.0*C2mS2*C2mS2)*cBetah5*sBetah, 0.25*(13.0 - 54.0*C2mS2 + 45.0*C2mS2*C2mS2)*cBetah6, (3.0/SQRT->sqrt2)*cBetah7*sBetah*(3.0-5.0*C2mS2), 3.0*SQRT->sqrt5*cBetah8*sBetah2};
        REAL8 d5m3[11] = {d53[10], -d53[9], d53[8], -d53[7], d53[6], -d53[5], d53[4], -d53[3], d53[2], -d53[1], d53[0]};

        REAL8 d52[11] = {2*SQRT->sqrt30*cBetah3*sBetah7, 2.0*SQRT->sqrt3*(2.0+5.0*C2mS2)*cBetah2*sBetah6, sqrt(1.5)*(1.0+12.0*C2mS2+15.0*C2mS2*C2mS2)*cBetah*sBetah5, (-1.0+3.0*C2mS2+18.0*C2mS2*C2mS2+15.0*C2mS2*C2mS2*C2mS2)*sBetah4,\
                        0.5*SQRT->sqrt7*(-1.0-3.0*C2mS2+9.0*C2mS2*C2mS2+15.0*C2mS2*C2mS2*C2mS2)*cBetah*sBetah3, 0.5*sqrt(52.5)*C2mS2*cBetah2*sBetah2*(1.0+3.0*cos2Beta), 0.5*SQRT->sqrt7*(1.0-3.0*C2mS2-9.0*C2mS2*C2mS2+15.0*C2mS2*C2mS2*C2mS2)*cBetah3*sBetah,\
                        (1.0+3.0*C2mS2-18.0*C2mS2*C2mS2+15.0*C2mS2*C2mS2*C2mS2)*cBetah4, -sqrt(1.5)*(1.0-12.0*C2mS2+15.0*C2mS2*C2mS2)*cBetah5*sBetah, 2.0*SQRT->sqrt3*(-2.0+5.0*C2mS2)*cBetah6*sBetah2, -2*SQRT->sqrt30*cBetah7*sBetah3};
        REAL8 d5m2[11] = {-d52[10], d52[9], -d52[8], d52[7], -d52[6], d52[5], -d52[4], d52[3], -d52[2], d52[1], -d52[0]};

        REAL8 d51[11] = {SQRT->sqrt210*cBetah4*sBetah6, SQRT->sqrt21*(1.0+5.0*C2mS2)*cBetah3*sBetah5, 0.5*sqrt(10.5)*(-1.0+6.0*C2mS2+15.0*C2mS2*C2mS2)*cBetah2*sBetah4,\
             0.5*SQRT->sqrt7*(-1.0-3.0*C2mS2+9.0*C2mS2*C2mS2+15.0*C2mS2*C2mS2*C2mS2)*cBetah*sBetah3, 0.125*(1.0-28.0*C2mS2-42.0*C2mS2*C2mS2+84.0*C2mS2*C2mS2*C2mS2+105.0*C2mS2*C2mS2*C2mS2*C2mS2)*sBetah2,\
             sqrt(7.5)/32.0*cBetah*sBetah*(15.0 + 28.0*cos2Beta + 21.0*cos4Beta), 0.125*(1.0+28.0*C2mS2-42.0*C2mS2*C2mS2-84.0*C2mS2*C2mS2*C2mS2+105.0*C2mS2*C2mS2*C2mS2*C2mS2)*cBetah2,\
             -0.5*SQRT->sqrt7*(1.0-3.0*C2mS2-9.0*C2mS2*C2mS2+15.0*C2mS2*C2mS2*C2mS2)*cBetah3*sBetah, 0.5*sqrt(10.5)*(-1.0-6.0*C2mS2+15.0*C2mS2*C2mS2)*cBetah4*sBetah2,\
             -SQRT->sqrt21*(-1.0+5.0*C2mS2)*cBetah5*sBetah3, SQRT->sqrt210*cBetah6*sBetah4};
        REAL8 d5m1[11] = {d51[10], -d51[9], d51[8], -d51[7], d51[6], -d51[5], d51[4], -d51[3], d51[2], -d51[1], d51[0]};

        REAL8 d50[11] = {6.0*SQRT->sqrt7*cBetah5*sBetah5, 3.0*SQRT->sqrt70*C2mS2*cBetah4*sBetah4, 0.25*SQRT->sqrt35*cBetah3*sBetah3*(7.0+9.0*cos2Beta), 0.5*sqrt(52.2)*C2mS2*cBetah2*sBetah2*(1.0+3.0*cos2Beta),\
            sqrt(7.5)/32.0*cBetah*sBetah*(15.0+28.0*cos2Beta+21.0*cos4Beta), 0.125*C2mS2*(15.0-70.0*C2mS2*C2mS2+63.0*C2mS2*C2mS2*C2mS2*C2mS2), -sqrt(7.5)/32.0*cBetah*sBetah*(15.0+28.0*cos2Beta+21.0*cos4Beta),\
            0.5*sqrt(52.2)*C2mS2*cBetah2*sBetah2*(1.0+3.0*cos2Beta), -0.25*SQRT->sqrt35*cBetah3*sBetah3*(7.0+9.0*cos2Beta), 3.0*SQRT->sqrt70*C2mS2*cBetah4*sBetah4, -6.0*SQRT->sqrt7*cBetah5*sBetah5};

        for(UINT4 i=0; i<5; i++){
          wS->wignerdL2[2][i] = d20[i];
        }
        for(UINT4 i=0; i<7; i++){
          wS->wignerdL3[1][i] = d3m2[i];
          wS->wignerdL3[2][i] = d3m1[i];
          wS->wignerdL3[3][i] = d30[i];
          wS->wignerdL3[4][i] = d31[i];
          wS->wignerdL3[5][i] = d32[i];
        }
        for(UINT4 i=0; i<9; i++){
          wS->wignerdL4[1][i] = d4m3[i];
          wS->wignerdL4[2][i] = d4m2[i];
          wS->wignerdL4[3][i] = d4m1[i];
          wS->wignerdL4[4][i] = d40[i];
          wS->wignerdL4[5][i] = d41[i];
          wS->wignerdL4[6][i] = d42[i];
          wS->wignerdL4[7][i] = d43[i];
        }
        for(UINT4 i=0; i<11; i++){
          wS->wignerdL5[1][i] = d5m4[i];
          wS->wignerdL5[2][i] = d5m3[i];
          wS->wignerdL5[3][i] = d5m2[i];
          wS->wignerdL5[4][i] = d5m1[i];
          wS->wignerdL5[5][i] = d50[i];
          wS->wignerdL5[6][i] = d51[i];
          wS->wignerdL5[7][i] = d52[i];
          wS->wignerdL5[8][i] = d53[i];
          wS->wignerdL5[9][i] = d54[i];
        }
        break;
      }
    }
  }

	/* expAlpha powers */

	COMPLEX16 expAlpha = cexp(I*alpha);
	COMPLEX16 expAlpham = 1./expAlpha;

	wS->expAlphaPowers[4] = expAlpha;
	wS->expAlphaPowers[3] = expAlpha*expAlpha;
	wS->expAlphaPowers[2] = expAlpha*wS->expAlphaPowers[3];
	wS->expAlphaPowers[1] = expAlpha*wS->expAlphaPowers[2];
	wS->expAlphaPowers[0] = expAlpha*wS->expAlphaPowers[1];
	wS->expAlphaPowers[5] = 1.0;
	wS->expAlphaPowers[6] = expAlpham;
	wS->expAlphaPowers[7] = expAlpham*expAlpham;
	wS->expAlphaPowers[8] = expAlpham*wS->expAlphaPowers[7];
	wS->expAlphaPowers[9] = expAlpham*wS->expAlphaPowers[8];
	wS->expAlphaPowers[10] = expAlpham*wS->expAlphaPowers[9];

	/* expGamma powers */

	COMPLEX16 expGamma = cexp(I*gamma);
	COMPLEX16 expGammam = 1./expGamma;

	wS->expGammaPowers[4] = expGamma;
	wS->expGammaPowers[3] = expGamma*expGamma;
	wS->expGammaPowers[2] = expGamma*wS->expGammaPowers[3];
	wS->expGammaPowers[1] = expGamma*wS->expGammaPowers[2];
	wS->expGammaPowers[0] = expGamma*wS->expGammaPowers[1];
	wS->expGammaPowers[5] = 1.0;
	wS->expGammaPowers[6] = expGammam;
	wS->expGammaPowers[7] = expGammam*expGammam;
	wS->expGammaPowers[8] = expGammam*wS->expGammaPowers[7];
	wS->expGammaPowers[9] = expGammam*wS->expGammaPowers[8];
	wS->expGammaPowers[10] = expGammam*wS->expGammaPowers[9];
}

/* Function to compute the WignerDMatrix. Analogous to XLALWignerDMatrix (which employs XLALWignerdMatrix), but Wignerd coefficients
here are read from the precomputed struct, in order to improve computational speed. */
COMPLEX16 PhenomTWignerDMatrix(
                                   int l,        /**< mode number l */
                                   int mp,        /**< mode number m' */
                                   int m,        /**< mode number m */
                                   PhenomTPWignerStruct *wS  /**< Wigner Struct, which contains Wignerd elements and powers of alpha and gamma */
    )
{
	COMPLEX16 wignerd = 0.0;

	switch(l) // For a given l, m, m' it returns the corresponding Wigner-d matrix coefficient, precomputed in the Wigner Struct
	{
		case 2:
			wignerd = (COMPLEX16)wS->wignerdL2[2+mp][2+m];
			break;
		case 3:
			wignerd = (COMPLEX16)wS->wignerdL3[3+mp][3+m];
			break;
		case 4:
			wignerd = (COMPLEX16)wS->wignerdL4[4+mp][4+m];
			break;
		case 5:
			wignerd = (COMPLEX16)wS->wignerdL5[5+mp][5+m];
			break;

	}

	// Construct the Wigner-D matrix coefficient adding the exponentials of alpha and gamma for the corresponding l, m, m'
	COMPLEX16 WignerD = wS->expAlphaPowers[mp+5]*wignerd*wS->expGammaPowers[m+5];

	return WignerD;
}

/* Function to rotate a given set of Spin-Weighted Spherical Harmonic mode time series with an time-dependent Euler rotation specified by three Euler angles.
Adaptation of XLALSimInspiralPrecessionRotateModes in SimInspiralPrecess.c , only difference is that at each time step the WignerDStruct is set and PhenomTWignerDMatrix
is then employed. */
int PhenomTPHM_RotateModes(
                SphHarmTimeSeries* h_lm, /**< spherical harmonic decomposed modes, modified in place */
                REAL8TimeSeries* alpha, /**< alpha Euler angle time series */
                REAL8TimeSeries* cosbeta, /**< beta Euler angle time series */
                REAL8TimeSeries* gam, /**< gamma Euler angle time series */
                PhenomT_precomputed_sqrt *SQRT /**< precomputed squared root factors */
){

	unsigned int i;
	int l, lmax, m, mp;
	lmax = XLALSphHarmTimeSeriesGetMaxL( h_lm );

	// Temporary holding variables
	complex double *x_lm = XLALCalloc( 2*lmax+1, sizeof(complex double) );
	COMPLEX16TimeSeries **h_xx = XLALCalloc( 2*lmax+1, sizeof(COMPLEX16TimeSeries) );

	for(i=0; i<alpha->data->length; i++){
		PhenomTPWignerStruct *wStruct;
		wStruct    = XLALMalloc(sizeof(PhenomTPWignerStruct));
		IMRPhenomTPHM_SetWignerDStruct(wStruct, SQRT, cosbeta->data->data[i], alpha->data->data[i], gam->data->data[i], lmax, 1, 0);
		for(l=2; l<=lmax; l++){
			for(m=0; m<2*l+1; m++){
				h_xx[m] = XLALSphHarmTimeSeriesGetMode(h_lm, l, m-l);
				if( !h_xx[m] ){
					x_lm[m] = 0;
				} else {
					x_lm[m] = h_xx[m]->data->data[i];
					h_xx[m]->data->data[i] = 0;
				}
			}

			for(m=0; m<2*l+1; m++){
				for(mp=0; mp<2*l+1; mp++){
					if( !h_xx[m] ) continue;
					if(!(creal(h_xx[m]->data->data[i])==0 && creal(x_lm[mp])==0)) {
					  h_xx[m]->data->data[i] +=
					    x_lm[mp] * PhenomTWignerDMatrix( l, mp-l, m-l, wStruct );
					  }
				}
			}
		}
		LALFree(wStruct);
	}


	XLALFree( x_lm );
	XLALFree( h_xx );
	return XLAL_SUCCESS;
}

/* Function to rotate a given set of Spin-Weighted Spherical Harmonic mode time series with an time-dependent Euler rotation specified by three Euler angles.
Main difference with previous function is that for a fixed rotation, Wigner-D matrix only have to be computed once. */
int PhenomTPHM_RotateModes_Global(
                SphHarmTimeSeries* h_lm, /**< spherical harmonic decomposed modes, modified in place */
                REAL8 alpha, /**< alpha Euler angle time series */
                REAL8 cosbeta, /**< beta Euler angle time series */
                REAL8 gam, /**< gamma Euler angle time series */
				        size_t length, /**< Length of waveform modes */
                PhenomT_precomputed_sqrt *SQRT /**< precomputed squared root factors */
){

	unsigned int i;
	int l, lmax, m, mp;
	lmax = XLALSphHarmTimeSeriesGetMaxL( h_lm );
	// Temporary holding variables
	complex double *x_lm = XLALCalloc( 2*lmax+1, sizeof(complex double) );
	COMPLEX16TimeSeries **h_xx = XLALCalloc( 2*lmax+1, sizeof(COMPLEX16TimeSeries) );

	PhenomTPWignerStruct *wStruct;
	wStruct    = XLALMalloc(sizeof(PhenomTPWignerStruct));
	IMRPhenomTPHM_SetWignerDStruct(wStruct, SQRT, cosbeta, alpha, gam, lmax, -1, 1);

	COMPLEX16 wigner[lmax-1][2*lmax+1][2*lmax+1];

	for(l=2; l<=lmax; l++){
		for(m=0; m<2*l+1; m++){
			for(mp=0; mp<2*l+1; mp++){
				wigner[l-2][mp][m] = PhenomTWignerDMatrix( l, mp-l, m-l, wStruct );
			}
		}
	}

	for(i=0; i<length; i++){
		for(l=2; l<=lmax; l++){
			for(m=0; m<2*l+1; m++){
				h_xx[m] = XLALSphHarmTimeSeriesGetMode(h_lm, l, m-l);
				if( !h_xx[m] ){
					x_lm[m] = 0;
				} else {
					x_lm[m] = h_xx[m]->data->data[i];
					h_xx[m]->data->data[i] = 0;
				}
			}

			for(m=0; m<2*l+1; m++){
				for(mp=0; mp<2*l+1; mp++){
					if( !h_xx[m] ) continue;
					if(!(creal(h_xx[m]->data->data[i])==0 && creal(x_lm[mp])==0)) {
					  h_xx[m]->data->data[i] += x_lm[mp] * wigner[l-2][mp][m];
					  }
				}
			}
		}
	}
	XLALFree(wStruct);
	XLALFree( x_lm );
	XLALFree( h_xx );
	return XLAL_SUCCESS;
}