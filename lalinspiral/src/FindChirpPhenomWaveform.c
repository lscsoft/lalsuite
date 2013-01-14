/*
*  Copyright (C) 2007 Santamaria L, Krishnan B, Dias M, Parameswaran A
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

/*-----------------------------------------------------------------------
 *
 * File Name: FindChirpSPData.c
 *
 * Author: Santamaria L, Krishnan B, Dias M, Parameswaran A
 *
 *
 *-----------------------------------------------------------------------
 */


#include <lal/LALDatatypes.h>
#include <lal/ComplexFFT.h>
#include <lal/LALInspiral.h>
#include <lal/LALInspiralBank.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/FindChirpDatatypes.h>
#include <lal/FindChirpChisq.h>
#include <lal/LALConfig.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/LALDatatypes.h>
#include <lal/AVFactories.h>
#include <lal/NRWaveIO.h>
#include <lal/NRWaveInject.h>
#include <lal/Inject.h>
#include <lal/FileIO.h>
#include <lal/Units.h>
#include <lal/FrequencySeries.h>
#include <lal/TimeSeries.h>
#include <lal/TimeFreqFFT.h>
#include <lal/VectorOps.h>
#include <lal/LALDetectors.h>
#include <lal/ComplexFFT.h>

#ifdef  __cplusplus
extern "C" {
#pragma }
#endif

typedef struct tagPhenomCoeffs{
    REAL8 fMerg_a;
    REAL8 fMerg_b;
    REAL8 fMerg_c;
    REAL8 fRing_a;
    REAL8 fRing_b;
    REAL8 fRing_c;
    REAL8 sigma_a;
    REAL8 sigma_b;
    REAL8 sigma_c;
    REAL8 fCut_a;
    REAL8 fCut_b;
    REAL8 fCut_c;
    REAL8 psi0_x;
    REAL8 psi0_y;
    REAL8 psi0_z;
    REAL8 psi2_x;
    REAL8 psi2_y;
    REAL8 psi2_z;
    REAL8 psi3_x;
    REAL8 psi3_y;
    REAL8 psi3_z;
    REAL8 psi4_x;
    REAL8 psi4_y;
    REAL8 psi4_z;
    REAL8 psi6_x;
    REAL8 psi6_y;
    REAL8 psi6_z;
    REAL8 psi7_x;
    REAL8 psi7_y;
    REAL8 psi7_z;
}
PhenomCoeffs;

typedef struct tagPhenomParams{
  REAL8 fMerger;
  REAL8 fRing;
  REAL8 fCut;
  REAL8 sigma;
  REAL8 psi0;
  REAL8 psi2;
  REAL8 psi3;
  REAL8 psi4;
  REAL8 psi6;
  REAL8 psi7;
}
PhenomParams;

void GetPhenomCoeffsLongJena(
			     PhenomCoeffs *co);

void ComputeParamsFromCoeffs(
			     PhenomParams *params,
			     PhenomCoeffs *coeffs,
			     REAL8        eta,
			     REAL8        M);

REAL4FrequencySeries *
XLALHybridP1Amplitude(
		      PhenomParams *params,
		      REAL8        fLow,
       		      REAL8        df,
		      REAL8        eta,
		      REAL8        M,
		      UINT4        len  );

REAL4FrequencySeries *
XLALHybridP1Phase(
		  PhenomParams  *params,
		  REAL8         fLow,
		  REAL8         df,
		  REAL8         eta,
		  REAL8         M,
		  UINT4         len );

REAL8
XLALLorentzian (
		REAL8 freq,
		REAL8 fRing,
		REAL8 sigma  );





/** This function contains the coeffs from the matching with the LONG
 * Jena waveforms (those are not the ones published in the paper)
 */
void
GetPhenomCoeffsLongJena(
			PhenomCoeffs *co)
{
  co->fMerg_a = 6.6389e-01; co->fMerg_b = -1.0321e-01; co->fMerg_c = 1.0979e-01;
  co->fRing_a = 1.3278e+00; co->fRing_b = -2.0642e-01; co->fRing_c = 2.1957e-01;
  co->sigma_a = 1.1383e+00; co->sigma_b = -1.7700e-01; co->sigma_c = 4.6834e-02;
  co->fCut_a = 1.7086e+00; co->fCut_b = -2.6592e-01; co->fCut_c = 2.8236e-01;

  co->psi0_x = -1.5829e-01; co->psi0_y = 8.7016e-02; co->psi0_z = -3.3382e-02;
  co->psi2_x = 3.2967e+01; co->psi2_y = -1.9000e+01; co->psi2_z = 2.1345e+00;
  co->psi3_x = -3.0849e+02; co->psi3_y = 1.8211e+02; co->psi3_z = -2.1727e+01;
  co->psi4_x = 1.1525e+03; co->psi4_y = -7.1477e+02; co->psi4_z = 9.9692e+01;
  co->psi6_x = 1.2057e+03; co->psi6_y = -8.4233e+02; co->psi6_z = 1.8046e+02;
  co->psi7_x = -0.0000e+00; co->psi7_y = 0.0000e+00; co->psi7_z = 0.0000e+00;
}


REAL8
XLALLorentzian (
		REAL8 freq,
		REAL8 fRing,
		REAL8 sigma  )
{
  REAL8 out;

  out = sigma / (2 * LAL_PI * ((freq - fRing)*(freq - fRing)
			       + sigma*sigma / 4.0));

  return(out);
}


/** Computes effective phase as in arXiv:0710.2335 [gr-qc] */
REAL4FrequencySeries *
XLALHybridP1Phase(
		  PhenomParams  *params,
		  REAL8         fLow,
		  REAL8         df,
		  REAL8         eta,
		  REAL8         M,
		  UINT4         n )
{
  UINT4 k;
  REAL8 piM;
  REAL8 f, psi0, psi2, psi3, psi4, psi6, psi7;
  REAL8 softfLow, softfCut;

  REAL4FrequencySeries *Phieff = NULL;

  LIGOTimeGPS epoch;

  psi0 = params -> psi0;
  psi2 = params -> psi2;
  psi3 = params -> psi3;
  psi4 = params -> psi4;
  psi6 = params -> psi6;
  psi7 = params -> psi7;

  piM = LAL_PI * M * LAL_MTSUN_SI;

  /* Allocate memory for the frequency series */
  epoch.gpsSeconds = 0;
  epoch.gpsNanoSeconds = 0;
  Phieff =
    XLALCreateREAL4FrequencySeries("", &epoch, 0, df, &lalDimensionlessUnit, n);

  /* We will soften the discontinuities of this function by multiplying it by
     the function (1/4)[1+tanh(k(f-fLow))][1-tanh(k(f-fCut))] with k=1.0.
     The step function is now a soft step. We estimate its width by requiring
     that the step function at the new boundaries takes the value 0.001
     Solving the eq. with Mathematica leads to a width = 3.45338, we take 3.5 */

  softfLow = fLow - 3.5;
  softfCut = params->fCut + 3.5;

  f = 0.0;

  for( k = 0 ; k < n ; k++ ) {

    if ( f <= softfLow || f > softfCut ) {
      Phieff->data->data[k] = 0.0;
    }
    else {

    /* QUESTION: what happens with the 2*pi*f*t_0 and psi_0 terms? */
    /* for the moment they're set to zero abd this seems to work   */
    /* (see Eq. (4.19) of the paper */

    Phieff->data->data[k] = psi0 * pow(f*piM , -5./3.) +
                            psi2 * pow(f*piM , -3./3.) +
                            psi3 * pow(f*piM , -2./3.) +
                            psi4 * pow(f*piM , -1./3.) +
                            psi6 * pow(f*piM , 1./3.) +
                            psi7 * pow(f*piM , 2./3.);

    Phieff->data->data[k] /= eta;
    }
    f += df;
  }

  return Phieff;
}


REAL4FrequencySeries *
XLALHybridP1Amplitude(
		      PhenomParams *params,
		      REAL8        fLow,
       		      REAL8        df,
		      REAL8        eta,
		      REAL8        M,
		      UINT4        n  )
{
  UINT4 k;
  REAL8 cConst;
  REAL8 f, fNorm, fMerg, fRing, fCut, sigma;
  REAL8 softfLow, softfCut, softFact;

  REAL4FrequencySeries *Aeff = NULL;

  LIGOTimeGPS epoch;

  INT4 sharpNess;

  fMerg = params->fMerger;
  fCut = params->fCut;
  fRing = params->fRing;
  sigma = params->sigma;

  /* Set amplitude of the wave (Ajith et al. Eq. 4.17) */
  cConst = pow(LAL_MTSUN_SI*M, 5./6.)*pow(fMerg,-7./6.)/pow(LAL_PI,2./3.);
  cConst *= pow(5.*eta/24., 1./2.);

  /* Allocate memory for the frequency series */
  epoch.gpsSeconds = 0;
  epoch.gpsNanoSeconds = 0;
  Aeff =
    XLALCreateREAL4FrequencySeries("", &epoch, 0, df, &lalDimensionlessUnit, n);

  f = 0.0;

  /* We will soften the discontinuities of this function by multiplying it by
     the function (1/4)[1+tanh(k(f-fLow))][1-tanh(k(f-fCut))] with k=1.0.
     The step function is now a soft step. We estimate its width by requiring
     that the step function at the new boundaries takes the value 0.001
     Solving the eq. with Mathematica leads to a width = 3.45338, we take 3.5 */

  softfLow = fLow - 3.5;
  softfCut = fCut + 3.5;

  sharpNess = 1;
  for( k = 0 ; k < n ; k++ ) {

    fNorm = f / fMerg;
    softFact = (1+tanh(sharpNess*(f-fLow)))*(1-tanh(sharpNess*(f-fCut)))/4.;

    if ( f <= softfLow || f > softfCut ) {
       Aeff->data->data[k] = 0.0;
    }
    else if ( f > softfLow && f <= fMerg ) {
      Aeff->data->data[k] = pow (fNorm, -7./6.);
      Aeff->data->data[k] *= softFact;
    }
    else if ( f > fMerg && f <= fRing ) {
      Aeff->data->data[k] = pow (fNorm, -2./3.);
      Aeff->data->data[k] *= softFact;
    }
    else if ( f > fRing && f <= softfCut ) {
      Aeff->data->data[k] = XLALLorentzian ( f, fRing, sigma);
      Aeff->data->data[k] *= LAL_PI_2*pow(fRing/fMerg,-2./3.)*sigma;
      Aeff->data->data[k] *= softFact;
    }
    Aeff->data->data[k] *= cConst;
    f += df;
  }

  return Aeff;
}


void
ComputeParamsFromCoeffs(
			PhenomParams *params,
			PhenomCoeffs *coeffs,
			REAL8        eta,
			REAL8        M)
{
  REAL8 piM;
  piM = LAL_PI * M * LAL_MTSUN_SI;

  params->fMerger = (coeffs->fMerg_a * eta * eta + coeffs->fMerg_b * eta +
		     coeffs->fMerg_c)/piM;
  params->fRing = (coeffs->fRing_a * eta * eta + coeffs->fRing_b * eta +
		   coeffs->fRing_c)/piM;
  params->fCut = (coeffs->fCut_a * eta * eta + coeffs->fCut_b * eta +
		  coeffs->fCut_c)/piM;
  params->sigma = (coeffs->sigma_a * eta * eta * coeffs->sigma_b * eta +
		   coeffs->sigma_c)/piM;

  params->psi0 = coeffs->psi0_x * eta * eta + coeffs->psi0_y * eta +
                   coeffs->psi0_z ;
  params->psi2 = coeffs->psi2_x * eta * eta + coeffs->psi2_y * eta +
                   coeffs->psi2_z ;
  params->psi3 = coeffs->psi3_x * eta * eta + coeffs->psi3_y * eta +
                   coeffs->psi3_z ;
  params->psi4 = coeffs->psi4_x * eta * eta + coeffs->psi4_y * eta +
                   coeffs->psi4_z ;
  params->psi6 = coeffs->psi6_x * eta * eta + coeffs->psi6_y * eta +
                   coeffs->psi6_z ;
  params->psi7 = coeffs->psi7_x * eta * eta + coeffs->psi7_y * eta +
                   coeffs->psi7_z ;

  return;

}
