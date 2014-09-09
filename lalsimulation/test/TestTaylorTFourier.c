/*
 *  Copyright (C) 2014 A. Klein
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

/**
 * \author Antoine Klein
 *
 * \file
 *
 * \brief Testing Fourier domain TaylorT4 and TaylorT2 waveforms.
 *
 * Returns matches between time-domain and frequency-domain waveforms,
 * for TaylorT4 and TaylorT2 approximants (in white noise for simplicity).
 *
 */



#include <lal/LALDatatypes.h>
#include <lal/LALConstants.h>
#include <lal/LALSimInspiral.h>
#include <lal/FrequencySeries.h>
#include <lal/LALStdio.h>

#include <math.h>
#include <fftw3.h>
#include <complex.h>
#include <time.h>







// random double with a flat distribution
REAL8 random_double(REAL8 min, REAL8 max);

// random double with a flat log-distribution
REAL8 random_double_log(REAL8 min, REAL8 max);

// random integer with a flat distribution
INT4 random_int(INT4 min, INT4 max);

// compute normalization factor for time domain waveform
REAL8 normalizeTD(COMPLEX16* hPlusTildeTD, COMPLEX16* hCrossTildeTD, REAL8 Fplus, REAL8 Fcross, INT4 iStart, INT4 jStart, INT4 nSkip, INT4 nMax);

// compute normalization factor for time domain waveform
REAL8 computeMatch(COMPLEX16FrequencySeries* hPlusTildeFD, COMPLEX16FrequencySeries* hCrossTildeFD, COMPLEX16* hPlusTildeTD, COMPLEX16* hCrossTildeTD, REAL8 Fplus, REAL8 Fcross, INT4 iStart, INT4 jStart, INT4 nSkip, INT4 nMax);

// randomizes a set of BHB parameters inside given ranges
void randomize(REAL8* phiRef, REAL8* m1, REAL8* m2, REAL8* r, REAL8* s1x, REAL8* s1y, REAL8* s1z, REAL8* s2x, REAL8* s2y, REAL8* s2z, REAL8* lnhatx, REAL8* lnhaty, REAL8* lnhatz, REAL8* e1x, REAL8* e1y, REAL8* e1z, REAL8* Fplus, REAL8* Fcross, REAL8 m1Min, REAL8 m1Max, REAL8 m2Min, REAL8 m2Max, REAL8 chi1Min, REAL8 chi1Max, REAL8 chi2Min, REAL8 chi2Max, REAL8 rMin, REAL8 rMax);

// compute mean time offset.
// the instantaneous time offset is given by deltaPhi/(2*pi*f)
// this needs the initial dephasing at low frequencies to be small (i.e. 2*pi*deltaT < 0.1)
REAL8 meanTimeOffset(COMPLEX16FrequencySeries* hPlusTildeFD, COMPLEX16FrequencySeries* hCrossTildeFD, COMPLEX16* hPlusTildeTD, COMPLEX16* hCrossTildeTD, REAL8 Fplus, REAL8 Fcross, INT4 iStart, INT4 jStart, INT4 nSkip, INT4 nMax, REAL8 fMin, REAL8 deltaF);

// random double with a flat distribution
REAL8 random_double(REAL8 min, REAL8 max)
{
  return min + ((max - min)*rand())/RAND_MAX;
}

// random double with a flat log-distribution
REAL8 random_double_log(REAL8 min, REAL8 max)
{
  return exp(log(min) + ((log(max) - log(min))*rand())/RAND_MAX);
}

// random integer with a flat distribution
INT4 random_int(INT4 min, INT4 max)
{
  INT4 n_results = max-min+1;
  INT4 mod = RAND_MAX/n_results;
  INT4 max_rand = n_results*mod - 1;
  INT4 rnd = rand();
  while(rnd > max_rand)
  {
    rnd = rand();
  }
  return min + rnd % n_results;
}

// compute normalization factor for time domain waveform
REAL8 normalizeTD(COMPLEX16* hPlusTildeTD, COMPLEX16* hCrossTildeTD, REAL8 Fplus, REAL8 Fcross, INT4 iStart, INT4 jStart, INT4 nSkip, INT4 nMax)
{
  INT4 i, j;
  REAL8 out = 0.;
  COMPLEX16 ht;
  for(i = iStart, j = jStart; j < nMax; i += nSkip, j++)
  {
    ht = Fplus*hPlusTildeTD[i] + Fcross*hCrossTildeTD[i];
    out += creal(ht*conj(ht));
  }
  return 1./sqrt(out);
}

// compute normalization factor for time domain waveform
REAL8 computeMatch(COMPLEX16FrequencySeries* hPlusTildeFD, COMPLEX16FrequencySeries* hCrossTildeFD, COMPLEX16* hPlusTildeTD, COMPLEX16* hCrossTildeTD, REAL8 Fplus, REAL8 Fcross, INT4 iStart, INT4 jStart, INT4 nSkip, INT4 nMax)
{
  INT4 i, j;
  REAL8 out = 0.;
  REAL8 normalizationFD = 0.;
  COMPLEX16 htTD, htFD;
  for(i = iStart, j = jStart; j < nMax; i += nSkip, j++)
  {
    htTD = Fplus*hPlusTildeTD[i] + Fcross*hCrossTildeTD[i];
    htFD = Fplus*hPlusTildeFD->data->data[j] + Fcross*hCrossTildeFD->data->data[j];
    out += creal(htTD*conj(htFD));
    normalizationFD += creal(htFD*conj(htFD));
  }
  return out/sqrt(normalizationFD);
}


// compute mean time offset.
// the instantaneous time offset is given by deltaPhi/(2*pi*f)
// this needs the initial dephasing at low frequencies to be small (i.e. 2*pi*deltaT < 0.1)
REAL8 meanTimeOffset(COMPLEX16FrequencySeries* hPlusTildeFD, COMPLEX16FrequencySeries* hCrossTildeFD, COMPLEX16* hPlusTildeTD, COMPLEX16* hCrossTildeTD, REAL8 Fplus, REAL8 Fcross, INT4 iStart, INT4 jStart, INT4 nSkip, INT4 nMax, REAL8 fMin, REAL8 deltaF)
{
  INT4 i, j;
  REAL8 out = 0.;
  REAL8 dPhi = 0.;
  REAL8 dPhiOld = 0.;
  REAL8 omega;
  COMPLEX16 htTD, htFD;
  for(i = iStart, j = jStart; j < nMax; i += nSkip, j++)
  {
    omega = LAL_TWOPI*(fMin + j*deltaF);
    htTD = Fplus*hPlusTildeTD[i] + Fcross*hCrossTildeTD[i];
    htFD = Fplus*hPlusTildeFD->data->data[j] + Fcross*hCrossTildeFD->data->data[j];
    dPhi = carg(htTD) - carg(htFD);
    while(dPhi - dPhiOld > LAL_PI)
    {
      dPhi -= LAL_TWOPI;
    }
    while(dPhi - dPhiOld < -LAL_PI)
    {
      dPhi += LAL_TWOPI;
    }
    dPhiOld = dPhi;
    out += dPhi/omega;
  }
  return out/(nMax-jStart);
}


// randomizes a set of BHB parameters inside given ranges
void randomize(REAL8* phiRef, REAL8* m1, REAL8* m2, REAL8* r, REAL8* s1x, REAL8* s1y, REAL8* s1z, REAL8* s2x, REAL8* s2y, REAL8* s2z, REAL8* lnhatx, REAL8* lnhaty, REAL8* lnhatz, REAL8* e1x, REAL8* e1y, REAL8* e1z, REAL8* Fplus, REAL8* Fcross, REAL8 m1Min, REAL8 m1Max, REAL8 m2Min, REAL8 m2Max, REAL8 chi1Min, REAL8 chi1Max, REAL8 chi2Min, REAL8 chi2Max, REAL8 rMin, REAL8 rMax)
{
  *phiRef = random_double(0., LAL_TWOPI);

  // random masses flat in log space
  REAL8 M1 = random_double_log(m1Min*LAL_MSUN_SI, m1Max*LAL_MSUN_SI);
  REAL8 M2 = random_double_log(m2Min*LAL_MSUN_SI, m2Max*LAL_MSUN_SI);
  if(M1 > M2)
  {
    *m1 = M1;
    *m2 = M2;
  }
  else
  {
    *m1 = M2;
    *m2 = M1;
  }

  // random distance flat in log space
  *r = random_double_log(rMin, rMax);

  // random spin magnitude flat, random orientation flat on sphere
  REAL8 chi1 = random_double(chi1Min, chi1Max);
  REAL8 mu1 = random_double(-1., 1.);
  REAL8 lambda1 = sqrt(1. - mu1*mu1);
  REAL8 phi1 = random_double(0., LAL_TWOPI);
  *s1x = chi1*lambda1*cos(phi1);
  *s1y = chi1*lambda1*sin(phi1);
  *s1z = chi1*mu1;


  // random spin magnitude flat, random orientation flat on sphere
  REAL8 chi2 = random_double(chi2Min, chi2Max);
  REAL8 mu2 = random_double(-1., 1.);
  REAL8 lambda2 = sqrt(1. - mu2*mu2);
  REAL8 phi2 = random_double(0., LAL_TWOPI);
  *s2x = chi2*lambda2*cos(phi2);
  *s2y = chi2*lambda2*sin(phi2);
  *s2z = chi2*mu2;


  // random orbital angular momentum orientation flat on sphere
  REAL8 mul = random_double(-1., 1.);
  REAL8 lambdal = sqrt(1. - mul*mul);
  REAL8 phil = random_double(0., LAL_TWOPI);
  REAL8 cl = cos(phil);
  REAL8 sl = sin(phil);
  *lnhatx = lambdal*cl;
  *lnhaty = lambdal*sl;
  *lnhatz = mul;

  // e1 has to be perpendicular to lnhat. So, rotate a horizontal vector perpendicular to lnhat around lnhat by a random angle.
  REAL8 alpha = random_double(0., LAL_TWOPI);
  REAL8 ca = cos(alpha);
  REAL8 sa = sin(alpha);

  *e1x = ca*mul*cl - sa*sl;
  *e1y = cl*sa + ca*mul*sl;
  *e1z = -ca*lambdal;

  // detector projection random angle flat
  double twopsi = random_double(0., LAL_TWOPI);
  *Fplus = cos(twopsi);
  *Fcross = sin(twopsi);
}






INT4 main (INT4 argc, char* argv[])
{
  UINT4 startAtZero = 1;
  // RNG seed
  UINT8 seed;

  REAL8 desiredDeltaF;

  // if passed parameters, assume first parameter is the seed
  if(argc > 1)
  {
    seed = atoi(argv[1]);
    if(seed == 0)
    {
      seed = time(0);
    }

    if(argc > 2)
    {
      desiredDeltaF = atof(argv[2]);
    }
    else
    {
      desiredDeltaF = 0.1;
    }
  }
  else
  {
    seed = time(0);
  }

  // seed the RNG
  srand(seed);
  printf("Starting TestTaylorTFourier with seed %" LAL_UINT8_FORMAT "\n" , seed);

  // define different orders
  LALSimInspiralSpinOrder spinO = 7;
  LALSimInspiralTidalOrder tideO = 0;
  INT4 phaseO = 7;
  INT4 amplitudeO = 3;

  // define parameters
  REAL8 phiRef, v0, m1, m2, r, s1x, s1y, s1z, s2x, s2y, s2z, lnhatx, lnhaty, lnhatz, e1x, e1y, e1z, lambda1, lambda2, quadparam1, quadparam2;

  v0 = 1.;

  // deformability parameters turned off
  lambda1 = 0.;
  lambda2 = 0.;
  quadparam1 = 0.;
  quadparam2 = 0.;

  // define instrument projection parameters
  REAL8 Fplus, Fcross;

  // define range for the parameters
  REAL8 m1Min, m1Max, m2Min, m2Max, chi1Min, chi1Max, chi2Min, chi2Max, rMin, rMax;

  m1Min = 1.; // in solar masses
  m1Max = 2.5; // in solar masses
  m2Min = 1.; // in solar masses
  m2Max = 2.5; // in solar masses
  chi1Min = 0.;
  chi1Max = 1.;
  chi2Min = 0.;
  chi2Max = 1.;
  rMin = 3.e24;
  rMax = 3.e24;

  // randomize parameters
  randomize(&phiRef, &m1, &m2, &r, &s1x, &s1y, &s1z, &s2x, &s2y, &s2z, &lnhatx, &lnhaty, &lnhatz, &e1x, &e1y, &e1z, &Fplus, &Fcross, m1Min, m1Max, m2Min, m2Max,  chi1Min, chi1Max, chi2Min, chi2Max, rMin, rMax);

  // misc parameters
  REAL8 tTot;
  REAL8 deltaF;
  INT4 nSamples;
  REAL8 fMin;
  REAL8 fMax;
  REAL8 tInit;
  REAL8 desiredFMin;
  UINT4 i;

  fMax = 2000.;
  desiredFMin = 10.;

  REAL8 M = m1 + m2;
  REAL8 eta = (m1*m2)/(M*M);
  // define time series for time domain waveform
  REAL8TimeSeries* hPlus = 0;
  REAL8TimeSeries* hCross = 0;

  // define extra parameters for TD WF
  REAL8 deltaT;

  REAL8 fStart = 9.;
  REAL8 fRef = 70.;
  INT4 phiRefAtEnd = 0;

  // more misc variables
  COMPLEX16* hPlusTildeTD;
  COMPLEX16* hCrossTildeTD;

  REAL8 factTukey, t1, t2, t3, t4;
  REAL8 sinfact, t;

  fftw_complex* ftilde_data;
  fftw_plan plan;

  INT4 iStart, jStart, nSkip;
  REAL8 rStop;

  INT4 nMax;


  // define frequency series for FD WFs
  COMPLEX16FrequencySeries* hPlusTildeFD = 0;
  COMPLEX16FrequencySeries* hCrossTildeFD = 0;

  INT4 kMax;
  REAL8 normalizationTD;
  REAL8 meanDeltaT;
  REAL8 match;


  printf("TaylorT4:\nStarting time domain...\n");

  // compute TD WF with deltaT=1 to get signal time duration
  deltaT = 1.;
  XLALSimInspiralSpinTaylorT4(&hPlus, &hCross, phiRef, v0, deltaT, m1, m2, fStart, fRef, r, s1x, s1y, s1z, s2x, s2y, s2z, lnhatx, lnhaty, lnhatz, e1x, e1y, e1z, lambda1, lambda2, quadparam1, quadparam2, spinO, tideO, phaseO, amplitudeO);

  // compute deltaT necessary for Nyquist frequency at fMax
  tTot = -(hPlus->epoch.gpsSeconds + 1.e-9*hPlus->epoch.gpsNanoSeconds);
  deltaF = 1./tTot;
  nSamples = 2.*fMax/deltaF;
  deltaT = tTot/nSamples;

  // compute TD WF with good deltaT
  XLALDestroyREAL8TimeSeries(hPlus);
  XLALDestroyREAL8TimeSeries(hCross);
  XLALSimInspiralSpinTaylorT4(&hPlus, &hCross, phiRef, v0, deltaT, m1, m2, fStart, fRef, r, s1x, s1y, s1z, s2x, s2y, s2z, lnhatx, lnhaty, lnhatz, e1x, e1y, e1z, lambda1, lambda2, quadparam1, quadparam2, spinO, tideO, phaseO, amplitudeO);

  // define and allocate DFT of TD WF
  hPlusTildeTD = (COMPLEX16*)malloc(sizeof(COMPLEX16)*(hPlus->data->length));
  hCrossTildeTD = (COMPLEX16*)malloc(sizeof(COMPLEX16)*(hCross->data->length));

  tInit = hPlus->epoch.gpsSeconds + 1.e-9*hPlus->epoch.gpsNanoSeconds;
  deltaF = 1./(deltaT*(hPlus->data->length));

  // define Tukey window parameters
  factTukey = (-5.*LAL_G_SI*M)/(((256.*eta*LAL_C_SI)*LAL_C_SI)*LAL_C_SI);
  t1 = tInit;
  t2 = factTukey*pow((LAL_PI*LAL_G_SI)*M*9.5/(((((REAL8)LAL_C_SI)*LAL_C_SI)*LAL_C_SI)), -8./3.); // 9.5 Hertz Newtonian
  t3 = factTukey*50625.; // 15M separation Newtonian
  t4 = 0.;

  // apply Tukey window
  i = 0;
  t = tInit;
  while(t < t2)
  {
    sinfact = sin(LAL_PI_2*(t - t1)/(t2 - t1));
    hPlus->data->data[i] *= sinfact*sinfact;
    hCross->data->data[i] *= sinfact*sinfact;
    i++;
    t = tInit + i*deltaT;
  }

  i = hPlus->data->length-1;
  t = tInit + i*deltaT;
  while(t > t3)
  {
    sinfact = sin(LAL_PI_2*(t - t4)/(t3 - t4));
    hPlus->data->data[i] *= sinfact*sinfact;
    hCross->data->data[i] *= sinfact*sinfact;
    i--;
    t = tInit + i*deltaT;
  }

  // compute the DFTs, applying the necessary time shift
  ftilde_data = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(hPlus->data->length));
  plan = fftw_plan_dft_r2c_1d(hPlus->data->length, hPlus->data->data, ftilde_data, FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan);

  for(i = 0; i <= hPlus->data->length >> 1; i++)
  {
    hPlusTildeTD[i] = deltaT*cexp(-I*LAL_TWOPI*i*deltaF*tInit)*ftilde_data[i];
  }

  plan = fftw_plan_dft_r2c_1d(hCross->data->length, hCross->data->data, ftilde_data, FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan);

  for(i = 0; i <= hCross->data->length >> 1; i++)
  {
    hCrossTildeTD[i] = deltaT*cexp(-I*LAL_TWOPI*i*deltaF*tInit)*ftilde_data[i];
  }

  fftw_free(ftilde_data);


  // compute fMin and deltaF for the frequency domain WFs, close to desiredFMin and desiredDeltaF but a multiple of deltaF for the TD WF
  iStart = (INT4)ceil(desiredFMin/deltaF);
  nSkip = (INT4)ceil(desiredDeltaF/deltaF);
  iStart -= iStart % nSkip;
  fMin = iStart*deltaF;
  deltaF *= nSkip;
  if(startAtZero)
  {
    jStart = iStart/nSkip;
  }
  else
  {
    jStart = 0;
  }

  // set maximum frequency for the comparison
  rStop = 20.; // frequency of r = 20M
  fMax = LAL_C_SI*(LAL_C_SI*(LAL_C_SI/(LAL_PI*LAL_G_SI*M*rStop*sqrt(rStop))));

  nMax = (INT4)floor((fMax - fMin)/deltaF) + jStart;


  // normalize TD WF
  normalizationTD = normalizeTD(hPlusTildeTD, hCrossTildeTD, Fplus, Fcross, iStart, jStart, nSkip, nMax);

  printf("Matches:\n");

  for(kMax = 0; kMax <= 10; kMax++) // loop over kMax
  {
    //compute FD WF
    XLALSimInspiralSpinTaylorT4Fourier(&hPlusTildeFD, &hCrossTildeFD, fMin, fMax, deltaF, kMax, phiRef, v0, m1, m2, fStart, fRef, r, s1x, s1y, s1z, s2x, s2y, s2z, lnhatx, lnhaty, lnhatz, e1x, e1y, e1z, lambda1, lambda2, quadparam1, quadparam2, spinO, tideO, phaseO, amplitudeO, phiRefAtEnd);

    // XLALSimInspiralSpinTaylorT4 and XLALSimInspiralSpinTaylorT4Fourier return with a slight time offset between the two. We get rid of that.
    if(startAtZero)
    {
      meanDeltaT = meanTimeOffset(hPlusTildeFD, hCrossTildeFD, hPlusTildeTD, hCrossTildeTD, Fplus, Fcross, iStart, jStart, nSkip, nMax, 0., deltaF);
      for(i = 0; i < hPlusTildeFD->data->length; i++)
      {
        hPlusTildeFD->data->data[i] *= cexp(I*LAL_TWOPI*(i*deltaF)*meanDeltaT);
        hCrossTildeFD->data->data[i] *= cexp(I*LAL_TWOPI*(i*deltaF)*meanDeltaT);
      }
    }
    else
    {
      meanDeltaT = meanTimeOffset(hPlusTildeFD, hCrossTildeFD, hPlusTildeTD, hCrossTildeTD, Fplus, Fcross, iStart, jStart, nSkip, nMax, fMin, deltaF);
      for(i = 0; i < hPlusTildeFD->data->length; i++)
      {
        hPlusTildeFD->data->data[i] *= cexp(I*LAL_TWOPI*(fMin + i*deltaF)*meanDeltaT);
        hCrossTildeFD->data->data[i] *= cexp(I*LAL_TWOPI*(fMin + i*deltaF)*meanDeltaT);
      }
    }

    // compute match
    match = normalizationTD*computeMatch(hPlusTildeFD, hCrossTildeFD, hPlusTildeTD, hCrossTildeTD, Fplus, Fcross, iStart, jStart, nSkip, nMax);
    printf("kMax = %2d: %.15f\n", kMax, match);

    XLALDestroyCOMPLEX16FrequencySeries(hPlusTildeFD);
    XLALDestroyCOMPLEX16FrequencySeries(hCrossTildeFD);
  }

  XLALDestroyREAL8TimeSeries(hPlus);
  XLALDestroyREAL8TimeSeries(hCross);
  free(hPlusTildeTD);
  free(hCrossTildeTD);








  printf("\nTaylorT2:\nStarting time domain...\n");

  // compute TD WF with deltaT=1 to get signal time duration
  deltaT = 1.;
  XLALSimInspiralSpinTaylorT2(&hPlus, &hCross, phiRef, v0, deltaT, m1, m2, fStart, fRef, r, s1x, s1y, s1z, s2x, s2y, s2z, lnhatx, lnhaty, lnhatz, e1x, e1y, e1z, lambda1, lambda2, quadparam1, quadparam2, spinO, tideO, phaseO, amplitudeO);

  // compute deltaT necessary for Nyquist frequency at fMax
  tTot = -(hPlus->epoch.gpsSeconds + 1.e-9*hPlus->epoch.gpsNanoSeconds);
  deltaF = 1./tTot;
  nSamples = 2.*fMax/deltaF;
  deltaT = tTot/nSamples;

  // compute TD WF with good deltaT
  XLALDestroyREAL8TimeSeries(hPlus);
  XLALDestroyREAL8TimeSeries(hCross);
  XLALSimInspiralSpinTaylorT2(&hPlus, &hCross, phiRef, v0, deltaT, m1, m2, fStart, fRef, r, s1x, s1y, s1z, s2x, s2y, s2z, lnhatx, lnhaty, lnhatz, e1x, e1y, e1z, lambda1, lambda2, quadparam1, quadparam2, spinO, tideO, phaseO, amplitudeO);

  // define and allocate DFT of TD WF
  hPlusTildeTD = (COMPLEX16*)malloc(sizeof(COMPLEX16)*(hPlus->data->length));
  hCrossTildeTD = (COMPLEX16*)malloc(sizeof(COMPLEX16)*(hCross->data->length));

  tInit = hPlus->epoch.gpsSeconds + 1.e-9*hPlus->epoch.gpsNanoSeconds;
  deltaF = 1./(deltaT*(hPlus->data->length));

  // define Tukey window parameters
  factTukey = (-5.*LAL_G_SI*M)/(((256.*eta*LAL_C_SI)*LAL_C_SI)*LAL_C_SI);
  t1 = tInit;
  t2 = factTukey*pow((LAL_PI*LAL_G_SI)*M*9.5/(((((REAL8)LAL_C_SI)*LAL_C_SI)*LAL_C_SI)), -8./3.); // 9.5 Hertz Newtonian
  t3 = factTukey*25.*25.*25.*25.; // 25M separation Newtonian. For TaylorT2 we need a huge half-Hann window at the end for some reason.
  t4 = 0.;

  // apply Tukey window
  i = 0;
  t = tInit;
  while(t < t2)
  {
    sinfact = sin(LAL_PI_2*(t - t1)/(t2 - t1));
    hPlus->data->data[i] *= sinfact*sinfact;
    hCross->data->data[i] *= sinfact*sinfact;
    i++;
    t = tInit + i*deltaT;
  }

  i = hPlus->data->length-1;
  t = tInit + i*deltaT;
  while(t > t3)
  {
    sinfact = sin(LAL_PI_2*(t - t4)/(t3 - t4));
    hPlus->data->data[i] *= sinfact*sinfact;
    hCross->data->data[i] *= sinfact*sinfact;
    i--;
    t = tInit + i*deltaT;
  }

  // compute the DFTs, applying the necessary time shift
  ftilde_data = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(hPlus->data->length));
  plan = fftw_plan_dft_r2c_1d(hPlus->data->length, hPlus->data->data, ftilde_data, FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan);

  for(i = 0; i <= hPlus->data->length >> 1; i++)
  {
    hPlusTildeTD[i] = deltaT*cexp(-I*LAL_TWOPI*i*deltaF*tInit)*ftilde_data[i];
  }

  plan = fftw_plan_dft_r2c_1d(hCross->data->length, hCross->data->data, ftilde_data, FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan);

  for(i = 0; i <= hCross->data->length >> 1; i++)
  {
    hCrossTildeTD[i] = deltaT*cexp(-I*LAL_TWOPI*i*deltaF*tInit)*ftilde_data[i];
  }

  fftw_free(ftilde_data);


  // compute fMin and deltaF for the frequency domain WFs, close to desiredFMin and desiredDeltaF but a multiple of deltaF for the TD WF
  iStart = (INT4)ceil(desiredFMin/deltaF);
  nSkip = (INT4)ceil(desiredDeltaF/deltaF);
  iStart -= iStart % nSkip;
  fMin = iStart*deltaF;
  deltaF *= nSkip;
  if(startAtZero)
  {
    jStart = iStart/nSkip;
  }
  else
  {
    jStart = 0;
  }

  // set maximum frequency for the comparison
  rStop = 25.; // frequency of r = 20M
  fMax = LAL_C_SI*(LAL_C_SI*(LAL_C_SI/(LAL_PI*LAL_G_SI*M*rStop*sqrt(rStop))));

  nMax = (INT4)floor((fMax - fMin)/deltaF) + jStart;


  // define frequency series for FD WFs
  hPlusTildeFD = 0;
  hCrossTildeFD = 0;

  // normalize TD WF
  normalizationTD = normalizeTD(hPlusTildeTD, hCrossTildeTD, Fplus, Fcross, iStart, jStart, nSkip, nMax);

  printf("Matches:\n");

  for(kMax = 0; kMax <= 10; kMax++) // loop over kMax
  {
    //compute FD WF
    XLALSimInspiralSpinTaylorT2Fourier(&hPlusTildeFD, &hCrossTildeFD, fMin, fMax, deltaF, kMax, phiRef, v0, m1, m2, fStart, fRef, r, s1x, s1y, s1z, s2x, s2y, s2z, lnhatx, lnhaty, lnhatz, e1x, e1y, e1z, lambda1, lambda2, quadparam1, quadparam2, spinO, tideO, phaseO, amplitudeO, phiRefAtEnd);

    // XLALSimInspiralSpinTaylorT2 and XLALSimInspiralSpinTaylorT2Fourier return with a slight time offset between the two. We get rid of that.
    if(startAtZero)
    {
      meanDeltaT = meanTimeOffset(hPlusTildeFD, hCrossTildeFD, hPlusTildeTD, hCrossTildeTD, Fplus, Fcross, iStart, jStart, nSkip, nMax, 0., deltaF);
      for(i = 0; i < hPlusTildeFD->data->length; i++)
      {
        hPlusTildeFD->data->data[i] *= cexp(I*LAL_TWOPI*(i*deltaF)*meanDeltaT);
        hCrossTildeFD->data->data[i] *= cexp(I*LAL_TWOPI*(i*deltaF)*meanDeltaT);
      }
    }
    else
    {
      meanDeltaT = meanTimeOffset(hPlusTildeFD, hCrossTildeFD, hPlusTildeTD, hCrossTildeTD, Fplus, Fcross, iStart, jStart, nSkip, nMax, fMin, deltaF);
      for(i = 0; i < hPlusTildeFD->data->length; i++)
      {
        hPlusTildeFD->data->data[i] *= cexp(I*LAL_TWOPI*(fMin + i*deltaF)*meanDeltaT);
        hCrossTildeFD->data->data[i] *= cexp(I*LAL_TWOPI*(fMin + i*deltaF)*meanDeltaT);
      }
    }

    // compute match
    match = normalizationTD*computeMatch(hPlusTildeFD, hCrossTildeFD, hPlusTildeTD, hCrossTildeTD, Fplus, Fcross, iStart, jStart, nSkip, nMax);
    printf("kMax = %2d: %.15f\n", kMax, match);

    XLALDestroyCOMPLEX16FrequencySeries(hPlusTildeFD);
    XLALDestroyCOMPLEX16FrequencySeries(hCrossTildeFD);
  }

  XLALDestroyREAL8TimeSeries(hPlus);
  XLALDestroyREAL8TimeSeries(hCross);
  free(hPlusTildeTD);
  free(hCrossTildeTD);

  fftw_cleanup();


  return XLAL_SUCCESS;
}
