/*
*  Copyright (C) 2008 Lucia Santamaria, Robert Adam Mercer, Badri Krishnan, P. Ajith
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

/*
$Id$
*/

/*
This code implements the phenomenological formula for hybrid PN-NR waves
given in P. Ajith et al. LIGO-P070111-00-Z, AEI-2007-143
e-Print: arXiv:0710.2335 [gr-qc]
 */

#include <math.h>
#include <stdio.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/LALInspiral.h>
#include <lal/PrintFTSeries.h>
#include <getopt.h>
#include <string.h>

#include <lal/LALConfig.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/LALDatatypes.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/AVFactories.h>
#include <lal/NRWaveIO.h>
#include <lal/NRWaveInject.h>
#include <lal/LIGOLwXMLRead.h>
#include <lal/Inject.h>
#include <lal/FileIO.h>
#include <lal/Units.h>
#include <lal/FrequencySeries.h>
#include <lal/TimeSeries.h>
#include <lal/TimeFreqFFT.h>
#include <lal/VectorOps.h>
#include <lal/LALDetectors.h>
#include <lal/LALFrameIO.h>
#include <lal/ComplexFFT.h>

/** DEFINE RCS ID STRING **/
NRCSID (XLALINSPIRALHYBRIDP1C, "$Id$");
#define CVS_ID_STRING "$Id$"
#define CVS_NAME_STRING "$Name$"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define PROGRAM_NAME "lalapps_phenom"

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

static void print_usage( CHAR *program );

void GetPhenomCoeffs(
		     PhenomCoeffs *co);

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

void
XLALComputeComplexVector(
			 COMPLEX8Vector **outPlus,
			 COMPLEX8Vector **outCross,
			 REAL4FrequencySeries *Ampl,
			 REAL4FrequencySeries *Phase);

REAL4Vector *
XLALComputeFreq(
		REAL4TimeSeries *hp,
		REAL4TimeSeries *hc);

REAL4TimeSeries *
XLALCutAtFreq(
	      REAL4TimeSeries *h,
	      REAL4Vector     *freq,
	      REAL8           cutFreq);

void LALPrintHPlusCross(
			REAL4TimeSeries *hp,
			REAL4TimeSeries *hc,
			CHAR            *out );


/* Main Program */
INT4 main ( INT4 argc, CHAR *argv[] ) {

  static LALStatus status;
  INT4 c;
  UINT4 i;
  REAL8 dt, totTime;
  REAL8 sampleRate = -1;
  REAL8 totalMass = -1, massRatio = -1;
  REAL8 lowFreq = -1, df, fLow;
  CHAR  *outFile = NULL, *outFileLong = NULL, tail[50];
  size_t optarg_len;
  REAL8 eta;
  REAL8 newtonianChirpTime, PN1ChirpTime, mergTime;
  UINT4 numPts;
  LIGOTimeGPS epoch;
  REAL8 offset;

  PhenomCoeffs coeffs;
  PhenomParams params;

  REAL4FrequencySeries     *Aeff   = NULL, *Phieff  = NULL;
  COMPLEX8Vector           *uFPlus = NULL, *uFCross = NULL;

  COMPLEX8 num;

  REAL4Vector      *hPlus = NULL, *hCross = NULL;
  REAL4TimeSeries  *hP = NULL, *hC = NULL;
  /*  REAL4TimeSeries  *hP = NULL, *hC = NULL;*/
  REAL4FFTPlan     *prevPlus = NULL, *prevCross = NULL;

  /*REAL4Vector      *Freq = NULL;*/

  UINT4 windowLength;
  INT4 hPLength;
  REAL8 linearWindow;

  /* getopt arguments */
  struct option long_options[] =
  {
    {"mass-ratio",              required_argument, 0,                'q'},
    {"low-freq (Hz)",           required_argument, 0,                'f'},
    {"total-mass (M_sun)",      required_argument, 0,                'm'},
    {"sample-rate",             required_argument, 0,                's'},
    {"output-file",             required_argument, 0,                'o'},
    {"help",                    no_argument,       0,                'h'},
    {"version",                 no_argument,       0,                'V'},
    {0, 0, 0, 0}
  };

  /* parse the arguments */
  while ( 1 )
  {
    /* getopt_long stores long option here */
    int option_index = 0;

    /* parse command line arguments */
    c = getopt_long_only( argc, argv, "q:t:d:hV",
        long_options, &option_index );

    /* detect the end of the options */
    if ( c == -1 )
      {
	break;
      }

    switch ( c )
    {
      case 0:
        fprintf( stderr, "Error parsing option '%s' with argument '%s'\n",
            long_options[option_index].name, optarg );
        exit( 1 );
        break;

      case 'h':
        /* help message */
        print_usage( argv[0] );
        exit( 0 );
        break;

      case 'V':
        /* print version information and exit */
        fprintf( stdout, "%s - Compute Ajith's Phenomenological Waveforms " \
		 "(arXiv:0710.2335) and output them to a plain text file\n" \
            "CVS Version: %s\nCVS Tag: %s\n", PROGRAM_NAME, CVS_ID_STRING, \
            CVS_NAME_STRING );
        exit( 0 );
        break;

      case 'q':
        /* set mass ratio */
        massRatio = atof( optarg );
        break;

      case 'f':
        /* set low freq */
        lowFreq = atof( optarg );
        break;

      case 'm':
        /* set total mass */
        totalMass = atof( optarg );
        break;

      case 's':
        /* set sample rate */
        sampleRate = atof( optarg );
        break;

      case 'o':
	/* set name of output file */
        optarg_len = strlen(optarg) + 1;
        outFile = (CHAR *)calloc(optarg_len, sizeof(CHAR));
        memcpy(outFile, optarg, optarg_len);
	break;

      case '?':
        print_usage( argv[0] );
        exit( 1 );
        break;

      default:
        fprintf( stderr, "ERROR: Unknown error while parsing options\n" );
        print_usage( argv[0] );
        exit( 1 );
    }
  }

  if ( optind < argc )
  {
    fprintf( stderr, "ERROR: Extraneous command line arguments:\n" );
    while ( optind < argc )
      {
	fprintf ( stderr, "%s\n", argv[optind++] );
      }
    exit( 1 );
  }


  /* * * * * * * * */
  /* Main Program  */
  /* * * * * * * * */

  eta = massRatio / pow(1. + massRatio, 2.);

  /* This freq low is the one used for the FFT */
  /* fLow = 2.E-3/(totalMass*LAL_MTSUN_SI); */
  fLow = lowFreq;       /* Changed by Ajith. 5 May 2008 */

  /* Phenomenological coefficients as in Ajith et. al */
  GetPhenomCoeffsLongJena( &coeffs );

  /* Compute phenomenologial parameters */
  ComputeParamsFromCoeffs( &params, &coeffs, eta, totalMass );


  /* Check validity of arguments */

  /* check we have freqs */
  if ( totalMass < 0 )
    {
      fprintf( stderr, "ERROR: --total-mass must be specified\n" );
      exit( 1 );
    }

  /* check we have mass ratio and delta t*/
  if ( massRatio < 0 )
    {
      fprintf( stderr, "ERROR: --mass-ratio must be specified\n" );
      exit( 1 );
    }
  if ( lowFreq < 0 )
    {
      fprintf( stderr, "ERROR: --low-freq must be specified\n" );
      exit( 1 );
    }
  if ( sampleRate < 0 )
    {
      fprintf( stderr, "ERROR: --sample-rate must be specified\n" );
      exit( 1 );
    }
  if ( lowFreq > params.fCut )
    {
      fprintf( stderr, "\nERROR in --low-freq\n"\
	       "The value chosen for the low frequency is larger "\
	       "than the frequency at the merger.\n"
	       "Frequency at the merger: %4.2f Hz\nPick either a lower value"\
	       " for --low-freq or a lower total mass\n\n", params.fCut);
      exit(1);
    }
  if ( lowFreq < fLow )
    {
      fprintf( stderr, "\nERROR in --low-freq\n"\
	       "The value chosen for the low frequency is lower "\
	       "than the lowest frequency computed\nby the implemented FFT.\n"
	       "Lowest frequency allowed: %4.2f Hz\nPick either a higher value"\
	       " for --low-freq or a higher total mass\n\n", fLow);
      exit(1);
    }
  if ( outFile == NULL )
    {
      fprintf( stderr, "ERROR: --output-file must be specified\n" );
      exit( 1 );
    }

  /* Complete file name with details of the input variables */
  sprintf(tail, "%s-Phenom_M%3.1f_R%2.1f.dat", outFile, totalMass, massRatio);
  optarg_len = strlen(tail) + strlen(outFile) + 1;
  outFileLong = (CHAR *)calloc(optarg_len, sizeof(CHAR));
  strcpy(outFileLong, tail);

  /* check sample rate is enough */
  if (sampleRate > 4.*params.fCut)   /* Changed by Ajith. 5 May 2008 */
    {
      dt = 1./sampleRate;
    }
  else
    {
      sampleRate = 4.*params.fCut;
      dt = 1./sampleRate;
    }

  /* Estimation of the time duration of the binary           */
  /* See Sathya (1994) for the Newtonian and PN1 chirp times */
  /* The merger time is overestimated                        */
  newtonianChirpTime =
    (5./(256.*eta))*pow(totalMass*LAL_MTSUN_SI,-5./3.)*pow(LAL_PI*fLow,-8./3.);
  PN1ChirpTime =
    5.*(743.+924.*eta)/(64512.*eta*totalMass*LAL_MTSUN_SI*pow(LAL_PI*fLow,2.));
  mergTime = 2000.*totalMass*LAL_MTSUN_SI;
  totTime = 1.2 * (newtonianChirpTime + PN1ChirpTime + mergTime);

  numPts = (UINT4) ceil(totTime/dt);
  df = 1/(numPts * dt);

  /* Compute Amplitude and Phase from the paper (Eq. 4.19) */
  Aeff = XLALHybridP1Amplitude(&params, fLow, df, eta, totalMass, numPts/2+1);
  Phieff = XLALHybridP1Phase(&params, fLow, df, eta, totalMass, numPts/2 +1);

  /* Construct u(f) = Aeff*e^(i*Phieff) */
  XLALComputeComplexVector(&uFPlus, &uFCross, Aeff, Phieff);

  /* Scale this to units of M */
  for (i = 0; i < numPts/2 + 1; i++) {
    num = uFPlus->data[i];
    num.re *= 1./(dt*totalMass*LAL_MTSUN_SI);
    num.im *= 1./(dt*totalMass*LAL_MTSUN_SI);
    uFPlus->data[i] = num;
    num = uFCross->data[i];
    num.re *= 1./(dt*totalMass*LAL_MTSUN_SI);
    num.im *= 1./(dt*totalMass*LAL_MTSUN_SI);
    uFCross->data[i] = num;
  }

  /* Inverse Fourier transform */
  LALCreateReverseREAL4FFTPlan( &status, &prevPlus, numPts, 0 );
  LALCreateReverseREAL4FFTPlan( &status, &prevCross, numPts, 0 );
  hPlus = XLALCreateREAL4Vector(numPts);
  hCross = XLALCreateREAL4Vector(numPts);

  LALReverseREAL4FFT( &status, hPlus, uFPlus, prevPlus );
  LALReverseREAL4FFT( &status, hCross, uFCross, prevCross );

  /* The LAL implementation of the FFT omits the factor 1/n */
  for (i = 0; i < numPts; i++) {
    hPlus->data[i] /= numPts;
    hCross->data[i] /= numPts;
  }

  /* Create TimeSeries to store more info about the waveforms */
  /* Note: it could be done easier using LALFreqTimeFFT instead of ReverseFFT */
  epoch.gpsSeconds = 0;
  epoch.gpsNanoSeconds = 0;

  hP =
    XLALCreateREAL4TimeSeries("", &epoch, 0, dt, &lalDimensionlessUnit, numPts);
  hP->data = hPlus;
  hC =
    XLALCreateREAL4TimeSeries("", &epoch, 0, dt, &lalDimensionlessUnit, numPts);
  hC->data = hCross;

  /* Cutting off the part of the waveform with f < fLow */
  /*  Freq = XLALComputeFreq( hP, hC);
      hP = XLALCutAtFreq( hP, Freq, lowFreq);
      hC = XLALCutAtFreq( hC, Freq, lowFreq); */

  /* multiply the last few samples of the time-series by a linearly
   * dropping window function in order to avid edges in the data
   * Added by Ajith 6 May 2008 */

  hPLength = hP->data->length;
  windowLength = (UINT4) (20.*totalMass * LAL_MTSUN_SI/dt);
  for (i=1; i<= windowLength; i++){
    linearWindow =  (i-1.)/windowLength;
    hP->data->data[hPLength-i] *= linearWindow;
    hC->data->data[hPLength-i] *= linearWindow;
  }

  /* Convert t column to units of (1/M) */
  /*  offset *= (1./(totalMass * LAL_MTSUN_SI));
      hP->deltaT *= (1./(totalMass * LAL_MTSUN_SI)); */

  /* Set t = 0 at the merger (defined as the max of the NR wave) */
  XLALFindNRCoalescenceTimeFromhoft( &offset, hP);
  XLALGPSAdd( &(hP->epoch), -offset);
  XLALGPSAdd( &(hC->epoch), -offset);

  /* Print waveforms to file */
  LALPrintHPlusCross( hP, hC, outFileLong );

  /* Free Memory */

  XLALDestroyREAL4FrequencySeries(Aeff);
  XLALDestroyREAL4FrequencySeries(Phieff);

  XLALDestroyREAL4FFTPlan(prevPlus);
  XLALDestroyREAL4FFTPlan(prevCross);

  XLALDestroyCOMPLEX8Vector(uFPlus);
  XLALDestroyCOMPLEX8Vector(uFCross);
  XLALDestroyREAL4TimeSeries(hP);
  XLALDestroyREAL4TimeSeries(hC);
  /*  XLALDestroyREAL4TimeSeries(hP); */
  /*  XLALDestroyREAL4TimeSeries(hC); */

  return(0);

}


/* My funtions */

void LALPrintHPlusCross(
			REAL4TimeSeries *hp,
			REAL4TimeSeries *hc,
			CHAR            *out )
{
  UINT4 i, n;
  FILE *file;
  REAL8 dt, off=0;

  n = hp->data->length;
  dt = hp->deltaT;
  XLALGPSSetREAL8( &(hp->epoch), off);

  file = LALFopen(out, "w");
  fprintf (file, "#   t[sec]\t    h_+\t\t   h_x\n");
  for (i=0; i < n; i++)
    {
      fprintf (file, "%e\t%e\t%e\n",
	       i*dt+off, hp->data->data[i], hc->data->data[i]);
    }
}


void
XLALComputeComplexVector(
			 COMPLEX8Vector **outPlus,
			 COMPLEX8Vector **outCross,
			 REAL4FrequencySeries *Ampl,
			 REAL4FrequencySeries *Phase)
{
  COMPLEX8Vector *uPlus = NULL, *uCross = NULL;
  COMPLEX8 num;
  REAL4 Re, Im;
  UINT4 k, n;

  n = Ampl->data->length;

  uPlus = LALCalloc(1, sizeof(*uPlus));
  uPlus = XLALCreateCOMPLEX8Vector(n);
  uCross = LALCalloc(1, sizeof(*uCross));
  uCross = XLALCreateCOMPLEX8Vector(n);

  for( k = 0 ; k < n ; k++ ) {
    Re = Ampl->data->data[k] * cos(Phase->data->data[k]);
    Im = Ampl->data->data[k] * sin(Phase->data->data[k]);
    num.re = Re;
    num.im = Im;
    uPlus->data[k] = num;

    Re = Ampl->data->data[k] * cos(Phase->data->data[k] - LAL_PI/2.);
    Im = Ampl->data->data[k] * sin(Phase->data->data[k] - LAL_PI/2.);
    num.re = Re;
    num.im = Im;
    uCross->data[k] = num;
  }

  *outPlus = uPlus;
  *outCross = uCross;

}


REAL4TimeSeries *
XLALCutAtFreq(
	      REAL4TimeSeries *h,
	      REAL4Vector     *freq,
	      REAL8           cutFreq)
{
  REAL8 dt;
  UINT4 k, k0, kMid, len;
  REAL4 currentFreq;
  UINT4 newLen;

  REAL4TimeSeries *newH = NULL;

  LIGOTimeGPS epoch;

  len = freq->length;
  dt = h->deltaT;

  /* Since the boundaries of this freq vector are likely to have   */
  /* FFT crap, let's scan the freq values starting from the middle */
  kMid = len/2;
  currentFreq = freq->data[kMid];
  k = kMid;

  /* freq is an increasing function of time */
  /* If we are above the cutFreq we move to the left; else to the right */
  if (currentFreq > cutFreq && k > 0)
    {
      while(currentFreq > cutFreq)
	{
	  currentFreq = freq->data[k];
	  k--;
	}
      k0 = k;
    }
  else
    {
      while(currentFreq < cutFreq && k < len)
	{
	  currentFreq = freq->data[k];
	  k++;
	}
      k0 = k;
    }

  newLen = len - k0;

  /* Allocate memory for the frequency series */
  epoch.gpsSeconds = 0;
  epoch.gpsNanoSeconds = 0;
  newH =
    XLALCreateREAL4TimeSeries("", &epoch, 0, dt, &lalDimensionlessUnit, newLen);

  for(k = 0; k < newLen; k++)
    {
      newH->data->data[k] = h->data->data[k0 + k];
    }

  return newH;

}


REAL4Vector *
XLALComputeFreq(
		REAL4TimeSeries *hp,
		REAL4TimeSeries *hc)
{
  REAL4Vector *Freq = NULL;
  REAL4Vector *hpDot = NULL, *hcDot = NULL;
  UINT4 k, len;
  REAL8 dt;

  len = hp->data->length;
  dt = hp->deltaT;
  Freq = LALCalloc(1, sizeof(*Freq));
  Freq= XLALCreateREAL4Vector(len);

  hpDot = LALCalloc(1, sizeof(*hpDot));
  hpDot= XLALCreateREAL4Vector(len);
  hcDot = LALCalloc(1, sizeof(*hcDot));
  hcDot= XLALCreateREAL4Vector(len);

  /* Construct the dot vectors (2nd order differencing) */
  hpDot->data[0] = 0.0;
  hpDot->data[len] = 0.0;
  hcDot->data[0] = 0.0;
  hcDot->data[len] = 0.0;
  for( k = 1; k < len-1; k++)
    {
      hpDot->data[k] = 1./(2.*dt) *
	(hp->data->data[k+1]-hp->data->data[k-1]);
      hcDot->data[k] = 1./(2.*dt) *
	(hc->data->data[k+1]-hc->data->data[k-1]);
    }

  /* Compute frequency using the fact that  */
  /*h(t) = A(t) e^(i Phi) = Re(h) + i Im(h) */
  for( k = 0; k < len; k++)
    {
      Freq->data[k] = hcDot->data[k] * hp->data->data[k] -
	hpDot->data[k] * hc->data->data[k];
      Freq->data[k] /= LAL_TWOPI;
      Freq->data[k] /= (pow(hp->data->data[k],2.) + pow(hc->data->data[k], 2.));
    }

  return Freq;

}


REAL4FrequencySeries *
XLALHybridP1Phase(
		  PhenomParams  *params,
		  REAL8         fLow,
		  REAL8         df,
		  REAL8         eta,
		  REAL8         M,
		  UINT4         n )
{
  /* Computes effective phase as in arXiv:0710.2335 [gr-qc] */

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
  REAL8 piM;
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

  piM = LAL_PI * M * LAL_MTSUN_SI;

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
}


/* Coeffs from the paper - matching with Jena short waveforms */
void
GetPhenomCoeffs(
		PhenomCoeffs *co)
{
  co->fMerg_a = 2.9740e-1; co->fMerg_b = 4.4810e-2; co->fMerg_c = 9.5560e-2;
  co->fRing_a = 5.9411e-1; co->fRing_b = 8.9794e-2;; co->fRing_c = 1.9111e-1;
  co->sigma_a = 5.0801e-1;; co->sigma_b = 7.7515e-2; co->sigma_c = 2.2369e-2;
  co->fCut_a = 8.4845e-1; co->fCut_b = 1.2848e-1; co->fCut_c = 2.7299e-1;

  co->psi0_x = 1.7516e-1; co->psi0_y = 7.9483e-2; co->psi0_z = -7.2390e-2;
  co->psi2_x = -5.1571e1; co->psi2_y = -1.7595e1; co->psi2_z = 1.3253e1;
  co->psi3_x = 6.5866e2; co->psi3_y = 1.7803e2; co->psi3_z = -1.5972e2;
  co->psi4_x = -3.9031e3; co->psi4_y = -7.7493e2; co->psi4_z = 8.8195e2;
  co->psi6_x = -2.4874e4; co->psi6_y = -1.4892e3; co->psi6_z = 4.4588e3;
  co->psi7_x = 2.5196e4; co->psi7_y = 3.3970e2; co->psi7_z = -3.9573e3;
}


/* This function contains the coeffs from the matching with the LONG */
/* Jena waveforms (those are not the ones published in the paper)         */
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


/* function to display program usage */
static void
print_usage( CHAR *program )
{
  fprintf( stderr,
      "\n Usage: %s [options]\n\n"\
      "  Options not surrounded in [] are required.\n\n"\
      "  [--help]                 display this message\n"\
      "  --mass-ratio     q       set mass ratio to q\n"\
      "  --total-mass     M       set total mass (in M_sun) of the binary "\
	   "system\n"\
      "  --low-freq (Hz)  fLow    set the lower frequency (in Hz) of the "\
	   "hybrid wave\n\t\t\t (see below for recommended values)\n"\
      "  --sample-rate (Hz)       set the sample rate of the output waveform\n" \
      "  --output-file            output file name\n\n", program);
  fprintf( stderr,
      " Recommended use: 5 < total-mass < 400 (M_sun) and 1 < mass-ratio < 4\n"\
      " Recommended values for low-freq:\n"\
        "\t\t35 Hz for total-mass <= 50 M_sun;\n"\
        "\t\t20 Hz for 50 M_sun < total-mass <= 100 M_sun;\n"\
        "\t\t15 Hz for 100 M_sun < total-mass <= 200 M_sun;\n"\
        "\t\t8 Hz for 200 M_sun < total-mass <= 400 M_sun;\n"\
        "\t\t5 Hz for total-mass > 400 M_sun.\n\n"\
      " NOTE: Output time-series is in units of seconds.\n"\
      "\n");
}
