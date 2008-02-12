/*
*  Copyright (C) 2008 Lucia Santamaria, Robert Adam Mercer, Badri Krishnan
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

REAL8FrequencySeries *
XLALHybridP1Amplitude( 
		      PhenomParams *params,
		      REAL8        fLow,
       		      REAL8        df,
		      REAL8        eta,
		      REAL8        M,
		      UINT4        len  );

REAL8FrequencySeries *
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

COMPLEX16Vector *
XLALConstructComplexVector(
			   REAL8FrequencySeries  *Ampl,
			   REAL8FrequencySeries  *Phase  );


void LALPrintComplex16Vec (
			   COMPLEX16Vector *hT,
			   REAL8 dt,
			   CHAR *outFile);

void  LALPrintReal8Vec(
		       REAL8Vector *h, 
		       REAL8       dt, 
		       CHAR        *outFile);


/* Main Program */
INT4 main ( INT4 argc, CHAR *argv[] ) {

  static LALStatus status;
  int c;
  REAL8 dt = -1, totTime = -1;
  REAL8 totalMass, massRatio = -1;
  REAL8 fLow, df; 
  CHAR  outFile[20];
  REAL8 eta;
  UINT4 numPoints;

  PhenomCoeffs coeffs;
  PhenomParams params;

  REAL8FrequencySeries  *Aeff = NULL;
  REAL8FrequencySeries  *Phieff = NULL;
  COMPLEX16Vector       *uF = NULL;

  REAL8Vector      *hPlus = NULL, *hCross = NULL;
  REAL8FFTPlan     *prev = NULL;

  CHAR fA[20], fPhi[20], hpFile[20], hcFile[20], fuF[20];

  strcpy(fA,"Aeff.dat");
  strcpy(fPhi,"Phieff.dat");
  strcpy(outFile, "phenom_out.dat");
  strcpy(hpFile, "hPlus.dat");
  strcpy(hcFile, "hCross.dat");
  strcpy(fuF, "uF.dat");
    
  /* getopt arguments */
  struct option long_options[] =
  {
    {"mass-ratio",              required_argument, 0,                'q'},
    {"total-time (units M_sun)",required_argument, 0,                't'},
    {"delta-t (units M_sun)",   required_argument, 0,                'd'},
    {"help",                    no_argument,       0,                'h'},
    {"version",                 no_argument,       0,                'V'},
    {0, 0, 0, 0}
  };

  /* parse the arguments */
  while ( 1 )
  {
    /* getopt_long stores long option here */
    int option_index = 0;
    size_t optarg_len;

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
        fprintf( stdout, "%s - Compute Ajith's Phenomenological Waveforms and output them to a plain text file\n" \
            "CVS Version: %s\nCVS Tag: %s\n", PROGRAM_NAME, CVS_ID_STRING, \
            CVS_NAME_STRING );
        exit( 0 );
        break;

      case 'q':
        /* set mass ratio */
        massRatio = atof( optarg );
        break;

      case 't':
        /* set low freq */
        totTime = atof( optarg );
        break;

      case 'd':
        /* set delta t */
        dt = atof( optarg );
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


  /* Check validity of arguments */

  /* check we have freqs */
  if ( totTime < 0 )
  {
    fprintf( stderr, "ERROR: --total-time must be specified\n" );
    exit( 1 );
  }

  /* check we have mass ratio and delta t*/
  if ( massRatio < 0 )
  {
    fprintf( stderr, "ERROR: --mass-ratio must be specified\n" );
    exit( 1 );
  }
  if ( dt < 0 )
  {
    fprintf( stderr, "ERROR: --delta-t must be specified\n" );
    exit( 1 );
  }

  /* * * * * * * * */
  /* Main Program  */
  /* * * * * * * * */

  eta = 1. / pow(1. + massRatio, 2.);
  totalMass = 1.0;
  fLow = 40.0;

  /* if dt and totTime are given in M_sun */
  numPoints = (UINT4) totTime/dt;
  
  /* df is given by the number of points and sampling ratio of h(t) */
  df = 1/(numPoints * dt * LAL_MTSUN_SI);

  /* Phenomenological coefficients as in Ajith et. al */
  GetPhenomCoeffs(&coeffs);

  /* Compute phenomenologial parameters */
  ComputeParamsFromCoeffs(&params, &coeffs, eta, totalMass);

  /* Compute and print Aeff */
  Aeff = XLALHybridP1Amplitude(&params, fLow, df, eta, totalMass, 1 + numPoints/2);
  LALDPrintFrequencySeries( Aeff, fA);

  /* Compute and print Phieff */
  Phieff = XLALHybridP1Phase(&params, fLow, df, eta, totalMass, 1 + numPoints/2);
  LALDPrintFrequencySeries( Phieff, fPhi); 

  /* Construct u(f) = Aeff*e^(i*Phieff) */
  uF = XLALConstructComplexVector( Aeff, Phieff);
  LALPrintComplex16Vec(uF,df,fuF);

  /* Inverse Fourier transform */
  LALCreateReverseREAL8FFTPlan( &status, &prev, numPoints, 0);
  hPlus = XLALCreateREAL8Vector(numPoints);
  LALReverseREAL8FFT( &status, hPlus, uF, prev);

  /* Print hplus to file */
  LALPrintReal8Vec(hPlus, dt, hpFile); 

  /* Free Memory */
  XLALDestroyREAL8FrequencySeries(Aeff);
  XLALDestroyREAL8FrequencySeries(Phieff);
  XLALDestroyCOMPLEX16Vector(uF);

  XLALDestroyREAL8FFTPlan(prev);
  XLALDestroyREAL8Vector(hPlus);

  return(0);

}


/* My funtions */

void  LALPrintReal8Vec(
		       REAL8Vector *h, 
		       REAL8       dt, 
		       CHAR        *outFile)  {

    UINT4 i, n;
    FILE *fp;
    n = h->length;
    
    fp = LALFopen(outFile, "w");
    for (i=0; i < n; i++) {
      fprintf (fp, "%e\t%e\t\n", i*dt, h->data[i]);
    }
}


void LALPrintComplex16Vec (
			   COMPLEX16Vector *hT,
			   REAL8 dt,
			   CHAR *outFile) {

    UINT4 i, n;
    FILE *fp;
    COMPLEX16 num;
    n = hT->length;
    
    fp = LALFopen(outFile, "w");
    for (i=0; i < n; i++) {
      num = hT->data[i];
      fprintf (fp, "%e\t%e\t%e\t%e\t%e\n", dt*i, num.re, num.im, sqrt( pow(num.re,2.) + pow(num.im,2.)), atan(num.im/num.re));
    }
}


COMPLEX16Vector *
XLALConstructComplexVector(
			   REAL8FrequencySeries  *Ampl,
			   REAL8FrequencySeries  *Phase  )
{
  COMPLEX16Vector *ComplVec = NULL;
  COMPLEX16 num;
  REAL8 Re, Im;
  UINT4 k, n;
  
  n = Ampl->data->length;

  ComplVec = LALCalloc(1, sizeof(*ComplVec)); 
  ComplVec = XLALCreateCOMPLEX16Vector(n);

  for( k = 0 ; k < n ; k++ ) {

    Re = Ampl->data->data[k] * cos(Phase->data->data[k]);
    Im = Ampl->data->data[k] * sin(Phase->data->data[k]);
    num.re = Re;
    num.im = Im;
    ComplVec->data[k] = num;
  }

  return ComplVec;

}

REAL8FrequencySeries *
XLALHybridP1Phase( 
		  PhenomParams  *params,
		  REAL8         fLow,
		  REAL8         df,
		  REAL8         eta,
		  REAL8         M,
		  UINT4         len )
{
  /* Computes phi_eff as in arXiv:0710.2335 [gr-qc] */

  INT4 k;
  REAL8 piM;
  REAL8 f, psi0, psi2, psi3, psi4, psi6, psi7;

  REAL8FrequencySeries *Phieff = NULL;

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
  Phieff = XLALCreateREAL8FrequencySeries("", &epoch, 0, df, &lalDimensionlessUnit, len);

  f = 0.0;

  for( k = 0 ; k < len ; k++ ) {

    if (f <= fLow ) {
      Phieff->data->data[k] = 0.0;
    }
    else {

    /* QUESTION: what happens with the 2*pi*f*t_0 and psi_0 terms? */
    /* for the moment they're set to zero */
    /* (see Eq. (4.19) of the paper */

    Phieff->data->data[k] = psi0 * pow(f*piM , -5./3.) + 
                            psi2 * pow(f*piM , -3./3.) +
                            psi3 * pow(f*piM , -3./3.) +
                            psi4 * pow(f*piM , -1./3.) +
                            psi6 * pow(f*piM , 1./3.) +
                            psi7 * pow(f*piM , 2./3.);
  
    Phieff->data->data[k] /= eta;
    }
    f += df; 
  }

  return Phieff;
}


REAL8FrequencySeries *
XLALHybridP1Amplitude( 
		      PhenomParams *params,
		      REAL8        fLow,
       		      REAL8        df,
		      REAL8        eta,
		      REAL8        M,
		      UINT4        len  )
{
  INT4 k;
  REAL8 piM;
  REAL8 cConst;
  REAL8 f, fNorm, fMerg, fRing, fCut, sigma;

  REAL8FrequencySeries *Aeff = NULL;

  LIGOTimeGPS epoch;

  fMerg = params->fMerger;
  fCut = params->fCut;
  fRing = params->fRing;
  sigma = params->sigma;

  piM = LAL_PI * M * LAL_MTSUN_SI;

  /* This constant has units of time */
  /* cConst = pow(LAL_PI,-1./2.) * LAL_MTSUN_SI * mass * 
    pow( piMfMerg, -1./6. ) /( piMfMerg/piM * effDist );
    cConst *= pow(5.*eta/24., 1./2.); */
  /* For the moment let's not care about this constant */
  cConst = 1.0;

  /* Allocate memory for the frequency series */
  epoch.gpsSeconds = 0;
  epoch.gpsNanoSeconds = 0;
  Aeff = XLALCreateREAL8FrequencySeries("", &epoch, 0, df, &lalDimensionlessUnit, len);

  f = 0.0;

  for( k = 0 ; k < len ; k++ ) {

    fNorm = f / fMerg;

    if (f <= fLow ) {
      Aeff->data->data[k] = 0.0;
    }
    else if ( f > fLow && f <= fMerg ) {
      Aeff->data->data[k] = pow (fNorm, -7./6.);
    }
    else if ( f > fMerg && f <= fRing ) {
      Aeff->data->data[k] = pow (fNorm, -2./3.);
    }
    else if ( f > fRing && f <= fCut ) {
      Aeff->data->data[k] = XLALLorentzian ( f, fRing, sigma);
      Aeff->data->data[k] *= LAL_PI_2*pow(fRing/fMerg,-2./3.)*sigma;
    }
    else if (f > fCut ) {
      Aeff->data->data[k] = 0.0;
    }
    Aeff->data->data[k] *= cConst; 
    f += df; 
  }

  return Aeff;
}

void ComputeParamsFromCoeffs(
			     PhenomParams *params,
			     PhenomCoeffs *coeffs,
			     REAL8        eta,
			     REAL8        M) {

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
void GetPhenomCoeffs(
		     PhenomCoeffs *co) {

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
/* Jena waveforms (those are not the published in the paper)         */
void GetPhenomCoeffsLongJena(
		     PhenomCoeffs *co) {

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

/* function to display program usgae */
static void print_usage( CHAR *program )
{
  fprintf( stderr,
      "Usage:  %s [options]\n"\
      "The following options are recognized.\n"\
	"   Options not surrounded in [] are required.\n"\
      "  [--help]                 display this message\n"\
      "  [--version]              print version information and exit\n"\
      "  --mass-ratio         q   set mass ratio to q\n"\
      "  --total-time (M_sun) T   set total length (in M_sun) of the hybrid waveform\n"\
      "  --delta-t (M_sun)    dt  set dt (in M_sun) of the hybrid waveform\n"\
      "  The output waveform will be written to hPlus.dat\n"\
      "\n", program );
}
