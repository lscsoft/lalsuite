/*
*  Copyright (C) 2007 Duncan Brown, Patrick Brady
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
 * File Name: chirplen.c
 *
 * Author: Brown, D. A. and Brady, P. R.
 * 
 * 
 *-----------------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <getopt.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <lal/LALError.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALStatusMacros.h>
#include <lal/LALConstants.h>
#include <lal/ReadNoiseSpectrum.h>
#include <lal/SeqFactories.h>
#include <lal/Units.h>
#include <lal/Date.h>
#include <lal/Inject.h>
#include <lal/SimulateCoherentGW.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/TimeFreqFFT.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOLwXMLRead.h>
#include <lal/LIGOMetadataUtils.h>
#include <lalapps.h>

#define CVS_ID_STRING "$Id$"
#define CVS_NAME_STRING "$Name$"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define PROGRAM_NAME "chirplen"

/* Usage format string. */
#define USAGE \
"Usage: %s [options] [LIGOLW XML input files]\n\n"\
"  --help                    display this message\n"\
"  --version                 print version information and exit\n"\
"  --machine                 prints space delimeted output\n"\
"  --m1 m1                   mass of 1st binary element in Msun\n"\
"  --m2 m2                   mass of 2nd binary element in Msun\n"\
"  --flow fstart             low frequency cutoff\n"\
"\n"

#define FSTART 30.0
#define DELTAT (1.0/4096.0)

int verbose = 0;
int machine = 0;
int chisqFlag = 0;
int printWaveform = 0;

int main ( int argc, char *argv[] )
{
  static LALStatus      status;
  UINT4 i;
  REAL4 m1 = 0.0;
  REAL4 m2 = 0.0;
  float mtot,eta,mchirp,c0,c2,c3,c4,x,x2,x3,x4,x8,chirpTime,f_max;
  REAL4 inc=0.0;
  REAL4 longit=0.0;
  REAL4 latit=0.0;

  /* for PN signals */
  REAL4                     fstart=FSTART;
  REAL4                     fstop=2000.0;
  PPNParamStruc             ppnParams;     /* wave generation parameters */
  REAL4Vector              *phiPPN=NULL;
  CoherentGW                waveform;      /* amplitude and phase structure */
  INT4                      phiOrder=0;
  CHAR                     *waveFile=NULL;

  /* power spectrum and chisqr information */
  CHAR                     *specFile=NULL;
  REAL4                      numChisqBins=0.0;

  /* getopt arguments */
  struct option long_options[] =
  {
    {"verbose",                 no_argument,       &verbose,          1 },
    {"machine",                 no_argument,       &machine,          1 },
    {"flow",			required_argument, 0,		     'f'},
    {"m1",			required_argument, 0,		     'm'},
    {"m2",			required_argument, 0,		     'n'},
    {"specfile",                required_argument, 0,                'o'},
    {"chisq-bins",              required_argument, 0,                'p' },
    {"wavefile",                required_argument, 0,                'q'},
    {"help",                    no_argument,       0,                'h'}, 
    {"version",                 no_argument,       0,                'V'},
    {0, 0, 0, 0}
  };
  int c;


  /*
   * 
   * initialize things
   *
   */

  lal_errhandler = LAL_ERR_EXIT;
  setvbuf( stdout, NULL, _IONBF, 0 );

  /* parse the arguments */
  while ( 1 )
  {
    /* getopt_long stores long option here */
    int option_index = 0;
    int optarg_len = 0;

    c = getopt_long_only( argc, argv, 
        "f:m:n:o:p:hV", long_options, 
	&option_index );

    /* detect the end of the options */
    if ( c == -1 )
    {
      break;
    }
    
    switch ( c )
    {
      case 0:
        /* if this option set a flag, do nothing else now */
        if ( long_options[option_index].flag != 0 )
        {
          break;
        }
        else
        {
          fprintf( stderr, "Error parsing option %s with argument %s\n",
              long_options[option_index].name, optarg );
          exit( 1 );
        }
        break;

      case 'f':
        /* flow */
        fstart=atof(optarg);
        break;

      case 'm':
        /* mass1 */
        m1=atof(optarg);
        break;

      case 'n':
        /* mass2 */
        m2=atof(optarg);
        break;

      case 'o':
        optarg_len = strlen( optarg ) + 1;
        specFile = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( specFile, optarg, optarg_len );
        break;
 
      case 'p':
        /* number of chisq bins */
        numChisqBins=atof(optarg);
        chisqFlag = 1;
        break;

      case 'q':
        optarg_len = strlen( optarg ) + 1;
        waveFile = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( waveFile, optarg, optarg_len );
        printWaveform = 1;
        break;

      case 'h':
        /* help message */
        fprintf( stderr, USAGE , argv[0]);
        exit( 1 );
        break;

      case 'V':
        /* print version information and exit */
        fprintf( stdout, "Compute some basic properties of inspiral signals\n" 
            "Patrick Brady and Duncan Brown\n"
            "CVS Version: " CVS_ID_STRING "\n"
            "CVS Tag: " CVS_NAME_STRING "\n" );
        exit( 0 );
        break;

      case '?':
        fprintf( stderr, USAGE , argv[0]);
        exit( 1 );
        break;

      default:
        fprintf( stderr, "Error: Unknown error while parsing options\n" );
        fprintf( stderr, USAGE, argv[0] );
        exit( 1 );
    }
  }

  /* if the chisq is to be computed,  check we have everything */
  if ( chisqFlag ){
    if ( numChisqBins <= 0.0 || !(specFile) ){
      fprintf(stderr, "Computing chisq bin boundaries:\n");
      fprintf(stderr, "Need both numChisqBins and specFile\n");
      exit(1);
    }
  }
  
  /* maul input and print out some information */
  if ( m1 <= 0.0 || m2 <= 0.0 ){
    fprintf(stderr, "Mass parameters m1 and m2 must be positive\n");
    exit(1);
  }

  mtot = m1 + m2;
  eta = ( m1 * m2 ) / ( mtot * mtot );
  mchirp = pow( eta, 0.6) * (mtot);
  fstop = f_max = 1.0 / (6.0 * sqrt(6.0) * LAL_PI * mtot * LAL_MTSUN_SI);

  if (verbose){
    fprintf( stdout, "m1 = %e\tm2 = %e\tfLow = %e\n", m1, m2, fstart );
    fprintf( stdout, "eta = %0.2f\tm = %0.2f\tmchirp = %0.2f\n", 
        eta, mtot, mchirp);
    fprintf( stdout, "isco freq = %e Hz\n", f_max );
  }

  
  /***************************************************************************
   * this is independent code to compute the duration of the chirp 
   **************************************************************************/
  c0 = 5*mtot*LAL_MTSUN_SI/(256*eta);
  c2 = 743.0/252.0 + eta*11.0/3.0;
  c3 = -32.*LAL_PI/5.;
  c4 = 3058673.0/508032.0 + eta*(5429.0/504.0 + eta*617.0/72.0);
  x  = pow(LAL_PI*mtot*LAL_MTSUN_SI*fstart, 1.0/3.0);
  x2 = x*x;
  x3 = x*x2;
  x4 = x2*x2;
  x8 = x4*x4;
  chirpTime = c0*(1 + c2*x2 + c3*x3 + c4*x4)/x8;

  if (verbose){
    fprintf( stdout, "length = %e seconds\n", chirpTime );
  }

  /***************************************************************************
   * Generate a PN waveform using codes in LAL to determine further
   * information about the waveforms themselves
   **************************************************************************/
  ppnParams.mTot = mtot;
  ppnParams.eta = eta;
  ppnParams.d = 1.0e6 * LAL_PC_SI ;
  ppnParams.phi = 0.0;
  ppnParams.inc = inc;

  /* Set up other parameter structures. */
  ppnParams.epoch.gpsSeconds = ppnParams.epoch.gpsNanoSeconds = 0;
  ppnParams.position.latitude = longit;
  ppnParams.position.longitude = latit;
  ppnParams.position.system = COORDINATESYSTEM_EQUATORIAL;
  ppnParams.psi = 0.0;
  ppnParams.fStartIn = fstart;
  ppnParams.fStopIn = -fstop;
  ppnParams.lengthIn = 0;
  ppnParams.ppn = NULL;
  ppnParams.deltaT = DELTAT;
  memset( &waveform, 0, sizeof(CoherentGW) );

  /*******************************************************************
   * Generate the waveform
   *******************************************************************/
  if (phiOrder){
    ppnParams.ppn = phiPPN;
  }
  LAL_CALL( LALGeneratePPNInspiral( &status, &waveform, &ppnParams ), &status );
  LALPrintError( "%d: %s\n", ppnParams.termCode, ppnParams.termDescription );
  if ( ppnParams.dfdt > 2.0 ) {
    LALPrintError( "Waveform sampling interval is too large:\n"
        "\tmaximum df*dt = %f", ppnParams.dfdt );
  }

  /*******************************************************************
   * Print out the information that's wanted for the inspiral group
   *******************************************************************/
  if (machine) {
    
    fprintf( stdout, "%e %e %e %e %e %e\n", m1, m2,
        ppnParams.fStart, ppnParams.fStop, ppnParams.tc, 
        (float)(0.5/LAL_PI) * (waveform.phi->data->data[waveform.phi->data->length-1] 
                               - waveform.phi->data->data[0]) );
  
  } else {

    fprintf( stdout, "fStart according to Tev = %e Hz\n", ppnParams.fStart );
    fprintf( stdout, "fStop  according to Tev = %e Hz\n", ppnParams.fStop );
    fprintf( stdout, "length according to Tev = %e seconds\n", ppnParams.tc );
    fprintf( stdout, "Ncycle according to Tev = %f \n", 
        (float)(0.5/LAL_PI) * (waveform.phi->data->data[waveform.phi->data->length-1] 
                               - waveform.phi->data->data[0]));
  
  }

  /***********************************************************************
   * Compute the chisq bin boundaries if requested
   ***********************************************************************/
  if ( chisqFlag ){
    REAL4FrequencySeries spectrum;
    REAL4 df = 0.125;
    REAL4 freq = 0.0;
    REAL4 sum = 0.0;
    REAL4 chisqSum = 0.0;
    REAL4 norm = 0.0;
    UINT4 k = 0;

    spectrum.f0 = ppnParams.fStart;
    spectrum.deltaF = df;
    spectrum.data = NULL;

    LAL_CALL( LALCreateVector( &status, &(spectrum.data), 
          1 + (INT4)( (ppnParams.fStop - ppnParams.fStart)/spectrum.deltaF ) ), 
        &status );

    LAL_CALL( LALReadNoiseSpectrum( &status, &spectrum, specFile), &status);

    sum = 0.0;
    norm = (spectrum.data->data[0] * spectrum.data->data[0]);
    for ( k=0 ; k<spectrum.data->length ; k++){
      freq = spectrum.f0 + k * df;
      sum += norm * pow(freq, -7.0/3.0) / 
        (spectrum.data->data[k] * spectrum.data->data[k]);
    }

    chisqSum = 0.0;
    for ( k=0 ; k<spectrum.data->length ; k++){
      freq = spectrum.f0 + k * df;
      chisqSum += norm * pow(freq, -7.0/3.0) 
        / (spectrum.data->data[k] * spectrum.data->data[k]);
      if ( chisqSum/sum >= 1.0/numChisqBins ){
        fprintf(stdout,"%f ",freq);
        chisqSum = 0.0;
      }
    }
    if ( chisqSum > 0.0 ){
      fprintf(stdout,"%f ",freq);
    }
  }

  /**********************************************************************
   * Print out the waveform information
   **********************************************************************/
  if ( printWaveform ){
    FILE *fpout=NULL;

    fpout = fopen(waveFile,"w");

    for(i = 0; i<waveform.phi->data->length ; i++)
    {
      float tmper =  0.5/LAL_PI;
      fprintf(fpout,"%e %e %e\n", i*ppnParams.deltaT, 
          waveform.f->data->data[i] ,
          tmper * waveform.phi->data->data[i] );
    }

    fclose(fpout);
  }

  return 0;
}
