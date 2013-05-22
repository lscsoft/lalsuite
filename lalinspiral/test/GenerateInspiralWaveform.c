/*
*  Copyright (C) 2007 Anand Sengupta, Thomas Cokelaer, Evan Ochsner
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
\author Sathyaprakash, B. S., Cokelaer T.
\file
\ingroup LALInspiral_h

\brief Test routine for wave generation codes.

To get some help just type the name of the executable and the option --h

Basically, you can provide all the arguments from the InspiralTemplate structure such as
--approximant, --order ....

*/


#define LALGENERATEINSPIRALWAVEFORMC_ENORM 0
#define LALGENERATEINSPIRALWAVEFORMC_ESUB  1
#define LALGENERATEINSPIRALWAVEFORMC_EARG  2
#define LALGENERATEINSPIRALWAVEFORMC_EVAL  3
#define LALGENERATEINSPIRALWAVEFORMC_EFILE 4
#define LALGENERATEINSPIRALWAVEFORMC_EMEM  5

#define LALGENERATEINSPIRALWAVEFORMC_MSGENORM "Normal exit"
#define LALGENERATEINSPIRALWAVEFORMC_MSGESUB  "Subroutine failed"
#define LALGENERATEINSPIRALWAVEFORMC_MSGEARG  "Error parsing arguments"
#define LALGENERATEINSPIRALWAVEFORMC_MSGEVAL "Input argument out of valid range"
#define LALGENERATEINSPIRALWAVEFORMC_MSGEFILE "Could not open file"
#define LALGENERATEINSPIRALWAVEFORMC_MSGEMEM  "Out of memory"

#include <stdio.h>
#include <lal/LALInspiral.h>
#include <lal/LALNoiseModels.h>
#include <lal/RealFFT.h>
#include <lal/AVFactories.h>
#include <lal/Random.h>
#include <lal/GenerateInspiral.h>

#define ERROR( code, msg, statement )                                \
do                                                                   \
if ( lalDebugLevel & LALERROR )                                      \
{                                                                    \
  LALPrintError( "Error[0] %d: program %s, file %s, line %d, %s\n"   \
         "        %s %s\n", (code), program, __FILE__,       \
         __LINE__, "$Id$", statement ? statement :  \
                 "", (msg) );                                        \
}                                                                    \
while (0)

#define WARNING( statement )                                         \
do                                                                   \
if ( lalDebugLevel & LALWARNING )                                    \
{                                                                    \
  LALPrintError( "Warning[0]: program %s, file %s, line %d, %s\n"    \
         "        %s\n", program, __FILE__, __LINE__,        \
         "$Id$", (statement) );                         \
}                                                                    \
while (0)

#define INFO( statement )                                            \
do                                                                   \
if ( lalDebugLevel & LALINFO )                                       \
{                                                                    \
  LALPrintError( "Info[0]: program %s, file %s, line %d, %s\n"       \
         "        %s\n", program, __FILE__, __LINE__,        \
         "$Id$", (statement) );                         \
}                                                                    \
while (0)

#define SUB( func, statusptr )                                       \
do                                                                   \
if ( (func), (statusptr)->statusCode )                               \
{                                                                    \
  ERROR( LALGENERATEINSPIRALWAVEFORMC_ESUB, LALGENERATEINSPIRALWAVEFORMC_MSGESUB,                      \
         "Function call \"" #func "\" failed:" );                    \
  exit( LALGENERATEINSPIRALWAVEFORMC_ESUB );                                          \
}                                                                    \
while (0)

typedef struct{
  INT4 PrintParameters;
  INT4 output;
  INT4 psd_data_file;
  INT4 taper;
  INT4 realImag;
  char tag[200];
  char waveformString[LIGOMETA_WAVEFORM_MAX];
} OtherParamIn;

char *program;

void printf_timeseries (FILE *f1, UINT4 n, REAL4 *signal1, REAL8 delta);
void printf_timeseries2 (UINT4 n, REAL4 *signal1, REAL4 *signal2, REAL8 delta);
void ParseParameters(UINT4 argc, CHAR **argv, OtherParamIn *otherIn);
void LALGenerateInspiralWaveformHelp(void);
void readPSD(REAL8 *psd, REAL4 Df, UINT4 length);
void buildhoft(LALStatus *status, REAL4Vector *wave, 
        InspiralTemplate *params, OtherParamIn *otherIn);


/* --- Main part --- */
int main (int argc , char **argv) {
  REAL4Vector *signal1 = NULL;   /* storing waveforms */
  REAL4Vector *signal2 = NULL;   /* storing waveforms */
  REAL4Vector *dummy = NULL;     /* stores h+ when hx desired */
  REAL8Vector *psd=NULL;
  static LALStatus status;
  InspiralTemplate params;       /* the parameters */
  REAL8 dt, cosI, ampFac, plusFac, crossFac, dist;
  UINT4 n,i;
  UINT4 nby2;
  InspiralInit paramsInit;
  expnCoeffs ak;

  RealFFTPlan *revp =NULL;
  RealFFTPlan *frwd =NULL;
  char name1[256];

  static OtherParamIn otherIn; /* some extra parameters to parse */
  FILE *f1=NULL, *f2=NULL, *f3=NULL, *f4=NULL; /* output file pointers */

  program = *argv;

  /* ---  we start real computation here --- */
  otherIn.PrintParameters = 0; /* by default we don't print the parameters */
  otherIn.output = 0;          /* by default we don't output the waveforms */
  otherIn.realImag = 0;        /* by default output FD waveforms as |h(f)| */
  strncpy(otherIn.tag, "1", sizeof(otherIn.tag));/*default tag for file names*/
  ParseParameters(argc, argv, &otherIn);/* let's parse user parameters */
  SUB( LALInspiralITStructureSetDefault(&status, &params),
       &status);

  SUB( LALInspiralITStructureParseParameters(&status, argc, argv, &params),
       &status); /* parse user inspiral template parameters */


  SUB( LALInspiralInit(&status, &params, &paramsInit), &status);

  /**************************************************************
   * The following are used only when the psd is read from a
   * file as, for example, the s5 psd
   *************************************************************/

  if (otherIn.psd_data_file)
  {
    n = params.tSampling * 8;
    if (paramsInit.nbins < n)
    {
      paramsInit.nbins = n;
    }
    else
    {
      fprintf(stderr, "Length is not enough\n");
      exit(0);
    }
  }


  if(otherIn.output)
  {
    sprintf(name1, "wave-TD-%s.dat", otherIn.tag); f1 = fopen(name1, "w");
    sprintf(name1, "wave-OT-%s.dat", otherIn.tag); f2 = fopen(name1, "w");
    sprintf(name1, "wave-NW-%s.dat", otherIn.tag); f3 = fopen(name1, "w");
    sprintf(name1, "wave-FD-%s.dat", otherIn.tag); f4 = fopen(name1, "w");
  }



  if (otherIn.PrintParameters)
  {
    fprintf(stderr, "the inspiral structure (your parameters) before the call to the waveform generation:\n");
    SUB( LALInspiralITStructurePrint(&status, params),  &status);
  }

  /* force those parameters */
  dt     = 1./params.tSampling;
  n      = paramsInit.nbins;
  nby2   = n/2;

  if( params.approximant  == AmpCorPPN )
  {
    n *= pow(1./(1.+params.ampOrder/2.), -8./3.);
  }

  if (otherIn.PrintParameters)
  {
    fprintf(stderr, "#Testing Inspiral Signal Generation Codes:\n");
    fprintf(stderr, "#Signal n=%d, t0=%e, t2=%e, t3=%e\n", 
        n, params.t0, params.t2, params.t3);
    fprintf(stderr,"#size in bins %d\n",n);
    fprintf(stderr,"#size in seconds %f\n",params.tC);
  }

  SUB( LALSCreateVector(&status, &(signal1), n), &status);
  SUB( LALSCreateVector(&status, &(signal2), n), &status);
  SUB( LALSCreateVector(&status, &(dummy), n), &status);
  SUB( LALCreateForwardRealFFTPlan(&status, &frwd, n, 0), &status);
  SUB( LALCreateReverseRealFFTPlan(&status, &revp, n, 0), &status);
  LALDCreateVector(&status, &psd, (nby2+1));

  params.ieta = 1;

  LALInspiralSetup (&status, &ak, &params);


  switch (params.approximant)
  {
    /* These FD approximants generate h(f) and IFFT to get h(t) */
    case PadeF1:
    case TaylorF1:
    case TaylorF2:
    case TaylorF2RedSpin:
    case FindChirpSP:
    case BCV:
    case BCVC:
    case BCVSpin:
    case IMRPhenomFA:
    case IMRPhenomFB:
      switch( otherIn.output )
      {
        default:
        case 0:
        case 1:
          cosI = cos(params.inclination);
          params.distance *= LAL_PC_SI*1e6;
          ampFac  = -2.0 * params.mu * LAL_MRSUN_SI/params.distance;
          plusFac = ampFac * (1.0 + cosI*cosI);
          params.signalAmplitude = plusFac;
          SUB(LALInspiralWave(&status, signal2, &params), &status);
          SUB(LALREAL4VectorFFT(&status, signal1, signal2, revp), &status);
          break;
        case 2:
          cosI = cos(params.inclination);
          params.distance *= LAL_PC_SI*1e6;
          ampFac  = -2.0 * params.mu * LAL_MRSUN_SI/params.distance;
          crossFac = ampFac * 2.0 * cosI;
          params.signalAmplitude = crossFac;
          SUB(LALInspiralWaveTemplates(&status, dummy, signal2, &params), 
              &status);
          SUB(LALREAL4VectorFFT(&status, signal1, signal2, revp), &status);
          break;
        case 3:
          buildhoft(&status, signal2, &params, &otherIn);
          SUB(LALREAL4VectorFFT(&status, signal1, signal2, revp), &status);
          break;
      }
      break;
    /* For default case, print a warning, but continue */
    /* If a new waveform was implemented, it may succeed */
    /* Otherwise, the call to LALInspiral, etc. will still fail */
    default:
      fprintf(stderr, "Warning! approximant appears to be unknown\n");
    /* These TD approximants generate h(t), h+(t), hx(t) directly */
    case TaylorT1:
    case TaylorT2:
    case TaylorT3:
    case TaylorT4:
    case TaylorEt:
    case TaylorN:
    case EOB:
    case EOBNR:
    case EOBNRv2:
    case EOBNRv2HM:
    case PadeT1:
    case GeneratePPN:
    case AmpCorPPN:
    case SpinTaylorT3:
    case SpinTaylor:
    case SpinTaylorFrameless:
    case SpinQuadTaylor:
    case PhenSpinTaylorRD:
    case IMRPhenomA:
    case IMRPhenomB:
    case Eccentricity:
      switch( otherIn.output )
      {
        default:
        case 0:
        case 1:
          cosI = cos(params.inclination);
          params.distance *= LAL_PC_SI*1e6;
          ampFac  = -2.0 * params.mu * LAL_MRSUN_SI/params.distance;
          plusFac = ampFac * (1.0 + cosI*cosI);
          params.signalAmplitude = plusFac;
          SUB(LALInspiralWave(&status, signal1, &params), &status);
          break;
        case 2:
          cosI = cos(params.inclination);
          params.distance *= LAL_PC_SI*1e6;
          ampFac  = -2.0 * params.mu * LAL_MRSUN_SI/params.distance;
          crossFac = ampFac * 2.0 * cosI;
          params.signalAmplitude = crossFac;
          SUB(LALInspiralWaveTemplates(&status, dummy, signal1, &params), 
              &status);
          break;
        case 3:
          buildhoft(&status, signal1, &params, &otherIn);
          break;
      }
      if(otherIn.taper) /* Taper if requested */
      {
        LALSimInspiralApplyTaper bookends;
        bookends = 0;
        if (otherIn.taper==1) bookends = LAL_SIM_INSPIRAL_TAPER_START;
        if (otherIn.taper==2) bookends = LAL_SIM_INSPIRAL_TAPER_END;
        if (otherIn.taper==3) bookends = LAL_SIM_INSPIRAL_TAPER_STARTEND;
        XLALSimInspiralREAL4WaveTaper(signal1, bookends);
      }
      break;
  }
  if(otherIn.output) printf_timeseries (f1, n, signal1->data, dt);
  {
    REAL8 df, f, hSq, sSq, hMag, sMag, sRe, sIm, rho, rhosq=0, rhoDet=8.;

    SUB( LALREAL4VectorFFT(&status, signal2, signal1, frwd), &status);
    df = 1./(n * dt);
    if (otherIn.psd_data_file)
    {
      readPSD(psd->data, df, nby2);
    }
    else
    {
      LALNoiseSpectralDensity(&status, psd, &LALLIGOIPsd, df);
    }
    for (i=1; i<nby2; i++)
    {
      UINT4 j = n-i;
      f = (double)i*df;
      if (psd->data[i])
      {
        sRe  = signal2->data[i] * dt;
        sIm  = signal2->data[j] * dt;
        hSq  = sRe*sRe + sIm*sIm;
        hMag = sqrt(hSq);
        sSq  = hSq / (psd->data[i]) ;
        sMag = hMag / (psd->data[i]) ;
        if (f>params.fLower)
        {
          if (params.approximant != EOBNR && params.approximant != EOBNRv2
               && params.approximant != EOBNRv2HM && params.fFinal > f)
            rhosq += sSq;
          else
            rhosq += sSq;
        }
        signal2->data[i] = sRe / sqrt(psd->data[i]);
        signal2->data[j] = sIm / sqrt(psd->data[i]);
        if( otherIn.realImag == 1 )
        {
          if(otherIn.output) fprintf(f3, "%e %e %e %e\n", f, 
            sRe/(psd->data[i]), sIm/(psd->data[i]), psd->data[i]);
          if(otherIn.output) fprintf(f4, "%e %e %e\n", f, sRe, sIm);
        }
        else
        {
          if(otherIn.output) fprintf(f3, "%e %e %e\n", 
            f, sMag, psd->data[i]);
          if(otherIn.output) fprintf(f4, "%e %e\n", f, hMag);
        }
      }
      else
      {
        signal2->data[i] = 0.;
        signal2->data[j] = 0.;
        sSq = 0.;
      }
    }
    /* Above, 'rhosq' = \Sum_i | h(f_i) |^2 / Sn(f_i)
     * Therefore, SNR^2 = 4 \int | h(f) |^2 / Sn(f) df ~= 4 * rhosq * df
     * so SNR = rho = sqrt(4 * rhosq * df)
     */
    rho = sqrt(4. * rhosq * df);
    /* Distance in Mpc at which the SNR is 8 */
    dist = (params.distance/(LAL_PC_SI*1e6))*(rho/rhoDet);
    if( otherIn.PrintParameters )
    {
      fprintf(stderr, "mass1: %e, mass2: %e, # samples: %d,\nSNR: %e, " 
          "distance at which SNR=8: %e\n", params.mass1, params.mass2, 
          signal2->length, rho, dist);
    }
    signal2->data[0] = 0.;
    signal2->data[nby2] = 0.;
    SUB( LALREAL4VectorFFT(&status, signal1, signal2, revp), &status);
    if(otherIn.output) printf_timeseries (f2, n, signal1->data, dt);
  }

  SUB( LALDestroyRealFFTPlan (&status, &revp), &status);
  SUB( LALDestroyRealFFTPlan (&status, &frwd), &status);
  SUB( LALSDestroyVector(&status, &signal1), &status);
  SUB( LALSDestroyVector(&status, &signal2), &status);
  SUB( LALDDestroyVector(&status, &psd), &status);

  if (otherIn.PrintParameters)
  {
    fprintf(stderr, "fFinal = %f Hz tFinal = %f seconds\n" , 
        params.fFinal, params.tC);
    fprintf(stderr, "the inspiral structure after the call "
        "to the waveform generation:\n");
    SUB( LALInspiralITStructurePrint(&status, params),  &status);
  }
  /*
  LALCheckMemoryLeaks();
  */
  if(otherIn.output)
  {
    fclose(f1);
    fclose(f2);
    fclose(f3);
    fclose(f4);
  }

  return 0;
}

void
ParseParameters( UINT4              argc,
                 CHAR             **argv,
                 OtherParamIn     *otherIn)
{
  UINT4 i = 1, i1 = 0;
  INT4 order = 0;

  while(i < argc)
    {
      if ( strcmp(argv[i],    "--verbose")     == 0 )
      {
        otherIn->PrintParameters = 1;
      }
      else if ( strcmp(argv[i],    "--output")     == 0 )
      {
        otherIn->output = atoi(argv[++i]);
      }
      else if ( strcmp(argv[i],    "--taper")     == 0 )
      {
        otherIn->taper = atoi(argv[++i]);
      }
      else if ( strcmp(argv[i],    "--psd-data-file") == 0 )
      {
        otherIn->psd_data_file = 1;
      }
      else if( strcmp(argv[i],    "--h")     == 0 )
      {
        LALGenerateInspiralWaveformHelp();
      }
      else if( strcmp(argv[i],"-h")     == 0 )
      {
        LALGenerateInspiralWaveformHelp();
      }
      else if( strcmp(argv[i],"-help")     == 0 )
      {
        LALGenerateInspiralWaveformHelp();
      }
      else if( strcmp(argv[i],"--help")    == 0 )
      {
        LALGenerateInspiralWaveformHelp();
      }
      else if( strcmp(argv[i],"--tag") == 0 )
      {
        strncpy(otherIn->tag, argv[++i], sizeof(otherIn->tag));
      }
      else if( strcmp(argv[i],"--approximant") == 0 )
      {
        i1 = i + 1;
      }
      else if( strcmp(argv[i],"--order") == 0 )
      {
        order = atoi(argv[++i]);
      }
      else if( strcmp(argv[i],"--real-imag") == 0 )
      {
        otherIn->realImag = atoi(argv[++i]);
      }
      i++;
    }
    if( otherIn->output == 3 )
    { /* concatenate approximant and order in a string for SimInspiralTable */
      strcpy( otherIn->waveformString, argv[i1] );
      switch( order )
      {
        case 8:
          strcat( otherIn->waveformString, "pseudoFourPN" );
          break;
        case 7:
          strcat( otherIn->waveformString, "threePointFivePN" );
          break;
        case 6:
          strcat( otherIn->waveformString, "threePN" );
          break;
        case 5:
          strcat( otherIn->waveformString, "twoPointFivePN" );
          break;
        case 4:
          strcat( otherIn->waveformString, "twoPN" );
          break;
        case 3:
          strcat( otherIn->waveformString, "onePointFivePN" );
          break;
        case 2:
          strcat( otherIn->waveformString, "onePN" );
          break;
        case 1:
          strcat( otherIn->waveformString, "oneHalfPN" );
          break;
        case 0:
	  printf(" WARNING: you have chose Newtonian order\n");
          strcat( otherIn->waveformString, "newtonian" );
          break;
        default:
          printf("Invalid PN order requested, using 3.5PN\n");
          strcat( otherIn->waveformString, "threePointFivePN" );
          break;
      }
    }
}


void LALGenerateInspiralWaveformHelp(void)
{

  fprintf(stderr,"GenerateInspiralWaveform Help\n");
  fprintf(stderr,"---------------------------------------------------------\n");
  fprintf(stderr,"All unspecified parameters use default values.\n");
  fprintf(stderr,"If you don't know what a parameter does it's");
  fprintf(stderr," (probably) safe to omit it.\n");
  fprintf(stderr,"---------------------------------------------------------\n");
  fprintf(stderr,"--h for help\n");
  fprintf(stderr,"--verbose to print Inspiral Template parameters\n");
  fprintf(stderr,"--tag TAG adds 'TAG' identifier to file names\n");
  fprintf(stderr,"--psd-data-file file.dat reads noise curve from file.dat (NOT WORKING YET!)\n");
  fprintf(stderr,"\t Initial LIGO design PSD used if no PSD file is given\n");
  fprintf(stderr,"--output=N to generate time and freq. domain waveforms:\n");
  fprintf(stderr,"         N = 0 - do not output any files (default)\n");
  fprintf(stderr,"         N = 1 - output plus polarization h+\n");
  fprintf(stderr,"         N = 2 - output cross polarization hx\n");
  fprintf(stderr,"         N = 3 - output measured strain h = F+ h+ + Fx hx\n");
  fprintf(stderr,"For N != 0, the waveform will be written to these files:\n");
  fprintf(stderr,"         wave-TD-TAG.dat (time domain)\n");
  fprintf(stderr,"         wave-OT-TAG.dat (TD (Wiener) optimal template)\n");
  fprintf(stderr,"         wave-FD-TAG.dat (frequency domain)\n");
  fprintf(stderr,"         wave-NW-TAG.dat (noise-weighted FD)\n");
  fprintf(stderr,"Note: Not all approximants support all types of output\n");
  fprintf(stderr,"         N = 1 calls LALInspiralWave()\n");
  fprintf(stderr,"         N = 2 calls LALInspiralWaveTemplates()\n");
  fprintf(stderr,"         N = 3 calls LALInspiralWaveForInjection() via LALGenerateInspiral()\n");
  fprintf(stderr,"If the approximant you want fails, make it callable from these functions\n");
  fprintf(stderr,"--real-imag=N controls output of complex FD and NW waveforms:\n");
  fprintf(stderr,"         N = 0 - output | h(f) | (default)\n");
  fprintf(stderr,"         N = 1 - output Re(h(f)) and Im(h(f))\n");
  fprintf(stderr,"---------------------------------------------------------\n");
  LALInspiralITStructureHelp();

}


void printf_timeseries (FILE *f1, UINT4 n, REAL4 *signal1, REAL8 delta)
{
  UINT4 i=0;

  /*
  FILE *outfile1;
  outfile1=fopen("wave1.dat","w");
  */
  do
     /* fprintf (outfile1,"%e %e\n", i*delta+t0, *(signal1+i)  );
      */
     fprintf (f1,"%e %e\n", i*delta, signal1[i]  );
  while (n-++i);

  /*
  fprintf(outfile1,"&\n");
  fclose(outfile1);
  */
}



void printf_timeseries2 (UINT4 n, REAL4 *signal1, REAL4 *signal2, REAL8 delta)
{
  UINT4 i=0;
  FILE *outfile1;

  outfile1=fopen("wave2.dat","w");

  do
     /* fprintf (outfile1,"%e %e\n", i*delta+t0, *(signal+i)  );
      */
     fprintf (stdout,"%e %e %e\n", i*delta, *(signal1+i), *(signal2+i)  );
  while (n-++i);

  fclose(outfile1);
}


void readPSD(REAL8 *psd, REAL4 Df, UINT4 length)
{
  FILE *f1;
  UINT4 i=1, j=0;
  REAL4 f;
  REAL4 x;

  psd[0] = 0;
  f1 = fopen("lho4k_070318_strain.txt", "r");

  if (fscanf(f1, "%e %e", &f, &x) == EOF)
  {
    fprintf(stderr, "Problem reading file lho4k_070318_strain.txt stopping\n");
    exit(0);
  }
  if (f!=Df)
  {
    fprintf(stderr, "Incompatible frequency resolution, exiting\n");
    exit(0);
  }

  psd[i] = x*x;

  while(fscanf(f1, "%e %e", &f, &x) != EOF && i<length)
  {
    ++i;
    psd[i] = x*x;
  }

  for (j=i; j<length; j++)
  {
    psd[j] = 0;
  }
  fclose(f1);

  return;
}

void buildhoft(LALStatus *status, REAL4Vector *wave, 
        InspiralTemplate *params, OtherParamIn *otherIn)
{
  REAL4 Fp, Fc, hp, hc, A1, A2, cosshift, sinshift;
  REAL4 cosphi, sinphi, theta, phi, psi;
  UINT4 i, len;
  CoherentGW waveform;
  memset( &waveform, 0, sizeof(CoherentGW) );
  PPNParamStruc ppnParams;
  memset( &ppnParams, 0, sizeof(PPNParamStruc) );
  ppnParams.deltaT   = 1./params->tSampling;
  ppnParams.lengthIn = 0;
  ppnParams.ppn      = NULL;

  /* Construct a SimInspiralTable... */
  SimInspiralTable simTable;
  memset( &simTable, 0, sizeof(SimInspiralTable) );
  /* ... and populate it with desired parameter values */
  /*simTable.waveform = otherIn->waveformString;*/
  memcpy( simTable.waveform, otherIn->waveformString, 
          sizeof(otherIn->waveformString) );

  if (strstr(simTable.waveform,"PhenSpinTaylorRD")) {
    if (params->axisChoice==LAL_SIM_INSPIRAL_FRAME_AXIS_ORBITAL_L) {
      strcat(simTable.waveform,"OrbitalL");}
    else if (params->axisChoice==LAL_SIM_INSPIRAL_FRAME_AXIS_TOTAL_J) {
      strcat(simTable.waveform,"TotalJ");
    }
    if (params->inspiralOnly==1) {
      strcat(simTable.waveform,"inspiralOnly");
    }
    if (params->fixedStep==1) {
      strcat(simTable.waveform,"fixedStep");
    }
  }
  switch (params->interaction) {
  case LAL_SIM_INSPIRAL_INTERACTION_NONE:
    strcat(simTable.waveform,"NO");
    break;
  case LAL_SIM_INSPIRAL_INTERACTION_SPIN_ORBIT_15PN:
    strcat(simTable.waveform,"SO15PN");
    break;
  case LAL_SIM_INSPIRAL_INTERACTION_SPIN_ORBIT_25PN:
    strcat(simTable.waveform,"SO25PN");
    break;
  case LAL_SIM_INSPIRAL_INTERACTION_SPIN_ORBIT_3PN:
    strcat(simTable.waveform,"SO");
    break;
  case LAL_SIM_INSPIRAL_INTERACTION_SPIN_SPIN_2PN:
    strcat(simTable.waveform,"SS");
    break;
  case LAL_SIM_INSPIRAL_INTERACTION_SPIN_SPIN_SELF_2PN:
    strcat(simTable.waveform,"SELF");
    break;
  case LAL_SIM_INSPIRAL_INTERACTION_QUAD_MONO_2PN:
    strcat(simTable.waveform,"QM");
    break;
  case LAL_SIM_INSPIRAL_INTERACTION_ALL_SPIN:
    strcat(simTable.waveform,"ALL_SPIN");
    break;
  case LAL_SIM_INSPIRAL_INTERACTION_TIDAL_5PN:
	strcat(simTable.waveform,"TIDAL5PN");
	break;
  case LAL_SIM_INSPIRAL_INTERACTION_TIDAL_6PN:
	strcat(simTable.waveform,"TIDAL");
	break;
  case LAL_SIM_INSPIRAL_INTERACTION_ALL:
	strcat(simTable.waveform,"ALL");
	break;
		  
  }

  simTable.mass1 = params->mass1;
  simTable.mass2 = params->mass2;
  simTable.eta = params->eta;
  simTable.mchirp = params->chirpMass;
  simTable.distance = params->distance;
  simTable.longitude = params->sourceTheta;
  simTable.latitude = params->sourcePhi;
  simTable.inclination = params->inclination;
  simTable.coa_phase = params->startPhase; /* Is this what I want?? */
  simTable.polarization = params->polarisationAngle;
  simTable.psi0 = params->psi0;
  simTable.psi3 = params->psi3;
  simTable.alpha = params->alpha;
  simTable.alpha1 = params->alpha1;
  simTable.alpha2 = params->alpha2;
  simTable.alpha3 = params->alpha3;
  simTable.alpha4 = params->alpha4;
  simTable.alpha5 = params->alpha5;
  simTable.alpha6 = params->alpha6;
  simTable.beta = params->beta;
  simTable.spin1x = params->spin1[0];
  simTable.spin1y = params->spin1[1];
  simTable.spin1z = params->spin1[2];
  simTable.spin2x = params->spin2[0];
  simTable.spin2y = params->spin2[1];
  simTable.spin2z = params->spin2[2];
  simTable.theta0 = params->orbitTheta0;
  simTable.phi0 = params->orbitPhi0;
  simTable.f_lower = params->fLower;
  simTable.f_final = params->fFinal;
  simTable.amp_order = params->ampOrder;
  sprintf(simTable.taper, "TAPER_NONE");
  /* We'll taper (or not) later in code, so just set TAPER_NONE here */
    
  theta = params->sourceTheta;
  phi   = params->sourcePhi;
  psi   = params->polarisationAngle;
  Fp    = 0.5 * (1. + cos(theta)*cos(theta) ) * cos(2.*phi) * cos(2.*psi)
              - cos(theta) * sin(2.*phi) * sin(2.*psi);
  Fc    = 0.5 * (1. + cos(theta)*cos(theta) ) * cos(2.*phi) * sin(2.*psi)
              + cos(theta) * sin(2.*phi) * cos(2.*psi);

  /* This function fills a CoherentGW with info to construct h(t) */
  LALGenerateInspiral(status, &waveform, &simTable, &ppnParams);

    
  if( waveform.h == NULL ) /* build h(t) from a, phi, shift */
  {
    len = waveform.phi->data->length;
    /* Some approximants do not set waveform.shift. We use shift(t) if set. */
    /* Otherwise, we assume the shift is a constant 0 */
    if( waveform.shift == NULL )
    {
      for(i = 0; i < len; i++)
      {
        A1 = waveform.a->data->data[2*i];
        A2 = waveform.a->data->data[2*i+1];
        cosshift = 1.;
        sinshift = 0.;
        cosphi = cos( waveform.phi->data->data[i] );
        sinphi = sin( waveform.phi->data->data[i] );
        hp = A1 * cosshift * cosphi - A2 * sinshift * sinphi;
        hc = A1 * sinshift * cosphi + A2 * cosshift * sinphi;
        wave->data[i] = Fp * hp + Fc * hc;
      }
    }
    else
    {
      for(i = 0; i < len; i++)
      {
        A1 = waveform.a->data->data[2*i];
        A2 = waveform.a->data->data[2*i+1];
        cosshift = cos( waveform.shift->data->data[i] );
        sinshift = sin( waveform.shift->data->data[i] );
        cosphi = cos( waveform.phi->data->data[i] );
        sinphi = sin( waveform.phi->data->data[i] );
        hp = A1 * cosshift * cosphi - A2 * sinshift * sinphi;
        hc = A1 * sinshift * cosphi + A2 * cosshift * sinphi;
        wave->data[i] = Fp * hp + Fc * hc;
      }
    }
  }
  else /* build h(t) from h+ and hx in waveform->h */
  {
    len=waveform.h->data->length < wave->length ? waveform.h->data->length : wave->length;
    for(i = 0; i < len; i++)
      {
	wave->data[i] = Fp * waveform.h->data->data[2*i] 
	  + Fc * waveform.h->data->data[2*i+1];
      }
    for (i=len;i<wave->length;i++) wave->data[i]=0.;
  }

  return;
}
