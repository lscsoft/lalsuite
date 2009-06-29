/*
*  Copyright (C) 2007 Anand Sengupta, Thomas Cokelaer
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

/*  <lalVerbatim file="LALGenerateInspiralWaveformCV">
Author: Sathyaprakash, B. S., Cokelaer T.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Test program \texttt{LALGenerateInspiralWaveform.c}}
Test routine for wave generation codes.

To get some help just type the name of the executable and the option --h

Basically, you can provide all the arguments from the InspiralTemplate structure such as
--approximant, --order ....

\vfill{\footnotesize\input{LALGenerateInspiralWaveformCV}}

</lalLaTeX> */


#define LALGENERATEINSPIRALWAVEFORMC_ENORM 0
#define LALGENERATEINSPIRALWAVEFORMC_ESUB  1
#define LALGENERATEINSPIRALWAVEFORMC_EARG  2
#define LALGENERATEINSPIRALWAVEFORMC_EVAL  3
#define LALGENERATEINSPIRALWAVEFORMC_EFILE 4
#define LALGENERATEINSPIRALWAVEFORMC_EMEM  5

#define LALGENERATEINSPIRALWAVEFORMC_MSGENORM "Normal exit"
#define LALGENERATEINSPIRALWAVEFORMC_MSGESUB  "Subroutine failed"
#define LALGENERATEINSPIRALWAVEFORMC_MSGEARG  "Error parsing arguments"
#define LALGENERATEINSPIRALWAVEFORMC_MSGEVAL  "Input argument out of valid range"
#define LALGENERATEINSPIRALWAVEFORMC_MSGEFILE "Could not open file"
#define LALGENERATEINSPIRALWAVEFORMC_MSGEMEM  "Out of memory"

#include <stdio.h>
#include <lal/LALInspiral.h>
#include <lal/LALNoiseModels.h>
#include <lal/RealFFT.h>
#include <lal/AVFactories.h>
#include <lal/Random.h>

NRCSID( LALGENERATEINSPIRALWAVEFORMC, "$Id$" );

INT4 lalDebugLevel=1;
#define ERROR( code, msg, statement )                                \
do                                                                   \
if ( lalDebugLevel & LALERROR )                                      \
{                                                                    \
  LALPrintError( "Error[0] %d: program %s, file %s, line %d, %s\n"   \
		 "        %s %s\n", (code), program, __FILE__,       \
		 __LINE__, LALGENERATEINSPIRALWAVEFORMC, statement ? statement :  \
                 "", (msg) );                                        \
}                                                                    \
while (0)

#define WARNING( statement )                                         \
do                                                                   \
if ( lalDebugLevel & LALWARNING )                                    \
{                                                                    \
  LALPrintError( "Warning[0]: program %s, file %s, line %d, %s\n"    \
		 "        %s\n", program, __FILE__, __LINE__,        \
		 LALGENERATEINSPIRALWAVEFORMC, (statement) );                         \
}                                                                    \
while (0)

#define INFO( statement )                                            \
do                                                                   \
if ( lalDebugLevel & LALINFO )                                       \
{                                                                    \
  LALPrintError( "Info[0]: program %s, file %s, line %d, %s\n"       \
		 "        %s\n", program, __FILE__, __LINE__,        \
		 LALGENERATEINSPIRALWAVEFORMC, (statement) );                         \
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
} OtherParamIn;

char *program;

void printf_timeseries (FILE *f1, UINT4 n, REAL4 *signal1, REAL8 delta) ;
void printf_timeseries2 (UINT4 n, REAL4 *signal1, REAL4 *signal2, REAL8 delta) ;
void ParseParameters(UINT4 argc, CHAR **argv, OtherParamIn *otherIn);
void LALGenerateInspiralWaveformHelp(void);
void readPSD(REAL8 *psd, REAL4 Df, UINT4 length);


/* --- Main part --- */
int main (int argc , char **argv) {
  REAL4Vector *signal1 = NULL;   /* storing waveforms */
  REAL4Vector *signal2 = NULL;   /* storing waveforms */
  REAL8Vector *psd=NULL;
  static LALStatus status;
  InspiralTemplate params;       /* the parameters */
  REAL8 dt, cosI, apFac, dist;   /* sampling interval, etc */
  UINT4 n,i;
  UINT4 nby2;
  InspiralInit paramsInit;
  expnCoeffs ak;

  RealFFTPlan *revp =NULL;
  RealFFTPlan *frwd =NULL;
  char name1[256];

  static OtherParamIn otherIn; /* some extra parameters to parse*/
  FILE *f1=NULL, *f2=NULL, *f3=NULL, *f4=NULL;

  program = *argv;

  /* ---  we start real computation here --- */
  otherIn.PrintParameters = 0; /* by default we don't print the parameters */
  otherIn.output = 0; /* by default we don't print the parameters */
  ParseParameters(argc, argv, &otherIn);/*let's parse user parameters     */
  SUB( LALInspiralITStructureSetDefault(&status, &params),
       &status); /*  */

  SUB( LALInspiralITStructureParseParameters(&status, argc, argv, &params),
       &status);/*parse user inspiral template parameters */


  SUB(  LALInspiralInit(&status, &params, &paramsInit), &status);
  /*  params.signalAmplitude =1;*/

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
    sprintf(name1, "wave-TD-%d.dat", params.approximant); f1 = fopen(name1, "w");
    sprintf(name1, "wave-OT-%d.dat", params.approximant); f2 = fopen(name1, "w");
    sprintf(name1, "wave-SS-%d.dat", params.approximant); f3 = fopen(name1, "w");
    sprintf(name1, "wave-SI-%d.dat", params.approximant); f4 = fopen(name1, "w");
  }



  if (otherIn.PrintParameters){
    fprintf(stderr, "the inspiral structure (your parameters) before the call to the waveform generation:\n");
      SUB( LALInspiralITStructurePrint(&status, params),  &status);
  }

  /* force those parameters */

  dt 	= 1./params.tSampling;
  n 	= paramsInit.nbins;
  nby2 = n/2;

  if( params.approximant  == AmpCorPPN )
  {
    n *= pow(1./(1.+params.ampOrder/2.), -8./3.);
  }

  if (otherIn.PrintParameters)
  {
    fprintf(stderr, "#Testing Inspiral Signal Generation Codes:\n");
    fprintf(stderr, "#Signal n=%d, t0=%e, t2=%e, t3=%e\n", n, params.t0, params.t2, params.t3);
    fprintf(stderr,"#size in bins %d\n",n);
    fprintf(stderr,"#size in seconds %f\n",params.tC);
  }

  /*
  params.ampOrder = 0;
   */
  SUB( LALSCreateVector(&status, &(signal1), n), &status);
  SUB( LALSCreateVector(&status, &(signal2), n), &status);
  SUB( LALCreateForwardRealFFTPlan(&status, &frwd, n, 0), &status);
  SUB( LALCreateReverseRealFFTPlan(&status, &revp, n, 0), &status);
  LALDCreateVector(&status, &psd, (nby2+1));

  params.ieta = 1;

  LALInspiralSetup (&status, &ak, &params);

  /*
  for (v=0.01; v<0.6; v+=0.01)
  {
	  double lum = paramsInit.func.flux(v, &paramsInit.ak)/(pow(v,10.)*paramsInit.ak.FTaN);
	  fprintf(stdout, "%e %e\n", v, lum);
  }
  exit(0);
  */

  switch (params.approximant)
  {
    case TaylorF1:
    case TaylorF2:
    case BCV:
    case BCVSpin:
      SUB( LALInspiralWave(&status, signal2, &params), &status);
      SUB( LALREAL4VectorFFT(&status, signal1, signal2, revp), &status);
      break;
    case SpinTaylorT3:
    case TaylorT1:
    case TaylorT2:
    case TaylorT3:
    case TaylorT4:
    case TaylorEt:
    case TaylorN:
    case EOB:
    case EOBNR:
    case PadeT1:
    case AmpCorPPN:
    case SpinTaylor:
    case Eccentricity:
      cosI = cos(params.inclination);
      /* We are expecting that the distance input to  this routine will be in Mpc */
      params.distance *= LAL_PC_SI*1e6;
      apFac  = -2.0 * params.mu * LAL_MRSUN_SI/params.distance;
      apFac *= (1.0 + cosI*cosI);
      params.signalAmplitude = apFac;
      SUB(LALInspiralWave(&status, signal1, &params), &status);
      if(otherIn.taper)
      {
	InspiralApplyTaper bookends;
	bookends = 0;
	if (otherIn.taper==1) bookends = INSPIRAL_TAPER_START;
	if (otherIn.taper==2) bookends = INSPIRAL_TAPER_END;
	if (otherIn.taper==3) bookends = INSPIRAL_TAPER_STARTEND;
        XLALInspiralWaveTaper(signal1, bookends);
      }
      break;
    case PadeF1:
    default:
      break;
  }
  if(otherIn.output) printf_timeseries (f1, n, signal1->data, dt);
  {
    REAL8 df, f, sSq, sRe, sIm, rho, rhosq=0, rhoDet=8.;

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
        double hSq;
	sRe = signal2->data[i];
	sIm = signal2->data[j];
	hSq = (sRe*sRe + sIm*sIm);
	sSq = hSq / (psd->data[i]) ;
	if (f>params.fLower)
	{
          if (params.approximant != EOBNR && params.fFinal > f)
            rhosq += sSq;
	  else
	    rhosq += sSq;
	}
	signal2->data[i] = sRe / sqrt(psd->data[i]) / (double)nby2;
	signal2->data[j] = sIm / sqrt(psd->data[i]) / (double)nby2;
	if(otherIn.output) fprintf(f3, "%e %e %e\n", f, sSq, sqrt(psd->data[i]));
	if(otherIn.output) fprintf(f4, "%e %e\n", f, hSq);
      }
      else
      {
        signal2->data[i] = 0.;
	signal2->data[j] = 0.;
	sSq = 0.;
      }
    }
    /* The normalization for rho^2 is dt^2 df. dt^2 is for two factors
     * of H(f), and df for the SNR integral.  A factor of 4 comes from
     * the definition of the scalar product.
     */
    rho = sqrt(4. * rhosq *dt*dt*df);
    /* Distance in Mpc at which the SNR is 10 */
    dist = (params.distance/(LAL_PC_SI*1e6))*(rho/rhoDet);
    fprintf(stdout, "%e %e %e %e %e %d\n", params.mass1, params.mass2,
      params.totalMass, rho, dist, signal2->length);
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

  fprintf(stderr, "fFinal = %f Hz tFinal = %f seconds\n" , params.fFinal, params.tC) ;
  if (otherIn.PrintParameters)
  {
    fprintf(stderr, "the inspiral structure after the call to the waveform generation:\n");
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
  if(otherIn.output)
  {
    unsigned int useed;
    REAL4Vector  *buff=NULL;
    static RandomParams *randomparams;

    useed = 1234567890;

    SUB( LALSCreateVector(&status, &(buff), n), &status);
    LALCreateRandomParams(&status, &randomparams, useed);
    LALNormalDeviates(&status, buff, randomparams);
    LALDestroyRandomParams(&status, &randomparams);
    sprintf(name1, "wave-Random-%d.dat", params.approximant);
    f4 = fopen(name1, "w");
    printf_timeseries (f4, n, buff->data, dt);
    SUB( LALSDestroyVector(&status, &(buff)), &status);
    fclose(f4);

  }

  return 0;
}

void
ParseParameters(	UINT4 			argc,
			CHAR 			**argv,
			OtherParamIn    	*otherIn)
{
  UINT4 		i = 1;

  while(i < argc)
    {
      if ( strcmp(argv[i],	"--verbose") 	== 0 )
      {
	otherIn->PrintParameters = 1;
      }
      else if ( strcmp(argv[i],	"--output") 	== 0 )
      {
	otherIn->output = 1;
      }
      else if ( strcmp(argv[i],	"--taper") 	== 0 )
      {
	otherIn->taper = atoi(argv[++i]);
      }
      else if ( strcmp(argv[i],	"--psd-data-file") == 0 )
      {
	otherIn->psd_data_file = 1;
      }
      else if( strcmp(argv[i],	"--h") 	== 0 )
      {
	LALGenerateInspiralWaveformHelp();
      }
      else if( strcmp(argv[i],"-h") 	== 0 )
      {
	LALGenerateInspiralWaveformHelp();
      }
      else if( strcmp(argv[i],"-help") 	== 0 )
      {
	LALGenerateInspiralWaveformHelp();
      }
      else if( strcmp(argv[i],"--help")	== 0 )
      {
	LALGenerateInspiralWaveformHelp();
      }
      i++;
    }
}


void LALGenerateInspiralWaveformHelp(void)
{

  fprintf(stderr,"LALGenerateInspiralWaveform Help\n");
  fprintf(stderr, "-----------------------------------------------\n");
  fprintf(stderr, "--h for help\n");
  fprintf(stderr, "--verbose to print Inspiral Template parameters\n");
  fprintf(stderr, "-----------------------------------------------\n");
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

