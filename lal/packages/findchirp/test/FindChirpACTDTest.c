/********************** <lalVerbatim file="COMPLETEFindChirpACTDTestCV">
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Program \texttt{FindChirpACTDlTest.c}}
\label{ss:FindChirpACTDTest.c}

Provides the necessary function to test the AmpCorPPN filter.


****************************************** </lalLaTeX><lalErrTable> */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>


#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/LALInspiral.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpSP.h>
#include <lal/FindChirpTD.h>
#include <lal/FindChirpACTD.h>
#include "FindChirpTDTest.h"
#include <lal/GeneratePPNInspiral.h>

#include <lal/LALDetectors.h>
#include <lal/DetResponse.h>
#include <lal/TimeDelay.h>
#include <lal/TimeFreqFFT.h>
#include <lal/SeqFactories.h>
#include <lal/RealFFT.h>


#include <lal/LALRCSID.h>
NRCSID (FINDCHIRPAMPCORTESTC,"$Id$");


static void print_usage( void );

#define TEST_STATUS( ps ) \
  ( ( ps )->statusCode && ( exit( 1 ), \
    fprintf( stderr, "LAL Routine failed near line %d\n", __LINE__ ), 1 ) )




/* Default parameter settings. */
#define M1       (15.0)
#define M2       (10.0)
#define DIST     (1000)

#define INC      (45.0)
#define PHIC     (1.8829478)
#define PHI      (3.8928495)
#define PSI      (1.8473762)
#define THETA    (2.6254967)


#define FMIN     (40.0)
#define FMAX     (1000.0)
#define SRATE    (2048)
#define ORDER    (4)
#define AMP      (1)

static void print_usage()
{
  fprintf( stderr, " --help                     : Print this message! \n");
  fprintf( stderr, " --overlap                  : Normalises input data \n");
  fprintf( stderr, " --signal                   : Replaces data with signal\n");
  fprintf( stderr, " --dynrange DYNRANGE        : set the dynamic range \n");
  fprintf( stderr, " --flatpsd                  : Use flat psd \n");
  fprintf( stderr, " --dominant                 : Inject only the domintnat harmonic \n");
  fprintf( stderr, " --h-plus                   : inject only h+\n");
  fprintf( stderr, " --enable-output            : Print output files \n");
  fprintf( stderr, " --tmplt-masses MASS1 MASS2 : Specify template masses \n");
  fprintf( stderr, " --sgnl-masses MASS1 MASS2  : Specify signal masses \n");
  fprintf( stderr, " --iota IOTA                : Specify Inclination \n");
  fprintf( stderr, " --phiC PHIC                : Specify coalescence phase \n");
  fprintf( stderr, " --phi PHI                  : Specify sky angle phi \n");
  fprintf( stderr, " --theta THETA              : Specify sky angle theta \n");
  fprintf( stderr, " --psi PSI                  : Specify polarisation psi \n");
  fprintf( stderr, " --dist DIST                : Specify signal distance \n");
  fprintf( stderr, " --amp-order AMP            : Specify signal amplitude order \n");
  fprintf( stderr, " --phase-order ORDER        : Specify signal phase order \n\n\n");
  return;
}

LALStatus status;
int lalDebugLevel = 1;

int main( int argc, char **argv )
{
  /*  command line parsing variables */
  int arg;
  REAL4 mass1 = M1;
  REAL4 mass2 = M2;

  REAL4 sigMass1 = M1;
  REAL4 sigMass2 = M2;

  UINT4 injSignal = 0;

  REAL4 inc   = INC   / 360 * LAL_TWOPI;
  REAL4 phiC  = PHIC  / 360.0 * LAL_TWOPI;
  REAL4 phi   = PHI   / 360.0 * LAL_TWOPI;
  REAL4 theta = THETA / 360.0 * LAL_TWOPI;
  REAL4 psi   = PSI   / 360.0 * LAL_TWOPI;

  REAL4 dist = DIST;
  REAL4 dt   = 1./SRATE;

  INT4 amp = AMP;
  INT4 order = ORDER;


  const UINT4 numSegments  = 1;
  UINT4 numPoints          = 262144;
  const UINT4 numChisqBins = 8;
  const UINT4 invSpecTrunc = 0;
  REAL4 srate              = SRATE;   /* Hz */
  REAL4 fmin               = FMIN;    /* Hz */
  REAL4 fmax               = FMAX;    /* Hz */
  REAL4 dynRange           = 1.0;

  int i, j;
  REAL4 ts;
  INT4 output = 0;
  INT4 overlap = 0;
  INT4 flatpsd = 0;
  INT4 dominant = 0;
  INT4 h_plus = 0;
  FILE *fp, *fpTwo;
  REAL4 segNormSum = 0;


  SnglInspiralTable *event = NULL;

  FindChirpInitParams initParams; /* need to populate this by hand */
  InspiralTemplate mytmplt;

  /* these are created by Start() */
  FindChirpFilterInput   *filterInput   = NULL;

  FindChirpFilterParams  *filterParams  = NULL;
  FindChirpSegmentVector *fcSegVec      = NULL;
  DataSegmentVector      *dataSegVec    = NULL;

  /* these are required for filtering; they are created by Init() */
  FindChirpTmpltParams *tmpltParams = NULL;
  FindChirpDataParams  *dataParams  = NULL;


  /* set initialization parameters */
  initParams.numSegments    = numSegments;
  initParams.numPoints      = numPoints;
  initParams.numChisqBins   = numChisqBins;
  initParams.approximant    = AmpCorPPN;
  initParams.createRhosqVec = 1;




  /*******************************************************************
   * ARGUMENT PARSING (arg stores the current position)              *
   *******************************************************************/

  arg = 1;
  while ( arg < argc )
  {
    /* Compute Overlap */
    if ( !strcmp( argv[arg], "--overlap" ) )
    {
      arg++;
      overlap = 1;
    }
    else if ( !strcmp( argv[arg], "--help") )
    {
      print_usage();
      return;
    }
    /* Set dynRange */
    else if ( !strcmp( argv[arg], "--dynrange" ) )
    {
      if ( argc > arg + 1 )
      {
        arg++;
        dynRange = atof( argv[arg++] );
      }
      else
      {
        print_usage();
        return;
      }
    }
    /* Use flat psd */
    else if ( !strcmp( argv[arg], "--flatpsd" ) )
    {
      arg++;
      flatpsd = 1;
    }
    /* Use dominant template */
    else if ( !strcmp( argv[arg], "--dominant" ) )
    {
      arg++;
      dominant = 1;
    }
    /* Print output files */
    else if ( !strcmp( argv[arg], "--h-plus" ) )
    {
      arg++;
      h_plus = 1;
    }
    /* Print output files */
    else if ( !strcmp( argv[arg], "--enable-output" ) )
    {
      arg++;
      output = 1;
    }
    /* Parse tmplt mass option. */
    else if ( !strcmp( argv[arg], "--tmplt-masses" ) )
    {
      if ( argc > arg + 2 )
      {
        arg++;
        mass1 = atof( argv[arg++] );
        mass2 = atof( argv[arg++] );
      }
      else
      {
        arg++;
        print_usage();
        return;
      }
    }
    /* Set signal mass */
    else if ( !strcmp( argv[arg], "--sgnl-masses" ) )
    {
      if ( argc > arg + 2 )
      {
        arg++;
        sigMass1 = atof( argv[arg++] );
        sigMass2 = atof( argv[arg++] );
      }
      else
      {
        arg++;
        print_usage();
        return;
      }
    }
    /* Replace data with injected signal */
    else if ( !strcmp( argv[arg], "--signal" ) )
    {
      arg++;
      injSignal = 1;
     }
    /* Parse iota option */
    else if ( !strcmp( argv[arg], "--iota" ) )
    {
      if (argc > arg + 1 )
      {
        arg++;
        inc = atof( argv[arg++] ) / 360.0 * LAL_TWOPI;
      }
      else
      {
        arg++;
        print_usage();
        return;
      }
    }
    /* parse phiC option */
    else if ( !strcmp( argv[arg], "--phiC" ) )
    {
      if (argc > arg + 1 )
      {
        arg++;
        phiC = atof( argv[arg++] ) / 360.0 * LAL_TWOPI;
      }
      else
      {
        arg++;
        print_usage();
        return;
      }
    }
    /* parse phi option */
    else if ( !strcmp( argv[arg], "--phi" ) )
    {
      if (argc > arg + 1 )
      {
        arg++;
        phi = atof( argv[arg++] ) / 360.0 * LAL_TWOPI;
      }
      else
      {
        arg++;
        print_usage();
        return;
      }
    }
    /* parse theta option */
    else if ( !strcmp( argv[arg], "--theta" ) )
    {
      if (argc > arg + 1 )
      {
        arg++;
        theta = atof( argv[arg++] ) / 360.0 * LAL_TWOPI;
      }
      else
      {
        arg++;
        print_usage();
        return;
      }
    }
    /* parse psi option */
    else if ( !strcmp( argv[arg], "--psi" ) )
    {
      if (argc > arg + 1 )
      {
        arg++;
        psi = atof( argv[arg++] ) / 360.0 * LAL_TWOPI;
      }
      else
      {
        arg++;
        print_usage();
        return;
      }
    }
    /* parse distance option */
    else if ( !strcmp( argv[arg], "--dist" ) )
    {
      if (argc > arg + 1 )
      {
        arg++;
        dist = atof( argv[arg++] ) ;
      }
      else
      {
        arg++;
        print_usage();
        return;
      }
    }
    /* parse amp order option */
    else if ( !strcmp( argv[arg], "--amp-order" ) )
    {
      if (argc > arg + 1 )
      {
        arg++;
        amp = atof( argv[arg++] ) ;
      }
      else
      {
        arg++;
        print_usage();
        return;
      }
    }
    /* parse phase order option */
    else if ( !strcmp( argv[arg], "--phase-order" ) )
    {
      if (argc > arg + 1 )
      {
        arg++;
        order = atof( argv[arg++] ) ;
      }
      else
      {
        arg++;
        print_usage();
        return;
      }
    }
    else
    {
      print_usage();
      return;
    }
  }



  /* create objects needed by  filters */
  fprintf( stderr, "Creating structures and variables...    " );
  Start( &dataSegVec, &filterInput, &filterParams, &fcSegVec, &initParams );
  fprintf( stderr, "      Done!\n" );


  /* set filter parameters, e.g., thresholds for events */
  filterParams->deltaT         = 1.0 / srate;
  filterParams->rhosqThresh    = 1e-6;
  filterParams->chisqThresh    = 1e+6;


  /* create some fake data */
  fprintf( stderr, "Making data segment...                  " );
  MakeData( dataSegVec, mass1, mass2, srate, fmin, fmax );
  fprintf( stderr, "      Done!\n" );

  for ( j = 0; j < dataSegVec->data->spec->data->length; ++j )
  {
    if( flatpsd == 1 )
    {
      dataSegVec->data->spec->data->data[j] = 1.0/dynRange;
    }
    else
    {
     dataSegVec->data->spec->data->data[j] *= 1.0/dynRange;
    }
  }

  /* Replace Data with Signal */
  if( injSignal )
  {
    PPNParamStruc   params;
    CoherentGW    waveform;
    REAL4TimeSeries  *hoft = NULL;
    COMPLEX8Vector   *Hoff;
    RealFFTPlan   *fwdPlan = NULL;
    LALDetector        det;
    double           fplus;
    double          fcross;
    double          tdelay;
    double            gmst;
    InterferometerNumber ifoNumber = LAL_IFO_H1;
    LIGOTimeGPS     time;
    time.gpsSeconds = 841000000;
    time.gpsNanoSeconds = 000000000;

    fprintf( stderr, "Replacing data with signal...           " );

    /* Fixed parameters. */
    params.position.latitude = theta;
    params.position.longitude = phi;
    params.position.system = COORDINATESYSTEM_EQUATORIAL;
    params.lengthIn = 0;

    /* Variable parameters. */
    params.deltaT = dt;
    params.mTot = sigMass1 + sigMass2;
    params.eta = sigMass1*sigMass2/( params.mTot*params.mTot );
    params.inc = inc;
    params.phi = phi;
    params.psi = psi;
    params.d = 1.0;
    params.fStartIn = fmin;
    params.fStopIn = - 1.0 /
                     ( 6.0 * sqrt(6.0) * LAL_PI * params.mTot * LAL_MTSUN_SI );

    /* Amplitude switches */
    params.ampOrder = amp;

    /* PPN parameter. */
    params.ppn = NULL;
    LALSCreateVector( &status, &(params.ppn), order + 1 );
    params.ppn->data[0] = 1.0;
    if ( order > 0 )
      params.ppn->data[1] = 0.0;
    for ( i = 2; i <= (UINT4)( order ); i++ )
      params.ppn->data[i] = 1.0;

    memset( &waveform, 0, sizeof( CoherentGW ) );
/*
    fprintf( stderr, " params.deltaT   = %e\n", params.deltaT );
    fprintf( stderr, " params.mTot     = %e\n", params.mTot );
    fprintf( stderr, " params.eta      = %e\n", params.eta );
    fprintf( stderr, " params.d        = %e\n", params.d );
    fprintf( stderr, " params.fStartIn = %e\n", params.fStartIn );
    fprintf( stderr, " params.fStopIn  = %e\n", params.fStopIn );
    fprintf( stderr, " params.inc      = %e\n", params.inc );
    fprintf( stderr, " params.amporder = %d\n", params.ampOrder );
    for( i = 0; i < order + 1; ++i )
    {
      fprintf( stderr, " params.ppn->data[%d] = %e\n", i, params.ppn->data[i]);
    }
*/
    /* Generate Signal */
    LALGeneratePPNAmpCorInspiral( &status, &waveform, &params );

    /* Compute h(t) */
    hoft = LALCalloc( 1, sizeof( *hoft ) );
    hoft->data = XLALCreateREAL4Vector( waveform.h->data->length );

    XLALReturnDetector( &det, ifoNumber );

    gmst = XLALGreenwichMeanSiderealTime( &time );

    XLALComputeDetAMResponse( &fplus, &fcross, det.response,
                        params.position.longitude, params.position.latitude,
                                                                  psi, gmst );
    tdelay = XLALTimeDelayFromEarthCenter( det.location,
                 params.position.longitude, params.position.latitude, &time );


    if ( h_plus == 1 )
    {
      for( j = 0; j < waveform.h->data->length; ++j )
      {
        hoft->data->data[j] = waveform.h->data->data[2*j];
      }
    }
    else
    {
      for( j = 0; j < waveform.h->data->length; ++j )
      {

        hoft->data->data[j] = (  fplus * waveform.h->data->data[2*j] +
                                fcross * waveform.h->data->data[2*j+1] );
      }
      /* Taper waveform */
      XLALInspiralWaveTaper( hoft->data, INSPIRAL_TAPER_STARTEND );
    }


    /* Replace Data */
    for( j = dataSegVec->data->chan->data->length-1; j > -1; --j )
    {
      INT4 diff = dataSegVec->data->chan->data->length - hoft->data->length;
      if( j >  diff )
      {
          dataSegVec->data->chan->data->data[j] =
               hoft->data->data[j - diff];
      }
      else
      {
        dataSegVec->data->chan->data->data[j] = 0.0;
      }
    }
    if( output == 1 )
    {
      fp = fopen( "tddata.dat", "w" );
      for( j = 0; j < numPoints - 1; ++j )
      {
		    fprintf( fp, "%e %e\n", j * dataSegVec->data->chan->deltaT,
                              dataSegVec->data->chan->data->data[j] );
      }
      fclose( fp );
    }
/*
      fp = fopen( "FTdata.dat", "w" );
*/
      /* Take FFT */
/*
       Hoff = XLALCreateCOMPLEX8Vector( dataSegVec->data->spec->data->length );
      LALCreateForwardRealFFTPlan( &status, &fwdPlan,
                                    dataSegVec->data->chan->data->length, 0 );
      XLALREAL4ForwardFFT( Hoff, dataSegVec->data->chan->data,
      fwdPlan ) ;

      for( j = 0; j < Hoff->length; ++j )
      {
        fprintf( fp, "%e %e %e\n", j * dataSegVec->data->spec->deltaF,
                     Hoff->data[j].re, Hoff->data[j].im );                }
      fclose( fp );
      XLALDestroyCOMPLEX8Vector( Hoff );
      LALDestroyRealFFTPlan( &status, &fwdPlan );
*/

    /* Clear Memory */
    LALSDestroyVector( &status, &(params.ppn) );

    XLALDestroyREAL4Vector( hoft->data );
    LALFree( hoft );


    LALSDestroyVectorSequence( &status, &(waveform.h->data) );
    LALSDestroyVectorSequence( &status, &(waveform.a->data) );
    LALSDestroyVector( &status, &(waveform.f->data) );
    LALDDestroyVector( &status, &(waveform.phi->data) );
    LALFree( waveform.h );
    LALFree( waveform.a );
    LALFree( waveform.f );
    LALFree( waveform.phi );

    fprintf( stderr, "      Done!\n" );
  }

  /*
   * setup the template
   */

  mytmplt.mass1           = mass1;
  mytmplt.mass2           = mass2;
  mytmplt.totalMass       = mytmplt.mass1 + mytmplt.mass2;
  mytmplt.mu              = mytmplt.mass1 * mytmplt.mass2 / mytmplt.totalMass;
  mytmplt.eta             = mytmplt.mu / mytmplt.totalMass;
  mytmplt.ieta            = 1;
  mytmplt.massChoice      = m1Andm2;
  mytmplt.startTime       = 0.;
  mytmplt.startPhase      = 0.;
  mytmplt.tSampling       = 1.0 / fcSegVec->data->deltaT;
  mytmplt.fLower          = fcSegVec->data->fLow;
  mytmplt.fCutoff         = fmax;
  mytmplt.signalAmplitude = 1;
  mytmplt.nStartPad       = 0;
  mytmplt.nEndPad         = 0;
  mytmplt.order           = twoPN;
  mytmplt.approximant     = AmpCorPPN;
  mytmplt.massChoice      = m1Andm2;
  mytmplt.OmegaS          = 0;
  mytmplt.Theta           = 0;
  mytmplt.ampOrder        = oneHalfPN;


  /*
   * initialize specific parameters
   */

  initParams.approximant = AmpCorPPN;
  initParams.numPoints = dataSegVec->data->chan->data->length;

  Init( &tmpltParams, &dataParams, &initParams, srate, fmin, dynRange,                     invSpecTrunc );


  tmpltParams->taperTmplt = INSPIRAL_TAPER_START;

  tmpltParams->bandPassTmplt = 0;

  fprintf( stderr, "Testing ACTDTemplate...                 " );

  LALFindChirpACTDTemplate( &status, filterInput->fcTmplt, &mytmplt,
                            tmpltParams  );

  ts = - (REAL4)(numPoints) * dt;
  ts = 0;
  if( output == 1 )
  {
    fp = fopen("tdtmplt.dat", "w");
    fpTwo = fopen("FTtmplt.dat", "w");
    for(i=0; i < numPoints; i++, ts+=dt )
    {
      fprintf( fp, "%e ", ts );
      for( j=0; j < NACTDVECS; ++j)
      {
        fprintf( fp, "%.4e ",
                             tmpltParams->ACTDVecs->data[i + j * numPoints] );
      }
      fprintf( fp, "\n" );

      if( i < numPoints / 2 + 1 )
      {
        fprintf( fpTwo, "%.5e ", i * srate /( numPoints ) );
        for( j=0; j < NACTDVECS; ++j )
        {
          fprintf( fpTwo, "%.5e %.5e ",
          filterInput->fcTmplt->ACTDtilde->data[i+j*(numPoints/2+1)].re,
          filterInput->fcTmplt->ACTDtilde->data[i+j*(numPoints/2+1)].im );
        }
       fprintf( fpTwo, "\n" );
      }
    }
    fclose( fp );
    fclose( fpTwo );
  }
  fprintf( stderr, "      Done!\n" );




  LALFindChirpTDData( &status, fcSegVec, dataSegVec, dataParams );



  fprintf( stderr, "Testing ACTDNormalize...                " );
  LALFindChirpACTDNormalize( &status, filterInput->fcTmplt, tmpltParams,
                               dataParams );
  if( output == 1 )
  {
    fp = fopen("FTtmpltNorm.dat","w");

    for( i = 0; i < ( numPoints / 2 + 1 ); ++i )
    {
      fprintf( fp, "%.5e ", i * srate /( numPoints ) );
      for( j = 0; j < NACTDVECS; ++j )
      {
        fprintf( fp, "%.5e %.5e ",
          filterInput->fcTmplt->ACTDtilde->data[i+j*(numPoints/2+1)].re,
          filterInput->fcTmplt->ACTDtilde->data[i+j*(numPoints/2+1)].im );
      }
      fprintf( fp, "\n");
    }
    fclose( fp );
  }
  fprintf( stderr, "      Done!\n" );


  /* Normalise data to compute overlap */
  if ( overlap == 1 )
  {
    REAL4 invRootData;
    REAL4 norm = 0.0, normTest;
    COMPLEX8Vector normTestVector;
    COMPLEX8Vector normTestVector2;

    normTestVector.length = fcSegVec->data->data->data->length;
    normTestVector.data = fcSegVec->data->data->data->data;
    normTestVector2.length = fcSegVec->data->data->data->length;
    normTestVector2.data = filterInput->fcTmplt->ACTDtilde->data +
                              (numPoints/2+1);

    memset( fcSegVec->data->segNorm->data, 0,
      fcSegVec->data->segNorm->length * sizeof(REAL4) );

    if ( dominant == 1 )
    {
      for ( j = 1; j < fcSegVec->data->data->data->length; ++j )
      {
        fcSegVec->data->data->data->data[j].re =
               filterInput->fcTmplt->ACTDtilde->data[j + (numPoints/2+1) ].re;
        fcSegVec->data->data->data->data[j].im =
               filterInput->fcTmplt->ACTDtilde->data[j + (numPoints/2+1) ].im;
      }
    }
    fprintf( stderr, "Normalising input data for overlap..." );

    for( j = 0; j < fcSegVec->data->data->data->length - 1; ++j )
    {
      if( j * fcSegVec->data->data->deltaF >= 40. )
      {
        REAL4 power;
        if ( dataParams->wtildeVec->data[j].re == 0.0 )
          printf(" We have a zero!!\n");
        power = fcSegVec->data->data->data->data[j].re *
                fcSegVec->data->data->data->data[j].re;
        power += fcSegVec->data->data->data->data[j].im *
                fcSegVec->data->data->data->data[j].im;
        norm +=  4.0 * dt * power / dataParams->wtildeVec->data[j].re
                                                            / (REAL4)numPoints;
      }
    }

    normTest = norm;

    /*
    normTest = XLALFindChirpACTDInnerProduct( &normTestVector,
                                              &normTestVector,
                                              dataParams->wtildeVec->data,
                                              40.,
                                              fcSegVec->data->data->deltaF );
    */
    invRootData = pow( normTest, -0.5 );

    for ( j = 0;  j < fcSegVec->data->data->data->length; ++j )
    {
      fcSegVec->data->data->data->data[j].re *= invRootData;
      fcSegVec->data->data->data->data[j].im *= invRootData;
    }
    fprintf( stderr, "         Done!\n");


/*
    fprintf( stderr, "   normTest  = %1.3e\n", normTest );
    normTest = XLALFindChirpACTDInnerProduct( &normTestVector,
                                              &normTestVector,
                                              dataParams->wtildeVec->data,
                                              40.,
                                              dt, numPoints );

    fprintf( stderr, "   < data, data>  = %1.3e\n", normTest );
    normTest = XLALFindChirpACTDInnerProduct( &normTestVector,
                                              &normTestVector2,
                                              dataParams->wtildeVec->data,
                                              40.,
                                              dt, numPoints );

    fprintf( stderr, "   < H2, data >   = %1.3e\n", normTest );

    fprintf( stderr, "                                              Done!\n" ); */
  }


  if( output == 1 )
  {
    fp = fopen("FTdata.dat","w");
    for( j = 0; j < fcSegVec->data->data->data->length; ++j )
    fprintf(fp, "%e %e %e \n",
             j * fcSegVec->data->data->deltaF,
             fcSegVec->data->data->data->data[j].re,
             fcSegVec->data->data->data->data[j].im );
    fclose( fp );
  }

  filterInput->segment = fcSegVec->data;



  fprintf( stderr, "Testing ACTDFilter...                   " );
  LALFindChirpACTDFilterSegment( &status, &event, filterInput, filterParams );
  fprintf( stderr, "      Done!\n" );

  if( output == 1 )
  {
    fp = fopen("rhosqVec.dat", "w");
    numPoints = filterParams->rhosqVec->data->length;

    for( i=0; i < numPoints; ++i  )
    {
      fprintf( fp, "%e \n", filterParams->rhosqVec->data->data[i] );
    }
    fclose( fp );
  }


  /* clean up memory and exit */
  Stop( &dataSegVec, &filterInput, &filterParams, &fcSegVec, numChisqBins );
  Fini( &tmpltParams, &dataParams );

  LALCheckMemoryLeaks();
  return 0;
}



/*
 *
 * Start(), Stop()
 *
 * Start() creates and initializes various structures needed by FindChirp.
 * Stop() destroys the allocated memory.
 *
 */

int Start(
    DataSegmentVector      **dataSegVec,
    FindChirpFilterInput   **filterInput,
    FindChirpFilterParams  **filterParams,
    FindChirpSegmentVector **fcSegVec,
    FindChirpInitParams     *initParams
    )
{
  LALCreateFindChirpSegmentVector( &status, fcSegVec, initParams );
  TEST_STATUS( &status );

  LALCreateFindChirpInput( &status, filterInput, initParams );
  TEST_STATUS( &status );

  LALFindChirpFilterInit( &status, filterParams, initParams );
  TEST_STATUS( &status );

  LALFindChirpChisqVetoInit( &status, (*filterParams)->chisqParams,
      initParams->numChisqBins, initParams->numPoints );
  TEST_STATUS( &status );

  LALCreateDataSegmentVector( &status, dataSegVec, initParams );
  TEST_STATUS( &status );

  return 0;
}

int Stop(
    DataSegmentVector      **dataSegVec,
    FindChirpFilterInput   **filterInput,
    FindChirpFilterParams  **filterParams,
    FindChirpSegmentVector **fcSegVec,
    UINT4 numChisqBins
    )
{
  LALDestroyFindChirpSegmentVector( &status, fcSegVec );
  TEST_STATUS( &status );

  LALDestroyFindChirpInput( &status, filterInput );
  TEST_STATUS( &status );

  LALFindChirpChisqVetoFinalize( &status, (*filterParams)->chisqParams,
     numChisqBins );
  TEST_STATUS( &status );

  LALFindChirpFilterFinalize( &status, filterParams );
  TEST_STATUS( &status );

  LALDestroyDataSegmentVector( &status, dataSegVec );
  TEST_STATUS( &status );

  return 0;
}


/*
 *
 * Init(), Fini()
 *
 * Init() creates and initializes various structures needed for filtering.
 * Fini() destroys the allocated memory.
 *
 */

int Init(
    FindChirpTmpltParams **tmpltParams,
    FindChirpDataParams  **dataParams,
    FindChirpInitParams     *initParams,
    REAL4 srate,
    REAL4 fmin,
    REAL4 dynRange,
    UINT4 trunc
    )
{
  LALFindChirpTemplateInit( &status, tmpltParams, initParams );
  TEST_STATUS( &status );

  (*tmpltParams)->deltaT   = 1.0 / srate;
  (*tmpltParams)->fLow     = fmin;
  (*tmpltParams)->dynRange = dynRange;

  LALFindChirpDataInit( &status, dataParams, initParams );
  TEST_STATUS( &status );

  (*dataParams)->fLow         = fmin;
  (*dataParams)->dynRange     = dynRange;
  (*dataParams)->invSpecTrunc = trunc;


  return 0;
}

int Fini(
  FindChirpTmpltParams **tmpltParams,
  FindChirpDataParams  **dataParams
        )
{

  LALFindChirpTemplateFinalize( &status, tmpltParams );
  TEST_STATUS( &status );


  LALFindChirpDataFinalize( &status, dataParams );
  TEST_STATUS( &status );

  return 0;
}


/*
 *
 * MakeData()
 *
 * Populate the dataSegVec with an injected inspiral and set the spectrum
 * and response vectors to unity.
 *
 */

int MakeData(
    DataSegmentVector *dataSegVec,
    REAL4 mass1,
    REAL4 mass2,
    REAL4 srate,
    REAL4 fmin,
    REAL4 fmax
    )
{
  InspiralTemplate tmplt;
  UINT4 i;
  UINT4 k;
  UINT4 n, nspec;
  FILE *fp = NULL;
  REAL8 fs, df, psdfs;

  memset( &tmplt, 0, sizeof( InspiralTemplate ) );

  tmplt.mass1           = mass1;
  tmplt.mass2           = mass2;
  tmplt.totalMass       = tmplt.mass1 + tmplt.mass2;
  tmplt.mu              = tmplt.mass1 * tmplt.mass2 / tmplt.totalMass;
  tmplt.eta             = tmplt.mu / tmplt.totalMass;
  tmplt.ieta            = 1;
  tmplt.massChoice      = m1Andm2;
  tmplt.startTime       = 0;
  tmplt.startPhase      = 0;
  tmplt.fLower          = fmin;
  tmplt.fCutoff         = fmax;
  tmplt.tSampling       = srate;
  tmplt.signalAmplitude = 1;
  tmplt.nStartPad       = dataSegVec->data->chan->data->length / 2;
  tmplt.nStartPad       = 0;
  tmplt.nEndPad         = 0;
  tmplt.approximant     = AmpCorPPN;
  tmplt.order           = twoPN;
  tmplt.massChoice      = m1Andm2;
  tmplt.OmegaS          = 0;
  tmplt.Theta           = 0;
  tmplt.ampOrder        = oneHalfPN;

  LALInspiralParameterCalc( &status, &tmplt );
  TEST_STATUS( &status );

  LALInspiralWaveLength( &status, &n, tmplt );
  TEST_STATUS( &status );
  if ( n > dataSegVec->data->chan->data->length )
  {
    fprintf( stderr, "Chirp is too long!\n" );
    exit( 1 );
  }

  n = dataSegVec->data->chan->data->length;
  nspec = dataSegVec->data->spec->data->length;


  memset( dataSegVec->data->chan->data->data, 0,
     n * sizeof( *dataSegVec->data->chan->data->data ) );
/*
 * Returns a memory error for AmpCorPPN. We will replace this vector anyway
 *
 */
/*
  LALInspiralWave( &status, dataSegVec->data->chan->data, &tmplt );
  TEST_STATUS( &status );
*/


  dataSegVec->data->chan->deltaT = 1.0 / srate;
  df = dataSegVec->data->spec->deltaF = srate / (double) n;

  dataSegVec->data->chan->epoch.gpsSeconds     = 0;
  dataSegVec->data->chan->epoch.gpsNanoSeconds = 0;

  fs = 40.;
  LALLIGOIPsd(&status, &psdfs, fs);

  for ( k = 0; k < nspec; ++k )
  {
    REAL8 f = k*df;
    REAL8 psd;
    if( f >= fs )
    {
      LALLIGOIPsd(&status, &psd, f);
      dataSegVec->data->spec->data->data[k] = psd;
    }
    else
    {
      dataSegVec->data->spec->data->data[k] = psdfs;
    }


    dataSegVec->data->resp->data->data[k].re = 1;
    dataSegVec->data->resp->data->data[k].im = 0;
  }

  return 0;
}



