/*
*  Copyright (C) 2007 Gareth Jones
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
 * File Name: FindChirpBCVTest.c
 *
 * Author: Brown, D. A., Jones, G
 *
 * Revision: $Id$
 *
 *-----------------------------------------------------------------------
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <lal/LALConfig.h>

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#ifdef HAVE_GETOPT_H
#include <getopt.h>
#endif

#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/Random.h>
#include <lal/LALInspiral.h>
#include <lal/FindChirp.h>
#include <lal/FindChirp.h>
#include <lal/TwoInterfFindChirp.h>

NRCSID (MAIN, "$Id$");

static void
TestStatus (LALStatus *status, const char *expectedCodes, int exitCode);

static void
ClearStatus (LALStatus *status);

static void
Usage (const char *program, int exitflag);

static void
ParseOptions (int argc, char *argv[]);

typedef enum
{
  impulse, gaussian, file
}
InputDataType;

extern char    *optarg;
extern int      optind;

int lalDebugLevel = 1;

static InputDataType    inputDataType   = gaussian;
static BOOLEAN          rhosqout        = 0;
static BOOLEAN          verbose         = 0;
static INT4             numPoints       = 32768; /*128; 32768;*/
static INT4             numSegments     = 8;
static INT4             numTmplts       = 1; /* 8; */
static INT4             numChisqBins    = 0;  /* changed default to zero */
static INT4             srate           = 8192;
static REAL4            sigmasq         = 64.0;
static INT4             invSpecTrunc    = 0;
static REAL4            fLow            = 40.0;
static REAL4            rhosqThresh     = 100.0;
static REAL4            chisqThresh     = 0.001;
static REAL4            mass            = 1.4;
static REAL4            dynRange        = 0.0; /* default = 1.0; */

static INT4		loopCount       = 0;



int
main (int argc, char *argv[])
{
  static LALStatus              status;

  INT4                          i;
  INT4                          j;
  INT4                          k;
  INT4                          l;

  INT4                          seed = 3;

  INT4                          flag = 0;
  FILE                         *fpData;
  FILE                         *fpSpec;
  FILE                         *fpResp;

  FILE                         *fpRhosq = NULL;
  FILE                         *fpRead = NULL;

  REAL4                         sigma;
  REAL4                         Sfk;
  REAL4                         respRe;
  REAL4                         respIm;
  REAL4                         deltaT;
  REAL4                         deltaF;

  REAL4                         temp;

  REAL4Vector                  *noiseVec = NULL;

  FindChirpInitParams          *initParams = NULL;

  RandomParams                 *randParams = NULL;

  FindChirpFilterParams        *filterParams = NULL;
  FindChirpDataParams        *dataParams   = NULL;
  FindChirpTmpltParams       *tmpltParams  = NULL;

  DataSegmentVector            *dataSegVec = NULL;
  DataSegment                  *dataSeg    = NULL;
  FindChirpSegmentVector       *fcSegVec   = NULL;

  FindChirpFilterInput         *filterInput = NULL;

  InspiralTemplate             *tmplt = NULL;
  InspiralEvent                *event = NULL;



  /*
   *
   * parse options, allocate memory, init params and set values
   *
   */


  ParseOptions (argc, argv);

  /* override numSegments if outputting rhosq */
  if ( rhosqout )
  {
    numSegments = 1;
    numTmplts = 1;
    fpRhosq = fopen ("rhosq.dat", "w");
  }

  initParams = (FindChirpInitParams *) LALMalloc (sizeof(FindChirpInitParams));

  initParams->numSegments       = numSegments;
  initParams->numPoints         = numPoints;
  initParams->numChisqBins      = numChisqBins;
  initParams->createRhosqVec    = rhosqout;


  LALCreateVector (&status, &noiseVec, numPoints);
  TestStatus (&status, "0", 1);
  ClearStatus (&status);


  /*
   *
   * allocate memory for segments
   *
   */


  LALCreateDataSegmentVector (&status, &dataSegVec, initParams);
  TestStatus (&status, "0", 1);
  ClearStatus (&status);

  /*  initParams->approximant = TaylorF2; */
  initParams->approximant = BCV;

  LALCreateFindChirpSegmentVector (&status, &fcSegVec, initParams);
  TestStatus (&status, "0", 1);
  ClearStatus (&status);

  LALCreateFindChirpInput (&status, &filterInput, initParams);
  TestStatus (&status, "0", 1);
  ClearStatus (&status);

  tmplt = filterInput->tmplt =
    (InspiralTemplate *) LALMalloc (sizeof(InspiralTemplate));
  memset( tmplt, 0, sizeof(InspiralTemplate) );


  /*
   *
   * initialize functions and create parameter structures
   *
   */


  LALCreateRandomParams (&status, &randParams, seed);
  TestStatus (&status, "0", 1);
  ClearStatus (&status);

  LALFindChirpFilterInit (&status, &filterParams, initParams);
  TestStatus (&status, "0", 1);
  ClearStatus (&status);

  LALFindChirpChisqVetoInit (&status, filterParams->chisqParams,
      numChisqBins, numPoints);
  TestStatus (&status, "0", 1);
  ClearStatus (&status);

  LALFindChirpDataInit (&status, &dataParams, initParams);
  TestStatus (&status, "0", 1);
  ClearStatus (&status);

  LALFindChirpTemplateInit (&status, &tmpltParams, initParams);
  TestStatus (&status, "0", 1);
  ClearStatus (&status);

  /*
   *
   * fill segments with data, spectrum and response
   *
   */


  /* analytical stuff */
  sigma  = sqrt( sigmasq );    /* standard deviation of gaussian noise*/
  deltaT = 1.0 / (float) srate;
  deltaF = 1.0 / (numPoints * deltaT);
  Sfk    = 2.0 * sigmasq * deltaT; /* used to calc gaussian noise spectrum*/
  respRe = 1.0;
  respIm = 0.0;

  /* print input parameters */
  fprintf( stdout, "\n              srate = %d\n          numPoints = %d\n",
      srate, numPoints );
  fprintf( stdout, "        numSegments = %d\n\n", numSegments );
  fprintf( stdout, "            sigma^2 = %5.2f\n", sigmasq );
  fprintf( stdout, "                Sfk = %10.8f\n", Sfk );
  fprintf( stdout, "             deltaT = %e\n             deltaF = %e\n\n",
      deltaT, deltaF);
  fprintf( stdout, "       invSpecTrunc = %d\n", invSpecTrunc );
  fprintf( stdout, "               fLow = %5.3f\n\n", fLow );
  fprintf( stdout, "     rhosqThreshold = %5.3f\n     chisqThreshold = %5.3f\n\n",
      rhosqThresh, chisqThresh);
  fprintf( stdout, "               mass = %5.2f\n\n", mass );
 fprintf( stdout, "        numSegments = %d\n\n", numSegments );
 fprintf( stdout, "        dataSegVec->length = %d\n\n",dataSegVec->length );



  for ( i = 0; (UINT4)i < dataSegVec->length; ++i )
  {
    REAL4 time;
    dataSeg = dataSegVec->data;

    dataSeg[i].chan->deltaT = (REAL8) deltaT;
    dataSeg[i].spec->deltaF = (REAL8) deltaF;

    time = numPoints * deltaT * (REAL4) i;
    dataSeg[i].chan->epoch.gpsSeconds     = (INT4) floor( time );
    time = (time - floor( time )) * 1.0E9;
    dataSeg[i].chan->epoch.gpsNanoSeconds = (INT4) floor( time );

    /* impulse */
    if ( inputDataType == impulse )
    {
      memset( dataSeg[i].chan->data->data, 0, numPoints * sizeof(REAL4) );

      dataSeg[i].chan->data->data[0] = 1.0;

      /* spectrum and response */
      for ( k = 0; k < numPoints/2 + 1; ++k )
      {
        dataSeg[i].spec->data->data[k]    = 1.0;
        dataSeg[i].resp->data->data[k].re = 1.0;
        dataSeg[i].resp->data->data[k].im = 0.0;
      }

    }
    /* gaussian noise */
    else if ( inputDataType == gaussian )
    {
      LALNormalDeviates (&status, noiseVec, randParams);
      TestStatus (&status, "0", 1);
      ClearStatus (&status);

      /* ifodmro */
      for ( j = 0; j < numPoints; ++j )
      {
        REAL4 noise;
        REAL4 ifodmro;

        noise =  0.5 + sigma * noiseVec->data[j];

        if ( noise > -2048 && noise < 2047 )
          ifodmro= noise;
        else if ( noise < -2048 )
          ifodmro = -2048.0;
        else
          ifodmro=2047.0;

        dataSeg[i].chan->data->data[j] = ifodmro;
      }

      /* spectrum and response */
      for ( k = 0; k < numPoints/2 + 1; ++k )
      {
        dataSeg[i].spec->data->data[k]    = Sfk;
        dataSeg[i].resp->data->data[k].re = 1.0;
        dataSeg[i].resp->data->data[k].im = 0.0;
      }
    }
    else if ( inputDataType == file )
    {
      fpRead = fopen ("Read.dat","w");

      /* open the input files */
      if ( !(fpData = fopen( "data.dat", "r" )) )
      {
        fprintf( stdout, "unable to open the data file for reading\n" );
        fflush( stdout );
        goto abort;
      }
      if ( !(fpSpec = fopen( "spectrum.dat", "r" )) )
      {
        fprintf( stdout, "unable to open the spectrum file for reading\n" );
        fflush( stdout );
        fclose( fpData );
        goto abort;
      }
      if ( !(fpResp = fopen( "response.dat", "r" )) )
      {
        fprintf( stdout, "unable to open the response file for reading\n" );
        fflush( stdout );
        fclose( fpData );
        fclose( fpSpec );
        goto abort;
      }

      /* read in ifodmro data */
      for ( j = 0; j < numPoints; ++j )
      {

       if (( (flag = fscanf( fpData, "%f\n",
                  &(dataSeg[i].chan->data->data[j]) )) != 1 || flag == EOF )
            && j <  numPoints )
        {
          fprintf( stdout, "error reading input data %f %d\n" , temp, j);
          fflush( stdout );
          fclose( fpData );
          fclose( fpSpec );
          fclose( fpResp );
          goto abort;
        }

      fprintf (fpRead, "%e\n",  dataSeg[i].chan->data->data[j] );

      }

      fclose (fpRead);

      /* read in spec and resp */
      for ( k = 0; k < numPoints/2 + 1; ++k )
      {
        if (( (flag = fscanf( fpSpec, "%f\n",
                  &(dataSeg[i].spec->data->data[k]) )) != 1 || flag == EOF )
            && k < numPoints/2 + 1 )
        {
          fprintf( stdout, "error reading input spectrum\n" );
          fflush( stdout );
          fclose( fpData );
          fclose( fpSpec );
          fclose( fpResp );
          goto abort;
        }

        if (( (flag = fscanf( fpResp, "%f %f\n",
                  &(dataSeg[i].resp->data->data[k].re),
                  &(dataSeg[i].resp->data->data[k].im) )) != 2 || flag == EOF )
            && k < numPoints/2 + 1 )
        {
          fprintf( stdout, "error reading input response\n" );
          fflush( stdout );
          fclose( fpData );
          fclose( fpSpec );
          fclose( fpResp );
          goto abort;
        }
      }

      /* close the files */
      fclose( fpData );
      fclose( fpSpec );
      fclose( fpResp );
    }
  }


  /*
   *
   * condition data for stationary phase chirps
   *
   */


  /* dynamic range */
  dynRange = pow( 2.0, dynRange );

  /* set parameters */
  dataParams->deltaT       = deltaT;
  dataParams->fLow         = fLow;
  dataParams->dynRange     = dynRange;
  dataParams->invSpecTrunc = invSpecTrunc;

  LALFindChirpBCVData (&status, fcSegVec, dataSegVec, dataParams);
  TestStatus (&status, "0", 1);
  ClearStatus (&status);

  /*
   *
   * loop over templates
   *
   */


  /* template generation parameters */
  tmpltParams->deltaT   = deltaT;
  tmpltParams->dynRange = dynRange;
  tmpltParams->fLow     = fLow;

  /* filter parameters */
  filterParams->deltaT            = deltaT;
  filterParams->rhosqThresh       = rhosqThresh;
  filterParams->chisqThresh       = chisqThresh;
  filterParams->computeNegFreq    = 0;
  filterParams->maximiseOverChirp = 1;

  for ( l = 0; l < numTmplts; ++l, mass +=0.01 )
  {


    /*
     *
     * generate stationary phase template
     *
     */


    /* generate dummy template parameters */
    {
      REAL4 m1 = mass+5;
      REAL4 m2 = mass+5;
      tmplt->mass1     = m1;
      tmplt->mass2     = m2;
      tmplt->totalMass = m1 + m2;
      tmplt->mu        = m1 * m2 / tmplt->totalMass;
      tmplt->eta       = tmplt->mu / tmplt->totalMass;
      tmplt->approximant = BCV;

      tmplt->fFinal 	= 1000;
      tmplt->psi0       = 205008;
      tmplt->psi3       = -1619;
      tmplt->beta       = 0;

      tmpltParams->fLow  = 40;
      tmpltParams->deltaT = deltaT;
}

    LALFindChirpBCVTemplate (&status, filterInput->fcTmplt, tmplt, tmpltParams);
    TestStatus (&status, "0", 1);
    ClearStatus (&status);


    /*
     *
     * filter the segments
     *
     */

    fprintf (stdout, "numSegments %d\n", numSegments);

    for ( i = 0; i < numSegments; ++i )
    {

      fprintf (stdout, "segment number %d\n", i);

      filterInput->segment = fcSegVec->data + i;

      event = NULL;

      LALFindChirpBCVFilterSegment (&status, &event, filterInput, filterParams);
      TestStatus (&status, "0", 1);


   /*   if ( event )
      {
        fprintf( stdout, "Events found in segment!\n" );
        while ( event )
        {
          InspiralEvent *thisEvent = event;
          event = thisEvent->next;
          fprintf( stdout, "event id       = %d\n", thisEvent->id );
          fprintf( stdout, "event GPS time = %d.%d\n",
              thisEvent->time.gpsSeconds, thisEvent->time.gpsNanoSeconds );
          fprintf( stdout, "timeIndex      = %d\n", thisEvent->timeIndex );
          fprintf( stdout, "m1             = %f\n", thisEvent->tmplt.mass1 );
          fprintf( stdout, "m2             = %f\n", thisEvent->tmplt.mass2 );
          fprintf( stdout, "snrsq          = %f\n", thisEvent->snrsq );
          fprintf( stdout, "chisq          = %f\n", thisEvent->chisq );
          fprintf( stdout, "sigma          = %e\n", thisEvent->sigma );
          fprintf( stdout, "effective dist = %f\n", thisEvent->effDist);
          fflush( stdout );

          LALFree( thisEvent );
	}
      }
    */

    }

loopCount = loopCount + 1;

  } /* end loop over templates */

 fprintf (stdout, "just after end of loop over templates \n" );
 fprintf (stdout, "no of templates: %d  \n", loopCount );


  if ( rhosqout )
  {
    for ( j = 0; (UINT4)j < filterParams->rhosqVec->data->length; ++j )
    {
       fprintf( fpRhosq, "%d\t%e\n", j, filterParams->rhosqVec->data->data[j] );
/*	 fprintf( fpRhosq, "%d\t%e\n", j, pow(filterParams->rhosqVec->data->data[j],0.5) ); */
    }
  }


  /*
   *
   * finalize functions and destroy parameter structures
   *
   */


abort:

  LALDestroyRandomParams (&status, &randParams);
  TestStatus (&status, "0", 1);
  ClearStatus (&status);

  LALFindChirpChisqVetoFinalize (&status, filterParams->chisqParams,
      numChisqBins);
  TestStatus (&status, "0", 1);
  ClearStatus (&status);

  LALFindChirpFilterFinalize (&status, &filterParams);
  TestStatus (&status, "0", 1);
  ClearStatus (&status);

  LALFindChirpDataFinalize (&status, &dataParams);
  TestStatus (&status, "0", 1);
  ClearStatus (&status);

  LALFindChirpTemplateFinalize (&status, &tmpltParams);
  TestStatus (&status, "0", 1);
  ClearStatus (&status);


  /*
   *
   * free memory for segments
   *
   */


  LALDestroyFindChirpInput (&status, &filterInput);
  TestStatus (&status, "0", 1);
  ClearStatus (&status);

  LALDestroyDataSegmentVector (&status, &dataSegVec);
  TestStatus (&status, "0", 1);
  ClearStatus (&status);

  LALDestroyFindChirpSegmentVector (&status, &fcSegVec);
  TestStatus (&status, "0", 1);
  ClearStatus (&status);


  /*
   *
   * free memory for init params
   *
   */


  if ( rhosqout )
  {
    fclose( fpRhosq );
  }

  LALDestroyVector (&status, &noiseVec);
  TestStatus (&status, "0", 1);
  ClearStatus (&status);

  LALFree (initParams);
  initParams = NULL;

  LALFree (tmplt);
  tmplt = NULL;

  LALCheckMemoryLeaks ();


  return 0;
}




/*
 * TestStatus ()
 *
 * Routine to check that the status code status->statusCode agrees with one of
 * the codes specified in the space-delimited string ignored; if not,
 * exit to the system with code exitcode.
 *
 */
static void
TestStatus (
    LALStatus  *status,
    const char *ignored,
    int         exitcode
           )
{
  char  str[64];
  char *tok;

  if (verbose)
  {
    REPORTSTATUS (status);
  }

  if (strncpy (str, ignored, sizeof (str)))
  {
    if ((tok = strtok (str, " ")))
    {
      do
      {
        if (status->statusCode == atoi (tok))
        {
          return;
        }
      }
      while ((tok = strtok (NULL, " ")));
    }
    else
    {
      if (status->statusCode == atoi (tok))
      {
        return;
      }
    }
  }

  fprintf (stderr, "\nExiting to system with code %d\n", exitcode);
  exit (exitcode);
}


/*
 *
 * ClearStatus ()
 *
 * Recursively applies DETATCHSTATUSPTR() to status structure to destroy
 * linked list of statuses.
 *
 */
void
ClearStatus (
    LALStatus   *status
            )
{
  if (status->statusPtr)
  {
    ClearStatus      (status->statusPtr);
    DETATCHSTATUSPTR (status);
  }
}


static void
Usage (
    const char *program,
    int         exitcode
      )
{
  fprintf (stderr, "Usage: %s [options]\n", program);
  fprintf (stderr, "Options (defauls shown in square brackets):\n");
  fprintf (stderr, "  General:\n");
  fprintf (stderr, "    -h                print this message\n");
  fprintf (stderr, "    -V                verbose\n");
  fprintf (stderr, "    -d debuglevel     LAL status debug level [1]\n");
  fprintf (stderr, "    -o                write output files\n");
  fprintf (stderr, "  Data Creation:\n");
  fprintf (stderr, "    -I(i|g|f)         input data is (i)mpulse, (g)aussian, (f)file [g]\n");
  fprintf (stderr, "    -n numPoints      number of points in a segment [32768]\n");
  fprintf (stderr, "    -s numSegments    number of data segments [8]\n");
  fprintf (stderr, "    -r srate          sampling rate of the data [8192]\n");
  fprintf (stderr, "    -v sigmasq        variance of the gaussian noise [64.0]\n");
  fprintf (stderr, "    -m mass           mass of thingies in solar masses [1.4]\n");
  fprintf (stderr, "  Data Conditioning:\n");
  fprintf (stderr, "    -i invSpecTrunc   number of points to truncate S^{-1}(f) [0]\n");
  fprintf (stderr, "    -f fLow           low frequency cutoff for S^{-1}(f) [150.0]\n");
  fprintf (stderr, "    -y dynRange       log_2( dynamicRange ) [1.0]\n");
  fprintf (stderr, "  Filtering:\n");
  fprintf (stderr, "    -N numTmplts      number templates to filter against [8]\n");
  fprintf (stderr, "    -b numChisqBins   number of bins for chi squared test [8]\n");
  fprintf (stderr, "    -t rhosqThresh    signal to noise squared threshold [100.0]\n");
  fprintf (stderr, "    -c chisqThresh    chi squared threshold [0.0001]\n");


  exit (exitcode);
}


static void
ParseOptions (
    int         argc,
    char       *argv[]
             )
{
  while (1)
  {
    int c = -1;

    c = getopt (argc, argv, "Vhd:""n:""s:""y:""r:""v:""i:""f:""R:""t:""b:""m:""N:""c:""oI:");
    if (c == -1)
    {
      break;
    }

    switch (c)
    {
      case 'd': /* set debuglevel */
        lalDebugLevel = atoi (optarg);
        break;
      case 'V': /* set verbosity */
        verbose = 1;
        break;
      case 'o': /* sets flag to write rhosq */
        rhosqout = 1;
        break;
      case 'n': /* set number of points in a segment */
        numPoints = atoi (optarg);
        break;
      case 's': /* set number of segments */
        numSegments = atoi (optarg);
        break;
      case 'r': /* set sampling rate */
        srate = atoi (optarg);
        break;
      case 'i': /* set invSpecTrunc */
        invSpecTrunc = atoi (optarg);
        break;
      case 'N': /* set numChisqBins */
        numTmplts = atoi (optarg);
        break;
      case 'b': /* set numChisqBins */
        numChisqBins = atoi (optarg);
        break;
      case 'f': /* set fLow */
        fLow = (REAL4) atof (optarg);
        break;
      case 'y': /* set dynRange */
        dynRange = (REAL4) atof (optarg);
        break;
      case 'm': /* mass */
        mass = (REAL4) atof (optarg);
        break;
      case 't': /* set rhosq threshold */
        rhosqThresh = (REAL4) atof (optarg);
        break;
      case 'c': /* set chisq threshold */
        chisqThresh = (REAL4) atof (optarg);
        break;
      case 'v': /* set variance */
        sigmasq = (float) atof (optarg);
        break;
      case 'I': /* input data */
        switch (*optarg)
        {
          case 'i':
            inputDataType = impulse;
            break;
          case 'g':
            inputDataType = gaussian;
            break;
          case 'f':
            inputDataType = file;
            break;
          default:
            inputDataType = gaussian;
        }
        break;
      case 'h':
        Usage (argv[0], 0);
        break;
      default:
        Usage (argv[0], 1);
    }
  }

  if (optind < argc)
  {
    Usage (argv[0], 1);
  }

  return;
}


/*
 *
 * graphREAL4 ()
 *
 * function to graph an array of REAL4's for debugging
 *
 */
static void
graphREAL4 (
    REAL4      *array,
    INT4        n,
    INT4        spacing
           )
{
  FILE *fp;
  INT4 i;

  /* open a file for writing */
  if ( !(fp = fopen( "temp.graph", "w" )) )
  {
    printf( "couldn't open file\n" );
  }

  /* print data into the file */
  for ( i = 0; i < n; i++ )
    fprintf( fp, "%d\t%e\n", i, array[i * spacing] );

  /* close the file */
  fclose( fp );

  /* start up graphing program with data in the file */
  /* system( "xmgr temp.graph 1>/dev/null 2>&1 &" ); */

  return;
}


/*
 *
 * graphINT2 ()
 *
 * function to graph an array of INT2's for debugging
 *
 */
static void
graphINT2 (
    INT2       *array,
    INT4        n,
    INT4        spacing
          )
{
  FILE *fp;
  INT4 i;

  /* open a file for writing */
  if ( !(fp = fopen( "temp.graph", "w" )) )
  {
    printf( "couldn't open file\n" );
  }

  /* print data into the file */
  for ( i = 0; i < n; i++ )
    fprintf( fp, "%d\t%d\n", i, array[i * spacing] );

  /* close the file */
  fclose( fp );

  /* start up graphing program with data in the file */
  /* system( "xmgr temp.graph 1>/dev/null 2>&1 &" ); */

  return;
}


/*
 *
 * graphINT4 ()
 *
 * function to graph an array of INT4's for debugging
 *
 */
static void
graphINT4 (
    INT4       *array,
    INT4        n,
    INT4        spacing
          )
{
  FILE *fp;
  INT4 i;

  /* open a file for writing */
  if ( !(fp = fopen( "temp.graph", "w" )) )
  {
    printf( "couldn't open file\n" );
  }

  /* print data into the file */
  for ( i = 0; i <n; i++ )
    fprintf( fp, "%d\t%d\n", i, array[i * spacing] );

  /* close the file */
  fclose( fp );

  /* start up graphing program with data in the file */
  /* system( "xmgr temp.graph 1>/dev/null 2>&1 &" ); */

return;
}
