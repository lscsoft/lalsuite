/*
*  Copyright (C) 2007 Sukanta Bose
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
 * File Name: CoherentInspiralFilterTest.c
 *
 * Author: Seader, S. E., Bose, S.
 *
 * Revision: $Id$
 *
 *-----------------------------------------------------------------------
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#ifdef HAVE_GETOPT_H
/* #include <getopt.h> */
#endif
#include <getopt.h>

#include <lal/LALRCSID.h>
#include <lal/LALConfig.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/AVFactories.h>
#include <lal/Random.h>
#include <lal/LALInspiral.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpSP.h>
#include <lal/DetectorSite.h>
#include <lal/LALStatusMacros.h>
#include <lal/StochasticCrossCorrelation.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/CoherentInspiral.h>


NRCSID(MAIN,"$Id$");

static void
ParseOptions (int argc, char *argv[]);

static void
TestStatus (LALStatus *status, const char *expectedCodes, int exitCode);

static void
ClearStatus (LALStatus *status);

static void
Usage (const char *program, int exitflag);

typedef enum{
  impulse, gaussian, file
    }
InputDataType;

extern char    *optarg;
extern int      optind;
int lalDebugLevel = 1;

static FILE *fp[4];
static FILE *fp2[4];

/* Parse Options - Hard code these and provide error checking */

char H1filename[256];
char H2filename[256];
char L1filename[256];
char GEOfilename[256];
char VIRGOfilename[256];
char TAMAfilename[256];

char H1Beam[256];
char H2Beam[256];
char L1Beam[256];
char GEOBeam[256];
char VIRGOBeam[256];
char TAMABeam[256];

UINT4 siteID[6] = {0,1,2,3,4,0}; /*H1 L V G T H2*/
UINT2 caseID[6] = {0,0,0,0,0,0};

static BOOLEAN          cohSNROut            = 1;
static BOOLEAN          verbose              = 0;

static INT4             numTmplts            = 1;
static UINT4            numSegments          = 1;
static UINT4            numPoints            = 262144;
static UINT4            numBeamPoints        = 2;
static UINT4            sampleRate           = 1024;
static REAL4            cohSNRThresh         = 10.0;
static REAL4            fLow                 = 70.0;
static REAL4            mass                 = 1.4;
static UINT4            startSeconds         = 10000;
static UINT4            maximiseOverChirp    = 1;

int
main (int argc, char *argv[])
{
  static LALStatus              status;


  INT4                          i,j,k,l;

  UINT4                         numDetectors=0;
  FILE                         *fpcohSNR;
  REAL4                         theta,phi,vPlus,vMinus;
  REAL4                         x,y;

  UINT2Vector                  *detIDVec;
  DetectorVector               *detectorVec;

  CoherentInspiralInitParams   *cohInspInitParams = NULL;
  CoherentInspiralFilterParams *cohInspFilterParams = NULL;
  CoherentInspiralFilterInput  *cohInspFilterInput = NULL;
  CoherentInspiralBeamVector   *cohInspBeamVec = NULL;
  CoherentInspiralZVector      *cohInspZVec = NULL;
  InspiralTemplate             *tmplt = NULL;
  CoherentInspiralEvent        *cohInspEvent = NULL;

  char namearray[6][256] = {"0","0","0","0","0","0"};
  char namearray2[6][256] = {"0","0","0","0","0","0"};

  /*
   *
   * parse options, allocate memory, init params and set values
   *
   */

  ParseOptions (argc, argv);


  /* override numSegments if outputting coherent SNR */
  if ( cohSNROut )
    {
      numSegments = 1;
      numTmplts = 1;
      fpcohSNR = fopen ("cohSNR.dat", "w");
    }


  /* read in the network detectors */
  for (l=0;l<6;l++)
    {
      if(caseID[l] == 1)
	numDetectors++;
    }

  fprintf(stdout, "You have specified %2d  detector(s).\n",numDetectors);
  fprintf(stdout, "The caseID is: %d %d %d %d %d %d (H1,H2,L1,GEO,VIRGO,TAMA) \n",caseID[0],caseID[1],caseID[2],caseID[3],caseID[4],caseID[5]);

  if (numDetectors > 4)
    {
    fprintf(stdout, "Too many detectors specified - exiting. \n");
    goto cleanexit;
    }
  if (numDetectors == 0)
    {
      fprintf(stdout, "You must specify data filename(s) for 1 to 4 detectors - exiting. \n");
      goto cleanexit;
    }


  /*
   *
   * allocate memory to structures
   *
   */

  /* fill the init params structure */
  cohInspInitParams = (CoherentInspiralInitParams *)
    LALMalloc (sizeof(CoherentInspiralInitParams));

  cohInspInitParams->numDetectors            = numDetectors;
  cohInspInitParams->numSegments             = numSegments;
  cohInspInitParams->numPoints               = numPoints;
  cohInspInitParams->numBeamPoints           = numBeamPoints;
  cohInspInitParams->cohSNROut               = cohSNROut;

  /* Create input structure for coherent filter code */
  LALCoherentInspiralFilterInputInit (&status, &cohInspFilterInput,
					cohInspInitParams);
  TestStatus (&status, "0", 1);
  ClearStatus (&status);



  /*
   * information for calculating chirp time
   */

  /* inspiral template structure */
  tmplt = cohInspFilterInput->tmplt = (InspiralTemplate *)
    LALMalloc (sizeof(InspiralTemplate));
  memset( tmplt, 0, sizeof(InspiralTemplate) );

  /* generate dummy template parameters */
  {
    REAL4 m1 = mass;
    REAL4 m2 = mass;
    tmplt->mass1     = m1;
    tmplt->mass2     = m2;
    tmplt->totalMass = m1 + m2;
    tmplt->mu        = m1 * m2 / tmplt->totalMass;
    tmplt->eta       = tmplt->mu / tmplt->totalMass;
  }

  /* fill the params structure */
  LALCoherentInspiralFilterParamsInit (&status, &cohInspFilterParams,
				       cohInspInitParams);
  TestStatus (&status, "0", 1);
  ClearStatus (&status);

  cohInspFilterParams->numDetectors            = numDetectors;
  cohInspFilterParams->numSegments             = numSegments;
  cohInspFilterParams->numPoints               = numPoints;
  cohInspFilterParams->numBeamPoints           = numBeamPoints;
  cohInspFilterParams->deltaT                  = 1/((REAL4) sampleRate);
  cohInspFilterParams->cohSNROut               = cohSNROut;
  cohInspFilterParams->cohSNRThresh            = cohSNRThresh;
  cohInspFilterParams->numTmplts               = numTmplts;
  cohInspFilterParams->fLow                    = fLow;
  cohInspFilterParams->maximiseOverChirp       = maximiseOverChirp;

  detIDVec = cohInspFilterParams->detIDVec;

  /*assign detIDs to the coincident detectors in the network */
  for ( i=0 ; i < 6 ; i++) {
    detIDVec->data[i] = caseID[i];
  }

 /* create and fill the DetectorVector structure of detector IDs*/

  detectorVec = cohInspFilterParams->detectorVec;

  i=0;
  for ( j=0 ; j < 6 ; j++ ) {
    /*    if (((j != 5) && caseID[j++])) { */
    if ( caseID[j] ) {
      detectorVec->detector[i++] = lalCachedDetectors[j];
    }
  }

  if (caseID[5]) {
    detectorVec->detector[numDetectors-1] = lalCachedDetectors[0];
  }


  /* Now read in all the filenames and store them in arrays */
  /* This will keep the files in the order:
     H1(0), L1(1), VIRGO(2), GEO(3), TAMA(4), H2(5)*/

  if(caseID[0])
    {
    strcpy(namearray[0],H1filename);
    strcpy(namearray2[0],H1Beam);
    }

  if(caseID[1])
    {
    strcpy(namearray[1],L1filename);
    strcpy(namearray2[1],L1Beam);
    }

  if(caseID[2])
    {
    strcpy(namearray[2],VIRGOfilename);
    strcpy(namearray2[2],VIRGOBeam);
    }

  if(caseID[3])
    {
    strcpy(namearray[3],GEOfilename);
    strcpy(namearray2[3],GEOBeam);
    }

  if(caseID[4])
    {
    strcpy(namearray[4],TAMAfilename);
    strcpy(namearray2[4],TAMABeam);
    }

  if (caseID[5])
    {
    strcpy(namearray[5],H2filename);
    strcpy(namearray2[5],H2Beam);
    }



  /* create and fill the CoherentInspiralBeamVector structure of beam-patterns*/
  cohInspBeamVec = cohInspFilterInput->beamVec;

  l=0;
  for ( j=0 ; j < 6 ; j++ ) {
    if ( caseID[j] ) {
      /*      for (l=0;l<numDetectors;l++) { */

      fp2[l] = fopen(namearray2[j], "r");
      if(!fp2[l])
	{
	  fprintf(stdout,"The file %s containing the coefficients could not be found - exiting...\n",namearray2[j]);
	  goto  cleanexit;
	}
      for ( k=0 ; k<numBeamPoints ; k++)
	{
	  fscanf(fp2[l],"%f, %f, %f, %f",&theta,&phi,&vPlus,&vMinus);
	  cohInspBeamVec->detBeamArray[l].thetaPhiVs[k].data->data[0] = theta;
	  cohInspBeamVec->detBeamArray[l].thetaPhiVs[k].data->data[1] = phi;
	  cohInspBeamVec->detBeamArray[l].thetaPhiVs[k].data->data[2] = vPlus;
	  cohInspBeamVec->detBeamArray[l].thetaPhiVs[k].data->data[3] = vMinus;
	}
      fclose(fp2[l++]);
    }
  }


  /*
   * CREATE the multi-z data structure;
   * the z=x+iy data from multiple detectors was read above along with
   * the beam-pattern functions
   */
  cohInspZVec = cohInspFilterInput->multiZData;

  /* First, the files will be tested for length consistency
     and then the z-data for multiple detectors will be read in */
  l=0;
  for ( j=0 ; j < 6 ; j++ ) {
    if ( caseID[j] ) {
      fp[l] = fopen(namearray[j], "r");
      if(!fp[l])
	{
	  fprintf(stdout,"The file %s does not exist - exiting...\n",
		  namearray[j]);
	  goto cleanexit;
	}
      for (k = 0; k < numPoints; k++)
	{
	  fscanf(fp[l],"%f %f", &x, &y);
	  cohInspZVec->zData[l].data->data[k].re = x;
	  cohInspZVec->zData[l].data->data[k].im = y;
	}
      fclose(fp[l++]);
    }
  }

  /*Do the filtering and output events */
  cohInspEvent = NULL;

  LALCoherentInspiralFilterSegment (&status, &cohInspEvent, cohInspFilterInput, cohInspFilterParams);
  TestStatus (&status, "0", 1);

  if ( cohInspEvent )
	{
	  fprintf( stdout, "\nEvents found in segment!\n\n" );
	  while (cohInspEvent  )
	    {
	      CoherentInspiralEvent *thisEvent = cohInspEvent;

	      cohInspEvent = thisEvent->next;

	      fprintf( stdout, "event id                              = %d\n\n", thisEvent->eventId+1 );

	      fprintf( stdout, "coherent SNR                         = %.2f\n", thisEvent->cohSNR );

	      fprintf( stdout, "'network' timeIndex                   = %d\n", thisEvent->timeIndex );

	      fflush( stdout );

	      LALFree( thisEvent );
	    }

	}
      else
	{
	  fprintf( stdout, "\nNo events found in segment for mass = %.4f solar mass template!\n", mass );
	}
  /* outputting coherent SNR */
  if ( cohSNROut )
    {
      for ( i = 0; i < cohInspFilterParams->cohSNRVec->data->length; ++i )
	{
	  fprintf( fpcohSNR, "%d\t%e\n", i, cohInspFilterParams->cohSNRVec->data->data[i] );
	}
    }

 cleanexit:

  fclose( fpcohSNR);

  /* Destroy params structure for coherent filter code */
  LALCoherentInspiralFilterParamsFinalize (&status, &cohInspFilterParams);
  TestStatus (&status, "0", 1);
  ClearStatus (&status);

  /* Destroy input structure for coherent filter code */
  LALCoherentInspiralFilterInputFinalize (&status, &cohInspFilterInput);
  TestStatus (&status, "0", 1);
  ClearStatus (&status);


  LALCheckMemoryLeaks ();

  return 0;

} /* end main */



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
  fprintf (stderr, "Options (defaults shown in square brackets):\n");
  fprintf (stderr, "  General:\n");
  fprintf (stderr, "    -w                verbose\n");
  fprintf (stderr, "    -h                print this message\n");
  fprintf (stderr, "    -d debuglevel     LAL status debug level [1]\n");
  fprintf (stderr, "    -o                write coherent SNR output file [1]\n");
  fprintf (stderr, "    -n numPoints      number of points in a segment [1048576]\n");
  fprintf (stderr, "    -b numBeamPoints  number of theta-phi template points \n");
  fprintf (stderr, "    -c maximiseOverChirp do clustering in chirp time [1]\n");
  fprintf (stderr, "    -N numTmplts      number of mass templates [1]\n");
  fprintf (stderr, "    -f fLow           low-frequency cut-off [70.0]\n");
  fprintf (stderr, "    -s numSegments    number of data segments [1]\n");
  fprintf (stderr, "    -r sampleRate     sampling rate of the z-data [8192]\n");
  fprintf (stderr, "    -p cohSNRThresh   coherent SNR threshold [10.0]\n");
  fprintf (stderr, "    -A site0          Use H1.dat and site number [0]\n");
  fprintf (stderr, "    -L site1          Use L1.dat and site number [1]\n");
  fprintf (stderr, "    -V site2          Use V.dat and site number [2]\n");
  fprintf (stderr, "    -G site3          Use G.dat and site number [3]\n");
  fprintf (stderr, "    -T site4          Use T.dat and site number [4]\n");
  fprintf (stderr, "    -Z H2site         Use H2.dat and site number [5]\n");
  fprintf (stderr, "    -a H1Beam         Use H1Beam.dat and site number [0]\n");
  fprintf (stderr, "    -l L1Beam         Use L1Beam.dat and site number [1]\n");
  fprintf (stderr, "    -v VBeam          Use VBeam.dat and site number [2]\n");
  fprintf (stderr, "    -g GBeam          Use GBeam.dat and site number [3]\n");
  fprintf (stderr, "    -t TBeam          Use TBeam.dat and site number [4]\n");
  fprintf (stderr, "    -z H2Beam         Use H2Beam.dat and site number [5]\n");
  exit (exitcode);
}


/* may want to add a verbose option flag to print out extra info here and there */
static void
ParseOptions (
	      int         argc,
	      char       *argv[]
	      )
{
  while (1)
    {
      int c = -1;

      c = getopt (argc, argv, "whd:""A:""L:""V:""G:""T:""Z:""a:""l:""v:""g:""t:""z:""o:""n:""b:""c:""s:""N:""f:""r:""p:");
      if (c == -1)
	{
	  break;
	}

      switch (c)
	{
	case 'w': /* set verbosity */
	  verbose = 1;
	  break;
	case 'h':
	  Usage (argv[0], 0);
	  break;
	case 'd': /* set debuglevel */
	  lalDebugLevel = atoi (optarg);
	  break;
  	case 'A':
	  strcpy(H1filename,optarg);
	  caseID[0] = 1;
	  break;
	case 'L':
	  strcpy(L1filename,optarg);
	  caseID[1] = 1;
	  break;
	case 'V':
	  strcpy(VIRGOfilename,optarg);
	  caseID[2] = 1;
	  break;
	case 'G':
	  strcpy(GEOfilename,optarg);
	  caseID[3] = 1;
	  break;
	case 'T':
	  strcpy(TAMAfilename,optarg);
	  caseID[4] = 1;
	  break;
	case 'Z':
	  strcpy(H2filename,optarg);
	  caseID[5] = 1;
	  break;
  	case 'a':
	  strcpy(H1Beam,optarg);
	  break;
	case 'l':
	  strcpy(L1Beam,optarg);
	  break;
	case 'v':
	  strcpy(VIRGOBeam,optarg);
	  break;
	case 'g':
	  strcpy(GEOBeam,optarg);
	  break;
	case 't':
	  strcpy(TAMABeam,optarg);
	  break;
	case 'z':
	  strcpy(H2Beam,optarg);
	  break;
	case 'o': /* sets flag to write coherent SNR */
	  cohSNROut = 1;
	  break;
	case 'n': /* set number of points in a segment */
	  numPoints = atoi (optarg);
	  break;
	case 'b': /* set number of theta-phi template points */
	  numBeamPoints = atoi (optarg);
	  break;
	case 'c': /* set maximiseOverChirp */
	  maximiseOverChirp  = 1;
	  break;
	case 's': /* set number of segments */
	  numSegments = atoi (optarg);
	  break;
	case 'N': /* set number of templates */
	  numTmplts = atoi (optarg);
	  break;
	case 'f': /* set fLow */
	  fLow = (REAL4) atof (optarg);
	  break;
	case 'r': /* set sampling rate */
	  sampleRate = atoi (optarg);
	  break;
	case 'p': /* set coherent SNR threshold */
	  cohSNRThresh = atof (optarg);
	  break;
        default:
	  Usage (argv[0], 1);
	  fprintf( stderr, "unknown error while parsing options\n" );
	  exit( 1 );

	}
    }

  if (optind < argc)
    {
      Usage (argv[0], 1);
    }

  return;
}
