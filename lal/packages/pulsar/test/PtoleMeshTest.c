/************************************** <lalVerbatim file="PtoleMeshTestCV">
Author: Owen, B. J.
$Id$
********************************************************** </lalVerbatim> */
 
/**************************************************************** <lalLaTeX>
 
\subsection{Program \texttt{PtoleMeshTest}}
\label{ss:PtoleMeshTest}
 
Tests and showcases the combination of \texttt{PtoleMetric} and
\texttt{TwoDMesh} modules.
 
\subsubsection*{Usage}
\begin{verbatim}
PtoleMeshTest
\end{verbatim}
 
\subsubsection*{Description}

The \texttt{-b} option sets the beginning time of integration to the option
argument. (Default is $7\times10^8$ seconds.)

The \texttt{-c} option determins the center of the patch. (Default is the
center of the globular cluster 47 Tuc.) This option is hardcoded to use
equatorial coordinates and the argument should be given in hh:mm:ss:dd:mm:ss
format.

The \texttt{-e} option sets \texttt{lalDebugLevel} to the option argument.
(Default is 1.)

The \texttt{-f} option sets the maximum frequency of integration (in Hz) to the
option argument. (The default value is 1000.)

The \texttt{-i} option does not function at this time.

The \texttt{-m} option sets the maximum mismatch of the mesh to the option
argument. (Default is 0.02.)

The \texttt{-n} option sets the maximum number of nodes in the mesh to the
option argument. (Default is $10^6$.)

The \texttt{-r} option sets the radius (in decimal degrees) of the circular sky
patch. (The default value is set for the globular cluster 47 Tuc.) At the
moment there is no option for another patch shape.

The \texttt{-t} option sets the duration of integration, in seconds. (The
default is $10^5$ seconds, which is of order one day but is not an integer
multiple of anything astronomically important.)

The \texttt{-x} option makes a plot of the mesh points on the sky patch using a
system call to \texttt{xmgrace}. If \texttt{xmgrace} is not installed on your
system, this option will not work. The plot goes to a file \texttt{mesh.agr}.

\subsubsection*{Exit Codes}
************************************************ </lalLaTeX><lalErrTable> */
#define PTOLEMESHTESTC_EMEM 1
#define PTOLEMESHTESTC_ERNG 2
#define PTOLEMESHTESTC_EFIO 3
#define PTOLEMESHTESTC_EOPT 4

#define PTOLEMESHTESTC_MSGEMEM "memory (de)allocation error"
#define PTOLEMESHTESTC_MSGERNG "value out of range"
#define PTOLEMESHTESTC_MSGEFIO "file I/O error"
#define PTOLEMESHTESTC_MSGEOPT "unknown command-line option"
/************************************************** </lalErrTable><lalLaTeX>
 
\subsubsection*{Algorithm}
 
\subsubsection*{Uses}

\begin{verbatim}
lalDebugLevel
LALCheckMemoryLeaks()
LALProjectMetric()
LALPtoleMetric()
LALXMGRPlotMesh()
\end{verbatim}
 
\subsubsection*{Notes}

\vfill{\footnotesize\input{PtoleMeshTestCV}}
 
************************************************************* </lalLaTeX> */


#include <math.h>
#include <stdio.h>
#include <lal/AVFactories.h>
#include <lal/LALXMGRInterface.h>
#include <lal/PtoleMetric.h>
#include <lal/StackMetric.h>
#include <lal/TwoDMesh.h>

NRCSID( PTOLEMESHTESTC, "$Id$" );

/* BEN: These aren't used right now, but should be. */
#define MIN_DURATION (86400./LAL_TWOPI) /* one radian of rotation */
#define MAX_DURATION (86400.*5)         /* not much orbital motion */
#define MIN_FREQ     1e2                /* more or less arbitrary */
#define MAX_FREQ     1e4                /* more or less arbitrary */

char *optarg = NULL;     /* option argument for getopt() */
int  lalDebugLevel = 1;  /* default value */

void getRange( LALStatus *, REAL4 [2], REAL4, void * );
void getMetric( LALStatus *, REAL4 [3], REAL4 [2], void * );

/* BEN: This is a cheat. Should make a search area structure. */
static SkyPosition center;  /* center of search */
REAL4              radius;  /* radius of search, in arcminutes */


int main( int argc, char **argv )
{
  static LALStatus     stat;      /* status structure */
  INT2                 opt;       /* command-line option character */
  BOOLEAN              errors;    /* whether or not to showcase error traps */
  BOOLEAN              grace;     /* whether or not to graph using xmgrace */
  TwoDMeshNode        *firstNode; /* head of linked list of nodes in mesh */
  static TwoDMeshParamStruc mesh; /* mesh parameters */
  static PtoleMetricIn search;    /* more mesh parameters */
  REAL4                mismatch;  /* mismatch threshold of mesh */
  UINT4                maxNodes;  /* maximum nodes in mesh */
  REAL4                begin;     /* start time of integration (seconds) */
  REAL4                duration;  /* duration of integration (seconds) */
  REAL4                fMax;      /* maximum frequency of search */
  LALFrDetector        site;      /* detector site */
  FILE                *fp;        /* where to write a plot */

  /* Set default values. */
  errors = 0; /* BEN: this is unused right now */
  grace = 0;
  begin = 7e8;
  duration = 1e5;
  fMax = 1e3;
  site = lalCachedDetectors[LALDetectorIndexGEO600DIFF].frDetector;
  mismatch = .02;
  maxNodes = 1e6;
  /* This is (roughly) the center of globular cluster 47 Tuc. */
  center.system = COORDINATESYSTEM_EQUATORIAL;
  center.longitude = (24.1/60)*LAL_PI_180*15;
  center.latitude = -(72+5./60)*LAL_PI_180;
  radius = 30.9/60*LAL_PI_180;

  /* Parse command-line options. */
  while( (opt = getopt( argc, argv, "b:c:e:f:i:m:n:t:x" )) != -1 )
  {
    switch( opt )
    {
      float a, b, c, d, e, f;
    case '?':
      return PTOLEMESHTESTC_EOPT;
    case 'b':
      begin = atof( optarg );
      break;
    case 'c':
      if( sscanf( optarg, "%f:%f:%f:%f:%f:%f", &a, &b, &c, &d, &e, &f ) != 6)
      {
        fprintf( stderr, "coordinates should be hh:mm:ss:dd:mm:ss\n" );
        return PTOLEMESHTESTC_EOPT;
      }
      center.longitude = (a+b/60+c/3600)*LAL_PI_180*15;
      center.latitude = (d+e/60+f/3600)*LAL_PI_180;
      break;
    case 'e':
      lalDebugLevel = atof( optarg );
      break;
    case 'f':
      fMax = atof( optarg );
      break;
    case 'i':
      break;
    case 'm':
      mismatch = atof( optarg );
      break;
    case 'n':
      maxNodes = atoi( optarg );
      break;
    case 'r':
      radius = LAL_PI_180/60*atof( optarg );
      break;
    case 't':
      duration = atof( optarg );
      break;
    case 'x':
      grace = 1;
      break;
    } /* switch( opt ) */
  } /* while( getopt... ) */

  /* Create the mesh. */
  mesh.mThresh = mismatch;
  mesh.nIn = maxNodes;
  mesh.getRange = getRange;
  mesh.rangeParams = NULL;
  mesh.getMetric = getMetric;
  mesh.metricParams = (void *) &search;
  mesh.domain[0] = center.longitude - radius;
  mesh.domain[1] = center.longitude + radius;
  search.position.system = COORDINATESYSTEM_EQUATORIAL;
  search.spindown = NULL;
  search.epoch.gpsSeconds = begin;
  search.epoch.gpsNanoSeconds = 0;
  search.duration = duration;
  search.maxFreq = fMax;
  search.site = site;
  firstNode = NULL;
  LALCreateTwoDMesh( &stat, &firstNode, &mesh );
  if( stat.statusCode )
    return stat.statusCode;
  printf( "created %d nodes\n", mesh.nOut );

  /* Plot what we've got. */
  if( grace )
  {
    fp = fopen( "mesh.agr", "w" );
    if( !fp )
      return PTOLEMESHTESTC_EFIO;
    LALXMGRPlotMesh( &stat, firstNode, fp, &mesh );
    if( stat.statusCode )
      return stat.statusCode;
    fclose( fp );
  }

  /* Clean up and leave. */
  LALDestroyTwoDMesh( &stat, &firstNode, &mesh.nOut );
  printf( "destroyed %d nodes\n", mesh.nOut );
  if( stat.statusCode )
    return PTOLEMESHTESTC_EMEM;
  LALCheckMemoryLeaks();
  return 0;
} /* main() */


/* This is the parameter range function as required by TwoDMesh. */

void getRange( LALStatus *stat, REAL4 y[2], REAL4 x, void *unused )
{
  /* Set up shop. */
  INITSTATUS( stat, "getRange", PTOLEMESHTESTC );
  ATTATCHSTATUSPTR( stat );
  unused = NULL;

  /* Search a circle. BEN: The 1.001 is a kludge. */
  y[0] = center.latitude - sqrt( pow( radius*1.001, 2 )
         - pow( x-center.longitude, 2 ) );
  y[1] = center.latitude + sqrt( pow( radius*1.001, 2 )
         - pow( x-center.longitude, 2 ) );

  /* Clean up and leave. */
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
} /* getRange() */


/* This is the wrapped metric function as required by TwoDMesh. */

void getMetric( LALStatus *stat,
                REAL4 g[3],
                REAL4 x[2],
                void *params )
{
  PtoleMetricIn *patch = params; /* avoid pesky gcc warnings */
  REAL8Vector   *metric = NULL;  /* for output of metric */

  /* Set up shop. */
  INITSTATUS( stat, "getMetric", PTOLEMESHTESTC );
  ATTATCHSTATUSPTR( stat );
  TRY( LALDCreateVector( stat->statusPtr, &metric, 10 ), stat );

  /* Translate input. */
  patch->position.longitude = x[0];
  patch->position.latitude = x[1];

  /* Call the real metric function. */
  LALPtoleMetric( stat->statusPtr, metric, patch );
  BEGINFAIL( stat )
    TRY( LALDDestroyVector( stat->statusPtr, &metric ), stat );
  ENDFAIL( stat );
  LALProjectMetric( stat->statusPtr, metric, 0 );
  BEGINFAIL( stat )
    TRY( LALDDestroyVector( stat->statusPtr, &metric ), stat );
  ENDFAIL( stat );

  /* Translate output. */
  g[0] = metric->data[2];
  g[1] = metric->data[5];
  g[2] = metric->data[4];

  /* Clean up and leave. */
  TRY( LALDDestroyVector( stat->statusPtr, &metric ), stat );
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
} /* getMetric() */
