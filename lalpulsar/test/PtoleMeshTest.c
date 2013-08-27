/*
*  Copyright (C) 2007 Jolien Creighton, Ian Jones, Benjamin Owen, Reinhard Prix
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
\author Owen, B. J.,   Jones, D. I.
\file
\ingroup PtoleMetric_h
\brief Tests and showcases the combination of \ref PtoleMetric_h and \ref TwoDMesh_h modules.

\heading{Program \c PtoleMeshTest}
\latexonly\tag{ss_PtoleMeshTest}\endlatexonly

\heading{Usage}
\code
PtoleMeshTest
\endcode

\heading{Description}

The <tt>-b</tt> option sets the beginning time of integration to the option
argument. (Default is \f$0\f$ seconds)

The <tt>-c</tt> option determins the center of the patch. (Default is the
center of the globular cluster 47 Tuc.) This option is hardcoded to use
equatorial coordinates and the argument should be given in hh:mm:ss:dd:mm:ss
format.

The <tt>-e</tt> option sets \c lalDebugLevel to the option argument.
(Default is 1.)

The <tt>-f</tt> option sets the maximum frequency of integration (in Hz) to the
option argument. (The default value is 1000.)

The <tt>-i</tt> option does not function at this time.

The <tt>-m</tt> option sets the maximum mismatch of the mesh to the option
argument. (Default is 0.02.)

The <tt>-n</tt> option sets the maximum number of nodes in the mesh to the
option argument. (Default is \f$10^6\f$.)

The texttt{-p} option causes the coordinates of the nodes to be written to
a file <tt>mesh.dat</tt>, for the benifit of users who don't have
\c xmgrace installed.  The format is one node per line, (RA, DEC),
with the angles in degrees.

The <tt>-r</tt> option sets the radius (in arcminutes) of the circular
sky patch. (The default value is set for the globular cluster 47 Tuc.)
At the moment there is no option for another patch shape, but if you
specify radius zero you will get a search over a rectangular region
whose limits in RA and dec are specified in the code.

The <tt>-t</tt> option sets the duration of integration, in seconds. (The
default is \f$10^5\f$ seconds, which is of order one day but is not an integer
multiple of anything astronomically important.)

The <tt>-x</tt> option makes a plot of the mesh points on the sky patch using a
system call to \c xmgrace. If \c xmgrace is not installed on your
system, this option will not work. The plot goes to a file <tt>mesh.agr</tt>.

\heading{Algorithm}

\heading{Uses}

\code
lalDebugLevel
LALCheckMemoryLeaks()
LALProjectMetric()
LALPtoleMetric()
LALXMGRPlotMesh()
\endcode

\heading{Notes}

*/

/** \name Error Codes */ /*@{*/
#define PTOLEMESHTESTC_EMEM 1
#define PTOLEMESHTESTC_ERNG 2
#define PTOLEMESHTESTC_EFIO 3
#define PTOLEMESHTESTC_EOPT 4

#define PTOLEMESHTESTC_MSGEMEM "memory (de)allocation error"
#define PTOLEMESHTESTC_MSGERNG "value out of range"
#define PTOLEMESHTESTC_MSGEFIO "file I/O error"
#define PTOLEMESHTESTC_MSGEOPT "unknown command-line option"
/*@}*/

#include <config.h>

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#ifdef HAVE_GETOPT_H
#include <getopt.h>
#endif

#include <math.h>
#include <stdio.h>
#include <lal/AVFactories.h>
#include <lal/LALXMGRInterface.h>
#include <lal/PtoleMetric.h>
#include <lal/StackMetric.h>
#include <lal/TwoDMesh.h>

/** \cond DONT_DOXYGEN */

/* BEN: These aren't used right now, but should be. */
#define MIN_DURATION (86400./LAL_TWOPI) /* one radian of rotation */
#define MAX_DURATION (86400.*5)         /* not much orbital motion */
#define MIN_FREQ     1e2                /* more or less arbitrary */
#define MAX_FREQ     1e4                /* more or less arbitrary */

/* IAN: Clumsy way of specifying rectangular search region if */
/* the r=0 option is invoked.                                 */

#define RA_min       0.0
#define RA_max       LAL_PI_2/2.0
#define dec_min      LAL_PI_2/2.0
#define dec_max      LAL_PI_2




extern char *optarg;     /* option argument for getopt() */

void getRange( LALStatus *, REAL4 [2], REAL4, void * );
void getMetric( LALStatus *, REAL4 [3], REAL4 [2], void * );

/* BEN: This is a cheat. Should make a search area structure. */
static SkyPosition center;  /* center of search */
REAL4              radius;  /* radius of search, in arcminutes */

int main( int argc, char **argv )
{
  static LALStatus     stat;      /* status structure */
  INT2                 opt;       /* command-line option character */
  BOOLEAN              grace;     /* whether or not to graph using xmgrace */
  BOOLEAN              nonGrace;  /* whether or not to write to data file */
  TwoDMeshNode        *firstNode; /* head of linked list of nodes in mesh */
  static TwoDMeshParamStruc mesh; /* mesh parameters */
  static PtoleMetricIn search;    /* more mesh parameters */
  REAL4                mismatch;  /* mismatch threshold of mesh */
  UINT4                maxNodes;  /* maximum nodes in mesh */
  REAL4                begin;     /* start time of integration (seconds) */
  REAL4                duration;  /* duration of integration (seconds) */
  REAL4                fMax;      /* maximum frequency of search */
  FILE                *fp;        /* where to write a plot */



  /* Set default values. */
  grace = 0;
  nonGrace = 0;
  begin = 0.0;
  duration = 1e5;
  fMax = 1e3;
  mismatch = .02;
  maxNodes = 1e6;
  /* This is (roughly) the center of globular cluster 47 Tuc. */
  center.system = COORDINATESYSTEM_EQUATORIAL;
  center.longitude = (24.1/60)*LAL_PI_180;
  center.latitude = -(72+5./60)*LAL_PI_180;
  radius = 24.0/60*LAL_PI_180;

  /* Parse and sanity-check the command-line options. */
  while( (opt = getopt( argc, argv, "b:c:e:f:im:n:pr:t:x" )) != -1 )
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
      center.longitude = (15*a+b/4+c/240)*LAL_PI_180;
      center.latitude = (d+e/60+f/3600)*LAL_PI_180;
      break;
    case 'e':
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
    case 'p':
      nonGrace = 1;
      break;
    case 'r':
      radius = LAL_PI_180/60*atof( optarg );
      if( radius < 0 ) {
        fprintf( stderr, "%s line %d: %s\n", __FILE__, __LINE__,
                 PTOLEMESHTESTC_MSGERNG );
        return PTOLEMESHTESTC_ERNG;
      }
      break;
    case 't':
      duration = atof( optarg );
      if( duration < MIN_DURATION || duration > MAX_DURATION ) {
	fprintf( stderr, "%s line %d: %s\n", __FILE__, __LINE__,
                 PTOLEMESHTESTC_MSGERNG );
        return PTOLEMESHTESTC_ERNG;
      }
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
  mesh.getMetric = getMetric;
  mesh.metricParams = (void *) &search;
  if( radius == 0 )
    {
      mesh.domain[0] = dec_min;
      mesh.domain[1] = dec_max;
      mesh.rangeParams = (void *) &search;
    }
  else
    {
      mesh.domain[0] = center.latitude - radius;
      mesh.domain[1] = center.latitude + radius;
      mesh.rangeParams = NULL;

    }
  search.position.system = COORDINATESYSTEM_EQUATORIAL;
  search.spindown = NULL;
  search.epoch.gpsSeconds = begin;
  search.epoch.gpsNanoSeconds = 0;
  search.duration = duration;
  search.maxFreq = fMax;
  /* BEN: make this flexible later */
  search.site = &lalCachedDetectors[LALDetectorIndexGEO600DIFF];
  firstNode = NULL;
  LALCreateTwoDMesh( &stat, &firstNode, &mesh );
  if( stat.statusCode )
    return stat.statusCode;
  printf( "created %d nodes\n", mesh.nOut );

  /* Plot what we've got, if asked. */
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

  /* Write what we've got to file mesh.dat */
  if( nonGrace )
  {
    TwoDMeshNode *node;
    fp = fopen( "mesh.dat", "w" );
    if( !fp )
      return PTOLEMESHTESTC_EFIO;

    for( node = firstNode; node; node = node->next )
      fprintf( fp, "%e %e\n", node->y, node->x);
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
  INITSTATUS(stat);
  ATTATCHSTATUSPTR( stat );

  /* Search a circle. BEN: The 1.001 is a kludge. */
  y[0] = center.longitude - sqrt( pow( radius*1.001, 2 )
				  - pow( x-center.latitude, 2 ) );
  y[1] = center.longitude + sqrt( pow( radius*1.001, 2 )
				  - pow( x-center.latitude, 2 ) );

  if( unused )
    {
      y[0] = RA_min;
      y[1] = RA_max;
    }

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
  INITSTATUS(stat);
  ATTATCHSTATUSPTR( stat );
  TRY( LALDCreateVector( stat->statusPtr, &metric, 6 ), stat );

  /* Translate input. */
  patch->position.longitude = x[1];
  patch->position.latitude =  x[0];


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
      g[1] = metric->data[2];
      g[0] = metric->data[5];
      g[2] = metric->data[4];

  /* Clean up and leave. */
  TRY( LALDDestroyVector( stat->statusPtr, &metric ), stat );
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
} /* getMetric() */

/** \endcond */
