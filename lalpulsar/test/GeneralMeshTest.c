/*
*  Copyright (C) 2007 Jolien Creighton, Ian Jones, Benjamin Owen
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
 * \author Jones, D. I.,   Owen, B. J.
 * \file
 * \ingroup PtoleMetric_h
 * \brief Tests and showcases the combination of a LAL metric function (of the
 * user's specification) and \ref TwoDMesh_h modules by producing a
 * template mesh.
 *
 * ### Program <tt>GeneralMeshTest.c</tt> ###
 *
 *
 * ### Usage ###
 *
 * \code
 * GeneralMeshTest
 * \endcode
 *
 * ### Description ###
 *
 * The <b>-b</b> option sets the beginning GPS time of integration to
 * the option argument. (Default is \f$731265908\f$ seconds, chosen to lie
 * within the S2 run).
 *
 * The <b>-c</b> option determins the center of the patch.  This option
 * is hardcoded to use equatorial coordinates and the argument should be
 * given in hh:mm:ss:dd:mm:ss format.  (Default is the center of the
 * globular cluster 47 Tuc).
 *
 * The <b>-f</b> option sets the maximum frequency of integration (in Hz) to the
 * option argument. (The default value is 1000.)
 *
 * The <b>-l</b> option determines, for a rectangular search region, the
 * limits in right ascension and declination of the grid.  The argument should
 * be given in degrees as RA(min):RA(max):dec(min):dec(max).  (The default is
 * the octant of the sky defined by \f$0 < \textrm{RA} < 90\f$ and \f$0< \textrm{dec} <
 * 85\f$; this avoids the coordinate singularity at the poles.) This option
 * automatically overrides whatever is specified by the <b>-r</b> option.
 *
 * The <b>-m</b> option sets the maximum mismatch of the mesh to the option
 * argument. (Default is 0.02.)
 *
 * The texttt{-p} option causes the coordinates of the nodes to be written to
 * a file <tt>mesh.dat</tt>, for the benifit of users who don't have
 * \c xmgrace installed.  The format is one node per line, (RA, DEC),
 * with the angles in degrees.
 *
 * The <b>-r</b> option sets the radius (in arcminutes) of the circular
 * sky patch.  If you specify radius zero you will get a search over a
 * rectangular region whose limits in RA and dec are specified by the
 * <b>-l</b> option.  (The default value is the radius of the globular
 * cluster 47 Tuc).
 *
 * The <b>-t</b> option sets the duration of integration in seconds. (The
 * default is \f$39600\f$ seconds \f$= 11\f$ hours, which is chosen because it is of
 * the right size for S2 analyses).
 *
 * ### Algorithm ###
 *
 *
 * ### Uses ###
 *
 * \code
 * lalDebugLevel                LALDCreateVector()
 * LALCheckMemoryLeaks()        LALDDestroyVector()
 * LALProjectMetric()           XLALGetEarthTimes()
 * LALPtoleMetric()             XLALInitBarycenter()
 * LALCreateTwoDMesh()          LALDestroyTwoDMesh()
 * \endcode
 *
 * ### Notes ###
 *
 * For most regions of parameter space the three metric codes seem to
 * agree well.  However, for short (less than one day) runs, they are all
 * capable of returning (unphysical) negative determinant metrics for
 * points very close to the equator.
 *
 */


/** \name Error Codes */ /*@{*/
#define GENERALMESHTESTC_EMEM 1
#define GENERALMESHTESTC_ERNG 2
#define GENERALMESHTESTC_EFIO 3
#define GENERALMESHTESTC_EOPT 4
#define GENERALMESHTESTC_EMET 5

#define GENERALMESHTESTC_MSGEMEM "memory (de)allocation error"
#define GENERALMESHTESTC_MSGERNG "value out of range"
#define GENERALMESHTESTC_MSGEFIO "file I/O error"
#define GENERALMESHTESTC_MSGEOPT "unknown command-line option"
#define GENERALMESHTESTC_MSGEMET "determinant of projected metric negative"
/*@}*/

#include <config.h>

#include <math.h>
#include <stdio.h>
#include <lal/AVFactories.h>
#include <lal/LALgetopt.h>
#include <lal/PtoleMetric.h>
#include <lal/TwoDMesh.h>
#include <lal/LALInitBarycenter.h>


/** \cond DONT_DOXYGEN */

#define MIN_DURATION (86400./LAL_TWOPI) /* one radian of rotation */
#define MAX_DURATION (3.16e7)           /* one year; arbitrary */
#define MIN_FREQ     1e2                /* more or less arbitrary */
#define MAX_FREQ     1e4                /* more or less arbitrary */
#define MAX_NODES    1e6                /* keep idiot users safe */

/* IAN: Clumsy way of specifying rectangular search region if */
/* the r=0 option is invoked.                                 */
REAL4 ra_min;
REAL4 ra_max;
REAL4 dec_min;
REAL4 dec_max;

void getRange( LALStatus *, REAL4 [2], REAL4, void * );
void getMetric( LALStatus *, REAL4 [3], REAL4 [2], void * );

/* BEN: This is a cheat. Should make a search area structure. */
static SkyPosition center;  /* center of search */
REAL4              radius;  /* radius of search, in arcminutes */

int main( int argc, char **argv )
{
  static LALStatus     stat;      /* status structure */
  INT2                 opt;       /* command-line option character */
  BOOLEAN              nonGrace;  /* whether or not to write to data file */
  TwoDMeshNode        *firstNode; /* head of linked list of nodes in mesh */
  static TwoDMeshParamStruc mesh; /* mesh parameters */
  static PtoleMetricIn search;    /* more mesh parameters */
  REAL4                mismatch;  /* mismatch threshold of mesh */
  int                  begin;     /* start time of integration (seconds) */
  REAL4                duration;  /* duration of integration (seconds) */
  REAL4                fMax;      /* maximum frequency of search */
  FILE                *fp;        /* where to write a plot */
  float a, b, c, d, e, f;         /* To specify center of search region */
  BOOLEAN rectangular;            /* is the search region rectangular? */

  /* Set default values. */
  nonGrace = 0;
  begin = 731265908;
  duration = 1e5;
  fMax = 1e3;
  mismatch = .02;
  /* This is (roughly) the center of globular cluster 47 Tuc. */
  center.system = COORDINATESYSTEM_EQUATORIAL;
  center.longitude = (24.1/60)*LAL_PI_180;
  center.latitude = -(72+5./60)*LAL_PI_180;
  radius = 24.0/60*LAL_PI_180;
  ra_min = 0.0;
  ra_max = LAL_PI_2;
  dec_min = 0.0;
  dec_max = LAL_PI_2;
  rectangular = 0;

  /* Parse and sanity-check the command-line options. */
  while( (opt = LALgetopt( argc, argv, "a:b:c:d:ef:l:m:pr:t:x" )) != -1 )
  {
    switch( opt )
    {
    case '?':
      return GENERALMESHTESTC_EOPT;
    case 'b':
      begin = atoi( LALoptarg );
      break;
    case 'c':
      if( sscanf( LALoptarg, "%f:%f:%f:%f:%f:%f", &a, &b, &c, &d, &e, &f ) != 6)
      {
        fprintf( stderr, "coordinates should be hh:mm:ss:dd:mm:ss\n" );
        return GENERALMESHTESTC_EOPT;
      }
      center.longitude = (15*a+b/4+c/240)*LAL_PI_180;
      center.latitude = (d+e/60+f/3600)*LAL_PI_180;
      break;
    case 'f':
      fMax = atof( LALoptarg );
      break;
    case 'l':
      if( sscanf( LALoptarg, "%f:%f:%f:%f",
		  &ra_min, &ra_max, &dec_min, &dec_max) != 4)
	{
	  fprintf( stderr, "coordinates should be ra_min, ra_max, dec_min, dec_max, all in degrees\n" );
          return GENERALMESHTESTC_EOPT;
	}
      ra_min  *= LAL_PI_180;
      ra_max  *= LAL_PI_180;
      dec_min *= LAL_PI_180;
      dec_max *= LAL_PI_180;
      rectangular = 1;
      radius = 0;
      break;
    case 'm':
      mismatch = atof( LALoptarg );
      break;
    case 'p':
      nonGrace = 1;
      break;
    case 'r':
      if( rectangular == 1 )
        break;
      radius = LAL_PI_180/60*atof( LALoptarg );
      if( radius < 0 ) {
        fprintf( stderr, "%s line %d: %s\n", __FILE__, __LINE__,
                 GENERALMESHTESTC_MSGERNG );
        return GENERALMESHTESTC_ERNG;
      }
      break;
    case 't':
      duration = atof( LALoptarg );
      if( duration < MIN_DURATION || duration > MAX_DURATION ) {
	fprintf( stderr, "%s line %d: %s\n", __FILE__, __LINE__,
                 GENERALMESHTESTC_MSGERNG );
        return GENERALMESHTESTC_ERNG;
      }
      break;
    } /* switch( opt ) */
  } /* while( LALgetopt... ) */


  /* Set the mesh input parameters. */
  mesh.mThresh = mismatch;
  mesh.nIn = MAX_NODES;
  mesh.getRange = getRange;
  mesh.getMetric = getMetric;
  mesh.metricParams = (void *) &search;
  if( radius == 0 )
    {
      mesh.domain[0] = dec_min;
      mesh.domain[1] = dec_max;
      /* I don't understand this line*/
      mesh.rangeParams = (void *) &search;
    }
  else
    {
      mesh.domain[0] = center.latitude - radius;
      mesh.domain[1] = center.latitude + radius;
      mesh.rangeParams = NULL;
    }


  /* Set detector location */
  search.site = &lalCachedDetectors[LALDetectorIndexGEO600DIFF];

  /* Ptolemetric constants */
  search.position.system = COORDINATESYSTEM_EQUATORIAL;
  search.spindown = NULL;
  search.epoch.gpsSeconds = begin;
  search.epoch.gpsNanoSeconds = 0;
  search.duration = duration;
  search.maxFreq = fMax;

  /* Create mesh */
  firstNode = NULL;
  LALCreateTwoDMesh( &stat, &firstNode, &mesh );
  if( stat.statusCode )
    return stat.statusCode;
  printf( "created %d nodes\n", mesh.nOut );
  if( mesh.nOut == MAX_NODES )
    printf( "This overflowed your limit. Try a smaller search.\n" );


  /* Write what we've got to file mesh.dat, if required */
  if( nonGrace )
  {
    TwoDMeshNode *node;
    fp = fopen( "mesh.dat", "w" );
    if( !fp )
      return GENERALMESHTESTC_EFIO;

    for( node = firstNode; node; node = node->next )
      fprintf( fp, "%e %e\n",
	       (double)((node->y)*180/LAL_PI), (double)((node->x)*180/LAL_PI));
    fclose( fp );
  }

  /* Clean up and leave. */
  LALDestroyTwoDMesh( &stat, &firstNode, &mesh.nOut );
  printf( "destroyed %d nodes\n", mesh.nOut );
  if( stat.statusCode )
    return GENERALMESHTESTC_EMEM;

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
      y[0] = ra_min;
      y[1] = ra_max;
    }

  /*
  printf("x = %le,  y[0] = %le,  y[1] = %le\n", (double)x, (double)y[0],
	 (double)y[1]);
  */

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

  REAL8Vector      *metric = NULL;   /* for output of metric */
  PtoleMetricIn    *Ppatch = NULL;   /* pointer for PtoleMetric params */
  REAL8             determinant;     /* Determinant of projected metric */


  Ppatch = params;

  /* Set up shop. */
  INITSTATUS(stat);
  ATTATCHSTATUSPTR( stat );
  TRY( LALDCreateVector( stat->statusPtr, &metric, 6 ), stat );

  /* Translate input. */
  Ppatch->position.longitude = x[1];
  Ppatch->position.latitude =  x[0];

  /* Call the real metric function. */
  LALPtoleMetric( stat->statusPtr, metric, Ppatch );

  BEGINFAIL( stat )
    TRY( LALDDestroyVector( stat->statusPtr, &metric ), stat );
  ENDFAIL( stat );
  LALProjectMetric( stat->statusPtr, metric, 0 );
  BEGINFAIL( stat )
    TRY( LALDDestroyVector( stat->statusPtr, &metric ), stat );
  ENDFAIL( stat );

  determinant = metric->data[5]*metric->data[2]-pow(metric->data[4],2.0);
  if(determinant < 0.0)
    {
      printf( "%s line %d: %s\n", __FILE__, __LINE__,
	      GENERALMESHTESTC_MSGEMET );
      return;
    }



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
