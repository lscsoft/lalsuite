/*
 * Copyright (C) 2007 Reinhard Prix
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
 * \author Reinhard Prix
 * \date 2006
 * \file
 * \brief Functions for optimal lattice-covering of Doppler-spaces.
 * NOTE: this module always uses *ECLIPTIC* coordinates for interal operations
 */

/*---------- INCLUDES ----------*/
#include <math.h>

#include <lal/LALStdlib.h>
#include <lal/DetectorSite.h>
#include <lal/LALError.h>
#include <lal/LatticeCovering.h>
#include <lal/LogPrintf.h>

#include "FlatPulsarMetric.h"

#include "DopplerFullScan.h"
#include "DopplerLatticeCovering.h"

/*---------- DEFINES ----------*/
NRCSID( DOPPLERLATTICECOVERING, "$Id$" );

#define TRUE (1==1)
#define FALSE (1==0)

#define MIN(x,y) (x < y ? x : y)
#define MAX(x,y) (x > y ? x : y)
#define SIGN(x) ( x < 0 ? -1 : ( x > 0 ? 1 : 0 ) )

#define SQUARE(x) ( (x) * (x) )
#define VECT_NORM(x) sqrt( SQUARE((x)[0]) + SQUARE((x)[1]) + SQUARE((x)[2]) )
#define VECT_ADD(x,y) do { (x)[0] += (y)[0]; (x)[1] += (y)[1]; (x)[2] += (y)[2]; } while(0)
#define VECT_MULT(x,k) do { (x)[0] *= k; (x)[1] *= k; (x)[2] *= k; } while(0)
#define VECT_COPY(dst,src) do { (dst)[0] = (src)[0]; (dst)[1] = (src)[1]; (dst)[2] = (src)[2]; } while(0)


/*---------- internal types ----------*/
typedef enum {
  HEMI_BOTH  =  0,	/**< points lie on both hemispheres */
  HEMI_NORTH =  1,	/**< all points on northern hemisphere */
  HEMI_SOUTH =  2	/**< all points on southern hemisphere */
} hemisphere_t;

typedef REAL8 vect2D_t[2];	/**< 2D vector */
typedef REAL8 vect3D_t[3];	/**< 3D vector */

/** 2D-polygon of points {nX, nY} on a single hemisphere of ecliptic sky-sphere */
typedef struct {
  UINT4 length;			/**< number of elements */
  vect2D_t *data;		/**< array of 2D vectors */
  hemisphere_t hemisphere;	/**< which sky-hemisphere the polygon lies on */
} hemiPolygon2D_t;

/** List of 3D vectors */
typedef struct {
  UINT4 length;			/**< number of elements */
  vect3D_t *data;		/**< array of 3D vectors */
} vect3Dlist_t;


struct tagDopplerLatticeScan {
  scan_state_t state;		/**< current state of the scan: idle, ready of finished */
  hemiPolygon2D_t skyRegionEcl; /**< (ecliptic) vector-polygon {nX, nY} defining a sky search-region */
  PulsarSpinRange wkRange;	/**< search-region in spin parameters ['canonical' units] */
  REAL8 Tspan;			/**< total observation time spanned */
  UINT4 dimSearch;		/**< dimension of search-space to cover [can be less than dim(latticeOrigin)!] */
  gsl_vector *latticeOrigin;	/**< 'origin' of the lattice {w0, kX, kY, w1, w2, ... } */
};

/*---------- empty initializers ---------- */

/*---------- Global variables ----------*/
extern INT4 lalDebugLevel;

/*---------- internal prototypes ----------*/
void skyRegionString2vect3D ( LALStatus *, vect3Dlist_t **skyRegionEcl, const CHAR *skyRegionString );
void setupSearchRegion ( LALStatus *status, DopplerLatticeScan *scan, const DopplerRegion *searchRegion );
hemisphere_t onWhichHemisphere ( const vect3Dlist_t *skypoints );
int skyposToVect3D ( vect3D_t *eclVect, const SkyPosition *skypos );
int vect3DToSkypos ( SkyPosition *skypos, const vect3D_t *vect );
int findCenterOfMass ( vect3D_t *center, const vect3Dlist_t *points );
int convertDoppler2Canonical ( gsl_vector **canonicalPoint, const vect3D_t vn, const PulsarSpins fkdot, REAL8 Tspan );
int convertSpins2Canonical ( PulsarSpins wk, const PulsarSpins fkdot, REAL8 Tspan );

/*==================== FUNCTION DEFINITIONS ====================*/

/** Initialized search-grid using optimal lattice-covering
 */
void
InitDopplerLatticeScan ( LALStatus *status, 
			 DopplerLatticeScan **scan, 		/**< [out] initialized scan-state for lattice-scan */
			 const DopplerLatticeInit *init		/**< [in] scan init parameters */
			 )
{
  DopplerLatticeScan *ret = NULL;
  gsl_matrix *gij;
  gsl_matrix *generator = NULL;

  INITSTATUS( status, "InitDopplerLatticeScan", DOPPLERLATTICECOVERING );
  ATTATCHSTATUSPTR ( status );

  ASSERT ( scan, status, DOPPLERSCANH_ENULL, DOPPLERSCANH_MSGENULL );  
  ASSERT ( *scan == NULL, status, DOPPLERSCANH_ENONULL, DOPPLERSCANH_MSGENONULL );  
  ASSERT ( init, status, DOPPLERSCANH_ENULL, DOPPLERSCANH_MSGENULL );  

  /* prepare scan-structure to return */
  if ( (ret = LALCalloc ( 1, sizeof(*ret) )) == NULL ) {
    ABORT (status, DOPPLERSCANH_EMEM, DOPPLERSCANH_MSGEMEM);
  }

  /* first store observation-time, required for conversion 'Doppler <--> canonical' */
  ret->Tspan = init->Tspan;

  /* ----- get search-Region ----- */
  TRY ( setupSearchRegion ( status->statusPtr, ret, &(init->searchRegion) ), status );

  /* ----- compute flat metric ----- */
  if ( (gij = gsl_matrix_calloc (ret->dimSearch, ret->dimSearch)) == NULL ) {
    ABORT (status, DOPPLERSCANH_EMEM, DOPPLERSCANH_MSGEMEM);
  }

  if ( XLALFlatMetricCW ( gij, init->searchRegion.refTime, init->startTime, init->Tspan, init->ephemeris ) )
    {
      LALPrintError ("\nCall to XLALFlatMetricCW() failed!\n\n");
      ABORT ( status, DOPPLERSCANH_EXLAL, DOPPLERSCANH_MSGEXLAL );
    }
  fprintf ( stderr, "gij = \\\n");
  XLALfprintfGSLmatrix ( stderr, "%.9g ", gij );

  /* ----- compute generating matrix for the lattice ----- */
  if ( XLALFindCoveringGenerator (&generator, LATTICE_TYPE_ANSTAR, sqrt(init->metricMismatch), gij ) < 0)
    {
      int code = xlalErrno;
      XLALClearErrno(); 
      LALPrintError ("\nERROR: XLALFindCoveringGenerator() failed (xlalErrno = %d)!\n\n", code);
      ABORT ( status, DOPPLERSCANH_EXLAL, DOPPLERSCANH_MSGEXLAL );
    }


  /* clean up memory */
  gsl_matrix_free ( gij );

  /* return final scan-state */
  ret->state = STATE_READY;
  (*scan) = ret;


  DETATCHSTATUSPTR ( status );
  RETURN ( status );

} /* InitDopplerLatticeScan() */



/** Free an allocated DopplerLatticeScan 
 */
void
FreeDopplerLatticeScan ( LALStatus *status, DopplerLatticeScan **scan )
{
  INITSTATUS( status, "FreeDopplerLatticeScan", DOPPLERLATTICECOVERING );

  ASSERT ( scan, status, DOPPLERSCANH_ENULL, DOPPLERSCANH_MSGENULL );  
  ASSERT ( *scan, status, DOPPLERSCANH_ENONULL, DOPPLERSCANH_MSGENONULL );  
  ASSERT ( (*scan)->state != STATE_IDLE, status, DOPPLERSCANH_EINPUT, DOPPLERSCANH_MSGEINPUT );

  if ( (*scan)->skyRegionEcl.data )
    LALFree ( (*scan)->skyRegionEcl.data );

  if ( (*scan)->latticeOrigin )
    gsl_vector_free ( (*scan)->latticeOrigin );

  /* FIXME: more freeing will be done here */


  
  LALFree ( (*scan) );
  (*scan) = NULL;

  RETURN(status);

} /* FreeDopplerLatticeScan() */

/** Translate the input 'DopplerRegion' into an internal representation for the scan.
 * 
 * \note DopplerLatticeScan->Tspan must have been set! This is used for conversion Doppler -> canonical 
 */
void
setupSearchRegion ( LALStatus *status, DopplerLatticeScan *scan, const DopplerRegion *searchRegion )
{
  UINT4 i;
  vect3Dlist_t *points3D = NULL;
  vect3D_t com = {0,0,0};
  PulsarSpins fkdotMid;
  UINT4 numSpins;

  INITSTATUS( status, "InitDopplerLatticeScan", DOPPLERLATTICECOVERING );
  ATTATCHSTATUSPTR ( status );

  ASSERT ( scan, status, DOPPLERSCANH_ENULL, DOPPLERSCANH_MSGENULL );  
  ASSERT ( searchRegion, status, DOPPLERSCANH_ENULL, DOPPLERSCANH_MSGENULL );  

  /* ----- sky ----- */
  TRY ( skyRegionString2vect3D ( status->statusPtr, &points3D, searchRegion->skyRegionString ), status );

  if ( (scan->skyRegionEcl.hemisphere = onWhichHemisphere ( points3D )) == HEMI_BOTH ) {
    LALPrintError ("\n\nSorry, currently only (ecliptic) single-hemisphere sky-regions are supported!\n");
    ABORT ( status, DOPPLERSCANH_EINPUT, DOPPLERSCANH_MSGEINPUT );
  }
  if ( (scan->skyRegionEcl.data = LALCalloc ( points3D->length, sizeof(scan->skyRegionEcl.data[0]) )) == NULL ) {
    ABORT (status, DOPPLERSCANH_EMEM, DOPPLERSCANH_MSGEMEM);
  }
  scan->skyRegionEcl.length = points3D->length;

  for ( i=0; i < points3D->length; i ++ ) {
    scan->skyRegionEcl.data[i][0] = points3D->data[i][0];
    scan->skyRegionEcl.data[i][1] = points3D->data[i][1];
  }

  findCenterOfMass ( &com, points3D );

  LALFree ( points3D->data );
  LALFree ( points3D );

  /* ----- spins ----- */
  scan->wkRange.refTime   = searchRegion->refTime;
  convertSpins2Canonical ( scan->wkRange.fkdot, searchRegion->fkdot, scan->Tspan );
  convertSpins2Canonical ( scan->wkRange.fkdotBand, searchRegion->fkdotBand, scan->Tspan );

  numSpins = sizeof(fkdotMid) / sizeof(fkdotMid[0]) ;	/* number of elements in PulsarSpin array */

  for ( i=0; i < numSpins; i ++ )
    fkdotMid[i] = searchRegion->fkdot[i] + 0.5 * searchRegion->fkdotBand[i];

  /* ----- use the center of the searchRegion as origin of the lattice ----- */
  if ( 0 != convertDoppler2Canonical ( &(scan->latticeOrigin), com, fkdotMid, scan->Tspan ) ) {
    LALPrintError ("\n\nconvertDoppler2Canonical() failed!\n");
    ABORT ( status, DOPPLERSCANH_EXLAL, DOPPLERSCANH_MSGEXLAL );
  }

  /* determine number of spins to compute metric for (at least 1) */
  while ( (numSpins > 1) && (searchRegion->fkdotBand[numSpins - 1] == 0) )
    numSpins --;

  scan->dimSearch = 2 + numSpins;	/* sky + spins (must be at least 3) */

  DETATCHSTATUSPTR ( status );
  RETURN ( status );

} /* setupSearchRegion() */


/** Convert Doppler-parameters from {nX, nY, nZ, fkdot} into internal 'canonical' form
 *  {w0, kX, kY, w1, w1, ... } 
 * \note The return-vector is allocated here and needs to be gsl_vector_free()'ed
 * 
 * Return: 0=OK, -1=ERROR
 */
int
convertDoppler2Canonical ( gsl_vector **canonicalPoint, const vect3D_t vn, const PulsarSpins fkdot, REAL8 Tspan )
{
  REAL8 kX, kY;
  REAL8 prefix;
  UINT4 s;
  PulsarSpins wk;
  gsl_vector *ret;
  UINT4 numSpins = sizeof(wk) / sizeof(wk[0]);

  if ( !canonicalPoint || (*canonicalPoint != NULL) || !vn )
    return -1;

  prefix = (LAL_TWOPI * LAL_AU_SI / LAL_C_SI ) * fkdot[0];
  kX = -prefix * vn[0];		/* vk = - 2*pi * Rorb/c * Freq * vn */
  kY = -prefix * vn[1];
  
  convertSpins2Canonical ( wk, fkdot, Tspan );
  
  if ( (ret = gsl_vector_calloc ( 2 + numSpins )) == NULL ) {
    LALPrintError("\n\nOut of memory!!\n");
    return -1;
  }
  gsl_vector_set (ret, 0, wk[0]);
  gsl_vector_set (ret, 1, kX );
  gsl_vector_set (ret, 2, kY );
  for ( s=1; s < numSpins; s ++ )
    gsl_vector_set (ret, 2 + s, wk[s] );

  /* return vector */
  (*canonicalPoint) = ret;

  return 0;

} /* convertDoppler2Canonical() */

/** Convert SI-unit spins 'fkdot' into canonical units w^(s) = 2*pi * f^(s) * T^(s+1) 
 */
int
convertSpins2Canonical ( PulsarSpins wk, const PulsarSpins fkdot, REAL8 Tspan )
{
  PulsarSpins dummy;
  UINT4 numSpins = sizeof(dummy) / sizeof(dummy[0]);
  REAL8 prefix = LAL_TWOPI * Tspan;
  UINT4 s;

  for ( s=0; s < numSpins; s ++ )
    {
      wk[s] = prefix * fkdot[s];	/* wk = 2*pi * T^(s+1) * fkdot */
      prefix *= Tspan;
    }

  return 0;

} /* convertSpins2Canonical() */


/** Convert a sky-region string into a list of vectors in ecliptic 2D coords {nX, nY} 
 */
void
skyRegionString2vect3D ( LALStatus *status, 
			 vect3Dlist_t **skyRegionEcl, 	/**< [out] list of skypoints in 3D-ecliptic coordinates {nX, nY, nZ} */
			 const CHAR *skyRegionString )	/**< [in] string of equatorial-coord. sky-positions */
{
  vect3Dlist_t *ret = NULL;
  SkyRegion region = empty_SkyRegion;
  UINT4 i;

  INITSTATUS( status, "skyRegionString2list", DOPPLERLATTICECOVERING );
  ATTATCHSTATUSPTR ( status );

  ASSERT ( skyRegionEcl, status, DOPPLERSCANH_ENULL, DOPPLERSCANH_MSGENULL );  
  ASSERT ( *skyRegionEcl == NULL, status, DOPPLERSCANH_ENONULL, DOPPLERSCANH_MSGENONULL );  
  ASSERT ( skyRegionString, status, DOPPLERSCANH_ENULL, DOPPLERSCANH_MSGENULL );  

  TRY ( ParseSkyRegionString (status->statusPtr, &region, skyRegionString), status );

  if ( ( ret = LALCalloc ( 1, sizeof(*ret)) ) == NULL ) {
    ABORT (status, DOPPLERSCANH_EMEM, DOPPLERSCANH_MSGEMEM);    
  }
  if ( (ret->data = LALCalloc ( region.numVertices, sizeof( ret->data[0] ) )) == NULL ) {
    LALFree ( ret );
    ABORT (status, DOPPLERSCANH_EMEM, DOPPLERSCANH_MSGEMEM);    
  }
  ret->length = region.numVertices;
  
  for ( i=0; i < region.numVertices; i ++ )
    skyposToVect3D ( &(ret->data[i]), &(region.vertices[i]) );

  /* return sky-region */
  (*skyRegionEcl) = ret;

  DETATCHSTATUSPTR ( status );
  RETURN ( status );

} /* skyRegionString2list() */

/** Check whether given list of skypoint lie on a single hemisphere 
 */
hemisphere_t
onWhichHemisphere ( const vect3Dlist_t *skypoints )
{
  UINT4 i;
  INT4 our_sign = 0;

  if ( !skypoints || (skypoints->length == 0) )
    return FALSE;

  for ( i=0; i < skypoints->length; i ++ )
    {
      vect3D_t *thisPoint = &(skypoints->data[i]);
      INT4 this_sign = SIGN ( (*thisPoint)[2] );
      if ( this_sign && (our_sign == 0) )	/* set our_sign to first non-zero sign */
	our_sign = this_sign;
      if ( this_sign && (this_sign != our_sign ) )
	return HEMI_BOTH;
    }
  
  if ( our_sign < 0 )
    return HEMI_SOUTH;
  else if ( our_sign > 0 )
    return HEMI_NORTH;
  else
    return HEMI_BOTH;		/* all points on equator.. */
  
} /* onWhichHemisphere() */

/** Convert a 'SkyPosition' {alpha, delta} into a 3D unit vector in ecliptic coordinates 
 * Return: 0=OK, -1=ERROR
 */
int 
skyposToVect3D ( vect3D_t *eclVect, const SkyPosition *skypos )
{
  REAL8 cosa, cosd, sina, sind, sineps, coseps;
  REAL8 nn[3];	/* unit-vector pointing to source */

  if ( !eclVect || !skypos )
    return -1;

  sina = sin(skypos->longitude);
  cosa = cos(skypos->longitude);
  sind = sin(skypos->latitude);
  cosd = cos(skypos->latitude);

  nn[0] = cosa * cosd;
  nn[1] = sina * cosd;
  nn[2] = sind;

  switch ( skypos->system )
    {
    case COORDINATESYSTEM_EQUATORIAL:
      sineps = SIN_EPS; coseps = COS_EPS;
      break;
    case COORDINATESYSTEM_ECLIPTIC:
      sineps = 0; coseps = 1;
      break;
    default:
      XLAL_ERROR ( "skyposToVect3D", XLAL_EINVAL );
      break;
    } /* switch(system) */

  (*eclVect)[0] =   nn[0];
  (*eclVect)[1] =   nn[1] * coseps + nn[2] * sineps;
  (*eclVect)[2] = - nn[1] * sineps + nn[2] * coseps;

  return 0;

} /* skyposToVect3D() */

/** Convert (ecliptic) vector pointing to a skypos back into (alpha, delta) 
 * The coordinate-system of the result is determined by the 'system'-setting in
 * the output skypos
 * return: 0=OK, -1=ERROR
 */
int
vect3DToSkypos ( SkyPosition *skypos, const vect3D_t *vect )
{
  REAL8 invnorm;
  vect3D_t nvect = {0,0,0};
  REAL8 longitude, latitude;
  REAL8 sineps, coseps;

  if ( !skypos || !vect )
    return -1;

  switch ( skypos->system )
    {
    case COORDINATESYSTEM_EQUATORIAL:
      sineps = SIN_EPS; coseps = COS_EPS;
      break;
    case COORDINATESYSTEM_ECLIPTIC:
      sineps = 0; coseps = 1;
      break;
    default:
      XLAL_ERROR ( "vect3DToSkypos", XLAL_EINVAL );
      break;
    } /* switch(system) */

  nvect[0] = (*vect)[0];
  nvect[1] = coseps * (*vect)[1] - sineps * (*vect)[2];
  nvect[2] = sineps * (*vect)[1] + coseps * (*vect)[2];

  invnorm = 1.0 / VECT_NORM ( nvect );
  
  VECT_MULT ( nvect, invnorm );

  longitude = atan2 ( nvect[1], nvect[0]);
  if ( longitude < 0 )
    longitude += LAL_TWOPI;

  latitude = asin ( nvect[2] );

  skypos->longitude = longitude;
  skypos->latitude = latitude;

  return 0;

} /* XLALSkyVector2Skypos() */


/** Find the "center of mass" of the given list of 3D points
 * Return: 0=OK, -1 on ERROR
 */
int
findCenterOfMass ( vect3D_t *center, const vect3Dlist_t *points )
{
  UINT4 i;
  vect3D_t com = {0, 0, 0};	/* center of mass */
  REAL8 norm;

  if ( !center || !points )
    return -1;

  for ( i=0; i < points->length ; i ++ )
    {
      VECT_ADD ( com, points->data[i] );
    }
  norm = 1.0 / points->length;
  VECT_MULT ( com, norm );

  VECT_COPY ( (*center), com );

  return 0;

} /* findCenterOfmass() */



