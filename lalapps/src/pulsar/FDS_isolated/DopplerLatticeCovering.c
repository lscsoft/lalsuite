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

#include "DopplerFullScan.h"
#include "FlatPulsarMetric.h"

/*---------- DEFINES ----------*/
NRCSID( DOPPLERLATTICECOVERING, "$Id$" );

#define TRUE (1==1)
#define FALSE (1==0)

#define MIN(x,y) (x < y ? x : y)
#define MAX(x,y) (x > y ? x : y)

#define SQUARE(x) ( (x) * (x) )
#define VECT_NORM(x) sqrt( SQUARE((x)[0]) + SQUARE((x)[1]) + SQUARE((x)[2]) )
#define VECT_ADD(x,y) do { (x)[0] += (y)[0]; (x)[1] += (y)[1]; (x)[2] += (y)[2]; } while(0)
#define VECT_MULT(x,k) do { (x)[0] *= k; (x)[1] *= k; (x)[2] *= k; } while(0)
#define VECT_COPY(dst,src) do { (dst)[0] = (src)[0]; (dst)[1] = (src)[1]; (dst)[2] = (src)[2]; } while(0)


/*---------- internal types ----------*/
typedef REAL8 vect2D_t[2];	/**< 2D vector */

/** List of 2D vectors */
typedef struct {
  UINT4 length;		/**< number of elements */
  vect2D_t *data;		/**< array of 2D vectors */
} vect2Dlist_t;

struct tagDopplerLatticeScan {
  scan_state_t state;		/**< current state of the scan: idle, ready of finished */

};

/* FIXME: remove */
typedef REAL8 SkyVector[3];
typedef struct {
  UINT4 length;			/**< number of list elements */
  SkyVector *data;		/**< array of sky-vectors */
} SkyVectorList;

int XLALSkypos2SkyVector ( SkyVector *vect, REAL8 longitude, REAL8 latitude );
int XLALSkyVector2Skypos ( REAL8 *longitude, REAL8 *latitude, const SkyVector *vect );
void FindSkyRegionCenter ( LALStatus *, SkyVector *center, const SkyRegion *region );


/*---------- empty initializers ---------- */

/*---------- Global variables ----------*/
extern INT4 lalDebugLevel;

/*---------- internal prototypes ----------*/
BOOLEAN searchBoundary ( const REAL8Vector *point );
void skyRegionString2list ( LALStatus *, vect2Dlist_t **skyRegionEcl, const CHAR *skyRegionString );
int XLALSkyPosition2EclVect ( vect2D_t *eclVect, const SkyPosition *skypos );

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
  vect2Dlist_t *skyRegionEcl = NULL;

  INITSTATUS( status, "InitDopplerLatticeScan", DOPPLERLATTICECOVERING );
  ATTATCHSTATUSPTR ( status );

  ASSERT ( scan, status, DOPPLERSCANH_ENULL, DOPPLERSCANH_MSGENULL );  
  ASSERT ( *scan == NULL, status, DOPPLERSCANH_ENONULL, DOPPLERSCANH_MSGENONULL );  
  ASSERT ( init, status, DOPPLERSCANH_ENULL, DOPPLERSCANH_MSGENULL );  

  if ( (ret = LALCalloc ( 1, sizeof(*ret) )) == NULL ) {
    ABORT (status, DOPPLERSCANH_EMEM, DOPPLERSCANH_MSGEMEM);
  }

  TRY ( skyRegionString2list ( status->statusPtr, &skyRegionEcl, init->searchRegion.skyRegionString ), status );

  

  DETATCHSTATUSPTR ( status );
  RETURN ( status );

} /* InitDopplerLatticeScan() */

/** Convert a sky-region string into a list of vectors in ecliptic 2D coords {nX, nY} 
 */
void
skyRegionString2list ( LALStatus *status, 
		       vect2Dlist_t **skyRegionEcl, 	/**< [out] list of skypoints in 2D-ecliptic coordinates {nX, nY} */
		       const CHAR *skyRegionString )	/**< [in] string of equatorial-coord. sky-positions */
{
  vect2Dlist_t *ret = NULL;
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
    XLALSkyPosition2EclVect ( &(ret->data[i]), &(region.vertices[i]) );

  DETATCHSTATUSPTR ( status );
  RETURN ( status );

} /* skyRegionString2list() */

int 
XLALSkyPosition2EclVect ( vect2D_t *eclVect, const SkyPosition *skypos )
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
    case   COORDINATESYSTEM_EQUATORIAL:
      sineps = SIN_EPS; coseps = COS_EPS;
      break;
    case COORDINATESYSTEM_ECLIPTIC:
      sineps = 0; coseps = 1;
      break;
    default:
      XLAL_ERROR ( "XLALSkyPosition2EclVect", XLAL_EINVAL );
      break;
    } /* switch(system) */

  (*eclVect)[0] = nn[0];
  (*eclVect)[1] = nn[1] * coseps + nn[2] * sineps;

  return 0;

} /* XLALSkypos2SkyVector() */


/** Free an allocated DopplerLatticeScan 
 */
void
FreeDopplerLatticeScan ( LALStatus *status, DopplerLatticeScan **scan )
{
  INITSTATUS( status, "FreeDopplerLatticeScan", DOPPLERLATTICECOVERING );

  ASSERT ( scan, status, DOPPLERSCANH_ENULL, DOPPLERSCANH_MSGENULL );  
  ASSERT ( *scan, status, DOPPLERSCANH_ENONULL, DOPPLERSCANH_MSGENONULL );  
  ASSERT ( (*scan)->state != STATE_IDLE, status, DOPPLERSCANH_EINPUT, DOPPLERSCANH_MSGEINPUT );

  /* FIXME: more freeing will be done here */
  
  LALFree ( (*scan) );
  (*scan) = NULL;

  RETURN(status);
} /* FreeDopplerLatticeScan() */


/** Initialized and construct an optimal lattice-covering for the given searchRegion.
 */
void
initLatticeCovering ( LALStatus *status, 
		      DopplerLatticeScan *scan, 
		      const DopplerFullScanInit *init)
{
  UINT4 numSpins, dim;
  int ret;
  REAL8VectorList *pts = NULL;
  REAL8Vector *p0 = NULL;
  gsl_matrix *gij;

  INITSTATUS( status, "initLatticeCovering", DOPPLERLATTICECOVERING );
  ATTATCHSTATUSPTR ( status );

  ASSERT ( scan, status, DOPPLERSCANH_EINPUT, DOPPLERSCANH_MSGEINPUT );
  ASSERT ( scan->state == STATE_IDLE, status, DOPPLERSCANH_EINPUT, DOPPLERSCANH_MSGEINPUT );
  ASSERT ( init, status, DOPPLERSCANH_EINPUT, DOPPLERSCANH_MSGEINPUT );

  /* determine number of spins to compute metric for (at least 1) */
  numSpins = PULSAR_MAX_SPINS;
  while ( (numSpins > 1) && (init->searchRegion.fkdotBand[numSpins - 1] == 0) )
    numSpins --;

  dim = 2 + numSpins;	/* sky + spins (must be at least 3) */

  /* compute flat metric */
  if ( (gij = gsl_matrix_calloc (dim, dim)) == NULL ) {
    ABORT (status, DOPPLERSCANH_EMEM, DOPPLERSCANH_MSGEMEM);
  }

  ret = XLALFlatMetricCW ( gij, init->searchRegion.refTime, init->startTime, init->Tspan, init->ephemeris );
  if ( ret != 0 )
    {
      LALPrintError ("\nCall to XLALFlatMetricCW() failed!\n\n");
      ABORT ( status, DOPPLERSCANH_EXLAL, DOPPLERSCANH_MSGEXLAL );
    }

  fprintf ( stderr, "gij = \\\n");
  XLALfprintfGSLmatrix ( stderr, "%.9g ", gij );

  if ( ( p0 = XLALCreateREAL8Vector ( dim ) ) == NULL ) {
    ABORT (status, DOPPLERSCANH_EMEM, DOPPLERSCANH_MSGEMEM);
  }

  {
    SkyVector skyCenter = {0,0,0};
    SkyRegion region;
    UINT4 s;
    TRY ( ParseSkyRegionString (status->statusPtr, &region, init->searchRegion.skyRegionString), status );

    TRY ( FindSkyRegionCenter ( status->statusPtr, &skyCenter, &region ), status );

    p0->data[0] = init->searchRegion.fkdot[0] + 0.5 * init->searchRegion.fkdotBand[0];
    p0->data[1] = 0;	/* FIXME */
    p0->data[2] = 0;
    for (s=1; s < numSpins; s ++ )
      p0->data[2 + s] = init->searchRegion.fkdot[s] + 0.5 * init->searchRegion.fkdotBand[s];
  }

  TRY ( LALLatticeCovering (status->statusPtr, &pts, sqrt( init->metricMismatch ), gij, p0, searchBoundary, LATTICE_TYPE_ANSTAR), status);

  XLALDestroyREAL8Vector ( p0 );

  DETATCHSTATUSPTR ( status );
  RETURN ( status );

} /* initLatticeCovering() */

/** implement the boundary of the search region 
 */
BOOLEAN 
searchBoundary ( const REAL8Vector *point )
{
  return TRUE;
} /* searchBoundary() */

/** Find the central point of the given sky-region: simply compute the center-of-mass
 */
void
FindSkyRegionCenter ( LALStatus *status, SkyVector *center, const SkyRegion *region )
{
  UINT4 i;
  SkyVector com = {0, 0, 0};	/* center of mass */
  REAL8 norm;

  INITSTATUS( status, "FindSkyRegionCenter", DOPPLERLATTICECOVERING ); 

  ASSERT ( center, status, DOPPLERSCANH_ENULL, DOPPLERSCANH_MSGENULL );  
  ASSERT ( region, status, DOPPLERSCANH_ENULL, DOPPLERSCANH_MSGENULL );  

  for ( i=0; i < region->numVertices ; i ++ )
    {
      SkyPosition *thisPos = &(region->vertices[i]);
      SkyVector vertex;
      XLALSkypos2SkyVector ( &vertex, thisPos->longitude, thisPos->latitude );
      VECT_ADD ( com, vertex );
    }
  norm = 1.0 / region->numVertices;
  VECT_MULT ( com, norm );

  VECT_COPY ( (*center), com );

  RETURN(status);

} /* FindSkyRegionCenter() */


/** Convert sky-pos (alpha, delta) into unit-vector pointing to skypos 
 */
int 
XLALSkypos2SkyVector ( SkyVector *vect, REAL8 longitude, REAL8 latitude )
{
  REAL8 cosa, cosd, sina, sind;

  if ( !vect )
    return -1;

  sina = sin(longitude);
  cosa = cos(longitude);
  sind = sin(latitude);
  cosd = cos(latitude);

  (*vect)[0] = cosa * cosd;
  (*vect)[1] = sina * cosd;
  (*vect)[2] = sind;

  return 0;

} /* XLALSkypos2SkyVector() */

/** Convert vector pointing to a skypos into (alpha, delta) 
 */
int
XLALSkyVector2Skypos ( REAL8 *longitude, REAL8 *latitude, const SkyVector *vect )
{
  REAL8 norm;
  SkyVector nvect;
  REAL8 alpha, delta;

  if ( !longitude || !latitude || !vect )
    return -1;

  norm = VECT_NORM ( (*vect) );

  nvect[0] = (*vect)[0] / norm;
  nvect[1] = (*vect)[1] / norm;
  nvect[2] = (*vect)[2] / norm;

  alpha = atan2 ( nvect[1], nvect[0]);
  if ( alpha < 0 )
    alpha += LAL_TWOPI;

  delta = asin ( nvect[2] );

  (*longitude) = alpha;
  (*latitude) = delta;

  return 0;

} /* XLALSkyVector2Skypos() */
