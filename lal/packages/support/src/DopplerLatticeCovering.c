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
#include <gsl/gsl_sys.h>

#include <lal/LALStdlib.h>
#include <lal/DetectorSite.h>
#include <lal/LALError.h>
#include <lal/LatticeCovering.h>
#include <lal/LogPrintf.h>

#include <lal/FlatPulsarMetric.h>
#include <lal/DopplerFullScan.h>
#include <lal/DopplerLatticeCovering.h>

/*---------- DEFINES ----------*/
NRCSID( DOPPLERLATTICECOVERING, "$Id$" );

/* turn off  gsl range-checking in non-debug compilations */
#ifdef LAL_NDEBUG
#define GSL_RANGE_CHECK_OFF 1
#endif

#define EPS_REAL8	1e-10	/**< relative error used in REAL8 comparisons */

#define TRUE (1==1)
#define FALSE (1==0)

#define NUM_SPINS	PULSAR_MAX_SPINS
#define DIM_CANONICAL	(2 + NUM_SPINS)		/**< dimension of 'canonical' Doppler-space */

#define MIN(x,y) (x < y ? x : y)
#define MAX(x,y) (x > y ? x : y)
#define SIGN(x) ( x < 0 ? -1 : ( x > 0 ? 1 : 0 ) )

#define SQUARE(x) ( (x) * (x) )
#define VECT_NORM(x) sqrt( SQUARE((x)[0]) + SQUARE((x)[1]) + SQUARE((x)[2]) )
#define VECT_ADD(x,y) do { (x)[0] += (y)[0]; (x)[1] += (y)[1]; (x)[2] += (y)[2]; } while(0)
#define VECT_MULT(x,k) do { (x)[0] *= k; (x)[1] *= k; (x)[2] *= k; } while(0)
#define VECT_COPY(dst,src) do { (dst)[0] = (src)[0]; (dst)[1] = (src)[1]; (dst)[2] = (src)[2]; } while(0)

#define VECT_HEMI(x) ( (x)[2] < 0 ? HEMI_SOUTH : ( (x)[2] > 0 ? HEMI_NORTH : 0 ) )

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
} vect2Dlist_t;

/** List of 3D vectors */
typedef struct {
  UINT4 length;			/**< number of elements */
  vect3D_t *data;		/**< array of 3D vectors */
} vect3Dlist_t;

/** 'standard' representation of Doppler-paremeters */
typedef struct {
  vect2D_t vn;			/**< X, Y, component of ecliptic unit-vector to sky-location */
  PulsarSpins fkdot;		/**< vector of spins f^(k) */
} dopplerParams_t;

/** boundary of (single-hemisphere) search region in Doppler-space */
typedef struct {
  vect2Dlist_t skyRegion; 	/**< (ecliptic) 2D vector-polygon {kX, kY} defining the sky search-region */
  hemisphere_t hemisphere;	/**< hemisphere of sky vector-polygon */
  LIGOTimeGPS refTime;		/**< SSB reference GPS-time at which spin-range is defined */
  PulsarSpins fkdot;		/**< Vector of canonical spin-values w^(k) */
  PulsarSpins fkdotBand;	/**< Vector of canonical spin-bands Delta w^(k) */
} dopplerBoundary_t;

struct tagDopplerLatticeScan {
  scan_state_t state;		/**< current state of the scan: idle, ready of finished */
  REAL8 Tspan;			/**< total observation time spanned */
  dopplerBoundary_t boundary;	/**< boundary of search-space to cover */
  UINT4 dimLattice;		/**< dimension of nonzero-band search-lattice (can be <= dimCanonical) */
  gsl_vector *canonicalOrigin;	/**< 'origin' of the lattice in canonical coords {lb0, lb1, lb2, ... } */
  gsl_vector *dopplerUnits;	/**< conversion-factors from 'doppler' {f0, nX, nY, f1, ...} to 'canonical'-coords */
  gsl_matrix *latticeGenerator;	/**< generating matrix for the lattice: rows are the lattice-vectors */
  gsl_vector_int *latticeIndex;	/**< "Index" = multi-dim index-counters of current lattice point */
  gsl_vector_int *prevIndexBoundaryUp; 	/**< last positive turning point */
  gsl_vector_int *prevIndexCenter; /**< where was our index-center previously */
  gsl_vector_int *map2canonical;/**< mapping of lattice-Index into (canonical) Index {i0, i1, i2, ... } */
};

/*---------- empty initializers ---------- */
dopplerParams_t empty_dopplerParams;

/*---------- Global variables ----------*/
extern INT4 lalDebugLevel;

/*---------- internal function prototypes ----------*/
void skyRegionString2vect3D ( LALStatus *, vect3Dlist_t **skyRegionEcl, const CHAR *skyRegionString );
void setupSearchRegion ( LALStatus *status, DopplerLatticeScan *scan, const DopplerRegion *searchRegion );

hemisphere_t onWhichHemisphere ( const vect3Dlist_t *skypoints );
int skyposToVect3D ( vect3D_t *eclVect, const SkyPosition *skypos );
int vect2DToSkypos ( SkyPosition *skypos, vect2D_t * const vect2D, hemisphere_t hemi );
int findCenterOfMass ( vect3D_t *center, const vect3Dlist_t *points );

int IndexToCanonical ( gsl_vector **canonicalOffset, const gsl_vector_int *Index, const DopplerLatticeScan *scan );
int XLALIndexToDoppler ( dopplerParams_t *doppler, const gsl_vector_int *Index, const DopplerLatticeScan *scan );
int isIndexInsideBoundary ( const gsl_vector_int *Index, const DopplerLatticeScan *scan );

int convertDoppler2Canonical ( gsl_vector **canonicalPoint, const dopplerParams_t *doppler, const gsl_vector *dopplerUnits );
int convertCanonical2Doppler ( dopplerParams_t *doppler, const gsl_vector *canonical, const gsl_vector *dopplerUnits );

int vect2DInPolygon ( const vect2D_t *point, const vect2Dlist_t *polygon );
int isDopplerInsideBoundary ( const dopplerParams_t *doppler,  const dopplerBoundary_t *boundary );

int fprintf_vect2D ( FILE *fp, vect2D_t * const vect, hemisphere_t hemi );
int fprintf_vect2Dlist ( FILE *fp, const vect2Dlist_t *list, hemisphere_t hemi );
DopplerLatticeScan *XLALDuplicateDopplerLatticeScan ( const DopplerLatticeScan *scan );

/*==================== FUNCTION DEFINITIONS ====================*/

/* ------------------------------------------------------ */
/* -------------------- EXTERNAL API -------------------- */
/* ------------------------------------------------------ */

/** Initialize search-grid using optimal lattice-covering
 */
void
InitDopplerLatticeScan ( LALStatus *status,
			 DopplerLatticeScan **scan, 		/**< [out] initialized scan-state for lattice-scan */
			 const DopplerLatticeInit *init		/**< [in] scan init parameters */
			 )
{
  DopplerLatticeScan *ret = NULL;
  gsl_matrix *gij, *gijLattice;
  UINT4 i, j;

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
  if ( (gij = gsl_matrix_calloc (DIM_CANONICAL, DIM_CANONICAL)) == NULL ) {
    ABORT (status, DOPPLERSCANH_EMEM, DOPPLERSCANH_MSGEMEM);
  }

  if ( XLALFlatMetricCW ( gij, init->searchRegion.refTime, init->startTime, init->Tspan, init->ephemeris ) ) {
    ABORT ( status, DOPPLERSCANH_EXLAL, DOPPLERSCANH_MSGEXLAL );
  }
#ifdef DEBUG_ANS
  fprintf ( stderr, "\ngij = ");
  XLALfprintfGSLmatrix ( stderr, "% 15.9e ", gij );
#endif

  /* ----- reduce metric to subspace with actual nonzero search-bands !! ----- */
  if ( ( gijLattice = gsl_matrix_calloc ( ret->dimLattice, ret->dimLattice )) == NULL ) {
    ABORT (status, DOPPLERSCANH_EMEM, DOPPLERSCANH_MSGEMEM);
  }
  for ( i=0; i < ret->dimLattice; i ++ )
    {
      for ( j=i; j < ret->dimLattice; j ++ )
	{
	  UINT4 iMap, jMap;
	  REAL8 gijMap;
	  iMap = gsl_vector_int_get( ret->map2canonical, i );
	  jMap = gsl_vector_int_get( ret->map2canonical, j );
	  gijMap = gsl_matrix_get ( gij, iMap, jMap );
	  gsl_matrix_set ( gijLattice, i, j, gijMap );
	  gsl_matrix_set ( gijLattice, j, i, gijMap );
	} /* for j < dimLattice */
    } /* for i < dimLattice */

#ifdef DEBUG_ANS
  fprintf ( stderr, "gijLattice = ");
  XLALfprintfGSLmatrix ( stderr, "% 15.9e ", gijLattice );
#endif

  /* ----- compute generating matrix for the lattice ----- */
  if ( XLALFindCoveringGenerator (&(ret->latticeGenerator), LATTICE_TYPE_ANSTAR, sqrt(init->metricMismatch), gijLattice ) ) {
    ABORT ( status, DOPPLERSCANH_EXLAL, DOPPLERSCANH_MSGEXLAL );
  }
#ifdef DEBUG_ANS
  fprintf ( stderr, "generator = ");
  XLALfprintfGSLmatrix ( stderr, "% 15.9e ", ret->latticeGenerator );

  fprintf ( stderr, "origin = ");
  XLALfprintfGSLvector ( stderr, "% .9e ", ret->canonicalOrigin );
  fprintf ( stderr, "\n");
#endif

  /* ----- prepare Index-counters to generate lattice-points ----- */
  if ( (ret->latticeIndex = gsl_vector_int_calloc ( ret->dimLattice )) == NULL ) {
    ABORT (status, DOPPLERSCANH_EMEM, DOPPLERSCANH_MSGEMEM);
  }
  if ( (ret->prevIndexBoundaryUp = gsl_vector_int_calloc ( ret->dimLattice )) == NULL ) {
    ABORT (status, DOPPLERSCANH_EMEM, DOPPLERSCANH_MSGEMEM);
  }
  if ( (ret->prevIndexCenter = gsl_vector_int_calloc ( ret->dimLattice )) == NULL ) {
    ABORT (status, DOPPLERSCANH_EMEM, DOPPLERSCANH_MSGEMEM);
  }

  /* clean up memory */
  gsl_matrix_free ( gij );
  gsl_matrix_free ( gijLattice );

  /* return final scan-state */
  ret->state = STATE_READY;
  (*scan) = ret;

  /* debug */
  {
    PulsarDopplerParams doppler = empty_PulsarDopplerParams;
    XLALgetCurrentDopplerPos ( &doppler, ret, COORDINATESYSTEM_EQUATORIAL );
    fprintf ( stderr, "dopplerOrigin: ");
    fprintfDopplerParams ( stderr, &doppler );
  }


  DETATCHSTATUSPTR ( status );
  RETURN ( status );

} /* InitDopplerLatticeScan() */


/** Free an allocated DopplerLatticeScan.
 * Return: 0=OK, -1=ERROR
 */
int
XLALFreeDopplerLatticeScan ( DopplerLatticeScan **scan )
{
  if ( !scan || !(*scan) )
    return -1;

  if ( (*scan)->boundary.skyRegion.data )
    LALFree ( (*scan)->boundary.skyRegion.data );

  if ( (*scan)->canonicalOrigin )
    gsl_vector_free ( (*scan)->canonicalOrigin );

  if ( (*scan)->latticeGenerator )
    gsl_matrix_free ( (*scan)->latticeGenerator );

  if ( (*scan)->latticeIndex )
    gsl_vector_int_free ( (*scan)->latticeIndex );

  if ( (*scan)->prevIndexBoundaryUp )
    gsl_vector_int_free ( (*scan)->prevIndexBoundaryUp );

  if ( (*scan)->prevIndexCenter )
    gsl_vector_int_free ( (*scan)->prevIndexCenter );

  if ( (*scan)->map2canonical )
    gsl_vector_int_free ( (*scan)->map2canonical );

  LALFree ( (*scan) );
  (*scan) = NULL;

  return 0;

} /* XLALFreeDopplerLatticeScan() */


/** Return current lattice Index of the scan:
 * allocate Index here if *Index == NULL, otherwise use given vector
 * [has to have right dimension]
 */
int
XLALgetCurrentLatticeIndex ( gsl_vector_int **Index, const DopplerLatticeScan *scan  )
{
  const CHAR *fn = "XLALgetCurrentLatticeIndex()";

  if ( !Index || !scan || (scan->state != STATE_READY ) ) {
    XLAL_ERROR (fn, XLAL_EINVAL );
  }

  if ( *Index == NULL )	/* allocate output-vector */
    {
      if ( ((*Index) = gsl_vector_int_calloc ( scan->dimLattice )) == NULL ) {
	XLAL_ERROR (fn, XLAL_ENOMEM );
      }
    }
  else if ( (*Index)->size != scan->dimLattice ) {
    LALPrintError ("\n\nOutput vector has wrong dimension %d instead of %d!\n\n", (*Index)->size, scan->dimLattice);
    XLAL_ERROR (fn, XLAL_EINVAL );
  }

  gsl_vector_int_memcpy ( (*Index), scan->latticeIndex );

  return 0;
} /* XLALgetCurrentLatticeIndex() */

/** Set the current index of the scan.
 */
int
XLALsetCurrentLatticeIndex ( DopplerLatticeScan *scan, const gsl_vector_int *Index )
{
  const CHAR *fn = "XLALsetCurrentLatticeIndex()";

  if ( !Index || !scan || (scan->state != STATE_READY) || (Index->size != scan->dimLattice) ) {
    XLAL_ERROR (fn, XLAL_EINVAL );
  }

  gsl_vector_int_memcpy ( scan->latticeIndex, Index );

  return 0;
} /* XLALsetCurrentLatticeIndex() */

/** Count number of templates in lattice-grid by stepping through the whole
 * grid. NOTE: original scan-state isn't modified.
 *
 * Return: > 0: OK, -1=ERROR
 */
REAL8
XLALCountLatticeTemplates ( const DopplerLatticeScan *scan )
{
  int ret;
  REAL8 counter;
  DopplerLatticeScan *dupl;

 if ( !scan || scan->state != STATE_READY )
    return -1;

 if ( (dupl = XLALDuplicateDopplerLatticeScan ( scan )) == NULL )
   return -1;

 counter = 0;
 while ( (ret = XLALadvanceLatticeIndex ( dupl )) == 0 )
   counter ++;

 XLALFreeDopplerLatticeScan ( &dupl );

 if ( ret < 0 )
   return -1;
 else
   return counter;

} /* XLALCountLatticeTemplates() */


/** The central "lattice-stepping" function: advance to the 'next' Index-point, taking care not
 * to leave the boundary, and to cover the whole (convex!) search-region eventually!
 *
 * \note Algorithm:
 o) start with first Index-dimension aI = 0
 1) if Index(aI) >= 0: Index(aI) ++; else   Index(aI) --;
    i.e. increase for pos Index, decrease for negative ones ==> always walk "outwards"
    from origin, towards boundaries)
 o) if resulting point lies inside boundary ==> keep & return

 o) if boundary was crossed: return aI to origin: Index(aI) = 0;
 o) step to next dimension: aI ++;
 o) if no more dimensions left: aI >= dim ==> no further lattice-points! RETURN
 o) else continue at 1)
 *
 *
 * Return: 0=OK, -1=ERROR,  +1=No more lattice points left
 */
int
XLALadvanceLatticeIndex ( DopplerLatticeScan *scan )
{
  const CHAR *fn = "XLALadvanceLatticeIndex()";
  UINT4 dim, aI;
  int ret;
  gsl_vector_int *next_Index;

  if ( !scan || scan->state != STATE_READY ) {
    XLAL_ERROR (fn, XLAL_EINVAL );
  }

  dim = scan->dimLattice;

  if ( (next_Index = gsl_vector_int_calloc ( dim )) == NULL ) {
    XLAL_ERROR (fn, XLAL_ENOMEM );
  }

  gsl_vector_int_memcpy ( next_Index, scan->latticeIndex );

  aI = 0;
  do
    {
      int *pindex = gsl_vector_int_ptr (next_Index, aI);	/* pointer to current Index-element */
      int prev_center = gsl_vector_int_get ( scan->prevIndexCenter, aI ) ;
      /* check sign of current Index-element  */
      if ( ((*pindex) - prev_center) >= 0 )	/* step "up" if "above" center */
	{
	  (*pindex) ++;
	  ret = isIndexInsideBoundary ( next_Index, scan );
	  if ( ret < 0 ) return -1;	/* ERROR */
	  else if ( ret > 0 )
	    {
	      gsl_vector_int_memcpy ( scan->latticeIndex, next_Index );
	      return 0;	/* OK: found point */
	    }
	  else	/* reached "upper" boundary on aI */
	    {
	      UINT4 bI;
	      gsl_vector_int_set ( scan->prevIndexBoundaryUp, aI, (*pindex) - 1 ); /* keep track of new turning point */
	      (*pindex) = prev_center;	/* return to (previous) center of the line */
	      for (bI=0; bI < aI; bI ++ )
		{
		  gsl_vector_int_set ( next_Index, bI, 0 ); /* return to origin for all lower dimensions */
		  gsl_vector_int_set ( scan->prevIndexBoundaryUp, bI, 0 );
		  gsl_vector_int_set ( scan->prevIndexCenter, bI, 0 );
		}
	    } /* reached 'upper' boundary */
	} /* try a step "up" */

      /* try a step "down" */
      (*pindex) --;
      ret = isIndexInsideBoundary ( next_Index, scan );
      if ( ret < 0 ) return -1;	/* ERROR */
      else if ( ret > 0 )
	{
	  gsl_vector_int_memcpy ( scan->latticeIndex, next_Index );
	  return 0;	/* OK: found point */
	}

      /* reached "lower" boundary on aI ==> bump next dimension*/
      {
	/* try advancing to the next dimension
	 * the key-trick is that we re-center the point in aI, so we should
	 * have better chances of success in advancing.
	 * Note: this will probably only work for sure [if at all?]
	 * for convex regions to be filled
	 */
	int thisBoundaryDown = (*pindex) + 1;
	int thisBoundaryUp   = gsl_vector_int_get ( scan->prevIndexBoundaryUp, aI );
	int new_center = (int)(0.5 * (thisBoundaryUp + thisBoundaryDown) );
	(*pindex) = new_center;
	gsl_vector_int_set ( scan->prevIndexCenter, aI, new_center );
      }

      aI ++;

    } while ( aI < dim );

  gsl_vector_int_free ( next_Index );

  return 1;	/* no more points could be found! */

} /* XLALadvanceLatticeIndex() */


/** Return the current doppler-position {Freq, Alpha, Delta, f1dot, f2dot, ... } of the lattice-scan
 *  NOTE: the skyposition coordinate-system is chosen via 'skyCoords' in [EQUATORIAL, ECLIPTIC]
 */
int
XLALgetCurrentDopplerPos ( PulsarDopplerParams *pos, const DopplerLatticeScan *scan, CoordinateSystem skyCoords )
{
  const CHAR *fn = "XLALgetCurrentDopplerPos()";
  dopplerParams_t doppler = empty_dopplerParams;
  SkyPosition skypos = empty_SkyPosition;

  if ( !pos || !scan || (scan->state != STATE_READY) ) {
    XLAL_ERROR (fn, XLAL_EINVAL );
  }

  if ( XLALIndexToDoppler ( &doppler, scan->latticeIndex, scan ) ) {
    XLAL_ERROR (fn, XLAL_EFUNC );
  }

  skypos.system = skyCoords;
  if ( vect2DToSkypos ( &skypos, &(doppler.vn), scan->boundary.hemisphere) ) {
    XLAL_ERROR (fn, XLAL_EFUNC );
  }

  /* convert into PulsarDopplerParams type */
  pos->refTime 		= scan->boundary.refTime;
  pos->Alpha		= skypos.longitude;
  pos->Delta		= skypos.latitude;
  memcpy ( pos->fkdot, doppler.fkdot, sizeof(pos->fkdot) );
  pos->orbit 		= NULL;		/* FIXME: not supported yet */

  return 0;
} /* XLALgetCurrentDopplerPos() */



/* ------------------------------------------------------------ */
/* -------------------- INTERNAL functions -------------------- */
/* ------------------------------------------------------------ */

/** Return a copy a full DopplerLatticeScan state, useful as a backup
 * while counting the lattice.
 *
 * Return: NULL=ERROR
 */
DopplerLatticeScan *
XLALDuplicateDopplerLatticeScan ( const DopplerLatticeScan *scan )
{
  DopplerLatticeScan *ret;
  size_t size;
  void *src, *dst;

  if (!scan || scan->state == STATE_IDLE )
    return NULL;

  if ( (ret = LALCalloc ( 1, sizeof(*ret) )) == NULL )
    return NULL;

  (*ret) = (*scan);	/* struct-copy everything non-alloc'ed */

  /* properly copy everything alloc'ed */

  /* 'char*' */
  if ( (src = scan->boundary.skyRegion.data) )
    {
      size = sizeof(*scan->boundary.skyRegion.data);
      if ( (dst = LALCalloc(1, size)) == NULL  )
	return NULL;
      memcpy ( dst, src, size );
      ret->boundary.skyRegion.data = dst;
    }

  /* gsl_vector */
  if ( (src = scan->canonicalOrigin) )
    {
      size = scan->canonicalOrigin->size;
      if ( (dst = gsl_vector_alloc(size)) == NULL  )
	return NULL;
      gsl_vector_memcpy ( dst, src );
      ret->canonicalOrigin = dst;
    }

  /* gsl_matrix */
  if ( (src = scan->latticeGenerator) )
    {
      size = scan->latticeGenerator->size1;
      if ( (dst = gsl_matrix_alloc(size,size)) == NULL  )
	return NULL;
      gsl_matrix_memcpy ( dst, src );
      ret->latticeGenerator = dst;
    }

  /* gsl_vector_int */
  if ( (src = scan->latticeIndex) )
    {
      size = scan->latticeIndex->size;
      if ( (dst = gsl_vector_int_alloc(size)) == NULL  )
	return NULL;
      gsl_vector_int_memcpy ( dst, src );
      ret->latticeIndex = dst;
    }
  if ( (src = scan->prevIndexBoundaryUp) )
    {
      size = scan->prevIndexBoundaryUp->size;
      if ( (dst = gsl_vector_int_alloc(size)) == NULL  )
	return NULL;
      gsl_vector_int_memcpy ( dst, src );
      ret->prevIndexBoundaryUp = dst;
    }
  if ( (src = scan->prevIndexCenter) )
    {
      size = scan->prevIndexCenter->size;
      if ( (dst = gsl_vector_int_alloc(size)) == NULL  )
	return NULL;
      gsl_vector_int_memcpy ( dst, src );
      ret->prevIndexCenter = dst;
    }
  if ( (src = scan->map2canonical) )
    {
      size = scan->map2canonical->size;
      if ( (dst = gsl_vector_int_alloc(size)) == NULL  )
	return NULL;
      gsl_vector_int_memcpy ( dst, src );
      ret->map2canonical = dst;
    }

  return ret;

} /* XLALDuplicateDopplerLatticeScan() */



/** Convert given Index into doppler-params
 * Return: 0=OK, -1=ERROR
 */
int
XLALIndexToDoppler ( dopplerParams_t *doppler, const gsl_vector_int *Index, const DopplerLatticeScan *scan )
{
  const CHAR *fn = "XLALIndexToDoppler()";
  gsl_vector *canonical = NULL;

  if ( !doppler || !Index || !scan ) {
    XLAL_ERROR (fn, XLAL_EINVAL );
  }

  if ( IndexToCanonical ( &canonical, Index, scan ) ) {
    XLAL_ERROR (fn, XLAL_EFUNC );
  }

  if ( convertCanonical2Doppler ( doppler, canonical, scan->dopplerUnits ) ) {
    gsl_vector_free ( canonical );
    XLAL_ERROR (fn, XLAL_EFUNC );
  }
  gsl_vector_free ( canonical );

  return 0;

} /* XLALIndexToDoppler() */



/** Determine whether the given lattice-Index corresponds to a Doppler-point
 * that lies within the search-boundary
 *  Return: TRUE, FALSE,  -1 = ERROR
 */
int
isIndexInsideBoundary ( const gsl_vector_int *Index, const DopplerLatticeScan *scan )
{
  dopplerParams_t doppler = empty_dopplerParams;

  if ( !Index || !scan || scan->state != STATE_READY )
    return -1;

  if ( XLALIndexToDoppler ( &doppler, Index, scan ) )
    return -1;

  return ( isDopplerInsideBoundary ( &doppler,  &scan->boundary ) );

} /* isIndexInsideBoundary() */


/** Determine whether the given Doppler-point lies within the search-boundary
 *  Return: TRUE, FALSE,  -1 = ERROR
 */
int
isDopplerInsideBoundary ( const dopplerParams_t *doppler,  const dopplerBoundary_t *boundary )
{
  vect2D_t skyPoint = {0, 0};
  int insideSky, insideSpins;
  UINT4 s;

  if ( !doppler || !boundary )
    return -1;

  skyPoint[0] = doppler->vn[0];
  skyPoint[1] = doppler->vn[1];

  if ( (insideSky = vect2DInPolygon ( (const vect2D_t*)skyPoint, &(boundary->skyRegion) )) < 0 )
    return -1;

#ifdef DEBUG_ANS
  /* ===== debug ===== */
  {
    FILE *fp = fopen ("gridpoints.dat", "ab");
    fprintf (fp, "%f, %f, ", skyPoint[0], skyPoint[1] );
    fprintf_vect2D ( fp, (const vect2D_t*)skyPoint, boundary->hemisphere );
    fclose(fp);
  }
#endif

  insideSpins = TRUE;
  for ( s=0; s < NUM_SPINS; s ++ )
    {
      double fkdotMax = boundary->fkdot[s] + boundary->fkdotBand[s];
      double fkdotMin = boundary->fkdot[s];
      if ( ( gsl_fcmp ( doppler->fkdot[s], fkdotMax, EPS_REAL8 ) > 0 ) ||
	   ( gsl_fcmp ( doppler->fkdot[s], fkdotMin, EPS_REAL8 ) < 0 ) )
	insideSpins = ( insideSpins && FALSE );
    } /* for s < NUM_SPINS */

  return ( insideSky && insideSpins );

} /* insideBoundary() */


/** Translate the input 'DopplerRegion' into an internal representation for the scan.
 *
 * \note DopplerLatticeScan->Tspan must have been set! This is used for conversion Doppler -> canonical
 */
void
setupSearchRegion ( LALStatus *status, DopplerLatticeScan *scan, const DopplerRegion *searchRegion )
{
  UINT4 i;
  vect3Dlist_t *points3D = NULL;

  INITSTATUS( status, "InitDopplerLatticeScan", DOPPLERLATTICECOVERING );
  ATTATCHSTATUSPTR ( status );

  ASSERT ( scan, status, DOPPLERSCANH_ENULL, DOPPLERSCANH_MSGENULL );
  ASSERT ( searchRegion, status, DOPPLERSCANH_ENULL, DOPPLERSCANH_MSGENULL );

  /* ----- convert sky-string to polygon ----- */
  TRY ( skyRegionString2vect3D ( status->statusPtr, &points3D, searchRegion->skyRegionString ), status );

  /* ----- map lattice-indices to Doppler-dimensions with nonzero search-bands ----- */
  {
    gsl_vector_int *mapTmp;
    UINT4 s, l;
    if ( (mapTmp = gsl_vector_int_calloc ( DIM_CANONICAL )) == NULL ) { /* <= dimLattice */
      ABORT (status, DOPPLERSCANH_EMEM, DOPPLERSCANH_MSGEMEM);
    }
    l=0;
    if ( searchRegion->fkdotBand[0] )
      gsl_vector_int_set ( mapTmp, l++, 0 );
    if ( points3D->length > 1 )		/* more than one point: must be a 2D area */
      {
	gsl_vector_int_set ( mapTmp, l++, 1 );
	gsl_vector_int_set ( mapTmp, l++, 2 );
      }
    for (s=1; s < NUM_SPINS; s ++ )
      {
	if ( searchRegion->fkdotBand[s] )
	  gsl_vector_int_set ( mapTmp, l++, 2 + s );
      } /* for s < NUM_SPINS */

    /* cut down to minimum size */
    if ( (l > 0) && (scan->map2canonical = gsl_vector_int_calloc ( l ) ) == NULL ) {
      ABORT (status, DOPPLERSCANH_EMEM, DOPPLERSCANH_MSGEMEM);
    }
    for ( i=0; i < l; i ++ )
      gsl_vector_int_set ( scan->map2canonical, i, gsl_vector_int_get ( mapTmp, i ) );
    gsl_vector_int_free ( mapTmp );

    scan->dimLattice = scan->map2canonical->size;	/* dimension of actual search-lattice */

  } /* find mapping lattice <--> dopplerSpace */


  /* ----- get unit-conversion factors from 'doppler' --> 'canonical' ----- */
  {
    REAL8 FreqMax = searchRegion->fkdot[0] + searchRegion->fkdotBand[0];
    REAL8 Tspan = scan->Tspan;
    REAL8 convFact;
    UINT4 s;

    if ( (scan->dopplerUnits = gsl_vector_calloc ( DIM_CANONICAL )) == NULL ) {
      ABORT (status, DOPPLERSCANH_EMEM, DOPPLERSCANH_MSGEMEM);
    }

    convFact = - (LAL_TWOPI * LAL_AU_SI / LAL_C_SI ) * FreqMax;	/* convert nX,Y --> kX,Y */
    gsl_vector_set ( scan->dopplerUnits, 1, convFact );
    gsl_vector_set ( scan->dopplerUnits, 2, convFact );

    convFact = LAL_TWOPI * Tspan;
    gsl_vector_set ( scan->dopplerUnits, 0, convFact );		/* convert Freq --> w0 */
    for ( s=1; s < NUM_SPINS; s ++ )
      {
	convFact *= Tspan;
	gsl_vector_set ( scan->dopplerUnits, s+2, convFact );	/* convert fkdot --> wk */
      }

  } /* get unit-conversions */

  if ( (scan->boundary.hemisphere = onWhichHemisphere ( points3D )) == HEMI_BOTH )
    {
      LALPrintError ("\n\nSorry, currently only (ecliptic) single-hemisphere sky-regions are supported!\n");
      ABORT ( status, DOPPLERSCANH_EINPUT, DOPPLERSCANH_MSGEINPUT );
    }
  if ( (scan->boundary.skyRegion.data = LALCalloc ( points3D->length, sizeof(scan->boundary.skyRegion.data[0]) )) == NULL ) {
    ABORT (status, DOPPLERSCANH_EMEM, DOPPLERSCANH_MSGEMEM);
  }
  scan->boundary.skyRegion.length = points3D->length;

  for ( i=0; i < points3D->length; i ++ ) {
    scan->boundary.skyRegion.data[i][0] = points3D->data[i][0];
    scan->boundary.skyRegion.data[i][1] = points3D->data[i][1];
  }

#ifdef DEBUG_ANS
  /* ===== debug ===== */
  {
    FILE *fp = fopen ( "boundary.dat", "wb" );

    for (i=0; i < scan->boundary.skyRegion.length; i ++ )
      fprintf (fp, "%f, %f\n", scan->boundary.skyRegion.data[i][0], scan->boundary.skyRegion.data[i][1] );
    fprintf (fp, "\n\n");

    fprintf_vect2Dlist ( fp, &(scan->boundary.skyRegion), scan->boundary.hemisphere );

    fclose ( fp );
  }
#endif


  /* ----- spins ----- */
  scan->boundary.refTime   = searchRegion->refTime;
  memcpy ( scan->boundary.fkdot, searchRegion->fkdot, sizeof(searchRegion->fkdot) );
  memcpy ( scan->boundary.fkdotBand, searchRegion->fkdotBand, sizeof(searchRegion->fkdotBand) );


  /* ----- use the center of the searchRegion as origin of the lattice ----- */
  {
    vect3D_t com = {0,0,0};		/* center-of-mass of points3D */
    dopplerParams_t midPoint = empty_dopplerParams;
    UINT4 s;
    findCenterOfMass ( &com, points3D );
    midPoint.vn[0] = com[0];
    midPoint.vn[1] = com[1];
    for ( s=0; s < NUM_SPINS; s ++ )
      midPoint.fkdot[s] = searchRegion->fkdot[s] + 0.5 * searchRegion->fkdotBand[s];
    if ( convertDoppler2Canonical ( &(scan->canonicalOrigin), &midPoint, scan->dopplerUnits) ) {
      ABORT ( status, DOPPLERSCANH_EXLAL, DOPPLERSCANH_MSGEXLAL );
    }
  }

  /* cleanup */
  LALFree ( points3D->data );
  LALFree ( points3D );

  DETATCHSTATUSPTR ( status );
  RETURN ( status );

} /* setupSearchRegion() */

/** Convert Index-vector {i0, i1, i2, ..} into canonical param-vector {w0, kX, kY, w1, w2, ...}
 *
 * \note the output-vector 'canonical' can be NULL, in which case it will be allocated here,
 * If allocated already, its dimension must be dimCanonical
 *
 * Return: 0=OK, -1=ERROR
 */
int
IndexToCanonical ( gsl_vector **canonical, const gsl_vector_int *Index, const DopplerLatticeScan *scan )
{
  UINT4 i, j;

  /* check input */
  if ( !canonical || !Index || !scan)
    return -1;

  if ( (*canonical != NULL) && ( (*canonical)->size != DIM_CANONICAL ) ) {
    LALPrintError ("\nIndexToCanonicalOffset(): output-vector has dim=%d instead of dimCanonical=%d!\n\n",
		   (*canonical)->size, DIM_CANONICAL );
    return -1;
  }

  /* allocate and initialized output-vector to zero */
  if ( (*canonical) == NULL ) {
    if ( ((*canonical) = gsl_vector_calloc ( DIM_CANONICAL )) == NULL )
      return -1;
  }
  gsl_vector_set_zero ( *canonical );

  /* lattice-vect^i = Index_j  generator^(ji) */
  for ( i=0; i < scan->dimLattice; i ++ )
    {
      REAL8 comp = 0;
      UINT4 iMap = gsl_vector_int_get ( scan->map2canonical, i );

      for (j=0; j < scan->dimLattice; j ++ )
	comp += 1.0 * gsl_vector_int_get ( Index, j ) * gsl_matrix_get ( scan->latticeGenerator, j, i );

      gsl_vector_set ( (*canonical), iMap, comp );

    } /* i < dim */

  if ( gsl_vector_add ( *canonical, scan->canonicalOrigin ) ) {
    return -1;
  }

  return 0;

} /* IndexToCanonical() */


/** Convert Doppler-parameters from {nX, nY, nZ, fkdot} into internal 'canonical' form
 *  {w0, kX, kY, w1, w1, ... }
 * \note If the return-vector is NULL, it is allocated here, otherwise it has to
 * have dimension DIM_CANONICAL
 *
 * Return: 0=OK, -1=ERROR
 */
int
convertDoppler2Canonical ( gsl_vector **canonical, const dopplerParams_t *doppler, const gsl_vector *dopplerUnits )
{
  UINT4 s;

  /* check input */
  if ( !canonical || !doppler || !dopplerUnits )
    return -1;

  if ( (*canonical) == NULL )
    {
      if ( ((*canonical) = gsl_vector_calloc ( DIM_CANONICAL )) == NULL )
	return -1;
    }
  else
    {
      if ( (*canonical)->size != DIM_CANONICAL )
	return -1;
    }

  gsl_vector_set ( (*canonical), 0, gsl_vector_get(dopplerUnits, 0) * doppler->fkdot[0] );

  gsl_vector_set ( (*canonical), 1, gsl_vector_get( dopplerUnits, 1) * doppler->vn[0] );
  gsl_vector_set ( (*canonical), 2, gsl_vector_get( dopplerUnits, 2) * doppler->vn[1] );

  for ( s=1; s < NUM_SPINS; s ++ )
    gsl_vector_set ( (*canonical), 2 + s, gsl_vector_get( dopplerUnits, 2 + s ) * doppler->fkdot[s] );

  return 0;

} /* convertDoppler2Canonical() */


/** Convert a 'canonical' parameter-space point back into a 'physical' Doppler units:
 * a unit (ecliptic) sky-vector 'vn' and spin-vector 'fkdot'
 *
 * Return: 0=OK, -1=ERROR
 */
int
convertCanonical2Doppler ( dopplerParams_t *doppler, const gsl_vector *canonical, const gsl_vector *dopplerUnits )
{
  UINT4 s;

  /* check input */
  if ( !canonical || !doppler || !dopplerUnits )
    return -1;
  if ( (canonical->size != DIM_CANONICAL) || (dopplerUnits->size != DIM_CANONICAL) )
    return -1;

  /* Freq */
  doppler->fkdot[0] = gsl_vector_get ( canonical, 0 ) / gsl_vector_get ( dopplerUnits, 0 );
  /* sky */
  doppler->vn[0]    = gsl_vector_get ( canonical, 1 ) / gsl_vector_get ( dopplerUnits, 1 );
  doppler->vn[1]    = gsl_vector_get ( canonical, 2 ) / gsl_vector_get ( dopplerUnits, 2 );
  /* fkdot */
  for ( s=1; s < NUM_SPINS; s ++ )
    doppler->fkdot[s] = gsl_vector_get ( canonical, s+2 ) / gsl_vector_get ( dopplerUnits, s+2 );

  return 0;

} /* convertCanonical2Doppler() */

/** Convert a sky-region string into a list of vectors in ecliptic 2D coords {nX, nY}
 */
void
skyRegionString2vect3D ( LALStatus *status,
			 vect3Dlist_t **skyRegion, 	/**< [out] list of skypoints in 3D-ecliptic coordinates {nX, nY, nZ} */
			 const CHAR *skyRegionString )	/**< [in] string of equatorial-coord. sky-positions */
{
  vect3Dlist_t *ret = NULL;
  SkyRegion region = empty_SkyRegion;
  UINT4 i, j, refineby, N0, N1;

  INITSTATUS( status, "skyRegionString2vect3D()", DOPPLERLATTICECOVERING );
  ATTATCHSTATUSPTR ( status );

  ASSERT ( skyRegion, status, DOPPLERSCANH_ENULL, DOPPLERSCANH_MSGENULL );
  ASSERT ( *skyRegion == NULL, status, DOPPLERSCANH_ENONULL, DOPPLERSCANH_MSGENONULL );
  ASSERT ( skyRegionString, status, DOPPLERSCANH_ENULL, DOPPLERSCANH_MSGENULL );

  TRY ( ParseSkyRegionString (status->statusPtr, &region, skyRegionString), status );

  N0 = region.numVertices;

  /* Note: the projection {alpha, delta} -> {nX, nY} distorts straight lines
   * we therefore try to 'refine' all connecting lines in {alpha, delta} by
   * 'refineby', in order to reduce the distortions to the intended polygon
   */
  if ( N0 >= 3)
    refineby = 10;
  else
    refineby = 1;

  N1 = N0 * refineby;

  if ( ( ret = LALCalloc ( 1, sizeof(*ret)) ) == NULL ) {
    ABORT (status, DOPPLERSCANH_EMEM, DOPPLERSCANH_MSGEMEM);
  }
  if ( (ret->data = LALCalloc ( N1, sizeof( ret->data[0] ) )) == NULL ) {
    LALFree ( ret );
    ABORT (status, DOPPLERSCANH_EMEM, DOPPLERSCANH_MSGEMEM);
  }
  ret->length = N1;

  for ( i=0; i < N0; i ++ )
    {
      SkyPosition inter;
      SkyPosition p0 = region.vertices[i];
      SkyPosition p1 = region.vertices[(i+1) % N0];
      REAL8 dAlpha, dDelta;

      inter.system = region.vertices[0].system;
      dAlpha = (p1.longitude - p0.longitude) / refineby;
      dDelta = (p1.latitude  - p0.latitude) / refineby;
      for ( j=0; j < refineby; j ++ )
	{
	  inter.longitude = p0.longitude + j * dAlpha;
	  inter.latitude  = p0.latitude  + j * dDelta;
	  skyposToVect3D ( &(ret->data[i*refineby + j]), &inter );
	} /* j < refineby */
    } /* i < N0 */

  /* free memory */
  LALFree ( region.vertices );

  /* return sky-region */
  (*skyRegion) = ret;

  DETATCHSTATUSPTR ( status );
  RETURN ( status );

} /* skyRegionString2list() */

/** Check whether given list of skypoint lie on a single hemisphere
 */
hemisphere_t
onWhichHemisphere ( const vect3Dlist_t *skypoints )
{
  UINT4 i;
  hemisphere_t our_hemi = HEMI_BOTH;

  if ( !skypoints || (skypoints->length == 0) )
    return FALSE;

  for ( i=0; i < skypoints->length; i ++ )
    {
      vect3D_t *thisPoint = &(skypoints->data[i]);
      hemisphere_t this_hemi = VECT_HEMI ( *thisPoint );
      if ( (our_hemi == HEMI_BOTH) && (this_hemi != HEMI_BOTH) ) /* set our_hemi to first non-zero hemisphere */
	our_hemi = this_hemi;
      if ( (this_hemi != HEMI_BOTH) && (this_hemi != our_hemi ) )
	return HEMI_BOTH;
    }

  return our_hemi;

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
 *
 * NOTE: The coordinate-system of the result is determined by the 'system'-setting in
 * the output 'skypos'
 *
 * return: 0=OK, -1=ERROR
 */
int
vect2DToSkypos ( SkyPosition *skypos, vect2D_t * const vect2D, hemisphere_t hemi )
{
  REAL8 invnorm;
  vect3D_t nvect = {0,0,0};
  vect3D_t vn3D = {0,0,0};
  REAL8 vn2;
  REAL8 longitude, latitude;
  REAL8 sineps, coseps;

  if ( !skypos || !vect2D )
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

  vn3D[0] = (*vect2D)[0];
  vn3D[1] = (*vect2D)[1];

  vn2 = SQUARE ( vn3D[0] ) + SQUARE ( vn3D[1] );
  if ( gsl_fcmp ( vn2, 1.0, EPS_REAL8 ) > 0 ) {
    LALPrintError ( "\n\nvect2DToSkypos(): Sky-vector has length > 1! (vn2 = %f)\n\n", vn2 );
    vn3D[2] = 0;
  }
  else
    vn3D[2] = sqrt ( fabs ( 1.0 - vn2 ) );		/* nZ = sqrt(1 - nX^2 - nY^2 ) */

  if ( hemi == HEMI_SOUTH )
    vn3D[2] *= -1;

  nvect[0] = vn3D[0];
  nvect[1] = coseps * vn3D[1] - sineps * vn3D[2];
  nvect[2] = sineps * vn3D[1] + coseps * vn3D[2];

  invnorm = 1.0 / VECT_NORM ( nvect );

  VECT_MULT ( nvect, invnorm );

  longitude = atan2 ( nvect[1], nvect[0]);
  if ( longitude < 0 )
    longitude += LAL_TWOPI;

  latitude = asin ( nvect[2] );

  skypos->longitude = longitude;
  skypos->latitude = latitude;

  return 0;

} /* vect2DToSkypos() */


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


/** Function for checking if a given vect2D-point lies inside or outside a given vect2D-polygon.
 * This is basically indentical to 'pointInPolygon()' only using different data-types.
 *
 * \par Note1:
 * 	The list of polygon-points must not close on itself, the last point
 * 	is automatically assumed to be connected to the first
 *
 * \par Alorithm:
 *     Count the number of intersections of rays emanating to the right
 *     from the point with the lines of the polygon: even=> outside, odd=> inside
 *
 * \par Note2:
 *     we try to get this algorith to count all boundary-points as 'inside'
 *     we do this by counting intersection to the left _AND_ to the right
 *     and consider the point inside if either of those says its inside...
 *
 * \par Note3:
 *     correctly handles the case of a 1-point 'polygon', in which the two
 *     points must agree within eps=1e-10 relative precision.
 *
 * \return : TRUE or FALSE, -1=ERROR
 *----------------------------------------------------------------------*/
int
vect2DInPolygon ( const vect2D_t *point, const vect2Dlist_t *polygon )
{
  UINT4 i;
  UINT4 N;
  UINT4 insideLeft, insideRight;
  BOOLEAN inside = 0;
  vect2D_t *vertex;
  REAL8 xinter, v1x, v1y, v2x, v2y, px, py;

  if (!point || !polygon || !polygon->data )
    return -1;

  /* convenience variables */
  vertex = polygon->data;
  N = polygon->length; 	/* num of vertices = num of edges */
  px = (*point)[0];
  py = (*point)[1];

  /* treat special case of 1-point 'polygon' ==> float-number comparisons */
  if ( N == 1 )
    {
      int diffx, diffy;
      double eps = 1e-10;	/* allow relative errors of up to 1e-10 */
      diffx = gsl_fcmp ( vertex[0][0], px, eps );
      diffy = gsl_fcmp ( vertex[0][1], py, eps );
      return ( (diffx == 0) && (diffy == 0) );
    }
  else if ( N < 3 )	/* need 3 points for an area */
    return -1;

  /* general case of 2D-polygon */

  insideLeft = insideRight = 0;

  for (i=0; i < N; i++)
    {
      v1x = vertex[i][0];
      v1y = vertex[i][1];
      v2x = vertex[(i+1) % N][0];
      v2y = vertex[(i+1) % N][1];

      /* pre-select candidate edges */
      if ( (py <  MIN(v1y,  v2y)) || (py >=  MAX(v1y, v2y) ) || (v1y == v2y) )
	continue;

      /* now calculate the actual intersection point of the horizontal ray with the edge in question*/
      xinter = v1x + (py - v1y) * (v2x - v1x) / (v2y - v1y);

      if (xinter > px)	      /* intersection lies to the right of point */
	insideLeft ++;

      if (xinter < px)       /* intersection lies to the left of point */
	insideRight ++;

    } /* for sides of polygon */

  inside = ( ((insideLeft %2) == 1) || (insideRight %2) == 1);

  return inside;

} /* vect2DInPolygon() */

int
fprintf_vect2D ( FILE *fp, vect2D_t * const vect, hemisphere_t hemi )
{
  SkyPosition skypos;

  if ( !vect || !fp )
    return -1;

  skypos.system = COORDINATESYSTEM_EQUATORIAL;

  vect2DToSkypos ( &skypos, vect, hemi );

  fprintf (fp, "%f, %f\n ", skypos.longitude, skypos.latitude );

  return 0;
} /* fprintf_vect2D() */

int
fprintf_vect2Dlist ( FILE *fp, const vect2Dlist_t *list, hemisphere_t hemi )
{
  UINT4 i;

  if ( !list || !fp )
    return -1;

  for ( i=0; i < list->length; i ++ )
    fprintf_vect2D ( fp, &(list->data[i]), hemi );

  return 0;
}


