/*----------------------------------------------------------------------- 
 * 
 * File Name: Overlap.c 
 * 
 * Author: J.D. Romano
 * 
 * Revision: $Id$
 * 
 *----------------------------------------------------------------------- 
 * 
 * NAME 
 * LALOverlap 
 * 
 * SYNOPSIS 
 * void LALOverlap (LALStatus *, REAL4Vector *, OverlapParameters *);
 *
 * typedef struct tagOverlapParameters {
 *   IFOsite   site1ID; 
 *   IFOsite   site2ID; 
 *   INT4      length;  
 *   REAL4     deltaF;  
 * } OverlapParameters;
 *
 * typedef enum {
 *   LHO, 
 *   LLO, 
 *   VIRGO,
 *   GEO600,
 *   TAMA300
 * } IFOsite;
 * 
 * DESCRIPTION 
 * Calculates the values of the overlap reduction function for a pair of 
 * interferometers.
 *
 * DIAGNOSTICS 
 * null pointer to input parameters
 * site ID having illegitimate enum type
 * specified length of output vector <= 0
 * frequency spacing <= 0
 * null pointer to output vector 
 * length of output vector not equal to length specified in input parameters
 * null pointer to data member of output vector
 *
 * CALLS
 * 
 * NOTES
 * 
 *-----------------------------------------------------------------------
 */

#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <math.h>
#include <lal/Overlap.h>

NRCSID (OVERLAPC, "$Id$");

static void
GetSiteParameters ( SiteParameters*, IFOsite );

static void
CalcSiteCoordinates ( SiteCoordinates*, SiteParameters* );

static const CHAR* siteNames[]=SITENAMELIST;

void
LALOverlap ( LALStatus            *status,
	  REAL4Vector       *vector,
	  OverlapParameters *parameters )
{
  INT4             site1ID;
  INT4             site2ID;
  INT4             length;
  REAL4            deltaF;
       
  SiteParameters   site1Parameters;
  SiteParameters   site2Parameters;

  SiteCoordinates  site1Coordinates;
  SiteCoordinates  site2Coordinates;

  INT4   i;
  REAL4  f;
  REAL4  sep[3];
  REAL4  distance;
  REAL4  x1[3], y_1[3];
  REAL4  x2[3], y2[3];
  REAL4  x1DOTx2,  x1DOTy2,  y1DOTy2,  y1DOTx2;
  REAL4  x1DOTsep, x2DOTsep, y1DOTsep, y2DOTsep;
  REAL4  c1, c2, c3;
  REAL8  alpha, alpha2, alpha4;
  REAL8  s, c;
  REAL8  besselJ0, besselJ1, besselJ2;

  /* initialize status structure */
  INITSTATUS( status, "LALOverlap", OVERLAPC );

  /* check that pointer to input parameters is not null */
  ASSERT(parameters!=NULL, status, OVERLAP_ENULLIP, OVERLAP_MSGENULLIP);

  /* check that site IDs have legitimate values */
  site1ID = parameters->site1ID;
  site2ID = parameters->site2ID;
  ASSERT(LHO<=site1ID && site1ID<NUMBEROFSITES, status,
         OVERLAP_ESITE, OVERLAP_MSGESITE);
  ASSERT(LHO<=site2ID && site2ID<NUMBEROFSITES, status,
         OVERLAP_ESITE, OVERLAP_MSGESITE);

  /* check that specified length of output vector is > 0 */ 
  length = parameters->length;
  ASSERT(length>0, status, OVERLAP_ESIZE, OVERLAP_MSGESIZE);

  /* check that frequency spacing is > 0 */ 
  deltaF = parameters->deltaF;
  ASSERT(deltaF>0, status, OVERLAP_EDELTAF, OVERLAP_MSGEDELTAF);

  /* check that pointer to output vector is not null */
  ASSERT(vector!=NULL, status, OVERLAP_ENULLOP, OVERLAP_MSGENULLOP);

  /* check that output vector length agrees with length specified in */
  /* input parameters */
  ASSERT((INT4)vector->length==length, status, OVERLAP_ESIZEMM, OVERLAP_MSGESIZEMM);

  /* check that pointer to data member of output vector is not null */
  ASSERT(vector->data!=NULL, status, OVERLAP_ENULLD, OVERLAP_MSGENULLD);

  /* everything okay here --------------------------------------------*/

  /* get parameters for each site */
  GetSiteParameters( &site1Parameters, site1ID );
  GetSiteParameters( &site2Parameters, site2ID );

  /* calculate coordinates for each site */
  CalcSiteCoordinates( &site1Coordinates, &site1Parameters );
  CalcSiteCoordinates( &site2Coordinates, &site2Parameters );

  /* calculate separation vector between sites */
  /* rename arm vectors */
  for ( i=0; i<3; i++) {
    sep[i] = site1Coordinates.vertex[i]-site2Coordinates.vertex[i];
    x1[i]  = site1Coordinates.arm1[i];
    y_1[i]  = site1Coordinates.arm2[i];
    x2[i]  = site2Coordinates.arm1[i];
    y2[i]  = site2Coordinates.arm2[i];
  }
  
  /* calculate distance between sites, in cm */
  distance = sqrt( sep[0]*sep[0] + sep[1]*sep[1] + sep[2]*sep[2] );

  /* calculate unit separation vector */
  if ( distance != 0 ) {
    for ( i=0; i<3; i++ ) {
      sep[i] /= distance;
    }
  }

  /* calculate dot products of unit arm vectors and separation vector */
  x1DOTx2   = x1[0]*x2[0]  + x1[1]*x2[1]  + x1[2]*x2[2];
  y1DOTy2   = y_1[0]*y2[0]  + y_1[1]*y2[1]  + y_1[2]*y2[2]; 
  x1DOTy2   = x1[0]*y2[0]  + x1[1]*y2[1]  + x1[2]*y2[2];
  y1DOTx2   = y_1[0]*x2[0]  + y_1[1]*x2[1]  + y_1[2]*x2[2];
  x1DOTsep  = x1[0]*sep[0] + x1[1]*sep[1] + x1[2]*sep[2];
  x2DOTsep  = x2[0]*sep[0] + x2[1]*sep[1] + x2[2]*sep[2];
  y1DOTsep  = y_1[0]*sep[0] + y_1[1]*sep[1] + y_1[2]*sep[2];
  y2DOTsep  = y2[0]*sep[0] + y2[1]*sep[1] + y2[2]*sep[2];
  
  /* calculate coefficients c1, c2, c3 for overlap reduction funtion */

  /* c1 = d1 : d2 */
  c1  = x1DOTx2*x1DOTx2 + y1DOTy2*y1DOTy2 - x1DOTy2*x1DOTy2 - y1DOTx2*y1DOTx2;
  c1 *= 0.25;

  /* c2 = s . d1 . d2 . s */
  c2  = x1DOTsep*x2DOTsep*x1DOTx2 + y1DOTsep*y2DOTsep*y1DOTy2 
      - y1DOTsep*x2DOTsep*y1DOTx2 - x1DOTsep*y2DOTsep*x1DOTy2;
  c2 *= 0.25;

  /* c3 = (s . d1 . s)(s . d2 . s) */
  c3  = (x1DOTsep*x1DOTsep - y1DOTsep*y1DOTsep)*
        (x2DOTsep*x2DOTsep - y2DOTsep*y2DOTsep);
  c3 *= 0.25;

  /* calculate overlap reduction function at each discrete frequency */
  for ( i=0; i<length; i++ ) {

    f = i*deltaF;
    alpha = 2*LAL_PI*f*distance/C_LIGHT;
    alpha2 = alpha*alpha;

    /* if the argument is close to zero, use power series */
    if ( alpha<0.01 ) {
      alpha4 = alpha2*alpha2;
      besselJ0 =  1.0 - alpha2/6.0  + alpha4/120.0;
      besselJ1 = (1.0 - alpha2/10.0 + alpha4/280.0)/3.0;
      besselJ2 = (1.0 - alpha2/14.0 + alpha4/504.0)/15.0;
    }
    else {
      s = sin(alpha);
      c = cos(alpha);

      /* the standard spherical bessel functions */
      besselJ0 = s/alpha;
      besselJ1 = (besselJ0 - c)/alpha;
      besselJ2 = (3.0*besselJ1 - s)/alpha;

      /* slightly modified spherical bessel functions */
      besselJ1 /= alpha;
      besselJ2 /= alpha2;
    }
    
    /* calculate values of overlap reduction function */
    vector->data[i] = 
      c1*(  5.0*besselJ0 - 10.0*besselJ1 +  5.0*besselJ2) +
      c2*(-10.0*besselJ0 + 40.0*besselJ1 - 50.0*besselJ2) +
      c3*(  2.5*besselJ0 - 25.0*besselJ1 + 87.5*besselJ2);

  }

  /* normal exit */
  RETURN (status);
}


/*-----------------------------------------------------------------------
 *
 * NOTE: the following function should probably replaced by one that reads 
 * in data from a file.
 * (See documentation for a detailed description of siteParameters)
 *
 *----------------------------------------------------------------------
 */

static void GetSiteParameters ( SiteParameters*  siteParameters,
				IFOsite          siteID )
{

 
  switch (siteID) {
    case LHO: /* LH0 - LIGO Hanford observatory */
      siteParameters->siteID            = siteID;
      siteParameters->siteName          = siteNames[siteID];
      siteParameters->vertexNorth       = 46.45236;
      siteParameters->vertexWest        = 119.40753;
      siteParameters->arm1Orientation   = 36.8;
      siteParameters->arm2Orientation   = 126.8;
      siteParameters->armLengthCM       = 4.e5;
      break;
    case LLO: /* LL0 - LIGO Livingston observatory */
      siteParameters->siteID            = siteID;
      siteParameters->siteName          = siteNames[siteID];
      siteParameters->vertexNorth       = 30.56277;
      siteParameters->vertexWest        = 90.77425;
      siteParameters->arm1Orientation   = 108.0;
      siteParameters->arm2Orientation   = 198.0;
      siteParameters->armLengthCM       = 4.0e5;
      break;
    case VIRGO:
      siteParameters->siteID            = siteID;
      siteParameters->siteName          = siteNames[siteID];
      siteParameters->vertexNorth       = 43.6333;
      siteParameters->vertexWest        = -10.5;
      siteParameters->arm1Orientation   = 71.5;
      siteParameters->arm2Orientation   = 341.5;
      siteParameters->armLengthCM       = 3.0e5;
      break;
    case GEO600: 
      siteParameters->siteID            = siteID;
      siteParameters->siteName          = siteNames[siteID];
      siteParameters->vertexNorth       = 52.2467;
      siteParameters->vertexWest        = -9.82167;
      siteParameters->arm1Orientation   = 25.94;
      siteParameters->arm2Orientation   = 291.61;
      siteParameters->armLengthCM       = 6.0e4;
      break;
    case TAMA300:
      siteParameters->siteID            = siteID;
      siteParameters->siteName          = siteNames[siteID];
      siteParameters->vertexNorth       = 35.6766;
      siteParameters->vertexWest        = -139.536;
      siteParameters->arm1Orientation   = 90.0;
      siteParameters->arm2Orientation   = 180.0;
      siteParameters->armLengthCM       = 3.0e4;
      break;
   default: /* this case will never occur since an error would have 
	       already been generated if siteID did not have a proper 
	       enum types */
      break;
  }  

  return;
}

/*----------------------------------------------------------------------
 *
 * The following static function calculates the siteCoordinates for 
 * a given interferometer site from the siteParameters (see documentation
 * for more details about siteParameters and siteCoordinates).
 *
 *----------------------------------------------------------------------
 */

#define EQUATORIAL (6.37814e+08) /* earth's equatorial radius, in cm */
#define POLAR      (6.35676e+08) /* earth's polar radius, in cm */

static void CalcSiteCoordinates ( SiteCoordinates*   siteCoordinates,
				  SiteParameters*    siteParameters )
{

  INT4     i;

  REAL8    psi, phi;
  REAL8    phi1, phi2;
  REAL4    a2, b2;
  REAL4    denom;
  REAL4    north[3], east[3];

  const REAL4  convFactor = LAL_PI/180.0;

  /* convert angles from degrees to radians */
  psi  =      convFactor*siteParameters->vertexNorth;
  phi  = -1.0*convFactor*siteParameters->vertexWest;
  phi1 =      convFactor*siteParameters->arm1Orientation;
  phi2 =      convFactor*siteParameters->arm2Orientation;

  /* 
   Correct for oblateness of earth, use reference spheroid with   
   EQUATORIAL and POLAR radii given in table 4.2 of
   "Spacecraft attitude determination and control",  
   Ed. James R. Wertz, D. Reidel Publishing Co., Boston, 1978.
  */

  a2 = EQUATORIAL*EQUATORIAL;
  b2 = POLAR*POLAR;
  denom = sqrt( a2*cos(psi)*cos(psi) + b2*sin(psi)*sin(psi) );

  /* calculate cartesian coordinates of central station */
  (siteCoordinates->vertex)[0] = a2*cos(phi)*cos(psi)/denom;
  (siteCoordinates->vertex)[1] = a2*sin(phi)*cos(psi)/denom;
  (siteCoordinates->vertex)[2] = b2*sin(psi)/denom;

  /* calculate unit vectors in the north and east directions */ 
  /* note: north[] and east[] are orthogonal */
  north[0] = -1.0*cos(phi)*sin(psi);
  north[1] = -1.0*sin(phi)*sin(psi);
  north[2] = cos(psi);

  east[0] = -1.0*sin(phi);
  east[1] = cos(phi);
  east[2] = 0.0;

  /* calculate unit vectors that point along arms 1, 2 */
  /* note: arm1[] and arm2[] are NOT necessarily orthogonal */ 
  for ( i=0; i<3; i++ ) {
    (siteCoordinates->arm1)[i] = ( cos(phi1)*north[i] - sin(phi1)*east[i] );
    (siteCoordinates->arm2)[i] = ( cos(phi2)*north[i] - sin(phi2)*east[i] );
  }

  return;
}
