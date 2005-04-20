/*!@file
   @author Itoh, Yousuke  
   $Id$
   @note 
   @brief A helper function for wide-sky-area search codes validation.

  [1] 
  A helper function for metric validations.
  At the same time, this is a wrapper for the Reinhard's InitDopplerScan().
  This routine refine a usre-specified skyRegion so that the number of the 
  grid points inside the region will be about a (hard-coded) desired value.
 
 
  [2] Purpose 
  When we search around a sky position, it may happen that we need to 
  use numbers, "targetNumGridPoints", of sky grid points around the  
  position. Taking the DopplerScanInit variable "scanInit", we try 
  "maxTrial" times to refine the initial skyRegion "scanInit->skyRegion" 
  so that the refined sky region contains "targetNumGridPoints +- tolerance"  
  grid points. 
  The output is a DopplerScanState variable "theScan".
 
  The input/output are same as InitDopplerScan to cope with it and NOT 
  to introduce additional command line arguments to ComputeFStatistic. 


 $Revision$

 @todo 
<ul>
   <li> Make this code take a set of parameters: 
        targetNumGridPoints, tolerance, maxTrial.
   <li> Currently the refining algorithm assumes flat geometry.
   Make it appropriate to the celestial sphere
</ul>  

*/

#include "DopplerScan.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>

#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALDatatypes.h>
#include <lal/LALConstants.h>
#include <lal/StringInput.h>

extern INT4 lalDebugLevel;
NRCSID( NEIGHBORINGGRIDPOINTSC, "$Id$" );

#define NEIGHBORINGGRIDPOINTSH_ENULL 	1
#define NEIGHBORINGGRIDPOINTSH_ENONULL	2
#define NEIGHBORINGGRIDPOINTSH_ESUB		3
#define NEIGHBORINGGRIDPOINTSH_EMEM		4

#define NEIGHBORINGGRIDPOINTSH_MSGENULL 		"Arguments contained an unexpected null pointer"
#define NEIGHBORINGGRIDPOINTSH_MSGENONULL		"Output pointer is not NULL"
#define NEIGHBORINGGRIDPOINTSH_MSGESUB		"Failure in a subroutine call."
#define NEIGHBORINGGRIDPOINTSH_MSGEMEM		"Memory allocation error."


/*! Computes ratio-scaled vertices  */
void RefineSkyRegion (LALStatus *status,  DopplerScanInit *scanInit, DopplerScanState *scanState,  SkyPosition *skyPosition, REAL8 *ratio);
/*! Compute the geometrical center of the given skyRegion */
void ComputeCenterOfMass (LALStatus *status,  SkyPosition *skyPosition, SkyRegion *skyRegion);
/*! 
refine the original skyRegion so that there are about "targetNumGridPoints" grid points within a new skyRegion. Return the new skyRegion 
*/
void InitDopplerScanOnRefinedGrid ( LALStatus *status, DopplerScanState *theScan, DopplerScanInit *scanInit);


/**
 * Refine the skyRegion to search only at neighboring grid points of the 
 * center of the original skyRegion. 
 */
void 
InitDopplerScanOnRefinedGrid ( LALStatus *status, 
			       DopplerScanState *theScan,  /*!<  Output DopplerScanState variable includes a refined skyRegion  */ 
			       DopplerScanInit *scanInit   /*!<  Input DopplerScanInit variable includes the original skyRegion */ )
{
  UINT4 ic;
  INT4 targetNumGridPoints = 40; 
  INT4 tolerance = 20;           
  INT4 trialCounter = 0, maxTrial = 10;
  REAL8 ratio = 0.75;
  REAL8 viscosity = 0.01;
  INT4 test = 0;
  CHAR *tmpptr;

  SkyRegion tmpSkyRegion;
  SkyPosition centerOfMass;
  DopplerScanInit *scanInitTmp;

  INITSTATUS( status, "InitDopplerScanOnRefinedGrid", NEIGHBORINGGRIDPOINTSC );
  ATTATCHSTATUSPTR( status );
  /* This traps coding errors in the calling routine. */
  ASSERT ( scanInit, status, NEIGHBORINGGRIDPOINTSH_ENULL , NEIGHBORINGGRIDPOINTSH_MSGENULL );  
  ASSERT ( theScan != NULL, status, NEIGHBORINGGRIDPOINTSH_ENONULL , NEIGHBORINGGRIDPOINTSH_MSGENONULL );  

  /* copy the origianl scanInit to a temporal scanInitTmp. */
  scanInitTmp = scanInit;
  tmpptr = (CHAR *) LALMalloc( sizeof(CHAR)*( strlen(scanInit->skyRegionString) + 1 ) );
  strcpy(tmpptr,scanInit->skyRegionString);
  scanInitTmp->skyRegionString = tmpptr; 

  TRY( InitDopplerScan ( status->statusPtr, theScan, scanInitTmp ), status );

  test = ( (INT4) theScan->numGridPoints ) - targetNumGridPoints;

  tmpSkyRegion.vertices = NULL; 
  TRY(ParseSkyRegionString(status->statusPtr,&tmpSkyRegion, scanInitTmp->skyRegionString),status);
  /* Compute a center of the mass, construct an initial sky polygon */
  TRY( ComputeCenterOfMass ( status->statusPtr, &centerOfMass, &tmpSkyRegion ), status );
  if ( status->statusCode ) {
    LALFree(tmpSkyRegion.vertices);
    REPORTSTATUS ( status );
    ABORT (status,  NEIGHBORINGGRIDPOINTSH_ESUB ,  NEIGHBORINGGRIDPOINTSH_MSGESUB );
  }
  LALFree(tmpSkyRegion.vertices);


  while ( ( abs(test) > tolerance ) && ( trialCounter < maxTrial ) )
    {
      /* Adaptive refinement of the ratio.
         If the initial polygon contains 0 grid points, 
	 we set a fixed ratio. */
      if ( theScan->numGridPoints == 0 ) {
	ratio = 2.0;
      } else {
	ratio = sqrt ( ((REAL8) targetNumGridPoints) / ((REAL8) theScan->numGridPoints) );
      }

	
      if( test < - tolerance ) {
	ratio += trialCounter * viscosity;
      } else {
	if ( ratio > trialCounter * viscosity )
	  ratio -= trialCounter * viscosity;
      }
	
      RefineSkyRegion ( status->statusPtr, scanInitTmp, theScan, &centerOfMass, &ratio );
      theScan->state = 2;
      BEGINFAIL( status )
	TRY( FreeDopplerScan ( status->statusPtr, theScan ), status );
      ENDFAIL( status );
      TRY( FreeDopplerScan ( status->statusPtr, theScan ), status );
      TRY( InitDopplerScan ( status->statusPtr, theScan, scanInitTmp), status );
      test = theScan->numGridPoints - targetNumGridPoints;
      trialCounter++;
    }
  

  if ( lalDebugLevel >= 3 ) {
    fprintf(stderr,"%10d\n",theScan->numGridPoints);
    for ( ic = 0; ic < (theScan->skyRegion.numVertices); ic++ ) {
      fprintf(stderr, "( %20.16f, %20.16f )\n",
	      theScan->skyRegion.vertices[ic].longitude,
	      theScan->skyRegion.vertices[ic].latitude);
    }
    fprintf(stderr, "\n");
  }


  /* clean up */
  LALFree( scanInitTmp->skyRegionString );
  DETATCHSTATUSPTR (status);
  RETURN( status );
} /* void InitDopplerScanOnRefinedGrid () */

/*--------------------------------------------------------------------
 * void ComputeCenterOfMass()
 *
 * This compute the center of the mass.
 * Find the center of the mass point by vectorial-averaging over 
 *    all the vertices. 
 *    $\vec g = \frac{1}{N}\sum_{i}^N\vec p_i$
 *--------------------------------------------------------------------*/
void 
ComputeCenterOfMass ( LALStatus *status, 
		     SkyPosition *skyposCM, /*!< output: center of mass of the vertices pointed by scanInit */
		     SkyRegion *skyRegion /*!< input: Initial vertices. */)
{
  REAL8 longitudeCM = 0.0, latitudeCM = 0.0;
  UINT4 ic;

  INITSTATUS( status, "ComputeCenterOfMass", NEIGHBORINGGRIDPOINTSC );
  /* This traps coding errors in the calling routine. */
  ASSERT ( skyRegion, status, NEIGHBORINGGRIDPOINTSH_ENULL , NEIGHBORINGGRIDPOINTSH_MSGENULL );  
  ASSERT ( skyposCM != NULL, status, NEIGHBORINGGRIDPOINTSH_ENONULL , NEIGHBORINGGRIDPOINTSH_MSGENONULL );  


  /* Find the center of mass */
  for ( ic = 0; ic < (skyRegion->numVertices); ic++ ) {
    longitudeCM += skyRegion->vertices[ic].longitude;
    latitudeCM  += skyRegion->vertices[ic].latitude;
  }

  longitudeCM /= skyRegion->numVertices;
  latitudeCM /= skyRegion->numVertices;

  skyposCM->system = COORDINATESYSTEM_EQUATORIAL;
  skyposCM->longitude = longitudeCM; 
  skyposCM->latitude = latitudeCM;

  RETURN( status );
} /* void ComputeCenterOfMass() */


/*--------------------------------------------------------------------
 *
 * void RefineSkyRegion()
 *
 * This routine computes ratio-scaled vertices.
 * (1)The cente of the mass is given and fixed. 
 * (2)Find the relative vectors of all the vertices w.r.t. 
 *    the center of the mass point.
 *    $\vec d_i = \vec p_i - \vec g$
 * (3)Multiply all the relative vectors by a user-specified ratio.
 *    $\vec l_i = ratio \vec d_i$
 * (4)Find a scaled positional vector by adding the scaled relative 
 *    vector to the center of the mass vector.
 *    $\vec q_i = \vec g + \vec l_i$.
 * (5)Check that all the vertices are within the nominal angular-range. 
 *--------------------------------------------------------------------*/
void 
RefineSkyRegion (LALStatus *status, 
		 DopplerScanInit *scanInit, /*!< output: includes scaled vertices */
		 DopplerScanState *scanState, /*!< input: includes skyPosition to be scaled */
		 SkyPosition *centerOfMass,  /*!< input: geometrical center of the skyRegion */
		 REAL8 *ratio /*!< input: ratio by which the skyRegion will be scaled */)
{
  REAL8 longitudeCM = 0.0, latitudeCM = 0.0;
  REAL8 longitudeIth = 0.0, latitudeIth = 0.0;
  CHAR *tmpSkyRegion;
  CHAR *tmpstr;
  CHAR *ptr;
  UINT4 ic;

  INITSTATUS( status, "RefineSkyRegion", NEIGHBORINGGRIDPOINTSC );

  longitudeCM = centerOfMass->longitude; 
  latitudeCM = centerOfMass->latitude; 

  tmpSkyRegion = (CHAR *) LALMalloc(strlen(scanInit->skyRegionString)+512);
  if ( tmpSkyRegion == NULL ) {
    ABORT (status,  NEIGHBORINGGRIDPOINTSH_EMEM ,  NEIGHBORINGGRIDPOINTSH_MSGEMEM );
  }
  tmpstr = (CHAR *) LALMalloc(512);
  if ( tmpstr == NULL ) {
    LALFree(tmpSkyRegion);
    ABORT (status,  NEIGHBORINGGRIDPOINTSH_EMEM ,  NEIGHBORINGGRIDPOINTSH_MSGEMEM );
  }
  memset(tmpSkyRegion,'\0',sizeof(tmpSkyRegion));
  memset(tmpstr,'\0',sizeof(tmpstr));

  /* Find the refined positional vector for the i-th verice. Then construct 
     CHAR STRING for the refined scanInit->skyRegion. */
  for ( ic = 0; ic < (scanState->skyRegion.numVertices); ic++ ) {
    longitudeIth = scanState->skyRegion.vertices[ic].longitude - longitudeCM;
    latitudeIth = scanState->skyRegion.vertices[ic].latitude - latitudeCM;

    longitudeIth *= (*ratio);
    latitudeIth *= (*ratio);

    longitudeIth += longitudeCM;
    latitudeIth += latitudeCM;


    sprintf ( tmpstr, "(%.16f, %.16f), ", longitudeIth, latitudeIth );
    strcat ( tmpSkyRegion, tmpstr );
  }
  /* ged rid of the last "," */
  ptr = strrchr ( tmpSkyRegion, ',' );
  *ptr = '\0';

  ptr =  LALRealloc(scanInit->skyRegionString, sizeof(CHAR) * (strlen(tmpSkyRegion)+1) );
  if ( ptr == NULL ) {
    fprintf(stderr,"Unable to reallocate memory.\n");
    LALFree(scanInit->skyRegionString);
    LALFree(tmpstr);
    LALFree(tmpSkyRegion);
    ABORT (status,  NEIGHBORINGGRIDPOINTSH_EMEM ,  NEIGHBORINGGRIDPOINTSH_MSGEMEM );
  }
  scanInit->skyRegionString = ptr;
  strcpy(scanInit->skyRegionString, tmpSkyRegion);

  LALFree(tmpstr);
  LALFree(tmpSkyRegion);

  RETURN( status );
} /* void RefineSkyRegion() */


