/*--------------------------------------------------------------------
 * void InitDopplerScanOnRefinedGrid()
 *
 *
 * [1] 
 * A helper function for metric validations.
 * At the same time, this is a wrapper for the Reinhard's InitDopplerScan().
 * This routine refine a usre-specified skyRegion so that the number of the 
 * grid points inside the region will be about a (hard-coded) desired value.
 *
 *
 * [2] Purpose 
 * When we search around a sky position, it may happen that we need to 
 * use numbers, "targetNumGridPoints", of sky grid points around the  
 * position. Taking the DopplerScanInit variable "scanInit", we try 
 * "maxTrial" times to refine the initial skyRegion "scanInit->skyRegion" 
 * so that the refined sky region contains "targetNumGridPoints +- tolerance"  
 * grid points. 
 * The output is a DopplerScanState variable "theScan".
 *
 * The input/output are same as InitDopplerScan to cope with it and NOT 
 * to introduce additional command line arguments to ComputeFStatistic. 
 * 
 * 
 * [3] Example usage: This is painful. 
 * #define NEARESTGRIDPOINTS_ON
 *
 * #ifdef NEARESTGRIDPOINTS_ON
 *     void InitDopplerScanOnRefinedGrid ( LALStatus *status, DopplerScanState *theScan, DopplerScanInit scanInit);
 * #endif
 *
 *
 *#ifdef NEARESTGRIDPOINTS_ON
 *   LAL_CALL ( InitDopplerScanOnRefinedGrid ( &status, &thisScan, &scanInit ), &status );
 *#else 
 *  if (lalDebugLevel) LALPrintError ("\nSetting up template grid ...");
 * LAL_CALL ( InitDopplerScan ( &status, &thisScan, &scanInit), &status); 
 *  if (lalDebugLevel) LALPrintError ("done.\n");
 * if ( uvar_outputSkyGrid ) {
 *   LALPrintError ("\nNow writing sky-grid into file '%s' ...", uvar_outputSkyGrid);
 *   LAL_CALL (writeSkyGridFile ( &status, thisScan.grid, uvar_outputSkyGrid, &scanInit), &status);
 *   LALPrintError (" done.\n\n");
 * }
 *#endif
 *
 *
 *
 *--------------------------------------------------------------------*/
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
NRCSID( NEARESTGRIDPOINTSC, "$Id$" );

#define NEARESTGRIDPOINTSH_ENULL 	1
#define NEARESTGRIDPOINTSH_ENONULL	2
#define NEARESTGRIDPOINTSH_ESUB		3
#define NEARESTGRIDPOINTSH_EMEM		4

#define NEARESTGRIDPOINTSH_MSGENULL 		"Arguments contained an unexpected null pointer"
#define NEARESTGRIDPOINTSH_MSGENONULL		"Output pointer is not NULL"
#define NEARESTGRIDPOINTSH_MSGESUB		"Failure in a subroutine call."
#define NEARESTGRIDPOINTSH_MSGEMEM		"Memory allocation error."


void RefineSkyRegion (LALStatus *status,  DopplerScanInit *scanInit, DopplerScanState *scanState,  SkyPosition *skyPosition, REAL8 *ratio);
void ComputeCenterOfMass (LALStatus *status,  SkyPosition *skyPosition, SkyRegion *skyRegion);
void InitDopplerScanOnRefinedGrid ( LALStatus *status, DopplerScanState *theScan, DopplerScanInit *scanInit);

void 
InitDopplerScanOnRefinedGrid ( LALStatus *status, 
			       DopplerScanState *theScan, /* output */ 
			       DopplerScanInit *scanInit   /* input */ )
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

  INITSTATUS( status, "InitDopplerScanOnRefinedGrid", NEARESTGRIDPOINTSC );
  ATTATCHSTATUSPTR( status );
  /* This traps coding errors in the calling routine. */
  ASSERT ( scanInit, status, NEARESTGRIDPOINTSH_ENULL , NEARESTGRIDPOINTSH_MSGENULL );  
  ASSERT ( theScan != NULL, status, NEARESTGRIDPOINTSH_ENONULL , NEARESTGRIDPOINTSH_MSGENONULL );  

  /* copy the origianl scanInit to a temporal scanInitTmp. */
  scanInitTmp = scanInit;
  tmpptr = (CHAR *) LALMalloc( sizeof(CHAR) * ( strlen(scanInit->skyRegion) + 1 ) );
  strcpy(tmpptr,scanInit->skyRegion);
  scanInitTmp->skyRegion = tmpptr; 

  TRY( InitDopplerScan ( status->statusPtr, theScan, scanInitTmp ), status );

  test = ( (INT4) theScan->numGridPoints ) - targetNumGridPoints;

  tmpSkyRegion.vertices = NULL; 
  TRY( ParseSkyRegion ( status->statusPtr, &tmpSkyRegion, scanInitTmp->skyRegion ), status );
  /* Compute a center of the mass, construct an initial sky polygon */
  TRY( ComputeCenterOfMass ( status->statusPtr, &centerOfMass, &tmpSkyRegion ), status );
  if ( status->statusCode ) {
    LALFree(tmpSkyRegion.vertices);
    REPORTSTATUS ( status );
    ABORT (status,  NEARESTGRIDPOINTSH_ESUB ,  NEARESTGRIDPOINTSH_MSGESUB );
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
  LALFree( scanInitTmp->skyRegion );
  DETATCHSTATUSPTR (status);
  RETURN( status );
} /* void InitDopplerScanOnRefinedGrid () */

/*--------------------------------------------------------------------
 * void ComputeCenterOfMass()
 *
 * This compute the center of the mass.
 * (1)Find the center of the mass point by vectorial-averaging over 
 *    all the vertices. 
 *    $\vec g = \frac{1}{N}\sum_{i}^N\vec p_i$
 *--------------------------------------------------------------------*/
void 
ComputeCenterOfMass ( LALStatus *status, 
		     SkyPosition *skyposCM, /* output center of mass of the vertices pointed by scanInit */
		     SkyRegion *skyRegion /* input Initial vertices. */)
{
  REAL8 longitudeCM = 0.0, latitudeCM = 0.0;
  UINT4 ic;

  INITSTATUS( status, "ComputeCenterOfMass", NEARESTGRIDPOINTSC );
  /* This traps coding errors in the calling routine. */
  ASSERT ( skyRegion, status, NEARESTGRIDPOINTSH_ENULL , NEARESTGRIDPOINTSH_MSGENULL );  
  ASSERT ( skyposCM != NULL, status, NEARESTGRIDPOINTSH_ENONULL , NEARESTGRIDPOINTSH_MSGENONULL );  


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
		 DopplerScanInit *scanInit, 
		 DopplerScanState *scanState, 
		 SkyPosition *centerOfMass,  
		 REAL8 *ratio )
{
  REAL8 longitudeCM = 0.0, latitudeCM = 0.0;
  REAL8 longitudeIth = 0.0, latitudeIth = 0.0;
  CHAR *tmpSkyRegion;
  CHAR *tmpstr;
  CHAR *ptr;
  UINT4 ic;
  /*
  REAL8 eps = sqrt(LAL_REAL4_EPS);
  */

  INITSTATUS( status, "RefineSkyRegion", NEARESTGRIDPOINTSC );

  longitudeCM = centerOfMass->longitude; 
  latitudeCM = centerOfMass->latitude; 

  tmpSkyRegion = (CHAR *) LALMalloc(strlen(scanInit->skyRegion)+512);
  if ( tmpSkyRegion == NULL ) {
    ABORT (status,  NEARESTGRIDPOINTSH_EMEM ,  NEARESTGRIDPOINTSH_MSGEMEM );
  }
  tmpstr = (CHAR *) LALMalloc(512);
  if ( tmpstr == NULL ) {
    LALFree(tmpSkyRegion);
    ABORT (status,  NEARESTGRIDPOINTSH_EMEM ,  NEARESTGRIDPOINTSH_MSGEMEM );
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


    /* This maps around the angular variable 
       modulo 0,2pi (alpha) and -pi/2,pi/2. Thanks to that (thanks 
       to Reinhard, I do not need the following restrictions!) */
    /*
    if( longitudeIth < 0.0 ) 
      longitudeIth = 0.0 + eps;
    if( longitudeIth > LAL_TWOPI ) 
      longitudeIth = LAL_TWOPI - eps;
    if( latitudeIth < - LAL_PI_2 ) 
      latitudeIth = - LAL_PI_2 + eps;
    if( latitudeIth > LAL_PI_2 ) 
      latitudeIth = LAL_PI_2 - eps;
    */

    sprintf ( tmpstr, "(%.16f, %.16f), ", longitudeIth, latitudeIth );
    strcat ( tmpSkyRegion, tmpstr );
  }
  /* ged rid of the last "," */
  ptr = strrchr ( tmpSkyRegion, ',' );
  *ptr = '\0';

  ptr =  LALRealloc(scanInit->skyRegion, sizeof(CHAR) * (strlen(tmpSkyRegion)+1) );
  if ( ptr == NULL ) {
    fprintf(stderr,"Unable to reallocate memory.\n");
    LALFree(scanInit->skyRegion);
    LALFree(tmpstr);
    LALFree(tmpSkyRegion);
    ABORT (status,  NEARESTGRIDPOINTSH_EMEM ,  NEARESTGRIDPOINTSH_MSGEMEM );
  }
  scanInit->skyRegion = ptr;
  /*
  memset(scanInit->skyRegion,'\0',strlen((scanInit->skyRegion)));
  */
  strcpy(scanInit->skyRegion,tmpSkyRegion);

  LALFree(tmpstr);
  LALFree(tmpSkyRegion);

  RETURN( status );
} /* void RefineSkyRegion() */


