/*----------------------------------------------------------------------- 
 * 
 * File Name: TwoInterfFindChirpSPData.c
 *
 * Author: Bose, S.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/DataBuffer.h>
#include <lal/LALInspiral.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpSP.h>
#include <lal/TwoInterfFindChirp.h>
#include <lal/FindChirpChisq.h>

NRCSID (TWOINTERFFINDCHIRPSPDATAC, "$Id$");

#pragma <lalVerbatim file="TwoInterfFindChirpSPDataCP">
void
LALTwoInterfFindChirpSPDataInit (
    LALStatus                                 *status,
    TwoInterfFindChirpSPDataParamsVector     **vector,
    TwoInterfFindChirpInitParams              *params
    )
#pragma </lalVerbatim>
{
  UINT4                                        j;
  UINT4                                        k;
  REAL4                                       *amp;
  FindChirpSPDataParams                       *dataParamPtr;
  const REAL4                                  exponent = -7.0/6.0;
  TwoInterfFindChirpSPDataParamsVector        *vectorPtr;


  
  INITSTATUS( status, "LALTwoInterfFindChirpSPDataInit", TWOINTERFFINDCHIRPSPDATAC );
  ATTATCHSTATUSPTR( status );


  /*
   *
   * check that the arguments are reasonable
   *
   */


  /* make sure the output handle exists, but points to a null pointer */
  ASSERT( vector, status, TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
  ASSERT( !*vector, status, TWOINTERFFINDCHIRPH_ENNUL, TWOINTERFFINDCHIRPH_MSGENNUL );

  /* make sure that the parameter structure exists */
  ASSERT (params, status, TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL);

  /* make sure that the number of detectors and points in a segment 
     is positive */
  ASSERT (params->numDetectors > 0, status, 
      TWOINTERFFINDCHIRPH_ENUMZ, TWOINTERFFINDCHIRPH_MSGENUMZ);
  ASSERT (params->numPoints > 0, status, 
      TWOINTERFFINDCHIRPH_ENUMZ, TWOINTERFFINDCHIRPH_MSGENUMZ);

  /*
   *
   * allocate memory for the FindChirpSPDataParams
   *
   */


  /* create the output structure */
  vectorPtr = *vector = (TwoInterfFindChirpSPDataParamsVector *)
    LALCalloc( 1, sizeof(TwoInterfFindChirpSPDataParamsVector) );
  if ( ! vectorPtr )
  {
    ABORT( status, TWOINTERFFINDCHIRPH_EALOC, TWOINTERFFINDCHIRPH_MSGEALOC );
  }

  /* set the number of detectors in the vector */
  vectorPtr->length = params->numDetectors;
  
  /* create vector sub-structures */
  dataParamPtr = vectorPtr->data = (FindChirpSPDataParams *)
    LALCalloc( 1, vectorPtr->length*sizeof(FindChirpSPDataParams) );
  if ( ! dataParamPtr )
  {
    ABORT( status, TWOINTERFFINDCHIRPH_EALOC, TWOINTERFFINDCHIRPH_MSGEALOC );
  }
  
  
  /*
   *
   * allocate and fill vector for exponent of amplitude
   *
   */
  for (j = 0; j < vectorPtr->length; ++j)
    {
      LALCreateVector( status->statusPtr, &dataParamPtr[j].ampVec, 
		       (params->numPoints)/2 + 1 );
      BEGINFAIL( status )
	{
	  LALFree( dataParamPtr );
	  *vector = NULL;
	}
      ENDFAIL( status );
      
      amp = dataParamPtr[j].ampVec->data;
      amp[0] = 0.0;
      
      for ( k = 1; k < dataParamPtr[j].ampVec->length; ++k )
	amp[k] = pow( ((REAL4) k / (REAL4)params->numPoints), exponent );
      
      
      /*
       *
       * create fft plans and workspace vectors
       *
       */
      
      
      /* foward fft plan */
      LALCreateForwardRealFFTPlan( status->statusPtr, &dataParamPtr[j].fwdPlan, 
				   params->numPoints, 0 );
      BEGINFAIL( status )
	{
	  TRY( LALDestroyVector( status->statusPtr, &dataParamPtr[j].ampVec ), status );
	  
	  LALFree( dataParamPtr );
	  *vector = NULL;
	}
      ENDFAIL( status );
      
      /* inverse fft plan */
      LALCreateReverseRealFFTPlan( status->statusPtr, &dataParamPtr[j].invPlan, 
				   params->numPoints, 0 );
      BEGINFAIL( status )
	{
	  TRY( LALDestroyRealFFTPlan( status->statusPtr, &dataParamPtr[j].fwdPlan ), status );
	  TRY( LALDestroyVector( status->statusPtr, &dataParamPtr[j].ampVec ), status );
	  
	  LALFree( dataParamPtr );
	  *vector = NULL;
	}
      ENDFAIL( status );
      
      /* workspace vector w: time domain */
      LALCreateVector( status->statusPtr, &dataParamPtr[j].wVec, 
		       params->numPoints );
      BEGINFAIL( status )
	{
	  TRY( LALDestroyRealFFTPlan( status->statusPtr, &dataParamPtr[j].invPlan ), status ); 
	  TRY( LALDestroyRealFFTPlan( status->statusPtr, &dataParamPtr[j].fwdPlan ), status );
	  TRY( LALDestroyVector( status->statusPtr, &dataParamPtr[j].ampVec ), status );
	  
	  LALFree( dataParamPtr );
	  *vector = NULL;
	}
      ENDFAIL( status );
      
      /* workspace vector w: freq domain */
      LALCCreateVector( status->statusPtr, &dataParamPtr[j].wtildeVec, 
			params->numPoints/2 + 1 );
      BEGINFAIL( status )
	{
	  TRY( LALDestroyVector( status->statusPtr, &dataParamPtr[j].wVec ), status ); 
	  TRY( LALDestroyRealFFTPlan( status->statusPtr, &dataParamPtr[j].invPlan ), status ); 
	  TRY( LALDestroyRealFFTPlan( status->statusPtr, &dataParamPtr[j].fwdPlan ), status );
	  TRY( LALDestroyVector( status->statusPtr, &dataParamPtr[j].ampVec ), status );
	  
	  LALFree( dataParamPtr );
	  *vector = NULL;
	}
      ENDFAIL( status );
      CHECKSTATUSPTR (status);
      
      /* template power vector */
      LALCreateVector( status->statusPtr, &dataParamPtr[j].tmpltPowerVec, 
		       params->numPoints/2 + 1 );
      BEGINFAIL( status )
	{
	  TRY( LALCDestroyVector( status->statusPtr, &dataParamPtr[j].wtildeVec), status );
	  TRY( LALDestroyVector( status->statusPtr, &dataParamPtr[j].wVec ), status ); 
	  TRY( LALDestroyRealFFTPlan( status->statusPtr, &dataParamPtr[j].invPlan ), status ); 
	  TRY( LALDestroyRealFFTPlan( status->statusPtr, &dataParamPtr[j].fwdPlan ), status );
	  TRY( LALDestroyVector( status->statusPtr, &dataParamPtr[j].ampVec ), status );
	  
	  LALFree( dataParamPtr );
	  *vector = NULL;
	}
      ENDFAIL( status );
      
    }
  
  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}



#pragma <lalVerbatim file="TwoInterfFindChirpSPDataCP">
void
LALTwoInterfFindChirpSPDataFinalize (
    LALStatus                                 *status,
    TwoInterfFindChirpSPDataParamsVector     **vector
    )
#pragma </lalVerbatim>
{
  UINT4                                        i;
  TwoInterfFindChirpSPDataParamsVector        *vectorPtr;
  FindChirpSPDataParams                       *dataParamPtr;

  INITSTATUS( status, "LALTwoInterfFindChirpSPDataFinalize", TWOINTERFFINDCHIRPSPDATAC );
  ATTATCHSTATUSPTR( status );


  /*
   *
   * check that the arguments are reasonable
   *
   */
  

  /* make sure handle is non-null and points to a non-null pointer */
  ASSERT( vector, status, TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
  ASSERT( *vector, status, TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );

  /* local pointer to structure */
  vectorPtr = *vector;

  
  /*
   *
   * destroy fft plans and workspace vectors
   *
   */

  dataParamPtr = (*vector)->data;

  for ( i = 0 ; i < (*vector)->length ; ++i )
    {
 
      LALDestroyRealFFTPlan (status->statusPtr, &dataParamPtr[i].fwdPlan);
      CHECKSTATUSPTR (status);
      
      LALDestroyRealFFTPlan (status->statusPtr, &dataParamPtr[i].invPlan);
      CHECKSTATUSPTR (status);
      
      LALDestroyVector (status->statusPtr, &dataParamPtr[i].wVec);
      CHECKSTATUSPTR (status);
      
      LALCDestroyVector (status->statusPtr, &dataParamPtr[i].wtildeVec);
      CHECKSTATUSPTR (status);
      
      LALDestroyVector (status->statusPtr, &dataParamPtr[i].tmpltPowerVec);
      CHECKSTATUSPTR (status);
      
      
      /*
       *
       * destroy vector for exponent of amplitude
       *
       */
      
      
      LALDestroyVector (status->statusPtr, &dataParamPtr[i].ampVec);
      CHECKSTATUSPTR (status);
      
    }
     
  /*
   *
   * free memory for the FindChirpSPDataParams
   *
   */
  
  LALFree (dataParamPtr);
  LALFree (vectorPtr);
  *vector = NULL;

  /* normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}


#pragma <lalVerbatim file="TwoInterfFindChirpSPDataCP">
void
LALTwoInterfFindChirpSPData (
    LALStatus                                    *status, 
    TwoInterfFindChirpSegmentVector              *twoInterfFcSegVec, 
    TwoInterfDataSegmentVector                   *twoInterfDataSegVec, 
    TwoInterfFindChirpSPDataParamsVector         *twoInterfDataParamsVec)
#pragma </lalVerbatim>
{
  UINT4                                  n; 
  FindChirpSPDataParams                 *params[2];
  FindChirpSegmentVector                *fcSegVecPtr[2];
  DataSegmentVector                     *dataSegVecPtr[2];

  INITSTATUS( status, "LALTwoInterfFindChirpSPData", TWOINTERFFINDCHIRPSPDATAC );
  ATTATCHSTATUSPTR( status );

  
  /*
   *
   * make sure that the arguments are reasonable
   *
   */

  ASSERT( twoInterfFcSegVec, status, 
      TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL 
      ": twoInterfFcSegVec" );

  for ( n = 0; n < 2 ; ++n ) /*CHECK: no. of detectors hardwired to 2*/
    {
      params[n] = &(twoInterfDataParamsVec->data[n]);
      fcSegVecPtr[n] = &(twoInterfFcSegVec->data[n]);
      dataSegVecPtr[n] = &(twoInterfDataSegVec->data[n]);
      
      ASSERT( params[n] , status, 
	      TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL 
	      ": twoInterfFcSegVec->data" );
      
      ASSERT( fcSegVecPtr[n] != NULL, status, 
	      TWOINTERFFINDCHIRPH_ENNUL, TWOINTERFFINDCHIRPH_MSGENNUL 
	      ": twoInterfDataSegVec->data" );
      
      ASSERT( dataSegVecPtr[n] != NULL , status, TWOINTERFFINDCHIRPH_ENNUL, 
	      TWOINTERFFINDCHIRPH_MSGENNUL ": twoInterfDataParamsVec->data" );

      LALFindChirpSPData (status->statusPtr, fcSegVecPtr[n], dataSegVecPtr[n], params[n]);
      CHECKSTATUSPTR(status);
    }    
  
  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}
