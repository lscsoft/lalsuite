/*----------------------------------------------------------------------- 
 * 
 * File Name: findchirp.c
 *
 * Author: Brown, D. A.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#include "findchirp.h"

RCSID( "$Id$" );

extern int vrbflg;
UINT4  inspiralDebugFlag = 0;

UINT4  numPoints        = 1048576;
UINT4  numSegments      = 1;
UINT4  ovrlap           = 0;

INT4   invSpecTrunc     = 262144;
REAL4  fLow             = 40.0;
REAL4  dynRange         = 69.0;

INT4   numChisqBins     = 8;
REAL4  rhosqThreshVec[] = { 10.0, 0.0 };
REAL4  chisqThreshVec[] = { 5.0, 0.0 };

int main( int argc, char *argv[] )
{
  LALStatus             status = blank_status;

  FindChirpInitParams  *initParamsPtr   = NULL;
  FindChirpSlaveParams *slaveParamsPtr  = NULL;
  DataSegmentVector    *dataSegVec      = NULL;


  /*
   *
   * initialisation
   *
   */

  
  lal_errhandler = LAL_ERR_EXIT;
  set_debug_level( "1" );

  findchirp_parse_options( argc, argv );
  if ( vrbflg )
    findchirp_print_options();


  /*
   *
   * initialize parameters
   *
   */


  if ( ! ( initParamsPtr = (FindChirpInitParams *) 
      LALCalloc( 1, sizeof(FindChirpInitParams) ) ) )
  {
    fprintf( stderr, "Could not allocate memory for init params.\n" );
    exit( 1 );
  }
  initParamsPtr->numPoints      = numPoints;
  initParamsPtr->numSegments    = numSegments;
  initParamsPtr->numChisqBins   = numChisqBins;
  initParamsPtr->createRhosqVec = 0;  

  if ( ! ( slaveParamsPtr = (FindChirpSlaveParams *) 
      LALCalloc( 1, sizeof(FindChirpSlaveParams) ) ) )
  {
    fprintf( stderr, "Could not allocate memory for slave params.\n" );
    exit( 1 );
  }

  slaveParamsPtr->dataConditioned = 0;

  slaveParamsPtr->inspiralDebugFlagPtr = &inspiralDebugFlag;

  slaveParamsPtr->rhosqThreshVec = rhosqThreshVec;

  slaveParamsPtr->chisqThreshVec = chisqThreshVec;

  LAL_CALL( 
      LALCreateFindChirpSegmentVector( &status, &(slaveParamsPtr->fcSegVec), 
        initParamsPtr ), 
      &status );

  LAL_CALL( 
      LALFindChirpSPDataInit( &status, &(slaveParamsPtr->dataParams), 
        initParamsPtr ), 
      &status );
  slaveParamsPtr->dataParams->invSpecTrunc = invSpecTrunc;
  slaveParamsPtr->dataParams->fLow = fLow;
  slaveParamsPtr->dataParams->dynRange = pow( 2.0, dynRange );

  LAL_CALL( 
      LALFindChirpSPTemplateInit( &status, &(slaveParamsPtr->tmpltParams), 
        initParamsPtr ), 
      &status );
  slaveParamsPtr->tmpltParams->fLow = fLow;
  slaveParamsPtr->tmpltParams->dynRange = pow( 2.0, dynRange );

  LAL_CALL( 
      LALFindChirpFilterInit( &status, &(slaveParamsPtr->filterParams), 
        initParamsPtr ), 
      &status );
  slaveParamsPtr->filterParams->computeNegFreq = 0;

  LAL_CALL( 
      LALCreateFindChirpInput( &status, &(slaveParamsPtr->filterInput), 
        initParamsPtr ), 
      &status );
  LAL_CALL(
      LALFindChirpChisqVetoInit( &status, 
        slaveParamsPtr->filterParams->chisqParams, initParamsPtr->numChisqBins,
        initParamsPtr->numPoints ), 
      &status );


  /* simulation params go here */

  slaveParamsPtr->useMPI = 0;

  slaveParamsPtr->mpiComm = NULL;


  /*
   *
   * create and fill data segment vector
   *
   */


  LAL_CALL( 
      LALCreateDataSegmentVector( &status, &dataSegVec, initParamsPtr ), 
      &status );



  /*
   *
   * engine
   *
   */
  


  /*
   *
   * finalize parameters, check for memory leaks and exit
   *
   */
  

  LAL_CALL( 
      LALDestroyDataSegmentVector( &status, &dataSegVec ), 
      &status );

  LAL_CALL( 
      LALFindChirpChisqVetoFinalize( &status, 
        slaveParamsPtr->filterParams->chisqParams, 
        initParamsPtr->numChisqBins ), 
      &status );
  LAL_CALL( 
      LALDestroyFindChirpInput( &status, &(slaveParamsPtr->filterInput) ), 
      &status );

  LAL_CALL( 
      LALFindChirpFilterFinalize( &status, &(slaveParamsPtr->filterParams) ), 
      &status );

  LAL_CALL( 
      LALFindChirpSPTemplateFinalize( &status, &(slaveParamsPtr->tmpltParams) ),
      &status );

  LAL_CALL( 
      LALFindChirpSPDataFinalize( &status, &(slaveParamsPtr->dataParams) ),
      &status );

  LAL_CALL( 
      LALDestroyFindChirpSegmentVector( &status, &(slaveParamsPtr->fcSegVec) ),
      &status );


  LALFree( slaveParamsPtr );
  LALFree( initParamsPtr );
  
  LALCheckMemoryLeaks();

  exit( 0 );
}
