/*----------------------------------------------------------------------- 
 * 
 * File Name: FindChirpTmplt.c
 *
 * Author: Brown, D. A.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#include <lal/FindChirpEngine.h>

NRCSID( FINDCHIRPTMPLTC, "$Id$" );


void
LALFindChirpCreateInspiralBank (
    LALStatus                  *status,
    InspiralCoarseBankIn       *coarseBankIn,
    InspiralTemplate          **bankHead,
    FindChirpCreateBankParams  *params
                               )
{

  InspiralTemplateList         *coarseList = NULL;
  InspiralTemplateList         *fineList = NULL;

  InspiralFineBankIn           *fineBankIn = NULL;

  InspiralTemplate             *tmpltPtr = NULL;
  InspiralTemplate             *fineTmpltPtr = NULL;

  INT4                          numFine;
  INT4                          tmpltCounter = 0;
  INT4                          i;

  INITSTATUS( status, "LALFindChirpCreateInspiralBank", FINDCHIRPTMPLTC );
  ATTATCHSTATUSPTR( status );


  /* 
   *
   * check arguments
   *
   */


  ASSERT( coarseBankIn, status, 
      FINDCHIRPENGINEH_ENULL, FINDCHIRPENGINEH_MSGENULL );
  ASSERT( bankHead, status, 
      FINDCHIRPENGINEH_ENULL, FINDCHIRPENGINEH_MSGENULL );
  ASSERT( params, status, 
      FINDCHIRPENGINEH_ENULL, FINDCHIRPENGINEH_MSGENULL );

  /* make sure we only request a flat or a tmplt bank one deep */
  if ( params->numLevel > 1 )
  {
    ABORT( status, FINDCHIRPENGINEH_ELVEL, FINDCHIRPENGINEH_MSGELVEL );
  }

  /* make sure that we are filtering a non zero number of segments */
  if ( params->numSegments <= 0 )
  {
    ABORT( status, FINDCHIRPENGINEH_ELVEL, FINDCHIRPENGINEH_MSGELVEL );
  }

  /* make sure this is always checked */
  if ( *bankHead ) 
  {
    ABORT( status, FINDCHIRPENGINEH_ENNUL, FINDCHIRPENGINEH_MSGENNUL );
  }


  /*
   *
   * create the coarse template bank and turn it into a linked list
   *
   */


  /* call the coarse bank generation package */
  LALInspiralCreateCoarseBank( status->statusPtr, &coarseList, 
      &(params->numCoarse), *coarseBankIn );
  BEGINFAIL( status )
  {
    LALFree( coarseList );
  }
  ENDFAIL( status );

  /* create the head of the bank... */
  tmpltPtr = *bankHead = (InspiralTemplate *) 
    LALCalloc( 1, sizeof(InspiralTemplate) );
  if ( ! bankHead )
  {
    LALFree( coarseList );
    ABORT( status, FINDCHIRPENGINEH_EALOC, FINDCHIRPENGINEH_MSGEALOC );
  }

  /* and copy the template data into it */
  memcpy( tmpltPtr, &(coarseList[0].params), 
      sizeof(InspiralTemplate) );
  tmpltPtr->number = tmpltCounter;
  tmpltPtr->level = 0;
  tmpltPtr->next = NULL;
  tmpltPtr->fine = NULL;

  /* create the list of segments to be filtered... */
  LALI4CreateVector( status->statusPtr, &(tmpltPtr->segmentIdVec), 
      params->numSegments );
  BEGINFAIL( status )
  {
    LALFree( coarseList );
    TRY( LALFindChirpDestroyInspiralBank( status->statusPtr, bankHead ),
        status );
  }
  ENDFAIL( status );

  /* ...and turn them all off */
  memset( tmpltPtr->segmentIdVec->data, 0, 
      tmpltPtr->segmentIdVec->length * sizeof(INT4) );

  /* ...and the rest of the bank */
  for ( tmpltCounter = 1; tmpltCounter < params->numCoarse; ++tmpltCounter )
  {
    tmpltPtr->next = (InspiralTemplate *) 
      LALCalloc( 1, sizeof(InspiralTemplate) );
    if ( ! tmpltPtr->next )
    {
      LALFree( coarseList );
      TRY( LALFindChirpDestroyInspiralBank( status->statusPtr, bankHead ),
          status );
      ABORT( status, FINDCHIRPENGINEH_EALOC, FINDCHIRPENGINEH_MSGEALOC );
    }
    tmpltPtr = tmpltPtr->next;

    /* copy the template data */
    memcpy( tmpltPtr, &(coarseList[tmpltCounter].params), 
        sizeof(InspiralTemplate) );
    tmpltPtr->number = tmpltCounter;
    tmpltPtr->level = 0;
    tmpltPtr->next = NULL;
    tmpltPtr->fine = NULL;

    /* create the list of segments to be filtered... */
    LALI4CreateVector( status->statusPtr, &(tmpltPtr->segmentIdVec), 
        params->numSegments );
    BEGINFAIL( status )
    {
      LALFree( coarseList );
      TRY( LALFindChirpDestroyInspiralBank( status->statusPtr, bankHead ),
          status );
    }
    ENDFAIL( status );

    /* ...and turn them all off */
    memset( tmpltPtr->segmentIdVec->data, 0, 
        tmpltPtr->segmentIdVec->length * sizeof(INT4) );
  }


  /*
   *
   * create a fine bank around each coarse template, if requested
   *
   */


  if ( params->numLevel > 0 )
  {
    /* create fine bank generation input and copy the coarse data into it */
    fineBankIn = (InspiralFineBankIn *)
      LALCalloc( 1, sizeof(InspiralFineBankIn) );
    if ( ! fineBankIn )
    {
      LALFree( coarseList );
      TRY( LALFindChirpDestroyInspiralBank( status->statusPtr, bankHead ),
          status );
      ABORT( status, FINDCHIRPENGINEH_EALOC, FINDCHIRPENGINEH_MSGEALOC );
    }
    memcpy( &(fineBankIn->coarseIn), coarseBankIn, sizeof(InspiralCoarseBankIn) );

    /* start with the first template */
    tmpltPtr = *bankHead;

    /* proceed through the linked list of coarse templates */
    while( tmpltPtr )
    {
      /* copy the template list structure into the input */
      memcpy( &(fineBankIn->templateList), coarseList + tmpltPtr->number,
          sizeof(InspiralTemplateList) );

      /* try and create the fine bank about this coarse template */
      LALInspiralCreateFineBank( status->statusPtr, 
          &fineList, &numFine, *fineBankIn );
      BEGINFAIL( status )
      {
        LALFree( coarseList );
        LALFree( fineList );
        LALFree( fineBankIn );
        TRY( LALFindChirpDestroyInspiralBank( status->statusPtr, bankHead ),
            status );
      }
      ENDFAIL( status );

      /* if any fine templates are returned, attach them to the bank */
      if ( numFine )
      {
        /* attach the first fine template to the coarse template */
        tmpltPtr->fine = fineTmpltPtr = (InspiralTemplate *)
          LALCalloc( 1, sizeof(InspiralTemplate) );
        if ( ! fineTmpltPtr )
        {
          LALFree( coarseList );
          LALFree( fineList );
          LALFree( fineBankIn );
          TRY( LALFindChirpDestroyInspiralBank( status->statusPtr, bankHead ),
              status );
          ABORT( status, FINDCHIRPENGINEH_EALOC, FINDCHIRPENGINEH_MSGEALOC );
        }

        /* and copy the template data into it */
        memcpy( fineTmpltPtr, &(fineList[0].params), sizeof(InspiralTemplate) );
        fineTmpltPtr->number = tmpltCounter++;
        fineTmpltPtr->level = 1;
        fineTmpltPtr->next = NULL;
        fineTmpltPtr->fine = NULL;

        /* create the list of segments to be filtered... */
        LALI4CreateVector( status->statusPtr, &(fineTmpltPtr->segmentIdVec), 
            params->numSegments );
        BEGINFAIL( status )
        {
          LALFree( coarseList );
          LALFree( fineList );
          LALFree( fineBankIn );
          TRY( LALFindChirpDestroyInspiralBank( status->statusPtr, bankHead ),
              status );
        }
        ENDFAIL( status );

        /* ...and turn them all off */
        memset( fineTmpltPtr->segmentIdVec->data, 0, 
            tmpltPtr->segmentIdVec->length * sizeof(INT4) );
         
        /* ...and the rest of the fine bank */
        for ( i = 1; i < numFine; ++i )
        {
          fineTmpltPtr->next = (InspiralTemplate *) 
            LALCalloc( 1, sizeof(InspiralTemplate) );
          if ( ! fineTmpltPtr->next )
          {
            LALFree( coarseList );
            LALFree( fineList );
            LALFree( fineBankIn );
            TRY( LALFindChirpDestroyInspiralBank( status->statusPtr, bankHead ),
                status );
            ABORT( status, FINDCHIRPENGINEH_EALOC, FINDCHIRPENGINEH_MSGEALOC );
          }
          fineTmpltPtr = fineTmpltPtr->next;

          /* and copy the template data into it */
          memcpy( fineTmpltPtr, &(fineList[i].params), sizeof(InspiralTemplate) );
          fineTmpltPtr->number = tmpltCounter++;
          fineTmpltPtr->level = 1;
          fineTmpltPtr->next = NULL;
          fineTmpltPtr->fine = NULL;

          /* create the list of segments to be filtered... */
          LALI4CreateVector( status->statusPtr, &(fineTmpltPtr->segmentIdVec), 
              params->numSegments );
          BEGINFAIL( status )
          {
            LALFree( coarseList );
            LALFree( fineList );
            LALFree( fineBankIn );
            TRY( LALFindChirpDestroyInspiralBank( status->statusPtr, bankHead ),
                status );
          }
          ENDFAIL( status );

          /* ...and turn them all off */
          memset( fineTmpltPtr->segmentIdVec->data, 0, 
              tmpltPtr->segmentIdVec->length * sizeof(INT4) );
        }

      } /* end if any fine templates returned */

      /* go to the next coarse template */
      tmpltPtr = tmpltPtr->next;
    }

    /* free fine bank generation */
    LALFree( fineBankIn );
    LALFree( fineList );
  }


  /* 
   *
   * destroy the allocated memory 
   *
   */


  LALFree( coarseList );


  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}


void
LALFindChirpDestroyInspiralBank (
    LALStatus           *status,
    InspiralTemplate   **bankHead
                                )
{
  InspiralTemplate    *current;

  INITSTATUS( status, "LALFindChirpDestroyInspiralBank", FINDCHIRPTMPLTC );
  ATTATCHSTATUSPTR( status );

  ASSERT( *bankHead, status, FINDCHIRPENGINEH_ENULL, FINDCHIRPENGINEH_MSGENULL );

  /* destroy the inspiral template parameter bank */
  while ( *bankHead )
  {
    current = *bankHead;
    if ( current->fine )
    {
      LALFindChirpDestroyInspiralBank( status->statusPtr, &(current->fine) );
      CHECKSTATUSPTR( status );
    }
    *bankHead = (*bankHead)->next;
    LALI4DestroyVector( status->statusPtr, &(current->segmentIdVec) );
    CHECKSTATUSPTR( status );

    LALFree( current );
  }

  DETATCHSTATUSPTR( status );
  RETURN( status );
}

