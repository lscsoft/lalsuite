/*----------------------------------------------------------------------- 
 * 
 * File Name: FindChripTmpltTest.c
 *
 * Author: Brown, D. A.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#include <stdio.h>
#include <lal/LALStdlib.h>
#include <lal/FindChirpEngine.h>

NRCSID (MAIN, "$Id$");

int lalDebugLevel = 1;

static
void PrintInspiralBank (
    LALStatus           *status,
    InspiralTemplate    *head,
    FILE                *fp
    );

int
main (int argc, char *argv[])
{
  static LALStatus  stat;

  InspiralTemplate             *tmplt = NULL;

  InspiralCoarseBankIn         *bankIn;
  FindChirpCreateBankParams    *createBankParams;

  
  /*
   *
   * create the input and parameters
   *
   */


  bankIn = (InspiralCoarseBankIn *)
    LALCalloc( 1, sizeof(InspiralCoarseBankIn) );
  if ( ! bankIn )
  {
    fprintf( stderr, "Unable to allocate memory for bankIn\n" );
    return 1;
  }

  /* bank generation parameters */
  bankIn->mMin          = 1.0;
  bankIn->MMax          = 10.0;
  bankIn->mmCoarse      = 0.90;
  bankIn->mmFine        = 0.98;
  bankIn->fLower        = 40.;
  bankIn->fUpper        = 2000;
  bankIn->iflso         = 0;
  bankIn->tSampling     = 4000.;
  bankIn->NoisePsd      = LALLIGOIPsd;
  bankIn->method        = one;
  bankIn->order         = twoPN;
  bankIn->approximant   = taylor;
  bankIn->domain        = TimeDomain;
  bankIn->space         = Tau0Tau3;
  bankIn->etamin        = bankIn->mMin * ( bankIn->MMax - bankIn->mMin) /
    ( bankIn->MMax * bankIn->MMax );

  createBankParams = (FindChirpCreateBankParams *)
    LALCalloc( 1, sizeof(FindChirpCreateBankParams) );
  if ( ! createBankParams )
  {
    fprintf( stderr, "Unable to allocate memory for createBankParams\n" );
    return 1;
  }

  /* request a flat template bank */
  createBankParams->numLevel = 0;


  /*
   *
   * generate the template bank, print it out and then free it
   *
   */


  LALFindChirpCreateInspiralBank( &stat, bankIn, &tmplt, createBankParams );
  REPORTSTATUS( &stat );
  if ( stat.statusCode ) return stat.statusCode;

  PrintInspiralBank( &stat, tmplt, stdout ); 
  REPORTSTATUS( &stat );
  if ( stat.statusCode ) return stat.statusCode;

  LALFindChirpDestroyInspiralBank( &stat, &tmplt );
  REPORTSTATUS( &stat );
  if ( stat.statusCode ) return stat.statusCode;


  /*
   *
   * free memory, check for leaks and exit
   *
   */


  LALFree( bankIn );
  LALFree( createBankParams );

  LALCheckMemoryLeaks ();

  return 0;
}

static
void
PrintInspiralBank (
    LALStatus           *status,
    InspiralTemplate    *head,
    FILE                *fp
    )
{
  InspiralTemplate    *current;

  INITSTATUS( status, "PrintInspiralBank", MAIN );
  ATTATCHSTATUSPTR( status );

  ASSERT( head, status, 1, "Null Pointer" );


  /*
   *
   * print the inspiral template parameter bank
   *
   */


  for ( current = head; current; current = current->next)
  {
    fprintf( fp, 
        "%p\t%d\t%f\t%f\t%f\t%f\t%f\t%p\t%p\n",
        current, current->number, 
        current->mass1, current->mass2, current->totalMass,
        current->eta, current->mu,
        current->next, current->fine );
    fflush( fp );
    if ( current->fine )
    {
      PrintInspiralBank( status->statusPtr, current->fine, fp );
      CHECKSTATUSPTR( status );
    }
  }

  DETATCHSTATUSPTR( status );
  RETURN( status );
}
