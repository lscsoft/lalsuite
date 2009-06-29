/*
*  Copyright (C) 2007 David Churches, Duncan Brown, Jolien Creighton
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

int lalDebugLevel = 0;

static
void PrintInspiralBank (
    LALStatus           *status,
    InspiralTemplate    *head,
    FILE                *fp
    );

int
main ( void )
{
  FILE                         *fp = NULL;
  static LALStatus              stat;

  InspiralTemplate             *tmplt = NULL;

  InspiralCoarseBankIn         *bankIn;
  FindChirpCreateBankParams    *createBankParams;

  const UINT4                   numSegments = 7;
  const UINT4                   numPts = 262144;
  UINT4                         k;

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
  bankIn->massRange     = MinMaxComponentMass;
  bankIn->mMin          = 0.3;
  bankIn->mMax          = 0.8;
  bankIn->MMax          = bankIn->mMax * 2.;
  bankIn->mmCoarse      = 0.97;
  bankIn->mmFine        = 0.99;
  bankIn->fLower        = 40.;
  bankIn->fUpper        = 1024L;
  bankIn->iflso         = 0;
  bankIn->tSampling     = 4096L;
  bankIn->order         = twoPN;
  bankIn->approximant   = TaylorT1;
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

  /* request a heirarchical template bank */
  createBankParams->numLevel = 0;
  createBankParams->numSegments = numSegments;


  /*
   *
   * create a power spectrum for the bank
   *
   */


  memset( &(bankIn->shf), 0, sizeof(REAL8FrequencySeries) );
  bankIn->shf.f0 = 0.0;
  bankIn->shf.deltaF = bankIn->tSampling / numPts;

  LALDCreateVector( &stat, &(bankIn->shf.data), numPts / 2 + 1 );
  REPORTSTATUS( &stat );
  if ( stat.statusCode ) return stat.statusCode;

  fp = fopen( "spec.dat", "r" );

  if ( ! fp )
  {
    for( k = 0; k < bankIn->shf.data->length ; ++k )
    {
      REAL8 freq = (REAL8) k * bankIn->shf.deltaF;
      if ( freq < bankIn->fLower )
      {
        LALLIGOIPsd( NULL, bankIn->shf.data->data + k, bankIn->fLower );
      }
      else
      {
        LALLIGOIPsd( NULL, bankIn->shf.data->data + k, freq );
      }
    }
  }
  else
  {
    for( k = 0; k < bankIn->shf.data->length ; ++k )
    {
      REAL8 specpt;
      fscanf( fp, "%le\n", &specpt );
      specpt /= 9.0e-46;
      bankIn->shf.data->data[k] = specpt;
    }
    fclose( fp );
  }



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


  LALDDestroyVector( &stat, &(bankIn->shf.data) );
  REPORTSTATUS( &stat );
  if ( stat.statusCode ) return stat.statusCode;
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
        "%p\t%d\t%.16e\t%.16e\t%p\t%p\n",
        current, current->number,
        current->mass1, current->mass2,
        current->next, current->fine );
    /*
     * fprintf( fp, "   > " );
     * for ( i = 0; i < current->segmentIdVec->length; ++i )
     * {
     *  fprintf( fp, "segment[%d] = %d  ",
     *      i, current->segmentIdVec->data[i] );
     * }
     * fprintf( fp, "\n\n" );
     */
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

static void
graphREAL4 (
    REAL4      *array,
    INT4        n,
    INT4        spacing
           )
{
  FILE *fp;
  INT4 i;

  /* open a file for writing */
  if ( !(fp = fopen( "temp.graph", "w" )) )
  {
    printf( "couldn't open file\n" );
  }

  /* print data into the file */
  for ( i = 0; i < n; i++ )
    fprintf( fp, "%d\t%e\n", i, array[i * spacing] );

  /* close the file */
  fclose( fp );

  /* start up graphing program with data in the file */
  /* system( "xmgr temp.graph 1>/dev/null 2>&1 &" ); */

  return;
}

static void
graphREAL8 (
    REAL8      *array,
    INT4        n,
    INT4        spacing
           )
{
  FILE *fp;
  INT4 i;

  /* open a file for writing */
  if ( !(fp = fopen( "temp.graph", "w" )) )
  {
    printf( "couldn't open file\n" );
  }

  /* print data into the file */
  for ( i = 0; i < n; i++ )
    fprintf( fp, "%d\t%e\n", i, array[i * spacing] );

  /* close the file */
  fclose( fp );

  /* start up graphing program with data in the file */
  /* system( "xmgr temp.graph 1>/dev/null 2>&1 &" ); */

  return;
}
