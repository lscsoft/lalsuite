/*
 * <lalVerbatim file="ExchangeTestCV">
 * $Id$
 * </lalVerbatim>
 */

/*
 * <lalLaTeX>
 *
 * \subsection{Program \texttt{ExchangeTest.c}}
 * \label{ss:ExchangeTest.c}
 *
 * Tests the LAL MPI exchange commands.
 *
 * \subsubsection*{Usage}
 *
 * Example: run five processes on local machine using LAM
 * \begin{verbatim}
 * /bin/sh
 * echo `hostname` > lamhosts
 * rm -f schema
 * i=0; while [ $i -lt 5 ]; do echo ExchangeTest >> schema; i=`expr $i + 1`; done
 * lamboot -v lamhosts
 * mpirun -v schema
 * wipe -v lamhosts
 * \end{verbatim}
 *
 * Example: run five processes on local machine using MPICH
 * \begin{verbatim}
 * /bin/sh
 * rm -f machines
 * i=0; while [ $i -lt 5 ]; do echo `hostname` >> schema; i=`expr $i + 1`; done
 * mpirun -np 5 -machinefile machines ExchangeTest
 * \end{verbatim}
 *
 * \subsubsection*{Exit codes}
 * \begin{tabular}{|c|l|}
 * \hline
 * Code  & Explanation           \\
 * \tt 0 & Success, normal exit. \\
 * \tt 1 & Failure.              \\
 * \hline
 * \end{tabular}
 *
 * \vfill{\footnotesize\input{ExchangeTestCV}}
 * </lalLaTeX>
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <mpi.h>
#include <lal/LALStdlib.h>
#include <lal/Comm.h>

/* the allowed exchange types and corresponding strings */
enum { ExchQuit, ExchData, ExchRslt };
const char *msgstr[] = { "QUIT", "DATA", "RSLT" };

INT4 lalDebugLevel = LALMEMDBG;

void Commissar( LALStatus *status, InitExchParams *params );
void Comrad( LALStatus *status, InitExchParams *params );

int main( int argc, char *argv[] )
{
  LALStatus *status;
  InitExchParams params;
  int rank;

  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  params.mpiComm   = MPI_COMM_WORLD;
  params.myProcNum = rank;

  status = LALCalloc( 1, sizeof( *status ) );
  if ( rank )
    Comrad( status, &params );
  else
    Commissar( status, &params );

  if ( status->statusCode )
  {
    REPORTSTATUS( status );
    return 1;
  }
  LALFree( status );

  LALCheckMemoryLeaks();
  MPI_Finalize();
  return 0;
}

void Commissar( LALStatus *status, InitExchParams *params )
{
  ExchParams *thisExch = NULL;
  int size;

  INITSTATUS( status, "Commissar", "$Id$" );
  ATTATCHSTATUSPTR( status );

  MPI_Comm_size( params->mpiComm, &size );
  while ( size > 1 )
  {
    INT4 obj;

    /* Commissar expects to be told the nature of the exchange */
    LALInitializeExchange( status->statusPtr, &thisExch, NULL, params );
    CHECKSTATUSPTR( status );
    fprintf( stderr, "  message %s from comrad %d... ",
        msgstr[thisExch->exchObjectType], thisExch->partnerProcNum );

    switch ( thisExch->exchObjectType ) /* hande the exchange type */
    {
      case ExchQuit: /* Comrad is quitting */
        fputs( "quitting\n", stderr );
        --size;
        break;

      case ExchData: /* Commissar requested to send data */
        srand( thisExch->partnerProcNum );
        for ( obj = 0; obj < thisExch->numObjects; ++obj )
        {
          INT4Vector *data = NULL;
          UINT4 i;

          /* create data for this Comrad */
          LALI4CreateVector( status->statusPtr, &data, thisExch->partnerProcNum );
          CHECKSTATUSPTR( status );
          for ( i = 0; i < data->length; ++i )
            data->data[i] = rand();

          /* send the data */
          LALExchangeINT4Vector( status->statusPtr, data, thisExch );
          CHECKSTATUSPTR( status );

          LALI4DestroyVector( status->statusPtr, &data );
          CHECKSTATUSPTR( status );
        }
        fputs( "data sent\n", stderr );
        break;

      case ExchRslt: /* Commissar receives results from Comrad */
        srand( thisExch->partnerProcNum );
        for ( obj = 0; obj < thisExch->numObjects; ++obj )
        {
          INT4 rslt;
          LALExchangeINT4( status->statusPtr, &rslt, thisExch );
          CHECKSTATUSPTR( status );
          if ( rslt != rand() / thisExch->partnerProcNum )
          {
            ABORT( status, 2, "Incorrect result recieved" );
          }
        }
        fputs( "results received\n", stderr );
        break;

      default: /* should never happen */
        fputs( "invalid message\n", stderr );
        ABORT( status, 1, "Invalid message received" );
    }

    LALFinalizeExchange( status->statusPtr, &thisExch );
  }

  DETATCHSTATUSPTR( status );
  RETURN( status );
}

void Comrad( LALStatus *status, InitExchParams *params )
{
  INT4Vector *data1    = NULL;
  INT4Vector *data2    = NULL;
  ExchParams *thisExch = NULL;
  ExchParams  exchQuit;
  ExchParams  exchData;
  ExchParams  exchRslt;
  UINT4 i;
  
  INITSTATUS( status, "Comrad", "$Id$" );
  ATTATCHSTATUSPTR( status );

  srand( params->myProcNum );
  sleep( rand() % 5 ); /* sleep for up to four seconds */

  /* all exchanges with Commissar */
  exchQuit.partnerProcNum=exchData.partnerProcNum=exchRslt.partnerProcNum = 0;

  /* set the exchange types */
  exchQuit.exchObjectType = ExchQuit;
  exchData.exchObjectType = ExchData;
  exchRslt.exchObjectType = ExchRslt;

  /* no objects exchanged with quit */
  exchQuit.numObjects = 0;
  exchQuit.send = 0; /* doesn't matter */

  /* Commissar will send two data objects */
  exchData.numObjects = 2;
  exchData.send = 0;

  /* Comrad will send a number of results equal to twice the Comrad's rank */
  exchRslt.numObjects = 2 * params->myProcNum;
  exchRslt.send = 1;

  /* create two data vectors with length equal to Comrad's rank */
  LALI4CreateVector( status->statusPtr, &data1, params->myProcNum );
  CHECKSTATUSPTR( status );
  LALI4CreateVector( status->statusPtr, &data2, params->myProcNum );
  CHECKSTATUSPTR( status );

  /* get the data */
  LALInitializeExchange( status->statusPtr, &thisExch, &exchData, params );
  CHECKSTATUSPTR( status );
  LALExchangeINT4Vector( status->statusPtr, data1, thisExch );
  CHECKSTATUSPTR( status );
  LALExchangeINT4Vector( status->statusPtr, data2, thisExch );
  CHECKSTATUSPTR( status );
  LALFinalizeExchange( status->statusPtr, &thisExch );
  CHECKSTATUSPTR( status );

  sleep( rand() % 5 ); /* sleep for up to four seconds */

  /* send the results */
  LALInitializeExchange( status->statusPtr, &thisExch, &exchRslt, params );
  CHECKSTATUSPTR( status );
  for ( i = 0; i < data1->length; ++i )
  {
    INT4 rslt = data1->data[i] / params->myProcNum;
    LALExchangeINT4( status->statusPtr, &rslt, thisExch );
    CHECKSTATUSPTR( status );
  }
  for ( i = 0; i < data2->length; ++i )
  {
    INT4 rslt = data2->data[i] / params->myProcNum;
    LALExchangeINT4( status->statusPtr, &rslt, thisExch );
    CHECKSTATUSPTR( status );
  }
  LALFinalizeExchange( status->statusPtr, &thisExch );
  CHECKSTATUSPTR( status );

  sleep( rand() % 5 ); /* sleep for up to four seconds */

  /* quit */
  LALInitializeExchange( status->statusPtr, &thisExch, &exchQuit, params );
  CHECKSTATUSPTR( status );
  LALFinalizeExchange( status->statusPtr, &thisExch );
  CHECKSTATUSPTR( status );

  LALI4DestroyVector( status->statusPtr, &data1 );
  CHECKSTATUSPTR( status );
  LALI4DestroyVector( status->statusPtr, &data2 );
  CHECKSTATUSPTR( status );

  DETATCHSTATUSPTR( status );
  RETURN( status );
}
