/*
 * <lalVerbatim file="SendRecvTestCV">
 * $Id$
 * </lalVerbatim>
 */

/*
 * <lalLaTeX>
 *
 * \subsection{Program \texttt{SendRecvTest.c}}
 * \label{ss:SendRecvTest.c}
 *
 * Tests the LAL MPI send and receive commands.
 *
 * \subsubsection*{Usage}
 *
 * Example: run five processes on local machine using LAM
 * \begin{verbatim}
 * /bin/sh
 * echo `hostname` > lamhosts
 * rm -f schema
 * i=0; while [ $i -lt 5 ]; do echo SendRecvTest >> schema; i=`expr $i + 1`; done
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
 * mpirun -np 5 -machinefile machines SendRecvTest
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
 * \vfill{\footnotesize\input{SendRecvTestCV}}
 * </lalLaTeX>
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <mpi.h>
#include <lal/LALStdlib.h>
#include <lal/Comm.h>

/* the allowed message codes and the message strings */
enum { QuitMsg, DataMsg, RsltMsg };
const char *msgstr[] = { "QUIT", "DATA", "RSLT" };

void Commissar( LALStatus *status, MPI_Comm comm );
void Comrad( LALStatus *status, INT4 rank, MPI_Comm comm );

INT4 lalDebugLevel = LALMEMDBG;

int main( int argc, char *argv[] )
{
  LALStatus *status;
  int rank;

  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  status = LALCalloc( 1, sizeof( *status ) );
  if ( rank )
    Comrad( status, rank, MPI_COMM_WORLD );
  else
    Commissar( status, MPI_COMM_WORLD );

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

void Commissar( LALStatus *status, MPI_Comm comm )
{
  MPIMessage msg;
  INT4 data;
  INT4 rslt;
  int  size;
  
  INITSTATUS( status, "Commissar", "$Id$" );
  ATTATCHSTATUSPTR( status );

  MPI_Comm_size( comm, &size );

  while ( size > 1 ) /* while there are still comrads */
  {
    LALMPIRecvMsg( status->statusPtr, &msg, comm ); /* get order from Comrad */
    CHECKSTATUSPTR( status );
    fprintf( stderr, "  message %s from comrad %d... ",
        msgstr[msg.msg], msg.source );
    
    switch ( msg.msg ) /* handle the message code */
    {
      case QuitMsg: /* Comrad is quitting */
        fputs( "quitting\n", stderr );
        --size;
        break;

      case DataMsg: /* Commissar sends data to Comrad */
        srand( msg.source );
        data = rand();
        LALMPISendINT4( status->statusPtr, &data, msg.source, comm );
        CHECKSTATUSPTR( status );
        fputs( "data sent\n", stderr );
        break;

      case RsltMsg: /* Commissar receives result from Comrad */
        srand( msg.source );
        LALMPIRecvINT4( status->statusPtr, &rslt, msg.source, comm );
        CHECKSTATUSPTR( status );
        fputs( "result received\n", stderr );
        if ( rslt != rand() / msg.source )
        {
          ABORT( status, 2, "Incorrect result received" );
        }
        break;

      default: /* should never happen */
        fputs( "invalid message\n", stderr );
        ABORT( status, 1, "Invalid message received" );
    }
  }

  DETATCHSTATUSPTR( status );
  RETURN( status );
}

void Comrad( LALStatus *status, INT4 rank, MPI_Comm comm )
{
  MPIMessage msg;
  INT4 data;
  INT4 rslt;

  INITSTATUS( status, "Comrad", "$Id$" );
  ATTATCHSTATUSPTR( status );

  srand( rank );
  sleep( rand() % 5 ); /* sleep for up to four seconds */

  msg.source = rank;

  /* send message to Commissar requesting data */
  msg.msg  = DataMsg;
  msg.send = 0; /* Commissar will send the data */
  LALMPISendMsg( status->statusPtr, &msg, 0, comm );
  CHECKSTATUSPTR( status );

  /* Commissar will now send the data */
  LALMPIRecvINT4( status->statusPtr, &data, 0, comm );
  CHECKSTATUSPTR( status );

  sleep( rand() % 5 ); /* sleep for up to four seconds */
  rslt = data / rank;  /* compute result */

  /* send message to Commissar indicating the result is to be communicated */
  msg.msg  = RsltMsg;
  msg.send = 1; /* Comrad will send the result */
  LALMPISendMsg( status->statusPtr, &msg, 0, comm );
  CHECKSTATUSPTR( status );

  /* Commissar will now receive the result */
  LALMPISendINT4( status->statusPtr, &rslt, 0, comm );
  CHECKSTATUSPTR( status );

  sleep( rand() % 5 ); /* sleep for up to four seconds */

  /* send message to Commissar indicating that Comrad is quitting */
  msg.msg = QuitMsg;
  LALMPISendMsg( status->statusPtr, &msg, 0, comm );
  CHECKSTATUSPTR( status );

  DETATCHSTATUSPTR( status );
  RETURN( status );
}
