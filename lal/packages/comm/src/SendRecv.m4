/*
 * <lalVerbatim file="SendRecvCV">
 * $Id$
 * </lalVerbatim>
 */

/*
 * <lalLaTeX>
 *
 * \subsection{Module \texttt{SendRecv.c}}
 *
 * Routines to perform basic MPI sending and receiving of LAL data types.
 *
 * \subsubsection*{Prototypes}
 * \vspace{0.1in}
 * \input{SendRecvCP}
 * \idx{LALMPISendMsg()}
 * \idx{LALMPIRecvMsg()}
 * \idx{LALMPISendCHAR()}
 * \idx{LALMPIRecvCHAR()}
 * \idx{LALMPISendINT2()}
 * \idx{LALMPIRecvINT2()}
 * \idx{LALMPISendINT4()}
 * \idx{LALMPIRecvINT4()}
 * \idx{LALMPISendINT8()}
 * \idx{LALMPIRecvINT8()}
 * \idx{LALMPISendUINT2()}
 * \idx{LALMPIRecvUINT2()}
 * \idx{LALMPISendUINT4()}
 * \idx{LALMPIRecvUINT4()}
 * \idx{LALMPISendUINT8()}
 * \idx{LALMPIRecvUINT8()}
 * \idx{LALMPISendREAL4()}
 * \idx{LALMPIRecvREAL4()}
 * \idx{LALMPISendREAL8()}
 * \idx{LALMPIRecvREAL8()}
 * \idx{LALMPISendCOMPLEX8()}
 * \idx{LALMPIRecvCOMPLEX8()}
 * \idx{LALMPISendCOMPLEX16()}
 * \idx{LALMPIRecvCOMPLEX16()}
 * \idx{LALMPISendCHARVector()}
 * \idx{LALMPIRecvCHARVector()}
 * \idx{LALMPISendINT2Vector()}
 * \idx{LALMPIRecvINT2Vector()}
 * \idx{LALMPISendINT4Vector()}
 * \idx{LALMPIRecvINT4Vector()}
 * \idx{LALMPISendINT8Vector()}
 * \idx{LALMPIRecvINT8Vector()}
 * \idx{LALMPISendUINT2Vector()}
 * \idx{LALMPIRecvUINT2Vector()}
 * \idx{LALMPISendUINT4Vector()}
 * \idx{LALMPIRecvUINT4Vector()}
 * \idx{LALMPISendUINT8Vector()}
 * \idx{LALMPIRecvUINT8Vector()}
 * \idx{LALMPISendREAL4Vector()}
 * \idx{LALMPIRecvREAL4Vector()}
 * \idx{LALMPISendREAL8Vector()}
 * \idx{LALMPIRecvREAL8Vector()}
 * \idx{LALMPISendCOMPLEX8Vector()}
 * \idx{LALMPIRecvCOMPLEX8Vector()}
 * \idx{LALMPISendCOMPLEX16Vector()}
 * \idx{LALMPIRecvCOMPLEX16Vector()}
 * \idx{LALMPISendINT2TimeSeries()}
 * \idx{LALMPIRecvINT2TimeSeries()}
 * \idx{LALMPISendINT4TimeSeries()}
 * \idx{LALMPIRecvINT4TimeSeries()}
 * \idx{LALMPISendINT8TimeSeries()}
 * \idx{LALMPIRecvINT8TimeSeries()}
 * \idx{LALMPISendUINT2TimeSeries()}
 * \idx{LALMPIRecvUINT2TimeSeries()}
 * \idx{LALMPISendUINT4TimeSeries()}
 * \idx{LALMPIRecvUINT4TimeSeries()}
 * \idx{LALMPISendUINT8TimeSeries()}
 * \idx{LALMPIRecvUINT8TimeSeries()}
 * \idx{LALMPISendREAL4TimeSeries()}
 * \idx{LALMPIRecvREAL4TimeSeries()}
 * \idx{LALMPISendREAL8TimeSeries()}
 * \idx{LALMPIRecvREAL8TimeSeries()}
 * \idx{LALMPISendCOMPLEX8TimeSeries()}
 * \idx{LALMPIRecvCOMPLEX8TimeSeries()}
 * \idx{LALMPISendCOMPLEX16TimeSeries()}
 * \idx{LALMPIRecvCOMPLEX16TimeSeries()}
 * \idx{LALMPISendINT2FrequencySeries()}
 * \idx{LALMPIRecvINT2FrequencySeries()}
 * \idx{LALMPISendINT4FrequencySeries()}
 * \idx{LALMPIRecvINT4FrequencySeries()}
 * \idx{LALMPISendINT8FrequencySeries()}
 * \idx{LALMPIRecvINT8FrequencySeries()}
 * \idx{LALMPISendUINT2FrequencySeries()}
 * \idx{LALMPIRecvUINT2FrequencySeries()}
 * \idx{LALMPISendUINT4FrequencySeries()}
 * \idx{LALMPIRecvUINT4FrequencySeries()}
 * \idx{LALMPISendUINT8FrequencySeries()}
 * \idx{LALMPIRecvUINT8FrequencySeries()}
 * \idx{LALMPISendREAL4FrequencySeries()}
 * \idx{LALMPIRecvREAL4FrequencySeries()}
 * \idx{LALMPISendREAL8FrequencySeries()}
 * \idx{LALMPIRecvREAL8FrequencySeries()}
 * \idx{LALMPISendCOMPLEX8FrequencySeries()}
 * \idx{LALMPIRecvCOMPLEX8FrequencySeries()}
 * \idx{LALMPISendCOMPLEX16FrequencySeries()}
 * \idx{LALMPIRecvCOMPLEX16FrequencySeries()}
 *
 * \subsubsection*{Description}
 *
 * The routines \verb+LALMPISend+$\langle\textit{type}\rangle$ and
 * \verb+LALMPIRecv+$\langle\textit{type}\rangle$ must be used in pairs:
 * the sender must specify the MPI process number of the receiving process
 * as the input variable \verb+dest+, while the recieving process must specify
 * the MPI process number of the sending process as its input variable
 * \verb+source+.  The object that is sent/recived is the output for both
 * classes of functions.  The parameter is the MPI communicator \verb+comm+.
 *
 * The exception is the function \verb+LALMPIRecvMsg()+.  This function does
 * not specify a source---it receives from any process.  The function
 * \verb+LALMPISendMsg()+ is used to inform the reciever, which should be
 * listening with \verb+LALMPIRecvMsg()+, that the sender wishes to exchange
 * information.
 *
 * </lalLaTeX>
 */

#if 0 /* autodoc block */
<lalLaTeX>
The following is a simple example of a commissar-comrad situation.
 
\begin{verbatim}
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <lal/LALStdlib.h>
#include <lal/Comm.h>

INT4 lalDebugLevel = 0;

int Commissar( INT4 size, MPI_Comm comm );
int Comrad( INT4 size, INT4 rank, MPI_Comm comm );

int main( int argc, char *argv[] )
{
  int size, rank, code;

  MPI_Init( &argc, &argv );
  MPI_Comm_size( MPI_COMM_WORLD, &size );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  code = rank ? Comrad( size, rank, MPI_COMM_WORLD )
    : Commissar( size, MPI_COMM_WORLD );

  MPI_Finalize();
  return code;
}

int Commissar( INT4 size, MPI_Comm comm )
{
  static LALStatus status;
  MPIMessage msg;
  REAL4 data;

  while ( size > 1 ) /* while there are still comrads */
  {
    /* receive order from comrad */
    LALMPIRecvMsg( &status, &msg, comm );

    /* handle the message code */
    switch ( msg.msg )
    {
      case 0: /* comrad is quitting */
        --size;
        break;

      case 1: /* comrad is sending data */
        LALMPIRecvREAL4( &status, &data, msg.source, comm );
        printf( "data from comrad %d: %f", msg.source, data );
        break;

      default: /* should never happen */
        return 1;
    }
  }

  return 0;
}

int Comrad( INT4 size, INT4 rank, MPI_Comm comm )
{
  static LALStatus status;
  MPIMessage msg;
  REAL4 data = fmod( rank, size ); /* some junk data */

  /* send message to commissar indicating that comrad is ready to send data */
  msg.msg    = 1;    /* indicates that data is to be communicated           */
  msg.send   = 1;    /* indicates that comrad (local process) will send     */
  msg.source = rank; /* rank of local process                               */
  LALMPISendMsg( &status, &msg, 0, comm ); /* rank zero is commissar */

  /* commissar now knows data is coming so send it */
  LALMPISendREAL4( &status, &data, 0, comm );

  /* send message to commissar indicating that comrad is quitting */
  msg.msg    = 0;    /* indicates that comrad is quitting         */
  msg.source = rank; /* rank of local process                     */
  LALMPISendMsg( &status, &msg, 0, comm );

  return 0;
}
\end{verbatim}
 
</lalLaTeX>
#endif /* autodoc block */

/*
 * <lalLaTeX>
 * \vfill{\footnotesize\input{SendRecvCV}}
 * </lalLaTeX>
 */

#include <lal/LALStdlib.h>
#include <lal/Comm.h>

NRCSID( SENDRECVC, "$Id$" );


/* <lalVerbatim file="SendRecvCP"> */
void LALMPISendMsg( LALStatus *status, MPIMessage *msg, INT4 dest, MPI_Comm comm )
{ /* </lalVerbatim> */
  INT4 code;

  INITSTATUS( status, "LALMPISendMsg", SENDRECVC );

  ASSERT( msg, status, COMMH_ENULL, COMMH_MSGENULL );

  code = MPI_Send( msg, sizeof( *msg ), MPI_BYTE, dest, MPIMsg, comm );
  if ( code != MPI_SUCCESS )
  {
    ABORT( status, COMMH_EMPIE, COMMH_MSGEMPIE );
  }

  RETURN( status );
}


/* <lalVerbatim file="SendRecvCP"> */
void LALMPIRecvMsg( LALStatus *status, MPIMessage *msg, MPI_Comm comm )
{ /* </lalVerbatim> */
  MPI_Status stat;
  INT4       code;

  INITSTATUS( status, "LALMPIRecvMsg", SENDRECVC );

  ASSERT( msg, status, COMMH_ENULL, COMMH_MSGENULL );

  code = MPI_Recv( msg, sizeof( *msg ), MPI_BYTE, MPI_ANY_SOURCE, MPIMsg,
                   comm, &stat );
  if ( code != MPI_SUCCESS || stat.MPI_ERROR != MPI_SUCCESS )
  {
    ABORT( status, COMMH_EMPIE, COMMH_MSGEMPIE );
  }

  RETURN( status );
}


define(`TYPECODE',`Z')
include(`SendRecvAtomic.m4')

define(`TYPECODE',`C')
include(`SendRecvAtomic.m4')

define(`TYPECODE',`D')
include(`SendRecvAtomic.m4')

define(`TYPECODE',`S')
include(`SendRecvAtomic.m4')

define(`TYPECODE',`I2')
include(`SendRecvAtomic.m4')

define(`TYPECODE',`I4')
include(`SendRecvAtomic.m4')

define(`TYPECODE',`I8')
include(`SendRecvAtomic.m4')

define(`TYPECODE',`U2')
include(`SendRecvAtomic.m4')

define(`TYPECODE',`U4')
include(`SendRecvAtomic.m4')

define(`TYPECODE',`U8')
include(`SendRecvAtomic.m4')

define(`TYPECODE',`CHAR')
include(`SendRecvAtomic.m4')


define(`TYPECODE',`Z')
include(`SendRecvVector.m4')

define(`TYPECODE',`C')
include(`SendRecvVector.m4')

define(`TYPECODE',`D')
include(`SendRecvVector.m4')

define(`TYPECODE',`S')
include(`SendRecvVector.m4')

define(`TYPECODE',`I2')
include(`SendRecvVector.m4')

define(`TYPECODE',`I4')
include(`SendRecvVector.m4')

define(`TYPECODE',`I8')
include(`SendRecvVector.m4')

define(`TYPECODE',`U2')
include(`SendRecvVector.m4')

define(`TYPECODE',`U4')
include(`SendRecvVector.m4')

define(`TYPECODE',`U8')
include(`SendRecvVector.m4')

define(`TYPECODE',`CHAR')
include(`SendRecvVector.m4')


define(`SERIESCODE',`T')

define(`TYPECODE',`Z')
include(`SendRecvSeries.m4')

define(`TYPECODE',`C')
include(`SendRecvSeries.m4')

define(`TYPECODE',`D')
include(`SendRecvSeries.m4')

define(`TYPECODE',`S')
include(`SendRecvSeries.m4')

define(`TYPECODE',`I2')
include(`SendRecvSeries.m4')

define(`TYPECODE',`I4')
include(`SendRecvSeries.m4')

define(`TYPECODE',`I8')
include(`SendRecvSeries.m4')

define(`TYPECODE',`U2')
include(`SendRecvSeries.m4')

define(`TYPECODE',`U4')
include(`SendRecvSeries.m4')

define(`TYPECODE',`U8')
include(`SendRecvSeries.m4')


define(`SERIESCODE',`F')

define(`TYPECODE',`Z')
include(`SendRecvSeries.m4')

define(`TYPECODE',`C')
include(`SendRecvSeries.m4')

define(`TYPECODE',`D')
include(`SendRecvSeries.m4')

define(`TYPECODE',`S')
include(`SendRecvSeries.m4')

define(`TYPECODE',`I2')
include(`SendRecvSeries.m4')

define(`TYPECODE',`I4')
include(`SendRecvSeries.m4')

define(`TYPECODE',`I8')
include(`SendRecvSeries.m4')

define(`TYPECODE',`U2')
include(`SendRecvSeries.m4')

define(`TYPECODE',`U4')
include(`SendRecvSeries.m4')

define(`TYPECODE',`U8')
include(`SendRecvSeries.m4')
