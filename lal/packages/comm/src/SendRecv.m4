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
 * \index{\texttt{LALMPISendMsg()}}
 * \index{\texttt{LALMPIRecvMsg()}}
 * \index{\texttt{LALMPISendCHAR()}}
 * \index{\texttt{LALMPIRecvCHAR()}}
 * \index{\texttt{LALMPISendINT2()}}
 * \index{\texttt{LALMPIRecvINT2()}}
 * \index{\texttt{LALMPISendINT4()}}
 * \index{\texttt{LALMPIRecvINT4()}}
 * \index{\texttt{LALMPISendINT8()}}
 * \index{\texttt{LALMPIRecvINT8()}}
 * \index{\texttt{LALMPISendUINT2()}}
 * \index{\texttt{LALMPIRecvUINT2()}}
 * \index{\texttt{LALMPISendUINT4()}}
 * \index{\texttt{LALMPIRecvUINT4()}}
 * \index{\texttt{LALMPISendUINT8()}}
 * \index{\texttt{LALMPIRecvUINT8()}}
 * \index{\texttt{LALMPISendREAL4()}}
 * \index{\texttt{LALMPIRecvREAL4()}}
 * \index{\texttt{LALMPISendREAL8()}}
 * \index{\texttt{LALMPIRecvREAL8()}}
 * \index{\texttt{LALMPISendCOMPLEX8()}}
 * \index{\texttt{LALMPIRecvCOMPLEX8()}}
 * \index{\texttt{LALMPISendCOMPLEX16()}}
 * \index{\texttt{LALMPIRecvCOMPLEX16()}}
 * \index{\texttt{LALMPISendCHARVector()}}
 * \index{\texttt{LALMPIRecvCHARVector()}}
 * \index{\texttt{LALMPISendINT2Vector()}}
 * \index{\texttt{LALMPIRecvINT2Vector()}}
 * \index{\texttt{LALMPISendINT4Vector()}}
 * \index{\texttt{LALMPIRecvINT4Vector()}}
 * \index{\texttt{LALMPISendINT8Vector()}}
 * \index{\texttt{LALMPIRecvINT8Vector()}}
 * \index{\texttt{LALMPISendUINT2Vector()}}
 * \index{\texttt{LALMPIRecvUINT2Vector()}}
 * \index{\texttt{LALMPISendUINT4Vector()}}
 * \index{\texttt{LALMPIRecvUINT4Vector()}}
 * \index{\texttt{LALMPISendUINT8Vector()}}
 * \index{\texttt{LALMPIRecvUINT8Vector()}}
 * \index{\texttt{LALMPISendREAL4Vector()}}
 * \index{\texttt{LALMPIRecvREAL4Vector()}}
 * \index{\texttt{LALMPISendREAL8Vector()}}
 * \index{\texttt{LALMPIRecvREAL8Vector()}}
 * \index{\texttt{LALMPISendCOMPLEX8Vector()}}
 * \index{\texttt{LALMPIRecvCOMPLEX8Vector()}}
 * \index{\texttt{LALMPISendCOMPLEX16Vector()}}
 * \index{\texttt{LALMPIRecvCOMPLEX16Vector()}}
 * \index{\texttt{LALMPISendINT2TimeSeries()}}
 * \index{\texttt{LALMPIRecvINT2TimeSeries()}}
 * \index{\texttt{LALMPISendINT4TimeSeries()}}
 * \index{\texttt{LALMPIRecvINT4TimeSeries()}}
 * \index{\texttt{LALMPISendINT8TimeSeries()}}
 * \index{\texttt{LALMPIRecvINT8TimeSeries()}}
 * \index{\texttt{LALMPISendUINT2TimeSeries()}}
 * \index{\texttt{LALMPIRecvUINT2TimeSeries()}}
 * \index{\texttt{LALMPISendUINT4TimeSeries()}}
 * \index{\texttt{LALMPIRecvUINT4TimeSeries()}}
 * \index{\texttt{LALMPISendUINT8TimeSeries()}}
 * \index{\texttt{LALMPIRecvUINT8TimeSeries()}}
 * \index{\texttt{LALMPISendREAL4TimeSeries()}}
 * \index{\texttt{LALMPIRecvREAL4TimeSeries()}}
 * \index{\texttt{LALMPISendREAL8TimeSeries()}}
 * \index{\texttt{LALMPIRecvREAL8TimeSeries()}}
 * \index{\texttt{LALMPISendCOMPLEX8TimeSeries()}}
 * \index{\texttt{LALMPIRecvCOMPLEX8TimeSeries()}}
 * \index{\texttt{LALMPISendCOMPLEX16TimeSeries()}}
 * \index{\texttt{LALMPIRecvCOMPLEX16TimeSeries()}}
 * \index{\texttt{LALMPISendINT2FrequencySeries()}}
 * \index{\texttt{LALMPIRecvINT2FrequencySeries()}}
 * \index{\texttt{LALMPISendINT4FrequencySeries()}}
 * \index{\texttt{LALMPIRecvINT4FrequencySeries()}}
 * \index{\texttt{LALMPISendINT8FrequencySeries()}}
 * \index{\texttt{LALMPIRecvINT8FrequencySeries()}}
 * \index{\texttt{LALMPISendUINT2FrequencySeries()}}
 * \index{\texttt{LALMPIRecvUINT2FrequencySeries()}}
 * \index{\texttt{LALMPISendUINT4FrequencySeries()}}
 * \index{\texttt{LALMPIRecvUINT4FrequencySeries()}}
 * \index{\texttt{LALMPISendUINT8FrequencySeries()}}
 * \index{\texttt{LALMPIRecvUINT8FrequencySeries()}}
 * \index{\texttt{LALMPISendREAL4FrequencySeries()}}
 * \index{\texttt{LALMPIRecvREAL4FrequencySeries()}}
 * \index{\texttt{LALMPISendREAL8FrequencySeries()}}
 * \index{\texttt{LALMPIRecvREAL8FrequencySeries()}}
 * \index{\texttt{LALMPISendCOMPLEX8FrequencySeries()}}
 * \index{\texttt{LALMPIRecvCOMPLEX8FrequencySeries()}}
 * \index{\texttt{LALMPISendCOMPLEX16FrequencySeries()}}
 * \index{\texttt{LALMPIRecvCOMPLEX16FrequencySeries()}}
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
