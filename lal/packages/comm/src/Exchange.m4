/*
 * <lalVerbatim file="ExchangeCV">
 * $Id$
 * </lalVerbatim>
 */

/*
 * <lalLaTeX>
 *
 * \subsection{Module \texttt{Exchange.c}}
 *
 * Routines to perform MPI exchanges of LAL data types.
 *
 * \subsubsection*{Prototypes}
 * \vspace{0.1in}
 * \input{ExchangeCP}
 * \idx{LALInitializeExchange()}
 * \idx{LALFinalizeExchange()}
 * \idx{LALExchangeCHAR()}
 * \idx{LALExchangeINT2()}
 * \idx{LALExchangeINT4()}
 * \idx{LALExchangeINT8()}
 * \idx{LALExchangeUINT2()}
 * \idx{LALExchangeUINT4()}
 * \idx{LALExchangeUINT8()}
 * \idx{LALExchangeREAL4()}
 * \idx{LALExchangeREAL8()}
 * \idx{LALExchangeCOMPLEX8()}
 * \idx{LALExchangeCOMPLEX16()}
 * \idx{LALExchangeCHARVector()}
 * \idx{LALExchangeINT2Vector()}
 * \idx{LALExchangeINT4Vector()}
 * \idx{LALExchangeINT8Vector()}
 * \idx{LALExchangeUINT2Vector()}
 * \idx{LALExchangeUINT4Vector()}
 * \idx{LALExchangeUINT8Vector()}
 * \idx{LALExchangeREAL4Vector()}
 * \idx{LALExchangeREAL8Vector()}
 * \idx{LALExchangeCOMPLEX8Vector()}
 * \idx{LALExchangeCOMPLEX16Vector()}
 * \idx{LALExchangeINT2TimeSeries()}
 * \idx{LALExchangeINT4TimeSeries()}
 * \idx{LALExchangeINT8TimeSeries()}
 * \idx{LALExchangeUINT2TimeSeries()}
 * \idx{LALExchangeUINT4TimeSeries()}
 * \idx{LALExchangeUINT8TimeSeries()}
 * \idx{LALExchangeREAL4TimeSeries()}
 * \idx{LALExchangeREAL8TimeSeries()}
 * \idx{LALExchangeCOMPLEX8TimeSeries()}
 * \idx{LALExchangeCOMPLEX16TimeSeries()}
 * \idx{LALExchangeINT2FrequencySeries()}
 * \idx{LALExchangeINT4FrequencySeries()}
 * \idx{LALExchangeINT8FrequencySeries()}
 * \idx{LALExchangeUINT2FrequencySeries()}
 * \idx{LALExchangeUINT4FrequencySeries()}
 * \idx{LALExchangeUINT8FrequencySeries()}
 * \idx{LALExchangeREAL4FrequencySeries()}
 * \idx{LALExchangeREAL8FrequencySeries()}
 * \idx{LALExchangeCOMPLEX8FrequencySeries()}
 * \idx{LALExchangeCOMPLEX16FrequencySeries()}
 *
 * \subsection*{Description}
 *
 * The routine \verb+LALInitializeExchange()+ is used to set up an exchange
 * protocol that contains the identities of the partners in the communication,
 * which one will be sending and which receiving, and the communicator.  The
 * protocol also contains additional information about the type of object that
 * is being exchanged and the number of objects to be exchanged.
 *
 * In initializing the exchange, there will be one active and one passive
 * partner.  The passive partner passes a \verb+NULL+ pointer as input to
 * \verb+LALInitializeExchange()+, which indicates that this partner is
 * listening for exchange requests.  The active partner sends the requested
 * protocol as input.  After \verb+LALInitializeExchange()+ returns, both
 * partners will have the agreed upon protocol allocated as output.  This
 * protocol is to be destroyed by \verb+LALFinalizeExchange()+ when all the
 * the communication for this exchange is completed.  The function
 * \verb+LALInitializeExchange()+ also requires a parameter structure that
 * contains the MPI communicator and the rank of the local process.
 *
 * The exchange protocol is used as the parameter to the various
 * \verb+LALExchange+$\langle\textit{type}\rangle$ functions.  Since the
 * exchange protocol knows who is the sender and who is the receiver, there
 * does not need to be a \verb+Send+-type and a \verb+Recv+-type function,
 * just the \verb+Exchange+-type.  If there are multiple objects to be
 * exchanged, the exchange function must be called the specified number of
 * times.
 *
 * </lalLaTeX>
 */

#if 0 /* autodoc block */
<lalLaTeX>
The following is an example with comrads which always initiate exchanges
with a commissar:

\begin{verbatim}
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <lal/LALStdlib.h>
#include <lal/Comm.h>

INT4 lalDebugLevel = 0;

int Commissar( INT4 size, InitExchParams params );
int Comrad( INT4 rank, InitExchParams params );

/* the types of exchanges */
enum { ExchData, ExchResults, ExchFinished };

int main( int argc, char *argv[] )
{
  InitExchParams params;
  int size, rank, code;

  MPI_Init( &argc, &argv );
  MPI_Comm_size( MPI_COMM_WORLD, &size );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  initParams.mpiComm   = comm;
  initParams.myProcNum = rank;

  code = rank ? Comrad( size, rank, &params ) : Commissar( size, &params );

  MPI_Finalize();
  return code;
}

int Commissar( INT4 size, InitExchParams *params )
{
  static LALStatus status;
  INT4 a = 0;
  INT4 b = 1;

  while ( size > 1 )
  {
    ExchParams *thisExch = NULL;
    INT4        obj;

    /* commissar expects to be told the nature of the exchange */
    LALInitializeExchange( &status, &thisExch, NULL, params );

    /* handle the exchange type */
    switch ( thisExch->exchObjectType )
    {
      case ExchData: /* commissar requested to send data */
        assert( thisExch->send == 1 );      /* commissar should send */
        for ( obj = 0; obj < thisExch->numObjects; ++obj )
        {
          INT4 fib = a + b; /* data is next Fibonacci number */
          LALExchangeINT4( &status, &fib, thisExch );
          a = b;
          b = fib;
        }
        break;

      case ExchResults: /* commissar receives results */
        assert( thisExch->send == 0 ); /* commissar should recieve */
        for ( obj = 0; obj < thisExch->numObjects; ++obj )
        {
          REAL4 result;
          LALExchangeREAL4( &status, &result, thisExch );
          printf( "result from slave %d: %f\n", thisExch->partnerProcNum, result );
        }
        break;

      case ExchFinished: /* comrad is finished */
        --size;
        break;

      default: /* this should not happen */
        return 1;
    }

    LALFinalizeExchange( &status, &thisExch );
  }

  return 0;
}

int Comrad( INT4 rank, InitExchParams *params )
{
  static LALStatus status;
  ExchParams *thisExch = NULL;
  ExchParams  exchData;
  ExchParams  exchResults;
  ExchParams  exchFinished;
  INT4 data[3];
  INT4 obj;

  /* define the various types of exchanges */

  exchData.exchObjectType = ExchData; /* data will be exchanged        */
  exchData.send = 0;                  /* commissar sends the data      */
  exchData.numObjects = 3;            /* always get three bits of data */
  exchData.partnerProcNum = 0;        /* partner is commissar          */

  exchResults.exchObjectType = ExchResults; /* results will be exchanged */
  exchResults.send = 1;                     /* comrad sends the data     */
  exchResults.numObjects = rank;            /* different for each comrad */
  exchResults.partnerProcNum = 0;           /* partner is commissar      */

  exchFinished.exchObjectType = ExchFinished; /* comrad is finished signal */
  exchFinished.send = 0;                      /* doesn't matter            */
  exchFinished.numObjects = 0;                /* no things are sent        */
  exchFinished.partnerProcNum = 0;            /* partner is commissar      */

  /* get some data */
  LALInitializeExchange( &status, &thisExch, &exchData, params );
  LALExchangeINT4( &status, &a, thisExch );
  LALExchangeINT4( &status, &b, thisExch );
  LALExchangeINT4( &status, &c, thisExch );
  LALFinalizeExchange( &status, &thisExch );

  /* send some results */
  LALInitializeExchange( &status, &thisExch, &exchResults, params );
  for ( obj = 0; obj < rank; ++obj )
  {
    REAL4 result = fabs( data[obj % 3], rank ); /* some stupid result */
    LALExchangeREAL4( &status, &result, thisExch );
  }
  LALFinalizeExchange( &status, &thisExch );

  /* tell commissar that comrad is finished */
  LALInitializeExchange( &status, &thisExch, &exchFinished, params );
  LALFinalizeExchange( &status, &thisExch );

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

NRCSID( EXCHANGEC, "$Id$" );

/* <lalVerbatim file="ExchangeCP"> */
void
LALInitializeExchange(
    LALStatus       *status,
    ExchParams     **exchParamsOut,
    ExchParams      *exchParamsInp,
    InitExchParams  *params
    )
{ /* </lalVerbatim> */
  MPIMessage hello; /* initialization message */

  INITSTATUS( status, "LALInitializeExchange", EXCHANGEC );
  ATTATCHSTATUSPTR( status );

  ASSERT( exchParamsOut,    status, COMMH_ENULL, COMMH_MSGENULL );
  ASSERT( ! *exchParamsOut, status, COMMH_ENNUL, COMMH_MSGENNUL );

  /* allocate memory for output exchange parameters */
  *exchParamsOut = (ExchParams *) LALCalloc( 1, sizeof( **exchParamsOut ) );
  if ( ! *exchParamsOut )
  {
    ABORT( status, COMMH_ENULL, COMMH_MSGENULL );
  }

  if ( exchParamsInp ) /* I am initializing the exchange */
  {
    INT4 dest = exchParamsInp->partnerProcNum;

    /* initialize communications */

    hello.msg    = exchParamsInp->exchObjectType;
    hello.source = params->myProcNum;

    /*
     * just use send field as a means to communicate number of objects:
     * set to negative to indicate that initializer (I) want to receive
     * those objects
     *
     * add one to the number of objects in case it is zero
     */

    if ( exchParamsInp->numObjects < 0 )
    {
      LALFree( *exchParamsOut );
      *exchParamsOut = NULL;
      ABORT( status, COMMH_ENOBJ, COMMH_MSGENOBJ );
    }

    if ( exchParamsInp->send )
    {
      hello.send = exchParamsInp->numObjects + 1;
    }
    else
    {
      hello.send = -(exchParamsInp->numObjects + 1);
    }

    /* send off the communications */
    LALMPISendMsg( status->statusPtr, &hello, dest, params->mpiComm );
    BEGINFAIL( status )
    {
      LALFree( *exchParamsOut );
      *exchParamsOut = NULL;
    }
    ENDFAIL( status );

    /* copy the input structure to the output structure */
    (*exchParamsOut)->exchObjectType = exchParamsInp->exchObjectType;
    (*exchParamsOut)->send           = exchParamsInp->send;
    (*exchParamsOut)->numObjects     = exchParamsInp->numObjects;
    (*exchParamsOut)->partnerProcNum = exchParamsInp->partnerProcNum;

    /* copy the communicator from the parameter structure */
    (*exchParamsOut)->myProcNum      = params->myProcNum;
    (*exchParamsOut)->mpiComm        = params->mpiComm;
  }
  else /* I am waiting for someone else to initialize the exchange */
  {
    /* wait for incoming message */
    LALMPIRecvMsg( status->statusPtr, &hello, params->mpiComm );
    BEGINFAIL( status )
    {
      LALFree( *exchParamsOut );
      *exchParamsOut = NULL;
    }
    ENDFAIL( status );

    /* the message contains all the information needed */
    (*exchParamsOut)->exchObjectType = hello.msg;
    (*exchParamsOut)->send           = (hello.send < 0);
    (*exchParamsOut)->numObjects     = abs(hello.send) - 1;
    (*exchParamsOut)->partnerProcNum = hello.source;

    /* copy the communicator from the parameter structure */
    (*exchParamsOut)->myProcNum      = params->myProcNum;
    (*exchParamsOut)->mpiComm        = params->mpiComm;
  }

  DETATCHSTATUSPTR( status );
  RETURN( status );
}

/* <lalVerbatim file="ExchangeCP"> */
void
LALFinalizeExchange(
    LALStatus   *status,
    ExchParams **exchParams
    )
{ /* </lalVerbatim> */
  INT2       magic = (INT2)0xA505; /* A SOS */
  INT2Vector goodbye;

  INITSTATUS( status, "LALFinalizeExchange", EXCHANGEC );
  ATTATCHSTATUSPTR( status );

  ASSERT( exchParams,  status, COMMH_ENULL, COMMH_MSGENULL );
  ASSERT( *exchParams, status, COMMH_ENULL, COMMH_MSGENULL );

  /* by convension the sending partner initializes the final handshake */

  if ( (*exchParams)->send )
  {
    INT4 dest = (*exchParams)->partnerProcNum;

    goodbye.length = 1;
    goodbye.data   = &magic;

    LALMPISendINT2Vector( status->statusPtr, &goodbye, dest, 
        (*exchParams)->mpiComm );
    CHECKSTATUSPTR( status );
  }
  else
  {
    INT4 source  = (*exchParams)->partnerProcNum;
    INT2 myMagic = 0;

    goodbye.length = 1;
    goodbye.data   = &myMagic;

    LALMPIRecvINT2Vector( status->statusPtr, &goodbye, source,
        (*exchParams)->mpiComm );
    CHECKSTATUSPTR( status );

    if ( goodbye.data[0] != magic )
    {
      ABORT( status, COMMH_EHAND, COMMH_MSGEHAND );
    }
  }

  /* empty memory */
  LALFree( *exchParams );
  *exchParams = NULL;

  DETATCHSTATUSPTR( status );
  RETURN( status );
}



define(`TYPE',`COMPLEX16')
include(`ExchangeObject.m4')

define(`TYPE',`COMPLEX8')
include(`ExchangeObject.m4')

define(`TYPE',`REAL8')
include(`ExchangeObject.m4')

define(`TYPE',`REAL4')
include(`ExchangeObject.m4')

define(`TYPE',`INT8')
include(`ExchangeObject.m4')

define(`TYPE',`INT4')
include(`ExchangeObject.m4')

define(`TYPE',`INT2')
include(`ExchangeObject.m4')

define(`TYPE',`UINT8')
include(`ExchangeObject.m4')

define(`TYPE',`UINT4')
include(`ExchangeObject.m4')

define(`TYPE',`UINT2')
include(`ExchangeObject.m4')

define(`TYPE',`CHAR')
include(`ExchangeObject.m4')



define(`TYPE',`COMPLEX16Vector')
include(`ExchangeObject.m4')

define(`TYPE',`COMPLEX8Vector')
include(`ExchangeObject.m4')

define(`TYPE',`REAL8Vector')
include(`ExchangeObject.m4')

define(`TYPE',`REAL4Vector')
include(`ExchangeObject.m4')

define(`TYPE',`INT8Vector')
include(`ExchangeObject.m4')

define(`TYPE',`INT4Vector')
include(`ExchangeObject.m4')

define(`TYPE',`INT2Vector')
include(`ExchangeObject.m4')

define(`TYPE',`UINT8Vector')
include(`ExchangeObject.m4')

define(`TYPE',`UINT4Vector')
include(`ExchangeObject.m4')

define(`TYPE',`UINT2Vector')
include(`ExchangeObject.m4')

define(`TYPE',`CHARVector')
include(`ExchangeObject.m4')



define(`TYPE',`COMPLEX16TimeSeries')
include(`ExchangeObject.m4')

define(`TYPE',`COMPLEX8TimeSeries')
include(`ExchangeObject.m4')

define(`TYPE',`REAL8TimeSeries')
include(`ExchangeObject.m4')

define(`TYPE',`REAL4TimeSeries')
include(`ExchangeObject.m4')

define(`TYPE',`INT8TimeSeries')
include(`ExchangeObject.m4')

define(`TYPE',`INT4TimeSeries')
include(`ExchangeObject.m4')

define(`TYPE',`INT2TimeSeries')
include(`ExchangeObject.m4')

define(`TYPE',`UINT8TimeSeries')
include(`ExchangeObject.m4')

define(`TYPE',`UINT4TimeSeries')
include(`ExchangeObject.m4')

define(`TYPE',`UINT2TimeSeries')
include(`ExchangeObject.m4')




define(`TYPE',`COMPLEX16FrequencySeries')
include(`ExchangeObject.m4')

define(`TYPE',`COMPLEX8FrequencySeries')
include(`ExchangeObject.m4')

define(`TYPE',`REAL8FrequencySeries')
include(`ExchangeObject.m4')

define(`TYPE',`REAL4FrequencySeries')
include(`ExchangeObject.m4')

define(`TYPE',`INT8FrequencySeries')
include(`ExchangeObject.m4')

define(`TYPE',`INT4FrequencySeries')
include(`ExchangeObject.m4')

define(`TYPE',`INT2FrequencySeries')
include(`ExchangeObject.m4')

define(`TYPE',`UINT8FrequencySeries')
include(`ExchangeObject.m4')

define(`TYPE',`UINT4FrequencySeries')
include(`ExchangeObject.m4')

define(`TYPE',`UINT2FrequencySeries')
include(`ExchangeObject.m4')
