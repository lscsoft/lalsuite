/*
 * <lalVerbatim file="CommHV">
 * $Id$
 * </lalVerbatim>
 */
   
/*
 * <lalLaTeX>
 *
 * \section{Header \texttt{Comm.h}}
 *
 * Provides routines for MPI communication.
 *
 * \subsection*{Synopsis}
 * \begin{verbatim}
 * #include <lal/Comm.h>
 * \end{verbatim}
 *
 * \noindent This header covers the routines for doing MPI communication.
 *
 * </lalLaTeX>
 */

#ifndef _COMM_H
#define _COMM_H

#include <mpi.h>
#include <lal/LALDatatypes.h>
#include <lal/AVFactories.h>

#ifdef  __cplusplus
extern "C" {
#endif

NRCSID (COMMH, "$Id$");

/*
 * <lalLaTeX>
 * \subsection*{Error conditions}
 * \input{CommHErrTab}
 * </lalLaTeX>
 *
 * <lalErrTable file="CommHErrTab">
 */
#define COMMH_ENULL 1
#define COMMH_ENNUL 2
#define COMMH_ESIZE 4
#define COMMH_ESZMM 8
#define COMMH_EMPIE 16
#define COMMH_EHAND 32
#define COMMH_ENOBJ 64

#define COMMH_MSGENULL "Null pointer"
#define COMMH_MSGENNUL "Non-Null pointer"
#define COMMH_MSGESIZE "Invalid size"
#define COMMH_MSGESZMM "Exchange size mismatch"
#define COMMH_MSGEMPIE "MPI error"
#define COMMH_MSGEHAND "Wrong handshake"
#define COMMH_MSGENOBJ "Invalid number of objects"
/*
 * </lalErrTable>
 */

/*
 * <lalLaTeX>
 * \subsection*{Types}
 * </lalLaTeX>
 */

typedef enum
{
  MPIDone,
  MPIMisc,
  MPIErr,
  MPIMsg,
  MPICHAR,
  MPII2,
  MPII4,
  MPII8,
  MPIU2,
  MPIU4,
  MPIU8,
  MPIS,
  MPID,
  MPIC,
  MPIZ,
  MPICHARVector,
  MPICHARVectorData,
  MPII2Vector,
  MPII2VectorData,
  MPII4Vector,
  MPII4VectorData,
  MPII8Vector,
  MPII8VectorData,
  MPIU2Vector,
  MPIU2VectorData,
  MPIU4Vector,
  MPIU4VectorData,
  MPIU8Vector,
  MPIU8VectorData,
  MPISVector,
  MPISVectorData,
  MPIDVector,
  MPIDVectorData,
  MPICVector,
  MPICVectorData,
  MPIZVector,
  MPIZVectorData,
  MPII2TimeSeries,
  MPII4TimeSeries,
  MPII8TimeSeries,
  MPIU2TimeSeries,
  MPIU4TimeSeries,
  MPIU8TimeSeries,
  MPISTimeSeries,
  MPIDTimeSeries,
  MPICTimeSeries,
  MPIZTimeSeries,
  MPII2FrequencySeries,
  MPII4FrequencySeries,
  MPII8FrequencySeries,
  MPIU2FrequencySeries,
  MPIU4FrequencySeries,
  MPIU8FrequencySeries,
  MPISFrequencySeries,
  MPIDFrequencySeries,
  MPICFrequencySeries,
  MPIZFrequencySeries
}
MPIMsgCode;

typedef struct
tagMPIMessage
{
  INT4 msg;
  INT4 send;
  INT4 source;
}
MPIMessage;

/*
 * <lalLaTeX>
 *
 * \subsubsection*{Structure \texttt{MPIMessage}}
 * \idx[Type]{MPIMessage}
 *
 * This structure is sent to a remote process, via \verb+LALMPISendMsg()+,
 * to alert that process that there is a message.  Note that
 * \verb+LALMPIRecvMsg()+ is the only \verb+Recv+-type function that does not 
 * require the source to be identified; the receiver can then identify the
 * source from the message received.
 *
 * Essentially, \verb+LALMPISendMsg()+ and \verb+LALMPIRecvMsg()+ form a
 * handshake for a subsequent transmission, and \verb+MPIMessage+ specifies
 * the protocol.  The local process uses \verb+LALMPISendMsg()+ to communicate
 * with a remote process, which is waiting to hear from \emph{any} process
 * using \verb+LALMPIRecvMsg()+.  The message the local process sends specifies
 * \begin{enumerate}
 *   \item An integer code telling the remote process what operation it should
 *     take (e.g., get ready to exchange some data, tell the remote process
 *     to terminate, etc.).
 *   \item A boolean integer that is zero if the remote process is expected
 *     to send something to the local process, or non-zero if the local process
 *     will send something to the remote process.
 *   \item An integer representing the MPI process number of the local process.
 * \end{enumerate}
 *
 * The fields are:
 * \begin{description}
 * \item[\texttt{INT4 msg}] An integer code specifying to the receiver what
 *   type of operation is to be taken.
 * \item[\texttt{INT4 send}] A boolean that is non-zero if the originator
 *   of the message will be sending something to the recipiant of the message
 *   in a subsequent communication, or zero if the recipiant of the message
 *   is expected to send something to the originator of the message.
 * \item[\texttt{INT4 source}] The MPI process number of the originator of the
 *   message.
 * \end{description}
 *
 * </lalLaTeX>
 */

typedef struct
tagExchParams
{
  INT4           send;
  INT4           numObjects;
  INT4           partnerProcNum;
  INT4           myProcNum;
  INT4           exchObjectType;
  MPI_Comm       mpiComm;
}
ExchParams;

typedef struct
tagInitExchParams
{
  INT4           myProcNum;
  MPI_Comm       mpiComm;
}
InitExchParams;

/*
 * <lalLaTeX>
 *
 * \subsubsection*{Structures \texttt{ExchParams} and \texttt{InitExchParams}}
 * \idx[Type]{ExchParams}
 * \idx[Type]{InitExchParams}
 *
 * These structures are used in the \verb+Exch+-type routines.  The structure
 * \verb+InitExchParams+ are the parameters used in initializing an exchange
 * protocol using \verb+LALInitializeExchange+.  The fields are:
 * \begin{description}
 *   \item[\texttt{INT4 myProcNum}] The MPI process number of the local process.
 *   \item[\texttt{MPI\_Comm mpiComm}] The MPI communicator.
 * \end{description}
 *
 * The structure \verb+ExchParams+ is created by \verb+LALInitializeExchange()+,
 * destroyed by \verb+LALFinalizeExchange()+, and serves as the parameter for
 * the various \verb+LALExchange+$\langle\textit{type}\rangle$ functions.  It
 * is also required as the input to \verb+LALInitializeExchange()+ for the
 * originator of the exchange request.  The fields are:
 * \begin{description}
 *   \item[\texttt{INT4 send}] A code that indicates whether this process is
 *     sending (non-zero) or recieving (zero).  The process that initializes
 *     the exchange chooses whether it will send in the subsequent exchanges
 *     (non-zero value for the \verb+send+ field of the input
 *     \verb+ExchParams+), or receive (zero value).
 *   \item[\texttt{INT4 numObjects}] The (maximum) number of objects to be
 *     exchanged.  (The partners in the exchange may have some mechanism to
 *     decide to terminate the exchange early, e.g., by exchanging a negative
 *     integer.)
 *   \item[\texttt{INT4 partnerProcNum}] The MPI process number of the partner
 *     in the exchange.
 *   \item[\texttt{INT4 myProcNum}] The MPI process number of the local process.
 *   \item[\texttt{INT4 exchObjectType}] An integer code representing the type
 *     of object that will be excanged.
 *   \item[\texttt{MPI\_Comm mpiComm}] The MPI communicator.
 * \end{description}
 *
 * </lalLaTeX>
 */


/*
 * <lalLaTeX>
 * \vfill{\footnotesize\input{CommHV}}
 * \newpage\input{SendRecvC}
 * \newpage\input{ExchangeC}
 * \newpage\input{SendRecvTestC}
 * \newpage\input{ExchangeTestC}
 * </lalLaTeX>
 */

void
LALMPISendMsg(
    LALStatus  *status,
    MPIMessage *msg,
    INT4        dest,
    MPI_Comm    mpiComm
    );

void
LALMPIRecvMsg(
    LALStatus  *status,
    MPIMessage *msg,
    MPI_Comm    mpiComm
    );



void
LALMPISendCHAR(
    LALStatus  *status,
    CHAR       *element,
    INT4        dest,
    MPI_Comm    mpiComm
    );

void
LALMPIRecvCHAR(
    LALStatus  *status,
    CHAR       *element,
    INT4        dest,
    MPI_Comm    mpiComm
    );

void
LALMPISendINT2(
    LALStatus  *status,
    INT2       *element,
    INT4        dest,
    MPI_Comm    mpiComm
    );

void
LALMPIRecvINT2(
    LALStatus  *status,
    INT2       *element,
    INT4        source,
    MPI_Comm    mpiComm
    );

void
LALMPISendINT4(
    LALStatus  *status,
    INT4       *element,
    INT4        dest,
    MPI_Comm    mpiComm
    );

void
LALMPIRecvINT4(
    LALStatus  *status,
    INT4       *element,
    INT4        source,
    MPI_Comm    mpiComm
    );

void
LALMPISendINT8(
    LALStatus  *status,
    INT8       *element,
    INT4        dest,
    MPI_Comm    mpiComm
    );

void
LALMPIRecvINT8(
    LALStatus  *status,
    INT8       *element,
    INT4        source,
    MPI_Comm    mpiComm
    );

void
LALMPISendUINT2(
    LALStatus   *status,
    UINT2       *element,
    INT4         dest,
    MPI_Comm     mpiComm
    );

void
LALMPIRecvUINT2(
    LALStatus   *status,
    UINT2       *element,
    INT4         source,
    MPI_Comm     mpiComm
    );

void
LALMPISendUINT4(
    LALStatus   *status,
    UINT4       *element,
    INT4         dest,
    MPI_Comm     mpiComm
    );

void
LALMPIRecvUINT4(
    LALStatus   *status,
    UINT4       *element,
    INT4         source,
    MPI_Comm     mpiComm
    );

void
LALMPISendUINT8(
    LALStatus   *status,
    UINT8       *element,
    INT4         dest,
    MPI_Comm     mpiComm
    );

void
LALMPIRecvUINT8(
    LALStatus   *status,
    UINT8       *element,
    INT4         source,
    MPI_Comm     mpiComm
    );

void
LALMPISendREAL4(
    LALStatus   *status,
    REAL4       *element,
    INT4         dest,
    MPI_Comm     mpiComm
    );

void
LALMPIRecvREAL4(
    LALStatus   *status,
    REAL4       *element,
    INT4         source,
    MPI_Comm     mpiComm
    );

void
LALMPISendREAL8(
    LALStatus   *status,
    REAL8       *element,
    INT4         dest,
    MPI_Comm     mpiComm
    );

void
LALMPIRecvREAL8(
    LALStatus   *status,
    REAL8       *element,
    INT4         source,
    MPI_Comm     mpiComm
    );

void
LALMPISendCOMPLEX8(
    LALStatus      *status,
    COMPLEX8       *element,
    INT4            dest,
    MPI_Comm        mpiComm
    );

void
LALMPIRecvCOMPLEX8(
    LALStatus      *status,
    COMPLEX8       *element,
    INT4            source,
    MPI_Comm        mpiComm
    );

void
LALMPISendCOMPLEX16(
    LALStatus       *status,
    COMPLEX16       *element,
    INT4             dest,
    MPI_Comm         mpiComm
    );

void
LALMPIRecvCOMPLEX16(
    LALStatus       *status,
    COMPLEX16       *element,
    INT4             source,
    MPI_Comm         mpiComm
    );



void
LALMPISendCHARVector(
    LALStatus  *status,
    CHARVector *vector,
    INT4        dest,
    MPI_Comm    mpiComm
    );

void
LALMPIRecvCHARVector(
    LALStatus  *status,
    CHARVector *vector,
    INT4        source,
    MPI_Comm    mpiComm
    );

void
LALMPISendINT2Vector(
    LALStatus  *status,
    INT2Vector *vector,
    INT4        dest,
    MPI_Comm    mpiComm
    );

void
LALMPIRecvINT2Vector(
    LALStatus  *status,
    INT2Vector *vector,
    INT4        source,
    MPI_Comm    mpiComm
    );

void
LALMPISendINT4Vector(
    LALStatus  *status,
    INT4Vector *vector,
    INT4        dest,
    MPI_Comm    mpiComm
    );

void
LALMPIRecvINT4Vector(
    LALStatus  *status,
    INT4Vector *vector,
    INT4        source,
    MPI_Comm    mpiComm
    );

void
LALMPISendINT8Vector(
    LALStatus  *status,
    INT8Vector *vector,
    INT4        dest,
    MPI_Comm    mpiComm
    );

void
LALMPIRecvINT8Vector(
    LALStatus  *status,
    INT8Vector *vector,
    INT4        source,
    MPI_Comm    mpiComm
    );

void
LALMPISendUINT2Vector(
    LALStatus   *status,
    UINT2Vector *vector,
    INT4         dest,
    MPI_Comm     mpiComm
    );

void
LALMPIRecvUINT2Vector(
    LALStatus   *status,
    UINT2Vector *vector,
    INT4         source,
    MPI_Comm     mpiComm
    );

void
LALMPISendUINT4Vector(
    LALStatus   *status,
    UINT4Vector *vector,
    INT4         dest,
    MPI_Comm     mpiComm
    );

void
LALMPIRecvUINT4Vector(
    LALStatus   *status,
    UINT4Vector *vector,
    INT4         source,
    MPI_Comm     mpiComm
    );

void
LALMPISendUINT8Vector(
    LALStatus   *status,
    UINT8Vector *vector,
    INT4         dest,
    MPI_Comm     mpiComm
    );

void
LALMPIRecvUINT8Vector(
    LALStatus   *status,
    UINT8Vector *vector,
    INT4         source,
    MPI_Comm     mpiComm
    );

void
LALMPISendREAL4Vector(
    LALStatus   *status,
    REAL4Vector *vector,
    INT4         dest,
    MPI_Comm     mpiComm
    );

void
LALMPIRecvREAL4Vector(
    LALStatus   *status,
    REAL4Vector *vector,
    INT4         source,
    MPI_Comm     mpiComm
    );

void
LALMPISendREAL8Vector(
    LALStatus   *status,
    REAL8Vector *vector,
    INT4         dest,
    MPI_Comm     mpiComm
    );

void
LALMPIRecvREAL8Vector(
    LALStatus   *status,
    REAL8Vector *vector,
    INT4         source,
    MPI_Comm     mpiComm
    );

void
LALMPISendCOMPLEX8Vector(
    LALStatus      *status,
    COMPLEX8Vector *vector,
    INT4            dest,
    MPI_Comm        mpiComm
    );

void
LALMPIRecvCOMPLEX8Vector(
    LALStatus      *status,
    COMPLEX8Vector *vector,
    INT4            source,
    MPI_Comm        mpiComm
    );

void
LALMPISendCOMPLEX16Vector(
    LALStatus       *status,
    COMPLEX16Vector *vector,
    INT4             dest,
    MPI_Comm         mpiComm
    );

void
LALMPIRecvCOMPLEX16Vector(
    LALStatus       *status,
    COMPLEX16Vector *vector,
    INT4             source,
    MPI_Comm         mpiComm
    );






void
LALMPISendINT2TimeSeries(
    LALStatus  *status,
    INT2TimeSeries *series,
    INT4        dest,
    MPI_Comm    mpiComm
    );

void
LALMPIRecvINT2TimeSeries(
    LALStatus  *status,
    INT2TimeSeries *series,
    INT4        source,
    MPI_Comm    mpiComm
    );

void
LALMPISendINT4TimeSeries(
    LALStatus  *status,
    INT4TimeSeries *series,
    INT4        dest,
    MPI_Comm    mpiComm
    );

void
LALMPIRecvINT4TimeSeries(
    LALStatus  *status,
    INT4TimeSeries *series,
    INT4        source,
    MPI_Comm    mpiComm
    );

void
LALMPISendINT8TimeSeries(
    LALStatus  *status,
    INT8TimeSeries *series,
    INT4        dest,
    MPI_Comm    mpiComm
    );

void
LALMPIRecvINT8TimeSeries(
    LALStatus  *status,
    INT8TimeSeries *series,
    INT4        source,
    MPI_Comm    mpiComm
    );

void
LALMPISendUINT2TimeSeries(
    LALStatus   *status,
    UINT2TimeSeries *series,
    INT4         dest,
    MPI_Comm     mpiComm
    );

void
LALMPIRecvUINT2TimeSeries(
    LALStatus   *status,
    UINT2TimeSeries *series,
    INT4         source,
    MPI_Comm     mpiComm
    );

void
LALMPISendUINT4TimeSeries(
    LALStatus   *status,
    UINT4TimeSeries *series,
    INT4         dest,
    MPI_Comm     mpiComm
    );

void
LALMPIRecvUINT4TimeSeries(
    LALStatus   *status,
    UINT4TimeSeries *series,
    INT4         source,
    MPI_Comm     mpiComm
    );

void
LALMPISendUINT8TimeSeries(
    LALStatus   *status,
    UINT8TimeSeries *series,
    INT4         dest,
    MPI_Comm     mpiComm
    );

void
LALMPIRecvUINT8TimeSeries(
    LALStatus   *status,
    UINT8TimeSeries *series,
    INT4         source,
    MPI_Comm     mpiComm
    );

void
LALMPISendREAL4TimeSeries(
    LALStatus   *status,
    REAL4TimeSeries *series,
    INT4         dest,
    MPI_Comm     mpiComm
    );

void
LALMPIRecvREAL4TimeSeries(
    LALStatus   *status,
    REAL4TimeSeries *series,
    INT4         source,
    MPI_Comm     mpiComm
    );

void
LALMPISendREAL8TimeSeries(
    LALStatus   *status,
    REAL8TimeSeries *series,
    INT4         dest,
    MPI_Comm     mpiComm
    );

void
LALMPIRecvREAL8TimeSeries(
    LALStatus   *status,
    REAL8TimeSeries *series,
    INT4         source,
    MPI_Comm     mpiComm
    );

void
LALMPISendCOMPLEX8TimeSeries(
    LALStatus      *status,
    COMPLEX8TimeSeries *series,
    INT4            dest,
    MPI_Comm        mpiComm
    );

void
LALMPIRecvCOMPLEX8TimeSeries(
    LALStatus      *status,
    COMPLEX8TimeSeries *series,
    INT4            source,
    MPI_Comm        mpiComm
    );

void
LALMPISendCOMPLEX16TimeSeries(
    LALStatus       *status,
    COMPLEX16TimeSeries *series,
    INT4             dest,
    MPI_Comm         mpiComm
    );

void
LALMPIRecvCOMPLEX16TimeSeries(
    LALStatus       *status,
    COMPLEX16TimeSeries *series,
    INT4             source,
    MPI_Comm         mpiComm
    );






void
LALMPISendINT2FrequencySeries(
    LALStatus  *status,
    INT2FrequencySeries *series,
    INT4        dest,
    MPI_Comm    mpiComm
    );

void
LALMPIRecvINT2FrequencySeries(
    LALStatus  *status,
    INT2FrequencySeries *series,
    INT4        source,
    MPI_Comm    mpiComm
    );

void
LALMPISendINT4FrequencySeries(
    LALStatus  *status,
    INT4FrequencySeries *series,
    INT4        dest,
    MPI_Comm    mpiComm
    );

void
LALMPIRecvINT4FrequencySeries(
    LALStatus  *status,
    INT4FrequencySeries *series,
    INT4        source,
    MPI_Comm    mpiComm
    );

void
LALMPISendINT8FrequencySeries(
    LALStatus  *status,
    INT8FrequencySeries *series,
    INT4        dest,
    MPI_Comm    mpiComm
    );

void
LALMPIRecvINT8FrequencySeries(
    LALStatus  *status,
    INT8FrequencySeries *series,
    INT4        source,
    MPI_Comm    mpiComm
    );

void
LALMPISendUINT2FrequencySeries(
    LALStatus   *status,
    UINT2FrequencySeries *series,
    INT4         dest,
    MPI_Comm     mpiComm
    );

void
LALMPIRecvUINT2FrequencySeries(
    LALStatus   *status,
    UINT2FrequencySeries *series,
    INT4         source,
    MPI_Comm     mpiComm
    );

void
LALMPISendUINT4FrequencySeries(
    LALStatus   *status,
    UINT4FrequencySeries *series,
    INT4         dest,
    MPI_Comm     mpiComm
    );

void
LALMPIRecvUINT4FrequencySeries(
    LALStatus   *status,
    UINT4FrequencySeries *series,
    INT4         source,
    MPI_Comm     mpiComm
    );

void
LALMPISendUINT8FrequencySeries(
    LALStatus   *status,
    UINT8FrequencySeries *series,
    INT4         dest,
    MPI_Comm     mpiComm
    );

void
LALMPIRecvUINT8FrequencySeries(
    LALStatus   *status,
    UINT8FrequencySeries *series,
    INT4         source,
    MPI_Comm     mpiComm
    );

void
LALMPISendREAL4FrequencySeries(
    LALStatus   *status,
    REAL4FrequencySeries *series,
    INT4         dest,
    MPI_Comm     mpiComm
    );

void
LALMPIRecvREAL4FrequencySeries(
    LALStatus   *status,
    REAL4FrequencySeries *series,
    INT4         source,
    MPI_Comm     mpiComm
    );

void
LALMPISendREAL8FrequencySeries(
    LALStatus   *status,
    REAL8FrequencySeries *series,
    INT4         dest,
    MPI_Comm     mpiComm
    );

void
LALMPIRecvREAL8FrequencySeries(
    LALStatus   *status,
    REAL8FrequencySeries *series,
    INT4         source,
    MPI_Comm     mpiComm
    );

void
LALMPISendCOMPLEX8FrequencySeries(
    LALStatus      *status,
    COMPLEX8FrequencySeries *series,
    INT4            dest,
    MPI_Comm        mpiComm
    );

void
LALMPIRecvCOMPLEX8FrequencySeries(
    LALStatus      *status,
    COMPLEX8FrequencySeries *series,
    INT4            source,
    MPI_Comm        mpiComm
    );

void
LALMPISendCOMPLEX16FrequencySeries(
    LALStatus       *status,
    COMPLEX16FrequencySeries *series,
    INT4             dest,
    MPI_Comm         mpiComm
    );

void
LALMPIRecvCOMPLEX16FrequencySeries(
    LALStatus       *status,
    COMPLEX16FrequencySeries *series,
    INT4             source,
    MPI_Comm         mpiComm
    );



void
LALInitializeExchange(
    LALStatus       *status,
    ExchParams     **exchParamsOut,
    ExchParams      *exchParamsInp,
    InitExchParams  *params
    );

void
LALFinalizeExchange(
    LALStatus   *status,
    ExchParams **exchParams
    );



void
LALExchangeCHAR(
    LALStatus  *status,
    CHAR       *object,
    ExchParams *exchParms
    );

void
LALExchangeINT2(
    LALStatus  *status,
    INT2       *object,
    ExchParams *exchParms
    );

void
LALExchangeINT4(
    LALStatus  *status,
    INT4       *object,
    ExchParams *exchParms
    );

void
LALExchangeINT8(
    LALStatus  *status,
    INT8       *object,
    ExchParams *exchParms
    );

void
LALExchangeUINT2(
    LALStatus  *status,
    UINT2      *object,
    ExchParams *exchParms
    );

void
LALExchangeUINT4(
    LALStatus  *status,
    UINT4      *object,
    ExchParams *exchParms
    );

void
LALExchangeUINT8(
    LALStatus  *status,
    UINT8      *object,
    ExchParams *exchParms
    );

void
LALExchangeREAL4(
    LALStatus  *status,
    REAL4      *object,
    ExchParams *exchParms
    );

void
LALExchangeREAL8(
    LALStatus  *status,
    REAL8      *object,
    ExchParams *exchParms
    );

void
LALExchangeCOMPLEX8(
    LALStatus  *status,
    COMPLEX8   *object,
    ExchParams *exchParms
    );

void
LALExchangeCOMPLEX16(
    LALStatus  *status,
    COMPLEX16  *object,
    ExchParams *exchParms
    );


void
LALExchangeCHARVector(
    LALStatus  *status,
    CHARVector       *object,
    ExchParams *exchParms
    );

void
LALExchangeINT2Vector(
    LALStatus  *status,
    INT2Vector       *object,
    ExchParams *exchParms
    );

void
LALExchangeINT4Vector(
    LALStatus  *status,
    INT4Vector       *object,
    ExchParams *exchParms
    );

void
LALExchangeINT8Vector(
    LALStatus  *status,
    INT8Vector       *object,
    ExchParams *exchParms
    );

void
LALExchangeUINT2Vector(
    LALStatus  *status,
    UINT2Vector      *object,
    ExchParams *exchParms
    );

void
LALExchangeUINT4Vector(
    LALStatus  *status,
    UINT4Vector      *object,
    ExchParams *exchParms
    );

void
LALExchangeUINT8Vector(
    LALStatus  *status,
    UINT8Vector      *object,
    ExchParams *exchParms
    );

void
LALExchangeREAL4Vector(
    LALStatus  *status,
    REAL4Vector      *object,
    ExchParams *exchParms
    );

void
LALExchangeREAL8Vector(
    LALStatus  *status,
    REAL8Vector      *object,
    ExchParams *exchParms
    );

void
LALExchangeCOMPLEX8Vector(
    LALStatus  *status,
    COMPLEX8Vector   *object,
    ExchParams *exchParms
    );

void
LALExchangeCOMPLEX16Vector(
    LALStatus  *status,
    COMPLEX16Vector  *object,
    ExchParams *exchParms
    );



void
LALExchangeINT2TimeSeries(
    LALStatus  *status,
    INT2TimeSeries       *object,
    ExchParams *exchParms
    );

void
LALExchangeINT4TimeSeries(
    LALStatus  *status,
    INT4TimeSeries       *object,
    ExchParams *exchParms
    );

void
LALExchangeINT8TimeSeries(
    LALStatus  *status,
    INT8TimeSeries       *object,
    ExchParams *exchParms
    );

void
LALExchangeUINT2TimeSeries(
    LALStatus  *status,
    UINT2TimeSeries      *object,
    ExchParams *exchParms
    );

void
LALExchangeUINT4TimeSeries(
    LALStatus  *status,
    UINT4TimeSeries      *object,
    ExchParams *exchParms
    );

void
LALExchangeUINT8TimeSeries(
    LALStatus  *status,
    UINT8TimeSeries      *object,
    ExchParams *exchParms
    );

void
LALExchangeREAL4TimeSeries(
    LALStatus  *status,
    REAL4TimeSeries      *object,
    ExchParams *exchParms
    );

void
LALExchangeREAL8TimeSeries(
    LALStatus  *status,
    REAL8TimeSeries      *object,
    ExchParams *exchParms
    );

void
LALExchangeCOMPLEX8TimeSeries(
    LALStatus  *status,
    COMPLEX8TimeSeries   *object,
    ExchParams *exchParms
    );

void
LALExchangeCOMPLEX16TimeSeries(
    LALStatus  *status,
    COMPLEX16TimeSeries  *object,
    ExchParams *exchParms
    );



void
LALExchangeINT2FrequencySeries(
    LALStatus  *status,
    INT2FrequencySeries       *object,
    ExchParams *exchParms
    );

void
LALExchangeINT4FrequencySeries(
    LALStatus  *status,
    INT4FrequencySeries       *object,
    ExchParams *exchParms
    );

void
LALExchangeINT8FrequencySeries(
    LALStatus  *status,
    INT8FrequencySeries       *object,
    ExchParams *exchParms
    );

void
LALExchangeUINT2FrequencySeries(
    LALStatus  *status,
    UINT2FrequencySeries      *object,
    ExchParams *exchParms
    );

void
LALExchangeUINT4FrequencySeries(
    LALStatus  *status,
    UINT4FrequencySeries      *object,
    ExchParams *exchParms
    );

void
LALExchangeUINT8FrequencySeries(
    LALStatus  *status,
    UINT8FrequencySeries      *object,
    ExchParams *exchParms
    );

void
LALExchangeREAL4FrequencySeries(
    LALStatus  *status,
    REAL4FrequencySeries      *object,
    ExchParams *exchParms
    );

void
LALExchangeREAL8FrequencySeries(
    LALStatus  *status,
    REAL8FrequencySeries      *object,
    ExchParams *exchParms
    );

void
LALExchangeCOMPLEX8FrequencySeries(
    LALStatus  *status,
    COMPLEX8FrequencySeries   *object,
    ExchParams *exchParms
    );

void
LALExchangeCOMPLEX16FrequencySeries(
    LALStatus  *status,
    COMPLEX16FrequencySeries  *object,
    ExchParams *exchParms
    );


#ifdef  __cplusplus
}
#endif

#endif /* _COMM_H */
