dnl $Id$

define(`XFUNC',`format(`LALExchange%s',TYPE)')
define(`SFUNC',`format(`LALMPISend%s',TYPE)')
define(`RFUNC',`format(`LALMPIRecv%s',TYPE)')

/* <lalVerbatim file="ExchangeCP"> */
void XFUNC ( LALStatus *status, TYPE *object, ExchParams *exchParams )
{ /* </lalVerbatim> */
  INITSTATUS( status, "XFUNC", EXCHANGEC );
  ATTATCHSTATUSPTR( status );

  ASSERT( exchParams, status, COMMH_ENULL, COMMH_MSGENULL );
  ASSERT( object, status, COMMH_ENULL, COMMH_MSGENULL );

  if ( exchParams->send ) /* I am sending */
  {
    SFUNC ( status->statusPtr, object, exchParams->partnerProcNum, exchParams->mpiComm );
    CHECKSTATUSPTR( status );
  }
  else /* I am receiving */
  {
    RFUNC ( status->statusPtr, object, exchParams->partnerProcNum, exchParams->mpiComm );
    CHECKSTATUSPTR( status );
  }

  DETATCHSTATUSPTR( status );
  RETURN( status );
}
