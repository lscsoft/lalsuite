dnl $Id$

ifelse(TYPECODE,`Z',`define(`TYPE',`COMPLEX16')')
ifelse(TYPECODE,`C',`define(`TYPE',`COMPLEX8')')
ifelse(TYPECODE,`D',`define(`TYPE',`REAL8')')
ifelse(TYPECODE,`S',`define(`TYPE',`REAL4')')
ifelse(TYPECODE,`I2',`define(`TYPE',`INT2')')
ifelse(TYPECODE,`I4',`define(`TYPE',`INT4')')
ifelse(TYPECODE,`I8',`define(`TYPE',`INT8')')
ifelse(TYPECODE,`U2',`define(`TYPE',`UINT2')')
ifelse(TYPECODE,`U4',`define(`TYPE',`UINT4')')
ifelse(TYPECODE,`U8',`define(`TYPE',`UINT8')')
ifelse(TYPECODE,`CHAR',`define(`TYPE',`CHAR')')

ifelse(SERIESCODE,`T',`define(`SERIES',`TimeSeries')')
ifelse(SERIESCODE,`F',`define(`SERIES',`FrequencySeries')')

define(`VTYPE',`format(`%sVector',TYPE)')
define(`STYPE',`format(`%s%s',TYPE,SERIES)')
define(`VSENDF',`format(`LALMPISend%s',VTYPE)')
define(`VRECVF',`format(`LALMPIRecv%s',VTYPE)')
define(`SSENDF',`format(`LALMPISend%s',STYPE)')
define(`SRECVF',`format(`LALMPIRecv%s',STYPE)')
define(`STAG',`format(`MPI%sTimeSeries',TYPECODE)')

/* <lalVerbatim file="SendRecvCP"> */
void SSENDF ( LALStatus *status, STYPE *series, INT4 dest, MPI_Comm comm )
{ /* </lalVerbatim> */
  INT4 code;

  INITSTATUS( status, "SSENDF", SENDRECVC );
  ATTATCHSTATUSPTR( status );

  ASSERT( series,       status, COMMH_ENULL, COMMH_MSGENULL );
  ASSERT( series->data, status, COMMH_ENULL, COMMH_MSGENULL );

  /* send the series structure */
  code = MPI_Send( series, sizeof( *series ), MPI_BYTE, dest, STAG, comm );
  if ( code != MPI_SUCCESS )
  {
    ABORT( status, COMMH_EMPIE, COMMH_MSGEMPIE );
  }

  /* send data vector */
  VSENDF ( status->statusPtr, series->data, dest, comm );
  CHECKSTATUSPTR( status );

  DETATCHSTATUSPTR( status );
  RETURN( status );
}

/* <lalVerbatim file="SendRecvCP"> */
void SRECVF ( LALStatus *status, STYPE *series, INT4 source, MPI_Comm comm )
{ /* </lalVerbatim> */
  MPI_Status stat;
  VTYPE *vtmp;
  INT4 code;

  INITSTATUS( status, "SRECVF", SENDRECVC );
  ATTATCHSTATUSPTR( status );

  ASSERT( series,       status, COMMH_ENULL, COMMH_MSGENULL );
  ASSERT( series->data, status, COMMH_ENULL, COMMH_MSGENULL );

  /* receive time series structure */
  vtmp = series->data;
  code = MPI_Recv( series, sizeof( *series ), MPI_BYTE, source, STAG,
      comm, &stat );
  if ( code != MPI_SUCCESS || stat.MPI_ERROR != MPI_SUCCESS )
  {
    ABORT( status, COMMH_EMPIE, COMMH_MSGEMPIE );
  }
  series->data = vtmp;

  /* receive data vector */
  VRECVF ( status->statusPtr, series->data, source, comm );
  CHECKSTATUSPTR( status );

  DETATCHSTATUSPTR( status );
  RETURN( status );
}
