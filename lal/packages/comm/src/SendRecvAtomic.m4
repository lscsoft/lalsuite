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

define(`ASENDF',`format(`LALMPISend%s',TYPE)')
define(`ARECVF',`format(`LALMPIRecv%s',TYPE)')
define(`TAG',`format(`MPI%s',TYPECODE)')

/* <lalVerbatim file="SendRecvCP"> */
void ASENDF ( LALStatus *status, TYPE *element, INT4 dest, MPI_Comm comm )
{ /* </lalVerbatim> */
  INT4 code;

  INITSTATUS( status, "ASENDF", SENDRECVC );

  ASSERT( element, status, COMMH_ENULL, COMMH_MSGENULL );

  /* send the element */
  code = MPI_Send( element, sizeof( *element ), MPI_BYTE, dest, TAG, comm );
  if ( code != MPI_SUCCESS )
  {
    ABORT( status, COMMH_EMPIE, COMMH_MSGEMPIE );
  }

  RETURN( status );
}

/* <lalVerbatim file="SendRecvCP"> */
void ARECVF ( LALStatus *status, TYPE *element, INT4 source, MPI_Comm comm )
{ /* </lalVerbatim> */
  MPI_Status stat;
  INT4 code;

  INITSTATUS( status, "ARECVF", SENDRECVC );

  ASSERT( element, status, COMMH_ENULL, COMMH_MSGENULL );

  code = MPI_Recv( element, sizeof( *element ), MPI_BYTE, source, TAG, comm, &stat );
  if ( code != MPI_SUCCESS || stat.MPI_ERROR != MPI_SUCCESS )
  {
    ABORT( status, COMMH_EMPIE, COMMH_MSGEMPIE );
  }

  RETURN( status );
}
