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

define(`VTYPE',`format(`%sVector',TYPE)')
define(`VSENDF',`format(`LALMPISend%s',VTYPE)')
define(`VRECVF',`format(`LALMPIRecv%s',VTYPE)')
define(`DTAG',`format(`MPI%sVectorData',TYPECODE)')
define(`VTAG',`format(`MPI%sVector',TYPECODE)')

/* <lalVerbatim file="SendRecvCP"> */
void VSENDF ( LALStatus *status, VTYPE *vector, INT4 dest, MPI_Comm comm )
{ /* </lalVerbatim> */
  INT4 code;

  INITSTATUS( status, "VSENDF", SENDRECVC );

  ASSERT( vector,         status, COMMH_ENULL, COMMH_MSGENULL );
  ASSERT( vector->data,   status, COMMH_ENULL, COMMH_MSGENULL );
  ASSERT( vector->length, status, COMMH_ESIZE, COMMH_MSGESIZE );

  /* send the vector structure */
  code = MPI_Send( vector, sizeof( *vector ), MPI_BYTE, dest, VTAG, comm );
  if ( code != MPI_SUCCESS )
  {
    ABORT( status, COMMH_EMPIE, COMMH_MSGEMPIE );
  }

  /* send the vector data */
  code = MPI_Send( vector->data, vector->length * sizeof( *vector->data ),
      MPI_BYTE, dest, DTAG, comm );
  if ( code != MPI_SUCCESS )
  {
    ABORT( status, COMMH_EMPIE, COMMH_MSGEMPIE );
  }

  RETURN( status );
}

/* <lalVerbatim file="SendRecvCP"> */
void VRECVF ( LALStatus *status, VTYPE *vector, INT4 source, MPI_Comm comm )
{ /* </lalVerbatim> */
  MPI_Status stat;
  VTYPE vtmp;
  INT4 code;

  INITSTATUS( status, "VRECVF", SENDRECVC );

  ASSERT( vector,         status, COMMH_ENULL, COMMH_MSGENULL );
  ASSERT( vector->data,   status, COMMH_ENULL, COMMH_MSGENULL );
  ASSERT( vector->length, status, COMMH_ESIZE, COMMH_MSGESIZE );

  /* receive temporary vector structure to check */
  code = MPI_Recv( &vtmp, sizeof( vtmp ), MPI_BYTE, source, VTAG, comm, &stat );
  if ( code != MPI_SUCCESS || stat.MPI_ERROR != MPI_SUCCESS )
  {
    ABORT( status, COMMH_EMPIE, COMMH_MSGEMPIE );
  }

  /* make sure that the vector lengths agree */
  if ( vector->length != vtmp.length )
  {
    ABORT( status, COMMH_ESZMM, COMMH_MSGESZMM );
  }

  /* receive the vector data */
  code = MPI_Recv( vector->data, vector->length * sizeof( *vector->data ),
      MPI_BYTE, source, DTAG, comm, &stat );

  if ( code != MPI_SUCCESS || stat.MPI_ERROR != MPI_SUCCESS )
  {
    ABORT( status, COMMH_EMPIE, COMMH_MSGEMPIE );
  }

  RETURN( status );
}
