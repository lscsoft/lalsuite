dnl $Id$
ifelse(TYPECODE,`Z',`define(`TYPE',`COMPLEX16')')dnl
ifelse(TYPECODE,`C',`define(`TYPE',`COMPLEX8')')dnl
ifelse(TYPECODE,`D',`define(`TYPE',`REAL8')')dnl
ifelse(TYPECODE,`S',`define(`TYPE',`REAL4')')dnl
ifelse(TYPECODE,`I2',`define(`TYPE',`INT2')')dnl
ifelse(TYPECODE,`I4',`define(`TYPE',`INT4')')dnl
ifelse(TYPECODE,`I8',`define(`TYPE',`INT8')')dnl
ifelse(TYPECODE,`U2',`define(`TYPE',`UINT2')')dnl
ifelse(TYPECODE,`U4',`define(`TYPE',`UINT4')')dnl
ifelse(TYPECODE,`U8',`define(`TYPE',`UINT8')')dnl
define(`GTYPE',`format(`%sGrid',TYPE)')dnl
define(`CFUNC',`format(`LAL%sCreateGrid',TYPECODE)')dnl
define(`DFUNC',`format(`LAL%sDestroyGrid',TYPECODE)')dnl
define(`CREATEFUNC',`format(`LAL%sCreateArray',TYPECODE)')dnl
define(`DESTROYFUNC',`format(`LAL%sDestroyArray',TYPECODE)')dnl

void
CFUNC ( LALStatus *stat, GTYPE **grid, UINT4Vector *dimLength, UINT4 dimension )
{
  INITSTATUS( stat, "CFUNC", GRIDC );
  ATTATCHSTATUSPTR( stat );

  /* Check for valid input arguments. */
  ASSERT( dimLength, stat, GRIDH_ENUL, GRIDH_MSGENUL );
  ASSERT( dimLength->data, stat, GRIDH_ENUL, GRIDH_MSGENUL );
  ASSERT( grid, stat, GRIDH_ENUL, GRIDH_MSGENUL );
  ASSERT( !*grid, stat, GRIDH_EOUT, GRIDH_MSGEOUT );
  ASSERT( dimension <= dimLength->length, stat, GRIDH_ENUL,
	  GRIDH_MSGENUL );

  /* Allocate the grid. */
  if ( !( *grid = ( GTYPE *)LALMalloc( sizeof( GTYPE ) ) ) ) {
    ABORT( stat, GRIDH_EMEM, GRIDH_MSGEMEM );
  }
  memset( *grid, 0, sizeof( GTYPE ) );

  /* Allocate the dimension units array. */
  if ( !( (*grid)->dimUnits = (LALUnit *)
	  LALMalloc( dimension*sizeof(LALUnit) ) ) ) {
    LALFree( *grid );
    *grid = NULL;
    ABORT( stat, GRIDH_EMEM, GRIDH_MSGEMEM );
  }
  memset( (*grid)->dimUnits, 0, dimension*sizeof(LALUnit) );

  /* Allocate the offset and interval vectors. */
  LALDCreateVector( stat->statusPtr, &((*grid)->offset), dimension );
  BEGINFAIL( stat ) {
    LALFree( (*grid)->dimUnits );
    LALFree( *grid );
    *grid = NULL;
  } ENDFAIL( stat );
  LALDCreateVector( stat->statusPtr, &((*grid)->interval), dimension );
  BEGINFAIL( stat ) {
    TRY( LALDDestroyVector( stat->statusPtr, &((*grid)->offset) ), stat );
    LALFree( (*grid)->dimUnits );
    LALFree( *grid );
    *grid = NULL;
  } ENDFAIL( stat );

  /* Allocate the data array. */
  CREATEFUNC ( stat->statusPtr, &((*grid)->data), dimLength );
  BEGINFAIL( stat ) {
    TRY( LALDDestroyVector( stat->statusPtr, &((*grid)->interval) ), stat );
    TRY( LALDDestroyVector( stat->statusPtr, &((*grid)->offset) ), stat );
    LALFree( (*grid)->dimUnits );
    LALFree( *grid );
    *grid = NULL;
  } ENDFAIL( stat );

  /* Done. */
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}


void
DFUNC ( LALStatus *stat, GTYPE **grid )
{
  INITSTATUS( stat, "DFUNC", GRIDC );
  ATTATCHSTATUSPTR( stat );

  /* Check for valid input argument. */
  ASSERT( grid, stat, GRIDH_ENUL, GRIDH_MSGENUL );
  ASSERT( *grid, stat, GRIDH_ENUL, GRIDH_MSGENUL );

  /* Destroy the internal arrays and vectors. */
  TRY( DESTROYFUNC ( stat->statusPtr, &((*grid)->data) ), stat );
  TRY( LALDDestroyVector( stat->statusPtr, &((*grid)->interval) ), stat );
  TRY( LALDDestroyVector( stat->statusPtr, &((*grid)->offset) ), stat );
  LALFree( (*grid)->dimUnits );

  /* Destroy the structure itself, and exit. */
  LALFree( *grid );
  *grid = NULL;
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}
