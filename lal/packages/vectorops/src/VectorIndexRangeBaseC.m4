/** -*- C -*- **/
/*
dnl $Id$
ifelse(TYPECODE,`CHAR',`define(`TYPE',`CHAR')')
ifelse(TYPECODE,`I2',`define(`TYPE',`INT2')')
ifelse(TYPECODE,`I4',`define(`TYPE',`INT4')')
ifelse(TYPECODE,`I8',`define(`TYPE',`INT8')')
ifelse(TYPECODE,`U2',`define(`TYPE',`UINT2')')
ifelse(TYPECODE,`U4',`define(`TYPE',`UINT4')')
ifelse(TYPECODE,`U8',`define(`TYPE',`UINT8')')
ifelse(TYPECODE,`S',`define(`TYPE',`REAL4')')
ifelse(TYPECODE,`D',`define(`TYPE',`REAL8')')
ifelse(TYPECODE,`C',`define(`TYPE',`COMPLEX8')')
ifelse(TYPECODE,`Z',`define(`TYPE',`COMPLEX16')')
define(`VTYPE',`format(`%sVector',TYPE)')
define(`F1',`format(`LAL%sVectorIndexRange',TYPECODE)')
*/

/******************************* <lalLaTeX file="VectorIndexRangeC">
\begin{verbatim}void F1 ( LALStatus *status, VTYPE **result,
       VTYPE *A, const UINT4Vector *indexVector )\end{verbatim}
*************************************************** </lalLaTeX> */

/* Take a vector, and remove a range entries by given indices, then stuff
   that truncated vector in to the result vector */
void F1 ( 
         LALStatus  *status,
         VTYPE     **p_result,
         VTYPE      *p_v,
         const UINT4Vector *p_indexVector
         )
{
  /*  Variable Declarations  */
  UINT4 iter;
  UINT4 rslt_length;

  INITSTATUS( status, "F1" , VECTORINDEXRANGEC );
  ATTATCHSTATUSPTR( status );

  /*
   * error-checking assertions
   */
  ASSERT(p_v, status, VECTORINDEXRANGEH_ENULL, VECTORINDEXRANGEH_MSGENULL);
  ASSERT(p_v->length > 1, status, VECTORINDEXRANGEH_ELNTH, 
         VECTORINDEXRANGEH_MSGELNTH);

  ASSERT(p_indexVector, status, VECTORINDEXRANGEH_ENULL, 
         VECTORINDEXRANGEH_MSGENULL );
  ASSERT(p_indexVector->length == 2, status, VECTORINDEXRANGEH_ELNTH,
         VECTORINDEXRANGEH_MSGELNTH );

  rslt_length = p_indexVector->data[1] - p_indexVector->data[0] + 1;

  if (lalDebugLevel & 8)
    {
      printf("rslt_length = %d\n", rslt_length);
      printf("BEFORE ALLOC:  p_result = %#x\n", p_result);
      printf("              *p_result = %#x\n", *p_result);
    }

  /* Handle memory */
  if ((*p_result != NULL) && ((*p_result)->length != rslt_length))
    {
      TRY( LAL`'TYPECODE`'DestroyVector(status->statusPtr, p_result),
           status );
    }
  
  /* allocate */
  TRY( LAL`'TYPECODE`'CreateVector(status->statusPtr, p_result, rslt_length),
       status );

  if (lalDebugLevel & 8)
    {
      printf("AFTER ALLOC: p_result = %#x\n", p_result);
      printf("            *p_result = %#x\n", *p_result);
      printf("            (*p_result)->length = %d\n", (*p_result)->length);
    }

  /* fill in the result vector */
  for (iter = p_indexVector->data[0]; iter <= p_indexVector->data[1]; ++iter)
    {
      (*p_result)->data[iter - p_indexVector->data[0]] = p_v->data[iter];
    }
  
  DETATCHSTATUSPTR( status );
  RETURN (status);
}

