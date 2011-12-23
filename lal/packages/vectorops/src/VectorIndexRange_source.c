/** -*- C -*- **/

#define CONCAT2x(a,b) a##b
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a##b##c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define STRING(a) #a

#define VTYPE CONCAT2(TYPE,Vector)
#define VPAIRTYPE CONCAT2(TYPE,VectorPair)

#define F1 CONCAT3(LAL,TYPECODE,VectorIndexRange)
#define F2 CONCAT3(LAL,TYPECODE,VectorIndexHole)
#define CFUNC CONCAT3(LAL,TYPECODE,CreateVector)
#define DFUNC CONCAT3(LAL,TYPECODE,DestroyVector)

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

  INITSTATUS(status);
  ATTATCHSTATUSPTR( status );

  /*
   * error-checking assertions
   */
  ASSERT(p_v, status, VECTORINDEXRANGEH_ENULL, VECTORINDEXRANGEH_MSGENULL);
  ASSERT(p_v->length >= 1, status, VECTORINDEXRANGEH_ELNTH,
         VECTORINDEXRANGEH_MSGELNTH);

  ASSERT(p_indexVector, status, VECTORINDEXRANGEH_ENULL,
         VECTORINDEXRANGEH_MSGENULL);
  ASSERT(p_indexVector->length == 2, status, VECTORINDEXRANGEH_ELNTH,
         VECTORINDEXRANGEH_MSGELNTH);

  rslt_length = p_indexVector->data[1] - p_indexVector->data[0] + 1;

  /* Handle memory */
  if ((*p_result != NULL) && ((*p_result)->length != rslt_length))
    {
      TRY(DFUNC(status->statusPtr, p_result),
          status);
    }

  /* allocate */
  TRY(CFUNC(status->statusPtr, p_result, rslt_length),
      status);

  /* fill in the result vector */
  for (iter = p_indexVector->data[0]; iter <= p_indexVector->data[1]; ++iter)
    {
      (*p_result)->data[iter - p_indexVector->data[0]] = p_v->data[iter];
    }

  DETATCHSTATUSPTR( status );
  RETURN (status);
}



void F2 ( LALStatus  *status,
          VPAIRTYPE  *p_result_pair,
          VTYPE      *p_v,
          const UINT4Vector *p_indexVector )
{
  /*  Variable Declarations  */
  UINT4Vector *head_index_range = NULL;
  UINT4Vector *tail_index_range = NULL;

  INITSTATUS(status);
  ATTATCHSTATUSPTR( status );



  /*
    assert indexVector range
    p_result_pair->head = VectorIndexRange(..., [0,indexVector->data[0]-1]);
    p_result_pair->tail = VectorIndexRange(..., [indexVector->data[1]+1,
                          p_v->length-1])
  */

  TRY(LALU4CreateVector(status->statusPtr, &head_index_range, 2),
      status);

  LALU4CreateVector(status->statusPtr, &tail_index_range, 2);

  BEGINFAIL(status)
    {
      TRY(LALU4DestroyVector(status->statusPtr, &head_index_range), status);
    }
  ENDFAIL(status);

  head_index_range->data[0] = 0;
  head_index_range->data[1] = p_indexVector->data[0] - 1;
  tail_index_range->data[0] = p_indexVector->data[1] + 1;
  tail_index_range->data[1] = p_v->length - 1;

  TRY(F1(status->statusPtr, p_result_pair->head,
                                      p_v, head_index_range), status);

  F1(status->statusPtr, p_result_pair->tail,
                                  p_v, tail_index_range);
  BEGINFAIL(status)
    TRY(DFUNC(status->statusPtr, p_result_pair->head),
        status);
  ENDFAIL(status);

  TRY(LALU4DestroyVector(status->statusPtr, &head_index_range), status);
  TRY(LALU4DestroyVector(status->statusPtr, &tail_index_range), status);

  DETATCHSTATUSPTR( status );
  RETURN( status );
}

#undef VTYPE
#undef VPAIRTYPE
#undef F1
#undef F2
#undef CFUNC
#undef DFUNC
