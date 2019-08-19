#define CONCAT2x(a,b) a##b
#define CONCAT2(a,b) CONCAT2x(a,b)
#define STRING(a) #a

#ifdef COMPLEX_DATA
#   define DBLDATATYPE COMPLEX16
#   ifdef SINGLE_PRECISION
#       define DATATYPE COMPLEX8
#   else
#       define DATATYPE COMPLEX16
#   endif
#else
#   define DBLDATATYPE REAL8
#   ifdef SINGLE_PRECISION
#       define DATATYPE REAL4
#   else
#       define DATATYPE REAL8
#   endif
#endif


#define VECTORTYPE CONCAT2(DATATYPE,Vector)
#define FILTERTYPE CONCAT2(DBLDATATYPE,IIRFilter)

#define FUNC CONCAT2(XLALIIRFilter,VECTORTYPE)

int FUNC(VECTORTYPE *vector, FILTERTYPE *filter)
{
  INT4 j;            /* Index for filter coeficients. */
  INT4 length;       /* Length of vector. */
  DATATYPE *data;    /* Vector data. */
  DBLDATATYPE w, datum;    /* Current auxiliary and output values. */
  INT4 directOrder;  /* Number of direct filter coefficients. */
  INT4 recursOrder;  /* Number of recursive filter coefficients. */
  INT4 numHist;      /* The number of history data. */
  REAL8 *directCoef; /* Direct filter coefficients. */
  REAL8 *recursCoef; /* Recursive filter coefficients. */
  DBLDATATYPE *temp=NULL;  /* Temporary storage for the filter history. */

  /* Make sure all the structures have been initialized. */
  if ( ! vector || ! filter )
    XLAL_ERROR( XLAL_EFAULT );
  if ( ! vector->data )
    XLAL_ERROR( XLAL_EINVAL );
  if ( ! filter->directCoef || ! filter->recursCoef || ! filter->history
      || !  filter->directCoef->data || ! filter->recursCoef->data
      || !  filter->history->data )
    XLAL_ERROR( XLAL_EINVAL );

  length=vector->length;
  data=vector->data;
  directOrder=filter->directCoef->length;
  recursOrder=filter->recursCoef->length;
  directCoef=filter->directCoef->data;
  recursCoef=filter->recursCoef->data;
  numHist=filter->history->length+1;
  temp = LALMalloc( numHist*sizeof(*temp) );
  if ( ! temp )
    XLAL_ERROR( XLAL_ENOMEM );
  memcpy(temp,filter->history->data,(numHist-1)*sizeof(*temp));

  /* Run through the vector. */
  while(length--){

    /* Compute the auxiliary variable. */
    for(j=numHist-1;j>=recursOrder;j--)
      temp[j]=temp[j-1];
    w=*data;
    for(;j;j--)
      w+=recursCoef[j]*(temp[j]=temp[j-1]);

    /* Compute filter output. */
    datum=*directCoef*(*temp=w);
    for(j=1;j<directOrder;j++)
      datum+=directCoef[j]*temp[j];
    *(data++)=datum;
  }

  /* Update the history. */
  memcpy(filter->history->data,temp,(numHist-1)*sizeof(*temp));
  LALFree(temp);

  /* Normal exit */
  return 0;
}

#undef FUNC
#undef VECTORTYPE
#undef FILTERTYPE
#undef DBLDATATYPE
#undef DATATYPE
#undef CONCAT2x
#undef CONCAT2
#undef STRING
