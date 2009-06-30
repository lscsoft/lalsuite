dnl $Id$
ifelse(TYPECODE,`D',`define(`TYPE',`REAL8')')
ifelse(TYPECODE,`S',`define(`TYPE',`REAL4')')
ifelse(TYPECODE,`I2',`define(`TYPE',`INT2')')
ifelse(TYPECODE,`I4',`define(`TYPE',`INT4')')
ifelse(TYPECODE,`I8',`define(`TYPE',`INT8')')
ifelse(TYPECODE2,`D',`define(`TYPE2',`REAL8')')
ifelse(TYPECODE2,`S',`define(`TYPE2',`REAL4')')
ifelse(TYPECODE2,`I2',`define(`TYPE2',`INT2')')
ifelse(TYPECODE2,`I4',`define(`TYPE2',`INT4')')
ifelse(TYPECODE2,`I8',`define(`TYPE2',`INT8')')
define(`VTYPE',`format(`%sVector',TYPE)')
define(`ATYPE',`format(`%sArray',TYPE)')
define(`VTYPE2',`format(`%sVector',TYPE2)')
define(`ATYPE2',`format(`%sArray',TYPE2)')
define(`F1',`format(`LAL%sDotPower%sVector',TYPECODE,TYPECODE2)')
define(`F2',`format(`LAL%sVectorDotPower%s',TYPECODE,TYPECODE2)')
define(`F3',`format(`LAL%sVectorDotPower%sVector',TYPECODE,TYPECODE2)')
define(`F4',`format(`LAL%sDotPower%sArray',TYPECODE,TYPECODE2)')
define(`F5',`format(`LAL%sArrayDotPower%s',TYPECODE,TYPECODE2)')
define(`F6',`format(`LAL%sArrayDotPower%sArray',TYPECODE,TYPECODE2)')

/******************************* <lalLaTeX file="MatrixPowerC">
\begin{verbatim}void F1 ( LALStatus *status, REAL8Vector **result,
		TYPE A, VTYPE2 *B )\end{verbatim}
 ************************************************** </lalLaTeX> */

void F1 (
        LALStatus		*status,
        REAL8Vector		**result,
        TYPE			A,
        VTYPE2		*B
)
{
        /*  Variable Declarations  */
        UINT4    iterator;
        UINT4   length;

        INITSTATUS( status, "F1" , MATLABMATRIXPOWC);
        ATTATCHSTATUSPTR( status );

        /*  Check input for existence.  */
        /*  Result should not come in Allocated  */
	ASSERT ( result, status, MATLABMATRIXH_ENULL, MATLABMATRIXH_MSGENULL);
        ASSERT ( !(*result), status, MATLABMATRIXH_ENNUL, MATLABMATRIXH_MSGENNUL);

        /*  data must be defined  */
        ASSERT ( B, status, MATLABMATRIXH_ENULL, MATLABMATRIXH_MSGENULL);

        length = ((VTYPE*)(B))->length;

        /*  length must be greater than one  */
        ASSERT ( length > 1, status, MATLABMATRIXH_ELNTH, MATLABMATRIXH_MSGELNTH);

        LALDCreateVector( status->statusPtr, result, length);

        for (iterator = 0; iterator < length; iterator++)
        {
                (*result)->data[iterator] = pow( B->data[iterator], A );
        }

        DETATCHSTATUSPTR( status );
        RETURN (status);
}

/******************************* <lalLaTeX file="MatrixPowerC">
\begin{verbatim}void F2 ( LALStatus *status, REAL8Vector **result,
		VTYPE *A, TYPE2 B )\end{verbatim}
 ************************************************** </lalLaTeX> */

void F2 (
        LALStatus               *status,
        REAL8Vector		**result,
        VTYPE		*A,
        TYPE2			B
)
{
        /*  Variable Declarations  */
        UINT4    iterator;
        UINT4   length;

        INITSTATUS( status, "F2" , MATLABMATRIXPOWC);
        ATTATCHSTATUSPTR( status );

        /*  Check input for existence.  */
        /*  Result should not come in Allocated  */
	ASSERT ( result, status, MATLABMATRIXH_ENULL, MATLABMATRIXH_MSGENULL);
        ASSERT ( !(*result), status, MATLABMATRIXH_ENNUL, MATLABMATRIXH_MSGENNUL);

        /*  data must be defined  */
        ASSERT ( A, status, MATLABMATRIXH_ENULL, MATLABMATRIXH_MSGENULL);

        length = ((VTYPE*)(A))->length;

        /*  length must be greater than one  */
        ASSERT ( length > 1, status, MATLABMATRIXH_ELNTH, MATLABMATRIXH_MSGELNTH);

        LALDCreateVector( status->statusPtr, result, length);

        for (iterator = 0; iterator < length; iterator++)
        {
                (*result)->data[iterator] = pow( A->data[iterator], B );
        }

        DETATCHSTATUSPTR( status );
        RETURN (status);
}

/******************************* <lalLaTeX file="MatrixPowerC">
\begin{verbatim}void F3 ( LALStatus *status, REAL8Vector **result,
		VTYPE *B, VTYPE2 *A )\end{verbatim}
 ************************************************** </lalLaTeX> */

void F3 (
        LALStatus		*status,
        REAL8Vector		**result,
        VTYPE		*B,
        VTYPE2		*A
)
{
        /*  Variable Declarations  */
        UINT4    iterator;
        UINT4   length;

        INITSTATUS( status, "F3" , MATLABMATRIXPOWC);
        ATTATCHSTATUSPTR( status );

        /*  Check input for existence.  */
        /*  Result should not come in Allocated  */
	ASSERT ( result, status, MATLABMATRIXH_ENULL, MATLABMATRIXH_MSGENULL);
        ASSERT ( !(*result), status, MATLABMATRIXH_ENNUL, MATLABMATRIXH_MSGENNUL);


        /*  data must be defined  */
        ASSERT ( A, status, MATLABMATRIXH_ENULL, MATLABMATRIXH_MSGENULL);

	ASSERT ( B, status, MATLABMATRIXH_ENULL, MATLABMATRIXH_MSGENULL);

	ASSERT ( A->length == B->length, status, MATLABMATRIXH_ELNTH, MATLABMATRIXH_MSGELNTH);

        length = ((VTYPE*)(A))->length;

        /*  length must be greater than one  */
        ASSERT ( length > 1, status, MATLABMATRIXH_ELNTH, MATLABMATRIXH_MSGELNTH);

        LALDCreateVector( status->statusPtr, result, length);

        for (iterator = 0; iterator < length; iterator++)
        {
                (*result)->data[iterator] = pow( A->data[iterator], B->data[iterator] );
        }

        DETATCHSTATUSPTR( status );
        RETURN (status);
}

/******************************* <lalLaTeX file="MatrixPowerC">
\begin{verbatim}void F4 ( LALStatus *status, REAL8Array **result,
		TYPE A, ATYPE2 *B )\end{verbatim}
 ************************************************** </lalLaTeX> */

void F4 (
        LALStatus               *status,
        REAL8Array		**result,
        TYPE			A,
        ATYPE2		*B
)
{
        /*  Variable Declarations  */
        UINT4Vector     *length;
        UINT4           ndims;
        UINT4		iterator, myindex;
	UINT4		row, column;

        INITSTATUS( status, "F4" , MATLABMATRIXPOWC);
        ATTATCHSTATUSPTR( status );

        /*  Check input for existence.  */
        /*  Result should not come in Allocated  */
	ASSERT ( result, status, MATLABMATRIXH_ENULL, MATLABMATRIXH_MSGENULL);
        ASSERT ( !(*result), status, MATLABMATRIXH_ENNUL, MATLABMATRIXH_MSGENNUL);

        /*  data must be defined  */
	ASSERT ( B, status, MATLABMATRIXH_ENULL, MATLABMATRIXH_MSGENULL);

	ASSERT ( B->dimLength, status, MATLABMATRIXH_ENULL, MATLABMATRIXH_MSGENULL);

	ndims = ((ATYPE2*)(B))->dimLength->length;
	length = NULL;

	LALU4CreateVector( status->statusPtr, &length, ndims);

	for ( iterator = 0; iterator < ndims; iterator++)
	{
		length->data[iterator] = ((ATYPE2*)(B))->dimLength->data[iterator];
	}

	/*  length must be greater than one  */
        ASSERT ( length->data[0] > 1, status, MATLABMATRIXH_ELNTH, MATLABMATRIXH_MSGELNTH);

        LALDCreateArray( status->statusPtr, result, length);

	if( ndims == 2 )
	{
	        for( row = 0; row < length->data[0]; row++)
		{
			for( column = 0; column < length->data[1]; column++)
			{
				myindex = (row * length->data[1]) + column;
				(*result)->data[myindex] = pow( A, B->data[myindex] );
			}
	        }
        }
        else
        {
		LALDDestroyArray( status->statusPtr, result);
		(*result) = NULL;
        }


	LALU4DestroyVector( status->statusPtr, &length );

        DETATCHSTATUSPTR( status );
        RETURN (status);
}

/******************************* <lalLaTeX file="MatrixPowerC">
\begin{verbatim}void F5 ( LALStatus *status, REAL8Array **result,
		ATYPE *A, TYPE2 B )\end{verbatim}
 ************************************************** </lalLaTeX> */

void F5 (
        LALStatus		*status,
        REAL8Array		**result,
        ATYPE		*A,
        TYPE2			B
)
{
        /*  Variable Declarations  */
        UINT4Vector     *length;
        UINT4           ndims;
        UINT4            iterator, myindex;
        UINT4            row, column;

        INITSTATUS( status, "F5" , MATLABMATRIXPOWC);
        ATTATCHSTATUSPTR( status );

        /*  Check input for existence.  */
        /*  Result should not come in Allocated  */
	ASSERT ( result, status, MATLABMATRIXH_ENULL, MATLABMATRIXH_MSGENULL);
        ASSERT ( !(*result), status, MATLABMATRIXH_ENNUL, MATLABMATRIXH_MSGENNUL);


        /*  data must be defined  */
        ASSERT ( A, status, MATLABMATRIXH_ENULL, MATLABMATRIXH_MSGENULL);

        ASSERT ( A->dimLength, status, MATLABMATRIXH_ENULL, MATLABMATRIXH_MSGENULL);

        ndims = ((ATYPE*)(A))->dimLength->length;
	length = NULL;

        LALU4CreateVector( status->statusPtr, &length, ndims);

        for ( iterator = 0; iterator < ndims; iterator++)
        {
                length->data[iterator] = ((ATYPE*)(A))->dimLength->data[iterator];
        }

        /*  length must be greater than one  */
        ASSERT ( length->data[0] > 1, status, MATLABMATRIXH_ELNTH, MATLABMATRIXH_MSGELNTH);

        LALDCreateArray( status->statusPtr, result, length);

	if( ndims == 2 )
	{
		for( row = 0; row < length->data[0]; row++)
		{
			for( column = 0; column < length->data[1]; column++)
			{
                                myindex = (row * length->data[1]) + column;
				(*result)->data[myindex] = pow( A->data[myindex], B );
			}
		}
        }
        else
        {
		LALDDestroyArray( status->statusPtr, result);
		(*result) = NULL;
        }


        LALU4DestroyVector( status->statusPtr, &length );

        DETATCHSTATUSPTR( status );
        RETURN (status);
}

/******************************* <lalLaTeX file="MatrixPowerC">
\begin{verbatim}void F6 ( LALStatus *status, REAL8Array **result,
		ATYPE *A, ATYPE2 *B )\end{verbatim}
 ************************************************** </lalLaTeX> */

void F6 (
        LALStatus		*status,
        REAL8Array		**result,
        ATYPE		*A,
        ATYPE2		*B
)
{
        /*  Variable Declarations  */
        UINT4Vector     *length;
        UINT4           ndims;
        UINT4           ndims2;
        UINT4            iterator, myindex;
        UINT4            row, column;

        INITSTATUS( status, "F6" , MATLABMATRIXPOWC);
        ATTATCHSTATUSPTR( status );

        /*  Check input for existence.  */
        /*  Result should not come in Allocated  */
	ASSERT ( result, status, MATLABMATRIXH_ENULL, MATLABMATRIXH_MSGENULL);
        ASSERT ( !(*result), status, MATLABMATRIXH_ENNUL, MATLABMATRIXH_MSGENNUL);


        /*  data must be defined  */
        ASSERT ( B, status, MATLABMATRIXH_ENULL, MATLABMATRIXH_MSGENULL);

        ASSERT ( B->dimLength, status, MATLABMATRIXH_ENULL, MATLABMATRIXH_MSGENULL);

        ndims = ((ATYPE*)(B))->dimLength->length;
	length = NULL;

        LALU4CreateVector( status->statusPtr, &length, ndims);

        for ( iterator = 0; iterator < ndims; iterator++)
        {
                length->data[iterator] = ((ATYPE*)(B))->dimLength->data[iterator];
        }

        /*  data must be defined  */
        ASSERT ( A, status, MATLABMATRIXH_ENULL, MATLABMATRIXH_MSGENULL);

        ASSERT ( A->dimLength, status, MATLABMATRIXH_ENULL, MATLABMATRIXH_MSGENULL);

        ndims2 = ((ATYPE2*)(A))->dimLength->length;

	ASSERT ( ndims == ndims2, status, MATLABMATRIXH_ELNTH, MATLABMATRIXH_MSGELNTH);

	for ( iterator = 0; iterator < ndims; iterator++)
	{
		ASSERT ( length->data[iterator] == ((ATYPE2*)(A))->dimLength->data[iterator], status, MATLABMATRIXH_ELNTH, MATLABMATRIXH_MSGELNTH);
	}

        LALDCreateArray( status->statusPtr, result, length);

	if ( ndims == 2 )
	{
		for( row = 0; row < length->data[0]; row++)
		{
			for( column = 0; column < length->data[1]; column++)
			{
                                myindex = (row * length->data[1]) + column;
				(*result)->data[myindex] = pow( (A->data[myindex]) , (B->data[myindex]));
			}
		}
	}
	else
	{
		LALDDestroyArray( status->statusPtr, result);
		(*result) = NULL;
	}

        LALU4DestroyVector( status->statusPtr, &length );

        DETATCHSTATUSPTR( status );
        RETURN (status);
}


