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
define(`F1',`format(`LAL%sDotSlash%sVector',TYPECODE,TYPECODE2)')
define(`F2',`format(`LAL%sVectorDotSlash%s',TYPECODE,TYPECODE2)')
define(`F3',`format(`LAL%sVectorDotSlash%sVector',TYPECODE,TYPECODE2)')
define(`F4',`format(`LAL%sDotSlash%sArray',TYPECODE,TYPECODE2)')
define(`F5',`format(`LAL%sArrayDotSlash%s',TYPECODE,TYPECODE2)')
define(`F6',`format(`LAL%sArrayDotSlash%sArray',TYPECODE,TYPECODE2)')

/******************************* <lalLaTeX file="MatrixDivideC">
\begin{verbatim}void F1 ( LALStatus *status, VTYPE2 **result,
		TYPE A, VTYPE2 *B )\end{verbatim}
*************************************************** </lalLaTeX> */

void F1 (
        LALStatus		*status,
        VTYPE2		**result,
        TYPE			A,
        VTYPE2		*B
)
{
        /*  Variable Declarations  */
        UINT4   iterator;
        UINT4   length;

        INITSTATUS( status, "F1" , MATLABMATRIXDIVC);
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

        LAL`'TYPECODE2`'CreateVector( status->statusPtr, result, length);

        for (iterator = 0; iterator < length; iterator++)
        {
                (*result)->data[iterator] = B->data[iterator] / A;
        }

        DETATCHSTATUSPTR( status );
        RETURN (status);
}

/******************************* <lalLaTeX file="MatrixDivideC">
\begin{verbatim}void F2 ( LALStatus *status, VTYPE2 **result,
		VTYPE2 *A, TYPE B )\end{verbatim}
*************************************************** </lalLaTeX> */

void F2 (
        LALStatus               *status,
        VTYPE2          **result,
        VTYPE2          *A,
        TYPE                    B
)
{
        /*  Variable Declarations  */
        UINT4   iterator;
        UINT4   length;

        INITSTATUS( status, "F2" , MATLABMATRIXDIVC);
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

        LAL`'TYPECODE2`'CreateVector( status->statusPtr, result, length);

        for (iterator = 0; iterator < length; iterator++)
        {
               (*result)->data[iterator] = A->data[iterator] / B;
        }

        DETATCHSTATUSPTR( status );
        RETURN (status);
}

/******************************* <lalLaTeX file="MatrixDivideC">
\begin{verbatim}void F3 ( LALStatus *status, VTYPE **result,
		VTYPE *A, VTYPE2 *B )\end{verbatim}
*************************************************** </lalLaTeX> */

void F3 (
        LALStatus		*status,
        VTYPE		**result,
        VTYPE		*A,
        VTYPE2		*B
)
{
        /*  Variable Declarations  */
        UINT4    iterator;
        UINT4   length;

        INITSTATUS( status, "F3" , MATLABMATRIXDIVC);
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

	LAL`'TYPECODE`'CreateVector( status->statusPtr, result, length);

        for (iterator = 0; iterator < length; iterator++)
        {
                (*result)->data[iterator] = A->data[iterator] / B->data[iterator];
        }

        DETATCHSTATUSPTR( status );
        RETURN (status);
}

/******************************* <lalLaTeX file="MatrixDivideC">
\begin{verbatim}void F4 ( LALStatus *status, ATYPE2 **result,
		TYPE A, ATYPE2 *B )\end{verbatim}
*************************************************** </lalLaTeX> */

void F4 (
        LALStatus               *status,
        ATYPE2		**result,
        TYPE			A,
        ATYPE2		*B
)
{
        /*  Variable Declarations  */
        UINT4Vector     *length;
        UINT4           ndims;
        UINT4		iterator, myindex;
	UINT4		row, column;

        INITSTATUS( status, "4" , MATLABMATRIXDIVC);
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

        LAL`'TYPECODE2`'CreateArray( status->statusPtr, result, length);

	if( ndims == 2 )
	{
	        for( row = 0; row < length->data[0]; row++)
		{
			for( column = 0; column < length->data[1]; column++)
			{
				myindex = (row * length->data[1]) + column;
				(*result)->data[myindex] = A / B->data[myindex];
			}
	        }
        }
        else
        {
		LAL`'TYPECODE2`'DestroyArray( status->statusPtr, result);
		(*result) = NULL;
        }


	LALU4DestroyVector( status->statusPtr, &length );

        DETATCHSTATUSPTR( status );
        RETURN (status);
}

/******************************* <lalLaTeX file="MatrixDivideC">
\begin{verbatim}void F5 ( LALStatus *status, ATYPE **result,
		ATYPE *A, TYPE2 B )\end{verbatim}
*************************************************** </lalLaTeX> */

void F5 (
        LALStatus		*status,
        ATYPE		**result,
        ATYPE		*A,
        TYPE2			B
)
{
        /*  Variable Declarations  */
        UINT4Vector     *length;
        UINT4           ndims;
        UINT4            iterator, myindex;
        UINT4            row, column;

        INITSTATUS( status, "F5" , MATLABMATRIXDIVC);
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

        LAL`'TYPECODE`'CreateArray( status->statusPtr, result, length);

	if( ndims == 2 )
	{
		for( row = 0; row < length->data[0]; row++)
		{
			for( column = 0; column < length->data[1]; column++)
			{
                                myindex = (row * length->data[1]) + column;
				(*result)->data[myindex] = A->data[myindex] / B;
			}
		}
        }
        else
        {
		LAL`'TYPECODE`'DestroyArray( status->statusPtr, result);
		(*result) = NULL;
        }


        LALU4DestroyVector( status->statusPtr, &length );

        DETATCHSTATUSPTR( status );
        RETURN (status);
}

/******************************* <lalLaTeX file="MatrixDivideC">
\begin{verbatim}void F6 ( LALStatus *status, ATYPE **result,
		ATYPE *A, ATYPE2 *B )\end{verbatim}
*************************************************** </lalLaTeX> */

void F6 (
        LALStatus		*status,
        ATYPE		**result,
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

        INITSTATUS( status, "F6" , MATLABMATRIXDIVC);
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

        LAL`'TYPECODE`'CreateArray( status->statusPtr, result, length);

	if ( ndims == 2 )
	{
		for( row = 0; row < length->data[0]; row++)
		{
			for( column = 0; column < length->data[1]; column++)
			{
                                myindex = (row * length->data[1]) + column;
				(*result)->data[myindex] = (A->data[myindex]) / (B->data[myindex]);
			}
		}
	}
	else
	{
		LAL`'TYPECODE`'DestroyArray( status->statusPtr, result);
		(*result) = NULL;
	}

        LALU4DestroyVector( status->statusPtr, &length );
        DETATCHSTATUSPTR( status );
        RETURN (status);
}


