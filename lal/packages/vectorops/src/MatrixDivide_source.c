#define CONCAT2x(a,b) a##b
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a##b##c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define CONCAT4x(a,b,c,d) a##b##c##d
#define CONCAT4(a,b,c,d) CONCAT4x(a,b,c,d)
#define CONCAT5x(a,b,c,d,e) a##b##c##d##e
#define CONCAT5(a,b,c,d,e) CONCAT5x(a,b,c,d,e)
#define STRING(a) #a

#define VTYPE CONCAT2(TYPE,Vector)
#define ATYPE CONCAT2(TYPE,Array)
#define VTYPE2 CONCAT2(TYPE2,Vector)
#define ATYPE2 CONCAT2(TYPE2,Array)

#define CAFUNC CONCAT3(LAL,TYPECODE,CreateArray)
#define CVFUNC CONCAT3(LAL,TYPECODE,CreateVector)
#define DAFUNC CONCAT3(LAL,TYPECODE,DestroyArray)
#define CAFUNC2 CONCAT3(LAL,TYPECODE2,CreateArray)
#define CVFUNC2 CONCAT3(LAL,TYPECODE2,CreateVector)
#define DAFUNC2 CONCAT3(LAL,TYPECODE2,DestroyArray)

#define F1 CONCAT5(LAL,TYPECODE,DotSlash,TYPECODE2,Vector)
#define F2 CONCAT4(LAL,TYPECODE,VectorDotSlash,TYPECODE2)
#define F3 CONCAT5(LAL,TYPECODE,VectorDotSlash,TYPECODE2,Vector)
#define F4 CONCAT5(LAL,TYPECODE,DotSlash,TYPECODE2,Array)
#define F5 CONCAT4(LAL,TYPECODE,ArrayDotSlash,TYPECODE2)
#define F6 CONCAT5(LAL,TYPECODE,ArrayDotSlash,TYPECODE2,Array)

void F1 (
         LALStatus	*status,
         VTYPE2		**result,
         TYPE		A,
         VTYPE2		*B
         )
{
        /*  Variable Declarations  */
        UINT4   iterator;
        UINT4   length;

        INITSTATUS( status, STRING(F1) , MATLABMATRIXDIVC);
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

        CVFUNC2( status->statusPtr, result, length);

        for (iterator = 0; iterator < length; iterator++)
        {
                (*result)->data[iterator] = B->data[iterator] / A;
        }

        DETATCHSTATUSPTR( status );
        RETURN (status);
}

void F2 (
         LALStatus       *status,
         VTYPE2          **result,
         VTYPE2          *A,
         TYPE            B
         )
{
        /*  Variable Declarations  */
        UINT4   iterator;
        UINT4   length;

        INITSTATUS( status, STRING(F2) , MATLABMATRIXDIVC);
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

        CVFUNC2( status->statusPtr, result, length);

        for (iterator = 0; iterator < length; iterator++)
        {
               (*result)->data[iterator] = A->data[iterator] / B;
        }

        DETATCHSTATUSPTR( status );
        RETURN (status);
}

void F3 (
         LALStatus	*status,
         VTYPE		**result,
         VTYPE		*A,
         VTYPE2		*B
         )
{
        /*  Variable Declarations  */
        UINT4    iterator;
        UINT4   length;

        INITSTATUS( status, STRING(F3) , MATLABMATRIXDIVC);
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

	CVFUNC( status->statusPtr, result, length);

        for (iterator = 0; iterator < length; iterator++)
        {
                (*result)->data[iterator] = A->data[iterator] / B->data[iterator];
        }

        DETATCHSTATUSPTR( status );
        RETURN (status);
}

void F4 (
         LALStatus      *status,
         ATYPE2		**result,
         TYPE		A,
         ATYPE2		*B
         )
{
        /*  Variable Declarations  */
        UINT4Vector     *length;
        UINT4           ndims;
        UINT4		iterator, myindex;
	UINT4		row, column;

        INITSTATUS( status, STRING(F4) , MATLABMATRIXDIVC);
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

        CAFUNC2( status->statusPtr, result, length);

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
		DAFUNC2( status->statusPtr, result);
		(*result) = NULL;
        }


	LALU4DestroyVector( status->statusPtr, &length );

        DETATCHSTATUSPTR( status );
        RETURN (status);
}

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

        INITSTATUS( status, STRING(F5) , MATLABMATRIXDIVC);
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

        CAFUNC( status->statusPtr, result, length);

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
		DAFUNC( status->statusPtr, result);
		(*result) = NULL;
        }


        LALU4DestroyVector( status->statusPtr, &length );

        DETATCHSTATUSPTR( status );
        RETURN (status);
}

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

        INITSTATUS( status, STRING(F6) , MATLABMATRIXDIVC);
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

        CAFUNC( status->statusPtr, result, length);

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
		DAFUNC( status->statusPtr, result);
		(*result) = NULL;
	}

        LALU4DestroyVector( status->statusPtr, &length );
        DETATCHSTATUSPTR( status );
        RETURN (status);
}
