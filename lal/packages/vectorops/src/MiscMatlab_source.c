#define CONCAT2x(a,b) a##b
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a##b##c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define STRING(a) #a

#define VTYPE CONCAT2(TYPE,Vector)

#define CVFUNC CONCAT3(LAL,TYPECODE,CreateVector)

#define F1 CONCAT3(LAL,TYPECODE,CumSum)
#define F2 CONCAT3(LAL,TYPECODE,Sum)
#define F3 CONCAT3(LAL,TYPECODE,Max)
#define F4 CONCAT3(LAL,TYPECODE,FlipVector)

void F1 (
	LALStatus		*status,
	VTYPE		**result,
	VTYPE		*data
)
{
	/*  Variable Declarations  */
	INT4	iterator;
	INT4	myindex;
	INT4	length;

	INITSTATUS( status, STRING(F1) , MATLABMATRIXSUMC);
        ATTATCHSTATUSPTR( status );

	/*  Check input for existence.  */
	/*  Result should come in NULL  */
        ASSERT ( result, status, MATLABMATRIXH_ENULL, MATLABMATRIXH_MSGENULL);
        ASSERT ( !(*result), status, MATLABMATRIXH_ENNUL, MATLABMATRIXH_MSGENNUL);

	/*  data must be defined  */
	ASSERT ( data, status, MATLABMATRIXH_ENULL, MATLABMATRIXH_MSGENULL);

	length = ((VTYPE*)(data))->length;

	/*  length must be greater than one  */
	ASSERT ( length >= 1, status, MATLABMATRIXH_ELNTH, MATLABMATRIXH_MSGELNTH);

	CVFUNC( status->statusPtr, result, length);
	CHECKSTATUSPTR( status );

	myindex = 0;

	while ( myindex < length )
	{
		(*result)->data[myindex] = 0;
		for (iterator = -1; iterator < myindex; iterator++)
		{
			(*result)->data[myindex] += ((VTYPE*)(data))->data[iterator+1];
		}

		myindex++;
	}

        DETATCHSTATUSPTR( status );
	RETURN( status );
}

void F2 (
        LALStatus		*status,
        TYPE			*result,
        VTYPE		*data
)
{
	/*  Variable Declarations  */
        INT4    iterator;
        INT4    length;

        INITSTATUS( status, STRING(F2) , MATLABMATRIXSUMC);

        /*  Check input for existence.  */
        /*  data must be defined  */
        ASSERT ( data, status, MATLABMATRIXH_ENULL, MATLABMATRIXH_MSGENULL);

        length = (( VTYPE* )(data))->length;

        /*  length must be greater than one  */
        ASSERT ( length >= 1, status, MATLABMATRIXH_ELNTH, MATLABMATRIXH_MSGELNTH);

	(*result) = 0;

	for (iterator = 0; iterator < length; iterator++)
	{
		(*result) += ((VTYPE*)(data))->data[iterator];
	}

        RETURN( status );
}

void F3 (
        LALStatus		*status,
        TYPE			*result,
        VTYPE		*data,
	INT4			*myindex
)
{
        /*  Variable Declarations  */
        INT4    iterator;
        INT4    length;

        INITSTATUS( status, STRING(F3) , MATLABMATRIXSUMC);

        /*  Check input for existence.  */
        /*  data must be defined  */
        ASSERT ( data, status, MATLABMATRIXH_ENULL, MATLABMATRIXH_MSGENULL);

	length = (( VTYPE* )(data))->length;

        /*  length must be greater than one  */
        ASSERT ( length >= 1, status, MATLABMATRIXH_ELNTH, MATLABMATRIXH_MSGELNTH);

	(*result) = ((VTYPE*)(data))->data[length - 1];
	(*myindex) = length - 1;

        for (iterator = 0; iterator < length; iterator++)
        {
                if ( ((VTYPE*)(data))->data[iterator] > (*result) )
		{
			(*result) = ((VTYPE*)(data))->data[iterator];
			(*myindex) = iterator;
		}
        }

        RETURN( status );
}

void F4 (
        LALStatus		*status,
        VTYPE		**result,
        VTYPE           *data
)
{
        /*  Variable Declarations  */
        INT4    iterator;
        INT4    length;

        INITSTATUS( status, STRING(F4) , MATLABMATRIXSUMC);
	ATTATCHSTATUSPTR( status );

        /*  Check input for existence.  */
        /*  Result should come in NULL  */
        ASSERT ( result, status, MATLABMATRIXH_ENULL, MATLABMATRIXH_MSGENULL);
        ASSERT ( !(*result), status, MATLABMATRIXH_ENNUL, MATLABMATRIXH_MSGENNUL);

        /*  data must be defined  */
        ASSERT ( data, status, MATLABMATRIXH_ENULL, MATLABMATRIXH_MSGENULL);

        length = ((VTYPE*)(data))->length;

        /*  length must be greater than one  */
        ASSERT ( length >= 1, status, MATLABMATRIXH_ELNTH, MATLABMATRIXH_MSGELNTH);

        CVFUNC( status->statusPtr, result, length);

	for (iterator = 0; iterator < length; iterator++)
	{
		(*result)->data[iterator] = ((VTYPE*)(data))->data[length - iterator - 1];
	}

	DETATCHSTATUSPTR( status );
        RETURN( status );
}


