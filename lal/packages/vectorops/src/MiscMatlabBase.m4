dnl $Id$
ifelse(TYPECODE,`D',`define(`TYPE',`REAL8')')
ifelse(TYPECODE,`S',`define(`TYPE',`REAL4')')
ifelse(TYPECODE,`I2',`define(`TYPE',`INT2')')
ifelse(TYPECODE,`I4',`define(`TYPE',`INT4')')
ifelse(TYPECODE,`I8',`define(`TYPE',`INT8')')
define(`VTYPE',`format(`%sVector',TYPE)')
define(`F1',`format(`LAL%sCumSum',TYPECODE)')
define(`F2',`format(`LAL%sSum',TYPECODE)')
define(`F3',`format(`LAL%sMax',TYPECODE)')
define(`F4',`format(`LAL%sFlipVector',TYPECODE)')

/******************************* <lalLaTeX file="MiscMatlabC">
\begin{verbatim}void F1 ( LALStatus *status, VTYPE **result, VTYPE *data )\end{verbatim}
*************************************************** </lalLaTeX> */

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

	INITSTATUS( status, "F1" , MATLABMATRIXSUMC);
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

	LAL`'TYPECODE`'CreateVector( status->statusPtr, result, length);
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

/******************************* <lalLaTeX file="MiscMatlabC">
\begin{verbatim}void F2 ( LALStatus *status, TYPE *result, VTYPE *data )\end{verbatim}
*************************************************** </lalLaTeX> */

void F2 (
        LALStatus		*status,
        TYPE			*result,
        VTYPE		*data
)
{
	/*  Variable Declarations  */
        INT4    iterator;
        INT4    length;

        INITSTATUS( status, "F2" , MATLABMATRIXSUMC);

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

/******************************* <lalLaTeX file="MiscMatlabC">
\begin{verbatim}void F3 (LALStatus *status, TYPE *result, VTYPE *data, INT4 *myindex )\end{verbatim}
*************************************************** </lalLaTeX> */

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

        INITSTATUS( status, "F3" , MATLABMATRIXSUMC);

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

/******************************* <lalLaTeX file="MiscMatlabC">
\begin{verbatim}void F4 (LALStatus *status, VTYPE **result, VTYPE *data )\end{verbatim}
*************************************************** </lalLaTeX> */

void F4 (
        LALStatus		*status,
        VTYPE		**result,
        VTYPE           *data
)
{
        /*  Variable Declarations  */
        INT4    iterator;
        INT4    length;

        INITSTATUS( status, "F4" , MATLABMATRIXSUMC);
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

        LAL`'TYPECODE`'CreateVector( status->statusPtr, result, length);

	for (iterator = 0; iterator < length; iterator++)
	{
		(*result)->data[iterator] = ((VTYPE*)(data))->data[length - iterator - 1];
	}

	DETATCHSTATUSPTR( status );
        RETURN( status );
}


