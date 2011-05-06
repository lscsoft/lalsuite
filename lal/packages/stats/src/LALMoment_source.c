#define CONCAT2x(a,b) a##b
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a##b##c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define STRING(a) #a

#define STYPE CONCAT2(TYPE,Sequence)
#define FUNC CONCAT3(LAL,TYPECODE,Moment)


void FUNC (
	LALStatus		*status,
	TYPE			*result,
	STYPE		*data,
	INT4			whichMoment
)
{	

	/*  Variable Declarations  */
	INT4	iterator;
	INT4	length;
	TYPE	ave 	= 0.0;
	TYPE	momentn	= 0.0;
	TYPE	sum	= 0.0;
	TYPE	base	= 0.0;

	INITSTATUS( status, STRING(FUNC) , LALMOMENTC);

	/*  Check input for existence.  */
	/*  Result should come in Allocated  */
	ASSERT ( result, status, LALMOMENTH_ENULL, LALMOMENTH_MSGENULL);

	/*  whichMoment must be greater than 1  */
	ASSERT ( whichMoment >= 1, status, LALMOMENTH_ENULL, LALMOMENTH_MSGENULL);

	/*  data must be defined  */
	ASSERT ( data, status, LALMOMENTH_ENULL, LALMOMENTH_MSGENULL);

	length = ((REAL8Sequence*)(data))->length;

        /*  length must be greater than one  */
        ASSERT ( length > 1, status, LALMOMENTH_ELNTH, LALMOMENTH_MSGELNTH);

	for (iterator = 0; iterator < length; iterator++)
	{
		sum += ((STYPE*)(data))->data[iterator];
	}

	ave = ( sum / length );

	/*  Return the Mean if whichMoment = 1  */
	if ( whichMoment == 1 )
	{
		*result = ave;
		RETURN(status);
	}

	for (iterator = 0; iterator < length; iterator++)
	{
		base = (data->data[iterator] - ave);
		sum = pow( base, whichMoment );
		momentn += sum;
	}

	momentn /= ((TYPE)(length - 1));


	if ( whichMoment > 2 )
	{
		for (iterator = 0; iterator < length; iterator++)
		{
			base = (data->data[iterator] - ave);
			sum = pow( base, whichMoment );
			momentn += sum;
		}

		momentn /= ((TYPE)(length));
	}


	*result = momentn;

	RETURN (status);
}
