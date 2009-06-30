dnl $Id$
ifelse(TYPECODE,`D',`define(`TYPE',`REAL8')')
ifelse(TYPECODE,`S',`define(`TYPE',`REAL4')')
ifelse(TYPECODE,`',`define(`TYPE',`REAL4')')
define(`STYPE',`format(`%sSequence',TYPE)')
define(`FUNC',`format(`LAL%sMoment',TYPECODE)')


/* <lalVerbatim file="LALMomentCP"> */
void FUNC (
	LALStatus		*status,
	TYPE			*result,
	STYPE		*data,
	INT4			whichMoment
)
{	/* </lalVerbatim> */

	/*  Variable Declarations  */
	INT4	iterator;
	INT4	length;
	TYPE	ave 	= 0.0;
	TYPE	momentn	= 0.0;
	TYPE	sum	= 0.0;
	TYPE	base	= 0.0;

	INITSTATUS( status, "FUNC" , LALMOMENTC);

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
