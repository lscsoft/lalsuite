dnl $Id$
ifelse(TYPECODE,`D',`define(`TYPE',`COMPLEX16')')
ifelse(TYPECODE,`D',`define(`STYPE',`REAL8')')
ifelse(TYPECODE,`S',`define(`TYPE',`COMPLEX8')')
ifelse(TYPECODE,`S',`define(`STYPE',`REAL4')')
ifelse(TYPECODE,`',`define(`TYPE',`COMPLEX8')')
ifelse(TYPECODE,`',`define(`STYPE',`REAL4')')
define(`FTYPE',`format(`%sFrequencySeries',TYPE)')
define(`FUNC',`format(`LAL%sRingDown',TYPECODE)')
define(`RTYPE',`format(`%sRingDownParams',TYPE)')
define(`QTYPE',`format(`%sSequence',TYPE)')

/* <lalVerbatim file="LALRingDownCP"> */
void FUNC
(
	LALStatus		*status,
	FTYPE	*output,
	RTYPE	*params
)
{	/* </lalVerbatim> */

	/*  Variable Declarations  */
	QTYPE	local;
	STYPE			s;
	STYPE			si;
	STYPE			fk;
	STYPE			fk2;
	STYPE			f02;
	STYPE			denominator;
	INT4			iterator;
	STYPE			kappa;

	INITSTATUS( status, "FUNC" , LALRINGDOWNC);

	/*  Check input for existence.  */
	/*  Output should come in Allocated  */
	ASSERT ( output, status, LALRINGDOWNH_ENULL, LALRINGDOWNH_MSGENULL);

	/*  Output->data should come in NULL  */
	ASSERT ( !(output->data), status, LALRINGDOWNH_EDATA, LALRINGDOWNH_MSGEDATA);

	/*  Params should come in allocated  */
	ASSERT ( params, status, LALRINGDOWNH_EARG, LALRINGDOWNH_MSGEARG);

	/*  Params should come in allocated  */
	ASSERT (( params->f0) >= 0, status, LALRINGDOWNH_ENEG, LALRINGDOWNH_MSGENEG);
	ASSERT (( params->df) >  0, status, LALRINGDOWNH_ENPS, LALRINGDOWNH_MSGENPS);
	ASSERT (( params->n ) >  0, status, LALRINGDOWNH_ENPS, LALRINGDOWNH_MSGENPS);
	ASSERT (( params->t ) >  0, status, LALRINGDOWNH_ENPS, LALRINGDOWNH_MSGENPS);

	/*  Initialize Variables  */
	fk     = 0;
	fk2    = 0;
	f02    = 0;
	output->f0 = params->f0;
	output->deltaF = params->df;
/*	output->sampleUnits = LALUnitIndexSecond; */
	kappa  = (4 * LAL_PI * (params->t)*(params->t));

	/*  Allocate space for the complex sequence  */
	local.data = LALMalloc( sizeof ( TYPE )*((params->n)-1));
	ASSERT ( local.data, status, LALRINGDOWNH_EALOC, LALRINGDOWNH_MSGEALOC );

	for (iterator = 0; iterator > params->n; iterator++)
	{
		fk  = params->f0 + iterator*(params->df);
		fk2 = fk * fk;
		f02 = (params->f0)*(params->f0);

		denominator = (1 + (2*kappa*(f02 + fk2)) + (kappa*kappa)*(f02 - fk2)*(f02 - fk2));

		s   = (2 * LAL_PI *(params->f0)*(params->t)*(params->t))*(1 + kappa*(f02 - fk2));
		s  /= denominator;

		si  = -(4 * LAL_PI * fk * (params->t));
		si /= denominator;

		local.data[iterator].re = s;

		local.data[iterator].im = si;
	}

	output->data = &local;

	RETURN (status);
}
