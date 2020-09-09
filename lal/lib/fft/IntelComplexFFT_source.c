#define CONCAT2x(a,b) a##b
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a##b##c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define STRING(a) #a

#ifdef SINGLE_PRECISION
#define COMPLEX_TYPE COMPLEX8
#define TYPESUFFIX f
#define DFTI_TYPE DFTI_SINGLE
#else
#define COMPLEX_TYPE COMPLEX16
#define TYPESUFFIX
#define DFTI_TYPE DFTI_DOUBLE
#endif

#define PLAN_TYPE			CONCAT2(COMPLEX_TYPE,FFTPlan)
#define COMPLEX_VECTOR_TYPE		CONCAT2(COMPLEX_TYPE,Vector)

#define CREATE_PLAN_FUNCTION		CONCAT2(XLALCreate,PLAN_TYPE)
#define CREATE_FORWARD_PLAN_FUNCTION	CONCAT2(XLALCreateForward,PLAN_TYPE)
#define CREATE_REVERSE_PLAN_FUNCTION	CONCAT2(XLALCreateReverse,PLAN_TYPE)
#define DESTROY_PLAN_FUNCTION		CONCAT2(XLALDestroy,PLAN_TYPE)
#define VECTOR_FFT_FUNCTION		CONCAT3(XLAL,COMPLEX_VECTOR_TYPE,FFT)

#define FFTWX				CONCAT2(fftw,TYPESUFFIX)
#define FFTWX_COMPLEX			CONCAT2(FFTWX,_complex)
#define FFTWX_PLAN_DFT_1D		CONCAT2(FFTWX,_plan_dft_1d)
#define FFTWX_DESTROY_PLAN		CONCAT2(FFTWX,_destroy_plan)
#define FFTWX_EXECUTE_DFT		CONCAT2(FFTWX,_execute_dft)

PLAN_TYPE *CREATE_PLAN_FUNCTION(UINT4 size, int fwdflg, __attribute__ ((unused)) int measurelvl)
{
    PLAN_TYPE *plan;
    INT8  fftStat;

    if (!size)
        XLAL_ERROR_NULL( XLAL_EBADLEN );

    /* allocate memory for the plan */
    plan = XLALMalloc( sizeof( *plan ) );

    if (!plan)
        XLAL_ERROR_NULL( XLAL_ENOMEM );

    /* make the intel fft descriptor */
    fftStat = DftiCreateDescriptor( &(plan->plan),
        DFTI_TYPE, DFTI_COMPLEX, 1, size);
    CHECKINTELFFTSTATUS_NULL( fftStat );

    /* configure intel fft descriptor */
    fftStat = DftiSetValue( plan->plan, DFTI_PLACEMENT,
        DFTI_NOT_INPLACE );
    CHECKINTELFFTSTATUS_NULL( fftStat );
    fftStat = DftiSetValue( plan->plan, DFTI_PACKED_FORMAT,
        DFTI_PACK_FORMAT );
    CHECKINTELFFTSTATUS_NULL( fftStat );

    /* commit the intel fft descriptor */
    fftStat = DftiCommitDescriptor( plan->plan );
    CHECKINTELFFTSTATUS_NULL( fftStat );

    /* now set remaining plan fields */
    plan->size = size;
    plan->sign = ( fwdflg ? -1 : 1 );

    return plan;
}

PLAN_TYPE *CREATE_FORWARD_PLAN_FUNCTION(UINT4 size, int measurelvl)
{
    PLAN_TYPE *plan;
    plan = CREATE_PLAN_FUNCTION(size, 1, measurelvl);
    if (!plan)
        XLAL_ERROR_NULL(XLAL_EFUNC);
    return plan;
}

PLAN_TYPE *CREATE_REVERSE_PLAN_FUNCTION(UINT4 size, int measurelvl)
{
    PLAN_TYPE *plan;
    plan = CREATE_PLAN_FUNCTION(size, 0, measurelvl);
    if (!plan)
        XLAL_ERROR_NULL(XLAL_EFUNC);
    return plan;
}

void DESTROY_PLAN_FUNCTION(PLAN_TYPE * plan)
{
    INT8  fftStat;

    if (!plan)
        XLAL_ERROR_VOID(XLAL_EFAULT);

    /* destroy intel fft descriptor */
    fftStat = DftiFreeDescriptor( &(plan->plan) );
    CHECKINTELFFTSTATUS_VOID( fftStat );

    XLALFree( plan );
}

int VECTOR_FFT_FUNCTION(COMPLEX_VECTOR_TYPE * _LAL_RESTRICT_ output, const COMPLEX_VECTOR_TYPE * _LAL_RESTRICT_ input,
    const PLAN_TYPE * plan)
{
    INT8  fftStat;

    /* sanity check on arguments */

    if (!output || !input || !plan)
        XLAL_ERROR(XLAL_EFAULT);
    if (!plan->plan || !plan->size)
        XLAL_ERROR(XLAL_EINVAL);
    if (!output->data || !input->data || output->data == input->data)
        XLAL_ERROR(XLAL_EINVAL);        /* note: must be out-of-place */
    if (output->length != plan->size || input->length != plan->size)
        XLAL_ERROR(XLAL_EBADLEN);

    /* complex intel fft */
    if ( plan->sign == -1 )
    {
        fftStat = DftiComputeForward( plan->plan, input->data, output->data );
        CHECKINTELFFTSTATUS( fftStat );
    }
    else if ( plan->sign == 1 )
    {
        fftStat = DftiComputeBackward( plan->plan, input->data, output->data );
        CHECKINTELFFTSTATUS( fftStat );
    }
    else
    {
        XLAL_ERROR( XLAL_EINVAL );
    }

    return 0;
}

#undef CONCAT2x
#undef CONCAT2
#undef CONCAT3x
#undef CONCAT3
#undef STRING

#undef COMPLEX_TYPE
#undef TYPESUFFIX
#undef DFTI_TYPE

#undef PLAN_TYPE
#undef COMPLEX_VECTOR_TYPE

#undef CREATE_PLAN_FUNCTION
#undef CREATE_FORWARD_PLAN_FUNCTION
#undef CREATE_REVERSE_PLAN_FUNCTION
#undef DESTROY_PLAN_FUNCTION
#undef VECTOR_FFT_FUNCTION

#undef FFTWX
#undef FFTWX_COMPLEX
#undef FFTWX_PLAN_DFT_1D
#undef FFTWX_DESTROY_PLAN
#undef FFTWX_EXECUTE_DFT
