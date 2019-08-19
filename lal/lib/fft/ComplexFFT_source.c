#define CONCAT2x(a,b) a##b
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a##b##c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define STRING(a) #a

#ifdef SINGLE_PRECISION
#define COMPLEX_TYPE COMPLEX8
#define TYPESUFFIX f
#else
#define COMPLEX_TYPE COMPLEX16
#define TYPESUFFIX
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

PLAN_TYPE *CREATE_PLAN_FUNCTION(UINT4 size, int fwdflg, int measurelvl)
{
    PLAN_TYPE *plan;
    COMPLEX_TYPE *tmp1;
    COMPLEX_TYPE *tmp2;
    size_t nbytes;
    int flags;

    if (!size)
        XLAL_ERROR_NULL(XLAL_EBADLEN);

    nbytes = size * sizeof(COMPLEX_TYPE);

    /* set fftw3 flags to perform requested degree of measurement */

#   ifdef LAL_FFTW3_MEMALIGN_ENABLED
    flags = 0;
#   else
    flags = FFTW_UNALIGNED;
#   endif

    switch (measurelvl) {
    case 0:    /* estimate */
        flags |= FFTW_ESTIMATE;
        break;
    default:   /* exhaustive measurement */
        flags |= FFTW_EXHAUSTIVE;
        /* fall-through */
    case 2:    /* lengthy measurement */
        flags |= FFTW_PATIENT;
        /* fall-through */
    case 1:    /* measure the best plan */
        flags |= FFTW_MEASURE;
        break;
    }

    /* allocate memory for the plan and the temporary arrays */

    plan = XLALMalloc(sizeof(*plan));
    if (!plan)
        XLAL_ERROR_NULL(XLAL_ENOMEM);

#   ifdef LAL_FFTW3_MEMALIGN_ENABLED
    tmp1 = XLALMallocAligned(nbytes);
    tmp2 = XLALMallocAligned(nbytes);
    if (!tmp1 || !tmp2) {
        XLALFreeAligned(tmp1);
        XLALFreeAligned(tmp2);
        XLAL_ERROR_NULL(XLAL_ENOMEM);
    }
#   else
    tmp1 = XLALMalloc(nbytes);
    tmp2 = XLALMalloc(nbytes);
    if (!tmp1 || !tmp2) {
        XLALFree(tmp1);
        XLALFree(tmp2);
        XLAL_ERROR_NULL(XLAL_ENOMEM);
    }
#   endif

    /* establish fftw mutex lock and create plan */

    LAL_FFTW_WISDOM_LOCK;
    plan->plan =
        FFTWX_PLAN_DFT_1D(size, (FFTWX_COMPLEX *) tmp1, (FFTWX_COMPLEX *) tmp2, fwdflg ? FFTW_FORWARD : FFTW_BACKWARD, flags);
    LAL_FFTW_WISDOM_UNLOCK;

    /* free the temporary arrays */

#   ifdef LAL_FFTW3_MEMALIGN_ENABLED
    XLALFreeAligned(tmp1);
    XLALFreeAligned(tmp2);
#   else
    XLALFree(tmp1);
    XLALFree(tmp2);
#   endif

    /* check to see success of plan creation */

    if (!plan->plan) {
        XLALFree(plan);
        XLAL_ERROR_NULL(XLAL_EFAILED);
    }

    /* set remaining plan fields */

    plan->size = size;
    plan->sign = (fwdflg ? -1 : 1);

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
    if (plan) {
        if (plan->plan) {
            LAL_FFTW_WISDOM_LOCK;
            FFTWX_DESTROY_PLAN(plan->plan);
            LAL_FFTW_WISDOM_UNLOCK;
        }
        memset(plan, 0, sizeof(*plan));
        XLALFree(plan);
    }
}

int VECTOR_FFT_FUNCTION(COMPLEX_VECTOR_TYPE * _LAL_RESTRICT_ output, const COMPLEX_VECTOR_TYPE * _LAL_RESTRICT_ input,
    const PLAN_TYPE * plan)
{
    COMPLEX_TYPE *input_data;
    COMPLEX_TYPE *output_data;

    /* sanity check on arguments */

    if (!output || !input || !plan)
        XLAL_ERROR(XLAL_EFAULT);
    if (!plan->plan || !plan->size)
        XLAL_ERROR(XLAL_EINVAL);
    if (!output->data || !input->data || output->data == input->data)
        XLAL_ERROR(XLAL_EINVAL);        /* note: must be out-of-place */
    if (output->length != plan->size || input->length != plan->size)
        XLAL_ERROR(XLAL_EBADLEN);

    input_data = input->data;
    output_data = output->data;

    /* if memory alignment is required, check memory alignment and create
     * temporary space if necessary */

#   ifdef LAL_FFTW3_MEMALIGN_ENABLED
    if (!LAL_IS_MEMORY_ALIGNED(input_data)) {
        input_data = XLALMallocAligned(plan->size * sizeof(COMPLEX_TYPE));
        if (!input_data)
            XLAL_ERROR(XLAL_ENOMEM);
        memcpy(input_data, input->data, plan->size * sizeof(COMPLEX_TYPE));
    }
    if (!LAL_IS_MEMORY_ALIGNED(output_data)) {
        output_data = XLALMallocAligned(plan->size * sizeof(COMPLEX_TYPE));
        if (!output_data) {
            if (input_data != input->data)
                XLALFreeAligned(input_data);
            XLAL_ERROR(XLAL_ENOMEM);
        }
    }
#   endif

    /* perform the fft */

    FFTWX_EXECUTE_DFT(plan->plan, (FFTWX_COMPLEX *)input_data, (FFTWX_COMPLEX *)output_data);

    /* cleanup aligned memory space if memory alignment is required;
     * copy data from temporary space to output vector */

#   ifdef LAL_FFTW3_MEMALIGN_ENABLED
    if (input_data != input->data)
        XLALFreeAligned(input_data);
    if (output_data != output->data) {
        memcpy(output->data, output_data, plan->size * sizeof(COMPLEX_TYPE));
        XLALFreeAligned(output_data);
    }
#   endif

    return 0;
}

#undef CONCAT2x
#undef CONCAT2
#undef CONCAT3x
#undef CONCAT3
#undef STRING

#undef COMPLEX_TYPE
#undef TYPESUFFIX

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
