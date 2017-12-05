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

/*
 * Legacy Routines
 */

#define OLD_CREATE_FORWARD_PLAN_FUNCTION	CONCAT2(LALCreateForward,PLAN_TYPE)
#define OLD_CREATE_REVERSE_PLAN_FUNCTION	CONCAT2(LALCreateReverse,PLAN_TYPE)
#define OLD_DESTROY_PLAN_FUNCTION		CONCAT2(LALDestroy,PLAN_TYPE)
#define OLD_VECTOR_FFT_FUNCTION			CONCAT3(LAL,COMPLEX_VECTOR_TYPE,FFT)

void OLD_CREATE_FORWARD_PLAN_FUNCTION(LALStatus * status, PLAN_TYPE ** plan, UINT4 size, INT4 measure)
{
    INITSTATUS(status);
    XLAL_PRINT_DEPRECATION_WARNING(STRING(CREATE_FORWARD_PLAN_FUNCTION));

    ASSERT(plan, status, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL);
    ASSERT(!*plan, status, COMPLEXFFTH_ENNUL, COMPLEXFFTH_MSGENNUL);
    ASSERT(size > 0, status, COMPLEXFFTH_ESIZE, COMPLEXFFTH_MSGESIZE);

    *plan = CREATE_FORWARD_PLAN_FUNCTION(size, measure);
    if (!*plan) {
        int code = xlalErrno;
        XLALClearErrno();
        switch (code) {
        case XLAL_EBADLEN:
            ABORT(status, COMPLEXFFTH_ESIZE, COMPLEXFFTH_MSGESIZE);
        case XLAL_ENOMEM:
            ABORT(status, COMPLEXFFTH_EALOC, COMPLEXFFTH_MSGEALOC);
        case XLAL_EFAILED:
            ABORT(status, COMPLEXFFTH_EFFTW, COMPLEXFFTH_MSGEFFTW);
        default:
            ABORTXLAL(status);
        }
    }

    RETURN(status);
}

void OLD_CREATE_REVERSE_PLAN_FUNCTION(LALStatus * status, PLAN_TYPE ** plan, UINT4 size, INT4 measure)
{
    INITSTATUS(status);
    XLAL_PRINT_DEPRECATION_WARNING(STRING(CREATE_REVERSE_PLAN_FUNCTION));

    ASSERT(plan, status, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL);
    ASSERT(!*plan, status, COMPLEXFFTH_ENNUL, COMPLEXFFTH_MSGENNUL);
    ASSERT(size > 0, status, COMPLEXFFTH_ESIZE, COMPLEXFFTH_MSGESIZE);

    *plan = CREATE_REVERSE_PLAN_FUNCTION(size, measure);
    if (!*plan) {
        int code = xlalErrno;
        XLALClearErrno();
        switch (code) {
        case XLAL_EBADLEN:
            ABORT(status, COMPLEXFFTH_ESIZE, COMPLEXFFTH_MSGESIZE);
        case XLAL_ENOMEM:
            ABORT(status, COMPLEXFFTH_EALOC, COMPLEXFFTH_MSGEALOC);
        case XLAL_EFAILED:
            ABORT(status, COMPLEXFFTH_EFFTW, COMPLEXFFTH_MSGEFFTW);
        default:
            ABORTXLAL(status);
        }
    }

    RETURN(status);
}

void OLD_DESTROY_PLAN_FUNCTION(LALStatus * status, PLAN_TYPE ** plan)
{
    INITSTATUS(status);
    XLAL_PRINT_DEPRECATION_WARNING(STRING(DESTROY_PLAN_FUNCTION));
    ASSERT(plan, status, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL);
    ASSERT(*plan, status, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL);
    DESTROY_PLAN_FUNCTION(*plan);
    if (xlalErrno) {
        int code = xlalErrno;
        XLALClearErrno();
        switch (code) {
        case XLAL_EINVAL:
            ABORT(status, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL);
        default:
            ABORTXLAL(status);
        }
    }
    RETURN(status);
}

void OLD_VECTOR_FFT_FUNCTION(LALStatus * status, COMPLEX_VECTOR_TYPE * output, COMPLEX_VECTOR_TYPE * input, PLAN_TYPE * plan)
{
    int code;
    INITSTATUS(status);
    XLAL_PRINT_DEPRECATION_WARNING(STRING(VECTOR_FFT_FUNCTION));

    ASSERT(output, status, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL);
    ASSERT(input, status, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL);
    ASSERT(plan, status, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL);

    ASSERT(output->data, status, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL);
    ASSERT(input->data, status, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL);

    /* make sure that it is not the same data! */
    ASSERT(output->data != input->data, status, COMPLEXFFTH_ESAME, COMPLEXFFTH_MSGESAME);

    /* make sure that the lengths agree */
    ASSERT(plan->size > 0, status, COMPLEXFFTH_ESIZE, COMPLEXFFTH_MSGESIZE);
    ASSERT(output->length == plan->size, status, COMPLEXFFTH_ESZMM, COMPLEXFFTH_MSGESZMM);
    ASSERT(input->length == plan->size, status, COMPLEXFFTH_ESZMM, COMPLEXFFTH_MSGESZMM);

    code = VECTOR_FFT_FUNCTION(output, input, plan);
    if (code) {
        code = xlalErrno;
        XLALClearErrno();
        switch (code) {
        case XLAL_EINVAL:
            if (!plan->size) {  /* plan size was invalid */
                ABORT(status, COMPLEXFFTH_ESIZE, COMPLEXFFTH_MSGESIZE);
            } else if (output->data == input->data) {   /* same data pointers */
                ABORT(status, COMPLEXFFTH_ESAME, COMPLEXFFTH_MSGESAME);
            } else {    /* one of the data pointers was NULL */

                ABORT(status, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL);
            }
        case XLAL_EBADLEN:     /* size mismatch */
            ABORT(status, COMPLEXFFTH_ESZMM, COMPLEXFFTH_MSGESZMM);
        default:
            ABORTXLAL(status);
        }
    }

    RETURN(status);
}

#undef OLD_CREATE_FORWARD_PLAN_FUNCTION
#undef OLD_CREATE_REVERSE_PLAN_FUNCTION
#undef OLD_DESTROY_PLAN_FUNCTION
#undef OLD_VECTOR_FFT_FUNCTION

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
