#define CONCAT2x(a,b) a##b
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a##b##c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define STRING(a) #a

#ifdef SINGLE_PRECISION
#define REAL_TYPE REAL4
#define COMPLEX_TYPE COMPLEX8
#define TYPESUFFIX f
#else
#define REAL_TYPE REAL8
#define COMPLEX_TYPE COMPLEX16
#define TYPESUFFIX
#endif

#define PLAN_TYPE			CONCAT2(REAL_TYPE,FFTPlan)
#define REAL_VECTOR_TYPE		CONCAT2(REAL_TYPE,Vector)
#define COMPLEX_VECTOR_TYPE		CONCAT2(COMPLEX_TYPE,Vector)

#define CREATE_PLAN_FUNCTION		CONCAT2(XLALCreate,PLAN_TYPE)
#define CREATE_FORWARD_PLAN_FUNCTION	CONCAT2(XLALCreateForward,PLAN_TYPE)
#define CREATE_REVERSE_PLAN_FUNCTION	CONCAT2(XLALCreateReverse,PLAN_TYPE)
#define DESTROY_PLAN_FUNCTION		CONCAT2(XLALDestroy,PLAN_TYPE)
#define FORWARD_FFT_FUNCTION		CONCAT3(XLAL,REAL_TYPE,ForwardFFT)
#define REVERSE_FFT_FUNCTION		CONCAT3(XLAL,REAL_TYPE,ReverseFFT)
#define VECTOR_FFT_FUNCTION		CONCAT3(XLAL,REAL_VECTOR_TYPE,FFT)
#define POWER_SPECTRUM_FUNCTION		CONCAT3(XLAL,REAL_TYPE,PowerSpectrum)

#define CREALX				CONCAT2(creal,TYPESUFFIX)
#define CIMAGX				CONCAT2(cimag,TYPESUFFIX)
#define FFTWX				CONCAT2(fftw,TYPESUFFIX)
#define FFTWX_PLAN_R2R_1D		CONCAT2(FFTWX,_plan_r2r_1d)
#define FFTWX_DESTROY_PLAN		CONCAT2(FFTWX,_destroy_plan)
#define FFTWX_EXECUTE_R2R		CONCAT2(FFTWX,_execute_r2r)

PLAN_TYPE *CREATE_PLAN_FUNCTION(UINT4 size, int fwdflg, int measurelvl)
{
    PLAN_TYPE *plan;
    REAL_TYPE *tmp1;
    REAL_TYPE *tmp2;
    size_t nbytes;
    int flags;

    if (!size)
        XLAL_ERROR_NULL(XLAL_EBADLEN);

    nbytes = size * sizeof(REAL_TYPE);

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
    if (fwdflg) /* forward */
        plan->plan = FFTWX_PLAN_R2R_1D(size, tmp1, tmp2, FFTW_R2HC, flags);
    else        /* reverse */
        plan->plan = FFTWX_PLAN_R2R_1D(size, tmp1, tmp2, FFTW_HC2R, flags);
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

int FORWARD_FFT_FUNCTION(COMPLEX_VECTOR_TYPE * output, const REAL_VECTOR_TYPE * input, const PLAN_TYPE * plan)
{
    REAL_TYPE *input_data;
    REAL_TYPE *tmp;
    UINT4 k;
    size_t nbytes;

    /* sanity checks on arguments */

    if (!output || !input || !plan)
        XLAL_ERROR(XLAL_EFAULT);
    if (!plan->plan || !plan->size || plan->sign != -1)
        XLAL_ERROR(XLAL_EINVAL);
    if (!output->data || !input->data)
        XLAL_ERROR(XLAL_EINVAL);
    if (input->length != plan->size || output->length != plan->size / 2 + 1)
        XLAL_ERROR(XLAL_EBADLEN);

    nbytes = plan->size * sizeof(REAL_TYPE);
    input_data = input->data;

    /* create temporary storage space; make sure that input data is
     * aligned, if memory alignment is required */

#   ifdef LAL_FFTW3_MEMALIGN_ENABLED
    tmp = XLALMallocAligned(nbytes);
    if (!tmp)
        XLAL_ERROR(XLAL_ENOMEM);
    if (!LAL_IS_MEMORY_ALIGNED(input_data)) {
        /* need to create temporary aligned space for input data */
        input_data = XLALMallocAligned(nbytes);
        if (!input_data) {
            XLALFreeAligned(tmp);
            XLAL_ERROR(XLAL_ENOMEM);
        }
        memcpy(input_data, input->data, nbytes);
    }
#   else
    tmp = XLALMalloc(nbytes);
    if (!tmp)
        XLAL_ERROR(XLAL_ENOMEM);
#   endif

    /* perform the fft */

    FFTWX_EXECUTE_R2R(plan->plan, input_data, tmp);

    /* unpack the results into the output vector */

    /* dc component */
    output->data[0] = tmp[0];

    /* other components */
    for (k = 1; k < (plan->size + 1) / 2; ++k)  /* k < size/2 rounded up */
        output->data[k] = tmp[k] + I * tmp[plan->size - k];

    /* Nyquist frequency */
    if (plan->size % 2 == 0)    /* n is even */
        output->data[plan->size / 2] = tmp[plan->size / 2];

    /* cleanup and return */

#   ifdef LAL_FFTW3_MEMALIGN_ENABLED
    if (input_data != input->data)
        XLALFreeAligned(input_data);
    XLALFreeAligned(tmp);
#   else
    XLALFree(tmp);
#   endif

    return 0;
}

int REVERSE_FFT_FUNCTION(REAL_VECTOR_TYPE * output, const COMPLEX_VECTOR_TYPE * input, const PLAN_TYPE * plan)
{
    REAL_TYPE *output_data;
    REAL_TYPE *tmp;
    UINT4 k;
    size_t nbytes;

    /* sanity checks on arguments */

    if (!output || !input || !plan)
        XLAL_ERROR(XLAL_EFAULT);
    if (!plan->plan || !plan->size || plan->sign != 1)
        XLAL_ERROR(XLAL_EINVAL);
    if (!output->data || !input->data)
        XLAL_ERROR(XLAL_EINVAL);
    if (output->length != plan->size || input->length != plan->size / 2 + 1)
        XLAL_ERROR(XLAL_EBADLEN);
    if (CIMAGX(input->data[0]) != 0.0)
        XLAL_ERROR(XLAL_EDOM);  /* imaginary part of DC must be zero */
    if (! plan->size % 2 && CIMAGX(input->data[plan->size / 2]) != 0.0)
        XLAL_ERROR(XLAL_EDOM);  /* imaginary part of Nyquist must be zero */

    output_data = output->data;
    nbytes = plan->size * sizeof(REAL_TYPE);

    /* create temporary storage space; make sure that output data is
     * aligned, if memory alignment is required */

#   ifdef LAL_FFTW3_MEMALIGN_ENABLED
    tmp = XLALMallocAligned(nbytes);
    if (!tmp)
        XLAL_ERROR(XLAL_ENOMEM);
    if (!LAL_IS_MEMORY_ALIGNED(output_data)) {
        output_data = XLALMallocAligned(nbytes);
        if (!output_data) {
            XLALFreeAligned(tmp);
            XLAL_ERROR(XLAL_ENOMEM);
        }
    }
#   else
    tmp = XLALMalloc(nbytes);
    if (!tmp)
        XLAL_ERROR(XLAL_ENOMEM);
#   endif

    /* unpack input into temporary array */

    /* dc component */
    tmp[0] = CREALX(input->data[0]);

    /* other components */
    for (k = 1; k < (plan->size + 1) / 2; ++k) {        /* k < size/2 rounded up */
        tmp[k] = CREALX(input->data[k]);
        tmp[plan->size - k] = CIMAGX(input->data[k]);
    }

    /* Nyquist component */
    if (plan->size % 2 == 0)    /* n is even */
        tmp[plan->size / 2] = CREALX(input->data[plan->size / 2]);

    /* perform the fft */

    FFTWX_EXECUTE_R2R(plan->plan, tmp, output_data);

    /* if temporary space for output data was created, copy data into
     * the output vector and free the temporary space */

#   ifdef LAL_FFTW3_MEMALIGN_ENABLED
    if (output_data != output->data) {
        memcpy(output->data, output_data, nbytes);
        XLALFreeAligned(output_data);
    }
#   endif

    /* cleanup and return */

#   ifdef LAL_FFTW3_MEMALIGN_ENABLED
    XLALFreeAligned(tmp);
#   else
    XLALFree(tmp);
#   endif

    return 0;
}

int VECTOR_FFT_FUNCTION(REAL_VECTOR_TYPE * _LAL_RESTRICT_ output, const REAL_VECTOR_TYPE * _LAL_RESTRICT_ input,
    const PLAN_TYPE * plan)
{
    REAL_TYPE *input_data;
    REAL_TYPE *output_data;

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
        input_data = XLALMallocAligned(plan->size * sizeof(REAL_TYPE));
        if (!input_data)
            XLAL_ERROR(XLAL_ENOMEM);
        memcpy(input_data, input->data, plan->size * sizeof(REAL_TYPE));
    }
    if (!LAL_IS_MEMORY_ALIGNED(output_data)) {
        output_data = XLALMallocAligned(plan->size * sizeof(REAL_TYPE));
        if (!output_data) {
            if (input_data != input->data)
                XLALFreeAligned(input_data);
            XLAL_ERROR(XLAL_ENOMEM);
        }
    }
#   endif

    /* perform the fft */

    FFTWX_EXECUTE_R2R(plan->plan, input_data, output_data);

    /* cleanup aligned memory space if memory alignment is required;
     * copy data from temporary space to output vector */

#   ifdef LAL_FFTW3_MEMALIGN_ENABLED
    if (input_data != input->data)
        XLALFreeAligned(input_data);
    if (output_data != output->data) {
        memcpy(output->data, output_data, plan->size * sizeof(REAL_TYPE));
        XLALFreeAligned(output_data);
    }
#   endif

    return 0;
}

int POWER_SPECTRUM_FUNCTION(REAL_VECTOR_TYPE * spec, const REAL_VECTOR_TYPE * data, const PLAN_TYPE * plan)
{
    REAL_TYPE *tmp;
    REAL_TYPE *input_data;
    UINT4 k;
    size_t nbytes;

    /* sanity check on arguments */

    if (!spec || !data || !plan)
        XLAL_ERROR(XLAL_EFAULT);
    if (!plan->plan || !plan->size)
        XLAL_ERROR(XLAL_EINVAL);
    if (!spec->data || !data->data)
        XLAL_ERROR(XLAL_EINVAL);
    if (data->length != plan->size || spec->length != plan->size / 2 + 1)
        XLAL_ERROR(XLAL_EBADLEN);

    input_data = data->data;
    nbytes = plan->size * sizeof(REAL_TYPE);

    /* create temporary storage space; make sure that input data is
     * aligned, if memory alignment is required */

#   ifdef LAL_FFTW3_MEMALIGN_ENABLED
    tmp = XLALMallocAligned(nbytes);
    if (!tmp)
        XLAL_ERROR(XLAL_ENOMEM);
    if (!LAL_IS_MEMORY_ALIGNED(input_data)) {
        /* need to create temporary aligned space for input data */
        input_data = XLALMallocAligned(nbytes);
        if (!input_data) {
            XLALFreeAligned(tmp);
            XLAL_ERROR(XLAL_ENOMEM);
        }
        memcpy(input_data, data->data, nbytes);
    }
#   else
    tmp = XLALMalloc(nbytes);
    if (!tmp)
        XLAL_ERROR(XLAL_ENOMEM);
#   endif

    /* perform the fft */

    FFTWX_EXECUTE_R2R(plan->plan, input_data, tmp);

    /* compute spectrum from the fft of the data */

    /* dc component */
    spec->data[0] = tmp[0] * tmp[0];

    /* other components */
    for (k = 1; k < (plan->size + 1) / 2; ++k) {        /* k < size/2 rounded up */
        REAL_TYPE re = tmp[k];
        REAL_TYPE im = tmp[plan->size - k];
        spec->data[k] = re * re + im * im;
        spec->data[k] *= 2.0;   /* accounts for negative frequency part */
    }

    /* Nyquist frequency */
    if (plan->size % 2 == 0)    /* size is even */
        spec->data[plan->size / 2] = tmp[plan->size / 2] * tmp[plan->size / 2];

    /* cleanup and return */

#   ifdef LAL_FFTW3_MEMALIGN_ENABLED
    if (input_data != data->data)
        XLALFreeAligned(input_data);
    XLALFreeAligned(tmp);
#   else
    XLALFree(tmp);
#   endif

    return 0;
}

/*
 * Legacy Routines
 */

#define OLD_CREATE_FORWARD_PLAN_FUNCTION	CONCAT2(LALCreateForward,PLAN_TYPE)
#define OLD_CREATE_REVERSE_PLAN_FUNCTION	CONCAT2(LALCreateReverse,PLAN_TYPE)
#define OLD_DESTROY_PLAN_FUNCTION		CONCAT2(LALDestroy,PLAN_TYPE)
#define OLD_FORWARD_FFT_FUNCTION		CONCAT3(LALForward,REAL_TYPE,FFT)
#define OLD_REVERSE_FFT_FUNCTION		CONCAT3(LALReverse,REAL_TYPE,FFT)
#define OLD_VECTOR_FFT_FUNCTION			CONCAT3(LAL,REAL_VECTOR_TYPE,FFT)
#define OLD_POWER_SPECTRUM_FUNCTION		CONCAT3(LAL,REAL_TYPE,PowerSpectrum)

void OLD_CREATE_FORWARD_PLAN_FUNCTION(LALStatus * status, PLAN_TYPE ** plan, UINT4 size, INT4 measure)
{
    INITSTATUS(status);
    XLAL_PRINT_DEPRECATION_WARNING(STRING(CREATE_FORWARD_FUNCTION));

    ASSERT(plan, status, REALFFTH_ENULL, REALFFTH_MSGENULL);
    ASSERT(!*plan, status, REALFFTH_ENNUL, REALFFTH_MSGENNUL);
    ASSERT(size > 0, status, REALFFTH_ESIZE, REALFFTH_MSGESIZE);

    *plan = CREATE_PLAN_FUNCTION(size, 1, measure);
    if (!*plan) {
        int code = xlalErrno;
        XLALClearErrno();
        switch (code) {
        case XLAL_EBADLEN:
            ABORT(status, REALFFTH_ESIZE, REALFFTH_MSGESIZE);
        case XLAL_ENOMEM:
            ABORT(status, REALFFTH_EALOC, REALFFTH_MSGEALOC);
        case XLAL_EFAILED:
            ABORT(status, REALFFTH_EFFTW, REALFFTH_MSGEFFTW);
        default:
            ABORTXLAL(status);
        }
    }

    RETURN(status);
}

void OLD_CREATE_REVERSE_PLAN_FUNCTION(LALStatus * status, PLAN_TYPE ** plan, UINT4 size, INT4 measure)
{
    INITSTATUS(status);
    XLAL_PRINT_DEPRECATION_WARNING(STRING(CREATE_REVERSE_FUNCTION));

    ASSERT(plan, status, REALFFTH_ENULL, REALFFTH_MSGENULL);
    ASSERT(!*plan, status, REALFFTH_ENNUL, REALFFTH_MSGENNUL);
    ASSERT(size > 0, status, REALFFTH_ESIZE, REALFFTH_MSGESIZE);

    *plan = CREATE_PLAN_FUNCTION(size, 0, measure);
    if (!*plan) {
        int code = xlalErrno;
        XLALClearErrno();
        switch (code) {
        case XLAL_EBADLEN:
            ABORT(status, REALFFTH_ESIZE, REALFFTH_MSGESIZE);
        case XLAL_ENOMEM:
            ABORT(status, REALFFTH_EALOC, REALFFTH_MSGEALOC);
        case XLAL_EFAILED:
            ABORT(status, REALFFTH_EFFTW, REALFFTH_MSGEFFTW);
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
    ASSERT(plan, status, REALFFTH_ENULL, REALFFTH_MSGENULL);
    ASSERT(*plan, status, REALFFTH_ENULL, REALFFTH_MSGENULL);
    DESTROY_PLAN_FUNCTION(*plan);
    if (xlalErrno) {
        int code = xlalErrno;
        XLALClearErrno();
        switch (code) {
        case XLAL_EINVAL:
            ABORT(status, REALFFTH_ENULL, REALFFTH_MSGENULL);
        default:
            ABORTXLAL(status);
        }
    }
    *plan = NULL;
    RETURN(status);
}

void OLD_FORWARD_FFT_FUNCTION(LALStatus * status, COMPLEX_VECTOR_TYPE * output, REAL_VECTOR_TYPE * input, PLAN_TYPE * plan)
{
    int code;
    UINT4 n;
    INITSTATUS(status);
    XLAL_PRINT_DEPRECATION_WARNING(STRING(FORWARD_FFT_FUNCTION));

    ASSERT(output, status, REALFFTH_ENULL, REALFFTH_MSGENULL);
    ASSERT(input, status, REALFFTH_ENULL, REALFFTH_MSGENULL);
    ASSERT(plan, status, REALFFTH_ENULL, REALFFTH_MSGENULL);

    ASSERT(output->data, status, REALFFTH_ENULL, REALFFTH_MSGENULL);
    ASSERT(input->data, status, REALFFTH_ENULL, REALFFTH_MSGENULL);
    ASSERT(plan->plan, status, REALFFTH_ENULL, REALFFTH_MSGENULL);

    n = plan->size;
    ASSERT(n > 0, status, REALFFTH_ESIZE, REALFFTH_MSGESIZE);
    ASSERT(input->length == n, status, REALFFTH_ESZMM, REALFFTH_MSGESZMM);
    ASSERT(output->length == n / 2 + 1, status, REALFFTH_ESZMM, REALFFTH_MSGESZMM);

    ASSERT(plan->sign == -1, status, REALFFTH_ESIGN, REALFFTH_MSGESIGN);

    code = FORWARD_FFT_FUNCTION(output, input, plan);
    if (code) {
        code = xlalErrno;
        XLALClearErrno();
        switch (code) {
        case XLAL_ENOMEM:
            ABORT(status, REALFFTH_EALOC, REALFFTH_MSGEALOC);
        case XLAL_EINVAL:
            if (!n) {   /* plan size was invalid */
                ABORT(status, REALFFTH_ESIZE, REALFFTH_MSGESIZE);
            } else if (plan->sign != -1) {      /* plan sign was wrong */
                ABORT(status, REALFFTH_ESIGN, REALFFTH_MSGESIGN);
            } else {    /* one of the data pointers was NULL */

                ABORT(status, REALFFTH_ENULL, REALFFTH_MSGENULL);
            }
        case XLAL_EBADLEN:     /* size mismatch */
            ABORT(status, REALFFTH_ESZMM, REALFFTH_MSGESZMM);
        default:
            ABORTXLAL(status);
        }
    }

    RETURN(status);
}

void OLD_REVERSE_FFT_FUNCTION(LALStatus * status, REAL_VECTOR_TYPE * output, COMPLEX_VECTOR_TYPE * input, PLAN_TYPE * plan)
{
    int code;
    UINT4 n;
    INITSTATUS(status);
    XLAL_PRINT_DEPRECATION_WARNING(STRING(REVERSE_FFT_FUNCTION));

    ASSERT(output, status, REALFFTH_ENULL, REALFFTH_MSGENULL);
    ASSERT(input, status, REALFFTH_ENULL, REALFFTH_MSGENULL);
    ASSERT(plan, status, REALFFTH_ENULL, REALFFTH_MSGENULL);

    ASSERT(output->data, status, REALFFTH_ENULL, REALFFTH_MSGENULL);
    ASSERT(input->data, status, REALFFTH_ENULL, REALFFTH_MSGENULL);
    ASSERT(plan->plan, status, REALFFTH_ENULL, REALFFTH_MSGENULL);

    n = plan->size;
    ASSERT(n > 0, status, REALFFTH_ESIZE, REALFFTH_MSGESIZE);
    ASSERT(output->length == n, status, REALFFTH_ESZMM, REALFFTH_MSGESZMM);
    ASSERT(input->length == n / 2 + 1, status, REALFFTH_ESZMM, REALFFTH_MSGESZMM);
    ASSERT(cimagf(input->data[0]) == 0, status, REALFFTH_EDATA, REALFFTH_MSGEDATA);
    ASSERT(n % 2 || cimagf(input->data[n / 2]) == 0, status, REALFFTH_EDATA, REALFFTH_MSGEDATA);

    ASSERT(plan->sign == 1, status, REALFFTH_ESIGN, REALFFTH_MSGESIGN);

    code = REVERSE_FFT_FUNCTION(output, input, plan);
    if (code) {
        code = xlalErrno;
        XLALClearErrno();
        switch (code) {
        case XLAL_ENOMEM:
            ABORT(status, REALFFTH_EALOC, REALFFTH_MSGEALOC);
        case XLAL_EINVAL:
            if (!n) {   /* plan size was invalid */
                ABORT(status, REALFFTH_ESIZE, REALFFTH_MSGESIZE);
            } else if (plan->sign != 1) {       /* plan sign was wrong */
                ABORT(status, REALFFTH_ESIGN, REALFFTH_MSGESIGN);
            } else {    /* one of the data pointers was NULL */

                ABORT(status, REALFFTH_ENULL, REALFFTH_MSGENULL);
            }
        case XLAL_EBADLEN:     /* size mismatch */
            ABORT(status, REALFFTH_ESZMM, REALFFTH_MSGESZMM);
        case XLAL_EDOM:        /* either DC or Nyquist imaginary part was non zero */
            ABORT(status, REALFFTH_EDATA, REALFFTH_MSGEDATA);
        default:
            ABORTXLAL(status);
        }
    }

    RETURN(status);
}

void OLD_VECTOR_FFT_FUNCTION(LALStatus * status, REAL_VECTOR_TYPE * output, REAL_VECTOR_TYPE * input, PLAN_TYPE * plan)
{
    int code;
    INITSTATUS(status);
    XLAL_PRINT_DEPRECATION_WARNING(STRING(VECTOR_FFT_FUNCTION));

    ASSERT(output, status, REALFFTH_ENULL, REALFFTH_MSGENULL);
    ASSERT(input, status, REALFFTH_ENULL, REALFFTH_MSGENULL);
    ASSERT(plan, status, REALFFTH_ENULL, REALFFTH_MSGENULL);

    ASSERT(output->data, status, REALFFTH_ENULL, REALFFTH_MSGENULL);
    ASSERT(input->data, status, REALFFTH_ENULL, REALFFTH_MSGENULL);
    ASSERT(plan->plan, status, REALFFTH_ENULL, REALFFTH_MSGENULL);

    /* make sure that it is not the same data! */
    ASSERT(output->data != input->data, status, REALFFTH_ESAME, REALFFTH_MSGESAME);

    ASSERT(plan->size > 0, status, REALFFTH_ESIZE, REALFFTH_MSGESIZE);
    ASSERT(output->length == plan->size, status, REALFFTH_ESZMM, REALFFTH_MSGESZMM);
    ASSERT(input->length == plan->size, status, REALFFTH_ESZMM, REALFFTH_MSGESZMM);

    code = VECTOR_FFT_FUNCTION(output, input, plan);
    if (code) {
        code = xlalErrno;
        XLALClearErrno();
        switch (code) {
        case XLAL_EINVAL:
            if (!plan->size) {  /* plan size was invalid */
                ABORT(status, REALFFTH_ESIZE, REALFFTH_MSGESIZE);
            } else if (output->data == input->data) {   /* same data pointers */
                ABORT(status, REALFFTH_ESAME, REALFFTH_MSGESAME);
            } else {    /* one of the data pointers was NULL */

                ABORT(status, REALFFTH_ENULL, REALFFTH_MSGENULL);
            }
        case XLAL_EBADLEN:     /* size mismatch */
            ABORT(status, REALFFTH_ESZMM, REALFFTH_MSGESZMM);
        default:
            ABORTXLAL(status);
        }
    }

    RETURN(status);
}

void OLD_POWER_SPECTRUM_FUNCTION(LALStatus * status, REAL_VECTOR_TYPE * spec, REAL_VECTOR_TYPE * data, PLAN_TYPE * plan)
{
    int code;
    UINT4 n;

    INITSTATUS(status);
    XLAL_PRINT_DEPRECATION_WARNING(STRING(POWER_SPECTRUM_FUNCTION));

    ASSERT(spec, status, REALFFTH_ENULL, REALFFTH_MSGENULL);
    ASSERT(data, status, REALFFTH_ENULL, REALFFTH_MSGENULL);
    ASSERT(plan, status, REALFFTH_ENULL, REALFFTH_MSGENULL);

    ASSERT(spec->data, status, REALFFTH_ENNUL, REALFFTH_MSGENNUL);
    ASSERT(data->data, status, REALFFTH_ENNUL, REALFFTH_MSGENNUL);
    ASSERT(plan->plan, status, REALFFTH_ENNUL, REALFFTH_MSGENNUL);

    n = plan->size;
    ASSERT(n > 0, status, REALFFTH_ESIZE, REALFFTH_MSGESIZE);
    ASSERT(data->length == n, status, REALFFTH_ESZMM, REALFFTH_MSGESZMM);
    ASSERT(spec->length == n / 2 + 1, status, REALFFTH_ESZMM, REALFFTH_MSGESZMM);

    code = POWER_SPECTRUM_FUNCTION(spec, data, plan);
    if (code) {
        code = xlalErrno;
        XLALClearErrno();
        switch (code) {
        case XLAL_ENOMEM:
            ABORT(status, REALFFTH_EALOC, REALFFTH_MSGEALOC);
        case XLAL_EINVAL:
            if (!n) {   /* plan size was invalid */
                ABORT(status, REALFFTH_ESIZE, REALFFTH_MSGESIZE);
            } else {    /* one of the data pointers was NULL */

                ABORT(status, REALFFTH_ENULL, REALFFTH_MSGENULL);
            }
        case XLAL_EBADLEN:     /* size mismatch */
            ABORT(status, REALFFTH_ESZMM, REALFFTH_MSGESZMM);
        default:
            ABORTXLAL(status);
        }
    }

    RETURN(status);
}

#undef OLD_CREATE_FORWARD_PLAN_FUNCTION
#undef OLD_CREATE_REVERSE_PLAN_FUNCTION
#undef OLD_DESTROY_PLAN_FUNCTION
#undef OLD_FORWARD_FFT_FUNCTION
#undef OLD_REVERSE_FFT_FUNCTION
#undef OLD_VECTOR_FFT_FUNCTION
#undef OLD_POWER_SPECTRUM_FUNCTION

#undef CONCAT2x
#undef CONCAT2
#undef CONCAT3x
#undef CONCAT3
#undef STRING

#undef REAL_TYPE
#undef COMPLEX_TYPE
#undef TYPESUFFIX

#undef PLAN_TYPE
#undef REAL_VECTOR_TYPE
#undef COMPLEX_VECTOR_TYPE

#undef CREATE_PLAN_FUNCTION
#undef CREATE_FORWARD_PLAN_FUNCTION
#undef CREATE_REVERSE_PLAN_FUNCTION
#undef DESTROY_PLAN_FUNCTION
#undef FORWARD_FFT_FUNCTION
#undef REVERSE_FFT_FUNCTION
#undef VECTOR_FFT_FUNCTION
#undef POWER_SPECTRUM_FUNCTION

#undef CREALX
#undef CIMAGX
#undef FFTWX
#undef FFTWX_PLAN_R2R_1D
#undef FFTWX_DESTROY_PLAN
#undef FFTWX_EXECUTE_R2R
