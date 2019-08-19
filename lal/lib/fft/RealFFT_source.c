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
