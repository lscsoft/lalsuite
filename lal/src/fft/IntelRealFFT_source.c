#define CONCAT2x(a,b) a##b
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a##b##c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define STRING(a) #a

#ifdef SINGLE_PRECISION
#define REAL_TYPE REAL4
#define DFTI_TYPE DFTI_SINGLE
#define COMPLEX_TYPE COMPLEX8
#define TYPESUFFIX f
#else
#define REAL_TYPE REAL8
#define DFTI_TYPE DFTI_DOUBLE
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

PLAN_TYPE *CREATE_PLAN_FUNCTION(UINT4 size, int fwdflg, __attribute__((unused)) int measurelvl)
{
    PLAN_TYPE *plan;
    INT8  fftStat;

    if (!size)
        XLAL_ERROR_NULL(XLAL_EBADLEN);

    /* allocate memory for the plan */
    plan = XLALMalloc(sizeof( *plan));

    if (!plan)
        XLAL_ERROR_NULL(XLAL_ENOMEM);

    /* make the intel fft descriptor */
    fftStat = DftiCreateDescriptor(&(plan->plan),
        DFTI_TYPE, DFTI_REAL, 1, size);
    CHECKINTELFFTSTATUS_NULL(fftStat);

    /* configure intel fft descriptor */
    fftStat = DftiSetValue(plan->plan, DFTI_PLACEMENT,
        DFTI_NOT_INPLACE);
    CHECKINTELFFTSTATUS_NULL(fftStat);
    fftStat = DftiSetValue(plan->plan, DFTI_PACKED_FORMAT,
        DFTI_PACK_FORMAT);
    CHECKINTELFFTSTATUS_NULL(fftStat);

    /* commit the intel fft descriptor */
    fftStat = DftiCommitDescriptor(plan->plan);
    CHECKINTELFFTSTATUS_NULL(fftStat);

    /* create workspace to do the fft into */
    plan->tmp = XLALMalloc(size * sizeof(REAL_TYPE));
    if (!plan->tmp)
    {
      XLALFree(plan);
      XLAL_ERROR_NULL(XLAL_ENOMEM);
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
    INT8  fftStat;

    if (!plan)
        XLAL_ERROR_VOID( XLAL_EFAULT );

    /* destroy intel fft descriptor */
    fftStat = DftiFreeDescriptor(&(plan->plan));
    CHECKINTELFFTSTATUS_VOID(fftStat);

    XLALFree(plan->tmp);
    XLALFree(plan);

    return;
}


int FORWARD_FFT_FUNCTION(COMPLEX_VECTOR_TYPE * output, const REAL_VECTOR_TYPE * input, const PLAN_TYPE * plan)
{
    INT8  fftStat;
    UINT4 n,k;

    if (!output || !input || !plan)
        XLAL_ERROR(XLAL_EFAULT);
    if (!plan->plan || !plan->size || plan->sign != -1)
        XLAL_ERROR(XLAL_EINVAL);
    if (!output->data || !input->data)
        XLAL_ERROR(XLAL_EINVAL);
    if (input->length != plan->size || output->length != plan->size/2 + 1)
        XLAL_ERROR(XLAL_EBADLEN);

    n = plan->size;

    /* execute intel fft */
    fftStat = DftiComputeForward(plan->plan, input->data, plan->tmp);
    CHECKINTELFFTSTATUS(fftStat);

    /* now unpack the results into the output vector */

    /* dc component */
    output->data[0] = plan->tmp[0] + 0.0 * _Complex_I;

    /* other components */
    for (k = 1; k < (n + 1) / 2; ++k) /* k < n/2 rounded up */
    {
        output->data[k] = plan->tmp[2 * k - 1] + (plan->tmp[2 * k] * _Complex_I);
    }

    /* Nyquist frequency */
    if (n % 2 == 0) /* n is even */
    {
        output->data[n / 2] = plan->tmp[n - 1] + 0.0 * _Complex_I;
    }

    return 0;
}


int REVERSE_FFT_FUNCTION(REAL_VECTOR_TYPE * output, const COMPLEX_VECTOR_TYPE * input, const PLAN_TYPE * plan)
{
    INT8  fftStat;
    UINT4 n,k;

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

    n = plan->size;

    /* dc component */
    plan->tmp[0] = CREALX(input->data[0]);

    /* other components */
    for ( k = 1; k < ( n + 1 ) / 2; ++k ) /* k < n / 2 rounded up */
    {
        plan->tmp[2 * k - 1] = CREALX(input->data[k]);
        plan->tmp[2 * k]     = CIMAGX(input->data[k]);
   }

    /* Nyquist component */
    if ( n % 2 == 0 ) /* n is even */
    {
        plan->tmp[n - 1] = CREALX(input->data[n / 2]);
    }

    /* execute intel fft */
    fftStat = DftiComputeBackward( plan->plan, plan->tmp, output->data );
    CHECKINTELFFTSTATUS( fftStat );

    return 0;
}


int VECTOR_FFT_FUNCTION(REAL_VECTOR_TYPE * _LAL_RESTRICT_ output, const REAL_VECTOR_TYPE * _LAL_RESTRICT_ input,
    const PLAN_TYPE * plan)
{
    INT8  fftStat;
    UINT4 n,k;

    /* sanity check on arguments */
    if (!output || !input || !plan)
        XLAL_ERROR(XLAL_EFAULT);
    if (!plan->plan || !plan->size)
        XLAL_ERROR(XLAL_EINVAL);
    if (!output->data || !input->data || output->data == input->data)
        XLAL_ERROR(XLAL_EINVAL);        /* note: must be out-of-place */
    if (output->length != plan->size || input->length != plan->size)
        XLAL_ERROR(XLAL_EBADLEN);

    n = plan->size;

    /* complex intel fft */
    if (plan->sign == -1)
    {
        fftStat = DftiComputeForward(plan->plan, input->data, plan->tmp);
        CHECKINTELFFTSTATUS(fftStat);

        output->data[0] = plan->tmp[0];
        for (k = 1 ; k < (n + 1) / 2 ; ++k) /* k < n/2 rounded up */
        {
            output->data[k]     = plan->tmp[2 * k - 1];
            output->data[n - k] = plan->tmp[2 * k];
        }
        if ( n % 2 == 0) /* n is even */
        {
            output->data[n / 2] = plan->tmp[n - 1];
        }
    }
    else if (plan->sign == 1)
    {
        plan->tmp[0] = input->data[0];
        for ( k = 1 ; k < (n + 1) / 2 ; ++k ) /* k < n/2 rounded up */
        {
            plan->tmp[2 * k - 1] = input->data[k];
            plan->tmp[2 * k]     = input->data[n - k];
        }
        if ( n % 2 == 0) /* n is even */
        {
            plan->tmp[n - 1] = input->data[n / 2];
        }

        fftStat = DftiComputeBackward( plan->plan, plan->tmp, output->data );
        CHECKINTELFFTSTATUS( fftStat );
    }
    else
    {
        XLAL_ERROR( XLAL_EINVAL );
    }

    return 0;
}

int POWER_SPECTRUM_FUNCTION(REAL_VECTOR_TYPE * spec, const REAL_VECTOR_TYPE * data, const PLAN_TYPE * plan)
{
    INT8  fftStat;
    UINT4 n,k;

    /* sanity check on arguments */
    if (!spec || !data || !plan)
        XLAL_ERROR(XLAL_EFAULT);
    if (!plan->plan || !plan->size)
        XLAL_ERROR(XLAL_EINVAL);
    if (!spec->data || !data->data)
        XLAL_ERROR(XLAL_EINVAL);
    if (data->length != plan->size || spec->length != plan->size / 2 + 1)
        XLAL_ERROR(XLAL_EBADLEN);

    n = plan->size;

    /* execute intel fft */
    fftStat = DftiComputeForward( plan->plan, data->data, plan->tmp );
    CHECKINTELFFTSTATUS( fftStat );

    /* dc component */
    spec->data[0] = plan->tmp[0] * plan->tmp[0];

    /* other components */
    for (k = 1; k < (n + 1)/2; ++k) /* k < n/2 rounded up */
    {
        REAL_TYPE re = plan->tmp[2 * k - 1];
        REAL_TYPE im = plan->tmp[2 * k];
        spec->data[k] = re * re + im * im;
    }

    /* Nyquist frequency */
    if (n % 2 == 0) /* n is even */
      spec->data[n / 2] = plan->tmp[n / 2] * plan->tmp[n / 2];

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
#undef DFTI_TYPE

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
