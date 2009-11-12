#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

#include <lal/LALConstants.h>
#include <lal/Skymap.h>
#include <lal/Random.h>
#include <lal/Sort.h>

#include <fftw3.h>

static void numericApply(
    XLALSkymapPlanType* plan,
    XLALSkymapDirectionPropertiesType* properties,
    double* wSw,
    XLALSkymapKernelType* kernel,
    double** xSw,
    double tau,
    double* logPosterior
    )
{
    // this function has almost the same interface as XLALSkymap2Apply,
    // but is implemented with numerical integration that should converge
    // on the result from XLALSkymap2Apply

    // parameters
    double a[2];

    // the step over them
    double da = 0.01;

    // accumulate probability in this
    double p = 0.0;

    // kernel is unused
    kernel = NULL;

    // range from +/- 5 sigma in prior on a[i] for adequate accuracy
    for (a[0] = -5.0; a[0] <= 5.0; a[0] += da)
    {
        for (a[1] = - 5.0; a[1] <= 5.0; a[1] += da)
        {
            double x[XLALSKYMAP_N];

            // start accumulating inner products with the priors
            double q = (a[0] * a[0] + a[1] * a[1]);

            // for each detector...
            int j;
            for (j = 0; j != plan->n; ++j)
            {
                int k;
                // get the time-shifted data
                x[j] = XLALSkymapInterpolate((tau + properties->delay[j]) * plan->sampleFrequency, xSw[j]);

				// subtract the inner product
                q -= x[j] * x[j] / wSw[j];
                // for each polarization
                for (k = 0; k != 2; ++k)
                {
                    // subtract the postulated signal out
                    x[j] -= a[k] * properties->f[j][k] * wSw[j];
                }
                // add the inner product after signal removed
                q += x[j] * x[j] / wSw[j];
            }
            // accumulate the piece of the integral
            p += da * da * exp(-0.5*q);
        }
    }

    // normalization from the priors on a
    p *= pow(LAL_TWOPI, -1);

    // return the log
    *logPosterior = log(p);

}

static void numerical(void)
{
    XLALSkymapPlanType plan;
    double direction[2];
    XLALSkymapDirectionPropertiesType properties;
    double wSw[5] = { 100., 100., 100., 100., 100. };
    XLALSkymapKernelType kernel;
    double *xSw[5];
    int siteNumbers[] = { LAL_LHO_4K_DETECTOR, LAL_LLO_4K_DETECTOR, LAL_VIRGO_DETECTOR, LAL_GEO_600_DETECTOR, LAL_LHO_2K_DETECTOR };
    RandomParams* rng;
    int n;

    rng = XLALCreateRandomParams(0);

    for (n = 1; n != 6; ++n)
    {

        XLALSkymapPlanConstruct(8192, n, siteNumbers, &plan);

        direction[0] = LAL_PI * XLALUniformDeviate(rng);
        direction[1] = LAL_TWOPI * XLALUniformDeviate(rng);

        XLALSkymapDirectionPropertiesConstruct(&plan, direction, &properties);

        XLALSkymapKernelConstruct(&plan, &properties, wSw, &kernel);

        {
            int i;

            for (i = 0; i != n; ++i)
            {
                int j;
                xSw[i] = malloc(sizeof(*xSw[i]) * plan.sampleFrequency);
                for (j = 0; j != plan.sampleFrequency; ++j)
                {
                    xSw[i][j] = XLALNormalDeviate(rng) * sqrt(wSw[i]);
                }
            }
        }

        {
            double logPosteriorAnalytic;
            double logPosteriorNumerical;

            XLALSkymapApply(&plan, &properties, &kernel, xSw, 0.5, &logPosteriorAnalytic);

            numericApply(&plan, &properties, wSw, &kernel, xSw, 0.5, & logPosteriorNumerical);

            if (abs(logPosteriorAnalytic - logPosteriorNumerical) > 1e-3)
            {
                fprintf(stderr, "Analytic expression does not match numerical result to expected accuracy\n");
                exit(1);
            }

        }

        {
            int i;
            for(i = 0; i != n; ++i)
                free(xSw[i]);
        }

    }

}

static void interpolation(void)
{
    double* x;
    int n;

    n = 1000;

    x = (double*) malloc(sizeof(double) * n);

    for (int i = 0; i != n; ++i)
    {
        x[i] = sin(i);
    }

    for (double t = n * 0.25; t < n * 0.75; t += 0.1)
    {
        if (abs(XLALSkymapInterpolate(t, x) - sin(t)) > 0.03)
        {
            fprintf(stderr, "Interpolation error larger than expected\n");
            exit(1);
        }
    }

    for (int t = n / 4; t < (n * 3) / 4; ++t)
    {
        if (abs(XLALSkymapInterpolate(t, x) - sin(t)) != 0)
        {
            fprintf(stderr, "Interpolation does not pass through data points\n");
            exit(1);
        }
    }
}

int main(void)
{

    // check the fast analytic bayesian statistic against simpler but
    // slower numerical integration

	numerical();

    interpolation();

    return 0;
}

