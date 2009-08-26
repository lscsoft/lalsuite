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
    XLALSkymap2PlanType* plan,
    XLALSkymap2DirectionPropertiesType* properties,
    double* wSw,
    XLALSkymap2KernelType* kernel,
    double** xSw,
    int tau,
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
            double x[XLALSKYMAP2_N];

            // start accumulating inner products with the priors
            double q = (a[0] * a[0] + a[1] * a[1]);

            // for each detector...
            int j;
            for (j = 0; j != plan->n; ++j)
            {
                int k;
                // get the time-shifted data
                x[j] = XLALSkymap2Interpolate(xSw[j] + tau + properties->delay[j] - 1, properties->weight[j]);

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
    XLALSkymap2PlanType plan;    
    double direction[2];
    XLALSkymap2DirectionPropertiesType properties;
    double wSw[5] = { 100., 100., 100., 100., 100. };
    XLALSkymap2KernelType kernel;
    double *xSw[5];
    int siteNumbers[] = { LAL_LHO_4K_DETECTOR, LAL_LLO_4K_DETECTOR, LAL_VIRGO_DETECTOR, LAL_GEO_600_DETECTOR, LAL_LHO_2K_DETECTOR };
    RandomParams* rng;
    int n;

    rng = XLALCreateRandomParams(0);
    
    for (n = 1; n != 6; ++n)
    {
    
        XLALSkymap2PlanConstruct(8192, n, siteNumbers, &plan);

        direction[0] = LAL_PI * XLALUniformDeviate(rng);
        direction[1] = LAL_TWOPI * XLALUniformDeviate(rng);

        XLALSkymap2DirectionPropertiesConstruct(&plan, direction, &properties);

        XLALSkymap2KernelConstruct(&plan, &properties, wSw, &kernel);

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
            
            XLALSkymap2Apply(&plan, &properties, &kernel, xSw, plan.sampleFrequency / 2, &logPosteriorAnalytic);        

            numericApply(&plan, &properties, wSw, &kernel, xSw, plan.sampleFrequency / 2, & logPosteriorNumerical);

            if (abs(logPosteriorAnalytic - logPosteriorNumerical) > 1e-3)
            {
                // test failed
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

int main(int argc, char** argv)
{
    
    // check the fast analytic bayesian statistic against simpler but 
    // slower numerical integration
    
	numerical();

    return 0;
}

