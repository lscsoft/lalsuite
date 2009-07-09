#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

#include <lal/LALConstants.h>
#include <lal/Skymap.h>
#include <lal/Random.h>
#include <lal/Sort.h>

#include <fftw3.h>

#define TEST(A) if (!(A)) { printf("FAIL: %s\n", #A); exit(1); }

#if 0
static void minimal(void)
{

    XLALSkymap2PlanType *plan;
    XLALSkymap2SphericalPolarType *direction;
    XLALSkymap2DirectionPropertiesType *properties;
    double wSw[3];
    XLALSkymap2KernelType *kernel;
    double *xSw[3];

    {
        plan = malloc(sizeof(*plan));
        XLALSkymap2PlanConstruct(8192, plan);

        TEST(plan->sampleFrequency == 8192);
    }

    {
        direction = malloc(sizeof(*direction));
        (*direction)[0] = LAL_PI_2;
        (*direction)[1] = 0.;
    }

    {
        properties = malloc(sizeof(*properties));
        XLALSkymap2DirectionPropertiesConstruct(plan, direction, properties);
    }

    {
        int i;
        for (i = 0; i != 3; ++i)
            wSw[i] = 1.;
    }

    {
        kernel = malloc(sizeof(*kernel));
        XLALSkymap2KernelConstruct(properties, wSw, kernel);
    }

    {
        int i;
        for (i = 0; i != 3; ++i)
        {
            int j;
            xSw[i] = malloc(sizeof(*xSw[i]) * plan->sampleFrequency);
            for (j = 0; j != plan->sampleFrequency; ++j)
            {
                xSw[i][j] = 0.;
            }
        }
    }

    {
        double logPosterior;
        XLALSkymap2Apply(properties, kernel, xSw, plan->sampleFrequency / 2, &logPosterior);

        printf("%g\n", exp(logPosterior));
    }

    {
        int i;
        for(i = 0; i != 3; ++i)
            free(xSw[i]);
    }

    free(kernel);
    free(properties);
    free(direction);
    free(plan);

}
#endif

static void directional(void)
{
    XLALSkymap2PlanType *plan;
    XLALSkymap2SphericalPolarType *directions;
    XLALSkymap2DirectionPropertiesType *properties;
    double wSw[3] = { 100., 100., 10. };
    XLALSkymap2KernelType *kernels;
    double *xSw[3];

    {
        plan = malloc(sizeof(*plan));
        XLALSkymap2PlanConstruct(512, plan);

        TEST(plan->sampleFrequency == 512);
    }

    {
        int i;
        directions = malloc(sizeof(*directions) * 180 * 360);
        for (i = 0; i != 180; ++i)
        {
            int j;
            for (j = 0; j != 360; ++j)
            {
                directions[i * 360 + j][0] = (i + 0.5) / 180. * LAL_PI;
                directions[i * 360 + j][1] = (j + 0.5) / 180. * LAL_PI;
            }
        }
    }

    {
        int i;
        properties = malloc(sizeof(*properties) * 180 * 360);
        for (i = 0; i != 180 * 360; ++i)
        {
            XLALSkymap2DirectionPropertiesConstruct(
                plan,
                directions + i,
                properties + i
                );
        }
    }

    {
        int i;
        kernels = malloc(sizeof(*kernels) * 180 * 360);
        for (i = 0; i != 180 * 360; ++i)
        {
            XLALSkymap2KernelConstruct(properties + i, wSw, kernels + i);
        }
    }

    {
        int i;
        RandomParams* rng;
        rng = XLALCreateRandomParams(0);
        for (i = 0; i != 3; ++i)
        {
            int j;
            xSw[i] = malloc(sizeof(*xSw[i]) * plan->sampleFrequency);
            for (j = 0; j != plan->sampleFrequency; ++j)
            {
                xSw[i][j] = XLALNormalDeviate(rng) * sqrt(wSw[i]);
            }
            //xSw[i][plan->sampleFrequency/2] += XLALNormalDeviate(rng) * wSw[i];
        }
    }

    {
        int i;
        for (i = 0; i != 180 * 360; ++i)
        {
            int t;
            double p = 0;
            for (t = plan->sampleFrequency / 4; t != plan->sampleFrequency / 4 * 3; ++t)
            {
                double logPosterior;
                XLALSkymap2Apply(properties + i, kernels + i, xSw, t, &logPosterior);
                p += exp(logPosterior) / (plan->sampleFrequency / 2);

            }
            printf("%g %g %g\n", directions[i][0], directions[i][1], log(p));
        }
    }

    {
        int i;
        for(i = 0; i != 3; ++i)
            free(xSw[i]);
    }

    free(kernels);
    free(properties);
    free(directions);
    free(plan);

}

int main(int argc, char** argv)
{
    // minimal();
    directional();

    return 0;
}

#if 0

int main(int argc, char** argv)
{

    int sampleFrequency = 8192;
    int m = 1024;
    int n = 2048;

    XLALSkymapPlanType* plan = 0;

    double* p = 0;
    double sigma = 1;
    double w[3] = {1, 1, 1};
    int begin[3] = { 4000, 4000, 4000 };
    int end[3] = {4200, 4200, 4200 };
    double* x[6] = { 0, 0, 0, 0, 0, 0 };
    int* counts = 0;
    int* modes = 0;
    int delay_limits[6] = { -8192, 8192, -8192, 8192, -8192, 8192 };

    RandomParams* rng_parameters;

    FILE* h;

    double* q = 0;

    int i;
    int j;

    lalDebugLevel = 0xFFFFFFFF;

    if (!(plan = XLALSkymapConstructPlanMN(sampleFrequency, m, n)))
    {
        printf("XLALSkymapConstructPlanMN failed\n");
        exit(1);
    }

    rng_parameters = XLALCreateRandomParams(0);

    for (i = 0; i != 6; ++i)
    {
        x[i] = XLALMalloc(sizeof(double) * sampleFrequency);
        for (j = 0; j != sampleFrequency; ++j)
        {
            x[i][j] = XLALNormalDeviate(rng_parameters);
        }
    }

    p = (double*) XLALMalloc(sizeof(double) * plan->pixelCount);
    counts = (int*) XLALMalloc(sizeof(int) * plan->pixelCount);
    modes = (int*) XLALMalloc(sizeof(int) * plan->pixelCount);

    if (XLALSkymapSignalHypothesisWithLimits(plan, p, sigma, w, begin, end, x, counts, modes, delay_limits))
    {
        printf("XLALSkymapSignalHypothesisWithLimits failed\n");
        exit(1);
    }

    h = fopen("p.bin", "w");
    fwrite(p, sizeof(double), plan->pixelCount, h);
    fclose(h);

    q = XLALMalloc(sizeof(double) * m * n);

    if (XLALSkymapRender(q, plan, p))
    {
        printf("XLALSkymapRender failed\n");
        exit(1);
    }

    h = fopen("q.bin", "w");
    fwrite(q, sizeof(double), m * n, h);
    fclose(h);

    return 0;
}

#endif

#if 0
#include <stdlib.h>
#include <math.h>

#include <fftw3.h>

typedef struct XLALSkymapFollowupInputStruct
{
    fftw_complex *data[3]; /* one-sided frequency domain double-whitened data */
    double *psd[3]; /* one-sided power spectrum */

    int duration; /* seconds */
    int inputRate; /* power-of-two integer */
    int analysisRate; /* power-of-two integer */

    double minFrequency;
    double maxFrequency;
    double minTime;
    double maxTime;

} XLALSkymapFollowupInput;

typedef struct XLALSkymapFollowupOutputStruct
{
} XLALSkymapFollowupOutput;

/* declare a C99 function */
double log2(double);

XLALSkymapFollowupOutput* XLALSkymapFollowup(XLALSkymapFollowupInput* in)
{
    /* validate arguments */

    /* now perform a time-freqiency-Q decomposition of the region of interest */

    printf("%f volume\n", (in->maxTime - in->minTime) * (in->maxFrequency - in->minFrequency));

    {
        int width;

        int minWidth = pow(2, ceil(log2(in->duration / (in->maxTime - in->minTime))));
        int maxWidth = pow(2, floor(log2(in->duration * (in->maxFrequency - in->minFrequency))));

        for (width = minWidth; width <= maxWidth; width *= 2)
        {
            int f;

            int minFrequency = ceil(4 * in->minFrequency * in->duration / width) * width / 4;
            int maxFrequency = floor(4 * in->maxFrequency * in->duration / width) * width / 4 - width;

            printf("%f Hz (%f...%f)\n", ((double) width) / in->duration, ((double) minFrequency) / in->duration, ((double) maxFrequency + width) / in->duration);

            for (f = minFrequency; f <= maxFrequency; f += (width / 4))
            {
                printf("    [%f, %f)\n", ((double) f) / in->duration, ((double) f + width) / in->duration);
            }
        }

    }

    return 0;
}

int main(int argc, char** argv)
{
    XLALSkymapFollowupInput* in;

    in = (XLALSkymapFollowupInput *) malloc(sizeof(XLALSkymapFollowupInput));

    in->duration  =   64;
    in->inputRate = 4096;
    in->analysisRate = 8192;

    in->minFrequency = 100;
    in->maxFrequency = 300;
    in->minTime = 31.9;
    in->maxTime = 32.1;

    XLALSkymapFollowup(in);

    return 0;
}

#endif

#if 0

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

#include <lal/LALConstants.h>
#include <lal/Skymap.h>
#include <lal/Random.h>
#include <lal/Sort.h>

#include <fftw3.h>

#define max(A,B) (((A) > (B)) ? (A) : (B))
#define min(A,B) (((A) < (B)) ? (A) : (B))

RandomParams* rng_parameters;

int duration;
int rate;

int samples;

double *h_t[3]; /* h(t) */

fftw_complex *h_f[3]; /* h(f) */

int f_min; /* Hz */
int f_max; /* Hz */
int *widths;

XLALSkymapPlanType* plan;

static double sq(double a)
{
    return a * a;
}

void make_h_t()
{
    int i, j;
    for (i = 0; i != 3; ++i)
    {
        h_t[i] = (double*) malloc(sizeof(double) * samples);
        for (j = 0; j != samples; ++j)
        {
            h_t[i][j] = XLALNormalDeviate(rng_parameters);
        }
    }
}

void free_h_t()
{
    int i;
    for (i = 0; i !=3; ++i)
    {
        free(h_t[i]);
    }
}

double thetaphi[2];

void add_injection(XLALSkymapPlanType* plan, double sigma)
{
    int i, j;
    double h[4];
    double direction[3];
    int delays[3];

    thetaphi[0] = acos(XLALUniformDeviate(rng_parameters) * 2 - 1);
    thetaphi[1] = XLALUniformDeviate(rng_parameters) * LAL_TWOPI;
    XLALSkymapCartesianFromSpherical(direction, thetaphi[0], thetaphi[1]);

    XLALSkymapDelaysFromDirection(plan, delays, direction);
    i = XLALSkymapIndexFromDirection(plan, direction);

    for (j = samples / 2 - rate / 128; j < samples /2 + rate / 128; ++j)
    {
        h[0] = sigma * XLALNormalDeviate(rng_parameters);
        h[1] = sigma * XLALNormalDeviate(rng_parameters);

        h_t[0][j            ] += plan->pixel[i].f[0][0] * h[0] + plan->pixel[i].f[0][1] * h[1];
        h_t[1][j + delays[0]] += plan->pixel[i].f[1][0] * h[0] + plan->pixel[i].f[1][1] * h[1];
        h_t[2][j + delays[1]] += plan->pixel[i].f[2][0] * h[0] + plan->pixel[i].f[2][1] * h[1];
    }
}

void add_glitch(double sigma)
{
    int j;
    int k[3];
    k[0] = (rand() % (samples/2)) - samples / 4;
    k[1] = (rand() % (samples/2)) - samples / 4;
    k[2] = (rand() % (samples/2)) - samples / 4;

    for (j = samples / 2 - rate / 256; j < samples /2 + rate / 256; ++j)
    {
        h_t[0][j] += sigma * XLALNormalDeviate(rng_parameters);
        h_t[1][j] += sigma * XLALNormalDeviate(rng_parameters);
        h_t[2][j] += sigma * XLALNormalDeviate(rng_parameters);
    }

}

void make_h_f()
{
    int i, j;

    fftw_complex *dft_in;
    fftw_complex *dft_out;
    fftw_plan dft_forward;

    dft_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * samples);
    dft_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * samples);

    dft_forward = fftw_plan_dft_1d(samples, dft_in, dft_out, FFTW_FORWARD, FFTW_ESTIMATE);

    for (i = 0; i != 3; ++i)
    {
        for (j = 0; j != samples; ++j)
        {
            dft_in[j][0] = h_t[i][j] / sqrt(samples);
            dft_in[j][1] = 0.;
        }
        fftw_execute(dft_forward);
        h_f[i] = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * samples);
        for (j = 0; j != samples; ++j)
        {
            h_f[i][j][0] = dft_out[j][0];
            h_f[i][j][1] = dft_out[j][1];
        }
    }

    fftw_destroy_plan(dft_forward);
    fftw_free(dft_out);
    fftw_free(dft_in);

}

void free_h_f()
{
    int i;
    for (i = 0; i != 3; ++i)
    {
        fftw_free(h_f[i]);
    }
}

int less_double(void *p, const void *a, const void *b)
{
    if (*((double*) a) < *((double*) b))
    {
        return -1;
    }
    else
    {
        if (*((double*) a) > *((double*) b))
        {
            return 1;
        }
        else
        {
            return 0;
        }
    }
}

int main(int argc, char **argv)
{
    FILE *h, *g;
    int width;

    double *probabilities, *scores;
    int *widths, *frequencies, *times[3], *used;
    int n_bands;
    int i_band;
    int* index;

    double* accumulated_skymap;
    double accumulated_glitch, accumulated_signal;
    double w[3];

    w[0] = 1.0; w[1] = w[0]; w[2] = w[0];

    rng_parameters = XLALCreateRandomParams(0);

    duration = 1;
    rate = 8192;

    samples = duration * rate;

    f_min = 0;
    f_max = 1024;

    printf("making skymap...\n");
    plan = XLALSkymapConstructPlan(rate);
    printf("    ...done\n");

    printf("making data...\n");
    make_h_t();
    printf("    ...done\n");

    add_injection(plan, 10);
    /*add_glitch(10);*/

    printf("transforming data...\n");
    make_h_f();
    printf("    ...done\n");

    g = fopen("bands.dat", "wb");


    n_bands = (1024 / 16) * 4 * 2;
    i_band = 0;
    probabilities = (double*) malloc(sizeof(double) * n_bands);
    scores        = (double*) malloc(sizeof(double) * n_bands);
    widths        = (int   *) malloc(sizeof(int   ) * n_bands);
    frequencies   = (int   *) malloc(sizeof(int   ) * n_bands);
    times[0]      = (int   *) malloc(sizeof(int   ) * n_bands);
    times[1]      = (int   *) malloc(sizeof(int   ) * n_bands);
    times[2]      = (int   *) malloc(sizeof(int   ) * n_bands);
    used          = (int   *) malloc(sizeof(int   ) * n_bands);

    /* loop through window bandwidths */
    for (width = 16; width <= 1024; width *= 2)
    {
        /* set up the transforms */
        int f, i, n;
        fftw_complex *dft_in;
        fftw_complex *dft_out;
        fftw_plan dft_backward;
        double *buffer;
        int n_sigmas = 3;
        double sigmas[] = { sqrt(10), 10, sqrt(1000) };
        /* sigmas[0] = 0.01; sigmas[1] = sigmas[0]; sigmas[2] = sigmas[0]; */

        n = width * duration * 4;

        dft_in  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
        dft_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);

        buffer = (double*) malloc(sizeof(double) * n);

        dft_backward = fftw_plan_dft_1d(n, dft_in, dft_out, FFTW_BACKWARD, FFTW_ESTIMATE);

        for (i = 0; i != n; ++i)
        {
            dft_in[i][0] = 0.0;
            dft_in[i][1] = 0.0;
        }
        /* loop through bands */
        for (f = 0; f + width <= 1024; f += (width / 4))
        {
            /* loop through detectors */
            int d;
            double p_detectors;
            double p_sigmas[3] = { log(0), log(0), log(0) };

            printf("    [%d, %d) Hz", f, f + width);
            p_detectors = 0.0;
            for (d = 0; d != 3; ++d)
            {
                double s;
                int j;

                double p_sigma;
                s = 1.0 / sqrt(width * duration);
                for (i = 0; i != width * duration; ++i)
                {
                    dft_in[i][0] = h_f[d][f * duration + i][0] * s;
                    dft_in[i][1] = h_f[d][f * duration + i][1] * s;
                }
                fftw_execute(dft_backward);
                fwrite(dft_out, sizeof(fftw_complex), n, g);


                {
                    /* find the peak power, which will be where the peak
                     * probability will be for any sigma */
                    int best_i = 0;
                    double best_power = 0;
                    for (i = 0; i != n; ++i)
                    {
                        double power = sq(dft_out[i][0]) + sq(dft_out[i][1]);
                        if (power > best_power)
                        {
                            best_power = power;
                            best_i = i;
                        }
                    }
                    times[d][i_band] = best_i * samples / n;
                }

                /* compute the marginalized likelihood ratio */
                for (j = 0; j != n_sigmas; ++j)
                {
                    double kernel;


                    kernel = 0.5 / (sq(w[d]) + 1.0 / sq(sigmas[j]));

                    for (i = 0; i != n; ++i)
                    {
                        buffer[i] = kernel * (sq(dft_out[i][0]) + sq(dft_out[i][1]));
                    }
                    p_sigma = XLALSkymapLogTotalExp(buffer, buffer + n) - log(sq(sigmas[j]) * sq(w[d]) + 1.0) - log(n);
                    p_sigmas[d] = XLALSkymapLogSumExp(p_sigmas[d], p_sigma);
                }
                /* the log odds that a waveform was added to a detector */
                p_sigmas[d] -= log(n_sigmas);
                /* compute the log odds that a waveform was added to all detectors */
                p_detectors += XLALSkymapLogSumExp(p_sigmas, 0);
            }

            p_detectors = XLALSkymapLogDifferenceExp(p_detectors, 0);
            p_detectors -= log(7);

            printf(" %f %f %f %f\n", p_detectors, ((double) times[0][i_band])/rate, ((double) times[1][i_band])/rate, ((double) times[2][i_band])/rate);

            probabilities[i_band] = p_detectors;
            widths[i_band] = width;
            frequencies[i_band] = f;
            used[i_band] = 0;
            scores[i_band] = p_detectors / width;

            if (p_sigmas[1] > max(p_sigmas[0], p_sigmas[2]))
            {
                times[0][i_band] = times[1][i_band];
            }
            if (p_sigmas[2] > max(p_sigmas[0], p_sigmas[1]))
            {
                times[0][i_band] = times[2][i_band];
            }

            ++i_band;
        }

        fftw_destroy_plan(dft_backward);
        free(buffer);
        fftw_free(dft_out);
        fftw_free(dft_in);

    }

    fclose(g);

    printf("\n    --------    \n\n");

    n_bands = i_band;

    g = fopen("used.dat", "wb");

    /* sort the basis waveforms by peak power */
    index = (INT4*) malloc(sizeof(INT4) * n_bands);
    XLALHeapIndex(index, scores, n_bands, sizeof(double), 0, less_double);

    /* work backwards from most powerful tile */
    accumulated_skymap = (double*) malloc(plan->pixelCount * sizeof(double));
    {
        int i;
        for (i = 0; i != plan->pixelCount; ++i)
        {
            accumulated_skymap[i] = 0.0;
        }
    }
    accumulated_glitch = 0;
    h = fopen("raw.dat", "wb");
    fwrite(thetaphi, sizeof(double), 2, h);
    for (i_band = n_bands - 1; i_band >= 0; --i_band)
    {
        double *z_t[6];
        int f, i, d, time;
        fftw_complex *dft_in;
        fftw_complex *dft_out;
        fftw_plan dft_backward;
        int useable;

        width = widths     [index[i_band]];
        f     = frequencies[index[i_band]];
        time  = times[0]   [index[i_band]];

        /* printf("with %d frequency %d time %d\n", width, f, time); */

        if (probabilities[index[i_band]] > log(1))
        {
            used[i_band] = 1;
            for (i = i_band + 1; i < n_bands; ++i)
            {
                if (used[i])
                {
                    if (!((frequencies[index[i]] >= f + width) || (frequencies[index[i]] + widths[index[i]] <= f)))
                    {
                        used[i_band] = 0;
                    }
                }
            }
        }

        if (used[i_band])
        {

            printf("[%d, %d) Hz time %f probability %f\n", f, f+width, ((double)time)/rate, probabilities[index[i_band]]);

            dft_in  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * samples);
            dft_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * samples);

            dft_backward = fftw_plan_dft_1d(samples, dft_in, dft_out, FFTW_BACKWARD, FFTW_ESTIMATE);

            for (i = 0; i != samples; ++i)
            {
                dft_in[i][0] = 0.0;
                dft_in[i][1] = 0.0;
            }

            for (d = 0; d != 3; ++d)
            {
                double power;
                double s;
                s = sqrt(2 * rate / width) / sqrt(samples);
                for (i = 0; i != width * duration; ++i)
                {
                    dft_in[f * duration + i][0] = h_f[d][f * duration + i][0] * s;
                    dft_in[f * duration + i][1] = h_f[d][f * duration + i][1] * s;
                }
                fftw_execute(dft_backward);
                z_t[d    ] = (double*) malloc(sizeof(double) * samples);
                z_t[d + 3] = (double*) malloc(sizeof(double) * samples);
                power = 0;
                for (i = 0; i != samples; ++i)
                {
                    z_t[d    ][i] = dft_out[i][0];
                    z_t[d + 3][i] = dft_out[i][1];
                    power += pow(dft_out[i][0], 2) + pow(dft_out[i][1], 2);
                }
                fwrite(z_t[d], sizeof(double), samples, g);
                fwrite(z_t[d + 3], sizeof(double), samples, g);
                power /= samples;
                /* printf("mean(power) = %f =?= 2\n", power); */
            }

            /* we now have band-limited time series */

            {

                int begin[3];
                int end[3];
                double* skymap[3];
                int n_sigmas = 3;
                double sigmas[] = { sqrt(10), 10, sqrt(1000) };
                double g[3][3];
                double psignal;
                double pglitch;
                int k;
                int *modes;
                int *counts;
                double c;

                modes = (int *) malloc(sizeof(int) * plan->pixelCount);
                counts = (int *) malloc(sizeof(int) * plan->pixelCount);

                begin[0] = max(times[0][i_band] - (rate / 32), 0);
                begin[1] = max(times[0][i_band] - (rate / 32), 0);
                begin[2] = max(times[0][i_band] - (rate / 32), 0);
                end[0] = min(times[0][i_band] + (rate / 32), samples);
                end[1] = min(times[0][i_band] + (rate / 32), samples);
                end[2] = min(times[0][i_band] + (rate / 32), samples);

                for (k = 0; k != n_sigmas; ++k)
                {
                    skymap[k] = (double*) malloc(plan->pixelCount * sizeof(double));
                    XLALSkymapSignalHypothesis(plan, skymap[k], sigmas[k], w, begin, end, z_t, counts, modes);
                }
                for (k = 1; k != n_sigmas; ++k)
                {
                    XLALSkymapSum(plan, skymap[0], skymap[0], skymap[k]);
                }
                for (i = 0; i != plan->pixelCount; ++i)
                {
                    skymap[0][i] -= log(n_sigmas);
                }

                {
                    double* image;
                    /* printf("rendering image...\n"); */
                    image = (double*) malloc(sizeof(double) * 1024 * 2048);
                    XLALSkymapRenderEqualArea(1024, 2048, image, plan, skymap[0]);

                    /* printf("writing image...\n"); */
                    fwrite(image, sizeof(double), 1024*2048, h);
                    free(image);
                    /* printf("    ...write complete\n"); */
                }

                for (i = 0; i != plan->pixelCount; ++i)
                {
                    accumulated_skymap[i] += skymap[0][i];
                }

                c = 0;
                for (i = 0; i != plan->pixelCount; ++i)
                {
                    if (plan->pixel[i].area > 0)
                    {
                        c += plan->pixel[i].area * counts[i];
                    }
                }

                accumulated_signal = log(0);
                for (i = 0; i != plan->pixelCount; ++i)
                {
                    if (plan->pixel[i].area > 0)
                    {
                        accumulated_signal = XLALSkymapLogSumExp(accumulated_signal,
                            accumulated_skymap[i] + log(plan->pixel[i].area * counts[i] / c));
                    }
                }

                psignal = log(0);
                for (i = 0; i != plan->pixelCount; ++i)
                {
                    if (plan->pixel[i].area > 0)
                    {
                        psignal = XLALSkymapLogSumExp(psignal,
                            skymap[0][i] + log(plan->pixel[i].area * counts[i] / c));
                    }
                }


                /* marginalize over glitch size too */
                for (k = 0; k != n_sigmas; ++k)
                {
                    XLALSkymapGlitchHypothesis(plan, g[k], sigmas[k], w, begin, end, z_t);
                }
                for (k = 0; k != 3; ++k)
                {
                    /* for each detector, marginalize over sigma */
                    g[0][k] = XLALSkymapLogSumExp(g[0][k], XLALSkymapLogSumExp(g[1][k], g[2][k]));
                }

                pglitch =
                    XLALSkymapLogDifferenceExp(XLALSkymapLogSumExp(0, g[0][0] - log(end[0] - begin[0])) +
                    XLALSkymapLogSumExp(0, g[0][1] - log(end[1] - begin[1])) +
                    XLALSkymapLogSumExp(0, g[0][2] - log(end[2] - begin[2])), 0) - log(7.0);
                accumulated_glitch += pglitch;

                printf("    exp(%f)        exp(%f)\n", (psignal - pglitch), (accumulated_signal - accumulated_glitch));

                free(skymap[0]);free(skymap[1]);free(skymap[2]);
                free(modes);
                free(counts);
            }

            /* free data */

            for (d = 0; d != 6; ++d)
            {
                free(z_t[d]);
            }
        } /*else { printf("    VETO\n"); } */
    }

    fclose(h);
    fclose(g);

    free_h_f();
    free_h_t();

    return 0;

}

#endif
#if 0

RandomParams* rng_parameters;

int duration;
int rate;

int samples;

double *strain[3];

XLALSkymapPlanType* plan;

double sigma;

fftw_complex *strain_tilde[3];

int limit;

int band;
int window;

double *filtered[6];

void make_data()
{
    int i, j;
    for (i = 0; i != 3; ++i)
    {
        strain[i] = (double*) malloc(sizeof(double) * samples);
        for (j = 0; j != samples; ++j)
        {
            strain[i][j] = XLALNormalDeviate(rng_parameters);
        }
    }
}

void free_data()
{
    int i;
    for (i = 0; i !=3; ++i)
    {
        free(strain[i]);
    }
}

int index_from_delays(XLALSkymapPlanType* plan, int hl, int hv, int hemisphere)
{
    return
	(hl + plan->hl) +
	(hv + plan->hv) * (plan->hl * 2 + 1) +
	(hemisphere   ) * (plan->hl * 2 + 1) * (plan->hv * 2 + 1);
}

void add_injection()
{
    int i, j;

    double h[2];

    /* for now, just inject perpendicular to network */

    i = index_from_delays(plan, 0, 0, 0);
    j = samples / 2;

    h[0] = sigma * XLALNormalDeviate(rng_parameters); // h+
    h[1] = sigma * XLALNormalDeviate(rng_parameters); // hx

    // hs = s+ + sx
    strain[0][j] += plan->pixel[i].f[0][0] * h[0] + plan->pixel[i].f[0][1] * h[1];
    strain[1][j] += plan->pixel[i].f[1][0] * h[0] + plan->pixel[i].f[1][1] * h[1];
    strain[2][j] += plan->pixel[i].f[2][0] * h[0] + plan->pixel[i].f[2][1] * h[1];

}

void make_strain_tilde()
{
    int i, j;

    fftw_complex *dft_in;
    fftw_complex *dft_out;
    fftw_plan dft_forward;

    dft_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * samples);
    dft_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * samples);

    dft_forward = fftw_plan_dft_1d(samples, dft_in, dft_out, FFTW_FORWARD, FFTW_ESTIMATE);

    for (i = 0; i != 3; ++i)
    {
        for (j = 0; j != samples; ++j)
        {
            dft_in[j][0] = strain[i][j];
            dft_in[j][1] = 0.;
        }
        fftw_execute(dft_forward);
        strain_tilde[i] = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * samples);
        for (j = 0; j != samples; ++j)
        {
            strain_tilde[i][j][0] = dft_out[j][0];
            strain_tilde[i][j][1] = dft_out[j][1];
        }
    }

    fftw_destroy_plan(dft_forward);
    fftw_free(dft_out);
    fftw_free(dft_in);

}

void free_strain_tilde()
{
    int i;
    for (i = 0; i != 3; ++i)
    {
        fftw_free(strain_tilde[i]);
    }
}

void make_filtered(int k)
{
    int i, j;
    fftw_complex *dft_in;
    fftw_complex *dft_out;
    fftw_plan dft_backward;

    printf("window %d...\n", k);

    dft_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * samples);
    dft_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * samples);

    dft_backward = fftw_plan_dft_1d(samples, dft_in, dft_out, FFTW_BACKWARD, FFTW_ESTIMATE);

    for (i = 0; i != 3; ++i)
    {
        double power;
        double total_power;
        double max_power;
        for (j = 0; j != samples; ++j)
        {
            dft_in[j][0] = 0.0;
            dft_in[j][1] = 0.0;
        }
        for (j = k * window; j != (k + 1) * window; ++j)
        {
            dft_in[j][0] = strain_tilde[i][j][0];
            dft_in[j][1] = strain_tilde[i][j][1];
        }
        fftw_execute(dft_backward);
        filtered[i] = (double*) malloc(sizeof(double) * samples);
        filtered[i + 3] = (double*) malloc(sizeof(double) * samples);
        total_power = 0;
        max_power = 0;
        for (j = 0; j != samples; ++j)
        {
            filtered[i][j] = dft_out[j][0] / samples * sqrt(rate/band);
            filtered[i + 3][j] = dft_out[j][1] / samples * sqrt(rate/band);
            power = pow(filtered[i][j], 2) + pow(filtered[i + 3][j], 2);
            total_power += power;
            if ((j > samples / 4) && (j < 3 * samples / 4) && (power > max_power))
            {
                max_power = power;
            }
        }
        total_power /= samples;
        printf("mean(power) = %f, sqrt(max(power)) = %f\n", total_power, sqrt(max_power));
    }

    fftw_destroy_plan(dft_backward);
    fftw_free(dft_out);
    fftw_free(dft_in);
}

void free_filtered()
{
    int i;
    for (i = 0; i != 6; ++i)
    {
        free(filtered[i]);
    }
}

int main(int argc, char** argv)
{
    FILE *h, *g;

    rng_parameters = XLALCreateRandomParams(0);

    duration = 1;
    rate = 4096;

    samples = duration * rate;

    make_data();

    plan = XLALSkymapConstructPlan(rate);

    /* inject a signal */

    sigma = 60.; /* get an SNR of 10 in each of 32 bands? */
    add_injection();
    sigma = 10;

    /* take to frequency domain */
    make_strain_tilde();

    /* process */
    limit = 1024; /* Hz */
    band  = 64; /* Hz */
    window = band * duration;

    h = fopen("raw.dat", "wb");
    g = fopen("bands.dat", "wb");

    {
        int k;
        for (k = 0; k != limit / band; ++k)
        {
            double *skymap, *image;
            int begin[3];
            int end[3];
            double w[3];

            make_filtered(k);

            {
                int gi;
                for (gi = 0; gi != 6; ++gi)
                {
                    fwrite(filtered[gi], sizeof(double), samples, g);
                }
            }

            w[0] = 1;
            w[1] = 1;
            w[2] = 1;
            begin[0] = samples / 2 - rate / 32;
            begin[1] = begin[0];
            begin[2] = begin[0];
            end[0] = samples / 2 + rate / 32;
            end[1] = end[0];
            end[2] = end[0];

            printf("samples [%d, %d)\n", begin[0], end[0]);

            skymap = (double*) malloc(plan->pixelCount * sizeof(double));
            XLALSkymapEllipticalHypothesis(plan, skymap, sigma, w, begin, end, filtered, 0);

            image = (double*) malloc(sizeof(double) * 1024 * 2048);
            XLALSkymapRenderEqualArea(1024, 2048, image, plan, skymap);

            printf("writing image...\n");
            fwrite(image, sizeof(double), 1024*2048, h);
            printf("    ...write complete\n");

            free(image);
            free(skymap);
            free_filtered();
        }
    }

    fclose(g);
    fclose(h);

    free_strain_tilde();
    XLALSkymapDestroyPlan(plan);
    free_data();
    XLALDestroyRandomParams(rng_parameters);

    return 0;
}

#endif

#if 0

static void make_data(double** z, int samples)
{
    int i;
    int j;
    RandomParams* params;
    params = XLALCreateRandomParams(0);
    for (i = 0; i != 6; ++i)
    {
        z[i] = (double*) malloc(samples * sizeof(double));
        for (j = 0; j != samples; ++j)
        {
            z[i][j] = XLALNormalDeviate(params);
        }
    }
    XLALDestroyRandomParams(params);
}

static void free_data(double** z)
{
    int i;
    for (i = 0; i != 6; ++i)
    {
        free(z[i]);
    }
}

static int index_from_delays(XLALSkymapPlanType* plan, int hl, int hv, int hemisphere)
{
    return
	(hl + plan->hl) +
	(hv + plan->hv) * (plan->hl * 2 + 1) +
	(hemisphere   ) * (plan->hl * 2 + 1) * (plan->hv * 2 + 1);
}

void make_injection(XLALSkymapPlanType* plan, double** z, int samples, double sigma)
{
    int i, j;
    RandomParams* params;
    double h[4];
    double theta, phi;
    params = XLALCreateRandomParams(0);

    theta = acos(XLALUniformDeviate(params) * 2 - 1);
    phi = LAL_TWOPI * XLALUniformDeviate(params);

    /* for now, just inject perpendicular to network */

    i = index_from_delays(plan, 0, 0, 0);
    j = samples / 2;

    h[0] = sigma * XLALNormalDeviate(params); // s+
    h[1] = sigma * XLALNormalDeviate(params); // sx
    h[2] = sigma * XLALNormalDeviate(params); // c+
    h[3] = sigma * XLALNormalDeviate(params); // cx

    // hs = s+ + sx
    z[0][j] += plan->pixel[i].f[0][0] * h[0] + plan->pixel[i].f[0][1] * h[1];
    z[1][j] += plan->pixel[i].f[1][0] * h[0] + plan->pixel[i].f[1][1] * h[1];
    z[2][j] += plan->pixel[i].f[2][0] * h[0] + plan->pixel[i].f[2][1] * h[1];
    z[3][j] += plan->pixel[i].f[0][0] * h[2] + plan->pixel[i].f[0][1] * h[3];
    z[4][j] += plan->pixel[i].f[1][0] * h[2] + plan->pixel[i].f[1][1] * h[3];
    z[5][j] += plan->pixel[i].f[2][0] * h[2] + plan->pixel[i].f[2][1] * h[3];

    XLALDestroyRandomParams(params);
}

void make_glitch(XLALSkymapPlanType* plan, double** z, int samples, double sigma)
{
    int i, t[3];
    RandomParams* params;
    params = XLALCreateRandomParams(0);

    i = (rand() % 7) + 1;

    t[0] = rand() % samples;
    t[1] = rand() % samples;
    t[2] = rand() % samples;

    if (i & 1)
    {
        z[0][t[0]] += sigma * XLALNormalDeviate(params);
        z[3][t[0]] += sigma * XLALNormalDeviate(params);
    }
    if (i & 2)
    {
        z[1][t[1]] += sigma * XLALNormalDeviate(params);
        z[4][t[1]] += sigma * XLALNormalDeviate(params);
    }
    if (i & 4)
    {
        z[2][t[2]] += sigma * XLALNormalDeviate(params);
        z[5][t[2]] += sigma * XLALNormalDeviate(params);
    }

    XLALDestroyRandomParams(params);
}

int lose_data(double **z)
{
    int i;
    i = rand() & 7;
    if (i & 1)
    {
        free(z[0]); free(z[3]);
        z[0] = 0; z[3] = 0;
    }
    if (i & 2)
    {
        free(z[1]); free(z[4]);
        z[1] = 0; z[4] = 0;
    }
    if (i & 4)
    {
        free(z[2]); free(z[5]);
        z[2] = 0; z[5] = 0;
    }
    return i;
}

int less_double(void *p, const void *a, const void *b)
{
    if (*((double*) a) < *((double*) b))
    {
        return -1;
    }
    else
    {
        if (*((double*) a) > *((double*) b))
        {
            return 1;
        }
        else
        {
            return 0;
        }
    }
}
void test_XLALHeapIndex()
{
    int n;
    double* p;
    INT4* i;
    int j;

    n = 10;
    p = (double*) malloc(sizeof(double) * n);
    for (j = 0; j != n; ++j)
        p[j] = -j;

    i = (INT4*) malloc(sizeof(INT4) * n);

    XLALHeapIndex(i, p, n, sizeof(double), 0, less_double);

    printf("Unsorted\t'Sorted'\n");
    for(j = 0; j != n; ++j)
        printf("%f\t%f\n", p[j], p[i[j]]);
}

void skymap_index(XLALSkymapPlanType* plan, double* skymap, INT4** pindex)
{
    *pindex = (INT4*) malloc(sizeof(INT4) * plan->pixelCount);
    XLALHeapIndex(*pindex, skymap, plan->pixelCount, sizeof(double), 0, less_double);
}

void cluster_test_raw()
{
    RandomParams* params;
    fftw_complex* dft_in;
    fftw_complex* dft_out;
    fftw_plan dft_forward;
    fftw_plan dft_backward;
    int f;
    int t;
    int i, j;
    int n;
    FILE* h;
    double* x;
    fftw_complex* y;
    int window;

    params = XLALCreateRandomParams(0);

    f = 8192;
    t = 1;
    n = t * f;

    x = (double*) malloc(sizeof(double) * n);
    y = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);

    dft_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
    dft_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);

    dft_forward = fftw_plan_dft_1d(n, dft_in, dft_out, FFTW_FORWARD, FFTW_ESTIMATE);
    dft_backward = fftw_plan_dft_1d(n, dft_in, dft_out, FFTW_BACKWARD, FFTW_ESTIMATE);

    /* construct the data */
    for (i = 0; i != n; ++i)
    {
        x[i] = 0.0; /*XLALNormalDeviate(params); */
    }
    x[n/2] = 1.0;

    /* transform the data to the frequency domain */
    for (i = 0; i != n; ++i)
    {
        dft_in[i][0] = x[i];
        dft_in[i][1] = 0.0;
    }
    fftw_execute(dft_forward);
    for (i = 0; i != n; ++i)
    {
        y[i][0] = dft_out[i][0];
        y[i][1] = dft_out[i][1];
    }
    h = fopen("dump.dat", "wb");
    /* window and inverse transform the data */
    window = 4096 / 16;
    for (j = 0; j < (n / 2); j += window)
    {
        printf("%d\n", j);
        /* zero the input array */
        for (i = 0; i != n; ++i)
        {
            dft_in[i][0] = 0.0;
            dft_in[i][1] = 0.0;
        }
        /* copy over a frequency window */
        for (i = j; i != j + window; ++i)
        {
            dft_in[i][0] = y[i][0];
            dft_in[i][1] = y[i][1];
        }
        /* inverse transform */
        fftw_execute(dft_backward);
        /* write the output */
        fwrite(dft_out, sizeof(fftw_complex), n, h);
    }

    fclose(h);

    fftw_destroy_plan(dft_forward);
    fftw_destroy_plan(dft_backward);
    fftw_free(dft_in);
    fftw_free(dft_out);
    fftw_free(y);
    free(x);

}

void cluster_test()
{
    int duration;
    int rate;
    int n;
    double *x[3];

    duration = 1;
    rate = 4096;
    n = duration * rate;

    {
        int i;
        for (i = 0; i != 3; ++i)
        {
            x[i] = (double*) malloc(sizeof(double) * n);
            for (j = 0; j != n; ++j)
            {
                x[i][j] =
            }
        }
    }


}

int main(int argc, char** argv)
{
    XLALSkymapPlanType* plan;
    double* raw;
    double* z[6];
    double w[3] = { 1, 1, 1};
    double sigma = 10;
    int samples = 512;
    FILE* h;
    int i;
    int begin[3] = { 0, 0, 0 };
    int end[3] = { 512, 512, 512 };
    double* image;
    int f = 8192; /* Hz */
    int bests[5];
    double g[3];
    double marginalizedSkymap;
    double marginalizedGlitch;
    int trial, trials;
    int activeDetectors;
    trials = 100;

    /* test_XLALHeapIndex();
    return 0; */

    cluster_test();
    return 0;

    printf("constructing plan...\n");
    plan = XLALSkymapConstructPlan(f);

    printf("Sigma = %f\n", sigma);
    printf("Injctn  p(Noise)+p(Gltch)+p(Signl)=1.00 log(p(Gltch),p(Signl))/p(Noise)\n");

    for (trial = 0; trial < trials; ++trial)
    {
        int injection_type;

        srand(time(0));

        make_data(z, samples);

        /* printf("injecting signal or glitch...\n"); */

        injection_type = rand() % 3;

        switch (injection_type)
        {
            case 0:
                printf("None  ");
                break;
            case 1:
                printf("Glitch");
                make_glitch(plan, z, samples, sigma);
                break;
            case 2:
                printf("Signal");
                make_injection(plan, z, samples, sigma);
                break;
        }

        /* drop virgo */
        /* z[2] = 0; z[5] = 0; */
        /* drop livingston */
        /* z[1] = 0; z[4] = 0; */
        /* printf("%d:", lose_data(z)); */
        printf(": ");

        /* printf("computing skymap...\n"); */
        raw = (double*) malloc(plan->pixelCount * sizeof(double));
        XLALSkymapEllipticalHypothesis(plan, raw, sigma, w, begin, end, z, bests);

        /* printf("marginalizing over skymap...\n"); */
        marginalizedSkymap = log(0);
        for (i = 0; i != plan->pixelCount; ++i)
        {
            marginalizedSkymap = XLALSkymapLogSumExp(marginalizedSkymap, raw[i]);
        }
        /* printf("Marginalized skymap = %f\n", marginalizedSkymap); */
        /* normalize by the hypothesis count */
        /*printf("Hypothesis count normalization %f - %f\n", marginalizedSkymap, log(bests[4]));
        marginalizedSkymap -= log(bests[4]);*/

        /* printf("computing glitch posterior...\n"); */
        XLALSkymapGlitchHypothesis(plan, g, sigma, w, begin, end, z);

        /* printf("Glitch %f %f %f %f\n", g[0], g[1], g[2], log(samples)); */

        marginalizedGlitch = 0;
        activeDetectors = 0;
        if (z[0])
        {
            marginalizedGlitch += XLALSkymapLogSumExp(0, g[0] - log(samples)) - log(2);
            ++activeDetectors;
        }
        if (z[1])
        {
            marginalizedGlitch += XLALSkymapLogSumExp(0, g[1] - log(samples)) - log(2);
            ++activeDetectors;
        }
        if (z[2])
        {
            marginalizedGlitch += XLALSkymapLogSumExp(0, g[2] - log(samples)) - log(2);
            ++activeDetectors;
        }
        /* remove the signal hypothesis */
        /* marginalizedGlitch = log(exp(marginalizedGlitch) - 0.125) + log(8) - log(7); */
        /* marginalizedGlitch = XLALSkymapLogDifferenceExp(marginalizedGlitch, log(0.125)) + log(8) - log(7); */
        marginalizedGlitch = XLALSkymapLogDifferenceExp(marginalizedGlitch, -log(pow(2,activeDetectors))) + log(pow(2,activeDetectors)) - log(pow(2,activeDetectors) - 1);

        /* printf("Signal %f - Glitch %f\n", marginalizedSkymap, marginalizedGlitch); */
        /* printf("%e %s %e", marginalizedSkymap, (marginalizedSkymap > marginalizedGlitch) ? ">" : "<", marginalizedGlitch); */

        {
            double p[3];
            p[0] = 1 / (1 + exp(marginalizedSkymap) + exp(marginalizedGlitch));
            p[1] = 1 / (exp(-marginalizedGlitch) + 1 + exp(marginalizedSkymap - marginalizedGlitch));
            p[2] = 1 / (exp(-marginalizedSkymap) + exp(marginalizedGlitch - marginalizedSkymap) + 1);
            printf("%f %f %f", p[0], p[1], p[2]);

            if (p[0] > max(p[1], p[2]))
            {
                if (p[0] > 0.95)
                {
                    printf(injection_type == 0 ? " PASS" : " FAIL");
                }
                else
                {
                    printf(" WEAK");
                }
            }
            if (p[1] > max(p[0], p[2]))
            {
                if (p[1] > 0.95)
                {
                    printf(injection_type == 1 ? " PASS" : " FAIL");
                }
                else
                {
                    printf(" WEAK");
                }
            }
            if (p[2] > max(p[1], p[0]))
            {
                if (p[2] > 0.95)
                {
                    printf(injection_type == 2 ? " PASS" : " FAIL");
                }
                else
                {
                    printf(" WEAK");
                }
            }
        }

        printf(" (%f, %f)\n", marginalizedGlitch, marginalizedSkymap);

        /*
        printf("writing directions...\n");
        h = fopen("out.txt", "wt");
        for (i = 0; i != plan->pixelCount; ++i)
        {
            if (plan->pixel[i].area > 0)
            {
                fprintf(h, "%f %f %f\n", plan->pixel[i].direction[0], plan->pixel[i].direction[1], plan->pixel[i].direction[2]);
            }
        }
        fclose(h);
        */

        {
            /* compute confidence region */
            INT4 *index;
            double total_p;
            double total_area;
            int i;

            skymap_index(plan, raw, &index);

            total_p = log(0);
            total_area = 0;
            i = plan->pixelCount - 1;
            while (total_p < marginalizedSkymap - log(2))
            {
                total_p = XLALSkymapLogSumExp(total_p, raw[index[i]]);
                total_area += plan->pixel[index[i]].area;
                --i;
            }

            printf("    %f confidence region area = %f\n", exp(total_p - marginalizedSkymap), total_area);

            free(index);
        }


        if (injection_type == 2)
        {
            printf("rendering image...\n");
            image = (double*) malloc(sizeof(double) * 1024 * 2048);
            XLALSkymapRenderEqualArea(1024, 2048, image, plan, raw);

            printf("writing image...\n");
            h = fopen("raw.dat", "wb");
            fwrite(image, sizeof(double), 1024*2048, h);
            fclose(h);
            free(image);
            printf("    ...write complete\n");
        }

        /* printf("cleanup...\n"); */
        free(raw);
        /* leak the data ? */
        free_data(z);

    }

    XLALSkymapDestroyPlan(plan);

    return 0;
}
#endif
