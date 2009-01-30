#include <malloc.h>
#include <stdio.h>
#include <math.h>

#include <lal/LALConstants.h>
#include <lal/Skymap.h>
#include <lal/Random.h>

#define max(A,B) (((A) > (B)) ? (A) : (B))

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

double logdifferenceexp(double a, double b)
{
    return a + log(1. - exp(b - a));
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

int main(int argc, char** argv)
{
    XLALSkymapPlanType* plan;
    double* raw;
    double* z[6];
    double w[3] = { 1, 1, 1};
    double sigma = 10.;
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
    trials = 100;
     
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

        marginalizedGlitch = 
                XLALSkymapLogSumExp(0, g[0] - log(samples)) + 
                XLALSkymapLogSumExp(0, g[1] - log(samples)) + 
                XLALSkymapLogSumExp(0, g[2] - log(samples)) -
                log(8);
        /* remove the signal hypothesis */
        /* marginalizedGlitch = log(exp(marginalizedGlitch) - 0.125) + log(8) - log(7); */
        marginalizedGlitch = logdifferenceexp(marginalizedGlitch, log(0.125)) + log(8) - log(7);

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
        /* free_data(z); */
        
    }
    
    XLALSkymapDestroyPlan(plan);
    
    return 0;
}
