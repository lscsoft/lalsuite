#include <malloc.h>
#include <stdio.h>
#include <math.h>

#include <lal/Skymap.h>
#include <lal/Random.h>

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
    params = XLALCreateRandomParams(0);
    
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
    int i, j;
    RandomParams* params;
    params = XLALCreateRandomParams(0);
    
    i = index_from_delays(plan, 0, 0, 0);
    j = samples / 2;
    
    z[0][j] += sigma * XLALNormalDeviate(params);
    z[1][j] += sigma * XLALNormalDeviate(params);
    z[2][j] += sigma * XLALNormalDeviate(params);
    z[3][j] += sigma * XLALNormalDeviate(params);
    z[4][j] += sigma * XLALNormalDeviate(params);
    z[5][j] += sigma * XLALNormalDeviate(params);
    
    XLALDestroyRandomParams(params);
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
    double marginalizedSkymap = log(0);
    double marginalizedGlitch;
    
    make_data(z, samples);    
    
    printf("constructing plan...\n");
    plan = XLALSkymapConstructPlan(f);
    
    printf("injecting signal or glitch...\n");
    make_injection(plan, z, samples, sigma);
    
    printf("computing skymap...\n");
    raw = (double*) malloc(plan->pixelCount * sizeof(double));
    XLALSkymapEllipticalHypothesis(plan, raw, sigma, w, begin, end, z, bests); 
    
    printf("marginalizing over skymap...\n");
    for (i = 0; i != plan->pixelCount; ++i)
    {
        marginalizedSkymap = XLALSkymapLogSumExp(marginalizedSkymap, raw[i]);
    }
    printf("Marginalized skymap = %f\n", marginalizedSkymap);
    /* normalize by the hypothesis count */
    /*printf("Hypothesis count normalization %f - %f\n", marginalizedSkymap, log(bests[4]));
    marginalizedSkymap -= log(bests[4]);*/
    
    printf("computing glitch posterior...\n");
    XLALSkymapGlitchHypothesis(plan, g, sigma, w, begin, end, z);
    
    printf("Glitch %f %f %f %f\n", g[0], g[1], g[2], log(samples));
    
    marginalizedGlitch = 
            XLALSkymapLogSumExp(0, g[0] - log(samples)) + 
            XLALSkymapLogSumExp(0, g[1] - log(samples)) + 
            XLALSkymapLogSumExp(0, g[2] - log(samples)) -
            log(8);
    /* remove the signal hypothesis */
    marginalizedGlitch = log(exp(marginalizedGlitch) - 0.125) + log(8) - log(7);
    
    printf("Signal %f - Glitch %f\n", marginalizedSkymap, marginalizedGlitch);
            
    
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
    
    printf("rendering image...\n");
    image = (double*) malloc(sizeof(double) * 512 * 1024);
    XLALSkymapRenderEqualArea(512, 1024, image, plan, raw);
    
    printf("writing image...\n");
    h = fopen("raw.dat", "wb");
    fwrite(image, sizeof(double), 512*1024, h);
    fclose(h);
    
    printf("cleanup...\n");
    free(raw);
    XLALSkymapDestroyPlan(plan);

    /* z leaks */
    
    return 0;
}
