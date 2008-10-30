#include <malloc.h>
#include <stdio.h>

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
	    z[i][j] = XLALNormalDeviate(params) * 1;
	}
    }
    XLALDestroyRandomParams(params);
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

    make_data(z, samples);
    
    plan = XLALSkymapConstructPlan(4096);
    raw = (double*) malloc(plan->pixelCount * sizeof(double));
    XLALSkymapAnalyzeElliptical(raw, plan, sigma, w, samples, z);
    
    /*
    h = fopen("out.txt", "wt");
    for (i = 0; i != plan->pixelCount; ++i)
    {
        if (plan->pixel[i].area > 0)
	{
	    fprintf(h, "raw[%d] = %f\n", i, raw[i]);
	}
    }

    fclose(h);
    */
    free(raw);
    XLALSkymapDestroyPlan(plan);
    
    return 0;
}
