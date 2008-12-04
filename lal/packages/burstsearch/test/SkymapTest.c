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
	    z[i][j] = XLALNormalDeviate(params);
	}
    }
    XLALDestroyRandomParams(params);
}

int main(int argc, char** argv)
{
    XLALSkymapPlanType* plan;
    double* raw;
    double* a[6];
    double* b[6];
    double w[3] = { 1, 1, 1};
    double sigma = 10;
    int samples = 512;
    FILE* h;
    int i;
    int begin[3] = { 0, 0, 0 };
    int end[3] = { 512, 512, 512 };
    double* image;
    
    make_data(a, samples);
    make_data(b, samples);
    
    printf("constructing plan...\n");
    plan = XLALSkymapConstructPlan(2048*4);
    printf("computing skymap...\n");
    raw = (double*) malloc(plan->pixelCount * sizeof(double));
    XLALSkymapEllipticalHypothesis(plan, raw, sigma, w, begin, end, a, 0); 
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

    {
	double *raw_b;
	raw_b = (double*) malloc(plan->pixelCount * sizeof(double));
	XLALSkymapEllipticalHypothesis(plan, raw_b, sigma, w, begin, end, b, 0);
	XLALSkymapSum(plan, raw, raw, raw_b);
	free(raw_b);
    }

    printf("rendering image...\n");
    image = (double*) malloc(sizeof(double) * 512 * 1024);
    XLALSkymapRenderEqualArea(512, 1024, image, plan, raw);
    printf("writing image...\n");
    h = fopen("raw.dat", "wb");
    fwrite(image, sizeof(double), 512*1024, h);
    fclose(h);
    
    printf("cleanup...\n");
    free(raw);
    /*free(raw);*/
    XLALSkymapDestroyPlan(plan);
    
    return 0;
}
