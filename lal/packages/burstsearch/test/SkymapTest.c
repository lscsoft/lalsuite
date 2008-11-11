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
    double* z[6];
    double* x[6];
    double w[3] = { 1, 1, 1};
    double sigma = 10;
    int samples = 512;
    FILE* h;
    int i;
    double* rbw;
    int begin[3] = { 200, 200, 200 };
    int end[3] = { 300, 300, 210 };
    double* image;
    
    make_data(z, samples);
    
    x[0] = z[0] + begin[0];
    x[1] = z[1] + begin[1];
    x[2] = z[2] + begin[2];
    x[3] = z[3] + begin[0];
    x[4] = z[4] + begin[1];
    x[5] = z[5] + begin[2];
    
    plan = XLALSkymapConstructPlan(2048);
    raw = (double*) malloc(plan->pixelCount * sizeof(double));
    XLALSkymapAnalyzeElliptical(raw, plan, sigma, w, samples, z);
    rbw = (double*) malloc(plan->pixelCount * sizeof(double));
    XLALSkymapEllipticalHypothesis(plan, rbw, sigma, w, begin, end, x); 
    
    h = fopen("out.txt", "wt");
    for (i = 0; i != plan->pixelCount; ++i)
    {
        if (plan->pixel[i].area > 0)
	{
	    /*fprintf(h, "raw[%d] = %f\n", i, raw[i]);*/
	    fprintf(h, "%f %f %f\n", plan->pixel[i].direction[0], plan->pixel[i].direction[1], plan->pixel[i].direction[2]);
	}
    }
    fclose(h);

    h = fopen("difference.txt", "wt");
    for (i = 0; i != plan->pixelCount; ++i)
    {
	if (plan->pixel[i].area > 0)
	{
	    fprintf(h, "%f\n", (raw[i] - rbw[i])/abs(raw[i]));
	}
    }
    fclose(h);

    image = (double*) malloc(sizeof(double) * 512 * 1024);
    XLALSkymapRenderEqualArea(512, 1024, image, plan, rbw);
    h = fopen("raw.dat", "wb");
    fwrite(image, sizeof(double), 512*1024, h);
    fclose(h);

    
    free(rbw);
    free(raw);
    XLALSkymapDestroyPlan(plan);
    
    return 0;
}
