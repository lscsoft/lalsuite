#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <math.h>

#ifdef SKYMAP_PNG
#include <png.h>
#endif

#include <lal/LALConstants.h>
#include <lal/Skymap.h>

#include <lal/FrameData.h>
#include <lal/FrameStream.h>
#include <lal/LALFrameIO.h>
#include <lal/LALDatatypes.h>
#include <lal/AVFactories.h>
#include <lal/Date.h>

#include <lal/LIGOLwXMLRead.h>
#include <lal/LIGOMetadataTables.h>

/*
 *  File names for input and output
 */
char* h1_frame_file = "none";
char* l1_frame_file = "none";
char* v1_frame_file = "none";

char* xml_file[3] = { "none", "none", "none"};

char* output_file = "skymap.txt";

char* event_id = 0;

/*
 *  Time of the trigger to be analyzed
 */
double trig_time[3] = {0, 0, 0};
double w[3] = { 1, 1, 1 };
/*
 *  Resolution of output sky map
 */
int ra_res = 1024;
int dec_res = 512;

/*
 *  Sampling frequency of analysis, which determines sky tiling
 */
int frequency = 4096;

/*
 *  Number of samples to analyze
 */
int samples = 512;

/*
 *  Array of pointers to the matched-filtered data z[t]
 * 
 *  x[0] = real(z_H)
 *  x[1] = real(z_L)
 *  x[2] = real(z_V)
 *  x[3] = imag(z_H)
 *  x[4] = imag(z_L)
 *  x[5] = imag(z_V)
 */

double* x[6] = { 0, 0, 0, 0, 0, 0 };

double greenwich;

void analyze();

double load_time(char* file, int detector);

void load_data(int detector, char* file, char* initial);

int main(int argc, char** argv)
{

    
    int c;
    while (1)
    {
        static struct option long_options[] =
        {
            {"h1-frame-file", required_argument, 0, 'h'},
            {"l1-frame-file", required_argument, 0, 'l'},
            {"v1-frame-file", required_argument, 0, 'v'},
            {"output-file", required_argument, 0, 'o'},
            {"ra-res", required_argument, 0, 'a'},
            {"dec-res", required_argument, 0, 'd'},
            {"h1-xml-file", required_argument, 0, 't'},
	    {"l1-xml-file", required_argument, 0, 's'},
	    {"v1-xml-file", required_argument, 0, 'r'},
            {"event-id", required_argument, 0, 'e'},
            {0, 0, 0, 0}
        };
        int option_index = 0;
        c = getopt_long_only(argc, argv, "h:l:v:o:a:d:t:s:r:e:", long_options, &option_index);
        if (c == -1)
            break;
        
        switch (c)
        {
            case 'h':
                h1_frame_file = optarg;
                break;
            case 'l':
                l1_frame_file = optarg;
                break;
            case 'v':
                v1_frame_file = optarg;
                break;
            case 'o':
                output_file = optarg;
                break;
            case 'a':
                ra_res = atoi(optarg);
                break;
            case 'd':
                dec_res = atoi(optarg);
                break;
            case 't':
                xml_file[0] = (optarg);
                break;
	    case 's':
	        xml_file[1] = (optarg);
		break;
            case 'r':
	        xml_file[2] = (optarg);
		break;
            case 'e':
                event_id = optarg;
                break;
            default:
                fprintf(stderr, "unknown error while parsing options\n");
                exit(1);
        }            
    }
    
    if (optind < argc)
    {
        fprintf(stderr, "extraneous command line arguments:\n");
        while (optind < argc)
        {
            fprintf(stderr, "\t%s\n", argv[optind++]);
        }
        exit(1);
    }
    
    trig_time[0] = load_time(xml_file[0],0);    
    trig_time[1] = load_time(xml_file[1],1);
    trig_time[2] = load_time(xml_file[2],2); 

    load_data(0, h1_frame_file, "H");
    fprintf(stderr, "H1 trigger:%f\n", trig_time[0]); 
    load_data(1, l1_frame_file, "L");
    load_data(2, v1_frame_file, "V");
    fprintf(stderr, "L1 trigger:%f\n", trig_time[1]);
    fprintf(stderr, "V1 trigger:%f\n", trig_time[2]);

    analyze(ra_res, dec_res, output_file);
 
    return 0;       
}

double load_time(char* file, int detector)
{
    if (strcmp(file, "none"))
    {
        SnglInspiralTable* a = 0;
        LALSnglInspiralTableFromLIGOLw(&a, file, 0, 1);
        w[detector] = sqrt(a->sigmasq);
        /*[detector] = a->snr;*/
        /*w[detector] = 1;*/
        fprintf(stderr, "w[%d] = %f\n", detector, w[detector]);
        greenwich = XLALGreenwichMeanSiderealTime(&(a->end_time));
        fprintf(stderr, "%d:greenwich hour %f\n", detector, fmod(greenwich, LAL_TWOPI));
        return XLALGPSGetREAL8(&(a->end_time));
    }
    else
    {
        fprintf(stderr, "warning: using unit template normalization");
        w[detector] = w[0]*0.1;
        return 0;
    }
}

void load_data(int detector, char* file, char* initial)
{
    char buffer[256];
    double mean_time = 0;
    {
       int i, j;
       j = 0;
       for (i = 0; i != 3; ++i)
       {
          if (trig_time[i] > 0)
          {
             mean_time += trig_time[i];
             ++j;
          }
       }
       if (j > 0)
       {
           mean_time /= j;
       }
       else
       {
           fprintf(stderr, "error: data from at least one interferometer must be supplied");
           exit(1);
       }
    }

    if (strcmp(file, "none"))
    {
        /* 
         *  Read the frame file here
         */

	/*
	 *  This is prototype code that needs to be made robust
	 */
  
        FrStream *stream = NULL;
        COMPLEX8TimeSeries H1series;
        int i;
        int inject_at;
	int start_index;
	double relative_time;

	FILE *FP = NULL;
	sprintf(buffer, "%s:%s.txt", initial, event_id);
	FP = fopen(buffer,"w");

        sprintf(H1series.name,"%s1:CBC-CData_%s", initial, event_id);
	fprintf(stderr, "series name: %s\n", H1series.name);
        stream = XLALFrOpen("./", file);
        H1series.data = XLALCreateCOMPLEX8Vector(4096);
        XLALFrGetCOMPLEX8TimeSeries(&H1series,stream);
        fprintf(stderr,"length %d\n",H1series.data->length);
        for (i = 0; i < H1series.data->length; i++)
        {
         /*   fprintf(FP,"%i %f %f\n",i, H1series.data->data[i].re, H1series.data->data[i].im);*/
        }
        XLALFrClose(stream);
        /*fclose(FP);*/
        /*
         *  Identify the index corresponding to the trigger
	 */
        
	fprintf(stderr, "%s:epoch: %f\n", initial, XLALGPSGetREAL8(&(H1series.epoch)));
        fprintf(stderr, "%s:mean_time: %f\n", initial, mean_time);
        /*relative_time = mean_time - trig_time[detector] + 0.5;*/ 
        relative_time = -trig_time[0] + trig_time[detector] + 0.5;/*- XLALGPSGetREAL8(&(H1series.epoch));*/
        /* FIXME CHAD DID THIS */
        relative_time = 0.5;
	fprintf(stderr,"relative time is %f\n",relative_time);
	start_index = (relative_time * frequency) - (samples / 2);

	/*
	 *  Allocate memory to repack the data in
	 */

        /* real, or waveform one */
	x[detector] = (double*) malloc(samples * sizeof(double));
	/* complex, or waveform two */
	x[detector + 3] = (double*) malloc(samples * sizeof(double));
        inject_at = samples/2 + ((int) ((trig_time[detector]-trig_time[0])*frequency));

	for (i = 0; i != samples; ++i)
	{
	    x[detector    ][i] = H1series.data->data[start_index + i].re;
	    x[detector + 3][i] = H1series.data->data[start_index + i].im;
	    fprintf(FP, "%d %f %f %d\n", i, x[detector][i], x[detector + 3][i], i == inject_at);
        }
        fclose(FP);
    }
    else
    {
        /*
         *  No frame given, generate white noise data
         */
        int i;
        fprintf(stdout, "Warning: generating white noise for detector %d\n", detector);
        for (i = detector; i < 6; i += 3)
        {
            x[i] = (double*) malloc(samples * sizeof(double));
            int j;
            for (j = 0; j != samples; ++j)
            {
                int k;
                x[i][j] = 0.;
                for (k = 0; k != 12; ++k)
                {
                    x[i][j] += rand();
                }
                x[i][j] /= RAND_MAX;
                x[i][j] -= 6;
            }
        }       
    }
}

#define max(A,B) (((A) > (B)) ? (A) : (B))

double logsumexp(double a, double b)
{
   double c;
   if (a > b)
   {
      c = log(1 + exp(b - a)) + a;
   }
   else
   {
      if (a < b)
      {
          c = log(exp(a - b) + 1) + b;
      }
      else
      {
          c = log(exp(a) + exp(b));
      }
   }
   /*if (c == (log(0) - log(0)))*/
   {
      fprintf(stderr, "warning: logsumexp(%f, %f) = %f\n", a, b, c);
   }
   return c;
}

void analyze()
{
    XLALSkymapPlanType* plan;
    double s[8];
    int i, j;
    double* raw;
    double* accumulator;
    double psigma;

    /* 
     *  The characteristic size of the signal 
     */
    /*s = 1e-20;*/
    /*s = 5 / w[0];*/

    s[0] =   1 / w[0];
    s[1] =   2 / w[0];
    s[2] =   4 / w[0];
    s[3] =   8 / w[0];
    s[4] =  16 / w[0];
    s[5] =  32 / w[0];
    s[6] =  64 / w[0];
    s[7] = 128 / w[0];
        
    /* 
     *  The inner product of the template with the noise in each detector 
     */

    /* 
     *  The filtered data z[t]:
     *  x[0] = real(z_H)
     *  x[1] = real(z_L)
     *  x[2] = real(z_V)
     *  x[3] = imag(z_H)
     *  x[4] = imag(z_L)
     *  x[5] = imag(z_V)
     */

     {
        int i, j;
        for (i = 0; i != 6; ++i)
        {
           for (j = 0; j != samples; ++j)
           {
               fprintf(stderr, "x[%d][%d] = %f\n", i, j, x[i][j]);
           }
        }
     }
    
   
    /*
     *  Construct an analysis plan (precomputes the area and direction of 
     *  the sky tiles implied by the frequency) 
     */
    plan = XLALSkymapConstructPlan(frequency);  

    /*
     *  Allocate a chunk of memory tto hold the sky map in the internal 
     *  timing format
     */
    raw = (double*) malloc(plan->pixelCount * sizeof(double));
    accumulator = (double*) malloc(plan->pixelCount * sizeof(double));
    /*
     *  Generate the skymap
     */
    XLALSkymapAnalyzeElliptical(accumulator, plan, s[0], w, samples, x);

    for (i = 1; i != 8; ++i)
    {
        XLALSkymapAnalyzeElliptical(raw, plan, s[i], w, samples, x);

        psigma = log(0);
        for (j = 0; j != plan->pixelCount; ++j)
        {
            if (plan->pixel[j].area != 0)
            {   
                fprintf(stderr, "sigma[%d] raw[%d] = %f\n", i, j, raw[j]);
                accumulator[j] = logsumexp(accumulator[j], raw[j]);
                psigma = logsumexp(psigma, raw[j]);
            }
        }
        fprintf(stderr, "sigma[%d] = %f with ln p(sigma) = %f\n", i, s[i], psigma);
        
    }

    free(raw);
    raw = accumulator;

    #if 0
    {   /*  
         *  Get the mode 
         */ 
        double thetaphi[2];
        XLALSkymapModeThetaPhi(plan, raw,thetaphi);
    }
    #endif
    
    {   /*  
         *  Render the timing-format sky map to a more friendly coordinate 
         *  system
         */
        double* render;
        int m = dec_res;
        int n = ra_res;
        double maximum;
        render = (double*) malloc(m * n * sizeof(double));
        XLALSkymapRenderEquirectangular(m, n, render, plan, raw);
        /*XLALSkymapRenderMollweide(m, n, render, plan, raw);*/

        maximum = render[0];
        for (j = 1; j != m * n; ++j)
        {
            if (render[j] > maximum)
            {
                maximum = render[j];
            }
        }
        for (j = 0; j != m * n; ++j)
        {
            render[j] -= maximum;
        }



#ifdef SKYMAP_PNG
        {
            FILE* fp;
            png_structp png_ptr;
            png_infop info_ptr;
            
            unsigned char* row_ptr;
            int i;
            int j;
            double minimum;
            double maximum;
            
            fp = fopen(output_file, "wb");
            
            png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, 0, 0, 0);
            info_ptr = png_create_info_struct(png_ptr);
            png_init_io(png_ptr, fp);
            
            png_set_IHDR(png_ptr, info_ptr, ra_res, dec_res, 8, 
                PNG_COLOR_TYPE_GRAY, PNG_INTERLACE_NONE, 
                PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
            
            png_write_info(png_ptr, info_ptr);
            
            i = 0;
            while (isinf(render[i]))
                ++i;
            minimum = render[i];
            maximum = render[i];
            for (++i; i != dec_res * ra_res; ++i)
            {
                if (!isinf(render[i]))
                {
                    if (render[i] < minimum)
                        minimum = render[i];
                    if (render[i] > maximum)
                        maximum = render[i];
                }
            }
            
            printf("[%e, %e]\n", minimum, maximum);
            
            row_ptr = malloc(ra_res);
            
            for (i = 0; i != dec_res; ++i)
            {
                for (j = 0; j != ra_res; ++j)
                {
                    row_ptr[j] = (int) (255.0 * (render[i + j * dec_res] - minimum) / (maximum - minimum));
                    /*printf("%d\n",row_ptr[j]);*/
                }
                png_write_row(png_ptr, row_ptr);            
            }
            
            free(row_ptr);
            
            png_write_end(png_ptr, 0);
            png_destroy_write_struct(&png_ptr, &info_ptr);
            
            fclose(fp);
        }
#else        
        {
            /*
             *  Write an ascii file describing the sky map
             */
            int i;
            int j;
            FILE* h;
            h = fopen(output_file, "wt");
            for (j = 0; j != n; ++j)
            {
                double phi, ra;
                phi = (LAL_TWOPI * (j + 0.5)) / n;
                /*ra = fmod(phi + greenwich, LAL_TWOPI);*/
                ra = fmod(phi+greenwich, LAL_TWOPI);
                while (ra < 0)
                    ra += LAL_TWOPI;
                while (ra >= LAL_TWOPI)
                    ra -= LAL_TWOPI;
                for (i = 0; i != m; ++i)
                {                    
                    double dec;
                    dec = LAL_PI_2 - (LAL_PI * (i + 0.5)) / m;
                    fprintf(h, "%.10e %.10e %.10e\n", ra, dec, (render[i + m * j]));
                }
            }
            fclose(h);
        }
#endif    
        free(render);
        
    }
        
    free(raw);
        
    /*
     *  The plan has allocated memory, so we have to destroy it
     */
    XLALSkymapDestroyPlan(plan);
    
    /*
     *  Free the data
     */
    {
        int i;
        for (i = 0; i != 6; ++i)
        {
            free(x[i]);
        }
    }

}



