#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <math.h>
#include <zlib.h>

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

#define max(A,B) (((A) > (B)) ? (A) : (B))

/* declare C99 function */
int (isnan)(double);

/*
 *  File names for input and output
 */
const char* h1_frame_file = 0;
const char* l1_frame_file = 0;
const char* v1_frame_file = 0;
const char* h2_frame_file = 0;

char* xml_file[4] = { 0, 0, 0, 0};

const char* output_file = "skymap.txt";

char* event_id = 0;
/*
 *  Time of the trigger to be analyzed
 */
double w[3] = { 1, 1, 1 };
/*
 *  Resolution of output sky map
 */
int ra_res = 512;
int dec_res = 256;

/*
 *  Sampling frequency of analysis, which determines sky tiling
 */
int frequency = 0;

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
/*
 * Orientation of the earth at the event time
 */
double greenwich = 0;
/*
 * Functions
 */
void load_metadata(char* file, int detector);
void load_data(int detector, const char* file, const char* initial);
void dump_data_to_file(FILE*);
void analyze(void);

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
                {"h2-frame-file", required_argument, 0, 'i'},             
                {"output-file", required_argument, 0, 'o'},
                {"ra-res", required_argument, 0, 'a'},
                {"dec-res", required_argument, 0, 'd'},
                {"h1-xml-file", required_argument, 0, 't'},
                {"l1-xml-file", required_argument, 0, 's'},
                {"v1-xml-file", required_argument, 0, 'r'},
                {"h2-xml-file", required_argument, 0, 'q'},
                {"event-id", required_argument, 0, 'e'},
                {"sample-rate", required_argument, 0, 'f'},
                {0, 0, 0, 0}
            };
        int option_index = 0;
        c = getopt_long_only(argc, argv, "h:l:v:i:o:a:d:t:s:r:q:e:f:", long_options, &option_index);
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
            case 'i':
                h2_frame_file = optarg;
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
                xml_file[0] = optarg;
                break;
            case 's':
                xml_file[1] = optarg;
                break;
            case 'r':
                xml_file[2] = optarg;
                break;
            case 'q':
                xml_file[3] = optarg;
                break;
            case 'e':
                event_id = optarg;
                break;
            case 'f':
                frequency = atoi(optarg);
                break;
            default:
                fprintf(stderr, "unknown error while parsing options\n");
                    exit(1);
        }
    }
    
    if (optind < argc)
    {
        fprintf(stderr, "error: extraneous command line argument ");
        while (optind < argc)
        {
            fprintf(stderr, "%s\n", argv[optind++]);
        }
        exit(1);
    }
    
    /* support "none" arguments */
    if (h1_frame_file && !strcmp("none", h1_frame_file)) { h1_frame_file = 0; }
    if (l1_frame_file && !strcmp("none", l1_frame_file)) { l1_frame_file = 0; }
    if (v1_frame_file && !strcmp("none", v1_frame_file)) { v1_frame_file = 0; }
    if (h2_frame_file && !strcmp("none", h2_frame_file)) { h2_frame_file = 0; }
    if (xml_file[0] && !strcmp("none", xml_file[0])) { xml_file[0] = 0; }
    if (xml_file[1] && !strcmp("none", xml_file[1])) { xml_file[1] = 0; }   
    if (xml_file[2] && !strcmp("none", xml_file[2])) { xml_file[2] = 0; } 
    if (xml_file[3] && !strcmp("none", xml_file[3])) { xml_file[3] = 0; }   
    
    /* check sanity */
    if (!(h1_frame_file || l1_frame_file || v1_frame_file || h2_frame_file))
    {
        fprintf(stderr, "error: Supply at least one of --h1-frame-file, --l1-frame-file or --v1-frame-file\n");
        exit(1);
    }
    if (h1_frame_file && !xml_file[0])
    {
        fprintf(stderr, "error: Supply --h1-xml-file to match --h1-frame-file\n");
        exit(1);
    }
    if (l1_frame_file && !xml_file[1])
    {
        fprintf(stderr, "error: Supply --l1-xml-file to match --l1-frame-file\n");
        exit(1);
    }
    if (v1_frame_file && !xml_file[2])
    {
        fprintf(stderr, "error: Supply --v1-xml-file to match --v1-frame-file\n");
        exit(1);
    }
    if (h2_frame_file && !xml_file[3])
    {
        fprintf(stderr, "error: Supply --h2-xml-file to match --h2-frame-file\n");
        exit(1);
    }
    if (!h1_frame_file && xml_file[0])
    {
        fprintf(stderr, "error: Supply --h1-frame-file to match --h1-xml-file\n");
        exit(1);
    }
    if (!l1_frame_file && xml_file[1])
    {
        fprintf(stderr, "error: Supply --l1-frame-file to match --l1-xml-file\n");
        exit(1);
    }
    if (!v1_frame_file && xml_file[2])
    {
        fprintf(stderr, "error: Supply --v1-frame-file to match --v1-xml-file\n");
        exit(1);
    }
    if (!h2_frame_file && xml_file[3])
    {
        fprintf(stderr, "error: Supply --h2-frame-file to match --h2-xml-file\n");
        exit(1);
    }
    
    if (frequency <= 0)
    {
        fprintf(stderr, "error: Supply positive integer Hertz --sample-rate\n");
        exit(1);
    }
    
    /* handle h2 */
    if (h2_frame_file)
    {
        if (h1_frame_file)
        {
            fprintf(stderr, "warning: H1 and H2 both supplied; H2 will be ignored\n");
        }
        else
        {
            /* rebadge H2 as H1 */
            h1_frame_file = h2_frame_file;
            xml_file[0] = xml_file[3];
        }
    }

    /* examine one second of data around the injection */
    samples = frequency;
    /*
     * Load metadata from XML files
     */
    load_metadata(xml_file[0], 0);
    load_metadata(xml_file[1], 1);
    load_metadata(xml_file[2], 2);
    /*
     * Load z[t] data from frame files
     */
    load_data(0, h1_frame_file, "H");
    load_data(1, l1_frame_file, "L");
    load_data(2, v1_frame_file, "V");
        
    /*
     * Analyze the data and save the skymap
     */

    analyze();
    
    return 0;
}

void load_metadata(char* file, int detector)
{
    if (file)
    {
        SnglInspiralTable* a = 0;
        LALSnglInspiralTableFromLIGOLw(&a, file, 0, 1);
        if (!a)
        {
            fprintf(stderr, "error: failed to read single inspiral table from file %s\n", file);
            exit(1);
        }
        w[detector] = sqrt(a->sigmasq);
        greenwich = fmod(XLALGreenwichMeanSiderealTime(&(a->end_time)), LAL_TWOPI);
        /*
        fprintf(stderr, "GPS %d -> GMS %e -> RAD %e \n", a->end_time.gpsSeconds, XLALGreenwichMeanSiderealTime(&(a->end_time)), greenwich);
         */
    }
}

void load_data(int detector, const char* file, const char* initial)
{
    if (file)
    {
        /* 
         *  Read the frame file here
         */

        FrStream *stream = NULL;
        COMPLEX8TimeSeries H1series;
        int i;

        sprintf(H1series.name,"%s1:CBC-CData_%s", initial, event_id);
        stream = XLALFrOpen("./", file);
        if (!stream)
        {
            fprintf(stderr, "error: failed to open FrStream from file %s\n", H1series.name);
            exit(1);
        }
        
        H1series.data = XLALCreateCOMPLEX8Vector(samples);
        if (!H1series.data)
        {
            fprintf(stderr, "error: failed to create COMPLEX8 vector\n");
            exit(1);
        }
        XLALFrGetCOMPLEX8TimeSeries(&H1series,stream);
        XLALFrClose(stream);

        /*
         *  Allocate memory to repack the data in
         */

        /* real, or waveform one */
        x[detector] = (double*) malloc(samples * sizeof(double));
        /* complex, or waveform two */
        x[detector + 3] = (double*) malloc(samples * sizeof(double));     

        for (i = 0; i != samples; ++i)
        {
            x[detector    ][i] = H1series.data->data[i].re;
            x[detector + 3][i] = H1series.data->data[i].im;
        }
    }
}

#define NSIGMA 3

void analyze(void)
{

    XLALSkymapPlanType* plan;
    double s[NSIGMA];
    int i, j;
    double* raw;
    double* accumulator;
    int begin[3], end[3];
    
    double min_w;

    begin[0] = 0;
    begin[1] = begin[0];
    begin[2] = begin[0];

    end[0] = samples;
    end[1] = end[0];
    end[2] = end[0];

    /* 
     *  The characteristic size of the signal 
     */

    /*fprintf(stderr, "w: %e %e %e\n", w[0], w[1], w[2]);*/
    
    {
        /* Validate the input */
        
        for (i = 0; i != 3; ++i)
        {
            if (x[i])
            {
                if (!((w[i] > 0) && (w[i] < -log(0))))
                {
                    fprintf(stderr, "error: w[%d] is %e\n", i, w[i]);
                    exit(1);
                }
                for (j = 0; j != samples; ++j)
                {
                    if (!((x[i][j] > log(0)) && (x[i][j] < -log(0))))
                    {
                        fprintf(stderr, "error: x[%d][%d] = %e\n", i, j, x[i][j]);
                        exit(1);
                    }
                    if (!((x[i + 3][j] > log(0)) && (x[i + 3][j] < -log(0))))
                    {
                        fprintf(stderr, "error: x[%d][%d] = %e\n", i + 3, j, x[i + 3][j]);
                        exit(1);
                    }
                }
            }
        }
        
    }
    
    {
        /*
         * Compute the glitch hypotheses
         */
        double pglitches[3];
        double pglitch;
        
        for (i = 0; i != 3; ++i)
        {
            pglitches[i] = 0;
            if (x[i])
            {
                for (j = 0; j != samples; ++j)
                {
                    int k;
                    for (k = 0; k != NSIGMA; ++k)
                    {
                        /*
                         * FIXME: exp may overflow
                         */
                        double p = pow(s[k], -2) * exp(0.5 * (1.0 - pow(s[k], -2)) * (pow(x[i][j], 2) + pow(x[i + 3][j], 2)));
                        pglitches[i] += (p / (NSIGMA * samples));
                    }
                }
            }            
        }

        pglitch = 1.0;
        for (i = 0; i != 3; ++i)
        {
            if (x[i])
            {
                pglitch *= pglitches[i];
            }
        }
        
        /*
         * FIXME: output pglitch somewhere
         */
        
        /*
         * FIXME: implications of glitch model
         */
    }
    
    {
        /* 
         * Use the most sensitive detector to set the bottom of the
         * logarithmic distribution of amplitude scales 
         */
        double min_w;
        
        min_w = -log(0);
        
        for (i = 0; i != 3; ++i)
        {
            if (x[i])
            {
                if (w[i] < min_w)
                {
                    min_w = w[i];
                }
            }
        }
        
        s[0] =   sqrt(10.0) / min_w;
        s[1] =  sqrt(100.0) / min_w;
        s[2] = sqrt(1000.0) / min_w;
    }
    
    {
        /* Convert the individually normalized matched filters to a common normalization */
        for (i = 0; i != 3; ++i)
        {
            if (x[i])
            {
                for (j = 0; j != samples; ++j)
                {
                    x[i][j] *= w[i];
                    x[i + 3][j] *= w[i];
                }
            }
        }
    }
     
    /*   
     *  the sky tiles implied by the frequency) 
     */
    plan = XLALSkymapConstructPlanMN(frequency, dec_res, ra_res);  
    /*
     *  Allocate a chunk of memory tto hold the sky map in the internal 
     *  timing format
     */
    raw = (double*) calloc(plan->pixelCount, sizeof(double));
    accumulator = (double*) calloc(plan->pixelCount, sizeof(double));
    /*
     *  Generate the skymap
     */
    XLALSkymapEllipticalHypothesis(plan, accumulator, s[0], w, begin, end, x, 0);
    for (i = 1; i != NSIGMA; ++i)
    {

        XLALSkymapEllipticalHypothesis(plan, raw, s[i], w, begin, end, x, 0);
        XLALSkymapSum(plan, accumulator, accumulator, raw);
    }

    free(raw);
    raw = accumulator;
    
    {
        /* validate the output */
        for (i = 0; i != plan->pixelCount; ++i)
        {
            if ((plan->pixel[i].area != 0) && isnan(raw[i]))
            {
                fprintf(stderr, "plan->pixel[%d].area = %e && raw[%d] = %e\n", i, plan->pixel[i].area, i, raw[i]);
                exit(1);
            }
        }
    }

    {   /*  
         *  Get the mode 
         */ 
        double thetaphi[2];
        XLALSkymapModeThetaPhi(plan, raw, thetaphi);
        while (thetaphi[1] < 0)
            thetaphi[1] += LAL_TWOPI;
        printf("%e %e\n", thetaphi[0], thetaphi[1]);
    }
    
    {   /*  
         *  Render the timing-format sky map to a more friendly coordinate 
         *  system
         */
        double* render;
        int m = dec_res;
        int n = ra_res;
        double maximum;
        render = (double*) malloc(m * n * sizeof(double));
        if (XLALSkymapRender(render, plan, raw))
        {
            fprintf(stderr, "XLALSkymapRender failed\n");
            exit(1);
        }
        
        for (j = 0; j != m * n; ++j)
        {
            if (isnan(render[j]))
            {
                fprintf(stderr, "WARNING: The rendered skymap contains a NaN at pixel %d\n", j);
            }
        }

        maximum = log(0);
        for (j = 0; j != m * n; ++j)
        {
            if (isnan(maximum))
            {
                fprintf(stderr, "ERROR: The rendered skymap's progressive maximum became NaN before pixel %d\n", j);
                exit(1);
            }               
            if (render[j] > maximum)
            {
                maximum = render[j];
            }
        }
        if (!(maximum < -log(0)))
        {
            fprintf(stderr, "ERROR: The skymap maximum was not finite, and normalization cannot be performed\n");
            exit(0);
        }
        for (j = 0; j != m * n; ++j)
        {
            render[j] -= maximum;
        }



        {
            /*
             *  Write a gzipped ascii file describing the sky map
             */
            gzFile* h = NULL;
            h = gzopen(output_file, "wb");
            if (h == NULL) 
            {
                fprintf(stderr,"cannot gzopen() %s\n",output_file);
                exit(1);
            }
            for (j = 0; j != n; ++j)
            {
                double phi, ra;
                phi = (LAL_TWOPI * (j + 0.5)) / n;
                ra = fmod(phi + greenwich, LAL_TWOPI);
                while (ra < 0)
                    ra += LAL_TWOPI;
                while (ra >= LAL_TWOPI)
                    ra -= LAL_TWOPI;
                for (i = 0; i != m; ++i)
                {
                    double dec;
                    dec = acos(1. - (i + 0.5) * (2. / m));
                    gzprintf(h, "%.10e %.10e %.10e\n", ra, dec, (render[i + m * j]));
                }
            }
            gzclose(h);
        }
        
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
        for (i = 0; i != 6; ++i)
        {
            free(x[i]);
        }
    }

}

