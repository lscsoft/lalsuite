/*----------------------------------------------------------------------- 
 * 
 * File Name: narrowBandExtract.c
 *
 * Author: Koranda, S.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <math.h>

#include <FrameL.h>

/* SFT header structure as of 7/21/2004 */
struct sftheader
{
  double endian;        /* endianess indicator, should be 1.0 */
  int    gpssec;        /* GPS start time seconds */
  int    gpsnan;        /* GPS start time residual nanoseconds */
  double tbase;         /* time base for the SFT or how many seconds of data represented */
  int    first;         /* index for first data point so that frequency is tbase/index */
  int    nsamples;      /* number of complex data points */
};

/* function prototypes */
int parse_command_line(int argc, char *argv[]);
int fileURL_to_localPath(char *url, char *path);
int read_sft( struct sftheader *header, float **dataF, const char *fname, float fmin, float deltaf );
void byte_swap( void *ptr, size_t size, size_t nobj );
int add_frame_to_frame_file(struct sftheader *header, float *dataF, int frnum, FrFile *outputFrameFile); 
int add_sft_to_mergedSFT_file(struct sftheader *header, float *dataF, FILE *outputSFTFile); 

/* global variables */
char inputFilePath[256] = "";   /* path to input file containing list of paths to wide-band SFT files */
char outputFilePath[256] = "";  /* path to output file, a many-frame frame file */
float startFrequency = 0.0;     /* desired starting frequency for the narrow band */
float bandwidth = 0.0;          /* width of the frequency band in the narrow band */
int mergedSFT = 0;              /* flag for generating merged SFT instead of frame */
float mysteryFactor = 1.0;      /* factor used to encrypt the data */

int main ( int argc, char *argv[] )
{
        char sftFilePath[256];          /* path to a wide-band SFT file */
        char url[256];                  /* a file:// URL */
        FILE *inputFilePathFILE = NULL; /* file object for the SFT file */
        FrFile *outputFrameFile = NULL; /* file object for the output frame file */
        FILE *outputSFTFile = NULL;     /* file object for the output merged SFT file */
        int dbglvl = 0;                 /* debug level for the Fr library */
        struct sftheader header;        /* header of the SFT under consideration */
        float *dataF = NULL;            /* pointer to frequency series data */
        int ret;                        
        int frnum = 0;                  /* frame number used internally in frame object */

        /* parse command line and populate global variables */
        parse_command_line(argc, argv);

        /* initialize the Fr library */
        if (!mergedSFT)
        {
                FrLibIni(NULL, stderr, dbglvl); 
        }
        

        /* open the input file containing paths to wide-band SFT files */
        inputFilePathFILE = fopen(inputFilePath, "r");
        if (inputFilePathFILE == NULL)
        {
                fprintf(stderr, "Error opening file %s\n", inputFilePath);
                exit(1);
        }

        /* open a new frame file for the frequency series stored in frames */
        if (!mergedSFT)
        {
                outputFrameFile = FrFileONew(outputFilePath, 0);
                if (outputFrameFile == NULL)
                {
                        fprintf(stderr, "Error opening file %s\n", outputFilePath);
                        exit(1);
                }
        }
        else
        {
                outputSFTFile = fopen(outputFilePath, "wb");
                if (outputSFTFile == NULL)
                {
                        fprintf(stderr, "Error opening file %s\n", outputFilePath);
                        exit(1);
                }
        }

        /* loop over wide-band SFT files, opening each, extracting the frequency
           series of interest, and writing it as a frame to the frame file */
        while ( 1 )
        {
                ret = fscanf(inputFilePathFILE, "%s\n", url);

                /* break when no more wide-band SFT files to process */
                if (ret == EOF)
                {
                        break;
                }

                /* convert file URLs to local path */
                fileURL_to_localPath(url, sftFilePath);
                
                /* fprintf(stdout, "%s\n", sftFilePath); */

                ret = read_sft(&header, &dataF, sftFilePath, startFrequency, bandwidth);
                if ( ret != 0)
                {
                        fprintf(stderr, "Error reading SFT file %s\n", sftFilePath);
                        exit(1);
                }

                if (!mergedSFT)
                {                
                        ret = add_frame_to_frame_file(&header, dataF, frnum++, outputFrameFile); 
                }
                else
                {
                        ret = add_sft_to_mergedSFT_file(&header, dataF, outputSFTFile); 
                }
                if ( ret != 0)
                {
                        fprintf(stderr, "Error writing frequency series to file %s\n", outputFilePath);
                        exit(1);
                }

                /* need to free memory allocated in read_sft */
                if (dataF)
                {
                        free(dataF);
                }
        } /* while */

        /* close the output file and exit */
        if (!mergedSFT)
        {
                FrFileOEnd(outputFrameFile);
        }
        else
        {
                fclose(outputSFTFile);
        }

        return 0;

}

#define USAGE \
"NAME                                                                       \n"\
"       lalapps_narrowBandExtract                                           \n"\
"                                                                           \n"\
"SYNOPSIS                                                                   \n"\
"       lalapps_narrowBandExtract --input=PATH --output=PATH                \n"\
"               --start-frequency=FREQUENCY --bandwidth=BANDWIDTH           \n"\
"               --mergedSFT  --mystery-factor=FACTOR                        \n"\
"                                                                           \n"\
"       lalapps_narrowBandExtract --help                                    \n"\
"                                                                           \n"\
"DESCRIPTION                                                                \n"\
"       Given an input file with a list of paths to SFT files, read         \n"\
"       each SFT and extract the appropriate narrow band. Output the        \n"\
"       collection of narrow band frequency series as a frame file          \n"\
"       containing a frame for each frequency series.                       \n"\
"                                                                           \n"\
"       -i, --input                                                         \n"\
"               path to file containing list of SFT file paths              \n"\
"                                                                           \n"\
"       -o, --output                                                        \n"\
"               path to output file into which to write the frame file      \n"\
"                                                                           \n"\
"       -f, --start-frequency                                               \n"\
"               starting frequency of the narrow band to extract            \n"\
"                                                                           \n"\
"       -b, --bandwidth                                                     \n"\
"               the bandwith of the narrow band to extract                  \n"\
"                                                                           \n"\
"       -m, --mergedSFT                                                     \n"\
"               generate a merged SFT file instead of a frame               \n"\
"               file                                                        \n"\
"                                                                           \n"\
"       -y, --mystery-factor                                                \n"\
"               secret factor that multiplies the data to make              \n"\
"               it hard to infer strain values                              \n"\
"                                                                           \n"\
"       -h, --help                                                          \n"\
"               print this usage message                                    \n"\
"                                                                           \n"\
"EXAMPLE                                                                    \n"\
"                                                                           \n"\
"$ lalapps_narrowBandExtract --input=SFTfilePaths                           \n"\
"       --output=H1-narrow_100.0_1.0-777777777-1800.gwf                     \n"\
"       --start-frequency=100.0 --bandwidth=1.0                             \n"\
"\n"

/*
Parse the command line arguments and populate necessary
global variables.
*/
int parse_command_line(int argc, char *argv[])
{
        struct option long_options[] =
        {
                {"input",               required_argument, 0,   'i'},
                {"output",              required_argument, 0,   'o'},
                {"start-frequency",     required_argument, 0,   'f'},
                {"bandwidth",           required_argument, 0,   'b'},
                {"mergedSFT",           no_argument,       0,   'm'},
		{"mystery-factor",      required_argument, 0,   'y'},
                {"help",                 no_argument,       0,   'h'},
                {0, 0, 0, 0}
        };

        char *short_options = "i:o:f:b:y:hm";
        int c;

        while ( 1 )
        {
                int option_index = 0;

                c = getopt_long( argc, argv, short_options, long_options, &option_index);

                if (c == -1)
                {
                        break;
                }


                switch (c)
                {
                        case 'i':
                                {
                                        strcpy(inputFilePath, optarg); 
                                        break;
                                }

                        case 'o':
                                {
                                        strcpy(outputFilePath, optarg);
                                        break;
                                }

                        case 'f':
                                {
                                        startFrequency = atof(optarg);
                                        break;
                                }

                        case 'y':
                                {
                                        mysteryFactor = atof(optarg);
                                        break;
                                }

                        case 'b':
                                {
                                        bandwidth = atof(optarg);
                                        break;
                                }
        
                        case 'm':
                                {
                                        mergedSFT = 1;
                                        break;
                                }

                        case 'h':
                                {
                                        fprintf(stderr, USAGE);
                                        exit(1);
                                }

                        default:
                                fprintf(stderr, "unknown error while parsing options\n");
                                fprintf(stderr, "Use --help for usage\n");
                                exit(1);

                } /* switch */

        } /* while */

        /* input file must be specified */
        if (strcmp(inputFilePath, "") == 0)
        {
                fprintf(stderr, "Input filename must be specified\n");
                fprintf(stderr, "Use --help for usage\n");
                exit(1);
        }

        /* output file must be specified */
        if (strcmp(outputFilePath, "") == 0)
        {
                fprintf(stderr, "Output filename must be specified\n");
                fprintf(stderr, "Use --help for usage\n");
                exit(1);
        }


        /* starting frequency must be specified */
        if (startFrequency == 0.0)
        {
                fprintf(stderr, "Starting frequency must be specified\n");
                fprintf(stderr, "Use --help for usage\n");
                exit(1);
        }

        /* bandwidth frequency must be specified */
        if (bandwidth == 0.0)
        {
                fprintf(stderr, "Bandwidth must be specified\n");
                fprintf(stderr, "Use --help for usage\n");
                exit(1);
        }


        return 0;
}


/*
Convert a file: URL to a local path. Note that this
is VERY simplistic and does little error checking.
*/
int fileURL_to_localPath(char *url, char *path)
{
        char *a;

        /* URL has form file://host/path */
        if (strncmp(url, "file://", 7) == 0)
        {
                /* find first slash after file:// */
                a = strchr(url + 7, '/');

                /* copy from first slash onward */
                strcpy(path, a);
                
        }

        /* URL has form file:/path */
        else if (strncmp(url, "file:/", 6) == 0)
        {
                strcpy(path, url + 5);

        }
        else
        {
                fprintf(stderr, "Bad file URL: %s", url);
                exit(1);
        }

        return 0;
}


/* read sft-format file fname; return header and allocate memory for data */
/*
Allocate memory and read in frequency series covering the range fmin to 
fmin + delta f. Also set the SFT header structure appropriately.

This code was taken almost completely from the file sft2frame.c written
by Jolien Creighton. Errors are owned by Scott Koranda.
*/
int read_sft( struct sftheader *header, float **dataF, const char *fname, float fmin, float deltaf )
{
        FILE *fp;               /* pointer to SFT wide-band file */
        int swapbytes = 0;      /* flag for whether indianess is wrong and need to swap */
        int bytes;              /* bytes necessary to hold comples freq series of interest */
        int code;               /* return code */
        int i;
        double tbase;           /* time period covered by the frequency series */
        float freqmin;          /* minimum frequency in the frequency series */
        float freqmax;          /* maximum frequency in the frequency series */
        int fminindex;          /* index associated with minimum frequency */
        int fmaxindex;          /* index associated with maximum frequency */
        long seekpos;           /* seek position in wide-band SFT file */
        float norm;             /* normalization applied per Xavi's instruction */

        if ( ! header || ! dataF || !fname )
        {
                fprintf( stderr, "Invalid input\n" );
                return 1;
        }
        fp = fopen( fname, "r" );
        if ( ! fp )
        {
                fprintf( stderr, "Could not read file %s\n", fname );
                return 1;
        }
        code = fread( header, sizeof( *header ), 1, fp );
        if ( code != 1 )
        {
                fprintf( stderr, "No SFT header in file %s\n", fname );
                fclose( fp );
                return 1;
        }
        if ( header->endian != 1.0 )
        {
                byte_swap( &header->endian, sizeof( header->endian ), 1 );
                if ( header->endian != 1.0 )
                {
                        fprintf( stderr, "Corrupt header or wrong endian in file %s\n", fname );
                        fclose( fp );
                        return 2;
                }
                else
                {
                        swapbytes = 1;
                        byte_swap( &header->gpssec, sizeof( header->gpssec ), 1 );
                        byte_swap( &header->gpsnan, sizeof( header->gpsnan ), 1 );
                        byte_swap( &header->tbase, sizeof( header->tbase ), 1 );
                        byte_swap( &header->first, sizeof( header->first ), 1 );
                        byte_swap( &header->nsamples, sizeof( header->nsamples ), 1 );
                }
        }
        if ( header->tbase <= 0 )
        {
                fprintf( stderr, "Timebase %f from data file %s must be positive\n",
                        header->tbase, fname );
                fclose( fp );
                return 3;
        }

        tbase = header->tbase;

        /* compute the minimum and maximum frequencies available in the wide-band series */
        freqmin = (float) (header->first / tbase);
        freqmax = (float) ((header->first + header->nsamples) / tbase);

        if ( fmin < freqmin )
        {
                fprintf(stderr, "Desired start frequency %f too small: SFT begins with %f\n", fmin, freqmin);
                return 4;
        }

        if ( (fmin + deltaf) > freqmax )
        {
                fprintf(stderr, "Desired end frequency %f too large: SFT ends with %f\n", fmin, freqmax);
                return 5;
        }

        /* compute the min and max index necessary to ensure fmin to fmin + deltaf
           is in extracted band */
        fminindex = (int) floor(fmin * tbase);
        fmaxindex = (int) ceil((fmin + deltaf) * tbase);

        /* allocate data to store the narrow band complex series
           this needs to be freed by the caller */
        bytes = 2 * (fmaxindex - fminindex + 1) * sizeof( float );
        *dataF = malloc( bytes );
        if ( ! *dataF )
        {
                fprintf( stderr, "Allocation error\n" );
                fclose( fp );
                return 5;
        }

        /* seek to the correct position to find the data point for fminindex
           assuming 4 bytes for a real and 2 for complex series */
        seekpos = (long) ((fminindex - (header -> first)) * 4 * 2);
        
        code = fseek(fp, seekpos, SEEK_CUR);
        if (code != 0)
        {
                fprintf(stderr, "Could not seek to position %ld in SFT file %s", seekpos, fname);
                fclose(fp);
                return 6;
        }

        /* read in the narrow band series and swap if necessary */
        code = fread( *dataF, bytes, 1, fp );
        if ( code != 1 )
        {
                fprintf( stderr, "Read error from file %s\n", fname );
                fclose( fp );
                return 7;
        }
        if ( swapbytes )
                byte_swap( *dataF, sizeof( **dataF ), 2 * header->nsamples );

        /* normalize the narrow band series
           this needs to be reviewed by an expert */
        norm = ((float)(fmaxindex - fminindex + 1))/( (float) header->nsamples);
        
        for(i = 0; i < (fmaxindex - fminindex + 1); i++)
        {
                (*dataF)[i] *= norm * mysteryFactor;
        }

        /* for the new narrow band series set the new first index
           and the new number of samples */
        header -> first = fminindex;
        header -> nsamples = fmaxindex -fminindex + 1;

        /* close the wide-band file and exit */
        fclose( fp );
        return 0;
}

/* byte swap nobj objects of size size stored in array ptr */
/* this converts between big- and little-endian objects */
/* directly from Jolien Creighton */
void byte_swap( void *ptr, size_t size, size_t nobj )
{
        while ( nobj-- > 0 )
        {
                size_t byte;
                for ( byte = 0; byte < size / 2; ++byte )
                {
                        char *b1 = ptr + byte;
                        char *b2 = ptr + size - byte - 1;
                        char tmp = *b1;
                        *b1 = *b2;
                        *b2 = tmp;
                }
                ptr += size;
        }
        return;
}


/*
Write SFT frequency series as a frame to a frame file.
With thanks to Jolien Creighton.
*/
int add_frame_to_frame_file(struct sftheader *header, float *dataF, int frnum, FrFile *outputFrameFile) 
{
        struct FrameH     *frame;       /* frame structure being created and populated */
        struct FrVect     *vect;        /* vector for series data */
        struct FrProcData *proc;        /* so-called processed data structure */
        char comment[] = "Generated from SFT data by narrowBandExtract";
        char hertz[] = "Hz";
        char units[] = "s";
        int ret;

        /* create a new frame */
        frame = FrameHNew( "LIGO" );
        if ( ! frame )
        {
            fprintf( stderr, "Could not create frame in frame file\n" );
            return 1;
        }

        /* -1 reserved for simulated data...see frame spec */
        frame->run    = -1;

        /* simple number to number frames inside this multi-frame frame file */
        frame->frame  = frnum;

        /* frame start time in GPS seconds */
        frame->GTimeS = header->gpssec;

        /* frame start time residual nanoseconds */
        frame->GTimeN = header->gpsnan;

        /* frame 'length' is the time-base of the SFT */
        frame->dt     = header->tbase;

        /* 
        name of vector is "SFT"
        type is complex 8
        number of samples is taken from SFT header
        step size is the size between frequency samples or 1/tbase
        step unit is Hertz
        the units of the series is seconds
        */
        vect = FrVectNew1D( "SFT", FR_VECT_8C, header->nsamples, 1.0 / header->tbase, hertz, units );
        if ( ! vect )
        {
            fprintf( stderr, "Could not create vector for frame file\n");
            FrameFree( frame );
            return 1;
        }

        /* the first frequency in the series */
        vect->startX[0] = header->first / header->tbase;

        /* create a processed data structure */ 
        proc = calloc( 1, sizeof( *proc ) );
        if ( ! proc )
        {
            fprintf( stderr, "Memory allocation error\n" );
            FrVectFree( vect );
            FrameFree( frame );
            return 1;
        }
        proc->classe    = FrProcDataDef();
        proc->type      = 2;                    /* 2 is for freq series */
        proc->subType   = 1;                    /* 1 is for DFT */
        proc->fShift    = 0;                    /* see frame spec */
        proc->tRange    = header -> tbase;      /* per agreement with Teviet */
        proc->data      = vect;                 
        proc->next      = frame->procData;
        frame->procData = proc;

        FrStrCpy( &proc->name, "SFT" );
        FrStrCpy( &proc->comment, comment );

        /* copy the series into the vector */
        memcpy( vect->dataF, dataF, 2 * header->nsamples * sizeof( float ) );

        /* write to frame file */
        ret = FrameWrite( frame, outputFrameFile );
        if ( ret != 0)
        {
                fprintf(stderr, "Could not write frame to file\n");
                FrameFree( frame );
                return 1;
        }

        /* done with this series so free frame */
        FrameFree(frame);

        return 0;
}  

int add_sft_to_mergedSFT_file(struct sftheader *header, float *dataF, FILE *outputSFTFile) 
{
        struct sftheader myHeader;
        int ret;

        myHeader.endian = header -> endian;
        myHeader.gpssec = header -> gpssec;
        myHeader.gpsnan = header -> gpsnan;
        myHeader.tbase  = header -> tbase; 
        myHeader.first  = header -> first;
        myHeader.nsamples = header -> nsamples;

        /* write header to the file */
        ret = fwrite((const void *) &myHeader, sizeof(myHeader), 1, outputSFTFile);

        if (ret != 1)
        {
                fprintf(stderr, "Error writing header to merged SFT file\n");
                exit(1);
        }


        /* write data to the file */
        ret = fwrite((const void *) dataF, 2 * sizeof(float), header -> nsamples, outputSFTFile);

        if (ret != header -> nsamples)
        {
                fprintf(stderr, "Error writing data to merged SFT file\n");
                exit(1);
        }

        return 0;

}
