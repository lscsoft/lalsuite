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

// Use detector names from lal/packages/tools/include/LALDetectors.h :
//      LAL_TAMA_300_DETECTOR   =       0,
//      LAL_VIRGO_DETECTOR      =       1,
//      LAL_GEO_600_DETECTOR    =       2,
//      LAL_LHO_2K_DETECTOR     =       3,
//      LAL_LHO_4K_DETECTOR     =       4,
//      LAL_LLO_4K_DETECTOR     =       5,

/*
 *  File names for input and output
 */
char* frame_file[6] = { 0, 0, 0, 0, 0, 0};

char* xml_file[6] = { 0, 0, 0, 0, 0, 0 };
typedef const char* cp;
cp channel_name[6] = { "T1" , "V1" , "G1" , "H2" , "H1" , "L1" };

const char* output_file = "skymap.txt";

char* event_id = 0;

/*
 *  Resolution of output sky map
 */
int ra_res = 512;
int dec_res = 256;

/*
 * Number of Detectors Used
 */
int numObs;

/*
 * Network Info...
 */
double *wSw;
int *detectors;

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
 *  ... etc ...
 *  x[ MAXOBS + 0 ] = imag(z_H)
 *  x[ MAXOBS + 1 ] = imag(z_L)
 *  x[ MAXOBS + 2 ] = imag(z_V)
 *  ... etc ...
 */
double** x;

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
  fprintf( stderr , "#INSIDE MAIN!!!\n");//<<-----------------FIXME
  int c;
  
  numObs = 0;
 
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
	  frame_file[LAL_LHO_4K_DETECTOR] = optarg;
	  break;
	case 'l':
	  frame_file[LAL_LLO_4K_DETECTOR] = optarg;
	  break;
	case 'v':
	  frame_file[LAL_VIRGO_DETECTOR] = optarg;
	  break;
	case 'i':
	  frame_file[LAL_LHO_2K_DETECTOR] = optarg;
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
	  xml_file[LAL_LHO_4K_DETECTOR] = optarg;
	  break;
	case 's':
	  xml_file[LAL_LLO_4K_DETECTOR] = optarg;
	  break;
	case 'r':
	  xml_file[LAL_VIRGO_DETECTOR] = optarg;
	  break;
	case 'q':
	  xml_file[LAL_LHO_2K_DETECTOR] = optarg;
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
        }// end switch
    }// end while
  
  if (optind < argc)
    {
      fprintf(stderr, "error: extraneous command line argument ");
      while (optind < argc)
        {
	  fprintf(stderr, "%s\n", argv[optind++]);
        }
      exit(1);
    }// end if
  
  /* support "none" arguments */
  {
    int i;
    for( i = 0 ; i < 6 ; ++i )
    {
      if( frame_file[i] && !strcmp("none", frame_file[i] ) ) { frame_file[i] = 0; }
      if( xml_file[i] && !strcmp("none" , xml_file[i] ) ) { xml_file[i] = 0; }
      if( ( frame_file[i] && !xml_file[i] ) || (xml_file[i] && !frame_file[i] ) )
        {
          fprintf( stderr , "error: Suppy matching pairs of frame & XML files\n");
          exit(1);
        }//end if
      if( xml_file[i] && frame_file[i] )
           ++numObs; 
    }// end i for
  }// end code block
  
  if (frequency <= 0)
    {
      fprintf(stderr, "error: Supply positive integer Hertz --sample-rate\n");
      exit(1);
    }
  
  fprintf( stderr , "#Cleared Sanity Check!!!\n");//<<-----------------FIXME
  
  /* figure out how many observatories and allocate memory for x and wSw */
  {
    int i,k=0;

    wSw = malloc( sizeof(double) * numObs ); 
    detectors = malloc( sizeof(int) * numObs );    
    
    for( i = 0 ; i < 6 ; ++i )
        if( frame_file[i] )
            detectors[k++] = i;

    x = malloc( sizeof( double* ) * numObs * 2 );
    
  }// end code block

  /* examine one second of data around the injection */
  samples = frequency;
  
  /*
   * Load metadata from XML files
   */
  {
    int i;
    for( i = 0 ; i < numObs ; ++i )
      {
        load_metadata( xml_file[ detectors[i] ] , i );
        load_data( i , frame_file[ detectors[i] ] , channel_name[ detectors[i] ] ); 
      }// end i for
   }
  
  fprintf( stderr , "#Entering Analsis...\n");//<<-----------------FIXME
  
  /*
   * Analyze the data and save the skymap
   */
  analyze();

  fprintf( stderr , "#Fully Analyzed... destorying last few arrays...\n");

  return 0;
}// end main

void load_metadata(char* file, int slot )
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
      wSw[slot] = a->sigmasq;
      greenwich = fmod(XLALGreenwichMeanSiderealTime(&(a->end_time)), LAL_TWOPI);
    }//end file if
}// end load_metadata

void load_data(int slot, const char* file, const char* initial)
{
  if (file)
    {
      /* 
       *  Read the frame file here
       */
      
      FrStream *stream = NULL;
      COMPLEX8TimeSeries H1series;
      int i;
      
      sprintf(H1series.name,"%s:CBC-CData_%s", initial, event_id);
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
      x[slot         ] = (double*) malloc(samples * sizeof(double));
      /* complex, or waveform two */
      x[slot + numObs] = (double*) malloc(samples * sizeof(double));     
      
      for (i = 0; i != samples; ++i)
        {
	  x[slot         ][i] = H1series.data->data[i].re;
	  x[slot + numObs][i] = H1series.data->data[i].im;
        }// end i for
    }// end file if
}// end load_data

#define NSIGMA 3

void analyze(void)
{
  int i,j;
   
  XLALSkymap2PlanType* plan;
  XLALSkymap2SphericalPolarType *directions;
  XLALSkymap2DirectionPropertiesType *properties;   
  XLALSkymap2KernelType *kernels;
  
  double* logPosteriors;
  double** xSw = malloc( sizeof( double* ) * numObs );
  double** xSw2 = malloc( sizeof( double* ) * numObs ); 
  
  /* 
   *  The characteristic size of the signal 
   */

  /* Validate the input */    
  for (i = 0; i != numObs; ++i)
    {
      if (x[i])
        {
          for (j = 0; j != samples; ++j)
            {
              if (!((x[i][j] > log(0)) && (x[i][j] < -log(0))))
	        {
	          fprintf(stderr, "error: x[%d][%d] = %e\n", i, j, x[i][j]);
	          exit(1);
	        }// end if
	      if (!((x[i + numObs][j] > log(0)) && (x[i + numObs][j] < -log(0))))
	        {
	          fprintf(stderr, "error: x[%d][%d] = %e\n", i + numObs, j, x[i + numObs][j]);
	          exit(1);
	        }//end if
	    }//end j for
	}// end if
    }//end i for
  
  fprintf( stderr , "#Input Validated...\n");//<<-----------------FIXME
  
  double log_pglitch;
  
  {
    
    /*
     * Compute the glitch hypotheses
     */
    double log_pglitches[numObs];
    
    for (i = 0; i != numObs; ++i)
      {
	log_pglitches[i] = log(0.0);
	if (x[i])
	  {
	    for (j = 0; j != samples; ++j)
	      {
		double log_p = -1.0 * log( 1 + pow( 1 , 2 )) + 
		.5 * pow( 1 + pow( 1 , -2 ) , -1 ) * ( pow( x[i][j] , 2 ) + pow( x[i+2][j] , 2 ) );
		log_pglitches[i] = XLALSkymapLogSumExp( log_pglitches[i] , log_p );
	      }// end j for
	    log_pglitches[i] -= log( samples );
	  }// end x if
      }// end i for
    
    log_pglitch = 0.0;
    for (i = 0; i != numObs; ++i)
      {
	if (x[i])
	  {
	    log_pglitch += log_pglitches[i];
	  }//end if
      }//end i for
  }// end code block
  
  fprintf( stderr , "#Glitch Hypothesis Found...\n");//<<-----------------FIXME
  
  /*   
   *  the sky tiles implied by the frequency) 
   */
  plan = malloc(sizeof(*plan));
  XLALSkymap2PlanConstruct( frequency , numObs , detectors , plan ); 
  
  /*
   *  Directions assigned for each pixel, using sine proj.
   */
  {
    directions = malloc(sizeof(*directions) * dec_res * ra_res );
    for ( i = 0 ; i != dec_res ; ++i )
      {
	for ( j = 0 ; j != ra_res ; ++j )
	  {
	    directions[i*ra_res + j][0] = acos((i - dec_res/2.0 + 0.5)/dec_res*2);
	    directions[i*ra_res + j][1] = ( j + .5 )/ra_res * 2 * LAL_PI;
	  }//end j for
      }// end i for
  }// end code block
  
  fprintf( stderr , "#Pixel Direction Assigned...\n");//<<-----------------FIXME
  
  /*
   *  Properties Constructed
   */    
  properties = malloc(sizeof(*properties) * dec_res * ra_res );
  for (i = 0 ; i != ra_res*dec_res ; ++i)
    {
      XLALSkymap2DirectionPropertiesConstruct(
 		                              plan, 
				              directions + i,
				              properties + i
				              );
    }// end i for
  
  
  /*
   *  Posterior probabilities found for 4 different signal amplitudes,
   *  at each pixel in the sky, over all arrival times...
   */
  {
    int n;
    {
      double maxwSw = 0;
      for( j = 0 ; j < numObs ; ++j )
	if( wSw[j] > maxwSw )
	  maxwSw = wSw[j];
      for( j = 0 ; j < numObs ; ++j )
	wSw[j] /= maxwSw;
    }// end code block
    
    for (i = 0; i != numObs; ++i)
    {
      xSw[i] = malloc(sizeof(*xSw[i]) * samples);
      xSw2[i] = malloc(sizeof(*xSw2[i]) * samples);
      for (j = 0; j != samples; ++j)
         xSw[i][j] = x[i][j] * sqrt( wSw[i] );
         xSw2[i][j] = x[i + numObs ][j] * sqrt( wSw[i] );
    }// end i for

    fprintf( stderr , "#Commencing kernel calculations...\n");//<<-----------------FIXME    

    kernels = malloc(sizeof(*kernels) * dec_res * ra_res);
    logPosteriors = malloc(sizeof(*logPosteriors) * ra_res * dec_res );
    for( j = 0 ; j < ra_res * dec_res ; ++j )
      logPosteriors[j] = log(0);
    double* buffer = malloc(sizeof(double) * samples);
    
    fprintf( stderr , "#Memory Allocated... btw's numObs == %d\n" , numObs ); //<<--------- FIXME
    
    {
      for( n = 0 ; n < 4 ; ++n ) /* over 4 waveForm amplitudes */
        {
  	  for (i = 0; i != ra_res * dec_res; ++i) /* over all pixels */
	    XLALSkymap2KernelConstruct( plan , properties + i, wSw, kernels + i );
          fprintf( stderr , "#samples = %d\n" , samples ); //<<--------------------- FIXME	
	  for (i = 0; i != ra_res * dec_res; ++i) /* over all pixels */
	    {
	      int t;
	      for( t = samples/4 ; t < 3 * samples / 4 ; ++t ) /* over all sample times */
	        {
                  double real , imag;
		  XLALSkymap2Apply( plan , 
                                    properties + i , 
                                    kernels + i , 
                                    xSw  , 
                                    t , 
                                    &real );       
		  XLALSkymap2Apply( plan , 
                                    properties + i , 
                                    kernels + i , 
                                    xSw2 , 
                                    t , 
                                    &imag );
		  buffer[t] = real + imag;  
	        }// end t for
	      logPosteriors[i] = XLALSkymapLogSumExp( logPosteriors[i] , 
	                                              XLALSkymapLogTotalExp( buffer + samples / 4 , buffer + 3 * samples / 4 ) - log(4) - log(samples/2) );
	    }//end i for
          fprintf( stderr , "#N = %d\n" , n );	
	  /*
	   *  For the next loop, look for a signal sqrt(10) larger...
	   */
	  for( i = 0 ; i < numObs ; ++i )
	  {
	    wSw[i] *= 10.0;
	    for( j = 0 ; j != samples; ++j )
	    {
	       xSw[i][j]  *= sqrt( 10.0 ); 
               xSw2[i][j] *= sqrt( 10.0 );
            }// end j for
	  }//end i for
        }// end n for
    
    }// end code block

    free(buffer);
  }// end posterior computing code block
 
  fprintf( stderr , "#Posterior Found, writing...\n");//<<-----------------FIXME 
 
  /*  Allocate a chunk of memory to hold the sky map in the internal 
   *  timing format
   */
  
  /* validate the output */
  for (i = 0; i != ra_res * dec_res; ++i)
    {
      if( isnan(logPosteriors[i]) )
	{
	  fprintf(stderr, "BAD DIRECTION :( THETA: %f PHI: %f i:%d\n", directions[i][0],directions[i][1],i);
	  exit(1);
	}
    }
  
  int iMode = 0;
  for( i = 0 ; i < ra_res * dec_res ; ++i )
    {
      if(logPosteriors[i] > logPosteriors[iMode])
	iMode = i;
    }
  
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
  
  for ( i = 0; i != ra_res * dec_res ; ++i )
    {
      double phi,ra,dec,longitude;
      phi = directions[i][1];
      ra  = fmod(phi + greenwich, LAL_TWOPI);
      longitude = fmod( phi , LAL_TWOPI);
      while (ra < 0)
     	ra += LAL_TWOPI;
      while ( longitude < 0 )
        longitude += LAL_TWOPI;
      while (ra >= LAL_TWOPI)
	ra -= LAL_TWOPI;
      while ( longitude >= LAL_TWOPI )
        longitude -= LAL_TWOPI;
      dec = LAL_PI_2 - directions[i][0];
      gzprintf(h, "%.10e %.10e %.10e %.10e %.10e\n", ra, dec, (logPosteriors[i]) , longitude , dec );
    }// end i for
  gzclose(h);
 
  fprintf( stderr , "#Written to file... writing summary...\n");//<<-----------------FIXME
   
  /*
   * Summary of Results printed to log file
   */
  {
    double log_totalProb = XLALSkymapLogTotalExp( logPosteriors , logPosteriors + ra_res*dec_res );
    log_totalProb -= log( ra_res * dec_res );
    double raMax = fmod( directions[iMode][0] + greenwich , LAL_TWOPI );
    double longMax = fmod( directions[iMode][0] , LAL_TWOPI );
    double decMax = directions[iMode][1];    
    printf("#EVENT ID \t\t FOUND RA \t\t FOUND DEC \t\t TOTAL PROB \t\t GLITCH PROB \t\t LONG \t\t LAT\n");
    printf("%s \t\t %f \t\t %f \t\t %e \t\t %e \t\t %f \t\t %f\n" , event_id , raMax , decMax , log_totalProb , log_pglitch , longMax , decMax ); 
  }// end summary print block       
  
  /*
   *  Free the data
   */
  free(plan);
  for (i = 0; i != 2*numObs; ++i)
    free(x[i]);
  free(logPosteriors);
  for(i = 0; i != numObs; ++i)
  {
     free(xSw2[i]);
     free(xSw[i]);    
  }// end i for
  free(kernels);
  free(properties);
  free(directions);    
  
}// end analyze(void)

