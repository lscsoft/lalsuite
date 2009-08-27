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
#define NSIGMA 3

typedef double XLALSkymap2SphericalPolarType[2];
        
typedef struct 
{
    XLALSkymap2SphericalPolarType *directions;
    double* logPosteriors;
    double total_logPosterior;
    int count;
} SkyMapProperties;

/* 
 * Use detector names from lal/packages/tools/include/LALDetectors.h :
 *
 *      LAL_TAMA_300_DETECTOR   =       0,
 *      LAL_VIRGO_DETECTOR      =       1,
 *      LAL_GEO_600_DETECTOR    =       2,
 *      LAL_LHO_2K_DETECTOR     =       3,
 *      LAL_LHO_4K_DETECTOR     =       4,
 *      LAL_LLO_4K_DETECTOR     =       5,
 *
 *-------------------------------------------------------
 *
 *  Array of pointers to the matched-filtered data z[t]
 * 
 *       x[0] = real(z_H)
 *       x[1] = real(z_L)
 *       x[2] = real(z_V)
 *          ... etc ...
 *       x[ MAXOBS + 0 ] = imag(z_H)
 *       x[ MAXOBS + 1 ] = imag(z_L)
 *       x[ MAXOBS + 2 ] = imag(z_V)
 *          ... etc ...
 */
typedef struct 
{
    int N;
    int* detectors;
    double* wSw;
    double** x;    
} NetworkProperties;

/*
 *  File names for input and output
 */
char* frame_file[6] = { 0, 0, 0, 0, 0, 0};
char* xml_file[6] = { 0, 0, 0, 0, 0, 0 };
typedef const char* cp;
cp channel_name[6] = { "T1" , "V1" , "G1" , "H2" , "H1" , "L1" };
const char* output_file = "skymap.txt";

/*
 *  Global Variables...
 */
char* event_id = 0;    // Event ID
int ra_res = 512;      // resolution of output skyMap
int dec_res = 256;     // ""
int frequency = 0;     // sampling freq. of analysis (det's sky tilting)
int samples = 512;     // Number of samples to analyze
double greenwich = 0;  // Orientation of Earth at event time

/*
 * Functions
 */
void  load_metadata( NetworkProperties* network,
                     int slot);
void      load_data( NetworkProperties* network,
                     int slot );
void writeSkymapOut( SkyMapProperties* skyMap );
void        analyze( NetworkProperties* network, 
                     SkyMapProperties* skyMap );

/*
 *  MAIN
 */
int main(int argc, char** argv)
{
  fprintf( stderr , "#INSIDE MAIN!!!\n");//<<-----------------FIXME
  int c,numObs = 0;
 
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
        {
           ++numObs;
           printf("NUM OBS: %d\n" , numObs ); 
        }// end if
    }// end i for
  }// end code block
  
  if (frequency <= 0)
    {
      fprintf(stderr, "error: Supply positive integer Hertz --sample-rate\n");
      exit(1);
    }
  
  /* examine one second of data around the injection */
  samples = frequency;

  
  fprintf( stderr , "#Cleared Sanity Check!!!\n");//<<-----------------FIXME
  
  /* 
   *  We create numObs + 1 networks... one for the whole system, and one network
   *  for each single detector:
   */
  NetworkProperties* network;
  NetworkProperties* singleObs;

  {
    int i,k=0;
    
    network = malloc(sizeof(*network));
    network->N = numObs; 
    
    network->wSw = malloc( sizeof(double) * numObs ); 
    network->detectors = malloc( sizeof(int) * numObs );    
    
    for( i = 0 ; i < 6 ; ++i )
        if( frame_file[i] )
            network->detectors[k++] = i;

    network->x = malloc( sizeof( double* ) * numObs * 2 );
  }// end code block

   
 
  /*
   * Load metadata from XML files
   */
  {
    int i;
    for( i = 0 ; i < numObs ; ++i )
      {
        load_metadata(network , i);
        load_data(network , i);
      }// end i for
   }

  {
      int j;
      double maxwSw = 0;
      for( j = 0 ; j < network->N ; ++j )
        if( network->wSw[j] > maxwSw )
          maxwSw = network->wSw[j];
      for( j = 0 ; j < network->N ; ++j )
          network->wSw[j] /= maxwSw;
  }// end code block
    
  if( numObs > 1 )
    {
      singleObs = malloc( sizeof(*network) * numObs );
      
      int i;
      for( i = 0 ; i != numObs ; ++i )
        {
          singleObs[i].N = 1;
          singleObs[i].wSw = malloc(sizeof(double));
          singleObs[i].wSw[0] = network->wSw[i];

          singleObs[i].detectors = malloc(sizeof(int));
          singleObs[i].detectors[0] = network->detectors[i];

          singleObs[i].x = malloc(sizeof(double*) * 2 );
          singleObs[i].x[0] = network->x[i];
          singleObs[i].x[1] = network->x[i + numObs]; 
        }//end for
    }//end if 

  fprintf( stderr , "#Entering Analsis...\n");//<<-----------------FIXME
  
  /*
   *  Analyze the data and save the skymap for the N-Detector Case
   */
  SkyMapProperties* skyMap;
  skyMap = malloc(sizeof(*skyMap));
  skyMap->count = ra_res * dec_res;
  skyMap->directions = malloc(sizeof(*(skyMap->directions))*skyMap->count);
  skyMap->logPosteriors = malloc(sizeof(*(skyMap->logPosteriors)) * skyMap->count );

  /*
   *  Analyze the data for the N single-detector cases, and store the net
   *  log posteriors...
   */
  double* total_logPosteriors;
  if( numObs > 1 )
    {
      total_logPosteriors = malloc(sizeof(double) * numObs );
      int i;
      for( i = 0 ; i != numObs ; ++i )
        {
          analyze(  singleObs + i , skyMap );
          total_logPosteriors[i] = skyMap->total_logPosterior;
          fprintf(stderr,"total_logPosterior[%d] = %e\n", i , total_logPosteriors[i]);
        }// end i for
    }//end if
  
  /*
   *   Analyze the full network, write out the skymap, and 
   *   then print the final summary...
   */

  analyze( network , skyMap );
  writeSkymapOut( skyMap );

  {
    int i,iMode = 0;
    for( i = 0 ; i < skyMap->count ; ++i )
      {
        if(skyMap->logPosteriors[i] > skyMap->logPosteriors[iMode])
          iMode = i;
      }

    double raMax = fmod( skyMap->directions[iMode][0] + greenwich , LAL_TWOPI );
    double longMax = fmod( skyMap->directions[iMode][0] , LAL_TWOPI );
    double decMax = skyMap->directions[iMode][1];

    printf("#EVENT ID \t FOUND RA \t FOUND DEC \t");
    printf(" TOTAL PROB \t LONG \t LAT\t");
    if( numObs > 1 )
        for( i = 0 ; i != numObs ; ++i )
            printf(" Detector_%d \t " , i );
    printf("\n");
    printf("%s  %f \t %f \t %e \t %f \t %f\t" ,
           event_id , raMax , decMax , skyMap->total_logPosterior ,
           longMax , decMax );
    if( numObs > 1 )
        for( i = 0 ; i != numObs ; ++i )
            printf(" %e \t",total_logPosteriors[i]);
    printf("\n");
  }// end summary block 
 
   
  /*
   *  Raw data free'd
   */
  {
    int i;
    for (i = 0; i != 2*numObs; ++i)
      free(network->x[i]); 
  }
  
  
  /*
   *  posterior data free'd
   */
  free( skyMap->logPosteriors );
  free( skyMap->directions );
  free( singleObs );

  return 0;
}// end main

void load_metadata(NetworkProperties* network , int slot)
{
  char* file = xml_file[ network->detectors[slot] ];
  
  if (file)
    {
      SnglInspiralTable* a = 0;
      LALSnglInspiralTableFromLIGOLw(&a, file, 0, 1);
      if (!a)
        {
	  fprintf(stderr, "error: failed to read single inspiral table from file %s\n", file);
	  exit(1);
        }
      network->wSw[slot] = a->sigmasq;
      greenwich = fmod(XLALGreenwichMeanSiderealTime(&(a->end_time)), LAL_TWOPI);
    }//end file if
}// end load_metadata

void load_data(NetworkProperties* network , int slot )
{
   const char* file = frame_file[ network->detectors[slot] ];
   const char* initial = channel_name[ network->detectors[slot] ];
  
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
       * real/complex or waveform one/two  
       */
      network->x[slot              ] = (double*) malloc(samples * sizeof(double));
      network->x[slot + network->N ] = (double*) malloc(samples * sizeof(double));     

      for (i = 0; i != samples; ++i)
        {
	  network->x[slot              ][i] = H1series.data->data[i].re;
	  network->x[slot + network->N ][i] = H1series.data->data[i].im;
        }// end i for
    }// end file if
}// end load_data

void writeSkymapOut( SkyMapProperties* skyMap )
{
    int i;
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

    for ( i = 0; i != skyMap->count ; ++i )
      {
        double phi,ra,dec,longitude;
        phi = skyMap->directions[i][1];
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
        dec = LAL_PI_2 - skyMap->directions[i][0];
        gzprintf( h, "%.10e %.10e %.10e %.10e %.10e\n", 
                  ra, dec, (skyMap->logPosteriors[i]) , longitude , dec );
      }// end i for
    gzclose(h);
}// end writeSkyMapOut



/*
 *  ANALYZE
 */
void analyze( NetworkProperties* network , SkyMapProperties* skyMap )
{
  int i,j;
  
  XLALSkymap2PlanType* plan;
  XLALSkymap2DirectionPropertiesType *properties;   
  XLALSkymap2KernelType *kernels;
  
  double** xSw = malloc( sizeof( double* ) * network->N );
  double** xSw2 = malloc( sizeof( double* ) * network->N ); 
  
  /* 
   *  The characteristic size of the signal 
   */

  /* Validate the input */    
  for (i = 0; i != network->N; ++i)
    {
      if (network->x[i])
        {
          for (j = 0; j != samples; ++j)
            {
              if (!( ( network->x[i][j] >  log(0) )  && 
                     ( network->x[i][j] < -log(0) )      ))
	        {
	          fprintf(stderr, "error: x[%d][%d] = %e\n", i, j, network->x[i][j]);
	          exit(1);
	        }// end if
	      if (!(  (network->x[i + network->N][j] > log(0)) && 
                      (network->x[i + network->N][j] < -log(0))   ))
	        {
	          fprintf( stderr , "error: x[%d][%d] = %e\n", 
                           i + network->N, j, network->x[i + network->N][j]);
	          exit(1);
	        }//end if
	    }//end j for
	}// end if
       else
         {
           fprintf(stderr , "Data not supplied for detector %d\n" , i );
           exit(1);
         }// end else
    }//end i for
  
  fprintf( stderr , "#Input Validated...\n");//<<-----------------FIXME
    
    /*
     * Compute the glitch hypotheses  //<----------- TEMP REMOVED FIXME
     *
  {
    double log_pglitches[network->N];
    
    for (i = 0; i != network->N; ++i)
      {
	log_pglitches[i] = log(0.0);
	if (network->x[i])
	  {
	    for (j = 0; j != samples; ++j)
	      {
		double log_p = -1.0 * log( 1 + pow( 1 , 2 )) + 
		.5 * pow(1 + pow(1,-2),-1) * (pow(x[i][j],2) + pow(x[i+network->N][j],2));
		log_pglitches[i] = XLALSkymapLogSumExp( log_pglitches[i] , log_p );
	      }// end j for
	    log_pglitches[i] -= log( samples );
	  }// end x if
      }// end i for
    
    skyMap->log_pGlitch = 0.0;
    for (i = 0; i != network->N ; ++i)
      {
	if (x[i])
	  {
	    skyMap->log_pGlitch += log_pglitches[i];
	  }//end if
      }//end i for
  }// end code block
  
  ----------- FIXME END TEMP REMOVED FIXME 
 
  */
  
  /*   
   *  the sky tiles implied by the frequency) 
   */
  plan = malloc(sizeof(*plan));
  XLALSkymap2PlanConstruct( frequency , network->N , network->detectors , plan ); 
  
  /*
   *  Directions assigned for each pixel, using sine proj.
   */
  {
    for ( i = 0 ; i != dec_res ; ++i )
      {
	for ( j = 0 ; j != ra_res ; ++j )
	  {
	    skyMap->directions[i*ra_res + j][0] = acos((i - dec_res/2.0 + 0.5)/dec_res*2);
	    skyMap->directions[i*ra_res + j][1] = (j + .5)/ra_res*2*LAL_PI;
	  }//end j for
      }// end i for
  }// end code block
  
  fprintf( stderr , "#Pixel Direction Assigned...\n");//<<-----------------FIXME
  
  /*
   *  Properties Constructed
   */    
  properties = malloc(sizeof(*properties) * skyMap->count );
  for (i = 0 ; i != skyMap->count ; ++i)
    {
      XLALSkymap2DirectionPropertiesConstruct(
 		                              plan, 
				              skyMap->directions + i,
				              properties + i
				              );
    }// end i for
  
  fprintf(stderr,"#Properties Constructed, N = %d \n" , network->N); 
  
  /*
   *  Posterior probabilities found for 4 different signal amplitudes,
   *  at each pixel in the sky, over all arrival times...
   */
  {
    for (i = 0; i != network->N; ++i)
    {
      xSw[i] = malloc(sizeof(double) * samples);
      xSw2[i] = malloc(sizeof(double) * samples);
      for (j = 0; j != samples; ++j)
        {
          xSw[i][j]  = network->x[i][j]               * sqrt( network->wSw[i] );
          xSw2[i][j] = network->x[i + network->N ][j] * sqrt( network->wSw[i] );
        }// end j for
    }// end i for

    fprintf( stderr , "#Commencing kernel calculations...\n");//<<-----------------FIXME    

    kernels = malloc( sizeof(*kernels) * skyMap->count );
    for( j = 0 ; j != skyMap->count ; ++j )
      skyMap->logPosteriors[j] = log(0);
    double* buffer = malloc(sizeof(double) * samples);
    
    fprintf( stderr , "#Memory Allocated... btw's numObs == %d\n" , network->N );//<<--- FIXME
    
    {
      int k;
      double* wSw;
      wSw = malloc( sizeof(double) * network->N );
      for( i = 0 ; i < network->N ; ++i )
        {
          wSw[i] = network->wSw[i];
        }// end i for

      for( k = 0 ; k < 4 ; ++k ) /* over 4 waveForm amplitudes */
        {
  	  for (i = 0; i != skyMap->count; ++i) /* over all pixels */
	    XLALSkymap2KernelConstruct( plan , properties + i, wSw, kernels + i );
          fprintf( stderr , "#samples = %d\n#count= %d\n" , samples , skyMap->count); //<< FIXME	
          
	  for (i = 0; i != skyMap->count; ++i) /* over all pixels */
	    {
	      int t;
              //fprintf(stderr, "%d: i = %d\n", __LINE__, i);
	      for( t = samples/4 ; t < 3 * samples / 4 ; ++t ) /* over all sample times */
	        {
                  double real , imag;
                  //fprintf(stderr, "%d: i = %d, t = %d\n", __LINE__, i, t);
		  XLALSkymap2Apply( plan , 
                                    properties + i , 
                                    kernels + i , 
                                    xSw  , 
                                    t , 
                                    &real );       
                  //fprintf(stderr, "%d: i = %d, t = %d\n", __LINE__, i, t);
		  XLALSkymap2Apply( plan , 
                                    properties + i , 
                                    kernels + i , 
                                    xSw2 , 
                                    t , 
                                    &imag );
		  buffer[t] = real + imag;  
	        }// end t for
                //if( isnan(skyMap->logPosteriors[i] ) )
                //  {
                //    fprintf( stderr , "NAN @ i = %d\n" , i );
                //    exit(1);
                //  }//end if
	      skyMap->logPosteriors[i] = XLALSkymapLogSumExp( skyMap->logPosteriors[i] , 
	                                 XLALSkymapLogTotalExp( buffer + samples / 4 , 
                                         buffer + 3*samples/4) - log(4) - log(samples/2) );
            }//end i for
          fprintf( stderr , "#N = %d\n" , k );	
	  /*
	   *  For the next loop, look for a signal sqrt(10) larger...
	   */
	  for( i = 0 ; i < network->N ; ++i )
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
 
  fprintf( stderr , "#Posterior Found, writing...\n");//<<--------FIXME 
 
  /* validate the output */
  for (i = 0; i != skyMap->count ; ++i)
    {
      if( isnan(skyMap->logPosteriors[i]) )
	{
	  fprintf(stderr, "BAD DIRECTION :( THETA: %f PHI: %f i:%d\n", 
                  skyMap->directions[i][0],skyMap->directions[i][1],i);
	  exit(1);
	}//end if
    }// end i for
  {
    double log_totalProb = XLALSkymapLogTotalExp( skyMap->logPosteriors ,
                                                  skyMap->logPosteriors + skyMap->count );
    skyMap->total_logPosterior = log_totalProb - log( skyMap->count );
    fprintf( stderr , "#TOTAL PROB: %e\n" , skyMap->total_logPosterior );
  }
  /*
   *  Free the data
   */
  free(plan);
  for(i = 0; i != network->N; ++i)
  {
     free(xSw2[i]);
     free(xSw[i]);    
  }// end i for
  free(kernels);
  free(properties);
}// end analyze(void)

