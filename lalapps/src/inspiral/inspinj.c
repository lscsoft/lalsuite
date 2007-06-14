/*----------------------------------------------------------------------- 
 * 
 * File Name: inspinj.c
 *
 * Author: Brown, D. A., Creighton, J. D. E. and Dietz A. 
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#include <ctype.h>
#include <getopt.h>
#include <lalapps.h>
#include <lal/Date.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/LIGOLwXMLRead.h>
#include <lal/LIGOLwXML.h>
#include <lal/Random.h>
#include <lal/AVFactories.h>
#include <lal/InspiralInjectionParams.h>

RCSID( "$Id$" );

#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define PROGRAM_NAME "inspinj"

#define ADD_PROCESS_PARAM( pptype, format, ppvalue ) \
  this_proc_param = this_proc_param->next = (ProcessParamsTable *) \
calloc( 1, sizeof(ProcessParamsTable) ); \
LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", \
    PROGRAM_NAME ); \
LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "--%s", \
    long_options[option_index].name ); \
LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "%s", pptype ); \
LALSnprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, format, ppvalue );


/* 
 *  *********************************
 *  Definition of the prototypes 
 *  *********************************
 */
ProcessParamsTable *next_process_param( const char *name, const char *type,
					const char *fmt, ... );
void printUsage(void);
void read_mass_data( char* filename );
void read_source_data( char* filename );
void drawDistanceFromSource( SimInspiralTable* table );
void drawLocationFromSource( SimInspiralTable* table );
void drawLocationFromExttrig( SimInspiralTable* table );
void drawMassFromSource( SimInspiralTable* table );

/* 
 *  *************************************
 *  Defining of the used global variables
 *  *************************************
 */

DistanceDistribution    dDistr;
SkyLocationDistribution lDistr;
MassDistribution        mDistr;
InclDistribution        iDistr;

SimInspiralTable *simTable;

char *massFileName = NULL;
char *sourceFileName = NULL;
char *outputFileName = NULL;
char *exttrigFileName = NULL;

float mwLuminosity = -1;
REAL4 dmin= -1;
REAL4 dmax= -1;
REAL4 minMass1=-1;
REAL4 maxMass1=-1;
REAL4 minMass2=-1;
REAL4 maxMass2=-1;
REAL4 minMtotal=-1;
REAL4 maxMtotal=-1;
REAL4 meanMass1=-1.0;
REAL4 meanMass2=-1.0;
REAL4 massStdev1=-1.0;
REAL4 massStdev2=-1.0;
REAL4 inclStd=-1.0;
int spinInjections=0;
REAL4 minSpin1=-1.0;
REAL4 maxSpin1=-1.0;
REAL4 minSpin2=-1.0;
REAL4 maxSpin2=-1.0;

static LALStatus status;
static RandomParams* randParams=NULL;
INT4 numExtTriggers = 0;
ExtTriggerTable   *exttrigHead = NULL;

int num_source;
struct {
  char   name[20];
  REAL8 ra;
  REAL8 dec;
  REAL8 dist;
  REAL8 lum;
  REAL8 fudge;
} *source_data;

REAL8* fracVec  =NULL;
REAL8* ratioVec = NULL;
REAL8 norm=0;

int num_mass;
struct {
  REAL8 mass1;
  REAL8 mass2;
} *mass_data;

/* 
 *  *********************************
 *  Implementation of the code pieces  
 *  *********************************
 */

/*
 *
 * code to step forward in the process table
 *
 */
ProcessParamsTable *next_process_param( const char *name, const char *type,
    const char *fmt, ... )
{
  ProcessParamsTable *pp;
  va_list ap;
  pp = calloc( 1, sizeof( *pp ) );
  if ( ! pp )
  {
    perror( "next_process_param" );
    exit( 1 );
  }
  strncpy( pp->program, PROGRAM_NAME, LIGOMETA_PROGRAM_MAX );
  LALSnprintf( pp->param, LIGOMETA_PARAM_MAX, "--%s", name );
  strncpy( pp->type, type, LIGOMETA_TYPE_MAX );
  va_start( ap, fmt );
  LALVsnprintf( pp->value, LIGOMETA_VALUE_MAX, fmt, ap );
  va_end( ap );
  return pp;
}

/*
 *
 * print-out of the usage
 *
 */
void printUsage(void)
{
  printf("lalapps_inspinj [options]\n" );
  printf("\nDefaults are shown in brackets\n\n" );
  printf("  --help                   display this message\n");
  printf("  --source-file FILE       read source parameters from FILE\n");
  printf("  --mass-file FILE         "
	 "read population mass parameters from FILE\n");
  printf("  --exttrig-file FILE      XML file containing external trigger\n");
  printf("  --f-lower FREQUENCY      lower cut-off frequency.\n");
  printf("  --gps-start-time TIME    "
	 "start injections at GPS time TIME\n");
  printf("  --gps-end-time TIME      "
	 "end injections at GPS time TIME \n");
  printf("  --time-step STEP         "
	 "space injections by ave of STEP sec (2630 / PI)\n");
  printf("  --time-interval TIME     "
	 "distribute injections in interval TIME (0)\n");
  printf("  --seed SEED              "
	 "seed random number generator with SEED (1)\n");
  printf("  --waveform NAME          set waveform type to NAME\n");
  printf("  --user-tag STRING        set the usertag to STRING\n");
  printf("  --enable-milkyway LUM    "
	 "enables Milky Way injections, set MW luminosity LUM\n");
  printf("  --disable-milkyway       disables Milky Way injections\n");
  printf("  --min-distance DMIN      "
	 "set the minimum distance to DMIN kpc\n");
  printf("  --max-distance DMAX      set the maximum distance to DMAX kpc\n");
  printf("  --d-distr DDISTR         "
	 "distribute injections uniformly over d\n");
  printf("                           (DDISTR = uniform), or over log10(d) \n");
  printf("                           (DDISTR = log10),\n");
  printf("                           or over volume (DDISTR = volume)\n");
  printf("                           or using source list (DDISTR = source)\n");
  printf("  --l-distr LDISTR         "
	 "sets the distribution for the source locations,\n");
  printf("                           "
	 "LDISTR=source uses the locations given in the \n");
  printf("                           "
         "source-file, LDISTR=exttrig uses external trigger file,\n");
  printf("                           "
         "LDISTR=random uses random locations\n");
  printf("  --i-distr INCDIST        "
	 "distribute inclinations uniformly ver arccos(i)\n");
  printf("                           "
	 "(INCDIST=uniform) or gaussian distributed in \n");
  printf("                           arccos(i)  (INCDIST=gaussian)\n");
  printf("  --inclStd INCLSTD        "
	 "sets the standard deviation for gaussian inclination distr.\n");
  printf("  --min-mass1 MIN          "
	 "set the minimum component mass to MIN\n");
  printf("  --max-mass1 MAX          "
	 "set the maximum component mass to MAX\n");
  printf("  --min-mass2 MIN          "
	 "set the minimum component mass2 to MIN (minMass1)\n");
  printf("  --max-mass2 MAX          "
	 "set the maximum component mass2 to MAX (maxMass1)\n"); 
  printf("  --min-mtotal MINTOTAL    "
	 "sets the minimum total mass to MINTOTAL\n");
  printf("  --max-mtotal MAXTOTAL    "
	 "sets the maximum total mass to MAXTOTAL\n");
  printf("  --m-distr MDISTR         distribute injections uniformly over\n");
  printf("                           "
	 "total mass (MDISTR = totalMass), or over mass1 and\n");
  printf("                           over mass2 (MDISTR = componentMass), \n");
  printf("                           or gaussian (MDISTR=gaussian), \n");
  printf("                           or log in comonent mass (MDISTR=log),\n");
  printf("                           or using mass file (MDISTR=source)\n");
  printf("  --mean-mass MASS         "
	 "set the mean value for both mass components if \n");
  printf("                           MDISTR=gaussian\n");
  printf("  --mean-mass2 MASS        "
	 "set the mean value for mass2 if MDISTR=gaussian\n");
  printf("  --stdev-mass MSTD        "
	 "set the standard deviation for both component\n");
  printf("                           masses if MDISTR=gaussian\n");
  printf("  --stdev-mass2 MSTD       "
	 "set the standard deviation for mass2 if \n");
  printf("                           MDISTR=gaussian\n");
  printf("  --disable-spin           disables spinning injections\n");
  printf("  --enable-spin            enables spinning injections\n");
  printf("  [--min-spin1 MIN]        set the minimum spin1 to MIN (0.0)\n");
  printf("  [--max-spin1 MAX]        set the maximum spin1 to MAX (0.0)\n");
  printf("  [--min-spin2 MIN]        set the minimum spin2 to MIN (0.0)\n");
  printf("  [--max-spin2 MAX]        set the maximum spin2 to MAX (0.0)\n");
  printf("  [--output NAME]          set the output filename \n");
  printf("\n");
}


/*
 *
 * functions to read source masses 
 *
 */

void 
read_mass_data( char* filename )
{
  /*const char *basename = massFileName ? massFileName : "BNSMasses.dat";*/
  char line[256];
  FILE   *fp;
  int n = 0;

  fp=fopen( filename, "r" );
  if ( ! fp )
  {
    perror( "read_mass_data" );
    fprintf( stderr, 
        "Error while trying to open file %s\n", 
        filename );
    exit( 1 );
  }

  /* count the number of lines in the file */
  num_mass=0;
  while ( fgets( line, sizeof( line ), fp ) )
    ++num_mass;

  /* alloc space for the data */
  mass_data = LALCalloc( num_mass, sizeof(*mass_data) );
  if ( !mass_data ) 
  {
    fprintf( stderr, "Allocation error for mass_data\n" );
    exit( 1 );
  }
  
  /* 'rewind' the file */
  rewind( fp );

  /* read the file finally */
  while ( fgets( line, sizeof( line ), fp ) )
  {
    sscanf( line, "%le %le", &mass_data[n].mass1, &mass_data[n].mass2 );
    n++;
  }

  /* close the file */
  fclose( fp );
}

/*
 *
 * functions to read source distribution
 *
 */

void 
read_source_data( char* filename )
{
  char line[256];
  FILE *fp;
  int i;

  fp = fopen (filename, "r" );
  if ( ! fp )
  {
    perror( "read_source_data" );
    fprintf( stderr, "Could not find file %s\n", filename );
    exit( 1 );
  }

  /* count the number of entries in this file */
  num_source = 0;
  while ( fgets( line, sizeof( line ), fp ) )
    if ( line[0] == '#' )
      continue;
    else 
      ++num_source;

  /* rewind the file */
  rewind( fp );

  /* allocate space */
  source_data = LALCalloc( num_source, sizeof( *source_data ) );
  if ( ! source_data )
  {
    fprintf( stderr, "Allocation error for source_data\n" );
    exit( 1 );
  }

  i = 0;
  while ( fgets( line, sizeof( line ), fp ) )
    if ( line[0] == '#' )
      continue;
    else
    {
      char ra_sgn, dec_sgn;
      REAL8 ra_h, ra_m, dec_d, dec_m;
      int c;

      c = sscanf( line, "%s %c%le:%le %c%le:%le %le %le %le",
          source_data[i].name, &ra_sgn, &ra_h, &ra_m, &dec_sgn, &dec_d, &dec_m,
          &source_data[i].dist, &source_data[i].lum, &source_data[i].fudge );
      if ( c != 10 )
      {
        fprintf( stderr, "error parsing source datafile %s\n", sourceFileName );
        exit( 1 );
      }

      /* by convention, overall sign is carried only on hours/degrees entry */
      source_data[i].ra  = ( ra_h + ra_m / 60.0 ) * LAL_PI / 12.0;
      source_data[i].dec = ( dec_d + dec_m / 60.0 ) * LAL_PI / 180.0;

      if ( ra_sgn == '-' )
        source_data[i].ra *= -1;
      if ( dec_sgn == '-' )
        source_data[i].dec *= -1;
      ++i;
    }

  /* close file */
  fclose( fp );


  /* generate ratio and fraction vectors */
  ratioVec = calloc( num_source, sizeof( REAL8 ) );
  fracVec  = calloc( num_source, sizeof( REAL8  ) );
  if ( !ratioVec || !fracVec )
  {
    fprintf( stderr, "Allocation error for ratioVec/fracVec\n" );
    exit( 1 );
  }
    
  /* MW luminosity might be zero */
  norm = mwLuminosity;
  
  /* calculate the fractions of the different sources */
  for ( i = 0; i < num_source; ++i )
    norm += ratioVec[i] = source_data[i].lum * source_data[i].fudge;
  fracVec[0] = ratioVec[0] / norm;
  for ( i = 1; i < num_source; ++i )
    fracVec[i] = fracVec[i-1] + ratioVec[i] / norm;
}


/*
 *
 * functions to draw masses from mass distribution
 *
 */
   
void drawMassFromSource( SimInspiralTable* table )
{ 
  REAL4 m1, m2, eta;
  int index=0;

  /* choose masses from the mass-list */  
  index = (int)( num_mass * XLALUniformDeviate( randParams ) );
  m1=mass_data[index].mass1;
  m2=mass_data[index].mass2;
  
  eta=m1 * m2 / ( ( m1 + m2 ) * ( m1 + m2 ) );
  table->mass1 = m1;
  table->mass2 = m2;
  table->eta = eta;
  table->mchirp = pow( eta, 0.6) * (m1 + m2); 
}

/*
 *
 * functions to draw sky location from source distribution
 *
 */
void drawLocationFromSource( SimInspiralTable* table )
{
  REAL4 u;
  int i;
 
  u=XLALUniformDeviate( randParams );

  /* draw from the source table */
  for ( i = 0; i < num_source; ++i )
  {
    if ( u < fracVec[i] )
    {
      /* put the parameters */
      table->longitude = source_data[i].ra;
      table->latitude  = source_data[i].dec;
      return;
    }
  }

  /* now then, draw from MilkyWay
     WARNING: This sets location AND distance */
  table=XLALRandomInspiralMilkywayLocation( table, randParams );
}

/*
 *
 * functions to draw sky location from exttrig source file
 *
 */
void drawLocationFromExttrig( SimInspiralTable* table )
{
  static LALStatus status;
  
  LIGOTimeGPS timeGRB;  /* real time of the GRB */
  LALMSTUnitsAndAcc unitsAndAcc; 
  REAL4 ra_rad, de_rad;
  REAL8 gmst1, gmst2;  

  /* convert the position (stored as degree) to radians first */
  ra_rad = exttrigHead->event_ra  * LAL_PI_180;
  de_rad = exttrigHead->event_dec * LAL_PI_180;

  /* set units and accuracy for GMST calculation*/
  unitsAndAcc.accuracy = LALLEAPSEC_STRICT;
  unitsAndAcc.units = MST_RAD;

  /* populate the time structures */
  timeGRB.gpsSeconds     = exttrigHead->start_time;
  timeGRB.gpsNanoSeconds = exttrigHead->start_time_ns;

  LALGPStoGMST1( &status, &gmst1, &timeGRB, &unitsAndAcc );
  LALGPStoGMST1( &status, &gmst2, &table->geocent_end_time, &unitsAndAcc );

  /* populate the table */
  table->longitude = ra_rad- gmst1 + gmst2;
  table->latitude  = de_rad;
}


/*
 *
 * functions to draw distance from source distribution
 *
 */
void drawDistanceFromSource( SimInspiralTable* table )
{ 
  REAL4 u;
  int i;
 
  u=XLALUniformDeviate( randParams );

  /* draw from the source table */
  for ( i = 0; i < num_source; ++i )
  {
    if ( u < fracVec[i] )
    {
      /* put the parameters, convert to Mpc */
      table->distance = source_data[i].dist/1000.0;
      return;
    }
  }

  /* now then, draw from MilkyWay
     WARNING: This sets location AND distance */
  table=XLALRandomInspiralMilkywayLocation( table, randParams );
}

/* 
 *
 * generate all parameters (sky position and angles) for a random inspiral 
 *
 */
int main( int argc, char *argv[] )
{ 
  LIGOTimeGPS gpsStartTime;
  LIGOTimeGPS gpsEndTime;
  long gpsDuration;

  REAL8 meanTimeStep = 2630 / LAL_PI; /* seconds between injections     */
  REAL8 timeInterval = 0;
  REAL4 fLower = -1;
  REAL4 eps=0.01;  /* needed for some awkward spinning injections */
  
  size_t ninj;
  int rand_seed = 1;

  /* waveform */
  CHAR waveform[LIGOMETA_WAVEFORM_MAX];
  CHAR dummy[256];
  /* xml output data */
  CHAR                  fname[256];
  CHAR                 *userTag = NULL;
  LALLeapSecAccuracy    accuracy = LALLEAPSEC_LOOSE;
  MetadataTable         proctable;
  MetadataTable         procparams;
  MetadataTable         injections;
  ProcessParamsTable   *this_proc_param;
  LIGOLwXMLStream       xmlfp;

  status=blank_status;
  gpsStartTime.gpsSeconds=-1;
  gpsEndTime.gpsSeconds=-1;

  /* getopt arguments */
  /* available letters: c q v x y z C H Q R S T V W X Y  */
  struct option long_options[] =
  {
    {"help",                          no_argument, 0,                'h'},
    {"source-file",             required_argument, 0,                'f'},
    {"mass-file",               required_argument, 0,                'm'},
    {"exttrig-file",            required_argument, 0,                'E'},
    {"f-lower",                 required_argument, 0,                'F'},
    {"gps-start-time",          required_argument, 0,                'a'},
    {"gps-end-time",            required_argument, 0,                'b'},
    {"time-step",               required_argument, 0,                't'},
    {"time-interval",           required_argument, 0,                'i'},
    {"seed",                    required_argument, 0,                's'},
    {"waveform",                required_argument, 0,                'w'},
    {"user-tag",                required_argument, 0,                'Z'},
    {"userTag",                 required_argument, 0,                'Z'},
    {"m-distr",                 required_argument, 0,                'd'},
    {"min-mass1",               required_argument, 0,                'j'},
    {"max-mass1",               required_argument, 0,                'k'},
    {"min-mass2",               required_argument, 0,                'J'},
    {"max-mass2",               required_argument, 0,                'K'},
    {"min-mtotal",              required_argument, 0,                'A'},
    {"max-mtotal",              required_argument, 0,                'L'},
    {"mean-mass1",              required_argument, 0,                'n'},
    {"mean-mass2",              required_argument, 0,                'N'},
    {"stdev-mass1",             required_argument, 0,                'o'},
    {"stdev-mass2",             required_argument, 0,                'O'},
    {"min-distance",            required_argument, 0,                'p'},
    {"max-distance",            required_argument, 0,                'r'},
    {"d-distr",                 required_argument, 0,                'e'},
    {"l-distr",                 required_argument, 0,                'l'},
    {"i-distr",                 required_argument, 0,                'I'},
    {"inclStd",                 required_argument, 0,                'B'},
    {"enable-milkyway",         required_argument, 0,                'M'},
    {"disable-milkyway",        no_argument,       0,                'D'},
    {"min-spin1",               required_argument, 0,                'g'},
    {"max-spin1",               required_argument, 0,                'G'},
    {"min-spin2",               required_argument, 0,                'u'},
    {"max-spin2",               required_argument, 0,                'U'},
    {"output",                  required_argument, 0,                'P'},
    {"enable-spin",             no_argument,       &spinInjections,    1},
    {"disable-spin",            no_argument,       &spinInjections,    0},
    {0, 0, 0, 0}
  };
  int c;

  /* set up inital debugging values */
  lal_errhandler = LAL_ERR_EXIT;
  set_debug_level( "LALMSGLVL2" );

  /* create the process and process params tables */
  proctable.processTable = (ProcessTable *) 
    calloc( 1, sizeof(ProcessTable) );
  LAL_CALL( LALGPSTimeNow ( &status, &(proctable.processTable->start_time),
        &accuracy ), &status );
  LAL_CALL( populate_process_table( &status, proctable.processTable, 
        PROGRAM_NAME, CVS_REVISION, CVS_SOURCE, CVS_DATE ), &status );
  LALSnprintf( proctable.processTable->comment, LIGOMETA_COMMENT_MAX, " " );
  this_proc_param = procparams.processParamsTable = (ProcessParamsTable *) 
    calloc( 1, sizeof(ProcessParamsTable) );

  /* clear the waveform field */
  memset( waveform, 0, LIGOMETA_WAVEFORM_MAX * sizeof(CHAR) );

  /* parse the arguments */
  while ( 1 )
  {
    /* getopt_long stores long option here */
    int option_index = 0;
    long int gpsinput;
    size_t optarg_len;

    c = getopt_long_only( argc, argv, 
        "hf:m:a:b:t:s:w:i:M:", long_options, &option_index );

    /* detect the end of the options */
    if ( c == - 1 )
    {
      break;
    }

    switch ( c )
    {
      case 0:
        /* if this option set a flag, do nothing else now */
        if ( long_options[option_index].flag != 0 )
        {
          break;
        }
        else
        {
          fprintf( stderr, "error parsing option %s with argument %s\n",
              long_options[option_index].name, optarg );
          exit( 1 );
        }
        break;

      case 'f':
        optarg_len = strlen( optarg ) + 1;
        sourceFileName = calloc( 1, optarg_len * sizeof(char) );
        memcpy( sourceFileName, optarg, optarg_len * sizeof(char) );
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, "string", 
              "%s", optarg );
        break;

      case 'm':
        optarg_len = strlen( optarg ) + 1;
        massFileName = calloc( 1, optarg_len * sizeof(char) );
        memcpy( massFileName, optarg, optarg_len * sizeof(char) );
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, "string", 
              "%s", optarg );
        break;

      case 'E':
        optarg_len = strlen( optarg ) + 1;
        exttrigFileName = calloc( 1, optarg_len * sizeof(char) );
        memcpy( exttrigFileName, optarg, optarg_len * sizeof(char) );
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, "string", 
              "%s", optarg );
        break;

      case 'F':
        fLower = atof( optarg );
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, "float", 
              "%f", fLower );
        break;

      case 'a':
        gpsinput = atol( optarg );
        if ( gpsinput < 441417609 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "GPS start time is prior to " 
              "Jan 01, 1994  00:00:00 UTC:\n"
              "(%ld specified)\n",
              long_options[option_index].name, gpsinput );
          exit( 1 );
        }
        if ( gpsinput > 999999999 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "GPS start time is after " 
              "Sep 14, 2011  01:46:26 UTC:\n"
              "(%ld specified)\n", 
              long_options[option_index].name, gpsinput );
          exit( 1 );
        }
        gpsStartTime.gpsSeconds = gpsinput;
        gpsStartTime.gpsNanoSeconds = 0;	
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, "int", 
              "%ld", gpsinput );
        break;

      case 'b':
        gpsinput = atol( optarg );
        if ( gpsinput < 441417609 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "GPS start time is prior to " 
              "Jan 01, 1994  00:00:00 UTC:\n"
              "(%ld specified)\n",
              long_options[option_index].name, gpsinput );
          exit( 1 );
        }
        if ( gpsinput > 999999999 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "GPS start time is after " 
              "Sep 14, 2011  01:46:26 UTC:\n"
              "(%ld specified)\n", 
              long_options[option_index].name, gpsinput );
          exit( 1 );
        }
	      gpsEndTime.gpsSeconds = gpsinput;
        gpsEndTime.gpsNanoSeconds = 0;
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, "int", 
              "%ld", gpsinput );
        break;

      case 's':
        rand_seed = atoi( optarg );
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, "int", 
              "%d", rand_seed );
        break;

      case 't':
        {
          meanTimeStep = atof( optarg );
          this_proc_param = this_proc_param->next = 
            next_process_param( long_options[option_index].name, "float", 
                "%le", meanTimeStep );
        }
        break;

      case 'i':
        timeInterval = atof( optarg );
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, 
              "float", "%le", timeInterval );
        break;

      case 'w':
        LALSnprintf( waveform, LIGOMETA_WAVEFORM_MAX * sizeof(CHAR), "%s",
            optarg );
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, "string", 
              "%s", optarg );
        break;

      case 'M':
        /* set the luminosity of the Milky Way */
        mwLuminosity = atof( optarg );
        if ( mwLuminosity < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "Milky Way luminosity must be positive" 
              "(%f specified)\n", 
              long_options[option_index].name, mwLuminosity );
          exit( 1 );
        }

        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, "float", 
              "%le", mwLuminosity );
        break;  

      case 'D':
        /* set the luminosity of the Milky Way */
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, "string", 
              "" );
        mwLuminosity = 0;
        break;

      case 'Z':
        /* create storage for the usertag */
        optarg_len = strlen( optarg ) + 1;
        userTag = (CHAR *) calloc( optarg_len, sizeof(CHAR) );
        memcpy( userTag, optarg, optarg_len );

        this_proc_param = this_proc_param->next = (ProcessParamsTable *)
          calloc( 1, sizeof(ProcessParamsTable) );
        LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", 
            PROGRAM_NAME );
        LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "--userTag" );
        LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
        LALSnprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, "%s",
            optarg );
        break;

      case 'd':
        optarg_len = strlen( optarg ) + 1;
        memcpy( dummy, optarg, optarg_len );
        this_proc_param = this_proc_param->next = (ProcessParamsTable *)
          calloc( 1, sizeof(ProcessParamsTable) );
        LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", 
            PROGRAM_NAME );
        LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "--m-distr" );
        LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
        LALSnprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, "%s",
            optarg );

        if (!strcmp(dummy, "source")) 
        {         
          mDistr=massFromSourceFile; 
        } 
        else if (!strcmp(dummy, "totalMass")) 
        {
          mDistr=uniformTotalMass;
        }  
        else if (!strcmp(dummy, "componentMass")) 
        {
          mDistr=uniformComponentMass;
        }  
        else if (!strcmp(dummy, "gaussian")) 
        {
          mDistr=gaussianMassDist;
        } 
	      else if (!strcmp(dummy, "log")) 
	      {
	        mDistr=logComponentMass;
        } 
        else
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "unknown mass distribution: "
              "%s must be one of (source,totalMass,componentMass,gaussian)\n", 
              long_options[option_index].name, optarg );
          exit( 1 );
        }
        break;

      case 'j':
        minMass1 = atof( optarg );
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, 
              "float", "%le", minMass1 );
        break; 

      case 'k':
        maxMass1 = atof( optarg );
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, 
              "float", "%le", maxMass1 );
        break;

      case 'J':
        minMass2 = atof( optarg );
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, 
              "float", "%le", minMass2 );
        break; 

      case 'K':
        maxMass2 = atof( optarg );
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, 
              "float", "%le", maxMass2 );
        break;

     case 'A':
        minMtotal = atof( optarg );
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, 
              "float", "%le", minMtotal );
        break;

      case 'L':
        maxMtotal = atof( optarg );
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, 
              "float", "%le", maxMtotal );
        break;

      case 'n':
        meanMass1 = atof( optarg );
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, 
              "float", "%le", meanMass1 );
        break;

      case 'N':
        meanMass2 = atof( optarg );
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, 
              "float", "%le", meanMass2 );
        break;

      case 'o':
        massStdev1 = atof( optarg );
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, 
              "float", "%le", massStdev1 );
        break;

      case 'O':
        massStdev2 = atof( optarg );
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, 
              "float", "%le", massStdev2 );
        break;

      case 'p':
        /* minimum distance from earth */
        dmin = (REAL4) atof( optarg );
        if ( dmin <= 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "minimum distance must be > 0: "
              "(%f kpc specified)\n",
              long_options[option_index].name, dmin );
          exit( 1 );
        }
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, 
              "float", "%e", dmin );
        break;

      case 'r':
        /* max distance from earth */
        dmax = (REAL4) atof( optarg );
        if ( dmax <= 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "maximum distance must be greater than 0: "
              "(%f kpc specified)\n",
              long_options[option_index].name, dmax );
          exit( 1 );
        }
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, 
              "float", "%e", dmax );
        break;

      case 'e':
        optarg_len = strlen( optarg ) + 1;
        memcpy( dummy, optarg, optarg_len );
        this_proc_param = this_proc_param->next = (ProcessParamsTable *)
          calloc( 1, sizeof(ProcessParamsTable) );
        LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", 
            PROGRAM_NAME );
        LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "--d-distr" );
        LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
        LALSnprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, "%s",
            optarg );

        if (!strcmp(dummy, "source")) 
        {         
          dDistr=distFromSourceFile; 
        } 
        else if (!strcmp(dummy, "uniform")) 
        {
          dDistr=uniformDistance;
        }
        else if (!strcmp(dummy, "log10")) 
        {
          dDistr=uniformLogDistance;
        } 
        else if (!strcmp(dummy, "volume")) 
        {
          dDistr=uniformVolume;
        } 
        else
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "unknown source distribution: "
              "%s, must be one of (uniform, log10, volume, source)\n", 
              long_options[option_index].name, optarg );
          exit( 1 );
        }

        break;
	
      case 'l':
        optarg_len = strlen( optarg ) + 1;
        memcpy( dummy, optarg, optarg_len );
        this_proc_param = this_proc_param->next = (ProcessParamsTable *)
          calloc( 1, sizeof(ProcessParamsTable) );
        LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", 
            PROGRAM_NAME );
        LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "--l-distr" );
        LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
        LALSnprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, "%s",
            optarg );

        if (!strcmp(dummy, "source")) 
        {         
          lDistr=locationFromSourceFile; 
        } 
        else if (!strcmp(dummy, "exttrig")) 
        {
          lDistr=locationFromExttrigFile;    
        } 
        else if (!strcmp(dummy, "random")) 
        {
          lDistr=uniformSkyLocation;        
        } 
        else
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "unknown location distribution: "
              "%s must be one of (source, random)\n", 
              long_options[option_index].name, optarg );
          exit( 1 );
        }

        break;
	
      case 'I':
        optarg_len = strlen( optarg ) + 1;
        memcpy( dummy, optarg, optarg_len );
        this_proc_param = this_proc_param->next = (ProcessParamsTable *)
          calloc( 1, sizeof(ProcessParamsTable) );
        LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", 
            PROGRAM_NAME );
        LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "--i-distr" );
        LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
        LALSnprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, "%s",
            optarg );

        if (!strcmp(dummy, "uniform")) 
        {         
          iDistr=uniformInclDist; 
        } 
        else if (!strcmp(dummy, "gaussian")) 
        {
          iDistr=gaussianInclDist;        
        } 
        else
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "unknown inclination distribution: "
              "%s must be one of (uniform, gaussian)\n", 
              long_options[option_index].name, optarg );
          exit( 1 );
        }

        break;


    case 'B':
      /* gaussian width for inclination */
      inclStd = (REAL4) atof( optarg );
      if ( inclStd <= 0 )
      {
	fprintf( stderr, "invalid argument to --%s:\n"
		 "inclination gaussian width must be greater than 0: "
		 "(%f specified)\n",
		 long_options[option_index].name, dmax );
	exit( 1 );
      }
      this_proc_param = this_proc_param->next = 
	next_process_param( long_options[option_index].name, 
			    "float", "%e", inclStd );
      break;

      case 'P':
        optarg_len = strlen( optarg ) + 1;
        outputFileName = calloc( 1, optarg_len * sizeof(char) );
        memcpy( outputFileName, optarg, optarg_len * sizeof(char) );
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, "string", 
              "%s", optarg );
        break;

      case 'g':
	minSpin1 = atof( optarg );
	this_proc_param = this_proc_param->next =
	  next_process_param( long_options[option_index].name,
			      "float", "%le", minSpin1 );
	break;

       case 'G':
         maxSpin1 = atof( optarg );
         this_proc_param = this_proc_param->next =
           next_process_param( long_options[option_index].name,
               "float", "%le", maxSpin1 );
         break;

       case 'u':
         minSpin2 = atof( optarg );
         this_proc_param = this_proc_param->next =
           next_process_param( long_options[option_index].name,
               "float", "%le", minSpin2 );
         break;

       case 'U':
         maxSpin2 = atof( optarg );
         this_proc_param = this_proc_param->next =
           next_process_param( long_options[option_index].name,
               "float", "%le", maxSpin2 );
         break;

      case 'h':
        printUsage();
        exit( 0 );
        break;

      case '?':
        printUsage();
        exit( 1 );
        break;

      default:
        fprintf( stderr, "unknown error while parsing options\n" );
        printUsage();
        exit( 1 );
    }
  }

  /* must set MW flag */
  if ( mwLuminosity < 0 ) 
  {
    fprintf( stderr, 
	     "Must specify either --enable-milkyway LUM or "
	     "--disable-milkyway\n" );
    printUsage();
    exit( 1 );
  }

  if (gpsStartTime.gpsSeconds==-1 || gpsEndTime.gpsSeconds==-1)
  {
    fprintf( stderr, 
	     "Must specify both --gps-start-time and --gps-end-time.\n");
    printUsage();
    exit( 1 );
  }

  gpsDuration=gpsEndTime.gpsSeconds-gpsStartTime.gpsSeconds;

  if ( dDistr == unknownDistanceDist )
  {
    printf("Must specify a distance distribution (--d-distr).\n");
    printUsage();
    exit( 1 );
  }

  if ( lDistr == unknownLocationDist )
  {
    printf("Must specify a location distribution (--l-distr).\n");
    printUsage();
    exit( 1 );
  }

  if ( mDistr == unknownMassDist )
  {
    printf("Must specify a mass distribution (--m-distr).\n");
    printUsage();
    exit( 1 );
  }

  if ( iDistr == unknownInclDist )
  {
    printf("Must specify an inclination distribution (--i-distr).\n");
    printUsage();
    exit( 1 );
  }

  /* if using source file, check that file and MW choice selected */
  if ( dDistr==distFromSourceFile || lDistr==locationFromSourceFile )
  {
    if ( ! sourceFileName )
    {
      fprintf( stderr, 
          "Must specify --source-file when using --d-distr source \n" );
      printUsage();
      exit( 1 );
    }

    /* read the source distribution here */
    read_source_data( sourceFileName );
  }

  /* check if the source file is specified for distance but NOT for 
     location */
  if ( dDistr==distFromSourceFile && lDistr!=locationFromSourceFile )
  {    
    fprintf( stderr, 
             "WARNING: source file specified for distance "
             "but NOT for location. This might give strange distributions\n" );
  }

  /* check if the location file is specified for location but NOT for 
     distances: GRB case */
  if ( dDistr!=distFromSourceFile && lDistr==locationFromSourceFile &&
       mwLuminosity>0.0 )
  {    
    fprintf( stderr, 
             "WARNING: source file specified for locations "
             "but NOT for distances, while Milky Way injections "
             "are allowed. This might give strange distributions\n" );
  }
  

  /* check selection of masses */
  if ( !massFileName && mDistr==massFromSourceFile )
  {
    fprintf( stderr, 
        "Must specify either a file contining the masses (--mass-file) "\
	     "or choose another mass-distribution (--m-distr).\n" );
    printUsage();
    exit( 1 );
  }

  /* read the masses from the mass file here */
  if ( massFileName && mDistr==massFromSourceFile )
  {
    read_mass_data( massFileName );
  } 

  /* read in the data from the external trigger file */
  if ( lDistr == locationFromExttrigFile && !exttrigFileName )
  {
    fprintf( stderr, 
             "If --l-distr exttrig is specified, must specify " \
             "external trigger XML file using --exttrig-file.\n");
    printUsage();
    exit( 1 );
  } 
  if ( lDistr == locationFromExttrigFile && exttrigFileName )
  {
    numExtTriggers=LALExtTriggerTableFromLIGOLw( &exttrigHead, exttrigFileName,
                                                 0, 1);
    printf("Number of triggers read from the external trigger file: %d\n",
           numExtTriggers);
    
    if (numExtTriggers>1)
    {
      printf("WARNING: Only 1 external trigger expected in the file '%s'",
             exttrigFileName );
    }
    if (numExtTriggers==0)
    {
      printf("ERROR: No external trigger found in file '%s'",exttrigFileName );
      exit(1);
    }
  }

  /* check inclination distribution */
  if ( iDistr==gaussianInclDist && inclStd<0.0 )
  {
    fprintf( stderr, 
	     "Must specify width for gaussian inclination distribution, "\
	     "use --inclStd.\n" );
    printUsage();
    exit( 1 );
  }

  /* require --f-lower be explicit */
  if ( fLower <= 0.0 )
  {
    fprintf( stderr, "--f-lower must be specified and non-zero\n" );
    printUsage();
    exit( 1 );
  }


  /* check for gaussian mass distribution parameters */
  if ( mDistr==gaussianMassDist && (meanMass1 < 0.0 || massStdev1 < 0.0 || 
				    meanMass2 < 0.0 || massStdev2 < 0.0))
  {
    fprintf( stderr, 
	     "Must specify --mean-mass1/2 and --stdev-mass1/2 if choosing"
	     " --m-distr=gaussian\n" );
    printUsage();
    exit( 1 );
  }

  /* check if the mass area is properly specified */
  if ( mDistr!=gaussianMassDist && (minMass1 <0.0 || minMass2 <0.0 || 
				    maxMass1 <0.0 || maxMass2 <0.0) )
  {
    fprintf( stderr, 
	     "Must specify --min-mass1/2 and --max-mass1/2 if choosing"
	     " --m-distr not gaussian\n" );
    printUsage();
    exit( 1 );
  }

   /* check if the maximum total mass is properly specified */
  if ( mDistr!=gaussianMassDist && maxMtotal<(minMass1 + minMass2 ))
  {
    fprintf( stderr, 
	     "Maximum total mass must be larger than minMass1+minMass2\n"); 
    printUsage();
    exit( 1 );
  }

  /* check if total mass is specified */
  if ( maxMtotal<0.0)
  {
    fprintf( stderr, 
	     "Must specify --max-mtotal.\n" );
    printUsage();
    exit( 1 );
  }
  if ( minMtotal<0.0)
  {
    fprintf( stderr, 
             "Must specify --min-mtotal.\n" );
    printUsage();
    exit( 1 );
  }

  if ( dDistr!=distFromSourceFile && (dmin<0.0 || dmax<0.0) )
  {
    fprintf( stderr, 
	     "Must specify --min-distance and --max-distance if "
	     "--d-distr is not source.\n" );
    printUsage();
    exit( 1 );

  } 
  /* check if waveform is specified */    
  if ( !*waveform )
  {
    fprintf( stderr, "No waveform specified (--waveform).\n" );
    printUsage();
    exit( 1 );
  }

  if ( spinInjections==-1 )
  {
    fprintf( stderr, 
	     "Must specify --disable-spin or --enable-spin.\n" );
    printUsage();
    exit( 1 );
  }

  if ( spinInjections==1 )
  {
    /* check that spins are in range 0 - 1 */
    if (minSpin1 < 0. || minSpin2 < 0. || maxSpin1 > 1. || maxSpin2 >1.)
    {
      fprintf( stderr,
               "Spins can only take values between 0 and 1.\n" );
      printUsage();
      exit( 1 );
    }

    /* check max and mins are the correct way around */
    if (minSpin1 > maxSpin1 || minSpin2 > maxSpin2 )
    {
      fprintf( stderr,
               "Minimal spins must be less than maximal spins.\n" );    
      printUsage();
      exit( 1 );
    }
  }
   
  if ( userTag )
  {
    LALSnprintf( fname, sizeof(fname), "HL-INJECTIONS_%d_%s-%d-%d.xml", 
		 rand_seed, userTag, gpsStartTime.gpsSeconds, 
		 gpsDuration);
  }
  else
  {
    LALSnprintf( fname, sizeof(fname), "HL-INJECTIONS_%d-%d-%ld.xml", 
        rand_seed, gpsStartTime, gpsDuration );
  }
  if ( outputFileName ) 
  {
    LALSnprintf( fname, sizeof(fname), "%s", 
        outputFileName);
  }

  /* set up the LAL random number generator */
  LALCreateRandomParams( &status, &randParams, rand_seed );

  this_proc_param = procparams.processParamsTable;
  procparams.processParamsTable = procparams.processParamsTable->next;
  free( this_proc_param );

  /* create the first injection */
  simTable = injections.simInspiralTable = (SimInspiralTable *)
    calloc( 1, sizeof(SimInspiralTable) );

  /* loop over parameter generation until end time is reached */
  ninj = 0;
  while ( 1 )
  {
    /* increase counter */
    ninj++;
		     
    /* store time in table */
    simTable=XLALRandomInspiralTime( simTable, randParams,
				     gpsStartTime, timeInterval );

    /* populate waveform and other parameters */
    memcpy( simTable->waveform, waveform, 
	    sizeof(CHAR) * LIGOMETA_WAVEFORM_MAX );
    simTable->f_lower = fLower;

    /* populate masses */
    if ( mDistr==massFromSourceFile )
    {
      drawMassFromSource( simTable );
    }
    else if ( mDistr==gaussianMassDist )
    { 
      simTable=XLALGaussianInspiralMasses( simTable, randParams,
                                           minMass1, maxMass1,
					   meanMass1, massStdev1,
                                           minMass2, maxMass2, 
					   meanMass2, massStdev2);
    } else {
      simTable=XLALRandomInspiralMasses( simTable, randParams, mDistr,
					 minMass1, maxMass1,
					 minMass2, maxMass2, 
                                         minMtotal, maxMtotal);
    }

    /* popualate distances */
    if ( dDistr == distFromSourceFile )
    {
      drawDistanceFromSource( simTable );
    }
    else
    {
      simTable=XLALRandomInspiralDistance(simTable, randParams, 
					  dDistr, dmin/1000.0, dmax/1000.0);
    }

    /* populate location */
    if ( lDistr == locationFromSourceFile )
    {
      drawLocationFromSource( simTable );
    }
    else if ( lDistr == locationFromExttrigFile )
    {
      drawLocationFromExttrig( simTable );
    }
    else
    {
      simTable=XLALRandomInspiralSkyLocation(simTable, randParams); 
    }
   
    /* populate orientations */
    do {
      simTable=XLALRandomInspiralOrientation(simTable, randParams,
					     iDistr, inclStd);

    } while ( ! strcmp(waveform, "SpinTaylorthreePointFivePN") &&
              ( simTable->inclination < eps ||
		simTable->inclination > LAL_PI-eps) );

    /* populate spins, if required */
    if (spinInjections)
    {
      simTable = XLALRandomInspiralSpins( simTable, randParams, 
					  minSpin1, maxSpin1,
					  minSpin2, maxSpin2);
    }

    /* populate the site specific information */
    LALPopulateSimInspiralSiteInfo( &status, simTable );

    /* increase time step, check if end of loop is reached */
    gpsStartTime=*XLALGPSAdd( &gpsStartTime, meanTimeStep );
    if ( gpsStartTime.gpsSeconds > gpsEndTime.gpsSeconds )
      break;

    /* allocate and go to next SimInspiralTable */
    simTable = simTable->next = (SimInspiralTable *)
      calloc( 1, sizeof(SimInspiralTable) );
  }


  /* destroy the structure containing the random params */
  LAL_CALL(  LALDestroyRandomParams( &status, &randParams ), &status);

  memset( &xmlfp, 0, sizeof(LIGOLwXMLStream) );

 
  LAL_CALL( LALOpenLIGOLwXMLFile( &status, &xmlfp, fname), &status );
  LAL_CALL( LALGPSTimeNow ( &status, &(proctable.processTable->end_time),
			    &accuracy ), &status );
  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlfp, process_table ), 
	    &status );
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlfp, proctable, 
				    process_table ), &status );
  LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlfp ), &status );
  
  if ( procparams.processParamsTable )
  {
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlfp, process_params_table ),
	      &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlfp, procparams, 
				      process_params_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlfp ), &status );
  }
  
  if ( injections.simInspiralTable )
  {
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlfp, sim_inspiral_table ), 
	      &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlfp, injections, 
				      sim_inspiral_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlfp ), &status );   
  }
  
  LAL_CALL( LALCloseLIGOLwXMLFile ( &status, &xmlfp ), &status );

  return 0;
}
