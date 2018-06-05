// Copyright Matthew Pitkin (2018)

// A code to compare the solar system barycenter delay calculated by the LALSuite
// routines to those calculated using TEMPO2

// Assuming TEMPO2 is installed compile with:
// gcc TEMPO2_vs_LALBarycenter.c -o TEMPO2_vs_LALBarycenter -L$LAL_PREFIX/lib -L$LALPULSAR_PREFIX/lib -L/usr/local/lib -I$LAL_PREFIX/include -I$LALPULSAR_PREFIX/include -I/usr/local/include -lm -lsofa -ltempo2 -lm -llalsupport -llal -llalpulsar

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// LAL headers
#include <lal/LALgetopt.h>
#include <lal/LALBarycenter.h>
#include <lal/LALInitBarycenter.h>
#include <lal/SFTutils.h>
#include <lal/LALString.h>

// TEMPO2 header
#include <tempo2.h>

#define USAGE \
"Usage: %s [options]\n\n"\
" --help (-h)              display this message\n"\
" --site (-s)              set the telescope site name\n"\
" --ra (-r)                the right ascension of the source (rads)\n"\
" --dec (-d)               the declination of the source (rads)\n"\
" --gps-start-time (-t)    the start GPS time-of-arrival at the detector\n"\
" --duration (-D)          the time-span over which to calculate delays (secs) (default: 0)\n"\
" --step (-z)              the time step between outputs (secs) (default: 60)\n"\
" --ephem (-e)             the JPL ephemeris to use (default: DE405)\n"\
" --units (-u)             the time system to use: TDB (default) or TCB\n"\
" --output (-o)            file to output time delays to (default: stdout)\n"\
"\n"

int main( int argc, char **argv ){

  REAL8 ra = 0.;      // right ascension (rads)
  REAL8 dec = 0.;     // declination (rads)
  char *site = NULL;  // telescope
  REAL8 gpstime = 0.; // GPS time for comparison 
  char *ephem = XLALStringDuplicate("DE405"); // JPL ephemeris file to use (e.g. DE405, the default)
  char *units = XLALStringDuplicate("TDB");   // the time system to use (TDB (default) or TCB)
  REAL8 step = 60.;      // time steps between outputs
  REAL8 duration = 0.;  // duration over which to produce outputs
  UINT4 nsteps = 0;
  char *output = NULL;

  REAL8 ssbdtLAL = 0., ssbdtTEMPO2 = 0.;

  struct LALoption long_options[] =
  {
    { "help",                     no_argument,        0, 'h' },
    { "site",                     required_argument,  0, 's' },
    { "ra",                       required_argument,  0, 'r' },
    { "dec",                      required_argument,  0, 'd' },
    { "gps-start-time",           required_argument,  0, 't' },
    { "duration",                 required_argument,  0, 'D' },
    { "step",                     required_argument,  0, 'z' },
    { "ephem",                    required_argument,  0, 'e' },
    { "units",                    required_argument,  0, 'u' },
    { "output",                   required_argument,  0, 'o' },
    { 0, 0, 0, 0 }
  };

  char args[] = "hs:r:d:t:D:z:e:u:o:";
  char *program = argv[0];

  /* get input arguments */
  while(1){
    int option_index = 0;
    int c;

    c = LALgetopt_long( argc, argv, args, long_options, &option_index );
    if ( c == -1 ) /* end of options */
      break;

    switch(c){
      case 0: /* if option set a flag, nothing else to do */
        if ( long_options[option_index].flag )
          break;
        else
          fprintf(stderr, "Error parsing option %s with argument %s\n", long_options[option_index].name, LALoptarg );
	  break;
      case 'h': /* help message */
        fprintf(stderr, USAGE, program);
        exit(0);
      case 's': /* verbose output */
        site = XLALStringDuplicate( LALoptarg );
        break;
      case 'r':
        ra = atof( LALoptarg );
        break;
      case 'd':
        dec = atof( LALoptarg );
        break;
      case 't':
        gpstime = atof( LALoptarg );
        break;
      case 'D':
        duration = atof( LALoptarg );
        if ( duration < 0. ){
          fprintf(stderr, "Error... duration must be positive");
          exit(1);
        }
        break;
      case 'z':
        step = atof( LALoptarg );
        if ( step <= 0. ){
          fprintf(stderr, "Error... step must be positive");
          exit(1);
        }
        break;
      case 'e':
        ephem = XLALStringDuplicate( LALoptarg );
        break;
      case 'u':
        units = XLALStringDuplicate( LALoptarg );
        break;
      case 'o':
        output = XLALStringDuplicate( LALoptarg );
        break;
      case '?':
        fprintf(stderr, "unknown error while parsing options\n" );
	break;
      default:
        fprintf(stderr, "unknown error while parsing options\n" );
	break;
    }
  }

  nsteps = (duration/step) > 1. ? (UINT4)floor(duration/step) : 1;

  // check ephemeris is valid
  if ( strcmp(ephem, "DE200") && strcmp(ephem, "DE405") && strcmp(ephem, "DE421") && strcmp(ephem, "DE430") ){
    fprintf(stderr, "Error... Invalid ephemeris '%s' given\n", ephem);
    exit(1);
  }

  // check units are valid
  if ( strcmp(units, "TDB") && strcmp(units, "TCB") ){
    fprintf(stderr, "Error... Invalid units '%s' given\n", units);
    exit(1);
  }

  // set site name to be consistent with those in TEMPO2's observatories.dat file
  char *temposite = NULL;
  if ( !strcmp( site, "LHO" ) || !strcmp( site, "H1" ) || !strcmp( site, "HANFORD" ) || !strcmp( site, "H2" ) ){
    temposite = XLALStringDuplicate( "HANFORD" );
  }
  else if ( !strcmp( site, "LLO" ) || !strcmp( site, "L1" ) || !strcmp( site, "LIVINGSTON" ) ){
    temposite = XLALStringDuplicate( "LIVINGSTON" );
  }
  else if ( !strcmp( site, "V1" ) || !strcmp( site, "VIRGO" ) ){
    temposite = XLALStringDuplicate( "VIRGO" );
  }
  else if ( !strcmp( site, "K1" ) || !strcmp( site, "KAGRA" ) ){
    temposite = XLALStringDuplicate( "KAGRA" );
  }
  else if ( !strcmp( site, "G1" ) || !strcmp( site, "GEO" ) || !strcmp( site, "GEO600" ) ){
    temposite = XLALStringDuplicate( "GEO600" );
  }
  else{
    fprintf(stderr, "Error... unrecognised observatory '%s'\n", site);
    exit(1);
  }

  // get the site
  LALDetector det;
  det = *XLALGetSiteInfo( site );

  // read in the LAL ephemeris files
  EphemerisData *edat=NULL;
  TimeCorrectionData *tdat=NULL;
  BarycenterInput baryinput;
  EarthState earth;
  EmissionTime emit;
  TimeCorrectionType ttype;
  
  char earthfile[1024], sunfile[1024], timefile[1024];
  snprintf(earthfile, 1024, "earth00-40-%s.dat.gz", ephem);
  snprintf(sunfile, 1024, "sun00-40-%s.dat.gz", ephem);
  if ( !strcmp( units, "TDB" ) ){
    snprintf(timefile, 1024, "tdb_2000-2040.dat.gz");
    ttype = TIMECORRECTION_TDB;
  }
  else{
    snprintf(timefile, 1024, "te405_2000-2040.dat.gz");
    ttype = TIMECORRECTION_TCB;
  }

  edat = XLALInitBarycenter( earthfile, sunfile );
  tdat = XLALInitTimeCorrections( timefile );

  // setup BarycenterInput structure
  baryinput.site.location[0] = det.location[0]/LAL_C_SI;
  baryinput.site.location[1] = det.location[1]/LAL_C_SI;
  baryinput.site.location[2] = det.location[2]/LAL_C_SI;
  baryinput.dInv = 0.; // no parallax
  baryinput.delta = dec;
  baryinput.alpha = ra;
  
  // TEMPO2 pulsar structure
  struct pulsar *psr = NULL;

  // initialise pulsar
  psr = (pulsar *)malloc(sizeof(pulsar)*1);
  MAX_OBSN = 10000; // maximum number of observations
  if ( nsteps > MAX_OBSN ){
    fprintf(stderr, "Number of steps is larger that 10000. Increase 'MAX_OBSN' value in file\n");
    exit(1);
  }
  initialiseOne(psr, 1, 1); // initialise pulsar

  char epfile[MAX_FILELEN]; // JPL ephemeris file for TEMPO
  snprintf(epfile, sizeof(char)*MAX_FILELEN, "%s/ephemeris/%s.1950.2050", getenv(TEMPO2_ENVIRON), ephem);
  strncpy(psr[0].ephemeris, ephem, sizeof(char)*MAX_FILELEN);
  strncpy(psr[0].JPL_EPHEMERIS, epfile, sizeof(char)*MAX_FILELEN);

  if ( !strcmp("TCB", units) ){ psr[0].units = SI_UNITS; }
  if ( !strcmp("TDB", units) ){ psr[0].units = TDB_UNITS; }
 
  // set the site
  strncpy(psr[0].obsn[0].telID, temposite, sizeof(psr[0].obsn[0].telID));
    
  REAL8 batdt = 0.;
    
  // set pulsar position
  psr[0].param[param_raj].val[0] = ra;
  psr[0].param[param_decj].val[0] = dec;
  psr[0].param[param_px].val[0] = 0.; // no parallax
  psr[0].correctTroposphere = 0;
  vectorPulsar(psr, 1);
 
  // set up time-of-arrivals for TEMPO2
  psr[0].nobs = nsteps;
  for ( UINT4 i=0; i < nsteps; i++ ){
    psr[0].obsn[i].delayCorr = 1;
    psr[0].obsn[i].clockCorr = 1;

    // set SAT value of observation
    REAL8 thistime = gpstime + i*step;
    REAL8 fracsec = thistime - floor(thistime);     // if not an integer get remaining fraction of seconds
    REAL8 mjd = 0.; // time as Modified Julian Date
      
    // this time is a GPS time, so this needs to be converted to UTC and then to an MJD
    // NOTE: this requires the times to be integer seconds
    struct tm utc;
    XLALGPSToUTC( &utc, (INT4)floor(thistime) ); // convert GPS to UTC
    mjd = XLALConvertCivilTimeToMJD( &utc );     // convert UTC into MJD format
    mjd += fracsec/86400.;                       // add on fractional seconds
  
    psr[0].obsn[i].sat = (long double)mjd;

    if ( i > 0 ){
      strncpy(psr[0].obsn[i].telID, psr[0].obsn[0].telID, sizeof(psr[0].obsn[0].telID));
    }
  }

  // get TEMPO2 time delays
  formBatsAll(psr, 1); // should contain everything that's required and in the right order

  FILE *fp = NULL;
  if ( output ){
    fp = fopen(output, "w");
  }
  else{
    fp = stdout;
  }

  // get LALBarycenter time delays and output differences
  for ( UINT4 i = 0; i < nsteps; i++ ){
    XLALGPSSetREAL8(&baryinput.tgps, gpstime + i*step);

    // get SSB time delay from LALBarycenter
    XLALBarycenterEarthNew( &earth, &baryinput.tgps, edat, tdat, ttype );
    XLALBarycenter( &emit, &baryinput, &earth );
    ssbdtLAL = emit.deltaT;

    // get TEMPO2 SSB time delay
    batdt = psr[0].obsn[i].correctionTT_TB + psr[0].obsn[i].roemer;
    REAL8 shapirodelay = psr[0].obsn[i].shapiroDelaySun; // only include Shapiro delay for the Sun
    ssbdtTEMPO2 = batdt - shapirodelay;

    //fprintf(stdout, "UT1 correction = %.9Lf\n", psr[0].obsn[i].correctionUT1*1e9/86400.);

    // output time delay difference in nanoseconds
    fprintf(fp, "%.9lf\n", fabs(ssbdtTEMPO2-ssbdtLAL)*1e9);
  }

  if ( output ) { fclose(fp); }

  destroyOne(psr); // free memory
  free(psr);

  return 0;
}

