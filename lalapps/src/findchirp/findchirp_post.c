#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <lalapps.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/FindChirp.h>

RCSID( "$Id$" );

#define LINEBUFSIZ 256

#define usgfmt \
  "Usage: %s [options]\n" \
  "Options [default in brackets]:\n" \
  "  -h            print this message\n" \
  "  -V            print version info\n" \
  "  -v            verbose\n" \
  "  -d dbglvl     set debug level to dbglvl [0]\n" \
  "  -m            maximise over inject lenght (100/pi)\n" \
  "  -s min-max    minimum and maximum of snr [7-14]\n" \
  "  -c min-max    minimum and maximum of chisq [0-50]\n" \
  "  -b bins       number of bins in histogram [10]\n" \
  "  -e eventfile  path to file containing events [events.dat]\n"

#define usage( program ) fprintf( stderr, usgfmt, program )

struct event { INT8 time; REAL4 snr; REAL4 chisq; struct event *next; };

extern char *optarg;
extern int optind, opterr, optopt;
extern int vrbflg;

int main( int argc, char *argv[] )
{
  INT4          i, j, k;
  
  const char   *program = argv[0];
  const char   *dbglvl  = NULL;
  
  REAL4         minsnr    = 7;
  REAL4         maxsnr    = 10;
  REAL4         minchisq  = 0;
  REAL4         maxchisq  = 50;
  INT4          numbins   = 10;
  INT4          max_over_inject = 0;
  CHAR          eventfile[LINEBUFSIZ] = "events.dat";

  REAL4         smin = 1e99, smax = 0, cmin = 1e99, cmax = 0;
  CHAR          linebuf[LINEBUFSIZ];

  struct event *elisthead = NULL;
  struct event *elistelem = NULL;

  UINT4       **histogram;

  REAL4        *snrbin, *chisqbin;
  REAL4         snrdiv, chisqdiv;

  INT8          first_event;
  INT8          stop_time;
  const INT8    window = 1000000000LL * (INT8) ( 1e2 / LAL_PI );

  int           opt;
  FILE         *event_fp  = NULL;
  
  /* parse options */
  while ( 0 < ( opt = getopt( argc, argv, "hVvmd:s:c:b:e:" ) ) )
    {
    switch ( opt )
    {
      case 'h':
        usage( program );
        return 0;
      case 'V':
        PRINT_VERSION( "hello" );
        return 0;
      case 'v':
        vrbflg = 1;
        break;
      case 'd':
        dbglvl = optarg;
        break;
      case 'm':
        max_over_inject = 1;
        break;
      case 's':
        minsnr = atof( optarg );
        maxsnr = atof( ( optarg = strchr( optarg, '-' ) ) ? ++optarg : "" );
        if ( maxsnr < minsnr || minsnr < 0 )
        {
          fprintf( stderr, "invalid range of snr %g-%g\n", 
              minsnr, maxsnr );
          return 1;
        }
        break;
      case 'c':
        minchisq = atof( optarg );
        maxchisq = atof( ( optarg = strchr( optarg, '-' ) ) ? ++optarg : "" );
        if ( maxchisq < minchisq || minchisq < 0 )
        {
          fprintf( stderr, "invalid range of chisq %g-%g\n", 
              minchisq, maxchisq );
          return 1;
        }
        break;
      case 'b':
        numbins = atoi( optarg );
        if ( numbins <= 0 )
        {
          fprintf( stderr, "invalid number of bins %d\n", numbins );
          return 1;
        }
        break;
      case 'e':
        strncpy( eventfile, optarg, LINEBUFSIZ );
        break;
      default:
        usage( program );
        return 1;
    }
  }
  if ( optind < argc )
  {
    usage( program );
    return 1;
  }
  fprintf( stderr, "reading events from %s\n", eventfile );
  fprintf( stderr, "generating 2D histogram with %d bins: ", numbins );
  fprintf( stderr, "snr range: [%g,%g] chisq range: [%g,%g]\n",
      minsnr, maxsnr, minchisq, maxchisq );

  /* open the event file and read in a list of events */
  event_fp = fopen( eventfile, "r" );
  if ( ! event_fp )
  {
    fprintf( stderr, "unable to open open data file\n" );
    return 1;
  }
  while ( fgets( linebuf, LINEBUFSIZ, event_fp ) )
  {
    UINT8 secbuf, nsbuf;

    if ( linebuf[0] == '#' ) continue;
    
    if ( ! elisthead )
    {
      elistelem = elisthead = calloc( 1, sizeof( *elistelem ) );
    }
    else
    {
      elistelem = elistelem->next = calloc( 1, sizeof( *elistelem ) );
    }
    
    if ( sscanf( linebuf, "%lld %lld %*f %*f %*f %*f %*f %f %f %*f\n",
          &secbuf, &nsbuf, &elistelem->snr, &elistelem->chisq ) != 4 )
    {
      fprintf( stderr, "error reading list of events\n" );
      return 1;
    }

    if ( elistelem->snr > smax ) smax = elistelem->snr;
    if ( elistelem->snr < smin ) smin = elistelem->snr;
    if ( elistelem->chisq > cmax ) cmax = elistelem->chisq;
    if ( elistelem->chisq < cmin ) cmin = elistelem->chisq;

    elistelem->time  = 1000000000LL * secbuf + nsbuf;
  }
  fclose( event_fp );
  
  fprintf( stderr, "snr: min = %f max = %f\nchi: min = %f max = %f\n",
      smin, smax, cmin, cmax );
  
  /* set up the bin boundaries for the histogram */
  snrbin     = calloc( numbins + 1, sizeof(REAL4) );
  chisqbin   = calloc( numbins + 1, sizeof(REAL4) );
  snrdiv     = ( maxsnr - minsnr ) / (REAL4) numbins;
  chisqdiv   = ( maxchisq - minchisq ) / (REAL4) numbins;
  for ( i = 0; i <= numbins; ++i )
  {
    snrbin[i] = minsnr + (REAL4) i * snrdiv;
    chisqbin[i] = minchisq + (REAL4) i * chisqdiv;
    fprintf( stderr, "snrbin[%d] = %f\tchisqbin[%d] = %f\n", 
        i, snrbin[i], i, chisqbin[i] );
  }
  histogram = calloc( numbins, sizeof(INT4*) );
  for ( i = 0; i < numbins; ++i )
  {
    histogram[i] = calloc( numbins, sizeof(INT4) );
  }

  /* record the time of the first event */
  first_event  = elisthead->time;

  for ( j = 0; j < numbins; ++j )
  {
    struct event loudest;
    stop_time    = first_event + window;
    memset ( &loudest, 0, sizeof( loudest ) );

    for ( elistelem = elisthead; elistelem; elistelem = elistelem->next )
    {
      if ( elistelem->chisq > chisqbin[j + 1] ) continue;

      if ( max_over_inject )
      {
        if ( elistelem->time < stop_time && elistelem->next )
        {
          if ( elistelem->snr > loudest.snr )
          {
            memcpy( &loudest, elistelem, sizeof( loudest ) );
          }
        }
        else
        {
          for ( k = 0; k < numbins; ++k )
          {
            if ( loudest.snr >= snrbin[k] ) ++histogram[k][j];
          }
          memcpy( &loudest, elistelem, sizeof( loudest ) );

          stop_time = first_event + (INT8) ceil( 
              (double) ( elistelem->time - first_event ) / (double) window )
            * window;
        }
      }
      else
      {
        for ( k = 0; k < numbins; ++k )
        {
          if ( elistelem->snr >= snrbin[k] ) ++histogram[k][j];
        }
      }        
    }
  }

  for ( j = 0; j < numbins; ++j )
  {
    for ( k = 0; k < numbins; ++k )
    {
      fprintf( stdout, "%d\t", histogram[k][j] );
    }
    fprintf( stdout, "\n" );
  }

#if 0
  for ( elistelem = elisthead; elistelem; elistelem = elistelem->next )
  {
    fprintf( stdout, "%lld\t%f\t%f\n", 
        elistelem->time, elistelem->snr, elistelem->chisq );
  }
#endif
  
  return 0;
}
