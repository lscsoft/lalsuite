#include <stdio.h>
#include <unistd.h>
#include <lalapps.h>
#include <lal/LALStdlib.h>
#include <lal/LALHello.h>

RCSID( "$Id$" );

#define usgfmt \
  "Usage: %s [options]\n" \
  "Options [default in brackets]:\n" \
  "  -h            print this message\n" \
  "  -v            verbose\n" \
  "  -d dbglvl     set debug level to dbglvl [0]\n" \
  "  -o outfile    use output file outfile [stdout]\n"

#define usage( program ) fprintf( stderr, usgfmt, program )

extern char *optarg;
extern int optind, opterr, optopt;
extern int vrbflg;

int main( int argc, char *argv[] )
{
  const char *program = argv[0];
  const char *outfile = NULL;
  const char *dbglvl  = NULL;
  lal_errhandler_t default_handler;
  LALStatus status = blank_status;
  int code;
  int opt;

  /* parse options */
  while ( 0 < ( opt = getopt( argc, argv, "hVvd:o:" ) ) )
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
      case 'o':
        outfile = optarg;
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

  /* set debug level */
  set_debug_level( dbglvl );

  /* try to call LALHello; catch error LALHELLOH_EOPEN */
  default_handler = lal_errhandler;
  lal_errhandler  = LAL_ERR_RTRN;
  code = LAL_CALL( LALHello( &status, outfile ), &status );
  if ( code == -1 && (status.statusPtr)->statusCode == LALHELLOH_EOPEN )
  {
    fprintf( stderr, "warning: couldn't open file %s for output"
        "(using stdout)\n", outfile );
    clear_status( &status );
    lal_errhandler = LAL_ERR_EXIT;
    LAL_CALL( LALHello( &status, NULL ), &status );
  }
  else if ( code )
  {
    exit( code );
  }
  lal_errhandler = default_handler; /* restore default handler */

  LALCheckMemoryLeaks();
  return 0;
}
