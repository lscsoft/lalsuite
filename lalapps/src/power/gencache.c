#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <lal/LALStdlib.h>
#include <lal/Date.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOLwXMLRead.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/FrameCache.h>
#include <lalapps.h>



RCSID("$Id$");

#define MAXSTR 2048

/* Usage format string. */
#define USAGE "Usage: %s --outcache cachefilename " \
        "[--inputdir dir name]" \
	"[--help]\n"

#define SNGLBURSTREADER_EARG   1
#define SNGLBURSTREADER_EROW   2
#define SNGLBURSTREADER_EFILE  3

#define SNGLBURSTREADER_MSGEARG   "Error parsing arguments"
#define SNGLBURSTREADER_MSGROW    "Error reading row from XML table"
#define SNGLBURSTREADER_MSGEFILE  "Could not open file"

#define TRUE  1
#define FALSE 0

/* Global var */
const char *cacheDir;

/*
 * Command-line options
 */


struct options_t {
	int verbose;
};


/*
 * Set defaults for command-line options.
 * (ANSI C structure initialization sucks)
 */

static void set_option_defaults(struct options_t *options)
{
	options->verbose = FALSE;
	cacheDir = NULL;
}


/*
 * Parse command line arguments.
 */

static void parse_command_line(int argc, char **argv, struct options_t *options, char **outfile)
{
	struct option long_options[] = {
		/* these options set a flag */
		{"verbose",         no_argument,        &options->verbose, TRUE},
		/* parameters which determine the output xml file */
		{"inputdir",        required_argument,  NULL,  'a'},
		{"outcache",        required_argument,  NULL,  'b'},
		{"help",            no_argument,        NULL,  'o'}, 
		{NULL, 0, NULL, 0}
	};
	int c;
	int option_index;

	*outfile = NULL;

	do {
		switch(c = getopt_long(argc, argv, "a:b:", long_options, &option_index)) {
			case -1:
			case 0:
			break;

			case 'a':
			/*
			 * file containing list of xml files to use
			 */
			cacheDir = optarg;
			break;

			case 'b':
			/*
			 * output cache file name
			 */
			*outfile = optarg;
			break;	

			case ':':
			case '?':
			case 'o':
			default:
			/*
			 * print usage
			 */
			LALPrintError(USAGE, *argv);
			exit(SNGLBURSTREADER_EARG);
		}
	} while(c != -1);

	if(optind < argc) {
		fprintf(stderr, "extraneous command line arguments:\n");
		while(optind < argc)
			fprintf(stderr, "%s\n", argv[optind++]);
		exit(SNGLBURSTREADER_EARG);
	}

	if(!*outfile) {
		LALPrintError( "Output cachefile name must be specified\n" );
		exit(SNGLBURSTREADER_EARG);
	}
}


/*
 * Entry Point
 */


int main(int argc, char **argv)
{
	static LALStatus  stat;
	struct options_t  options;

	FrCache *trigCache = NULL;
	CHAR              *outfile;

	/*
	 * Initialize things
	 */

	set_option_defaults(&options);
	parse_command_line(argc, argv, &options, &outfile);


	lal_errhandler = LAL_ERR_EXIT;
	set_debug_level("1");


	/* generate the internal cache */
	LAL_CALL( LALFrCacheGenerate(&stat, &trigCache, cacheDir, "*.xml" ), &stat);
	if ( ! trigCache->numFrameFiles )
	  {
	    fprintf( stderr, "error: no LwXml trigger files found\n");
	    exit( 1 );
	  }

	/* output in a cache file */
	LAL_CALL(LALFrCacheExport(&stat, trigCache, outfile), &stat);
	LAL_CALL( LALDestroyFrCache( &stat, &trigCache ), &stat );

	exit(0);
}
