#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <lal/FrameCache.h>
#include <lal/LIGOLwXMLRead.h>
#include <lal/LIGOMetadataTables.h>


/*
 ******************************************************************************
 *
 *                            Command-Line Parsing
 *
 ******************************************************************************
 */

struct options {
	char *input_cache_name;
	int verbose;
};


/*
 * Usage message.
 */

static void print_usage(const char *prog)
{
	fprintf(stderr,
"Usage:  %s [options ...] [filename ...]\n" \
"\n" \
"Valid options:\n" \
"	-i,--input-cache filename\n" \
"		Get trigger files from LAL cache file filename.\n",
		prog
	);
}


/*
 * Default values.
 */

static void options_set_defaults(struct options *options)
{
	options->input_cache_name = NULL;
	options->verbose = 0;
}


/*
 * Validate the command line options.
 */

static int options_validate(struct options *options, int argc)
{
	int good = 1;

	good &= (options->input_cache_name != NULL) || argc;

	return good;
}


/*
 * Parse the command line.
 */

static struct options *options_parse(int *argc, char **argv[])
{
	struct options *options = malloc(sizeof(*options));
	struct option long_options[] = {
		{"help", no_argument, NULL, 'h'},
		{"input-cache", required_argument, NULL, 'i'},
		{"verbose", no_argument, &options->verbose, 1},
		{NULL, 0, NULL, 0}
	};
	int c, index;

	options_set_defaults(options);

	do switch(c = getopt_long(*argc, *argv, "hi:", long_options, &index)) {
		/* --help */
		case 'h':
			print_usage(*argv[0]);
			exit(0);

		/* --input */
		case 'i':
			options->input_cache_name = optarg;
			break;

		/* end of options */
		case -1:
			*argc -= optind;
			*argv = &(*argv)[optind];
			break;

		/* option sets a flag */
		case 0:
			break;

		/* missing an argument for an option */
		case ':':
			print_usage(*argv[0]);
			exit(-1);

		/* unrecognized option */
		case '?':
			print_usage(*argv[0]);
			exit(-1);

		/* who knows */
		default:
			print_usage(*argv[0]);
			exit(-1);
	} while(c != -1);

	return options;
}


/*
 * Clean up.
 */

static void options_destroy(struct options *options)
{
	free(options);
}


/*
 ******************************************************************************
 *
 *                        Input File List Construction
 *
 ******************************************************************************
 */

struct inputlist {
	int length;
	char **file;
};


/*
 * Append a file to the input list, and assign it a process id string.
 */

static void inputlist_append(struct inputlist *inputlist, const char *filename)
{
	/* allocate new memory */
	inputlist->file = realloc(inputlist->file, (inputlist->length + 1) * sizeof(*inputlist->file));

	/* copy the file name */
	inputlist->file[inputlist->length] = malloc(strlen(filename) + 1);
	strcpy(inputlist->file[inputlist->length], filename);

	/* increment length */
	inputlist->length++;
}


/*
 * Use LAL's frame cache code to read a cache of trigger files.
 */

static void inputlist_append_from_cache(struct inputlist *inputlist, const char *filename)
{
	FrCache *cache = XLALFrImportCache(filename);
	int offset;
	unsigned i;

	if(!cache)
		exit(-1);

	for(i = 0; i < cache->numFrameFiles; i++) {
		/* find path in URL and append to list*/
		sscanf(cache->frameFiles[i].url, "%*[^:]://%*[^/]%n", &offset);
		inputlist_append(inputlist, &cache->frameFiles[i].url[offset]);
	}
	XLALFrDestroyCache(cache);
}


/*
 * Build input list from cache file and command line.
 */

static struct inputlist *inputlist_build(struct options *options, int argc, char *argv[])
{
	struct inputlist *inputlist = calloc(1, sizeof(*inputlist));

	/* append file names from cache file */
	if(options->input_cache_name)
		inputlist_append_from_cache(inputlist, options->input_cache_name);

	/* append file names from command line */
	while(argc--)
		inputlist_append(inputlist, *argv++);

	return inputlist;
}


/*
 * Clean up.
 */

static void inputlist_destroy(struct inputlist *inputlist)
{
	if(inputlist)
		while(inputlist->length--)
			free(inputlist->file[inputlist->length]);
	free(inputlist);
}


/*
 ******************************************************************************
 *
 *                       Assemble Trigger File Metadata
 *
 ******************************************************************************
 */

struct metadata {
	int length;
	ProcessTable *process;
	ProcessParamsTable *params;
	SearchSummaryTable *summary;
	char **new_id;
};


/*
 * Clean up.
 */

#define LISTFREE(list) do { __typeof(list) tmp; while(list) { tmp = list; list = list->next; free(tmp); } } while(0)

static void metadata_destroy(struct metadata *metadata)
{
	if(metadata) {
		while(metadata->length--)
			free(metadata->new_id[metadata->length]);
		free(metadata->new_id);
		LISTFREE(metadata->process);
		LISTFREE(metadata->params);
		LISTFREE(metadata->summary);
	}
	free(metadata);
}


/*
 * Add a process id string to the new_id array
 */

static void metadata_append_id(struct metadata *metadata, int id)
{
	int length, i;

	/* allocate memory */
	metadata->new_id = realloc(metadata->new_id, (metadata->length + 1) * sizeof(*metadata->new_id));

	/* compute the length of the process id string */
	for(length = 1, i = id; i > 9; i /= 10, length++);
	length += 19 + 1; /* length of "process:process_id:" + \0 */

	/* construct the process id string */
	metadata->new_id[metadata->length] = malloc(length);
	sprintf(metadata->new_id[metadata->length], "process:process_id:%d", id);
}


/*
 * Load the metadata from the tables in one file
 */

static void metadata_from_file(struct metadata *metadata, char *filename, int id)
{
	ProcessTable *process;
	SearchSummaryTable *summary;

	/* read the metadata tables */
	metadata->process = XLALProcessTableFromLIGOLw(filename);
	metadata->params = XLALProcessParamsTableFromLIGOLw(filename);
	metadata->summary = XLALSearchSummaryTableFromLIGOLw(filename);

	/* loop over rows */
	metadata->length = 0;
	process = metadata->process;
	summary = metadata->summary;
	while(process && summary) {
		process = process->next;
		summary = summary->next;
		metadata_append_id(metadata, id++);
		metadata->length++;
	}

	/* check that the process and search summary tables have the same
	 * number of rows */
	if(process || summary) {
		fprintf(stderr, "error: process and search summary tables different lengths in %s\n", filename);
		exit(-1);
	}
}


/*
 * Turn list of input files into array of meta data info.
 */

static struct metadata *metadata_from_inputlist(struct inputlist *inputlist)
{
	struct metadata *metadata = malloc(inputlist->length * sizeof(*metadata));
	int i;
	int id;

	for(i = id = 0; i < inputlist->length; i++) {
		metadata_from_file(&metadata[i], inputlist->file[i], id);
		id += metadata[i].length;
	}

	return metadata;
}


/*
 ******************************************************************************
 *
 *                                Diagnostics
 *
 ******************************************************************************
 */

static void print_process_id_mapping(struct inputlist *inputlist, struct metadata *metadata)
{
	ProcessTable *process;
	int i, j;

	fprintf(stderr, "Process ID mappings:\n");
	for(i = 0; i < inputlist->length; i++) {
		fprintf(stderr, "\t%s:\n", inputlist->file[i]);
		process = metadata->process;
		for(j = 0; j < metadata->length; j++) {
			fprintf(stderr, "\t\t%s --> %s\n", process->node, metadata->new_id[j]);
			process = process->next;
		}
		metadata++;
	}
}


/*
 ******************************************************************************
 *
 *                                Entry Point
 *
 ******************************************************************************
 */

int main(int argc, char *argv[])
{
	struct options *options;
	struct inputlist *inputlist;
	struct metadata *metadata;
	char *prog = argv[0];

	/* parse and validate the command line */
	options = options_parse(&argc, &argv);
	if(!options_validate(options, argc)) {
		fprintf(stderr, "%s: invalid arguments\n", prog);
		print_usage(prog);
		exit(-1);
	}

	/* construct list of input files */
	inputlist = inputlist_build(options, argc, argv);

	/* load all metadata */
	metadata = metadata_from_inputlist(inputlist);

	if(options->verbose)
		print_process_id_mapping(inputlist, metadata);

	/* success */
	metadata_destroy(metadata);
	inputlist_destroy(inputlist);
	options_destroy(options);
	exit(0);
}
