/*
 * Copyright (C) 2007 Duncan Brown, Stephen Fairhurst
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with with program; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 * MA  02111-1307  USA
 */


#include <config.h>
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <FrameL.h>


struct options {
	char *outfile;
	char **infiles;
	int n_infiles;
	int verbose;
};


static int print_usage(const char *prog)
{
	return fprintf(stderr,
"Usage:  %s [option ...] --output filename inputfilename ...\n" \
"\n" \
"The following options are recognized:\n" \
"	--help\n" \
"	--output (required)\n" \
"	--verbose\n", prog);
}


static struct options parse_command_line(int argc, char *argv[])
{
	int c;
	int option_index;
	struct options options = {
		NULL,
		NULL,
		0,
		0
	};
	struct option long_options[] = {
		{"help", no_argument, NULL, 'h'},
		{"output", required_argument, NULL, 'o'},
		{"verbose", no_argument, NULL, 'v'},
		{NULL, 0, NULL, 0}
	};

	do switch(c = getopt_long(argc, argv, "", long_options, &option_index)) {
	/* --output */
	case 'o':
		options.outfile = optarg;
		break;

	/* --help */
	case 'h':
		print_usage(argv[0]);
		exit(0);

	/* --verbose */
	case 'v':
		options.verbose = 1;
		break;

	/* option sets a flag */
	case 0:
		break;

	/* end of arguments */
	case -1:
		if(optind >= argc) {
			/* nothing on command line after options */
			print_usage(argv[0]);
			exit(1);
		}
		options.infiles = &argv[optind];
		options.n_infiles = argc - optind;
		break;

	/* unrecognized option */
	case '?':
		print_usage(argv[0]);
		exit(1);

	/* missing argument for an option */
	case ':':
		print_usage(argv[0]);
		exit(1);
	} while(c != -1);

	if(!options.outfile) {
		/* --output not among command line options, use first
		 * filename on command line as output (ugh). emulates
		 * legacy behaviour */
		if(optind + 1 >= argc) {
			/* not enough file names on command line */
			print_usage(argv[0]);
			exit(1);
		}
		options.outfile = options.infiles[0];
		options.infiles++;
		options.n_infiles--;
	}

	return options;
}


int main(int argc, char *argv[])
{
	static const char history[] = "Created by " PACKAGE "-" VERSION ".";
	FILE *devnull;
	struct options options;
	int i;
	struct FrFile *frfileout;

	options = parse_command_line(argc, argv);

	/* note the hack to silence libframe if --verbose is not given */
	devnull = fopen("/dev/null", "w");
	FrLibIni(NULL, options.verbose || !devnull ? stderr : devnull, 0);
	if(devnull)
		fclose(devnull);

	if(options.verbose)
		fprintf(stderr, "%s: writing to \"%s\"\n", argv[0], options.outfile);

	frfileout = FrFileONewH(options.outfile, 0, history);
	if(!frfileout) {
		fprintf(stderr, "%s: error: could not open output file \"%s\"\n", argv[0], options.outfile);
		exit(1);
	}

	for(i = 0; i < options.n_infiles; i++) {
		struct FrFile *frfilein;
		struct FrameH *frame;

		if(options.verbose)
			fprintf(stderr, "%s: reading \"%s\"\n", argv[0], options.infiles[i]);

		frfilein = FrFileINew(options.infiles[i]);
		if(!frfilein) {
			FrFileOEnd(frfileout);
			fprintf(stderr, "%s: error: input file \"%s\" not found\n", argv[0], options.infiles[i]);
			exit(1);
		}

		while((frame = FrameRead(frfilein))) {
			FrameWrite(frame, frfileout);
			FrameFree(frame);
		}
		FrFileIEnd(frfilein);
	}

	FrFileOEnd(frfileout);

	exit(0);
}
