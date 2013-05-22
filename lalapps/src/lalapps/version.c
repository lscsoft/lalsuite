/*
*  Copyright (C) 2007 Jolien Creighton
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

#include <lal/LALMalloc.h>

#include <lalapps.h>
#include <LALAppsVCSInfo.h>

/* program information */
#define PROGRAM_NAME "lalapps_version"

/* verbose flag */
extern int vrbflg;

/* function prototypes */
static void parse_options(int argc, char **argv);
static void print_usage(FILE *ptr, const char *program_name);

/*
 * main program entry point
 */

int main(int argc, char **argv)
{
  /* setup error handler */
  lal_errhandler = LAL_ERR_EXIT;

  /* parse command line options */
  parse_options(argc, argv);

  /* print version information */
  XLALOutputVersionString(stderr, vrbflg);

  /* clean up and exit */
  LALCheckMemoryLeaks();
  return 0;
}

/*
 * local helper functions
 */

/* parse command line options */
static void parse_options(int argc, char **argv)
{
  /* counters */
  int c;

  /* getopt arguments */
  struct option long_options[] =
  {
    /* options that set a flag */
    {"verbose", no_argument, &vrbflg, 1},
    {"help", no_argument, 0, 'h'},
    {0, 0, 0, 0}
  };

  /* initialise vrbflg */
  vrbflg = 0;

  /* parse the arguments */
  while(1)
  {
    /* getopt_long stores option here */
    int option_index = 0;

    /* parse options */
    c = getopt_long_only(argc, argv, "h", long_options, &option_index);

    /* detect the end of the options */
    if (c == -1)
      break;

    switch(c)
    {
      case 0:
        /* if this option sets a flag, do nothing else for now */
        if (long_options[option_index].flag != 0)
          break;
        else
        {
          fprintf(stderr, "Error parsing option %s with argument %s\n",
              long_options[option_index].name, optarg);
          exit(1);
        }
        break;

      case 'h':
        /* print help message, and exit */
        print_usage(stderr, PROGRAM_NAME);
        exit(0);
        break;

      case '?':
        print_usage(stderr, PROGRAM_NAME);
        exit(1);
        break;

      default:
        fprintf(stderr, "Unknown error while parsing arguments\n");
        print_usage(stderr, PROGRAM_NAME);
        exit(1);
        break;
    }
  }

  /* check for extraneous command line arguments */
  if (optind < argc)
  {
    fprintf(stderr, "Extraneous command line arguments:\n");
    while(optind < argc)
      fprintf(stderr, "\t%s\n", argv[optind++]);
    exit(1);
  }
}

/* function to display program usage information */
static void print_usage(FILE *ptr, const char *program_name)
{
  /* print usage information */
  fprintf(ptr,
      "Usage: %s [options]\n"
      " --help display this messgage and exit\n"
      " --verbose display verbose version information\n", program_name);
}
