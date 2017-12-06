/*
 * Copyright (c) 2011 Joshua L. Willis
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
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

/**
 * \file
 * \author Josh Willis
 *
 * \brief Utility to create single-precision FFTW wisdom files
 *
 * This utility is designed to serve as a replacement for the utility fftwf-wisdom that
 * comes with FFTW.  It will generate only the kind of plans that LAL supports, as well
 * as a detailed, human readable description of those plans that is sent to stdout.  This
 * description can then be saved so that the user knows what transforms are included in
 * the wisdom file (which is not itself human-readable).
 *
 * The transforms to be planned must be specified in an input file, with one transform
 * per line.  The format is:
 *
 * \<type\>\<direc\>\<size\>
 *
 * where:
 *
 * - \<type\>     may be 'c' for a complex transform, or 'r' for a real transform
 * - \<direc\>    may be 'f' for a forward transform, or either 'r' or 'b' for a reverse or
 * backward transform. The latter two are equivalent; the dual notation is
 * designed to support both the 'backward' terminology of FFTW and the
 * 'reverse' terminology of LAL.
 * - \<size\>     is the size of the desired transform.
 *
 * Both the \<type\> and \<direc\> specifiers are insensitive to case.
 *
 * The invocation must also specify an output file, where the wisdom itself will be
 * written.  You may redirect stdout to save the human-readable description. This should
 * always be APPENDED to the human readable description of any existing wisdom files
 * (system or otherwise) that were read in as part of the execution; there is no way
 * for lalapps_fftwf_wisdom to print out the human readable description of any existing
 * wisdom read in.
 *
 * So, as an example, suppose that there is an existing system-wide file in
 * /etc/fftw/wisdomf, and that its human readable description is
 * /etc/fftw/wisdomf_description.  Then if there is presently no plan for a complex,
 * reverse transform of size 1048576 (or 2^20), then put the line:
 *
 * cb1048576
 *
 * into the file input_plans (say) and run:
 *
 * cp /etc/fftw/wisdomf_description new_description
 *
 * lalapps_fftwf_wisdom -i input_plans -o new_wisdom -l 3 \>\> new_description
 *
 * When this has finished you can (as root) move new_wisdom to /etc/fftw/wisdomf and
 * new_description to /etc/fftw/wisdomf_description.  The same results could also be
 * achieved by specifying any of cr1048576, CB1048576, etc, in input_plans.
 *
 * Aside from the options specifying input and output files (which are mandatory) there
 * are three other optional arguments.  They are:
 *
 * - -n or --no-system-wisdom  This disables reading in /etc/fftw/wisdomf.  Use this only
 * if you do NOT want the wisdom there to be incorporated into the wisdom file you generate
 * (for instance, if it is outdated because of hardware or FFTW library change). In such a
 * case you should also not append to the existing human readable description, but replace it.
 * - -w \<FILE\> or --wisdom=\<FILE\>  Read in existing wisdom from \<FILE\>.  Program exits if this
 * option is given but the corresponding file is not readable. As with the system wisdom, you
 * should append to a human readable description of the wisdom in this file; what is there
 * already will not be descirbed in the output of lalapps_fftwf_wisdom.
 * - -l \<int\> or --measurelvl=\<int\>  The planning measure level used in creating the plans.
 * Defaults to 3 (exhaustive) if not given.  Levels 0 or 1 (estimate and measure, respectively)
 * are possible, but not appropriate for a system-wide wisdom file, as plans of those levels
 * can be efficiently generated on the fly in application code. Level 3 (exhaustive) is the
 * highest level and would normally be expected to give the best performance, though for large
 * transform sizes it can take hours or even days to create the plans.  The human readable
 * description of the generated wisdom will note the measure level at which it was created.
 *
 * \note All of the human readable descriptions will also indicate that the wisdom was created
 * with FFTW_UNALIGNED.  Users unfamiliar with this may safely ignore it; all LAL plans
 * are presently created with this flag.  The notation is there in case that should change
 * at some point in the future, so that the record of the created wisdom preserves the
 * distinction.
 */

#include <stdlib.h>
#include <stdio.h>
#include <fftw3.h>
#include <getopt.h>
#include <limits.h>  /* For LINE_MAX */

#include <lal/FileIO.h>
#include <lal/StringInput.h>

/* prototypes */
void print_help(void);
int plan_problem(char type, char direc, UINT4 transform_size, int measurelvl);

/** Print basic usage information about the program and exit */
void print_help(void)
{
  fprintf(stderr,"\nlalapps_fftwf_wisdom -i INPUT_FILE -o OUTPUT_FILE [OPTIONS]:\n");
  fprintf(stderr,"This program creates a wisdom file based on an input file specifying transform\n");
  fprintf(stderr,"kinds and sizes.  This is much like a simplified version of the fftwf-wisdom\n");
  fprintf(stderr,"utility provided with FFTW. The transforms to be planned for must be specified\n");
  fprintf(stderr,"through the input file, one transform per line, using the format:\n");
  fprintf(stderr,"   <type><direction><size>\n");
  fprintf(stderr,"where:\n");
  fprintf(stderr,"    <type> is 'r' or 'c'                      for real or complex\n");
  fprintf(stderr,"    <direction> is 'f' or  either 'b' or 'r'  for forward or backward/reverse\n");
  fprintf(stderr,"    <size>                                    is the size of the transform\n");
  fprintf(stderr,"Backward and reverse are synonymous, and <type> and <direction> are not case-\n");
  fprintf(stderr,"sensitive.  For further details and warnings, see the doxygen documentation.\n");
  fprintf(stderr,"\n");
  fprintf(stderr,"Options and their behavior are:\n");
  fprintf(stderr,"\n");
  fprintf(stderr,"   --measurelvl=, -l <int>  The measurelvl argument to plan creation calls. The\n");
  fprintf(stderr,"                            lowest level of zero corresponds to no planning,\n");
  fprintf(stderr,"                            whereas the highest level of 3 could take days, de-\n");
  fprintf(stderr,"                            pending on the transform size. Defaults to 3.\n");
  fprintf(stderr,"     --input=, -i  <file>   Input file containing transforms to plan. Mandatory.\n");
  fprintf(stderr,"     --output=, -o <file>   File in which to save generated wisdom. Mandatory.\n");
  fprintf(stderr,"     --wisdom=, -w <file>   Existing wisdom file to import. Program exits if this\n");
  fprintf(stderr,"                            option is given but the file is unreadable. Optional.\n");
  fprintf(stderr,"     --no-system-wisdom, -n By default, system wisdom file /etc/fftw/wisdomf is \n");
  fprintf(stderr,"                            imported.  If this option is specified, that import is\n");
  fprintf(stderr,"                            disabled.  If the system file is not present or not \n");
  fprintf(stderr,"                            readable, a warning is printed but execution continues.\n");
  fprintf(stderr,"     --help, -h             Print this help message and exit.\n");
  fprintf(stderr,"In addition to writing the accumulated wisdom to the specified output file, the\n");
  fprintf(stderr,"program also prints to stdout a human-readable description of the wisdom success-\n");
  fprintf(stderr,"fully created during this invocation.  That description should be appended to the\n");
  fprintf(stderr,"corresponding description for any wisdom files (including system) read in.\n");
  exit(EXIT_SUCCESS);
}


/**
 * Function used only internally, to create an FFTW plan for a specified problem (thereby adding to wisdom)
 */
int plan_problem(char type,            /**< 'r' for real or 'c' for complex transform */
		 char direc,           /**< 'f' for forward or 'b'/'r' for backward/reverse transform */
		 UINT4 transform_size, /**< Size of transform to plan */
		 int measurelvl)       /**< Level of patience in planning (0 least, 3 most) */
{
  fftwf_plan genericPlan;
  void *indata, *outdata;
  int fwdflag, planning_flags;

  fwdflag = ( (direc=='f') || (direc=='F') );

  /* We call FFTW routines directly, rather than through LAL, so that if XLAL planning routines
     are changed to always read in wisdom, we can still toggle that behavior through the command line.
     In case we ever allow for aligned memory, we allocate everything with fftwf_malloc().
  */

  /* If we ever allow for aligned memory, this will have to toggle depending on input: */
  planning_flags = FFTW_UNALIGNED;

  switch(measurelvl)
    {
    case 0:
      planning_flags |= FFTW_ESTIMATE;
      break;
    default:
    case 3:
      planning_flags |= FFTW_EXHAUSTIVE;
      /* Fall through: */
    case 2:
      planning_flags |= FFTW_PATIENT;
      /* Fall through */
    case 1:
      planning_flags |= FFTW_MEASURE;
      break;
    }

  /* Ugly, but makes us 'locale' independent */

  if ( (type=='r') || (type=='R') )
    {

      indata  = (float *) fftwf_malloc(transform_size*sizeof(float));
      outdata = (float *) fftwf_malloc(transform_size*sizeof(float));

      if ( (!indata) || (!outdata) )
	{
	  if (indata) fftwf_free(indata);
	  if (outdata) fftwf_free(outdata);
	  return 1;
	}

      genericPlan = fftwf_plan_r2r_1d(transform_size,indata,outdata,
				      (fwdflag ? FFTW_R2HC : FFTW_HC2R),
				      planning_flags);
      if (!genericPlan)
	{
	  fftwf_free(indata);
	  fftwf_free(outdata);
	  return 1;
	}
      else
	{
	  fftwf_free(indata);
	  fftwf_free(outdata);
	  fftwf_destroy_plan(genericPlan);
	  return 0;
	}
    }
  else
    {  /* type == 'c' */

      indata  = (fftwf_complex *) fftwf_malloc(transform_size*sizeof(fftwf_complex));
      outdata = (fftwf_complex *) fftwf_malloc(transform_size*sizeof(fftwf_complex));

      if ( (!indata) || (!outdata) )
	{
	  if (indata) fftwf_free(indata);
	  if (outdata) fftwf_free(outdata);
	  return 1;
	}

      genericPlan = fftwf_plan_dft_1d(transform_size,indata,outdata,
				      (fwdflag ? FFTW_FORWARD : FFTW_BACKWARD),
				      planning_flags);

      if (!genericPlan)
	{
	  fftwf_free(indata);
	  fftwf_free(outdata);
	  return 1;
	}
      else
	{
	  fftwf_free(indata);
	  fftwf_free(outdata);
	  fftwf_destroy_plan(genericPlan);
	  return 0;
	}
    }
}

/**
 * Main function
 *
 * Reads command line specifying input and output files, and optionally wisdom file,
 * measure level, and flag preventing import of system wide wisdom, and then creates
 * plans for each problem specified in the input file.  Accumulated wisdom is written
 * to the output file, and a human readable description of the plans successfully
 * created is written to stdout.  Any warnings or errors are written to stderr.
 */
int main(int argc, char **argv)
{
  static int measurelvl=3;
  static int nosys=0;
  UINT4 transform_size;
  char input_line[LINE_MAX];
  char type;
  char direc;
  FILE *infp=NULL, *outfp=NULL, *wisfp=NULL;
  int optindex, optreturn, retval;

  static struct option long_options[] =
    {
      /* Options setting flags */
      {"no-system-wisdom",no_argument,&nosys,1},
      /* Options specifying input/output  */
      {"input",required_argument,NULL,'i'},
      {"output",required_argument,NULL,'o'},
      {"wisdom",required_argument,NULL,'w'},
      {"measurelvl",required_argument,NULL,'l'},
      {"help",no_argument,NULL,'h'},
      {0,0,0,0}
    };

  while ( (optreturn = getopt_long(argc,argv,"ni:o:w:l:h",long_options,&optindex)) != -1)
    {
      switch(optreturn)
	{
	case 0:
	  break;  /* Everything done in setting flag */
	case 'n':
	  nosys=1;
	  break;
	case 'i':
	  infp = LALFopen(optarg,"r+");
	  if (!infp)
	    {
	      fprintf(stderr,"Error: Could not open input file %s\n",optarg);
	      if (outfp) LALFclose(outfp);
	      exit(EXIT_FAILURE);
	    }
	  break;
	case 'o':
	  outfp = LALFopen(optarg,"w+");
	  if (!outfp)
	    {
	      fprintf(stderr,"Error: Could not open output file %s\n",optarg);
	      if (infp) LALFclose(infp);
	      exit(EXIT_FAILURE);
	    }
	  break;
	case 'w':
	  wisfp = LALFopen(optarg,"r+");
	  if (!wisfp)
	    {
	      fprintf(stderr,"Error: Could not open input wisdom file %s for reading\n",optarg);
	      if (infp)  LALFclose(infp);
	      if (outfp) LALFclose(outfp);
	      exit(EXIT_FAILURE);
	    }
	  else
	    {
	      retval = fftwf_import_wisdom_from_file(wisfp);
	      if (!retval)
		{
		  /* Retval is zero if UNsuccessful */
		  fprintf(stderr,"Error: Could not read wisdom from input wisdom file %s\n",optarg);
		  if (infp)  LALFclose(infp);
		  if (outfp) LALFclose(outfp);
		  LALFclose(wisfp);
		  exit(EXIT_FAILURE);
		}
	      LALFclose(wisfp);
	    }
	  fprintf(stderr,"Read in existing wisdom from file %s\n",optarg);
	  break;
	case 'l':
	  if ( sscanf(optarg,"%d",&measurelvl) != 1)
	    {
	      fprintf(stderr,"Error: invalid measure level %s.\n",optarg);
	      if (infp)  LALFclose(infp);
	      if (outfp) LALFclose(outfp);
	      exit(EXIT_FAILURE);
	    }
	  if ( (measurelvl<0) || (measurelvl>3) )
	    {
	      fprintf(stderr,"Error: invalid measure level %d.\n",measurelvl);
	      if (infp)  LALFclose(infp);
	      if (outfp) LALFclose(outfp);
	    }
	  break;
	case 'h': /* Fall through */
	case '?':
	  print_help();
	  break;
	default:
	  exit(EXIT_FAILURE);
	} /* switch(optreturn) */
    } /* while(optreturn != -1) */

  /* Check to make sure mandatory options were given */

  if (!infp)
    {
      fprintf(stderr,"Error: You must specify an input file with -i <FILE> or --input=<FILE>\n");
      if (outfp) LALFclose(outfp);
      exit(EXIT_FAILURE);
    }

  if (!outfp)
    {
      fprintf(stderr,"Error: You must specify an output file with -o <FILE> or --output=<FILE>\n");
      if (infp) LALFclose(infp);
      exit(EXIT_FAILURE);
    }

  /* Only after processing all options do we know if we should read in system wisdom file */

  if (!nosys)
    {
      retval = fftwf_import_system_wisdom();
      if (!retval)
	{
	  /* Retval is zero if UNsuccessful */
	  fprintf(stderr,"Warning: Could not import system wisdom file /etc/fftw/wisdomf\n");
	}
    }
  else
    {
      fprintf(stderr,"Skipped import of system wisdom file /etc/fftw/wisdomf\n");
    }


  /* Process the input file */

  while ( (fgets(input_line,LINE_MAX,infp) != NULL) )
    {
      if (sscanf(input_line,"%c%c%" LAL_UINT4_FORMAT, &type, &direc, &transform_size) == 3)
	{
	  /* Yes, it's ugly, but we don't have to worry about locales: */
	  if ( !( (type=='r') || (type=='R') || (type=='c') || (type=='C') ) )
	    {
	      fprintf(stderr,"Error: Invalid type specifier %c; must be 'r' (real) or 'c' (complex). ",type);
	      fprintf(stderr,"Problem %c%c%" LAL_UINT4_FORMAT " will be skipped!\n", type, direc, transform_size);
	    }
	  else if ( !( (direc=='f') || (direc=='b') || (direc=='r') || (direc=='F') || (direc=='B') || (direc=='R') ) )
	    {
	      fprintf(stderr,"Error: Invalid direction specifier %c; must be 'f' (forward) or 'b'/'r' (backward/reverse). ",
		      direc);
	      fprintf(stderr,"Problem %c%c%" LAL_UINT4_FORMAT " will be skipped!\n",type,direc,transform_size);
	    }
	  else
	    {
	      retval = plan_problem(type,direc,transform_size,measurelvl);
	      if (retval)
		{
		  fprintf(stderr,"Unable to create plan %c%c%" LAL_UINT4_FORMAT "; skipping!\n",
			  type,direc,transform_size);
		}
	      else
		{
		  fprintf(stdout,"Created single-precision %s %s plan, size %" LAL_UINT4_FORMAT
			  " with measure level %d and FFTW_UNALIGNED\n",
			  ( (type=='r') || (type=='R') ) ? "REAL4" : "COMPLEX8",
			  ( (direc=='f') || (direc=='F') ) ? "forward" : "reverse",
			  transform_size, measurelvl);
		}
	    }
	}
      else
	{
	  fprintf(stderr,"Error: Invalid problem specifier. Problem: %s will be skipped\n",input_line);
	}
    }

  fftwf_export_wisdom_to_file(outfp);
  LALFclose(infp);
  LALFclose(outfp);

  exit(EXIT_SUCCESS);
}
