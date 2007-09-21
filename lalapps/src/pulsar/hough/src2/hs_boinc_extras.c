/*
*  Copyright (C) 2007 Bernd Machenschalk, Reinhard Prix
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

/* Extras for building an Einstein@Home BOINC App from HierarchicalSearch
*/

/* TODO:
   - behavior when boinc_is_standlone()?
   - check for critical sections
*/


/** INCLUDES **/

/* BOINC - needs to be before the #defines in hs_boinc_extras.h */
#include "boinc_api.h"
#include "diagnostics.h"
#include "boinc_zip.h"
#if BOINC_GRAPHICS
#include "graphics_api.h"
#include "graphics_lib.h"
#endif

/* our own win_lib includes patches for chdir() and sleep() */
#ifdef _WIN32
#include "win_lib.h"
#endif

/* probably already included by previous headers, but make sure they are included */
#include <stdlib.h>
#include <string.h>
#if (BOINC_GRAPHICS == 2) && !defined(_MSC_VER)
#include <dlfcn.h>
#endif

/* headers of our own code */
#include "HierarchicalSearch.h"
#include <lal/LogPrintf.h>
#include "hs_boinc_extras.h"

NRCSID(HSBOINCEXTRASCRCSID,"$Id$");


/*^* MACROS *^*/

#define MAX_PATH_LEN 512

/** don't want to include LAL headers just for PI */
#define LAL_PI 3.1415926535897932384626433832795029  /**< pi */


/** compare strings s1 and s2 up to the length of s1 (w/o the '\0'!!)
    and set l to the length */
#define MATCH_START(s1,s2,l) (0 == strncmp(s1,s2,(l=strlen(s1))-1))


/*^* global VARIABLES *^*/

/** The program might have multiple output file(s) that are to be zipped into the archive
    to be returned. A program can "register" files to be finally sent back by calling
    register_output_file(). The function stores the information it gets in the following
    variables (global within this module)
*/
static char **outfiles = NULL;        /**< the names  of the output files */
static int  noutfiles  = 0;           /**< the number of the output files */
static char resultfile[MAX_PATH_LEN]; /**< the name of the file / zip archive to return */

/** FLOPS estimation - may be set by command line option --WUfpops=.
    When set, ((skypoint_counter / total_skypoints) * estimated_flops) is periodically
    reported to the BOINC Client as the number of flops, so that together with information
    from the Workunit Genrator, Scheduler and Validator leads to claiming the Credit that
    the system intends to grant for a Workunit
 **/
static double estimated_flops = -1;

/** hooks for communication with the graphics thread */
int (*boinc_init_graphics_hook)(void (*worker)(void)) = NULL; /**< boinc_init_graphics hook -
								 no graphics if this can't be loaded */
void (*set_search_pos_hook)(float,float) = NULL; /**< updates the search position on the starsphere */
double *fraction_done_hook = NULL; /**< hooks the "fraction done" counter of the graphics */

/** if we don't get these symbols from a dynamic library (BOINC_GRAPHICS == 2) we declare them here */
#if (BOINC_GRAPHICS == 1) || ((BOINC_GRAPHICS == 2) && defined(_MSC_VER))
extern double fraction_done;
extern void set_search_pos(float RAdeg, float DEdeg);
extern int boinc_init_graphics(void (*worker)(void));
#endif
/** allow for telling apps with "dynamic graphics" to not use graphics */
static int no_graphics = 0;

/** worker() doesn't take arguments, so we have to pass it argv/c as global vars :-( */
static int global_argc;
static char **global_argv;


/** variables for checkpointing */
static char* cptfilename;                 /**< name of the checkpoint file */
static FStatCheckpointFile* cptf = NULL;  /**< FStatCheckpointFile structure */
static UINT4 bufsize = 8*1024;            /**< size of output file buffer */
static UINT4 maxsize = 1024*1024;         /**< maximal size of the output file */
static double last_rac, last_dec;         /**< last sky position, set by show_progress(),
					       used by set_checkpoint() */
static UINT4 last_count, last_total;      /**< last template count, see last_rac */




/*^* LOCAL FUNCTION PROTOTYPES *^*/
static void worker (void);
static void sighandler(int);
static int write_checkpoint(char*);
static int is_zipped(const char *);
static int resolve_and_unzip(const char*, char*, const size_t);
static int load_graphics_dll(void);



/*^* FUNCTIONS *^*/

#ifdef _MSC_VER
/** Attempt to load the dlls that are required to display graphics.
   If any of them fail do not start the application in graphics mode.
*/
int load_graphics_dll(void) {
  if (FAILED(__HrLoadAllImportsForDll("GDI32.dll"))) {
    LogPrintf(LOG_NORMAL, "WARNING: Failed to load GDI32.DLL - running w/o graphics\n" );
    return(-1);
  } else
    LogPrintf(LOG_NORMAL, "INFO: GDI32.DLL loaded\n" );
  if (FAILED(__HrLoadAllImportsForDll("OPENGL32.dll"))) {
    LogPrintf(LOG_NORMAL, "WARNING: Failed to load OPENGL32.DLL - running w/o graphics\n" );
    return(-1);
  } else
    LogPrintf(LOG_NORMAL, "INFO: OPENGL32.DLL loaded\n" );
  if (FAILED(__HrLoadAllImportsForDll("GLU32.dll"))) {
    LogPrintf(LOG_NORMAL, "WARNING: Failed to load GLU32.DLL - running w/o graphics\n" );
    return(-1);
  } else
    LogPrintf(LOG_NORMAL, "INFO: GLU32.DLL loaded\n" );
  return(0);
}
#endif

/** freaking LAL's REPORTSTATUS just won't work with any of NDEBUG or 
    LAL_NDEBUG set, so we write our own function that dumps the LALStatus
    based on LogPrintf()
 */
void ReportStatus(LALStatus *status)
{ /* </lalVerbatim> */
  LALStatus *ptr;
  for ( ptr = status; ptr ; ptr = ptr->statusPtr ) {                                         
    LogPrintf ( LOG_NORMAL, "\nLevel %i: %s\n", ptr->level, ptr->Id );
    if ( ptr->statusCode ) {
      LogPrintf ( LOG_NORMAL, "\tStatus code %i: %s\n", ptr->statusCode,
		  ptr->statusDescription );
    } else {
      LogPrintf ( LOG_NORMAL, "\tStatus code 0: Nominal\n" );
    }
    LogPrintf ( LOG_NORMAL, "\tfunction %s, file %s, line %i\n",
                   ptr->function, ptr->file, ptr->line );
  }
  return;
}

/** BOINC-compatible LAL(Apps) error handler */
int BOINC_LAL_ErrHand (LALStatus  *stat,
		       const char *func,
		       const char *file,
		       const int line,
		       volatile const char *id) {
  if (stat->statusCode) {
    fprintf(stderr,
            "Level 0: %s\n"
            "\tFunction call `%s' failed.\n"
            "\tfile %s, line %d\n",
            id, func, file, line );
    ReportStatus(stat);
    LogPrintf (LOG_CRITICAL, "BOINC_LAL_ErrHand(): now calling boinc_finish()\n");
    boinc_finish( COMPUTEFSTAT_EXIT_LALCALLERROR+stat->statusCode );
  }
  /* should this call boinc_finish too?? */
  return 0;
}


/**
  our own signal handler
*/
#ifdef __GLIBC__
  /* needed to define backtrace() which is glibc specific*/
#include <execinfo.h>
#endif /* __GLIBC__ */

/* signal handlers */
static void sighandler(int sig){
  LALStatus *mystat = global_status;
  static int killcounter = 0;
#ifdef __GLIBC__
  /* for glibc stacktrace */
  void *stackframes[64];
  size_t nostackframes;
#endif

  /* lets start by ignoring ANY further occurences of this signal
     (hopefully just in THIS thread, if truly implementing POSIX threads */
  LogPrintfVerbatim(LOG_CRITICAL, "\n");
  LogPrintf (LOG_CRITICAL, "APP DEBUG: Application caught signal %d.\n\n", sig );

  /* ignore TERM interrupts once  */
  if ( sig == SIGTERM || sig == SIGINT ) {
    killcounter ++;
    if ( killcounter >= 4 ) {
      LogPrintf (LOG_CRITICAL, "APP DEBUG: got 4th kill-signal, guess you mean it. Exiting now\n\n");
      boinc_finish(COMPUTEFSTAT_EXIT_USER);
    }
    else
      return;
  } /* termination signals */

  if (mystat)
    LogPrintf (LOG_CRITICAL,   "Stack trace of LAL functions in worker thread:\n");
  while (mystat) {
    LogPrintf (LOG_CRITICAL,   "%s at line %d of file %s\n", mystat->function, mystat->line, mystat->file);
    if (!(mystat->statusPtr)) {
      const char *p=mystat->statusDescription;
      LogPrintf (LOG_CRITICAL,   "At lowest level status code = %d, description: %s\n", mystat->statusCode, p?p:"NO LAL ERROR REGISTERED");
    }
    mystat=mystat->statusPtr;
  }
  
#ifdef __GLIBC__
  /* now get TRUE stacktrace */
  nostackframes = backtrace (stackframes, 64);
  LogPrintf (LOG_CRITICAL,   "Obtained %zd stack frames for this thread.\n", nostackframes);
  LogPrintf (LOG_CRITICAL,   "Use gdb command: 'info line *0xADDRESS' to print corresponding line numbers.\n");
  backtrace_symbols_fd(stackframes, nostackframes, fileno(stderr));
#endif /* __GLIBC__ */
  /* sleep a few seconds to let the OTHER thread(s) catch the signal too... */
  sleep(5);
  boinc_finish(COMPUTEFSTAT_EXIT_SIGNAL);
  return;
} /* sighandler */




/**
  show_progress() just sets some variables,
  so should be pretty fast and can be called several times a second
 */
void show_progress(double rac,  /**< right ascension */
		   double dec,  /**< declination */
		   UINT4 count, /**< current skypoint counter */
		   UINT4 total  /**< total number of skypoints */
		   ) {
  double fraction = (double)count / (double)total;

  /* set globals to be written into next checkpoint */
  last_rac = rac;
  last_dec = dec;
  last_count = count;
  last_total = total;

  /* tell graphics thread about fraction done and sky position */
  if (fraction_done_hook)
    *fraction_done_hook = fraction;
  if (set_search_pos_hook)
    set_search_pos_hook(rac * 180.0/LAL_PI, dec * 180.0/LAL_PI);

  /* tell BOINC client about fraction done and flops so far (faked from estimation) */
  boinc_fraction_done(fraction);
  if (estimated_flops >= 0)
    boinc_ops_cumulative( estimated_flops * fraction, 0 /*ignore IOPS*/ );
}





/**
  this registers a new output file to be zipped into the archive that is returned
  to the server as a result file
 */
void register_output_file(char*filename /**< name of the output file to 'register' */
			  ) {
  int len = strlen(filename)+1;
  outfiles = (char**)realloc(outfiles,(noutfiles+1)*sizeof(char*));
  if (outfiles == NULL) {
    LogPrintf (LOG_CRITICAL, "ERROR: Can't allocate output filename '%s'\n", filename);
    noutfiles = 0;
    return;
  }
  outfiles[noutfiles] = calloc(len,sizeof(char));
  if (outfiles[noutfiles] == NULL) {
    LogPrintf (LOG_CRITICAL, "ERROR: Can't allocate output filename '%s'\n", filename);
    return;
  }
  strncpy(outfiles[noutfiles],filename,len);
  noutfiles++;
}



/**
  check if given file is a zip archive by looking for the zip-magic header 'PK\003\044'
  returns 1 if true, 0 if false, -1 if an error occurred
 */
static int is_zipped ( const char *fname /**< name of the file to check for being zipped */
		       ) {
  FILE *fp;
  char zip_magic[] = {'P', 'K', 3, 4 };
  char file_header[4];

  if ( (fp = fopen( fname, "rb")) == NULL ) {
    LogPrintf (LOG_CRITICAL, "Failed to open '%s' for reading: %d: %s\n", fname,errno,strerror(errno));
#ifdef _MSC_VER
      LogPrintf (LOG_CRITICAL, "Windows system call returned: %d\n", _doserrno);
#endif
    return -1;
  }
  if ( 4 != fread ( file_header, sizeof(char), 4, fp ) ) {
    LogPrintf (LOG_CRITICAL, "Failed to read first 4 bytes from '%s'.\n", fname);
    return -1;
  }
  fclose(fp);

  if ( memcmp ( file_header, zip_magic, 4 ) )
    return 0;	/* false: no zip-file */
  else
    return 1;	/* yep, found magic zip-header */
} /* is_zipped() */



/**
  prepare an input file for the program, i.e. boinc_resolve and/or unzip it if necessary
 */
/* better: if the file is a BOINC symlink to a zipped file, (boinc_resolve succeeds),
   first rename the link, then unzip the file, then remove the renamed link.
   Thus, at the beginning, if the file couldn't be found (i.e. resolved), try to resolve
   the renamed link, and upon success, unzip the file and remove the link.
*/
#define ZIPPED_EXT ".zip"
#define LINKED_EXT ".lnk"
static int resolve_and_unzip(const char*filename, /**< filename to resolve */
			     char*resfilename,    /**< resolved filename */
			     const size_t size    /**< size of the buffer for resolved name */
			     ) {
  char buf[MAX_PATH_LEN]; /**< buffer for filename modifications */
  int zipped; /**< flag: is the file zipped? */
  int ret; /** keep return values */

  ret = boinc_resolve_filename(filename,resfilename,size);
  if (ret) {
    LogPrintf(LOG_CRITICAL,"ERROR %d boinc_resolving file '%s'\n", ret, filename);
    return(-1);
  }
  if (strncmp(filename,resfilename,size) == 0) {
    /* filename wasn't a symbolic link */
    strncpy(buf,filename,sizeof(buf));
    strncat(buf,LINKED_EXT,sizeof(buf));
    if (!boinc_resolve_filename(buf,resfilename,size)) {
      /* this could only be the remainder of a previous interrupted unzip */
      LogPrintf (LOG_NORMAL, "WARNING: found old link file '%s'\n", buf);

      /* unzip */
      if (boinc_zip(UNZIP_IT,resfilename,".") ) {
	LogPrintf (LOG_CRITICAL, "ERROR: Couldn't unzip '%s'\n", resfilename);
	return(-1);
      }

      /* delete the link to avoid later confusion */
      if(boinc_delete_file(buf)) {
	LogPrintf (LOG_CRITICAL, "WARNING: Couldn't delete link '%s'\n", buf);
      }

      /* the new resolved filename is the unzipped file */
      strncpy(resfilename,filename,size);
      return(0);
    }

    zipped = is_zipped (filename);

    if (zipped<0) {
      LogPrintf (LOG_DEBUG, "ERROR: Couldn't open '%s'\n", filename);
      return(-1);

    } else if (zipped) { 

      /** unzip in-place: rename file to file.zip, then unzip it */
      LogPrintf (LOG_NORMAL, "WARNING: Unzipping '%s' in-place\n", filename);
      strncpy(resfilename,filename,size);
      strncat(resfilename,ZIPPED_EXT,size);
      if( boinc_rename(filename,resfilename) ) {
	LogPrintf (LOG_CRITICAL, "ERROR: Couldn't rename '%s' to '%s'\n", filename, resfilename );
	return(-1);
      }
      if( boinc_zip(UNZIP_IT,resfilename,".") ) {
	LogPrintf (LOG_CRITICAL, "ERROR: Couldn't unzip '%s'\n", resfilename);
	return(-1);
      }
    }

    /* copy the filename into resfile as if boinc_resove() had succeeded */
    strncpy(resfilename,filename,size);
    return(0);
  }

  /** we end up here if boinc_resolve found the filename to be a symlink */
  zipped = is_zipped (resfilename);

  /** return if not zipped or couldn't find out because of an error */
  if (zipped <= 0)
    return(zipped);

  /** rename the local link so we can unzip to that name */
  strncpy(buf,filename,sizeof(buf));
  strncat(buf,LINKED_EXT,sizeof(buf));
  if( boinc_rename(filename,buf) ) {
    LogPrintf (LOG_CRITICAL, "ERROR: Couldn't rename '%s' to '%s'\n", filename, buf);
    return(-1);
  }

  /* unzip */
  if ( boinc_zip(UNZIP_IT,resfilename,".") ) {
    LogPrintf (LOG_CRITICAL, "ERROR: Couldn't unzip '%s'\n", resfilename);
    return(-1);
  }

  /* delete the link to avoid later confusion */
  if(boinc_delete_file(buf)) {
    LogPrintf (LOG_CRITICAL, "WARNING: Couldn't delete link '%s'\n", buf);
  }

  /* the new resolved filename is the unzipped file */
  strncpy(resfilename,filename,size);
  return(0);
}



/**
  The worker() ist called either from main() directly or from boinc_init_graphics
  (in a separate thread). It does some funny things to the command line (mostly
  boinc-resolving filenames), then calls MAIN() (from HierarchicalSearch.c), and
  finally handles the output / result file(s) before exiting with boinc_finish().
*/
static void worker (void) {
  int argc    = global_argc;   /**< as worker is defined void worker(void), ... */
  char**argv  = global_argv;   /**< ...  take argc and argv from global variables */
  char**rargv = NULL;          /**< argv and ... */
  int rargc   = global_argc;   /**< ... argc values for calling the MAIN() function of
				    HierarchicalSearch.c. Until we know better, we expect to
				    pass the same number of arguments / options than we got */
  int arg, rarg;               /**< current command-line argument */
  int i;                       /**< loop counter */
  int l;                       /**< length of matched string */
  int res = 0;                 /**< return value of a function call */
  char *startc,*endc,*appc;    /**< pointers for parsing a command-line argument */
  int output_help = 0;         /**< flag: should we write out an additional help string?
				    describing additional command-line arguments handled
			            only by this BOINC-wrapper? */
  int breakpoint = 0;          /**< stop at breakpoint? (for testing the Windows Runtime Debugger) */

  /* init BOINC diagnostics */
  boinc_init_diagnostics(BOINC_DIAG_DUMPCALLSTACKENABLED |
                         BOINC_DIAG_HEAPCHECKENABLED |
                         BOINC_DIAG_ARCHIVESTDERR |
                         BOINC_DIAG_REDIRECTSTDERR |
                         BOINC_DIAG_TRACETOSTDERR);

#ifdef _WIN32
  /* point the Windows Runtime Debugger to the Symbol Store on einstein */
  diagnostics_set_symstore("http://einstein.phys.uwm.edu/symstore");
#endif

  /* try to load the graphics shared object and, if succeeded, hook the symbols */
#if (BOINC_GRAPHICS == 2) && !defined(_MSC_VER)
  if (graphics_lib_handle) {
    if (!(set_search_pos_hook = dlsym(graphics_lib_handle,"set_search_pos"))) {
      LogPrintf (LOG_CRITICAL,   "unable to resolve set_search_pos(): %s\n", dlerror());
      boinc_finish(HIERARCHICALSEARCH_EDLOPEN);
    }
    if (!(fraction_done_hook = dlsym(graphics_lib_handle,"fraction_done"))) {
      LogPrintf (LOG_CRITICAL,   "unable to resolve fraction_done(): %s\n", dlerror());
      boinc_finish(HIERARCHICALSEARCH_EDLOPEN);
    }
  }
  else
    LogPrintf (LOG_CRITICAL,  "graphics_lib_handle NULL: running without graphics\n");
#endif

  /* PATCH THE COMMAND LINE

     The actual parsing of the command line will be left to the
     MAIN() of HierarchicalSearch.c. However, command line arguments
     that can be identified as filenames must be boinc_resolved
     before passing them to the main function.
     We will also look if input files are possibly zipped and unzip
     them as needed. Output filename(s) will be recorded (resolved
     and unresolved) and the flops estimation, if present,
     will be stored for later use.
  */

  /* allocate space for the vectorof arguments passed to MAIN() of
     HierarchicalSearch. None of the operations below _adds_ an argument,
     so it's safe to allocate space for as many arguments as we got */
  rargv = (char**)calloc(1,argc*sizeof(char*));
  if(!rargv){
    LogPrintf(LOG_CRITICAL, "Out of memory\n");
    boinc_finish(HIERARCHICALSEARCH_EMEM);
  }

  /* the program name (argv[0]) remains the same in any case */
  rargv[0] = argv[0];
  rarg = 1;

  /* for all args in the command line (except argv[0]) */
  for (arg=1; arg<argc; arg++) {
    
    /* a possible config file is boinc_resolved, but filenames contained in it are not! */
    if (argv[arg][0] == '@') {
      rargv[rarg] = (char*)calloc(MAX_PATH_LEN,sizeof(char));
      if(!rargv[rarg]){
	LogPrintf(LOG_CRITICAL, "Out of memory\n");
	boinc_finish(HIERARCHICALSEARCH_EMEM);
      }
      rargv[rarg][0] = '@';
      if (boinc_resolve_filename(argv[arg]+1,rargv[rarg]+1,MAX_PATH_LEN-1)) {
        LogPrintf (LOG_NORMAL, "WARNING: Can't boinc-resolve config file '%s'\n", argv[arg]+1);
      }
    }

    /* boinc_resolve and unzip skygrid file */
    else if (MATCH_START("--skyGridFile=",argv[arg],l)) {
      rargv[rarg] = (char*)calloc(MAX_PATH_LEN,sizeof(char));
      if(!rargv[rarg]){
	LogPrintf(LOG_CRITICAL, "Out of memory\n");
	boinc_finish(HIERARCHICALSEARCH_EMEM);
      }
      strncpy(rargv[rarg],argv[arg],l);
      if (resolve_and_unzip(argv[arg]+l, rargv[rarg]+l, MAX_PATH_LEN-l) < 0)
	res = HIERARCHICALSEARCH_EFILE;
    }

    /* boinc_resolve and unzip ephermeris files */
    else if (MATCH_START("--ephemE=",argv[arg],l)) {
      rargv[rarg] = (char*)calloc(MAX_PATH_LEN,sizeof(char));
      if(!rargv[rarg]){
	LogPrintf(LOG_CRITICAL, "Out of memory\n");
	boinc_finish(HIERARCHICALSEARCH_EMEM);
      }
      strncpy(rargv[rarg],argv[arg],l);
      if (resolve_and_unzip(argv[arg]+l, rargv[rarg]+l, MAX_PATH_LEN-l) < 0)
	res = HIERARCHICALSEARCH_EFILE;
    }
    else if (MATCH_START("--ephemS=",argv[arg],l)) {
      rargv[rarg] = (char*)calloc(MAX_PATH_LEN,sizeof(char));
      if(!rargv[rarg]){
	LogPrintf(LOG_CRITICAL, "Out of memory\n");
	boinc_finish(HIERARCHICALSEARCH_EMEM);
      }
      strncpy(rargv[rarg],argv[arg],l);
      if (resolve_and_unzip(argv[arg]+l, rargv[rarg]+l, MAX_PATH_LEN-l) < 0)
	res = HIERARCHICALSEARCH_EFILE;
    }


    /* boinc_resolve SFT files (no unzipping, but dealing with multiple files separated by ';' */
    else if (0 == strncmp("--DataFiles",argv[arg],11)) {
      rargv[rarg] = (char*)calloc(1024,sizeof(char));
      if(!rargv[rarg]){
	LogPrintf(LOG_CRITICAL, "Out of memory\n");
	boinc_finish(HIERARCHICALSEARCH_EMEM);
      }

      /* copy & skip the "[1|2]=" characters, too */
      strncpy(rargv[rarg],argv[arg],13);
      appc = rargv[rarg]+13;
      startc = argv[arg]+13;

      /* skip one set of single quotes if and only if they are surrounding the complete path-string */
      if ((*startc == '\'') && (*(startc+(strlen(startc)-1)) == '\'')) {
        LogPrintf (LOG_DEBUG, "DEBUG: removing quotes from path %s\n", argv[arg]);
	*(startc+strlen(startc)-1) = '\0';
	startc++;
      }

      /* look for multiple paths separated by ';' */
      while((endc = strchr(startc,';'))) {
	*endc = '\0';
	if (boinc_resolve_filename(startc,appc,255)) {
	  LogPrintf (LOG_NORMAL, "WARNING: Can't boinc-resolve input file '%s'\n", startc);
	}

#ifdef _WIN32
	/* for Windows, we have to translate the path separator '/' to '\' */
	{
	  char *c = appc;
	  while((*c != '\0') && (c < appc+255)) {
	    if(*c == '/') *c = '\\';
	    c++;
	  }
	}
#endif
	/* append a ';' to resolved string */
	appc = appc + strlen(appc) + 1;
	*(appc-1) = ';';
	*appc = '\0';

	/* skip the ';' in the original string */
	startc = endc+1;
      }

      /* handle last (or only) filename (comments see above) */
      if (boinc_resolve_filename(startc,appc,255)) {
	LogPrintf (LOG_NORMAL, "WARNING: Can't boinc-resolve input file '%s'\n", startc);
      }

#ifdef _WIN32
      {
	char *c = appc;
	while((*c != '\0') && (c < appc+255)) {
	  if(*c == '/') *c = '\\';
	  c++;
	}
      }
#endif
    }

    /* handle output file:
       these are two similar but not equal cases (long and short option name) */
#define OUTPUT_EXT ".zip"
    else if (MATCH_START("--fnameout=",argv[arg],l)) {
      int s;
      if (boinc_resolve_filename(argv[arg]+l,resultfile,sizeof(resultfile))) {
        LogPrintf (LOG_NORMAL, "WARNING: Can't boinc-resolve result file '%s'\n", argv[arg]+l);
      }
      /* derive the name of the local output file from the boinc-resolved output file */
      startc = strrchr(resultfile,'/');
      if(startc == NULL)
	startc = strrchr(resultfile,'\\');
      if(startc == NULL) {
	/* boinc_resolve() doesn't give us a file outside the current directory, so we can't
	   use the same name for the zip archive and the uncompressed file. So we apend the
	   OUTPUT_EXT to the archive filename */
        s = strlen(argv[arg])+1;
        rargv[rarg] = (char*)calloc(s,sizeof(char));
	if(!rargv[rarg]){
	  LogPrintf(LOG_CRITICAL, "Out of memory\n");
	  boinc_finish(HIERARCHICALSEARCH_EMEM);
	}
        strncpy(rargv[rarg],argv[arg],s);
        strncat(resultfile,OUTPUT_EXT,sizeof(resultfile));
        register_output_file(rargv[rarg]+l);
	LogPrintf (LOG_NORMAL, "WARNING: boinc-resolved result file \"%s\" in local directory - will zip into \"%s\"\n",
		   argv[arg]+l,resultfile);
      } else {
	/* boinc_resolve() points us to a file outside the local directory. We will derive that
	   filename from the returned string, write the output to a local file with that name
	   and at the end zip the output file into an archive boinc_resolve() pointed us to */
	startc++;
	s = l+strlen(startc)+1;
        rargv[rarg] = (char*)calloc(s,sizeof(char));
	if(!rargv[rarg]){
	  LogPrintf(LOG_CRITICAL, "Out of memory\n");
	  boinc_finish(HIERARCHICALSEARCH_EMEM);
	}
	strncpy(rargv[rarg],argv[arg],l);
        strncat(rargv[rarg],startc,s);
	register_output_file(startc);
      }
    }
    else if (0 == strncmp("-o",argv[arg],strlen("-o"))) {
      int s;
      rargv[rarg] = argv[arg]; /* copy the "-o" */
      arg++;                   /* grab next argument */
      rarg++;
      if(arg >= argc) {
	LogPrintf(LOG_CRITICAL,"ERROR in command line: no argument following '-o' option\n");
	res = HIERARCHICALSEARCH_EFILE;
      } else {
	if (boinc_resolve_filename(argv[arg],resultfile,sizeof(resultfile))) {
	  LogPrintf (LOG_NORMAL, "WARNING: Can't boinc-resolve result file '%s'\n", argv[arg]);
	}
	/* derive the name of the local output file from the boinc-resolved output file */
	startc = strrchr(resultfile,'/');
	if(startc == NULL)
	  startc = strrchr(resultfile,'\\');
	if(startc == NULL) {
	  /* see previous case - local filename, add OUTPUT_EXT  */
	  s = strlen(argv[arg])+1;
	  rargv[rarg] = (char*)calloc(s,sizeof(char));
	  if(!rargv[rarg]){
	    LogPrintf(LOG_CRITICAL, "Out of memory\n");
	    boinc_finish(HIERARCHICALSEARCH_EMEM);
	  }
	  strncpy(rargv[rarg],argv[arg],s);
	  strncat(resultfile,OUTPUT_EXT,sizeof(resultfile));
	  register_output_file(rargv[rarg]);
	  LogPrintf (LOG_NORMAL, "WARNING: boinc-resolved result file \"%s\" in local directory - will zip into \"%s\"\n",
		     argv[arg],resultfile);
	} else {
	  /* see previous case - different directory - derive local filename */
	  startc++;
	  s = strlen(startc)+1;
	  rargv[rarg] = (char*)calloc(s,sizeof(char));
	  if(!rargv[rarg]){
	    LogPrintf(LOG_CRITICAL, "Out of memory\n");
	    boinc_finish(HIERARCHICALSEARCH_EMEM);
	  }
	  strncpy(rargv[rarg],startc,s);
	  register_output_file(startc);
	}
      }
    }

    /* set the "flops estimation" */
    else if (MATCH_START("--WUfpops=",argv[arg],l)) {
      estimated_flops = atof(argv[arg]+l);
      rarg--; rargc--; /* this argument is not passed to the main worker function */
    }

    /* set maximal output filesize (roughly - can grow beyond this until the next checkpoint) */
    else if (MATCH_START("--MaxFileSize=",argv[arg],l)) {
      maxsize = 1024*atoi(argv[arg]+l);
      rarg--; rargc--; /* this argument is not passed to the main worker function */
    }

    /* set size of output file buffer */
    else if (MATCH_START("--OutputBufSize=",argv[arg],l)) {
      bufsize = 1024*atoi(argv[arg]+l);
      rarg--; rargc--; /* this argument is not passed to the main worker function */
    }

    /* fire up debugger at breakpoint, solely for testing the debugger (and symbols) */
    else if (MATCH_START("--BreakPoint",argv[arg],l)) {
      breakpoint = -1;
      rarg--; rargc--; /* this argument is not passed to the main worker function */
    }

    /* allow for telling apps with "dynamic graphics" to not use graphics */
    else if (MATCH_START("--NoGraphics",argv[arg],l)) {
      no_graphics = -1;
      rarg--; rargc--; /* this argument is not passed to the main worker function */
    }

    /* record a help otion (to later write help for additional command-line options) */
    else if ((0 == strncmp("--help",argv[arg],strlen("--help"))) ||
	     (0 == strncmp("-h",argv[arg],strlen("--help")))) {
      output_help = 1;
      rargv[rarg] = argv[arg];
    }

    /* any other argument - simply pass unchanged */
    else 
      rargv[rarg] = argv[arg];

    /* next argument */
    rarg++;
  } /* for all command line arguments */


  /* sanity check */
  if (!resultfile) {
      LogPrintf (LOG_CRITICAL, "ERROR: no result file has been specified\n");
      res = HIERARCHICALSEARCH_EFILE;
  }

  /* debug: dump the modified command line */
  /*
  fprintf(stderr,"command line:");
  for(i=0;i<rargc;i++)
    fprintf(stderr," %s",rargv[i]);
  fprintf(stderr,"\n");
  */

  /* if there already was an error, there is no use in continuing */
  if (res) {
    LogPrintf (LOG_CRITICAL, "ERROR: error %d in command-line parsing\n", res);
    boinc_finish(res);
  }

  /* test the debugger (and symbol loading) here if we were told to */
#ifdef _MSC_VER
  /* break on file present */
#define DEBUG_BREAKPOINT_FNAME "EAH_MSC_BREAKPOINT"
  {
    FILE*fp_debug;
    if ((fp_debug=fopen("..\\..\\" DEBUG_BREAKPOINT_FNAME, "r")) || (fp_debug=fopen(DEBUG_BREAKPOINT_FNAME, "r")) ) 
      DebugBreak();
  }

  /* break on command-line option present */
  if (breakpoint)
    DebugBreak();
#elif defined(__GNUC__)
  if (breakpoint)
    attach_gdb();
#endif


  /* CALL WORKER's MAIN()
   */
  res = MAIN(rargc,rargv);
  if (res) {
    LogPrintf (LOG_CRITICAL, "ERROR: MAIN() returned with error '%d'\n",res);
  }

  /* if the program was called for help, we write out usage for command-line options this wrapper adds to it and exit */
  if(output_help) {
    printf("Additional options the BOINC version understands:\n");
    printf("      --WUfpops         REAL     \"flops estimation\", passed to the BOINC client as the number of Flops\n");
    printf("      --MaxFileSize     INT      maximum size the outpufile may grow to befor compacted (in 1k)\n");
    printf("      --OutputBufSize   INT      size of the output file buffer (in 1k)\n");
    printf("      --BreakPoint       -       if present fire up the Windows Runtime Debugger at internal breakpoint (WIN32 only)\n");
    printf("      --NoGraphics       -       if present Apps that dynamically look for graphics capabilities will run without graphics in any case\n");
    boinc_finish(0);
  }


  /* HANDLE OUTPUT FILES
   */

  /* we'll still try to zip and send back what's left from an output file for diagnostics */
  /* in case of an error before any output was written the result will contain the link file */
  if(noutfiles == 0)
    LogPrintf (LOG_CRITICAL, "ERROR: no output file has been specified\n");
/* critical> */
  for(i=0;i<noutfiles;i++)
    if ( 0 == strncmp(resultfile, outfiles[i], sizeof(resultfile)) )
      LogPrintf (LOG_NORMAL, "WARNING: output (%d) and result file are identical (%s) - output not zipped\n", i, resultfile);
    else if ( boinc_zip(ZIP_IT, resultfile, outfiles[i]) ) {
      LogPrintf (LOG_NORMAL, "WARNING: Can't zip output file '%s'\n", outfiles[i]);
    }
/* <critical */

  /* finally set (fl)ops count if given */
  if (estimated_flops >= 0)
    boinc_ops_cumulative( estimated_flops, 0 /*ignore IOPS*/ );

  boinc_finish(res);
} /* worker() */



/**
  the main function of the BOINC App
  deals with boinc_init(_graphics) and calls the worker
*/

int main(int argc, char**argv) {
  FILE* fp_debug;
  int skipsighandler = 0;

  LogSetLevel(LOG_DETAIL); /* as long as we are debugging */

  /* dummy for keeping the RCSIDs */
  if(skipsighandler)
    printf(stderr,"%s\n%s\n",HSBOINCEXTRASHRCSID,HSBOINCEXTRASCRCSID);

  LogPrintf(LOG_NORMAL, "Built at: " __DATE__ " " __TIME__ "\n");

  /* pass argc/v to the worker via global vars */
  global_argc = argc;
  global_argv = argv;

  /* setup windows diagnostics (e.g. redirect stderr into a file!) */
#ifdef _WIN32
  boinc_init_diagnostics(BOINC_DIAG_DUMPCALLSTACKENABLED |
                         BOINC_DIAG_HEAPCHECKENABLED |
                         BOINC_DIAG_ARCHIVESTDERR |
                         BOINC_DIAG_REDIRECTSTDERR |
                         BOINC_DIAG_TRACETOSTDERR);
#endif



  /* debugging support by files */

#define DEBUG_LEVEL_FNAME "EAH_DEBUG_LEVEL"
#define DEBUG_DDD_FNAME   "EAH_DEBUG_DDD"

  LogPrintfVerbatim (LOG_NORMAL, "\n");
  LogPrintf (LOG_NORMAL, "Start of BOINC application '%s'.\n", argv[0]);
  
  /* see if user has a DEBUG_LEVEL_FNAME file: read integer and set lalDebugLevel */
  if ((fp_debug=fopen("../../" DEBUG_LEVEL_FNAME, "r")) || (fp_debug=fopen("./" DEBUG_LEVEL_FNAME, "r")))
    {
      int read_int;

      LogPrintf (LOG_NORMAL, "Found '%s' file\n", DEBUG_LEVEL_FNAME);
      if ( 1 == fscanf(fp_debug, "%d", &read_int ) ) 
	{
	  LogPrintf (LOG_NORMAL, "...containing int: Setting lalDebugLevel -> %d\n", read_int );
	  lalDebugLevel = read_int;
	}
      else
	{
	  LogPrintf (LOG_NORMAL, "...with no parsable int: Setting lalDebugLevel -> 1\n");
	  lalDebugLevel = 1;
	}
      fclose (fp_debug);

    } /* if DEBUG_LEVEL_FNAME file found */

  
#if defined(__GNUC__)
  /* see if user has created a DEBUG_DDD_FNAME file: turn on debuggin using 'ddd' */
  if ((fp_debug=fopen("../../" DEBUG_DDD_FNAME, "r")) || (fp_debug=fopen("./" DEBUG_DDD_FNAME, "r")) ) 
    {
      char commandstring[256];
      char resolved_name[MAXFILENAMELENGTH];
      char *ptr;
      pid_t process_id=getpid();
      
      fclose(fp_debug);
      LogPrintf ( LOG_NORMAL, "Found '%s' file, trying debugging with 'ddd'\n", DEBUG_DDD_FNAME);
      
      /* see if the path is absolute or has slashes.  If it has
	 slashes, take tail name */
      if ((ptr = strrchr(argv[0], '/'))) {
	ptr++;
      } else {
	ptr = argv[0];
      }
      
      /* if file name is an XML soft link, resolve it */
      if (boinc_resolve_filename(ptr, resolved_name, sizeof(resolved_name)))
	LogPrintf (LOG_NORMAL,  "Unable to boinc_resolve_filename(%s), so no debugging\n", ptr);
      else {
	skipsighandler = 1;
	LALSnprintf(commandstring,sizeof(commandstring),"ddd %s %d &", resolved_name ,process_id);
	system(commandstring);
	sleep(20);
      }
    } /* DEBUGGING */
#endif // GNUC



  /* install signal handler */
#ifndef _WIN32
  /* install signal-handler for SIGTERM, SIGINT and SIGABRT(?) 
   * NOTE: it is critical to catch SIGINT, because a user
   * pressing Ctrl-C under boinc should not directly kill the
   * app (which is attached to the same terminal), but the app 
   * should wait for the client to send <quit/> and cleanly exit. 
   */
  boinc_set_signal_handler(SIGTERM, sighandler);
  boinc_set_signal_handler(SIGABRT, sighandler);

  /* install signal handler (for ALL threads) for catching
   * Segmentation violations, floating point exceptions, Bus
   * violations and Illegal instructions */
  if ( !skipsighandler )
    {
      boinc_set_signal_handler(SIGINT,  sighandler);
      boinc_set_signal_handler(SIGSEGV, sighandler);
      boinc_set_signal_handler(SIGFPE,  sighandler);
      boinc_set_signal_handler(SIGILL,  sighandler);
      boinc_set_signal_handler(SIGBUS,  sighandler);
    } /* if !skipsighandler */
#else /* WIN32 */
  signal(SIGTERM, sighandler);
  signal(SIGABRT, sighandler);
  if ( !skipsighandler ) {
    signal(SIGINT,  sighandler);
    signal(SIGSEGV, sighandler);
    signal(SIGFPE,  sighandler);
    signal(SIGILL,  sighandler);
  }
#endif /* WIN32 */



  /* boinc_init variations */
#if (BOINC_GRAPHICS == 2) && defined(_MSC_VER)
  /* We don't load an own DLL on Windows, but we check if we can (manually)
     load the system DLLs necessary to do graphics on Windows, and will run
     without graphics if this fails */
  if((!no_graphics) || (!load_graphics_dll())) {
    int retval;
    set_search_pos_hook = set_search_pos;
    fraction_done_hook = &fraction_done;
    retval = boinc_init_graphics(worker);
    LogPrintf (LOG_CRITICAL, "boinc_init_graphics() returned %d.\n", retval);
    boinc_finish(HIERARCHICALSEARCH_EWORKER );
  }
#elif BOINC_GRAPHICS == 2
  /* Try loading screensaver-graphics as a dynamic library.  If this
     succeeds then extern void* graphics_lib_handle is set, and can
     be used with dlsym() to resolve symbols from that library as
     needed. */
  if (!no_graphics) {
    int retval;
    retval = boinc_init_graphics_lib(worker, argv[0]);
    LogPrintf (LOG_CRITICAL, "ERROR: boinc_init_graphics_lib() returned %d.\n", retval);
    boinc_finish(HIERARCHICALSEARCH_EWORKER );
  }
#elif BOINC_GRAPHICS == 1
  {
    int retval;
    /* if we don't get them from the shared library, use variables local to here */
    set_search_pos_hook = set_search_pos;
    fraction_done_hook = &fraction_done;
    /* no dynamic library, just call boinc_init_graphics() */
    retval = boinc_init_graphics(worker);
    LogPrintf (LOG_CRITICAL, "ERROR: boinc_init_graphics() returned %d\n", retval);
    boinc_finish(HIERARCHICALSEARCH_EWORKER );
  }
#endif /* BOINC_GRAPHICS== */

  /* we end up here only if BOINC_GRAPHICS == 0 or a call to boinc_init_graphics failed */
  boinc_init();
  worker();
  boinc_finish(HIERARCHICALSEARCH_ENORM);
  /* boinc_init_graphics() or boinc_finish() ends the program, we never get here */
  return(0);
}



/* CHECKPOINTING FUNCTIONS */

/** log an I/O error, i.e. source code line no., ferror, errno and strerror, and doserrno on Windows, too */
#ifdef _MSC_VER
#define LOGIOERROR(mess,filename) \
    LogPrintf(LOG_CRITICAL, "ERROR: %s %s: line:%d, doserr:%d, ferr:%d, errno:%d: %s\n",\
	      mess,filename,__LINE__,_doserrno,ferror(fp),errno,strerror(errno))
#else
#define LOGIOERROR(mess,filename) \
    LogPrintf(LOG_CRITICAL, "ERROR: %s %s: line:%d, ferr:%d, errno:%d: %s\n",\
	      mess,filename,__LINE__,ferror(fp),errno,strerror(errno))
#endif

/** init checkpointing and read a checkpoint if already there */
int init_and_read_checkpoint(toplist_t*toplist, /**< the toplist to checkpoint */
			     UINT4*count,       /**< returns the skypoint counter if a checkpoint was found */
			     UINT4 total,       /**< total number of skypoints */
			     char*outputname,   /**< name of checkpointed output file */
			     char*cptname       /**< name of checkpoint file */
			     ) {

  FILE*fp;
  unsigned int  checksum, bytes;
  unsigned long count_read, total_read;
  int scanned;
  int ret;

  /* create data needed for checkpointing */
  if (fstat_cpt_file_create (&cptf, outputname, bufsize, maxsize, toplist))
    return(-2);

  /* open the "checkpointed" (output) file for possible writing */
  if (fstat_cpt_file_open(cptf))
    return(-3);

  /* store the name of the checkpoint file in global cptfilename */
  if(cptname) {
    int s = strlen(cptname)+1;
    cptfilename = (char*)calloc(s,sizeof(char));
    if(!cptfilename){
      LogPrintf(LOG_CRITICAL, "Out of memory\n");
      return(-2);
    }
    strncpy(cptfilename,cptname,s);
  } else {
    /* create an own checkpoint file name if we didn't get passed one */
#define CHECKPOINT_EXT ".cpt"
    int s = strlen(outputname)+strlen(CHECKPOINT_EXT)+1;
    cptfilename = (char*)calloc(s,sizeof(char));
    if(!cptfilename){
      LogPrintf(LOG_CRITICAL, "Out of memory\n");
      return(-2);
    }
    strncpy(cptfilename,outputname,s);
    strncat(cptfilename,CHECKPOINT_EXT,s);
  }

  /* now try to open an old checkpoint file with that name */
  fp = fopen(cptfilename,"r");
  
  /* found no previous checkpoint - start from beginning, nothing more to do here */
  if (!fp) {
    LogPrintf (LOG_DEBUG,  "Couldn't open checkpoint (%d) - starting from beginning\n", errno);
    return(0);
  }

  LogPrintf (LOG_DEBUG,  "Found checkpoint - reading...\n");

  errno = 0;

  /* try to read and parse the checkpoint if we found one */
  scanned = fscanf(fp,"%lf,%lf,%lu,%lu,%u,%u\n",
		   &last_rac, &last_dec,
		   &count_read, &total_read,
		   &checksum, &bytes);

  /* that didn't work as expected - deal with various types of errors */
  if (scanned != 6) {
    if (scanned == EOF) {
      /* file was empty */
      LOGIOERROR("ERROR: EOF encountered reading checkpoint", cptfilename);
    } else {
      int c;
      LOGIOERROR("Could't parse checkpoint", cptfilename);
      if(errno){
	/* file-I/O error */
	LogPrintf(LOG_CRITICAL,"scanned %d/6\n", scanned);
      } else {
	char buf[256];
	/* unexpected file content - write it to stderr */
	LogPrintf(LOG_CRITICAL,"scanned %d/6, first 256 bytes were:'", scanned);
	rewind(fp);
	memset(buf,0,sizeof(buf));
	fwrite(buf,sizeof(char),fread(buf,sizeof(char),sizeof(buf)-1,fp),stderr);
	fprintf(stderr,"'\n");
      }
    }
 
    /* try to remove a broken checkpoint */
    if(fclose(fp))
      LOGIOERROR("Couldn't close checkpoint",cptfilename);
    if(remove(cptfilename))
      LOGIOERROR("Couldn't remove broken checkpoint",cptfilename);
    return(-1);
  }

  /* close the checkpoint file after reading */
  if(fclose(fp))
    LOGIOERROR("Couldn't close checkpoint",cptfilename);

  /* compare the number of skypoints read from the checkpoint with that determined from the program input */
  if (total_read != total) {
    LogPrintf (LOG_DEBUG,  "ERROR reading checkpoint: n.o. skypoints inconsistent (%ul != %ul)\n", total_read, total);
    if(remove(cptfilename))
      LOGIOERROR("Couldn't remove broken checkpoint",cptfilename);
    return(-1);
  }

  /* checkpoint seems sensible - try to read the previous output back into the toplist */
  LogPrintf (LOG_DEBUG,  "Read checkpoint - reading previous output...\n");

  ret = fstat_cpt_file_read (cptf, checksum, bytes);

  /* check if something went wrong */
  if (ret < 0) {
    LogPrintf (LOG_DEBUG,  "ERROR reading previous output\n");
    return(-1);
  } else if (ret == 1) {
    return (2); /* "end marker" was found */
  }

  /* make sure the point of next writing is where we stopped reading -
     apparently this isn't necessarily the case with a buffered file on BSD (MacOS) */
  if(fseek(cptf->fp,cptf->bytes,SEEK_SET)) {
    LOGIOERROR("Could't seek to point of last reading", cptf->filename);
    return(-1);
  }
  *count = count_read;

  return(1);
}


/** adds a candidate to the toplist and to the checkpointing file, too, if it was actually inserted
    into the "toplist".
    Compacting, if necessary, is NOT done here, but in set_checkpoint() - doing it here would lead to
    inconsistent state on disk until the next set_checkpoint call.
*/
int add_checkpoint_candidate (toplist_t*toplist,    /**< the toplist */
			      FstatOutputEntry cand /**< the candidate to insert into the toplist */
			      ) {
  if(toplist != cptf->list) {
    LogPrintf (LOG_CRITICAL,  "ERROR: wrong toplist passed to add_checkpoint_candidate()\n", cptfilename);
    return(-2);
  }

  return(fstat_cpt_file_add (cptf, cand));
}


/** actually writes a checkpoint.
    single point to contain the checkpoint format string.
    called only from 2 places in set_checkpoint()
*/
static int write_checkpoint (char*filename) {
  int ret;

  /* not much done here: open ... */
  FILE* fp = fopen(filename,"w");

  if (fp == NULL) {
    LOGIOERROR("Couldn't open checkpoint file",filename);
    return(-1);
  }

  /* ... write ... */
  ret = fprintf(fp,"%lf,%lf,%u,%u,%u,%u\n",
		last_rac, last_dec, last_count, last_total,
		cptf->checksum, cptf->bytes);
  if (ret <= 0) {
    LOGIOERROR("Couldn't write checkpoint",filename);
    if(fclose(fp))
      LOGIOERROR("Couldn't close checkpoint",filename);
    return(ret);
  }

  /* ... close */
  ret = fclose(fp);
  if (ret != 0) {
    LOGIOERROR("Couldn't close checkpoint",filename);
    return(ret);
  }
  
  return(0);
}


/** sets a checkpoint.
    It also "compacts" the output file, i.e. completely rewrites it from
    the toplist in memory, when it has reached the maximum size. When doing
    so it also writes another checkpoint for consistency, regardles of whether
    it's time to checkpoint or not.
*/
void set_checkpoint (void) {
  /* checkpoint every time (i.e. sky position) if FORCE_CHECKPOINTING */
#ifndef FORCE_CHECKPOINTING
  if (boinc_time_to_checkpoint())
#endif
    {
      /* It's time to checkpoint - should we compact, too? */
      if ((cptf->maxsize > 0) && (cptf->bytes >= cptf->maxsize)) {
	/* compact the file */
	fstat_cpt_file_compact(cptf);
	/* write a fresh checkpoint with the new checksum&length */
	if (write_checkpoint(cptfilename))
	  LogPrintf (LOG_CRITICAL, "ERROR: Couldn't write checkpoint file\n", cptfilename);
	fprintf(stderr,"C\n");
      } else {
	/* just checkpoint w/o compacting */
	fstat_cpt_file_flush(cptf);
	/* write a temporary checkpoint, then (atomically) rename */
#define TEMPCHECKPOINT "checkpoint.tmp"
	if (write_checkpoint(TEMPCHECKPOINT)) {
	  LogPrintf (LOG_CRITICAL, "ERROR: Couldn't write checkpoint file\n", cptfilename);
	} else {
	  if(boinc_rename(TEMPCHECKPOINT,cptfilename) ) {
	    LogPrintf (LOG_CRITICAL, "ERROR: Couldn't rename checkpoint file\n", cptfilename);
	  }
	}
	fprintf(stderr,"c\n");
      }
      /* checkpoint done */
      boinc_checkpoint_completed();
      LogPrintf(LOG_DETAIL,"\nset_checkpt(): bytes: %u, file: %d\n", cptf->bytes, ftell(cptf->fp));
    }
#ifndef FORCE_CHECKPOINTING
  else if (cptf->bytes >= cptf->maxsize)
    /* BOINC says it's not yet time for a checkpoint, but we need to compact the file */
    {
      /* don't let the Client interrupt us here */
      boinc_begin_critical_section();
      /* compact and write checkpoint (see above) */
      if (cptf->maxsize > 0)
	fstat_cpt_file_compact(cptf);
      if (write_checkpoint(cptfilename))
	LogPrintf (LOG_CRITICAL, "ERROR: Couldn't write checkpoint file\n", cptfilename);
      boinc_end_critical_section();
      fprintf(stderr,"C\n");
      LogPrintf(LOG_DETAIL,"\nset_checkpt(): bytes: %u, file: %d\n", cptf->bytes, ftell(cptf->fp));
    }
#endif
}

/** finally writes a minimal (compacted) version of the toplist and cleans up
    all structures related to the toplist. After that, the toplist is invalid.
 */
void write_and_close_checkpointed_file (void) {
  fstat_cpt_file_close(cptf);
  fstat_cpt_file_destroy(&cptf);
}

/** attach gdb to the running process. */
void attach_gdb() {
#ifdef __GLIBC__
  char cmd[256];
  pid_t pid=getpid();
  snprintf(cmd, sizeof(cmd), "gdb -pid %d -ex gcore --args %s", pid, global_argv[0]); 
  system(cmd);
  sleep(20);
#endif
}
