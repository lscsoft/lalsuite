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

/* BOINC includes - need to be before the #defines in hs_boinc_extras.h */
#include "boinc_api.h"
#include "diagnostics.h"
#include "boinc_zip.h"

/* our own win_lib includes patches for chdir() and sleep() */
#ifdef _WIN32
#include "win_lib.h"
#endif

/* probably already included by previous headers, but make sure they are included */
#include <stdlib.h>
#include <string.h>

/* for finding out and logging the glibc version */
#ifdef __GLIBC__
#define _GNU_SOURCE
#include <gnu/libc-version.h>
#endif

/* our own exception handler / runtime debugger */
#ifdef __MINGW32__
#include "exchndl.h"
#endif

/* headers of our own code */
#include <lal/LogPrintf.h>
#include "HierarchicalSearch.h"
#include "hs_boinc_extras.h"
#include "hs_boinc_options.h"

/* for Linux extended backtrace */
#if defined(__GLIBC__) && defined(__i386__)
#include "erp_execinfo_plus.h"
#endif

/* FIXME: we should probably eliminate the references to this */
#include "ComputeFStatistic.h"

char* HSBOINCEXTRASCRCSID = "$Id$";

#ifdef __APPLE__
#include "EaH_Mac_Icon.h" 
#endif


/*^* MACROS *^*/

#define MAX_PATH_LEN 512

/** don't want to include LAL headers just for PI */
#define LAL_PI 3.1415926535897932384626433832795029  /**< pi */

#ifndef MIN
#define MIN(a,b) ((a)<(b)?(a):(b))
#endif

#define DEBUG_COMMAND_LINE_MANGLING 1

/** compare strings s1 and s2 up to the length of s1 (w/o the '\0'!!)
    and set l to the length */
#define MATCH_START(s1,s2,l) (0 == strncmp(s1,s2,(l=strlen(s1))-1))

/** write the FPU status flags / exception mask bits to stderr */
#define PRINT_FPU_EXCEPTION_MASK(fpstat) \
  if (fpstat & FPU_STATUS_PRECISION)	 \
    fprintf(stderr," PRECISION");	 \
  if (fpstat & FPU_STATUS_UNDERFLOW)	 \
    fprintf(stderr," UNDERFLOW");	 \
  if (fpstat & FPU_STATUS_OVERFLOW)	 \
    fprintf(stderr," OVERFLOW");	 \
  if (fpstat & FPU_STATUS_ZERO_DIVIDE)	 \
    fprintf(stderr," ZERO_DIVIDE");	 \
  if (fpstat & FPU_STATUS_DENORMALIZED)	 \
    fprintf(stderr," DENORMALIZED");	 \
  if (fpstat & FPU_STATUS_INVALID)	 \
    fprintf(stderr," INVALID")

#define PRINT_FPU_STATUS_FLAGS(fpstat) \
  if (fpstat & FPU_STATUS_COND_3)      \
    fprintf(stderr," COND_3");	       \
  if (fpstat & FPU_STATUS_COND_2)      \
    fprintf(stderr," COND_2");	       \
  if (fpstat & FPU_STATUS_COND_1)      \
    fprintf(stderr," COND_1");	       \
  if (fpstat & FPU_STATUS_COND_0)      \
    fprintf(stderr," COND_0");	       \
  if (fpstat & FPU_STATUS_ERROR_SUMM)  \
    fprintf(stderr," ERR_SUMM");       \
  if (fpstat & FPU_STATUS_STACK_FAULT) \
    fprintf(stderr," STACK_FAULT");    \
  PRINT_FPU_EXCEPTION_MASK(fpstat)


/*^* global VARIABLES *^*/

/** the cpu type, see cpu_type_features.h */
int global_cpu_type;

/** output filename - probably not needed to be public anymore */
static char resultfile[MAX_PATH_LEN]; /**< the name of the file / zip archive to return */


/** FLOPS estimation - may be set by command line option --WUfpops=.
    When set, ((skypoint_counter / total_skypoints) * estimated_flops) is periodically
    reported to the BOINC Client as the number of flops, so that together with information
    from the Workunit Genrator, Scheduler and Validator leads to claiming the Credit that
    the system intends to grant for a Workunit
 **/
static double estimated_flops = -1.0;


/** worker() doesn't take arguments, so we have to pass it argv/c as global vars :-( */
static int global_argc;
static char **global_argv;


/** variables for checkpointing */
static char* cptfilename;                 /**< name of the checkpoint file */
static char* outfilename;                 /**< name of the output file */
static toplist_t* toplist;                /**< the toplist we're checkpointing */
static UINT4 last_count, last_total;      /**< last template count, see last_rac */
static BOOLEAN do_sync = -1;              /**< sync checkpoint file to disk, default: yes */



/*^* LOCAL FUNCTION PROTOTYPES *^*/
#ifdef __GLIBC__
static void sighandler(int, siginfo_t*, void*);
#else
static void sighandler(int);
#endif
static void worker (void);
static int is_zipped(const char *);
static int resolve_and_unzip(const char*, char*, const size_t);
static void drain_fpu_stack(void);
static REAL4 get_nan(void);
#ifdef _NO_MSC_VER
#include "graphics_dlls.h"
#include "delayload_dlls.h"
static int try_load_dlls(const char*, const char*);
#endif

typedef UINT2 fpuw_t;
typedef UINT4 ssew_t;
static void   set_fpu_control_word(const fpuw_t word);
static fpuw_t get_fpu_control_word(void);
static fpuw_t get_fpu_status(void);
static ssew_t get_sse_control_status(void);
static void   set_sse_control_status(const ssew_t cword);

/* constants in FPU status word and control word mask */
#define FPU_STATUS_INVALID      (1<<0)
#define FPU_STATUS_DENORMALIZED (1<<1)
#define FPU_STATUS_ZERO_DIVIDE  (1<<2)
#define FPU_STATUS_OVERFLOW     (1<<3)
#define FPU_STATUS_UNDERFLOW    (1<<4)
#define FPU_STATUS_PRECISION    (1<<5)
#define FPU_STATUS_STACK_FAULT  (1<<6)
#define FPU_STATUS_ERROR_SUMM   (1<<7)
#define FPU_STATUS_COND_0       (1<<8)
#define FPU_STATUS_COND_1       (1<<9)
#define FPU_STATUS_COND_2       (1<<10)
#define FPU_STATUS_COND_3       (1<<14)
/* for SSE, status and control information is in the same register
   the status bits 0-5 are identical to the FPU status bits,
   the exception mask bits follow */
#define SSE_MASK_INVALID        (1<<7)
#define SSE_MASK_DENORMALIZED   (1<<8)
#define SSE_MASK_ZERO_DIVIDE    (1<<9)
#define SSE_MASK_OVERFLOW       (1<<10)
#define SSE_MASK_UNDERFLOW      (1<<11)
#define SSE_MASK_PRECISION      (1<<12)


/*^* FUNCTIONS *^*/

#ifdef _NO_MSC_VER
/** Attempt to load the dlls that are required to display graphics.
    returns 0 if successful, -1 in case of a failure.
*/
int try_load_dlls(const char*dlls, const char*mess) {
  char *startc = dlls, *endc = dlls;
  char dll_name[13]; /* DLLs should have 8.3 names */
  
  while((endc = strchr(startc,' '))) {
    memset(dll_name,'\0',sizeof(dll_name));
    strncpy(dll_name, startc, MIN( (endc - startc), sizeof(dll_name) ) );
    if (FAILED(__HrLoadAllImportsForDll(dll_name))) {
      LogPrintf(LOG_NORMAL, mess, dll_name );
      return(-1);
    } else
      LogPrintf(LOG_NORMAL, "INFO: %s loaded\n", dll_name );
    startc = endc + 1;
  }
  return(0);
}
#endif


/** LAL's REPORTSTATUS just won't work with any of NDEBUG or 
    LAL_NDEBUG set, so we write our own function that dumps the LALStatus
    based on LogPrintf()
 */
void ReportStatus(LALStatus *status)
{ /* </lalVerbatim> */
  LALStatus *ptr;
  for ( ptr = status; ptr ; ptr = ptr->statusPtr ) {                                         
    fprintf(stderr, "\nLevel %i: %s\n", ptr->level, ptr->Id );
    if ( ptr->statusCode ) {
      fprintf(stderr, "\tStatus code %i: %s\n", ptr->statusCode,
	      ptr->statusDescription );
    } else {
      fprintf(stderr, "\tStatus code 0: Nominal\n" );
    }
    fprintf(stderr, "\tfunction %s, file %s, line %i\n",
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
#include <signal.h>
#include <execinfo.h>
/* get REG_EIP from ucontext.h, see
   http://www.linuxjournal.com/article/6391 */
#define __USE_GNU
#include <ucontext.h>
#endif /* __GLIBC__ */


#ifdef __GLIBC__
static void sighandler(int sig,
		       siginfo_t *info,
		       void *secret)
#else
static void sighandler(int sig)
#endif /* __GLIBC__ */
{
  static int killcounter = 0;
#ifdef __GLIBC__
  /* for glibc stacktrace */
  static void *stackframes[64];
  static size_t nostackframes;
  ucontext_t *uc = (ucontext_t *)secret;
#endif

  /* lets start by ignoring ANY further occurences of this signal
     (hopefully just in THIS thread, if truly implementing POSIX threads */
  fprintf(stderr, "\n");
  fprintf(stderr, "APP DEBUG: Application caught signal %d.\n\n", sig );

  /* ignore TERM interrupts once  */
  if ( sig == SIGTERM || sig == SIGINT ) {
    killcounter ++;
    if ( killcounter >= 4 ) {
      fprintf(stderr, "APP DEBUG: got 4th kill-signal, guess you mean it. Exiting now\n\n");
      boinc_finish(COMPUTEFSTAT_EXIT_USER);
    }
    else
      return;
  } /* termination signals */

#ifdef __GLIBC__
#ifdef __i386__
  /* in case of a floating-point exception write out the FPU status */
  if ( sig == SIGFPE ) {
    fpuw_t fpstat = uc->uc_mcontext.fpregs->sw;
    fprintf(stderr,"FPU status word %lx, flags: ", uc->uc_mcontext.fpregs->sw);
    PRINT_FPU_STATUS_FLAGS(fpstat);
    fprintf(stderr,"\n");
  }
#endif /* __i386__ */
  /* now get TRUE stacktrace */
  nostackframes = backtrace (stackframes, 64);
  fprintf(stderr,   "Obtained %zd stack frames for this thread.\n", nostackframes);
  fprintf(stderr,   "Use gdb command: 'info line *0xADDRESS' to print corresponding line numbers.\n");
  /* overwrite sigaction with caller's address */
#ifdef __i386__
  stackframes[1] = (void *) uc->uc_mcontext.gregs[REG_EIP];
  backtrace_symbols_fd_plus(stackframes, nostackframes, fileno(stderr));
#else /* __i386__ */
  backtrace_symbols_fd(stackframes, nostackframes, fileno(stderr));
#endif /* __i386__ */
#endif /* __GLIBC__ */

  if (global_status)
    fprintf(stderr,   "Stack trace of LAL functions in worker thread:\n");
  while (global_status) {
    fprintf(stderr,   "%s at line %d of file %s\n", global_status->function, global_status->line, global_status->file);
    if (!(global_status->statusPtr)) {
      const char *p=global_status->statusDescription;
      fprintf(stderr,   "At lowest level status code = %d, description: %s\n", global_status->statusCode, p?p:"NO LAL ERROR REGISTERED");
    }
    global_status=global_status->statusPtr;
  }
  
  /* sleep a few seconds to let the OTHER thread(s) catch the signal too... */
  sleep(5);
  boinc_finish(COMPUTEFSTAT_EXIT_SIGNAL + sig);
  return;
} /* sighandler */




/**
  show_progress() just sets some variables,
  so should be pretty fast and can be called several times a second
 */
void show_progress(REAL8 rac,   /**< right ascension */
		   REAL8 dec,   /**< declination */
		   REAL8 count, /**< current skypoint counter */
		   REAL8 total, /**< total number of skypoints */
		   REAL8 freq,  /**< base frequency */
		   REAL8 fband  /**< frequency bandwidth */
		   ) {
  double fraction = count / total;

  /* set globals to be written into next checkpoint */
  last_count = count;
  last_total = total;

  /* tell BOINC client about fraction done and flops so far (faked from estimation) */
  boinc_fraction_done(fraction);
  if (estimated_flops >= 0)
    boinc_ops_cumulative( estimated_flops * fraction, 0 /*ignore IOPS*/ );

  /* tell APIv6 graphics about status */
  boincv6_progress.skypos_rac = rac;
  boincv6_progress.skypos_dec = dec;
  if(toplist->elems > 0) {
    /* take the last (rightmost) leaf of the heap tree - might not be the
       "best" candidate, but for the graphics it should be good enough */
    HoughFStatOutputEntry *line = (HoughFStatOutputEntry*)(toplist->heap[toplist->elems - 1]);

    boincv6_progress.cand_frequency  = line->Freq;
    boincv6_progress.cand_spindown   = line->f1dot;
    boincv6_progress.cand_rac        = line->Alpha;
    boincv6_progress.cand_dec        = line->Delta;
    boincv6_progress.cand_hough_sign = line->HoughFStat;
    boincv6_progress.frequency       = freq;
    boincv6_progress.bandwidth       = fband;
  } else {
    boincv6_progress.cand_frequency  = 0.0;
    boincv6_progress.cand_spindown   = 0.0;
    boincv6_progress.cand_rac        = 0.0;
    boincv6_progress.cand_dec        = 0.0;
    boincv6_progress.cand_hough_sign = 0.0;
    boincv6_progress.frequency       = 0.0;
    boincv6_progress.bandwidth       = 0.0;
  }
}





/**
  check if given file is a zip archive by looking for the zip-magic header 'PK\003\044'
  returns 1 if a zip file, 0 if not, -1 if an error occurred
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
/* better: if the file is a BOINC softlink to a zipped file, (boinc_resolve succeeds),
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
  FILE*fp;

  ret = boinc_resolve_filename(filename,resfilename,size);
  if (ret) {
    LogPrintf(LOG_CRITICAL,"ERROR %d boinc_resolving file '%s'\n", ret, filename);
    return(-1);
  }
  if (strncmp(filename,resfilename,size) == 0) {
    /* boinc_resove() returned the same filename, so filename wasn't a softlink */

    strncpy(buf,filename,sizeof(buf));
    strncat(buf,LINKED_EXT,sizeof(buf));
    /* f**king BOINC's new symlink behavior returns no error if the link file doesn't,
       exist, so we need to check it manually */
    if((fp=fopen(buf,"r"))) {
      fclose(fp);
      /* this could only be the remainder of a previous interrupted unzip */
      LogPrintf (LOG_NORMAL, "WARNING: found old link file '%s'\n", buf);

      /* try to resolve it (again) */
      if (boinc_resolve_filename(buf,resfilename,size))
	LogPrintf (LOG_NORMAL, "WARNING: Couldn't boinc_resolve '%s' though should be a softlink\n", buf);
	
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

  /** we end up here if boinc_resolve found the filename to be a softlink */
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
  int l;                       /**< length of matched string */
  int res = 0;                 /**< return value of a function call */
  char *startc,*endc;          /**< pointers for parsing a command-line argument */
  int output_help = 0;         /**< flag: should we write out an additional help string?
				    describing additional command-line arguments handled
			            only by this BOINC-wrapper? */
  FILE*fp;                     /**< file pointer to check only if a file can be opened */
  int breakpoint = 0;          /**< stop at breakpoint? (for testing the Windows Runtime Debugger) */
  int crash_fpu = 0;
  int test_nan  = 0;
  int test_snan = 0;
  int test_sqrt = 0;

  int resultfile_present = 0;

#ifdef _WIN32
  /* point the Windows Runtime Debugger to the Symbol Store on einstein */
  diagnostics_set_symstore("http://einstein.phys.uwm.edu/symstore");
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
    else if (0 == strncmp("--DataFiles",argv[arg],strlen("--DataFiles"))) {
      int chars = strlen("--DataFiles1=");

      /* initially allocate a buffer for "--DataFiles1=" plus MAX_PATH_LEN chars */
      rargv[rarg] = (char*)calloc(MAX_PATH_LEN + chars, sizeof(char));
      if(!rargv[rarg]){
	LogPrintf(LOG_CRITICAL, "Out of memory\n");
	boinc_finish(HIERARCHICALSEARCH_EMEM);
      }

      /* copy & skip the "[1|2]=" characters, too */
      strncpy(rargv[rarg],argv[arg],chars);
      startc = argv[arg]+chars;

      /* skip one set of single quotes if and only if they are surrounding the complete path-string */
      if ((*startc == '\'') && (*(startc+(strlen(startc)-1)) == '\'')) {
        LogPrintf (LOG_DEBUG, "DEBUG: removing quotes from path %s\n", argv[arg]);
	*(startc+strlen(startc)-1) = '\0';
	startc++;
      }

      /* look for multiple paths separated by ';' */
      while((endc = strchr(startc,';'))) {
	*endc = '\0';
	if (boinc_resolve_filename(startc,&(rargv[rarg][chars]),MAX_PATH_LEN)) {
	  LogPrintf (LOG_NORMAL, "WARNING: Can't boinc-resolve input file '%s'\n", startc);
	}

	/* append a ';' to resolved string */
	chars = strlen(rargv[rarg]) + 1;
	rargv[rarg][chars-1] =  ';';
	rargv[rarg][chars]   = '\0';

	/* make sure the next boinc_resolve() has a buffer of MAX_PATH_LEN */
	rargv[rarg] = (char*)realloc(rargv[rarg], (MAX_PATH_LEN + chars) * sizeof(char));
	if(!rargv[rarg]){
	  LogPrintf(LOG_CRITICAL, "Out of memory\n");
	  boinc_finish(HIERARCHICALSEARCH_EMEM);
	}

	/* skip the ';' in the original string */
	startc = endc+1;
      }

      /* handle last (or only) filename (comments see above) */
      if (boinc_resolve_filename(startc,&(rargv[rarg][chars]),MAX_PATH_LEN)) {
	LogPrintf (LOG_NORMAL, "WARNING: Can't boinc-resolve input file '%s'\n", startc);
      }

      /* include the terminating '\0' here */
      chars = strlen(rargv[rarg]) + 1;

#ifdef _WIN32
      /* for Windows, we have to translate the path separator '/' to '\' */
      {
	int c;
	for(c=0; c < chars; c++)
	  if(rargv[rarg][c] == '/')
	    rargv[rarg][c] = '\\';
      }
#endif
      /* truncate to the memory actually needed */
      rargv[rarg] = (char*)realloc(rargv[rarg], chars * sizeof(char));
    }

    /* handle output file */
    else if ((0 == strncmp("-o",argv[arg],strlen("-o"))) ||
	     (0 == strncmp("--fnameout",argv[arg],strlen("--fnameout")))) {

      char targetpath[MAX_PATH_LEN];
      memset(targetpath,0,MAX_PATH_LEN);

      /* these are two similar but not equal cases ("option file" and "option=file") */
      if ((0 == strncmp("-o=",argv[arg],strlen("-o="))) ||
	  (0 == strncmp("--fnameout=",argv[arg],strlen("--fnameout="))))
	{
	  int s;
	  startc = strchr(argv[arg],'=');
	  startc++; /* filename begins _after_ '=' */
	  if (boinc_resolve_filename(startc,resultfile,sizeof(resultfile))) {
	    LogPrintf (LOG_NORMAL, "WARNING: Can't boinc-resolve result file '%s'\n", startc);
	  }
#ifndef _WIN32
	  /* if boinc_resolve() returns a symbolic link, resolve it outself */
	  if (readlink(resultfile,targetpath,sizeof(targetpath)) != -1)
	    strncpy(resultfile,targetpath,sizeof(resultfile));
#endif
	  s = (startc - argv[arg]) + strlen(resultfile) + 1;
	  rargv[rarg] = (char*)calloc(s,sizeof(char));
	  if(!rargv[rarg]){
	    LogPrintf(LOG_CRITICAL, "Out of memory\n");
	    boinc_finish(HIERARCHICALSEARCH_EMEM);
	  }
	  strncpy(rargv[rarg],argv[arg], (startc - argv[arg]));
	  strncat(rargv[rarg],resultfile,s);
	}
      else /* option and argument, not joined with "=" */
	{
	  if(arg + 1 >= argc) {
	    LogPrintf(LOG_CRITICAL,"ERROR in command line: no argument following %s option\n",argv[arg]);
	    res = HIERARCHICALSEARCH_EFILE;
	  } else {
	    rargv[rarg] = argv[arg]; /* copy the "-o" */
	    arg++;                   /* grab next argument */
	    rarg++;
	    if (boinc_resolve_filename(argv[arg],resultfile,sizeof(resultfile))) {
	      LogPrintf (LOG_NORMAL, "WARNING: Can't boinc-resolve result file '%s'\n", argv[arg]);
	    }
#ifndef _WIN32
	    /* if boinc_resolve() returns a symbolic link, resolve it outself */
	    if (readlink(resultfile,targetpath,sizeof(targetpath)) != -1)
	      strncpy(resultfile,targetpath,sizeof(resultfile));
#endif
	    rargv[rarg] = resultfile;
	  }
	}
    }

    /* set the "flops estimation" */
    else if (MATCH_START("--WUfpops=",argv[arg],l)) {
      estimated_flops = atof(argv[arg]+l);
      rarg--; rargc--; /* this argument is not passed to the main worker function */
    }

    /* fire up debugger at breakpoint, solely for testing the debugger (and symbols) */
    else if (MATCH_START("--BreakPoint",argv[arg],l)) {
      breakpoint = -1;
      rarg--; rargc--; /* this argument is not passed to the main worker function */
    }

    /* drain fpu stack, solely for testing FPU exceptions */
    else if (MATCH_START("--CrashFPU",argv[arg],l)) {
      crash_fpu = -1;
      rarg--; rargc--; /* this argument is not passed to the main worker function */
    }

    /* produce a NaN, solely for testing FPU exceptions */
    else if (MATCH_START("--TestNaN",argv[arg],l)) {
      test_nan = -1;
      rarg--; rargc--; /* this argument is not passed to the main worker function */
    }

    /* produce a NaN, solely for testing FPU exceptions */
    else if (MATCH_START("--TestSNaN",argv[arg],l)) {
      test_snan = -1;
      rarg--; rargc--; /* this argument is not passed to the main worker function */
    }

    else if (MATCH_START("--TestSQRT",argv[arg],l)) {
      test_sqrt = -1;
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


#if DEBUG_COMMAND_LINE_MANGLING
  /* debug: dump the modified command line */
  {
    int i;
    fprintf(stderr,"command line:");
    for(i=0;i<rargc;i++)
      fprintf(stderr," %s",rargv[i]);
    fprintf(stderr,"\n");
  }
#endif

#ifdef __GLIBC__
  /* log the glibc version */
  fprintf(stderr, "glibc version/release: %s/%s\n", gnu_get_libc_version(), gnu_get_libc_release());
#endif

  /* if there already was an error, there is no use in continuing */
  if (res) {
    LogPrintf (LOG_CRITICAL, "ERROR: error %d in command-line parsing\n", res);
    boinc_finish(res);
  }

  /* test the debugger (and symbol loading) here if we were told to */
#ifdef _MSC_VER
  /* break on file present */
#define DEBUG_BREAKPOINT_FNAME "EAH_MSC_BREAKPOINT"
  if ((fp=fopen("..\\..\\" DEBUG_BREAKPOINT_FNAME, "r")) ||
      (fp=fopen(DEBUG_BREAKPOINT_FNAME, "r")) ) {
    fclose(fp);
    DebugBreak();
  }

  /* break on command-line option present */
  if (breakpoint)
    DebugBreak();
#elif defined(__GNUC__)
  if (breakpoint)
    attach_gdb();
#endif

  enable_floating_point_exceptions();

  if(crash_fpu)
    drain_fpu_stack();

  if(test_nan)
    fprintf(stderr,"NaN:%f\n", get_nan());

#ifdef _MSC_VER
  if(test_snan)
    fprintf(stderr,"sNaN:%f\n", get_float_snan());
#endif

  if(test_sqrt)
    fprintf(stderr,"NaN:%f\n", sqrt(-1));

  if((fp = boinc_fopen(resultfile,"r"))) {
    fclose(fp);
    resultfile_present = 1;
  }

#ifdef BOINC_APIV6
  if(setup_shmem())
    LogPrintf (LOG_NORMAL, "WARNING: Couldn't set up communication with graphics process\n");
  else
    LogPrintf (LOG_DEBUG, "Set up communication with graphics process.\n");
#endif

#ifdef EAH_DETECT_CPU_FEATURES
  detect_cpu_type();
  LogPrintf (LOG_DEBUG, "CPU type:%d\n",global_cpu_type);
#endif

  if (output_help || !resultfile_present) {
    /* CALL WORKER's MAIN()
     */
    res = MAIN(rargc,rargv);
    if (res) {
      LogPrintf (LOG_CRITICAL, "ERROR: MAIN() returned with error '%d'\n",res);
    }

#if defined(__GNUC__) && defined(__i386__)
    {
      fpuw_t fpstat = get_fpu_status();
      fprintf(stderr,"FPU status flags: ");
      PRINT_FPU_STATUS_FLAGS(fpstat);
      fprintf(stderr,"\n");
    }
#endif

    /* if the program was called for help, we write out usage for command-line options this wrapper adds to it and exit */
    if(output_help) {
      printf("Additional options the BOINC version understands:\n");
      printf("      --WUfpops         REAL     \"flops estimation\", passed to the BOINC client as the number of Flops\n");
      printf("      --BreakPoint       -       if present fire up the Windows Runtime Debugger at internal breakpoint (WIN32 only)\n");
      printf("      --CrashFPU         -       if present drain the FPU stack to test FPE\n");
      printf("      --TestNaN          -       if present raise a NaN to test FPE\n");
      printf("      --TestSQRT         -       if present try to calculate sqrt(-1) to test FPE\n");
      boinc_finish(0);
    }
  }


  /* HANDLE OUTPUT FILES
   */

  /* we'll still try to zip and send back what's left from an output file for diagnostics */
  /* in case of an error before any output was written the result will contain the link file */
  {
#define OUTPUT_EXT ".zip"
    int zipped = is_zipped(resultfile);
    if(zipped == 0) {
      int s = strlen(resultfile) + strlen(OUTPUT_EXT) + 1;
      char*zipfile = (char*)calloc(s,sizeof(char));
      strncpy(zipfile,resultfile,s);
      strncat(zipfile,OUTPUT_EXT,s);
      
      if ( boinc_zip(ZIP_IT, zipfile, resultfile) ) {
	LogPrintf (LOG_NORMAL, "WARNING: Can't zip output file '%s'\n", resultfile);
      } else if( boinc_rename(zipfile, resultfile) ) {
	LogPrintf (LOG_NORMAL, "WARNING: Couldn't rename '%s' to '%s'\n", zipfile, resultfile);
      }
    } else if (zipped < 0) {
      LogPrintf (LOG_NORMAL, "WARNING: is_zipped() couldn't open output file '%s' (%d)\n", resultfile,res);
      /* if there wasn't an error before, a problem opening the output file is serious, so signal it */
      if(res == 0)
	res = zipped;
    }
  }

  /* finally set (fl)ops count if given */
  if (estimated_flops >= 0)
    boinc_ops_cumulative( estimated_flops, 0 /*ignore IOPS*/ );

  LogPrintf (LOG_NORMAL, "done. calling boinc_finish(%d).\n",res);
  boinc_finish(res);
} /* worker() */



/**
  the main function of the BOINC App
  deals with boinc_init(_graphics) and calls the worker
*/

int main(int argc, char**argv) {
  FILE* fp_debug;
  int skipsighandler = 0;

  /* init BOINC diagnostics */
  boinc_init_diagnostics(BOINC_DIAG_DUMPCALLSTACKENABLED |
                         BOINC_DIAG_HEAPCHECKENABLED |
                         BOINC_DIAG_ARCHIVESTDERR |
                         BOINC_DIAG_REDIRECTSTDERR |
                         BOINC_DIAG_TRACETOSTDERR);

  LogSetLevel(LOG_DETAIL); /* as long as we are debugging */

  /* dummy for keeping the RCSIDs */
  HSBOINCEXTRASCRCSID = HSBOINCEXTRASHRCSID;

  LogPrintf(LOG_NORMAL, "This program is published under the GNU General Public License, version 2\n");
  LogPrintf(LOG_NORMAL, "For details see http://einstein.phys.uwm.edu/license.php\n");
  LogPrintf(LOG_NORMAL, "This Einstein@home App was built at: " __DATE__ " " __TIME__ "\n");

  /* pass argc/v to the worker via global vars */
  global_argc = argc;
  global_argv = argv;



  /* debugging support by files */

#define DEBUG_LEVEL_FNAME "EAH_DEBUG_LEVEL"
#define NO_SYNC_FNAME     "EAH_NO_SYNC"
#define DEBUG_DDD_FNAME   "EAH_DEBUG_DDD"
#define DEBUG_GDB_FNAME   "EAH_DEBUG_GDB"

  LogPrintfVerbatim (LOG_NORMAL, "\n");
  LogPrintf (LOG_NORMAL, "Start of BOINC application '%s'.\n", argv[0]);
  
  /* don't force syncing the checkpoint file if demanded */
  if ((fp_debug=fopen("../../" NO_SYNC_FNAME, "r")) || (fp_debug=fopen("./" NO_SYNC_FNAME, "r")))
    do_sync = 0;

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
    } /* DDD DEBUGGING */

  /* see if user has created a DEBUG_GDB_FNAME file: turn on debuggin using 'gdb' */
  if ((fp_debug=fopen("../../" DEBUG_GDB_FNAME, "r")) || (fp_debug=fopen("./" DEBUG_GDB_FNAME, "r")) ) 
    {
      char commandstring[256];
      char resolved_name[MAXFILENAMELENGTH];
      char *ptr;
      pid_t process_id=getpid();
      
      fclose(fp_debug);
      LogPrintf ( LOG_NORMAL, "Found '%s' file, trying debugging with 'gdb'\n", DEBUG_GDB_FNAME);
      
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
	LALSnprintf(commandstring,sizeof(commandstring),"gdb %s %d &", resolved_name ,process_id);
	system(commandstring);
	sleep(20);
      }
    } /* GDB DEBUGGING */
#endif // GNUC



  /* install signal handler */

  /* the previous boinc_init_diagnostics() call should have installed boinc_catch_signal() for
     SIGILL
     SIGABRT
     SIGBUS
     SIGSEGV
     SIGSYS
     SIGPIPE
     With the current debugging stuff now in boinc/diagnostic.C (for Windows & MacOS)
     it's probably best to leave it that way on everything else but Linux (glibc), where
     bactrace() would give messed up stacframes in the signal handler and we are
     interested in the FPU status word, too.

     NOTE: it is critical to catch SIGINT with our own handler, because a user
     pressing Ctrl-C under boinc should not directly kill the app (which is attached to the
     same terminal), but the app should wait for the client to send <quit/> and cleanly exit. 
   */

#ifdef __MINGW32__
  ExchndlSetup();
#elif _WIN32
  signal(SIGTERM, sighandler);
  if ( !skipsighandler ) {
    signal(SIGINT, sighandler);
    signal(SIGFPE, sighandler);
  }
#elif __GLIBC__
  /* this uses unsupported features of the glibc, so don't
     use the (rather portable) boinc_set_signal_handler() here */
  {
    struct sigaction sa;

    sa.sa_sigaction = (void *)sighandler;
    sigemptyset (&sa.sa_mask);
    sa.sa_flags = SA_RESTART | SA_SIGINFO;

    sigaction(SIGTERM, &sa, NULL);
    sigaction(SIGABRT, &sa, NULL);

    if ( !skipsighandler ) {
      sigaction(SIGINT,  &sa, NULL);
      sigaction(SIGSEGV, &sa, NULL);
      sigaction(SIGFPE,  &sa, NULL);
      sigaction(SIGILL,  &sa, NULL);
      sigaction(SIGBUS,  &sa, NULL);
    }
  }
#else
  /* install signal handler (generic unix) */
  boinc_set_signal_handler(SIGTERM, sighandler);
  if ( !skipsighandler ) {
      boinc_set_signal_handler(SIGINT, sighandler);
      boinc_set_signal_handler(SIGFPE, boinc_catch_signal);
  } /* if !skipsighandler */
#endif /* WIN32 */


#ifdef _NO_MSC_VER
  if (try_load_dlls(delayload_dlls, "ERROR: Failed to load %s - terminating\n")) {
    LogPrintf(LOG_NORMAL,"ERROR: Loading of mandantory DLLs failed\n");
    boinc_init();
    boinc_finish(29);
  }
#endif

#if defined(__APPLE__) && ! defined(BOINC_APIV6)
  setMacIcon(argv[0], MacAppIconData, sizeof(MacAppIconData));
#endif

  /* boinc_init */
  set_boinc_options();
  boinc_init();
  worker();
  LogPrintf (LOG_NORMAL, "done. calling boinc_finish(%d).\n",0);
  boinc_finish(0);
  /* boinc_finish() ends the program, we never get here */
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
int init_and_read_checkpoint(toplist_t*tl     , /**< the toplist to checkpoint */
			     UINT4*count,       /**< returns the skypoint counter if a checkpoint was found */
			     UINT4 total,       /**< total number of skypoints */
			     char*outputname,   /**< name of checkpointed output file */
			     char*cptname       /**< name of checkpoint file */
			     ) {
  FILE*fp;

  /* remember the toplist pointer */
  toplist = tl;

  /* store the name of the output file in global outfilename */
  {
    int s = strlen(outputname)+1;
    outfilename = (char*)calloc(s,sizeof(char));
    if(!outfilename){
      LogPrintf(LOG_CRITICAL, "Out of memory\n");
      return(-2);
    }
    strncpy(outfilename,outputname,s);
  }

  /* nothing to do if the output file already contains an end marker
     (it always exists when running under BOINC) */
  {
    int alldone = 0;
    fp=fopen(outputname,"rb");
    if(fp) {
      int len = strlen("%DONE\n");
      if(!fseek(fp,-len,SEEK_END)) {
	char *buf;
	if((buf=((char*)LALCalloc(len+1,sizeof(char))))) {
	  if(len == fread(buf,sizeof(char),len,fp))
	    if (0 == strcmp(buf,"%DONE\n"))
		alldone = -1;
	  LALFree(buf);
	}
      }
      fclose(fp);
    }
    if(alldone)
      return(2);
  }

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
    int s;
    char *c;
    c = strrchr(outputname,'/');
    if (c) {
      c++;
    } else {
      c = strrchr(outputname,'\\');
      if (c) {
	c++;
      } else {
	c = outputname;
      }
    }
    s = strlen(c)+strlen(CHECKPOINT_EXT)+1;
    cptfilename = (char*)calloc(s,sizeof(char));
    if(!cptfilename){
      LogPrintf(LOG_CRITICAL, "Out of memory\n");
      return(-2);
    }
    strncpy(cptfilename,c,s);
    strncat(cptfilename,CHECKPOINT_EXT,s);
  }
  
  return(read_hfs_checkpoint(cptfilename,toplist,count));
}


/** sets a checkpoint.
*/
void set_checkpoint (void) {
  /* make sure the exception mask isn't messed up by a badly written device driver etc.,
     so restore it periodically */
  enable_floating_point_exceptions();
  /* checkpoint every time (i.e. sky position) if FORCE_CHECKPOINTING */
#ifndef FORCE_CHECKPOINTING
  if (boinc_time_to_checkpoint())
#endif
    {
      write_hfs_checkpoint(cptfilename,toplist,last_count, do_sync);
      fprintf(stderr,"c\n");
      boinc_checkpoint_completed();
    }
}


/** finally writes a minimal (compacted) version of the toplist and cleans up
    all structures related to the toplist. After that, the toplist is invalid.
 */
void write_and_close_checkpointed_file (void) {
  write_hfs_oputput(outfilename,toplist);
}

/* Experimental and / or debugging stuff */

/** attach gdb to the running process; for debugging. */
void attach_gdb() {
#ifdef __GLIBC__
  char cmd[256];
  pid_t pid=getpid();
  snprintf(cmd, sizeof(cmd), "gdb -batch -pid %d -ex gcore --args %s", pid, global_argv[0]); 
  system(cmd);
  sleep(20);
#endif
}


/** sets the FPU control word.
    The argument should be a (possibly modified) 
    fpuw_t gotten from get_fpu_control_word() */
void set_fpu_control_word(const fpuw_t cword) {
  static fpuw_t fpucw;
  fpucw = cword;
#ifdef _MSC_VER
  __asm fldcw fpucw;
#elif defined(__GNUC__) && defined(__i386__)
  __asm("fldcw %0\n\t" : : "m" (fpucw));
#endif
}

/** returns the fpu control word */
fpuw_t get_fpu_control_word(void) {
  static fpuw_t fpucw = 0;
#ifdef _MSC_VER
  __asm fstcw fpucw;
#elif defined(__GNUC__) && defined(__i386__)
  __asm("fstcw %0\n\t" : "=m" (fpucw));
#endif
  return(fpucw);
}

/** returns the fpu status word */
fpuw_t get_fpu_status(void) {
  static fpuw_t fpusw = 0;
#ifdef _MSC_VER
  __asm fnstsw fpusw;
#elif defined(__GNUC__) && defined(__i386__)
  __asm("fnstsw %0\n\t" : "=m" (fpusw));
#endif
  return(fpusw);
}

/** sets the sse control/status word */
void set_sse_control_status(const ssew_t cword) {
  static ssew_t ssecw;
  ssecw = cword;
#ifdef _MSC_VER
  __asm ldmxcsr ssecw;
#elif defined(__GNUC__) && defined(__i386__)
  __asm("ldmxcsr %0\n\t" : : "m" (ssecw));
#endif
}

/** returns the sse control/status word */
ssew_t get_sse_control_status(void) {
  static ssew_t ssesw = 0;
#ifdef _MSC_VER
  __asm stmxcsr ssesw;
#elif defined(__GNUC__) && defined(__i386__)
  __asm("stmxcsr %0\n\t" : "=m" (ssesw));
#endif
  return(ssesw);
}

static void drain_fpu_stack(void) {
  static double dummy;
#if defined(__GNUC__) && defined(__i386__)
  __asm(
	"fstpl %0\n\t"
	"fstpl %0\n\t"
	"fstpl %0\n\t"
	"fstpl %0\n\t"
	"fstpl %0\n\t"
	"fstpl %0\n\t"
	"fstpl %0\n\t"
	"fstpl %0\n\t"
	"fstpl %0\n\t"
	: "=m" (dummy));
#elif defined(_MSC_VER)
  __asm {
    fstp dummy;
    fstp dummy;
    fstp dummy;
    fstp dummy;
    fstp dummy;
    fstp dummy;
    fstp dummy;
    fstp dummy;
    fstp dummy;
  }
#endif
}

static REAL4 get_nan(void) {
  static const UINT4 inan =
    /* 0xFFFFFFFF; /* quiet NaN */
       0xFF8001FF; /* signaling NaN palindrome */
  return((*((REAL4*)&inan)) * ((REAL4)estimated_flops));
}


void enable_floating_point_exceptions(void) {
#if defined(_MSC_VER) && 0
#define MY_INVALID 1 /* _EM_INVALID /**/
  /*
    _controlfp(MY_INVALID,_MCW_EM);
    _controlfp_s(NULL,MY_INVALID,_MCW_EM);
  */
  {
    unsigned int cw87, cwSSE;
    __control87_2(MY_INVALID,_MCW_EM,&cw87,&cwSSE);
  }
#elif defined(_MSC_VER) || defined(__GNUC__) && defined(__i386__)
  /* write out the masked FPU exceptions */
  /*
  {
    fpuw_t fpstat = get_fpu_status();
    fprintf(stderr,"FPU status flags: ");
    PRINT_FPU_STATUS_FLAGS(fpstat);
    fprintf(stderr,"\n");
  }
  */

  {
    fpuw_t fpstat;
    
    fpstat = get_fpu_control_word();
    /*
    fprintf(stderr,"FPU masked exceptions now: %4x:",fpstat);
    PRINT_FPU_EXCEPTION_MASK(fpstat);
    fprintf(stderr,"\n");
    */

    /* throw an exception at an invalid operation */
    fpstat &= ~FPU_STATUS_INVALID;
    set_fpu_control_word(fpstat);    

    /*
    fprintf(stderr,"FPU masked exceptions set: %4x:",fpstat);
    PRINT_FPU_EXCEPTION_MASK(fpstat);
    fprintf(stderr,"\n");
    */

    /* this is weird - somtimes gcc seems to cache the fpstat value
       (e.g. on MacOS Intel), so the second reading of the control_word
       doesn't seem to do what is expected
    fpstat = 0;

    fpstat = get_fpu_control_word();
    fprintf(stderr,"FPU exception mask set to:  %4x:",fpstat);
    PRINT_FPU_EXCEPTION_MASK(fpstat);
    fprintf(stderr,"\n");
    */
#ifdef ENABLE_SSE_EXCEPTIONS
    set_sse_control_status(get_sse_control_status() & ~SSE_MASK_INVALID);
#endif
  }
#endif
}

int segfault (void) {
  volatile int i = *((int*)1);
  return(i);
}


#ifdef EAH_DETECT_CPU_FEATURES

/***********************
 * cpu feature detection
 ***********************/

/* Todo:
   - check handling of AMD extensions
 */

/* ALLOW_FILE_OVERRIDE by default */
#ifndef ALLOW_FILE_OVERRIDE
#define ALLOW_FILE_OVERRIDE 1
#endif

/* for stat() */
#if ALLOW_FILE_OVERRIDE
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#endif

#include "cpu_type_features.h"

#ifndef FEAT_CPU_TYPE_1
#define FEAT_CPU_TYPE_1 0
#endif
#ifndef FEAT_CPU_TYPE_2
#define FEAT_CPU_TYPE_2 0
#endif
#ifndef FEAT_CPU_TYPE_3
#define FEAT_CPU_TYPE_3 0
#endif
#ifndef FEAT_CPU_TYPE_4
#define FEAT_CPU_TYPE_4 0
#endif
#ifndef FEAT_CPU_TYPE_5
#define FEAT_CPU_TYPE_5 0
#endif
#ifndef FEAT_CPU_TYPE_6
#define FEAT_CPU_TYPE_6 0
#endif
#ifndef FEAT_CPU_TYPE_7
#define FEAT_CPU_TYPE_7 0
#endif
#ifndef FEAT_CPU_TYPE_8
#define FEAT_CPU_TYPE_8 0
#endif
#ifndef FEAT_CPU_TYPE_9
#define FEAT_CPU_TYPE_9 0
#endif


#if defined(__POWERPC__) && defined(__APPLE__)

#include <CoreServices/CoreServices.h>

#ifndef gestaltPowerPCHas64BitSupport
#define gestaltPowerPCHas64BitSupport 6
#endif
  
/* feature flags */
#define FEAT_ALTIVEC (1 << gestaltPowerPCHasVectorInstructions)
#define FEAT_G5      (1 << gestaltPowerPCHas64BitSupport)

long get_cpu_features(void) {
  long cpuAttributes;
  OSErr oserr = Gestalt(gestaltPowerPCProcessorFeatures, &cpuAttributes);
  if (oserr)
    return(0);
  else
    return(cpuAttributes);
}


#else /* ifdef __APPLE__ ; x86 follows here */


/* feature flags */
#define FEAT_MMX    (1<<23)
#define FEAT_3DNow  (1<<31)
#define FEAT_3DNow2 (1<<30)
#define FEAT_SSE    (1<<25)
#define FEAT_SSE2   (1<<26)
#define FEAT_SSE3   (1<<10)

#ifdef _MSC_VER /* MSC, Intel assembler syntax */

#include <windows.h> 

/* this was taken from the VC++ CPUID example - I don't like this method... */
int has_CPUID() {
        __try { 
                _asm { 
                        xor eax, eax 
                        cpuid 
                } 
        } 
        __except (EXCEPTION_EXECUTE_HANDLER) { 
                return 0; 
        } 
        return 1; 
}


long get_cpu_features() {

        DWORD dwStandard = 0; 
        DWORD dwFeature = 0; 
        DWORD dwExt = 0; 

        if(!has_CPUID())
                return(0);

        _asm { 
                push ebx 
                push ecx 
                push edx 

                // get the Standard bits 
                mov eax, 1 
                cpuid 
                mov dwStandard, eax 
                mov dwFeature, edx 

                // get AMD-specials 
                mov eax, 80000000h 
                cpuid 
                cmp eax, 80000000h 
                jc notamd 
                mov eax, 80000001h 
                cpuid 
                mov dwExt, edx 

                notamd: 
                pop ecx 
                pop ebx 
                pop edx 
        } 

        return(dwFeature);
}


#else /* ifdef _MSC_VER -> gcc, AT&T assembler syntax */


/* This was mainly inspired by the CPU type detection of ffmpeg */

/* macro to get CPUID
   ebx saving is necessary for PIC. gcc seems unable to see it alone */

#ifdef SAVE_EBX_REG

#ifdef ARCH_X86_64
#define REG_b "rbx"
#define REG_S "rsi"
#else
#define REG_b "ebx"
#define REG_S "esi"
#endif
#define cpuid(index,eax,ebx,ecx,edx)\
    __asm __volatile\
        ("mov %%"REG_b", %%"REG_S"\n\t"\
         "cpuid\n\t"\
         "xchg %%"REG_b", %%"REG_S\
         : "=a" (eax), "=S" (ebx),\
           "=c" (ecx), "=d" (edx)\
         : "0" (index));

#else /* ifdef SAVE_EBX_REG */

#define cpuid(func,ax,bx,cx,dx)\
         __asm__ __volatile__\
       (\
         "cpuid\n\t"\
         :\
         "=a" (ax),\
         "=b" (bx),\
         "=c" (cx),\
         "=d" (dx)\
         :\
         "a" (func));

#endif /* ifdef SAVE_EBX_REG */


/* See if CPUID instruction is supported */
int cpu_has_cpuid() {
  long a,c;
 
  __asm__ __volatile__
    (
     /* ... Get copies of EFLAGS into eax and ecx */
     "pushf\n\t"
     "pop %0\n\t"
     "mov %0, %1\n\t"
     
     /* ... Toggle the ID bit in one copy and store */
     /*     to the EFLAGS reg */
     "xor $0x200000, %0\n\t"
     "push %0\n\t"
     "popf\n\t"
     
     /* ... Get the (hopefully modified) EFLAGS */
     "pushf\n\t"
     "pop %0\n\t"
     : "=a" (a), "=c" (c)
     :
     : "cc" 
     );

  return((a != c));
}


/* Check for SSE support of the OS.
   The OS can disable SSE use by not setting Bit 9 of CR4,
   e.g. if it doesn't save XMM registers in a context switch.
   In this case FXSAVE and FXRSTOR won't save/restore XMM registers.
   As CR4 can only be read in privileged mode, we check FXSAVE/FXRSTOR here.
   The code (including comments) is taken from Agner Fog's asmlib.
*/
int sse_disabled() {
  int ret;
  __asm__ __volatile__
    (
     "MOVL    $1,%%EAX                   \n\t" 
     "CPUID                              \n\t"
     "TEST    $0x2000000,%%EDX           \n\t" /* test SSE support */
     "JZ      noSSE%=                    \n\t"
     "TEST    $0x1000000,%%EDX           \n\t" /* test FXSAVE/FXRSTOR support */
     "JZ      noSSE%=                    \n\t"
     "SMSW    %%AX                       \n\t" /* read lower part of CR0 without violating privilege */
     "TEST    $4,%%AL                    \n\t" /* test if FXSAVE is emulated */
     "JNZ     noSSE%=                    \n\t"
     "PUSH    %%EBP                      \n\t"
     "MOV     %%ESP,%%ESI                \n\t" /* save stack pointer */
     "SUB     $0x200,%%ESP               \n\t" /* allocate space for FXSAVE */
     "AND     $-0x10,%%ESP               \n\t" /* align by 16 */
     "MOV     %%ESP,%%EBP                \n\t"
     "FXSAVE  (%%EBP)                    \n\t" /* save FP/MMX and XMM registers */
     "MOV     0x10C(%%EBP),%%EAX         \n\t" /* read part of XMM6 register */
     "XORL    $0x0D95A34BE,0x10C(%%EBP)  \n\t" /* change value to random test value */
     "FXRSTOR (%%EBP)                    \n\t" /* load changed value into XMM6 */
     "MOV     %%EAX,0x10C(%%EBP)         \n\t" /* restore old value in buffer */
     "FXSAVE  (%%EBP)                    \n\t" /* save again */
     "MOV     0x10C(%%EBP),%%EBX         \n\t" /* read changed XMM6 register */
     "MOV     %%EAX,0x10C(%%EBP)         \n\t" /* restore old value */
     "FXRSTOR (%%EBP)                    \n\t" /* load restored value into XMM6 */
     "XOR     %%EBX,%%EAX                \n\t" /* get difference between old and new value */
     "CMP     $0x0D95A34BE,%%EAX         \n\t" /* test if XMM6 was changed correctly */
     "SETE    %%AL                       \n\t" /* 1 if test successful */
     "AND     $1,%%EAX                   \n\t" /* clear rest of EAX */
     "SHL     $11,%%EAX                  \n\t" /* move to bit 11 */
     "OR      %%EAX,%%EDI                \n\t" /* combine bits in EDI */
     "MOV     %%ESI,%%ESP                \n\t" /* restore ESP */
     "POP     %%EBP                      \n\t" /* restore EBP */
     "SUB     %%EAX,%%EAX                \n\t" /* mark SSE support */
     "JMP     okSSE%=                    \n\t" 
     "noSSE%=:                           \n\t" 
     "MOV     $-1,%%EAX                  \n\t" 
     "okSSE%=:                           \n\t" 
     "MOVL    %%EAX,%0                   \n\t" /* store the result */
     : "=m" (ret)
     : 
     : "eax","ebx","ecx","edx","esi","edi"
     );
    return(ret);
}


/* Get feature flags of the cpu, return 0 (none) if not detectable
   SSE features will be reset if the OS disabled SSE */
long get_cpu_features() {
  long ax,bx,cx,dx;
  long ret;
  if (cpu_has_cpuid()) {
    cpuid(1,ax,bx,cx,dx);
    ret = dx;

    /* put the SSE3 information in an unused bit */
    if (cx & 1)
      ret |= FEAT_SSE3;
    else
      ret &= ~FEAT_SSE3;

    if (FEAT_SSE & ret)
      if (sse_disabled()) {
	fprintf(stderr,"SSE not supported by OS\n");
	ret &= ~(FEAT_SSE | FEAT_SSE2 | FEAT_SSE3);
      }
    return(ret);
  } else
    return(0);
}

#endif /* ifdef _MSC_VER */
#endif /* ifdef __POWERPC__ && __APPLE__ */



/* set a cpu type depending on features,
   caching the detection result in a static variable */
int detect_cpu_type() {
  static int t = -1;
  long features = 0;
#if ALLOW_FILE_OVERRIDE
  struct stat statb;
#endif
  if (t<0)
#if ALLOW_FILE_OVERRIDE
    if( !stat("CPU_TYPE_0",&statb) || !stat("../../CPU_TYPE_0",&statb) ) {
      t=0;
      fprintf(stderr,"Overridden CPU type %d\n",t);
    } else
#endif
    {
      features = get_cpu_features();
      if      (features & FEAT_CPU_TYPE_9)
        t=9;
      else if (features & FEAT_CPU_TYPE_8)
        t=8;
      else if (features & FEAT_CPU_TYPE_7)
        t=7;
      else if (features & FEAT_CPU_TYPE_6)
        t=6;
      else if (features & FEAT_CPU_TYPE_5)
        t=5;
      else if (features & FEAT_CPU_TYPE_4)
        t=4;
      else if (features & FEAT_CPU_TYPE_3)
        t=3;
      else if (features & FEAT_CPU_TYPE_2)
        t=2;
      else if (features & FEAT_CPU_TYPE_1)
        t=1;
      else
        t=0;
      fprintf(stderr,"Detected CPU type %d\n",t);
    }

  global_cpu_type = t;

  return(t);
}

#endif /* EAH_DETECT_CPU_FEATURES */
