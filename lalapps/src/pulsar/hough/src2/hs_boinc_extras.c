/* Extras for building an Einstein@Home BOINC App from HierarchicalSearch
   Bernd Machenschalk for Einstein@Home
*/

/* TODO:
   - cleanup of toplist if something goes wrong reading it from checkpoint
   - error handling in checkpointing
   - catch malloc errors in worker()
   - behavior when boinc_is_standlone()?
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

/* this includes patches for chdir() and sleep() */
#ifdef _WIN32
#include "win_lib.h"
#endif

#include "hs_boinc_extras.h"
#include "HierarchicalSearch.h"
#include <lal/LogPrintf.h>

NRCSID(HSBOINCEXTRASCRCSID,"$Id$");

/* probably already included by previous headers, but anyway */
#include <stdlib.h>
#include <string.h>
#if (BOINC_GRAPHICS == 2)
#include <dlfcn.h>
#endif



/** MACROS **/

#define MAX_PATH_LEN 512

/* don't want to include LAL headers for PI */
#define LAL_PI 3.1415926535897932384626433832795029  /**< pi */

#define BOINC_TRY(test,code,mess) \
  if (test) { \
    LogPrintf (LOG_CRITICAL, "ERROR %d: %s\n", code, mess); \
    boinc_finish(code); \
  }
    

/* compare strings s1 and s2 up to the length of s1 (w/o the '\0'!!)
   and set l to the length */
#define MATCH_START(s1,s2,l) (0 == strncmp(s1,s2,(l=strlen(s1))-1))




/** global VARIABLES **/

/* information about the output files */
static char **outfiles = NULL;        /* the names  of the output files of the program */
static int  noutfiles  = 0;           /* the number of the output files of the program */
static char resultfile[MAX_PATH_LEN]; /* the name of the file / zip archive to return */
static double estimated_flops = -1;

/* hooks for communication with the graphics thread */
void (*set_search_pos_hook)(float,float) = NULL;
int (*boinc_init_graphics_hook)(void (*worker)(void)) = NULL;
double *fraction_done_hook = NULL;

/* declare graphics stuff here if we don't get it from a dynamic library */
#if (BOINC_GRAPHICS == 1)
extern double fraction_done;
extern void set_search_pos(float RAdeg, float DEdeg);
extern int boinc_init_graphics(void (*worker)(void));
#endif

/* worker doesn't take arguments, so we have to pass them as (mol) global vars :-( */
static int global_argc;
static char **global_argv;

/* variables for checkpointing */
static char* cptfilename;
static FStatCheckpointFile* cptf = NULL;
static UINT4 bufsize = 8*1024;
static UINT4 maxsize = 1024*1024;
static double last_rac, last_dec;
static UINT4 last_count, last_total;




/** PROTOTYPES **/
static void worker (void);
static void sighandler(int);




/** FUNCTIONS **/


/*
  sighandler()
*/
#ifdef __GLIBC__
  /* needed to define backtrace() which is glibc specific*/
#include <execinfo.h>
#endif /* __GLIBC__ */

/* signal handlers */
static void sighandler(int sig){
  LALStatus *mystat = global_status;

#ifdef __GLIBC__
  void *array[64];
  size_t size;
#endif
  static int killcounter = 0;

  /* RP: not sure what this is for. FIXME: better remove?
#ifndef _WIN32
  sigset_t signalset;
  sigemptyset(&signalset);
  sigaddset(&signalset, sig);
  pthread_sigmask(SIG_BLOCK, &signalset, NULL);
#endif
  */


  /* lets start by ignoring ANY further occurences of this signal
     (hopefully just in THIS thread, if truly implementing POSIX threads */
  LogPrintfVerbatim(LOG_CRITICAL, "\n");
  LogPrintf (LOG_CRITICAL, "APP DEBUG: Application caught signal %d.\n\n", sig );

  /* ignore TERM interrupts once  */
  if ( sig == SIGTERM || sig == SIGINT )
    {
      killcounter ++;

      if ( killcounter >= 4 )
        {
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
  size = backtrace (array, 64);
  LogPrintf (LOG_CRITICAL,   "Obtained %zd stack frames for this thread.\n", size);
  LogPrintf (LOG_CRITICAL,   "Use gdb command: 'info line *0xADDRESS' to print corresponding line numbers.\n");
  backtrace_symbols_fd(array, size, fileno(stderr));
#endif /* __GLIBC__ */
  /* sleep a few seconds to let the OTHER thread(s) catch the signal too... */
  sleep(5);
  boinc_finish(COMPUTEFSTAT_EXIT_SIGNAL);
  return;
} /* sighandler */




/*
  show_progress()
 */
void show_progress(double rac, double dec, UINT4 count, UINT4 total) {
  double fraction = (double)count / (double)total;

  /* set globals */
  last_rac = rac;
  last_dec = dec;
  last_count = count;
  last_total = total;

  if (fraction_done_hook)
    *fraction_done_hook = fraction;
  if (set_search_pos_hook)
    set_search_pos_hook(rac * 180.0/LAL_PI, dec * 180.0/LAL_PI);
  boinc_fraction_done(fraction);
  if (estimated_flops >= 0)
    boinc_ops_cumulative( estimated_flops * fraction, 0 /*ignore IOPS*/ );
}





/*
  register_output_file()

  this registers a new output file to be zipped into the archive that is returned
  to the server as a result file
 */
void register_output_file(char*filename) {
  outfiles = (char**)realloc(outfiles,(noutfiles+1)*sizeof(char*));
  if (outfiles == NULL) {
    LogPrintf (LOG_CRITICAL, "ERROR: Can't allocate output filename '%s'\n", filename);
    noutfiles = 0;
    return;
  }
  outfiles[noutfiles] = malloc(strlen(filename));
  if (outfiles[noutfiles] == NULL) {
    LogPrintf (LOG_CRITICAL, "ERROR: Can't allocate output filename '%s'\n", filename);
    return;
  }
  strcpy(outfiles[noutfiles],filename);
  noutfiles++;
}



/* to be inmplemented... */
void set_checkpoint(char*filename,double rac,double dec, long tpl_count, long tpl_total);
void get_checkpoint(char*filename);





/*
  is_zipped()
  
  check if given file is a zip archive by looking for the zip-magic header 'PK\003\044'
  RETURN: 1 = true, 0 = false, -1 = error
 */
static int is_zipped ( const char *fname ) {
  FILE *fp;
  char zip_magic[] = {'P', 'K', 3, 4 };
  char file_header[4];

  if ( (fp = fopen( fname, "rb")) == NULL ) {
    LogPrintf (LOG_CRITICAL, "Failed to open '%s' for reading.\n", fname);
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



/*
  resolve_and_unzip()

  prepare an input file for the program, i.e. boinc_resolve and/or unzip it
 */
static int resolve_and_unzip(const char*filename, char*resfilename, const size_t size) {
  int zipped;

  if (boinc_resolve_filename(filename,resfilename,size)) {
    LogPrintf (LOG_NORMAL, "WARNING: Can't boinc-resolve file '%s'\n", filename);

    zipped = is_zipped (filename);

    if (zipped<0) {
      return(-1);

    } else if (zipped) { 

      /* unzip in-place:
	 unzip file to "filename.uz", then replace original file with "filename.uz" */
      LogPrintf (LOG_NORMAL, "WARNING: Unzipping '%s' in-place\n", filename);
      strncpy(resfilename,filename,size);
      strncat(resfilename,".uz",size);
      boinc_delete_file(resfilename);
      if( boinc_zip(UNZIP_IT,filename,resfilename) ) {
	LogPrintf (LOG_CRITICAL, "ERROR: Couldn't unzip '%s'\n", filename);
	return(-1);
      }
      if( boinc_delete_file(filename) ) {
	LogPrintf (LOG_CRITICAL, "ERROR: Couldn't delete '%s'\n", filename);
	return(-1);
      }
      if( boinc_rename(resfilename,filename) ) {
	LogPrintf (LOG_CRITICAL, "ERROR: Couldn't rename to '%s'\n", filename);
	return(-1);
      }
    }

    /* copy the filename into resfile as if boinc_resove() had succeeded */
    strncpy(resfilename,filename,size);
    return(0);
  }

  /* we end up here if boinc_resolve was successful */
  zipped = is_zipped (resfilename);

  /* return if not zipped or couldn't find out */
  if (zipped <= 0)
    return(zipped);

  /* delete the local link so we can unzip to that name */
  if( boinc_delete_file(filename) ) {
    LogPrintf (LOG_CRITICAL, "ERROR: Couldn't delete '%s'\n", filename);
    return(-1);
  }

  /* unzip */
  if ( boinc_zip(UNZIP_IT,resfilename,filename) ) {
    LogPrintf (LOG_CRITICAL, "ERROR: Couldn't unzip '%s'\n", resfilename);
    return(-1);
  }

  /* the new resolved filename is the unzipped file */
  strncpy(resfilename,filename,size);
  return(0);
}



/*
  worker()

  The worker() ist called either from main() directly or from boinc_init_graphics
  (in a separate thread). It does some funny things to the command line (mostly
  boinc-resolving filenames), then calls MAIN() (from HierarchicalSearch.c), and
  finally handles the output / result file before properly exiting with boinc_finish().
*/
static void worker (void) {
  int argc    = global_argc;   /* as worker is defined void worker(void), ... */
  char**argv  = global_argv;   /* ...  take argc and argv from global variables */
  int rargc   = global_argc;   /* argc and ... */
  char**rargv = NULL;          /* ... argv values for calling the MAIN() function of the worker */
  int arg, rarg;               /* current command-line argument */
  int i;                       /* loop counter */
  int l;                       /* length of matched string */
  int res = 0;                 /* return value of function call */
  char *startc,*endc,*appc;

  boinc_init_diagnostics(BOINC_DIAG_DUMPCALLSTACKENABLED |
                         BOINC_DIAG_HEAPCHECKENABLED |
                         BOINC_DIAG_ARCHIVESTDERR |
                         BOINC_DIAG_REDIRECTSTDERR |
                         BOINC_DIAG_TRACETOSTDERR);

  /* try to load the graphics library and set the hooks if successful */
#if (BOINC_GRAPHICS == 2) 
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
     and unresolved) and the flops estimation will be dealt with.
  */
  rargv = (char**)malloc(argc*sizeof(char*));
  rargv[0] = argv[0];
  rarg = 1;

  /* for all args in the command line (except argv[0]) */
  for (arg=1; arg<argc; arg++) {
    
    /* config file */
    if (argv[arg][0] == '@') {
      rargv[rarg] = (char*)malloc(MAX_PATH_LEN);
      rargv[rarg][0] = '@';
      if (boinc_resolve_filename(argv[arg]+1,rargv[rarg]+1,MAX_PATH_LEN-1)) {
        LogPrintf (LOG_NORMAL, "WARNING: Can't boinc-resolve config file '%s'\n", argv[arg]+1);
      }
    }

    /* skygrid file */
    else if (MATCH_START("--skyGridFile=",argv[arg],l)) {
      rargv[rarg] = (char*)malloc(MAX_PATH_LEN);
      strncpy(rargv[rarg],argv[arg],l);
      if (resolve_and_unzip(argv[arg]+l, rargv[rarg]+l, MAX_PATH_LEN-l) < 0)
	res = HIERARCHICALSEARCH_EFILE;
    }

    /* ephermeris files */
    else if (MATCH_START("--ephemE=",argv[arg],l)) {
      rargv[rarg] = (char*)malloc(MAX_PATH_LEN);
      strncpy(rargv[rarg],argv[arg],l);
      if (resolve_and_unzip(argv[arg]+l, rargv[rarg]+l, MAX_PATH_LEN-l) < 0)
	res = HIERARCHICALSEARCH_EFILE;
    }
    else if (MATCH_START("--ephemS=",argv[arg],l)) {
      rargv[rarg] = (char*)malloc(MAX_PATH_LEN);
      strncpy(rargv[rarg],argv[arg],l);
      if (resolve_and_unzip(argv[arg]+l, rargv[rarg]+l, MAX_PATH_LEN-l) < 0)
	res = HIERARCHICALSEARCH_EFILE;
    }


    /* SFT files (no unzipping, but dealing with multiple files separated by ';' */
    else if (0 == strncmp("--DataFiles",argv[arg],11)) {
      rargv[rarg] = (char*)malloc(1024);

      /* copy & skip the "[1|2]=" characters, too */
      strncpy(rargv[rarg],argv[arg],13);
      appc = rargv[rarg]+13;
      startc = argv[arg]+13;

      /* skip single quotes if and only if they are surrounding the complete path-string */
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

	/* append a ';' to resolved string */
	appc = appc + strlen(appc) + 1;
	*(appc-1) = ';';
	*appc = '\0';

	/* skip the ';' in the original string */
	startc = endc+1;
      }

      /* handle last (or only) filename */
      if (boinc_resolve_filename(startc,appc,255)) {
	LogPrintf (LOG_NORMAL, "WARNING: Can't boinc-resolve input file '%s'\n", startc);
      }
    }

    /* output file */
#define OUTPUT_EXT ".res"
    else if (MATCH_START("--fnameout=",argv[arg],l)) {
      int s;
      if (boinc_resolve_filename(argv[arg]+l,resultfile,sizeof(resultfile))) {
        LogPrintf (LOG_NORMAL, "WARNING: Can't boinc-resolve result file '%s'\n", argv[arg]+l);
      }
      s = strlen(argv[arg])+strlen(OUTPUT_EXT)+1;
      rargv[rarg] = (char*)malloc(s);
      strncpy(rargv[rarg],argv[arg],s);
      strncat(rargv[rarg],OUTPUT_EXT,s);
      register_output_file(rargv[rarg]+l);
    }
    else if (0 == strncmp("-o",argv[arg],strlen("-o"))) {
      int s;
      rargv[rarg] = argv[arg]; /* copy the "-o" */
      arg++;                /* grab next argument */
      rarg++;
      if(arg >= argc) {
	LogPrintf(LOG_CRITICAL,"ERROR in command line: no argument following '-o' option\n");
	res = HIERARCHICALSEARCH_EFILE;
      } else {
	if (boinc_resolve_filename(argv[arg],resultfile,sizeof(resultfile))) {
	  LogPrintf (LOG_NORMAL, "WARNING: Can't boinc-resolve result file '%s'\n", argv[arg]);
	}
	s = strlen(argv[arg])+strlen(OUTPUT_EXT)+1;
	rargv[rarg] = (char*)malloc(s);
	strncpy(rargv[rarg],argv[arg],s);
	strncat(rargv[rarg],OUTPUT_EXT,s);
	register_output_file(rargv[rarg]);
      }
    }

    /* flops estimation */
    else if (MATCH_START("--WUfpops=",argv[arg],l)) {
      estimated_flops = atof(argv[arg]+l);
      rarg--; rargc--; /* this argument is not passed to the main worker function */
    }

    /* any other argument */
    else 
      rargv[rarg] = argv[arg];

    rarg++;
  } /* for all command line arguments */



  /* sanity check */
  if (!resultfile) {
      LogPrintf (LOG_CRITICAL, "ERROR: no result file has been specified\n");
      res = HIERARCHICALSEARCH_EFILE;
  }

  /* debug */
  for(i=0;i<rargc;i++)
    LogPrintf (LOG_DETAIL, "DETAIL: command-line argument %d: %s\n", i,rargv[i]);

  /* if there already was an error, there is no use in continuing */
  if (res) {
    LogPrintf (LOG_CRITICAL, "ERROR: error %d in command-line parsing\n", res);
    boinc_finish(res);
  }


  /* CALL WORKER's MAIN()
   */

  res = MAIN(rargc,rargv);
  if (res) {
    LogPrintf (LOG_CRITICAL, "ERROR: MAIN() returned with error '%d'\n",res);
  }

  /* we'll still try to zip and send back what's left from an output file for diagnostics */
  /* in case of an error before any output was written the result will contain the link file */


  /* HANDLE OUTPUT FILES
   */
  if(noutfiles == 0)
    LogPrintf (LOG_CRITICAL, "ERROR: no output file has been specified\n");

  for(i=0;i<noutfiles;i++)
    if ( 0 == strncmp(resultfile, outfiles[i], sizeof(resultfile)) )
      LogPrintf (LOG_NORMAL, "WARNING: output (%d) and result file are identical (%s) - output not zipped\n", i, resultfile);
    else if ( boinc_zip(ZIP_IT, resultfile, outfiles[i]) ) {
      LogPrintf (LOG_NORMAL, "WARNING: Can't zip output file '%s'\n", outfiles[i]);
    }

  /* finally set (fl)ops count if given */
  if (estimated_flops >= 0)
    boinc_ops_cumulative( estimated_flops, 0 /*ignore IOPS*/ );

  boinc_finish(res);
} /* worker() */





/*
  main()

  the main function of the BOINC App
  deals with boinc_init(_graphics) and calls the worker
*/

int main(int argc, char**argv) {
  FILE* fp_debug;
  int skipsighandler = 0;

  LogSetLevel(LOG_DETAIL); /* as long as we are debugging */

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

#if BOINC_GRAPHICS>0
  {
    int retval;
#if BOINC_GRAPHICS==2
    /* Try loading screensaver-graphics as a dynamic library.  If this
       succeeds then extern void* graphics_lib_handle is set, and can
       be used with dlsym() to resolve symbols from that library as
       needed. */
    retval = boinc_init_graphics_lib(worker, argv[0]);
#endif /* BOINC_GRAPHICS==2 */
#if BOINC_GRAPHICS==1
    /* if we don't get them from the shared library, use variables local to here */
    set_search_pos_hook = set_search_pos;
    fraction_done_hook = &fraction_done;
    /* no dynamic library, just call boinc_init_graphics() */
    retval = boinc_init_graphics(worker);
#endif /* BOINC_GRAPHICS==1 */
    LogPrintf (LOG_CRITICAL,  "boinc_init_graphics[_lib]() returned %d. This indicates an error...\n", retval);
    boinc_finish(HIERARCHICALSEARCH_EWORKER );
  }
#endif /*  BOINC_GRAPHICS>0 */
    
  /* we end up hereo only if BOINC_GRAPHICS == 0 or a call to boinc_init_graphics failed */
  boinc_init();
  worker();
  boinc_finish( HIERARCHICALSEARCH_ENORM );
  /* boinc_init_graphics() or boinc_finish() ends the program, we never get here */
  return(0);
}



/* inits checkpointing and read a checkpoint if already there */
void init_and_read_checkpoint(toplist_t*toplist, UINT4*count,
			      UINT4 total, char*outputname, char*cptname) {
  FILE*fp;
  unsigned int  checksum, bytes;
  unsigned long count_read, total_read;
  int scanned;

  fstat_cpt_file_create (&cptf, outputname, bufsize, maxsize, toplist);

  fstat_cpt_file_open (cptf);

  /* create checkpoint file name if necc. */
  if(!cptname) {
#define CHECKPOINT_EXT ".cpt"
    int s = strlen(outputname)+strlen(CHECKPOINT_EXT)+1;
    cptfilename = (char*)malloc(s);
    strncpy(cptfilename,outputname,s);
    strncat(cptfilename,CHECKPOINT_EXT,s);
  }

  fp = fopen(cptfilename,"r");
  
  if (!fp) {
    LogPrintf (LOG_DEBUG,  "Couldn't open checkpoint - startng over\n");
    return;
  }

  LogPrintf (LOG_DEBUG,  "Found checkpoint - reading...\n");

  scanned = fscanf(fp,"%lf,%lf,%lu,%lu,%u,%u\n",
		  &last_rac, &last_dec,
		  &count_read, &total_read,
		  &checksum, &bytes);

  fclose(fp);

  if (scanned != 6) {
    LogPrintf (LOG_DEBUG,  "Error reading checkpoint - startng over\n");
    return;
  }

  if (total_read != total) {
    LogPrintf (LOG_DEBUG,  "Error reading checkpoint - startng over\n");
    return;
  }

  LogPrintf (LOG_DEBUG,  "Read checkpoint - reading previous output...\n");

  if (0 > fstat_cpt_file_read (cptf, checksum, bytes)) {
    LogPrintf (LOG_DEBUG,  "Error reading previous output - startng over\n");
    return;
  }

  *count = count_read;
}



/* set_checkpoint() */
int add_candidate_and_checkpoint (toplist_t*toplist, FstatOutputEntry cand) {
  int ret = fstat_cpt_file_add (cptf, cand);
  if (boinc_time_to_checkpoint()) {
    FILE* fp;
    fstat_cpt_file_flush (cptf);
    fp = fopen(cptfilename,"w");
    fprintf(fp,"%f,%f,%d,%d,%d,%d\n",
	    last_rac, last_dec, last_count, last_total,
	    cptf->checksum, cptf->bytes);
    fclose(fp);
    boinc_checkpoint_completed();
  }
  return (ret);
}


void write_and_close_checkpointed_file (void) {
  fstat_cpt_file_close(cptf);
  fstat_cpt_file_destroy(&cptf);
}
