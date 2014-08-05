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

/**
 * \file
 * \ingroup pulsarApps
 * \author Bernd Machenschalk, Reinhard Prix
 */

/* Extras for building an Einstein@Home BOINC App from HierarchicalSearch
*/

/* TODO:
   - behavior when boinc_is_standlone()?
   - check for critical sections
*/


#if defined(__i386__) || defined(__x86_64__)
#define __X86__ 1
#ifdef __GNUC__
#define __GNUX86__ 1
#endif
#else
#undef __X86__
#undef __GNUC__
#endif

/** INCLUDES **/

/* BOINC includes - need to be before the #defines in hs_boinc_extras.h */
#include "boinc/boinc_api.h"
#include "boinc/diagnostics.h"
#ifdef HAVE_BOINC_ZIP
#include "boinc/boinc_zip.h"
#endif
#include "boinc/svn_version.h"
/* this ultimately needs to be fixed in boinc_api.h,
   #include "app_ipc.h" must be moved outside the C++ section */
extern int boinc_resolve_filename(const char*, char*, int len);

/* our own win_lib includes patches for chdir() and sleep() */
#ifdef _WIN32
#include "win_lib.h"
#endif

/* probably already included by previous headers, but make sure they are included */
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>

/* for finding out and logging the glibc version */
#ifdef __GLIBC__
#define _GNU_SOURCE
#include <gnu/libc-version.h>
#endif

#ifdef HAVE_BUILD_INFO_H
#include "build_info.h"
#endif

/* try to dlopen("libgcc_s.so.1") */
#ifdef DLOPEN_LIBGCC
#include <dlfcn.h>
#endif

/* our own exception handler / runtime debugger */
#if HAVE_EXCHNDL
#include "exchndl.h"
#endif

/* headers of our own code */
#include <lal/LogPrintf.h>
#include "../HierarchicalSearch.h"
#include "hs_boinc_extras.h"
#include "hs_boinc_options.h"

/* for Linux extended backtrace */
#if defined(__GLIBC__) && defined(__i386__) && defined(EXT_STACKTRACE)
#include "erp_execinfo_plus.h"
#endif

/* FIXME: we should probably eliminate the references to this */
#include "ComputeFStatistic.h"

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

typedef enum gdbcmd { gdb_dump_core, gdb_attach } gdb_cmd;

/**
 * compare strings s1 and s2 up to the length of s1 (without the trailing 0!!)
 * and set l to the length
 */
#define MATCH_START(s1,s2,l) (0 == strncmp(s1,s2,(l=strlen(s1))-1))

/** write the FPU status flags / exception mask bits to stderr */
#define PRINT_FPU_EXCEPTION_MASK(fpstat) \
  if (fpstat & FPU_STATUS_PRECISION)	 \
    fputs(" PRECISION",stderr);		 \
  if (fpstat & FPU_STATUS_UNDERFLOW)	 \
    fputs(" UNDERFLOW",stderr);		 \
  if (fpstat & FPU_STATUS_OVERFLOW)	 \
    fputs(" OVERFLOW",stderr);		 \
  if (fpstat & FPU_STATUS_ZERO_DIVIDE)	 \
    fputs(" ZERO_DIVIDE",stderr);	 \
  if (fpstat & FPU_STATUS_DENORMALIZED)	 \
    fputs(" DENORMALIZED",stderr);	 \
  if (fpstat & FPU_STATUS_INVALID)	 \
    fputs(" INVALID",stderr)

#define PRINT_FPU_STATUS_FLAGS(fpstat) \
  if (fpstat & FPU_STATUS_COND_3)      \
    fputs(" COND_3",stderr);	       \
  if (fpstat & FPU_STATUS_COND_2)      \
    fputs(" COND_2",stderr);	       \
  if (fpstat & FPU_STATUS_COND_1)      \
    fputs(" COND_1",stderr);	       \
  if (fpstat & FPU_STATUS_COND_0)      \
    fputs(" COND_0",stderr);	       \
  if (fpstat & FPU_STATUS_ERROR_SUMM)  \
    fputs(" ERR_SUMM",stderr);	       \
  if (fpstat & FPU_STATUS_STACK_FAULT) \
    fputs(" STACK_FAULT",stderr);      \
  PRINT_FPU_EXCEPTION_MASK(fpstat)

#ifndef MIN
#define MIN(a,b) (((a)<(b))?(a):(b))
#endif

static char* myultoa2(unsigned long n, char*buf, size_t size, size_t mindig) {
  int i;
  memset(buf,'0',size);
  buf[size-1] = '\0';
  for(i=size-2; i>=0; i--) {
    buf[i] = n % 10 + '0';
    n /= 10;
    if (!n)
      break;
  }
  if (i > 0)
    return MIN(buf + i, buf + (size - mindig - 1));
  return buf;
}

static char* myultoa(unsigned long n, char*buf, size_t size) {
  return myultoa2(n, buf, size, 1);
}

static char* myltoa(long n, char*buf, size_t size) {
  int i, m=0;
  memset(buf,'\0',size);
  if (n<0) {
    n = -n;
    m = -1;
  }
  for(i=size-2; i>=0; i--) {
    buf[i] = n % 10 + '0';
    n /= 10;
    if (!n)
      break;
  }
  if(m && i>0) {
    i--;
    buf[i] = '-';
  }
  if (i > 0)
    return buf + i;
  return buf;
}

/* my own time function - not using vnprintf() */
void mytime(void) {
  char buf[64];
  struct timeval tv;
  struct tm *tmv;
  if(gettimeofday(&tv,NULL)) {
    fputs("Couldn't gettimeofday()", stderr);
    return;
  }
  tmv=localtime(&tv.tv_sec);
  if (!tmv) {
    fputs("Couldn't get localtime(gettimeofday))", stderr);
    return;
  }
  fputs(myultoa(tmv->tm_year+1900, buf, sizeof(buf)), stderr);
  fputc('-', stderr);
  fputs(myultoa2(tmv->tm_mon+1, buf, sizeof(buf), 2), stderr);
  fputc('-', stderr);
  fputs(myultoa2(tmv->tm_mday, buf, sizeof(buf), 2), stderr);
  fputc(' ', stderr);
  fputs(myultoa2(tmv->tm_hour, buf, sizeof(buf), 2), stderr);
  fputc(':', stderr);
  fputs(myultoa2(tmv->tm_min, buf, sizeof(buf), 2), stderr);
  fputc(':', stderr);
  fputs(myultoa2(tmv->tm_sec, buf, sizeof(buf), 2), stderr);
  fputc('.', stderr);
  fputs(myultoa(tv.tv_usec, buf, sizeof(buf)), stderr);
}

/*^* global VARIABLES *^*/

/** the cpu type, see cpu_type_features.h */
int global_cpu_type;


/** output filename - probably not needed to be public anymore */
static char resultfile[MAX_PATH_LEN]; /**< the name of the file / zip archive to return */


/** GPU (CUDA or OpenCL) device id */
#if USE_OPENCL_KERNEL || defined(USE_CUDA)
extern int gpu_device_id;
#endif


/**
 * FLOPS estimation - may be set by command line option --WUfpops=.
 * When set, ((skypoint_counter / total_skypoints) * estimated_flops) is periodically
 * reported to the BOINC Client as the number of flops, so that together with information
 * from the Workunit Genrator, Scheduler and Validator leads to claiming the Credit that
 * the system intends to grant for a Workunit
 */
static double estimated_flops = -1.0;


/** worker() doesn't take arguments, so we have to pass it argv/c as global vars :-( */
int global_argc;
char **global_argv;


/** variables for checkpointing */
static char* cptfilename;                 /**< name of the checkpoint file */
static char* outfilename;                 /**< name of the output file */
static toplist_t* toplist;                /**< the toplist we're checkpointing */
static UINT4 last_count, last_total;      /**< last template count, see last_rac */
static BOOLEAN do_sync = -1;              /**< sync checkpoint file to disk, default: yes */


/** record whether loading libgcc_s.so.1 succeeded */
static int libgcc_s_loaded = 0;


/** record the status of the last boinc_finish call in case boinc_finish() throws a signal */
static int boinc_finish_status = 0;


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

#ifdef __GNUC__
void run_gdb(gdb_cmd command);
#endif

void ReportStatus(LALStatus *status);

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
/**
 * Attempt to load the dlls that are required to display graphics.
 * returns 0 if successful, -1 in case of a failure.
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


/**
 * LAL's REPORTSTATUS just won't work with any of NDEBUG or
 * LAL_NDEBUG set, so we write our own function that dumps the LALStatus
 * based on LogPrintf()
 */
void ReportStatus(LALStatus *status)
{ 
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
int BOINC_LAL_ErrHand (LALStatus  *status,
		       const char *func,
		       const char *file,
		       const int line,
		       volatile const char *id) {
  if (status->statusCode) {
    fprintf(stderr,
            "Level 0: %s\n"
            "\tFunction call `%s' failed.\n"
            "\tfile %s, line %d\n",
            id, func, file, line );
    ReportStatus(status);
    LogPrintf (LOG_CRITICAL, "BOINC_LAL_ErrHand(): now calling boinc_finish()\n");
    boinc_finish(boinc_finish_status= COMPUTEFSTAT_EXIT_LALCALLERROR+status->statusCode );
  }
  /* should this call boinc_finish too?? */
  return 0;
}


/**
 * our own signal handler
 */
#ifdef __GLIBC__
  /* needed to define backtrace() which is glibc specific*/
#include <signal.h>
#include <execinfo.h>
/* get REG_EIP from ucontext.h, see
   http://www.linuxjournal.com/article/6391 */
#define __USE_GNU
#include <ucontext.h>
/* see https://bugzilla.redhat.com/show_bug.cgi?id=234560 */
#ifdef __x86_64__
#ifdef  REG_RIP
#define REG_EIP REG_RIP
#endif
#define FP_SW swd
#else  /*  __x86_64__ */
#define FP_SW sw
#endif /*  __x86_64__ */
#endif /* __GLIBC__ */

/* Use thread-local storage only where needed & supported */
#if defined(__GLIBC__) && defined(__GNUC__)
#define THREAD_LOCAL static __thread
#else
#define THREAD_LOCAL static
#endif

#ifdef __GLIBC__
static void sighandler(int sig,
		       siginfo_t *info,
		       void *secret)
#else
static void sighandler(int sig)
#endif /* __GLIBC__ */
{
  static char buf[32];
  static int killcounter = 0;
#ifdef __GLIBC__
  /* for glibc stacktrace */
  static void *stackframes[64];
  static size_t nostackframes;
  static char **backtracesymbols = NULL;
  ucontext_t *uc = (ucontext_t *)secret;
#endif
  THREAD_LOCAL int boinc_finish_in_sighabdler = 0;

  fputc('\n',stderr);
  mytime();
  fputs("\n-- signal handler called: signal ",stderr);
  fputs(myultoa(sig, buf, sizeof(buf)), stderr);
  fputc('\n',stderr);

  /* ignore INT interrupts once  */
  if ( sig == SIGINT ) {
    killcounter ++;
    if ( killcounter >= 4 ) {
      fputs("App got 4th SIGINT, guess you mean it.\nCalling boinc_finish().\n",stderr);
      boinc_finish(boinc_finish_status=COMPUTEFSTAT_EXIT_USER);
    }
    else
      return;
  } /* termination signals */

  /* A SIGABRT most likely came from a failure to load libgcc_s.so.1,
     which is required for boinc_finish() (calling pthread_exit() calling
     pthread_cancel()) to work properly. In this case take the "emergency
     exit" with exit status 0 - the worst that can happen is that
     the tasks ends up with "too many exits" error. */
  if ( ( libgcc_s_loaded == -1 ) && ( sig == 6 ) ) {
    fputs("Program received SIGABRT probably because libgcc_s.so.1 wasn't loaded - trying exit(0)\n", stderr);
    /* sleep a few seconds to let the OTHER thread(s) catch the signal too... */
    sleep(5);
    exit(boinc_finish_status);
  }

#ifdef __GLIBC__
#ifdef __X86__
  /* in case of a floating-point exception write out the FPU status */
  if ( sig == SIGFPE ) {
    fpuw_t fpstat = uc->uc_mcontext.fpregs->FP_SW;
    fputs("FPU status word: ",stderr);
    fputs(myultoa(fpstat, buf, sizeof(buf)), stderr);
    fputs(", flags: ",stderr);
    PRINT_FPU_STATUS_FLAGS(fpstat);
    fputs("\n",stderr);
  }
#endif /* __X86__ */
  /* now get TRUE stacktrace */
  nostackframes = backtrace (stackframes, 64);
  if (nostackframes == 0) {
    fputs("no stack frames obtained for this thread:\n", stderr);
  } else {
    fputs(myltoa(nostackframes, buf, sizeof(buf)), stderr);
    fputs(" stack frames obtained for this thread:\n", stderr);
    /* overwrite sigaction with caller's address */
#ifndef REWRITE_STACK_SIGACTION
#define REWRITE_STACK_SIGACTION 1
#endif
#if REWRITE_STACK_SIGACTION && defined(REG_EIP)
    stackframes[1] = (void *) uc->uc_mcontext.gregs[REG_EIP];
#endif
    fputs("Use gdb command: 'info line *0xADDRESS' to print corresponding line numbers.\n",stderr);
    backtrace_symbols_fd(stackframes, nostackframes, fileno(stderr));
#if defined(__i386__) && defined(EXT_STACKTRACE)
    fputs("Trying extended stacktrace:\n", stderr);
    backtracesymbols = backtrace_symbols(stackframes, nostackframes);
    if(backtracesymbols != NULL) {
      backtrace_symbols_fd_plus((const char *const *)(void**)backtracesymbols, nostackframes, fileno(stderr));
      free(backtracesymbols);
    }
#endif /* EXT_STACKTRACE */
    fputs("\nEnd of stacktrace\n",stderr);
  }
#endif /* __GLIBC__ */

  if (global_status)
    fputs("Stack trace of LAL functions in worker thread:\n", stderr);
  while (global_status) {
    if(global_status->function)
      fputs(global_status->function, stderr);
    else
      fputs("global_status->function=NULL", stderr);
    fputs(" at ", stderr);
    if(global_status->file)
      fputs(global_status->file, stderr);
    else
      fputs("global_status->file=NULL", stderr);
    fputc(':', stderr);
    fputs(myultoa(global_status->line, buf, sizeof(buf)), stderr);
    fputc('\n', stderr);
    if (!(global_status->statusPtr)) {
      const char *p=global_status->statusDescription;
      fputs("At lowest level status code = ", stderr);
      fputs(myltoa(global_status->statusCode, buf, sizeof(buf)), stderr);
      fputs(": ", stderr);
      fputs(p?p:"NO LAL ERROR REGISTERED", stderr);
      fputc('\n', stderr);
    }
    global_status=global_status->statusPtr;
  }
  
  /* sleep a few seconds to let the OTHER thread(s) catch the signal too... */
  sleep(5);
  /* if boinc_finish() was already called from the signal handler and we end up
     here again, then boinc_finish() itself causes a signal. Prevent an endless
     loop by callin plain exit() */
  if(boinc_finish_in_sighabdler) {
    fputs("boinc_finish already called from signal handler, trying exit()\n", stderr);
    exit(boinc_finish_status);
  }
  boinc_finish_in_sighabdler = 1;
  boinc_finish(boinc_finish_status=COMPUTEFSTAT_EXIT_SIGNAL + sig);
  return;
} /* sighandler */




/**
 * show_progress() just sets some variables,
 * so should be pretty fast and can be called several times a second
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

  /* tell APIv6 graphics about status */
  boincv6_progress.skypos_rac      = rac;
  boincv6_progress.skypos_dec      = dec;
  boincv6_progress.frequency       = freq;
  boincv6_progress.bandwidth       = fband;
#ifndef HIERARCHSEARCHGCT /* used for Hough HierarchicalSearch, not GCT */
  if(toplist->elems > 0) {
    /* take the last (rightmost) leaf of the heap tree - might not be the
       "best" candidate, but for the graphics it should be good enough */
    HoughFStatOutputEntry *line = (HoughFStatOutputEntry*)(toplist->heap[toplist->elems - 1]);

    boincv6_progress.cand_frequency  = line->Freq;
    boincv6_progress.cand_spindown   = line->f1dot;
    boincv6_progress.cand_rac        = line->Alpha;
    boincv6_progress.cand_dec        = line->Delta;
    boincv6_progress.cand_hough_sign = line->HoughFStat;
  }
#endif /* used for Hough HierarchicalSearch, not GCT */
}





/**
 * check if given file is a zip archive by looking for the zip-magic header 'PK\003\044'
 * returns 1 if a zip file, 0 if not, -1 if an error occurred
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
    LogPrintf (LOG_DEBUG, "Failed to read first 4 bytes from '%s'.\n", fname);
    fclose(fp);
    return 0;
  }
  fclose(fp);

  if ( memcmp ( file_header, zip_magic, 4 ) )
    return 0;	/* false: no zip-file */
  else
    return 1;	/* yep, found magic zip-header */
} /* is_zipped() */



/**
 * prepare an input file for the program, i.e. boinc_resolve and/or unzip it if necessary
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

#ifndef HAVE_BOINC_ZIP
      LogPrintf (LOG_CRITICAL, "ERROR: would unzip '%s' if had boinc_zip\n", resfilename);
      return(-1);
#else
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
#endif
    }

    zipped = is_zipped (filename);

    if (zipped<0) {
      LogPrintf (LOG_DEBUG, "ERROR: Couldn't open '%s'\n", filename);
      return(-1);

    } else if (zipped) {

#ifndef HAVE_BOINC_ZIP
      LogPrintf (LOG_CRITICAL, "ERROR: would unzip '%s' if had boinc_zip\n", filename);
      return(-1);
#else
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
#endif

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

#ifndef HAVE_BOINC_ZIP
  LogPrintf (LOG_CRITICAL, "ERROR: would unzip '%s' if had boinc_zip\n", resfilename);
  return(-1);
#else
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
#endif
}



/**
 * The worker() ist called either from main() directly or from boinc_init_graphics
 * (in a separate thread). It does some funny things to the command line (mostly
 * boinc-resolving filenames), then calls MAIN() (from HierarchicalSearch.c), and
 * finally handles the output / result file(s) before exiting with boinc_finish().
 */
/**
 * rules for "bundled" workunits
 * - the command-line needs to begin with --BundleSize=<no_wus>
 * - each WU within a bundle must have its own config file specified on the command line
 *   the number of config files on the command-line must equal the bundle_size
 * - the config file should not refer to any files, as filenames in the config files are noy 'boinc_resolve'd
 *   in particular skygrid- and SFT files need to be specified on the command-line, thus being the
 *   same for all WUs of a bundle
 * - only a single output file (toplist) per WU is allowed
 *   the output file must be specified as a separate argument to the '-o' option,
 *   e.g. '-o outputfile', NOT '--OutputFile=outputfile'
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
  int output_version = 0;      /**< flag: version requested? This skips a check for an output file option */
  FILE*fp;                     /**< file pointer to check only if a file can be opened */
  int breakpoint = 0;          /**< stop at breakpoint? (for testing the Windows Runtime Debugger) */
  int crash_fpu = 0;
  int test_nan  = 0;
  int test_snan = 0;
  int test_sqrt = 0;

  int bundle_size = 0;           /**< how many workunits are bundled */
  int current_config_file = 0;   /**< number of current WU in bundle (or config file on command line during parsing) */
  char**config_files = NULL;     /**< list of (boinc-resolved) config files, omne for each WU in bundle */
  char**config_file_arg = NULL;  /**< points to a placeholder for the config file in the command-line paased to MAIN() */
  char wu_result_file[MAX_PATH_LEN];

  int second_outfile = 0;        /**< flag: is there a second output file, i.e. --SortToplist=3 */
  int resultfile_present = 0;

  resultfile[0] = '\0';

#ifdef _MSC_VER
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
    boinc_finish(boinc_finish_status=HIERARCHICALSEARCH_EMEM);
  }

  /* the program name (argv[0]) remains the same in any case */
  rargv[0] = argv[0];
  rarg = 1;

  /* for all args in the command line (except argv[0]) */
  for (arg=1; arg<argc; arg++) {
    
    /* if this is a "workunit bundle", --BundleSize must be the first command-line option */
    if (MATCH_START("--BundleSize=",argv[arg],l)) {
      bundle_size = atoi(argv[arg]+l);
      rarg--; rargc--; /* this argument is not passed to the main worker function */
    }

    /* a possible config file is boinc_resolved, but filenames contained in it are not! */
    else if (argv[arg][0] == '@') {
      rargv[rarg] = (char*)calloc(MAX_PATH_LEN,sizeof(char));
      if(!rargv[rarg]){
	LogPrintf(LOG_CRITICAL, "Out of memory\n");
	boinc_finish(boinc_finish_status=HIERARCHICALSEARCH_EMEM);
      }
      rargv[rarg][0] = '@';
      if (boinc_resolve_filename(argv[arg]+1,rargv[rarg]+1,MAX_PATH_LEN-1)) {
        LogPrintf (LOG_NORMAL, "WARNING: Can't boinc-resolve config file '%s'\n", argv[arg]+1);
      }
      if (bundle_size) {
	config_files = realloc(config_files, sizeof(char*) * (current_config_file + 1));
	if(!config_files){
	  LogPrintf(LOG_CRITICAL, "Out of memory\n");
	  boinc_finish(boinc_finish_status=HIERARCHICALSEARCH_EMEM);
	}
	config_files[current_config_file] = rargv[rarg];
	if (current_config_file) {
	  rarg--; rargc--; /* skip that config file */
	} else {
	  config_file_arg = &rargv[rarg]; /* save the place of the first config file in the command-line */
	}
	current_config_file ++;
      }
    }

    /* boinc_resolve and unzip skygrid file */
    else if (MATCH_START("--skyGridFile=",argv[arg],l)) {
      rargv[rarg] = (char*)calloc(MAX_PATH_LEN,sizeof(char));
      if(!rargv[rarg]){
	LogPrintf(LOG_CRITICAL, "Out of memory\n");
	boinc_finish(boinc_finish_status=HIERARCHICALSEARCH_EMEM);
      }
      strncpy(rargv[rarg],argv[arg],l);
      if (resolve_and_unzip(argv[arg]+l, rargv[rarg]+l, MAX_PATH_LEN-l) < 0)
	res = HIERARCHICALSEARCH_EFILE;
    }

    /* boinc_resolve and unzip segment list */
    else if (MATCH_START("--segmentList=",argv[arg],l)) {
      rargv[rarg] = (char*)calloc(MAX_PATH_LEN,sizeof(char));
      if(!rargv[rarg]){
	LogPrintf(LOG_CRITICAL, "Out of memory\n");
	boinc_finish(boinc_finish_status=HIERARCHICALSEARCH_EMEM);
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
	boinc_finish(boinc_finish_status=HIERARCHICALSEARCH_EMEM);
      }
      strncpy(rargv[rarg],argv[arg],l);
      if (resolve_and_unzip(argv[arg]+l, rargv[rarg]+l, MAX_PATH_LEN-l) < 0)
	res = HIERARCHICALSEARCH_EFILE;
    }
    else if (MATCH_START("--ephemS=",argv[arg],l)) {
      rargv[rarg] = (char*)calloc(MAX_PATH_LEN,sizeof(char));
      if(!rargv[rarg]){
	LogPrintf(LOG_CRITICAL, "Out of memory\n");
	boinc_finish(boinc_finish_status=HIERARCHICALSEARCH_EMEM);
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
	boinc_finish(boinc_finish_status=HIERARCHICALSEARCH_EMEM);
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
	  boinc_finish(boinc_finish_status=HIERARCHICALSEARCH_EMEM);
	}

	/* put back the ';' in the original string and skip it for next iteration */
	*endc = ';';
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
	    boinc_finish(boinc_finish_status=HIERARCHICALSEARCH_EMEM);
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
	    if (bundle_size) {
	      rargv[rarg] = wu_result_file; /* this will get copied from resultfile later */
	    } else {
	      rargv[rarg] = resultfile;
	    }
	  }
	}
    }

    /* record if there will be a second output file */
    else if (!strcmp("--SortToplist=3",argv[arg])) {
      second_outfile = -1;
    }

    /* set the "flops estimation" */
    else if (MATCH_START("--WUfpops=",argv[arg],l)) {
      estimated_flops = atof(argv[arg]+l);
      rarg--; rargc--; /* this argument is not passed to the main worker function */
    }

#if USE_OPENCL_KERNEL || defined(USE_CUDA)
    /* if building a CUDA App, handle --device option given by BOINC client */
    else if (MATCH_START("--device",argv[arg],l)) {
      arg++; /* next argument */
      gpu_device_id = atoi(argv[arg]);
      rarg-=2; rargc-=2; /* these arguments are not passed to the main worker function */
    }
#endif

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
	     (0 == strncmp("-h",argv[arg],strlen("-h")))) {
      output_help = 1;
      rargv[rarg] = argv[arg];
    }

    else if (0 == strncmp("--version",argv[arg],strlen("--version"))) {
      output_version = 1;
      rargv[rarg] = argv[arg];
    }

    /* any other argument - simply pass unchanged */
    else 
      rargv[rarg] = argv[arg];

    /* next argument */
    rarg++;
  } /* for all command line arguments */

  /* sanity checks */
  if (!resultfile[0] && !output_help && !output_version) {
      LogPrintf (LOG_CRITICAL, "ERROR: no result file has been specified\n");
      res = HIERARCHICALSEARCH_EFILE;
  }
  if (bundle_size != current_config_file) {
    LogPrintf (LOG_CRITICAL, "ERROR: bundle size %d doesn't match number of config files %d\n",
	       bundle_size, current_config_file);
    res = HIERARCHICALSEARCH_EFILE;
  }


#if DEBUG_COMMAND_LINE_MANGLING
  /* debug: dump the modified command line */
  {
    int i;
    fputs("command line:",stderr);
    for(i=0;i<rargc;i++)
      fprintf(stderr," %s",rargv[i]);
    fprintf(stderr,"\n");
  }
#endif


  LogPrintf (LOG_DEBUG, "Flags: "
#ifdef LAL_NDEBUG
	  "LAL_NDEBUG"
#else
	  "LAL_DEBUG"
#endif
#if __OPTIMIZE__
	  ", OPTIMIZE"
#endif
#ifdef HS_OPTIMIZATION
	  ", HS_OPTIMIZATION"
#endif
#ifdef GC_SSE2_OPT
	  ", GC_SSE2_OPT"
#endif
#ifdef USE_CUDA
	  ", CUDA"
#endif
#ifdef AUTOVECT_HOTLOOP
	  ", AUTOVECT"
#endif
#if _ARCH_PPC
	  ", PPC"
#elif __x86_64__
	  ", X64"
#elif __i386__
	  ", i386"
#endif
#if __ALTIVEC__
	  ", ALTIVEC"
#endif
#if __SSE__
	  ", SSE"
#endif
#if __SSE2__
	  ", SSE2"
#endif
#ifdef __BIG_ENDIAN__
	  ", BIG_ENDIAN"
#endif
#if __GNUC__
	  ", GNUC"
#endif
#if __X86__
          " X86"
#endif
#if __GNUX86__
          " GNUX86"
#endif
#ifdef _MSC_VER
	  ", _MSC_VER:%d\n",_MSC_VER
#else
	  "\n"
#endif
	  );

#ifdef __GLIBC__
  /* log the glibc version */
  LogPrintf (LOG_DEBUG, "glibc version/release: %s/%s\n", gnu_get_libc_version(), gnu_get_libc_release());
  /* test mytime() */
  mytime();
  fputs(" - mytime()\n",stderr);
#endif

  /* if there already was an error, there is no use in continuing */
  if (res) {
    LogPrintf (LOG_CRITICAL, "ERROR: error %d in command-line parsing\n", res);
    boinc_finish(boinc_finish_status=res);
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
    run_gdb(gdb_dump_core);
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

  if(!bundle_size && (fp = boinc_fopen(resultfile,"r"))) {
    fclose(fp);
    LogPrintf (LOG_NORMAL, "WARNING: Resultfile '%s' present - doing nothing\n", resultfile);
    resultfile_present = 1;
  }

#ifdef BOINC_APIV6
  if(setup_shmem())
    LogPrintf (LOG_NORMAL, "WARNING: Couldn't set up communication with graphics process\n");
  else
    LogPrintf (LOG_DEBUG, "Set up communication with graphics process.\n");
#endif

  /* if the program was called to output the version, output the BOINC revision, too */
  if(output_version)
#ifdef BUILD_INFO
    printf("%%%% " BUILD_INFO "\n");
#endif
    printf("%%%% BOINC: " SVN_VERSION "\n");

  if (output_help || output_version || !resultfile_present) {

    /* loop over WUs in bundles */
    current_config_file = 0;
    do {

      if (bundle_size) {
	unsigned int rlen = strlen(resultfile);
	strcpy(wu_result_file, resultfile);
	myltoa(second_outfile ? current_config_file*2 : current_config_file, &wu_result_file[rlen-1], MAX_PATH_LEN-rlen);
	*config_file_arg = config_files[current_config_file];
      }

      /* CALL WORKER's MAIN()
       */
      res = MAIN(rargc,rargv);
      if (res) {
	LogPrintf (LOG_CRITICAL, "ERROR: MAIN() returned with error '%d'\n",res);
      }

      /* if there is a file <wuname>_<instance>_0-LV, rename it to <wuname>_<instance>_1 */
      if (second_outfile) {
	unsigned int len = strlen(resultfile);
	char*lv_file = (char*)malloc(len+6);
	if (lv_file) {
	  strcpy(lv_file, resultfile);
	  strcat(lv_file, "-LV");
	  if (boinc_file_exists(lv_file) && resultfile[len-1]=='0') {
	    myltoa(current_config_file*2+1, &resultfile[len-1], 5);
	    boinc_rename(lv_file,resultfile);
	  }
	} else {
	  LogPrintf(LOG_CRITICAL,"ERROR: out of memory, can't allocate lv_file\n");
	  res = HIERARCHICALSEARCH_EMEM;
	}
      }

      current_config_file ++;
    } while (!res && bundle_size && current_config_file < bundle_size);

#ifdef __GNUX86__
    {
      fpuw_t fpstat = get_fpu_status();
      fputs("FPU status flags: ",stderr);
      PRINT_FPU_STATUS_FLAGS(fpstat);
      fputs("\n",stderr);
    }
#endif

    /* if the program was called for help, we write out usage for command-line options this wrapper adds to it and exit */
    if(output_help) {
      printf("Additional options the BOINC version understands:\n");
      printf("      --WUfpops         REAL     \"flops estimation\", passed to the BOINC client as the number of Flops\n");
#ifdef _MSC_VER
      printf("      --BreakPoint       -       fire up the Windows Runtime Debugger at internal breakpoint\n");
#elif defined(__GNUC__)
      printf("      --BreakPoint       -       attach gdb to dump a corefile, then continue\n");
#endif
      printf("      --CrashFPU         -       drain the FPU stack to test FPE\n");
      printf("      --TestNaN          -       raise a NaN to test FPE\n");
      printf("      --TestSQRT         -       try to calculate sqrt(-1) to test FPE\n");
      boinc_finish(boinc_finish_status=0);
    }

  }

  /* FIXME: HANDLE OUTPUT FILES
   */

#ifdef HAVE_BOINC_ZIP
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
#endif

  LogPrintf (LOG_NORMAL, "done. calling boinc_finish(%d).\n",res);
  boinc_finish(boinc_finish_status=res);
} /* worker() */



/**
 * the main function of the BOINC App
 * deals with boinc_init(_graphics) and calls the worker
 */

int main(int argc, char**argv) {
  FILE* fp_debug;
  char* name_debug;
  int skipsighandler = 0;

  /* init BOINC diagnostics */
  boinc_init_diagnostics(BOINC_DIAG_DUMPCALLSTACKENABLED |
                         BOINC_DIAG_HEAPCHECKENABLED |
                         BOINC_DIAG_ARCHIVESTDERR |
                         BOINC_DIAG_REDIRECTSTDERR |
                         BOINC_DIAG_TRACETOSTDERR);

  LogSetLevel(LOG_DETAIL); /* as long as we are debugging */

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
	}
      else
	{
	  LogPrintf (LOG_NORMAL, "...with no parsable int: Setting lalDebugLevel -> 1\n");
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
	snprintf(commandstring,sizeof(commandstring),"ddd %s %d &", resolved_name ,process_id);
	system(commandstring);
	sleep(20);
      }
    } /* DDD DEBUGGING */

  /* see if user has created a DEBUG_GDB_FNAME file: turn on debuggin using 'gdb' */
  /* record the debugging filename found in name_debug to use it as a command-file for gdb*/
  name_debug="../../" DEBUG_GDB_FNAME;
  if ((fp_debug=fopen(name_debug, "r"))) {
    fclose(fp_debug);
  } else {
    name_debug="." DEBUG_GDB_FNAME;
    if ((fp_debug=fopen(name_debug, "r"))) {
      fclose(fp_debug);
    } else {
      name_debug=NULL;
    }
  }
  if(name_debug)
    {
      char commandstring[256];
      char resolved_name[MAXFILENAMELENGTH];
      char *ptr;
      pid_t process_id=getpid();
      
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
	snprintf(commandstring,sizeof(commandstring),
		 "gdb -n -x %s %s %d >&2&",
		 name_debug, resolved_name, process_id);
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

#if HAVE_EXCHNDL
  ExchndlSetup();
#elif defined(_WIN32)
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

#ifdef DLOPEN_LIBGCC
  {
    void *lib_handle = dlopen("libgcc_s.so.1", RTLD_LAZY);
    if(lib_handle) {
      LogPrintf (LOG_DEBUG, "Successfully loaded libgcc_s.so.1\n");
      libgcc_s_loaded = 1;
    } else {
      LogPrintf (LOG_DEBUG, "Couldn't load libgcc_s.so.1: %s\n", dlerror());
      libgcc_s_loaded = -1;
    }
  }
#endif

#ifdef _NO_MSC_VER
  if (try_load_dlls(delayload_dlls, "ERROR: Failed to load %s - terminating\n")) {
    LogPrintf(LOG_NORMAL,"ERROR: Loading of mandantory DLLs failed\n");
    boinc_init();
    boinc_finish(boinc_finish_status=29);
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
  boinc_finish(boinc_finish_status=0);
  /* boinc_finish() ends the program, we never get here */
  return(0);
}

/* CHECKPOINTING FUNCTIONS */

/** log an I/O error, i.e. source code line no., ferror, errno and strerror, and doserrno on Windows, too */
#ifdef _WIN32
#define LOGIOERROR(mess,filename) \
    LogPrintf(LOG_CRITICAL, "ERROR: %s %s: line:%d, doserr:%d, ferr:%d, errno:%d: %s\n",\
	      mess,filename,__LINE__,_doserrno,ferror(fp),errno,strerror(errno))
#else
#define LOGIOERROR(mess,filename) \
    LogPrintf(LOG_CRITICAL, "ERROR: %s %s: line:%d, ferr:%d, errno:%d: %s\n",\
	      mess,filename,__LINE__,ferror(fp),errno,strerror(errno))
#endif

#ifdef HIERARCHSEARCHGCT

/**
 * sets a checkpoint.
 */
int write_boinc_gct_checkpoint(const char*filename, toplist_t*tl, toplist_t*t2, UINT4 counter, BOOLEAN do_sync) {
  int ret;
  /* make sure the exception mask isn't messed up by a badly written device driver etc.,
     so restore it periodically */
  enable_floating_point_exceptions();
  /* checkpoint every time (i.e. sky position) if FORCE_CHECKPOINTING */
#ifndef FORCE_CHECKPOINTING
  if (!(boinc_is_standalone() || boinc_time_to_checkpoint()))
    return 1; /* >0, no checkpoint written, no error */
#endif
  ret = write_gct_checkpoint(filename, tl, t2, counter, do_sync);
  fprintf(stderr,"c\n");
  boinc_checkpoint_completed();
  return(ret);
}

#else /* #ifdef HIERARCHSEARCHGCT */

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
    UINT4 s = total;
    s = strlen(outputname)+1;
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
	  if((unsigned int)len == fread(buf,sizeof(char),len,fp))
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

void set_checkpoint (void) {
  /* make sure the exception mask isn't messed up by a badly written device driver etc.,
     so restore it periodically */
  enable_floating_point_exceptions();
  /* checkpoint every time (i.e. sky position) if FORCE_CHECKPOINTING */
#ifndef FORCE_CHECKPOINTING
  if (boinc_time_to_checkpoint())
#endif
    {
      write_hfs_checkpoint(cptfilename, toplist, last_count, do_sync);
      fprintf(stderr,"c\n");
      boinc_checkpoint_completed();
    }
}

/**
 * finally writes a minimal (compacted) version of the toplist and cleans up
 * all structures related to the toplist. After that, the toplist is invalid.
 */
void write_and_close_checkpointed_file (void) {
  write_hfs_oputput(outfilename,toplist);
}

#endif /* #ifdef HIERARCHSEARCHGCT */


/* Experimental and / or debugging stuff */

/** attach gdb to the running process and do something; for debugging. */
void run_gdb(gdb_cmd command) {
#ifdef __GLIBC__
  fputs("attaching gdb...\n",stderr);
  char cmd[256];
  pid_t pid=getpid();
  switch (command) {
  case gdb_attach:
    /* FIXME: that should write a stackdump when somthing (signal) happens, not right away */
    snprintf(cmd, sizeof(cmd), "echo bt | gdb -batch %s %d -ex cont >&2", global_argv[0], pid);
    break;
  case gdb_dump_core:
    snprintf(cmd, sizeof(cmd), "echo y | gdb -batch %s %d -ex gcore -ex quit >&2", global_argv[0], pid);
    break;
  }
  fprintf(stderr,"executing '%s'\n",cmd);
  system(cmd);
  sleep(20);
#endif
}


/**
 * sets the FPU control word.
 * The argument should be a (possibly modified)
 * fpuw_t gotten from get_fpu_control_word()
 */
void set_fpu_control_word(const fpuw_t cword) {
  static fpuw_t fpucw;
  fpucw = cword;
#ifdef _MSC_VER
  __asm fldcw fpucw;
#elif defined(__GNUX86__)
  __asm("fldcw %0\n\t" : : "m" (fpucw));
#endif
}

/** returns the fpu control word */
fpuw_t get_fpu_control_word(void) {
  static fpuw_t fpucw = 0;
#ifdef _MSC_VER
  __asm fstcw fpucw;
#elif defined(__GNUX86__)
  __asm("fstcw %0\n\t" : "=m" (fpucw));
#endif
  return(fpucw);
}

/** returns the fpu status word */
fpuw_t get_fpu_status(void) {
  static fpuw_t fpusw = 0;
#ifdef _MSC_VER
  __asm fnstsw fpusw;
#elif defined(__GNUX86__)
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
#elif defined(__GNUX86__)
  __asm("ldmxcsr %0\n\t" : : "m" (ssecw));
#endif
}

/** returns the sse control/status word */
ssew_t get_sse_control_status(void) {
  static ssew_t ssesw = 0;
#ifdef _MSC_VER
  __asm stmxcsr ssesw;
#elif defined(__GNUX86__)
  __asm("stmxcsr %0\n\t" : "=m" (ssesw));
#endif
  return(ssesw);
}

static void drain_fpu_stack(void) {
  static double dummy;
#ifdef __GNUX86__
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
#if defined(_MSC_VER) || defined(__GNUX86__)
  /* write out the masked FPU exceptions */
  /*
  {
    fpuw_t fpstat = get_fpu_status();
    fputs("FPU status flags: ",stderr);
    PRINT_FPU_STATUS_FLAGS(fpstat);
    fputs("\n",stderr);
  }
  */

  {
    fpuw_t fpstat;
    
    fpstat = get_fpu_control_word();
    /*
    fprintf(stderr,"FPU masked exceptions now: %4x:",fpstat);
    PRINT_FPU_EXCEPTION_MASK(fpstat);
    fputs("\n",stderr);
    */

    /* throw an exception at an invalid operation */
    fpstat &= ~FPU_STATUS_INVALID;
    set_fpu_control_word(fpstat);    

    /*
    fprintf(stderr,"FPU masked exceptions set: %4x:",fpstat);
    PRINT_FPU_EXCEPTION_MASK(fpstat);
    fputs("\n",stderr);
    */

    /* this is weird - somtimes gcc seems to cache the fpstat value
       (e.g. on MacOS Intel), so the second reading of the control_word
       doesn't seem to do what is expected
    fpstat = 0;

    fpstat = get_fpu_control_word();
    fprintf(stderr,"FPU exception mask set to:  %4x:",fpstat);
    PRINT_FPU_EXCEPTION_MASK(fpstat);
    fputs("\n",stderr);
    */
#if __SSE__
    set_sse_control_status(get_sse_control_status() & ~SSE_MASK_INVALID);
#endif
  }
#endif
}

int segfault (void) {
  volatile int i = *((int*)1);
  return(i);
}
