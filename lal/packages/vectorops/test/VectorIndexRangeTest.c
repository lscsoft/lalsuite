/*************** <lalVerbatim file="VectorIndexRangeTestCV"> **************
$Id$
**************** </lalVerbatim> ***********************************/

/* <lalLaTeX>
   
\subsection{Program \texttt{VectorIndexRangeTest.c}}
\label{ss:VectorIndexRangeTest.c}

Tests the routines in \verb@VectorIndexRange.h@.  Exercises some of the error
conditions and makes sure that they work.

\subsubsection*{Usage}
\begin{verbatim}
VectorIndexRangeTest [options]
Options:
  -h         print help
  -q         quiet: run silently
  -v         verbose: print extra information
  -d level   set lalDebugLevel to level
\end{verbatim}

\subsubsection*{Description}
\subsubsection*{Exit codes}
\begin{tabular}{|c|l|}
\hline
 Code & Explanation                   \\
\hline
\tt 0 & Success, normal exit.         \\
\tt 1 & Subroutine failed.            \\
\hline
\end{tabular}

\vfill{\footnotesize\input{VectorIndexRangeTestCV}}

</lalLaTeX> */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <config.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/VectorIndexRange.h>

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#ifdef HAVE_GETOPT_H
#include <getopt.h>
#endif


NRCSID (MAIN, "$Id$");

#define CODES_(x) #x
#define CODES(X)  CODES_(x)

#define FALSE 0
#define TRUE  1

extern char *optarg;
extern int   optind;

int lalDebugLevel = 0;
BOOLEAN verbose_p = FALSE;

static void
usage(char **argv);


static void
test_status(LALStatus *status, const char *expected_codes, int exit_code);

int
main(int argc, char **argv)
{
  static LALStatus status;
  const UINT4      vec_length = 64;
  UINT4            iter;
  REAL4Vector     *x = NULL;
  REAL4Vector     *y = NULL;
  UINT4Vector     *index_range = NULL;
  BOOLEAN          result_ok_p = TRUE;
  int              retval = 1;

  LALU4CreateVector(&status, &index_range, 2);
  test_status(&status, CODES(0), 1);

  index_range->data[0] = 16;
  index_range->data[1] = 48;
  if (lalDebugLevel & 8)
    {
      printf("index_range->length = %d\n", index_range->length);
      for (iter = 0; iter < 2; ++iter)
        {
          printf("index_range->data[%.2d] = %.3d\n", iter,
                 index_range->data[iter]);
        }
    }

  LALSCreateVector(&status, &x, vec_length);
  test_status(&status, CODES(0), 1);

  for (iter = 0; iter < x->length; ++iter)
    {
      x->data[iter] = (REAL4)iter;
      if (lalDebugLevel & 8)
        printf("x->data[%.2d] = % 7.13e\n", iter, x->data[iter]);
    }

  if (verbose_p)
    printf("* * * * * * * * * * * *\n");
  /*
  LALSCreateVector(&status, &y, 1); 
  test_status(&status, CODES(0), 1); 
  printf("y->length = %d\n", y->length);
  */
  
  LALSVectorIndexRange(&status, &y, x, index_range);
  test_status(&status, CODES(0), 1);

  if (verbose_p)
    {
      printf(" - - - - - - - -\n");
      printf("y = %#x\n", y);
      printf("y->length = %d\n", y->length);
    }

  for (iter = 0; iter < y->length; ++iter)
    {
      if (lalDebugLevel & 8)
        printf("y->data[%.2d] = % 7.13e\n", iter, y->data[iter]);
      result_ok_p = result_ok_p &&
        (y->data[iter] == x->data[iter + index_range->data[0]]);
    }

  LALSDestroyVector(&status, &x);
  test_status(&status, CODES(0), 1);

  LALSDestroyVector(&status, &y);
  test_status(&status, CODES(0), 1);

  LALU4DestroyVector(&status, &index_range);
  test_status(&status, CODES(0), 1);

  LALCheckMemoryLeaks();
  
  if (result_ok_p)
    retval = 0;
  
  return retval;
} /* END: main() */

/*
 * test_status() -- copied from VectorOpsTest.c -- TestStatus()
 *
 * Routine to check that the status code status->statusCode agrees with one of
 * the codes specified in the space-delimited string expected_codes; if not,
 * exit to the system with code exit_code.
 *
 */
static void
test_status( LALStatus *status, const char *expected_codes, int exit_code )
{
  char  str[64];
  char *tok;
  BOOLEAN status_ok_p = FALSE;

  if ( verbose_p )
  {
    REPORTSTATUS( status );
  }

  if ( strncpy( str, expected_codes, sizeof( str ) ) )
  {
    if ( ( tok = strtok( str, " " ) ) )
    {
      do
      {
        if ( status->statusCode == atoi( tok ) )
        {
          status_ok_p = TRUE;
        }
      }
      while ( ( tok = strtok( NULL, " " ) ) );
    }
    else
    {
      if ( status->statusCode == atoi( tok ) )
      {
        status_ok_p = TRUE;
      }
    }
  }

  if (!status_ok_p)
    {
      fprintf( stderr, "\nExiting to system with code %d\n", exit_code );
      exit( exit_code );
    }
  else
    {
      return;
    }
}


static void
usage(char **argv)
{
  fprintf(stderr, "%s: read the code\n", argv[0]);
  exit(0);
}
