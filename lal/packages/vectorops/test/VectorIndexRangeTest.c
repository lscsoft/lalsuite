/*
*  Copyright (C) 2007 David Chin
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

static void
trail_status_maybe(LALStatus *status);

int
main(int argc, char **argv)
{
  static LALStatus status;
  const UINT4      vec_length = 64;
  UINT4            iter;
  REAL4Vector     *x = NULL;
  REAL4Vector     *y = NULL;
  REAL4Vector     *z = NULL;
  REAL4VectorPair  vector_pair;
  UINT4Vector     *index_range = NULL;
  BOOLEAN          result_ok_p = TRUE;
  int              retval = 1;

  if (argc > 1)
    {
      lalDebugLevel = atoi(argv[1]);
      verbose_p = TRUE;
    }

  trail_status_maybe(&status);

  LALU4CreateVector(&status, &index_range, 2);
  test_status(&status, CODES(0), 1);
  trail_status_maybe(&status);

  index_range->data[0] = 32;
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

  if (lalDebugLevel & 8)
    printf("* * * * * * * * * * * *\n");

  trail_status_maybe(&status);

  /*
  LALSCreateVector(&status, &y, 1);
  test_status(&status, CODES(0), 1);
  printf("y->length = %d\n", y->length);
  */

  LALSVectorIndexRange(&status, &y, x, index_range);
  test_status(&status, CODES(0), 1);
  trail_status_maybe(&status);

  if (verbose_p)
    {
      printf(" - - - - - - - -\n");
      printf("y = %#x\n", (unsigned int)y);
      printf("y->length = %d\n", y->length);
    }

  for (iter = 0; iter < y->length; ++iter)
    {
      if (lalDebugLevel & 8)
        printf("y->data[%.2d] = % 7.13e; x->data[%.2d] = % 7.13e\n",
               iter, y->data[iter], iter, x->data[iter]);
      result_ok_p = result_ok_p &&
        (y->data[iter] == x->data[iter + index_range->data[0]]);
    }

  trail_status_maybe(&status);

  if (!result_ok_p)
    {
      if (verbose_p)
        fprintf(stderr, "%s: LAL*VectorIndexRange() failed: line %d\n",
                argv[0], __LINE__);
      goto end;
    }

  if (verbose_p)
    {
      printf("*-*-*-*-*-*-*\n");
      printf("trivial range:\n");
    }

  /* test trivial range */
  index_range->data[0] = 32;
  index_range->data[1] = 32;
  LALSVectorIndexRange(&status, &y, x, index_range);
  for (iter = 0; iter < y->length; ++iter)
    {
      if (lalDebugLevel & 8)
        printf("y->data[%.2d] = % 7.13e; x->data[%.2d] = % 7.13e\n",
               iter, y->data[iter], iter, x->data[iter]);
      result_ok_p = result_ok_p &&
        (y->data[iter] == x->data[iter + index_range->data[0]]);
    }

  if (!result_ok_p)
    {
      if (verbose_p)
        fprintf(stderr, "%s: LAL*VectorIndexRange() failed: line %d\n",
                argv[0], __LINE__);
      goto end;
    }

  if (verbose_p)
    {
      printf("\n");
      printf("*********************\n");
    }

  /*
   * VectorIndexHole
   */
  LALSDestroyVector(&status, &y);
  test_status(&status, CODES(0), 1);

  vector_pair.head = &y;
  vector_pair.tail = &z;

  LALSVectorIndexHole(&status, &vector_pair, x, index_range);

  if (lalDebugLevel & 8)
    {
      printf("y->length = %d\n", y->length);
      printf("z->length = %d\n", z->length);
      printf("\n");
    }


  for (iter = 0; iter < y->length; ++iter)
    {
      if (lalDebugLevel & 8)
        printf("y->data[%.2d] = % 7.13e; x->data[%.2d] = % 7.13e\n",
               iter, y->data[iter], iter, x->data[iter]);
      result_ok_p = result_ok_p && (y->data[iter] == x->data[iter]);
    }

  if (!result_ok_p)
    {
      if (verbose_p)
        fprintf(stderr, "%s: LAL*VectorIndexHole() failed: line %d\n",
                argv[0], __LINE__);
      goto end;
    }

  for (iter = 0; iter < z->length; ++iter)
    {
      if (lalDebugLevel & 8)
        printf("z->data[%.2d] = % 7.13e; x->data[%.2d] = % 7.13e\n",
               iter, z->data[iter], iter, x->data[iter]);
      result_ok_p = result_ok_p &&
        (z->data[iter] == x->data[iter + index_range->data[1] + 1]);
    }

  if (!result_ok_p)
    {
      if (verbose_p)
        fprintf(stderr, "%s: LAL*VectorIndexHole() failed: line %d\n",
                argv[0], __LINE__);
      goto end;
    }

  if (lalDebugLevel & 8)
    {
      printf("y->length = %d\n", y->length);
      printf("z->length = %d\n", z->length);
      printf("\n");
    }

  /*
   * Housekeeping
   */
 end:
  LALSDestroyVector(&status, &x);
  test_status(&status, CODES(0), 1);

  LALSDestroyVector(&status, &y);
  test_status(&status, CODES(0), 1);

  LALSDestroyVector(&status, &z);
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
trail_status_maybe(LALStatus *status)
{
  if (lalDebugLevel & 15)
    {
      printf("TRAIL STATUS:\n");
      printf("  status = %#x\n", (unsigned int)status);
      if (status)
        printf("  status->statusPtr = %#x\n",
               (unsigned int)(status->statusPtr));
      printf("\n");

      if (status)
        {
          REPORTSTATUS(status);
          trail_status_maybe(status->statusPtr);
        }
    }
}


static void
usage(char **argv)
{
  fprintf(stderr, "%s: read the code\n", argv[0]);
  exit(0);
}
