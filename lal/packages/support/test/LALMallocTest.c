/*************** <lalVerbatim file="LALMallocTestCV"> *************
$Id$
**************** </lalVerbatim> ***********************************/

/* <lalLaTeX>

\subsection{Program \texttt{LALMallocTest.c}}
\label{s:LALMallocTest.c}

Tests the routines in \verb@LALMalloc.h@.

\subsubsection*{Usage}
\begin{verbatim}
LALMallocTest
\end{verbatim}

\subsubsection*{Description}

This program has ugly code for testing the LAL memory allocation and freeing
routines.  It always returns success, so it is not very useful.  The actual
success should be interpreted from the output messages.  This should be
changed someday.

\subsubsection*{Exit codes}
\begin{tabular}{|c|l|}
\hline
 Code & Explanation                   \\
\hline
\tt 0 & Always returned.              \\
\hline
\end{tabular}

\subsubsection*{Uses}
\begin{verbatim}
lalDebugLevel
LALMalloc()
LALFree()
LALCheckMemoryLeaks()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALMallocTestCV}}

</lalLaTeX> */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <lal/LALStdlib.h>

int lalDebugLevel = 3;

int raise (int sig);

int main( void )
{
  int    *p;
  int    *q;
  int     i;
  size_t  n;

  fprintf (stderr, "Test 0: allocate and free an int\n");
  p = (int *) LALMalloc (sizeof(int));
  LALFree (p);
  fprintf (stderr, "\n");

  fprintf (stderr, "Test 1: zero size alloc\n");
  LALMalloc (0);
  fprintf (stderr, "\n");

  fprintf (stderr, "Test 2: NULL ptr free\n");
  LALFree (NULL);
  fprintf (stderr, "\n");

  fprintf (stderr, "Test 3: corrupt magic\n");
  p = (int *) LALMalloc (2*sizeof(int));
  n = ((size_t *)p)[-1];
  ((size_t *)p)[-1] = 17;
  LALFree (p);
  ((size_t *)p)[-1] = n;
  LALFree (p);
  fprintf (stderr, "\n");

  fprintf (stderr, "Test 4: corrupt size\n");
  p = (int *) LALMalloc (4*sizeof(int));
  n = ((size_t *)p)[-2];
  ((size_t *)p)[-2] = -1;
  LALFree (p);
  ((size_t *)p)[-2] = n;
  LALFree (p);
  fprintf (stderr, "\n");

  fprintf (stderr, "Test 5: overwrite bounds\n");
  p = (int *) LALMalloc (8*sizeof(int));
  i = p[10];
  p[10] = 17;
  LALFree (p);
  p[10] = i;
  LALFree (p);
  fprintf (stderr, "\n");

  fprintf (stderr, "Test 6: free freed memory\n");
  q = p = (int *) LALMalloc (16*sizeof(int));
  LALFree (p);
  LALFree (q);
  fprintf (stderr, "\n");

  fprintf (stderr, "Test 7: should be no memory leaks so far\n");
  LALCheckMemoryLeaks ();
  fprintf (stderr, "\n");

  fprintf (stderr, "Test 8: forget to free some memory\n");
  p = (int *) LALMalloc (32*sizeof(int));
  LALCheckMemoryLeaks ();
  LALFree (p);
  LALCheckMemoryLeaks ();
  fprintf (stderr, "\n");

  fprintf (stderr, "Test 9: free the `wrong' amount of memory\n");
  p = (int *) LALMalloc (64*sizeof(int));
  q = (int *) malloc (2*sizeof(size_t) + 2*64*sizeof(int));
  memcpy ( q, (size_t *)p - 2, 2*sizeof(size_t) + 2*64*sizeof(int));
  ((size_t *)p)[-2] = 32*sizeof(int);
  LALFree (p);
  LALCheckMemoryLeaks ();
  p = (int *) ((size_t *) malloc (2*sizeof(size_t) + 2*64*sizeof(int)) + 2);
  memcpy ((size_t *)p - 2, q, 2*sizeof(size_t) + 2*64*sizeof(int));
  free (q);
  q = LALMalloc (32*sizeof(int));
  free ((size_t *)q - 2);
  LALFree (p);
  LALCheckMemoryLeaks ();
  fprintf (stderr, "\n");

  fprintf (stderr, "Test 10: try calloc\n");
  p = (int *) LALCalloc (2, 16*sizeof(int));
  LALFree (p);
  LALCheckMemoryLeaks ();
  fprintf (stderr, "\n");

  fprintf (stderr, "Test 11: try several reallocs\n");
  p = (int *) LALRealloc (NULL, 4*sizeof(int));
  p = (int *) LALRealloc (p, 4*sizeof(int));
  p = (int *) LALRealloc (p, 2*sizeof(int));
  p = (int *) LALRealloc (p, 8*sizeof(int));
  p = (int *) LALRealloc (p, 0);
  LALCheckMemoryLeaks ();
  fprintf (stderr, "\n");

  fprintf (stderr, "Test 12: try to free too much memory\n");
  p = (int *) LALMalloc (128*sizeof(int));
  free ( (size_t *)p - 2 );
  p = (int *) ((size_t *) malloc (2*sizeof(size_t) + 2*256*sizeof(int)) + 2 );
  q = (int *) LALMalloc (256*sizeof(int));
  memcpy ((size_t *)p - 2, (size_t *)q - 2, 2*sizeof(size_t) + 2*256*sizeof(int));
  LALFree (q);
  ((size_t *)p)[-2] = 256*sizeof(int);
  LALFree (p);

  return 0;
}

int raise (int sig)
{
  fprintf (stderr, "Ignored signal %d\n", sig);
  return 1;
}
