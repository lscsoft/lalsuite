/************************************ <lalVerbatim file="PrintVectorCV">
Author: Allen, B.
$Id$
************************************* </lalVerbatim> */

/* <lalLaTeX>

\subsection{Module \texttt{PrintVector.c}}
\label{ss:PrintVector.c}

Print a REAL4Vector object into a file.  For use in non-production and test
code only.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{PrintVectorCP}
\index{\verb&PrintVector()&}

\subsubsection*{Description}

This function prints the elements of \verb+vector+ into a file.  Note: the
file names are \verb+PrintVector.000+, \verb+PrintVector.001+, and so on.  The
file numbers are incremented with each additional call.  This function is for
debugging use only: it uses a static internal variable to keep track of the
file number so it should not be used in any real analysis codes.

\subsubsection*{Algorithm}

\subsubsection*{Uses}

\subsubsection*{Notes}

This function uses an internal static variable to keep track of file numbers.
For this reason it should only be used for debugging purposes in test
functions, not in any production code.

\vfill{\footnotesize\input{PrintVectorCV}}

</lalLaTeX> */


#include "LALStdlib.h"
#include <stdio.h>
#include "LALDatatypes.h"
#include "PrintVector.h"

NRCSID( PRINTVECTORC, "$Id$" );

/* <lalVerbatim file="PrintVectorCP"> */
void
PrintVector( REAL4Vector *vector )
{ /* </lalVerbatim> */
  int i;
  static int fileno=0;
  FILE *fp;
  char fname[256];


  if (vector==NULL) return;

  /* open output file */
  sprintf(fname,"PrintVector.%03d",fileno++);
  fp=fopen(fname,"w");

  for (i=0;i<vector->length;i++)
    fprintf(fp,"%i %f\n",i,vector->data[i]);

  fclose(fp);

  return;
}
