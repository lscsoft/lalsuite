#if 0 /* autodoc block */
<lalVerbatim file="LALVersionTestCV">
$Id$
</lalVerbatim>

<lalLaTeX>
\subsection{Program \texttt{LALVersionTest.c}}
\label{s:LALVersionTest.c}

Prints the version and configure options of the LAL library being used.

\subsubsection*{Usage}
\begin{verbatim}
LALVersionTest
\end{verbatim}

\subsubsection*{Description}

This program prints the current version of LAL.\@  If the version information
in the library differs from the version information in the header file, this
program prints the two versions and exits with code 1.  This is useful for
determining which version of the LAL library and header files you are linking
to.

\subsubsection*{Exit codes}
\begin{tabular}{|c|l|}
\hline
 Code & Explanation                   \\
\hline
\tt 0 & Success, normal exit.         \\
\tt 1 & Version info in library disagrees with header file. \\
\tt 2 & Subroutine failed.            \\
\hline
\end{tabular}


\vfill{\footnotesize\input{LALVersionTestCV}}


</lalLaTeX>
#endif /* autodoc block */

#include <stdio.h>
#include <string.h>
#include <lal/LALStdlib.h>
#include <lal/LALVersion.h>

int lalDebugLevel = 0;

int main( void )
{
  static LALStatus status;
  char msg[1024];
  int verbose = 1;

  if ( strcmp( LAL_VERSION, lalVersion ) ||
       strcmp( LAL_CONFIGURE_ARGS, lalConfigureArgs ) ||
       strcmp( LAL_CONFIGURE_DATE, lalConfigureDate ) )
  {
    fputs( "LAL Version Mismatch!\n\n", stderr );
    fputs( "Header Version ",           stderr );
    fputs( LAL_VERSION,                 stderr );
    fputs( "\nCompiled on ",            stderr );
    fputs( LAL_CONFIGURE_DATE,          stderr );
    fputs( "\nWith arguments ",         stderr );
    fputs( LAL_CONFIGURE_ARGS,          stderr );
    fputs( "\n\n",                      stderr );
    fputs( "Library Version ",          stderr );
    fputs( lalVersion,                  stderr );
    fputs( "\nCompiled on ",            stderr );
    fputs( lalConfigureDate,            stderr );
    fputs( "\nWith arguments ",         stderr );
    fputs( lalConfigureArgs,            stderr );
    fputs( "\n",                        stderr );
    return 1;
  }

  LALVersion( &status, msg, sizeof( msg ), verbose );

  if ( status.statusCode )
  {
    LALStatus *next = &status;
    do
    {
      fputs( next->statusDescription, stderr );
      fputs( "\n", stderr );
      next = next->statusPtr;
    }
    while ( next );
    return 2;
  }

  puts( msg );

  return 0;
}
