#if 0 /* autodoc block */
<lalVerbatim file="FileIOCV">
$Id$
</lalVerbatim>

<lalLaTeX>
\subsection{Module \texttt{FileIO.c}}

File IO routines for LAL.  These should \emph{not} be used in routines within
the LAL library---only within test programs.

\subsection*{Prototypes}
\input{FileIOCP}
\begin{verbatim}
FILE *LALFopen( const char *path, const char *mode );
int LALFclose( FILE *stream );
\end{verbatim}
\idx{LALOpenDataFile()}
\idx{LALFopen()}
\idx{LALFclose()}

\subsection*{Description}

The routines \verb+LALFopen()+ and \verb+LALFclose()+ are macro defined to be
the same as the standard C routines \verb+fopen()+ and \verb+fclose()+.  These
should only be used in test programs.

The routine \verb+LALOpenDataFile()+ is used to open a data file for reading.
This routine is also to be used in test programs only.  Unless the data file
is specified with an absolute path (beginning with a \verb+/+), the directory
that the file is in is obtained from the environment variable
\verb+LAL_DATA_PATH+, which must be set at run-time.  (If the environment
variable is not set, the default path is \verb+.+ --- i.e., the current
directory.)  \verb+LAL_DATA_PATH+ should typically be set to
\verb+/usr/local/share/lal+, or wherever LAL data is installed in your system
(which may be different if you used a \verb+--prefix+ argument when
configuring LAL), but when the test suite is run with \verb+make check+, the
variable \verb+LAL_DATA_PATH+ is set to the current source directory.  If the
filename (including the directory path) is too long (more than 256
characters), \verb+LALOpenDataFile()+ returns \verb+NULL+ and sets
\verb+errno+ to \verb+ENAMETOOLONG+.  It is strongly recommended that
\verb+LALOpenDataFile()+ be used when writing test code.

\vfill{\footnotesize\input{FileIOCV}}

</lalLaTeX>
#endif /* autodoc block */

#include <stdlib.h>
#include <errno.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/FileIO.h>

/*
 * This needs to be in some source code somewhere,
 * so may as well put it here.
 */
#if defined(NDEBUG) || defined(LAL_NDEBUG)
const int lalNoDebug = 1;
#else
const int lalNoDebug = 0;
#endif

/* <lalVerbatim file="FileIOCP"> */
FILE *
LALOpenDataFile( const char *fname )
{ /* </lalVerbatim> */
  char *path;
  char  fdata[265];
  int   n;

  if ( ! fname )
    return NULL;

  if ( *fname == '/' ) /* absolute path is given */
    return LALFopen( fname, "r" );

  path = getenv( "LAL_DATA_PATH" );

  n = LALSnprintf( fdata, sizeof( fdata ), "%s/%s", path ? path : ".", fname );
  if ( n > (int) sizeof( fdata ) ) /* data file name too long */
  {
    errno = ENAMETOOLONG;
    return NULL;
  }

  return LALFopen( fdata, "r" );
}
