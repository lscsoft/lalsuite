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
is specified with an absolute path (beginning with a \verb+/+), or a specific
path (beginning with a \verb+./+ or a \verb+../+), the directory
that the file is in is obtained from the environment variable
\verb+LAL_DATA_PATH+, which must be set at run-time.  (If the environment
variable is not set, the default path is \verb+.+ --- i.e., the current
directory.)

\verb+LAL_DATA_PATH+ should typically be set to
\verb+/usr/local/share/lal+, or wherever LAL data is installed in your system
(which may be different if you used a \verb+--prefix+ argument when
configuring LAL), but when the test suite is run with \verb+make check+, the
variable \verb+LAL_DATA_PATH+ is set to the current source directory.  If the
filename (including the directory path) is too long (more than 256
characters), \verb+LALOpenDataFile()+ returns \verb+NULL+ and sets
\verb+errno+ to \verb+ENAMETOOLONG+.

\verb+LAL_DATA_PATH+ can be any colon-delimeted list of directories, which
are searched in order (just like the \verb+PATH+ environment variable).
An extra colon inserts the default data directory
($\langle$prefix$\rangle$\verb+/share/lal+) into the search path at that
point.  E.g., a leading/trailing colon will look for the default data
directory at the start/end of the list of directories.

It is strongly recommended that
\verb+LALOpenDataFile()+ be used when writing test code.

\vfill{\footnotesize\input{FileIOCV}}

</lalLaTeX>
#endif /* autodoc block */

#include <stdlib.h>
#include <errno.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/FileIO.h>

#define STR( x ) #x
#define XSTR( x ) STR( x )
#define INFOMSG( msg, file ) ( ( lalDebugLevel & LALINFO ) ? \
    LALPrintError( "Info: function LALOpenDataFile, file " __FILE__ ", line " \
      XSTR( __LINE__ ) ", $Id$\n\t%s %s\n", msg, file ) : 0 )
#define ERRORMSG( file ) ( ( lalDebugLevel & LALERROR ) ? \
    LALPrintError( "Info: function LALOpenDataFile, file " __FILE__ ", line " \
      XSTR( __LINE__ ) ", $Id$\n\tCould not open data file %s\n", file ) : 0 )

#ifndef LAL_PREFIX
#define LAL_PREFIX "/usr/local"
#endif


/* <lalVerbatim file="FileIOCP"> */
FILE *
LALOpenDataFile( const char *fname )
{ /* </lalVerbatim> */
  FILE *fp;
  const char *path;
  char *datapath;	/* locally allocated copy of env-var LAL_DATA_PATH */
  const char *p0; 	/* pointer to current sub-path of datapath*/
  char *p1;		/* pointer to next sub-path */
  char  fdata[265];
  int   n;

  if ( ! fname )
    return NULL;

  if ( *fname == '/' ) /* absolute path is given */
  {
    fp = LALFopen( fname, "r" );
    if ( ! fp )
      ERRORMSG( fname );
    else
      INFOMSG( "Opening data file", fname );
    return fp;
  }

  n = strlen( fname );
  if ( *fname == '.' && n > 0 && ( fname[1] == '/' || ( n > 1 && fname[1] == '.'
          && fname[2] == '/' ) ) ) /* specific path is given */
  {
    fp = LALFopen( fname, "r" );
    if ( ! fp )
      ERRORMSG( fname );
    else
      INFOMSG( "Opening data file", fname );
    return fp;
  }

  path = getenv( "LAL_DATA_PATH" );

  if ( ! path || ! strlen( path ) ) /* path is NULL or empty */
  {
    fp = LALFopen( fname, "r" );
    if ( ! fp )
      ERRORMSG( fname );
    else
      INFOMSG( "Opening data file", fname );
    return fp;
  }

  /* scan through all directories in colon-delmited list of directories */
  if ( (datapath = LALCalloc (strlen(path)+1, 1)) == NULL)	/* we need local copy */
    {
      ERRORMSG( fname );
      return NULL;
    }
  strcpy (datapath, path);
  p0 = datapath;
  do {
    p1 = strchr( p0, ':' ); /* look for additional directories */
    if ( p1 ) /* there are more things in the list */
      *p1++ = 0; /* NUL-terminate current directory */
    if ( ! strlen( p0 ) ) /* this directory is empty */
      p0 = LAL_PREFIX "/share/lal"; /* default data directory */

    n = LALSnprintf( fdata, sizeof(fdata), "%s/%s", p0 ? p0 : ".", fname );
    if ( n > (int) sizeof( fdata ) ) /* data file name too long */
    {
      errno = ENAMETOOLONG;
      LALFree (datapath);
      return NULL;
    }

    INFOMSG( "Looking for file", fdata );
    fp = LALFopen( fdata, "r" );
    if ( fp ) /* we've found it! */
    {
      INFOMSG( "Opening data file", fdata );
      LALFree (datapath);
      return fp;
    }

    p0 = p1;
  }
  while ( p0 );

  LALFree (datapath);
  ERRORMSG( fname );
  return NULL;
}
