#if 0 /* autodoc block */
<lalVerbatim file="LALVersionCV">
$Id$
</lalVerbatim>

<lalLaTeX>
\subsection*{Module \texttt{LALVersion.c}}

Routine that returns the version of LAL.

\subsubsection*{Prototypes}
\input{LALVersionCP}
\idx{LALVersion()}

\subsubsection*{Description}

This function writes a version message into the string buffer of specified
size (and is truncated if the buffer is too small).  Configuration information
is also provided if the verbose flag is set.

\vfill{\footnotesize\input{LALVersionCV}}

</lalLaTeX>
#endif /* autodoc block */

/*
 *
 * Use LAL's config.h file rather than LALConfig.h which may be from some
 * previous installation.
 *
 */
#include <config.h>

#include <stdio.h>
#include <lal/LALStatusMacros.h>
#include <lal/LALStdio.h>
#include <lal/LALVersion.h>

NRCSID( LALVERSIONC, "$Id$" );

const char *lalVersion       = LAL_VERSION;
const int   lalVersionMajor  = LAL_VERSION_MAJOR;
const int   lalVersionMinor  = LAL_VERSION_MINOR;
const int   lalVersionMicro  = LAL_VERSION_MICRO;
const char *lalConfigureArgs = LAL_CONFIGURE_ARGS;
const char *lalConfigureDate = LAL_CONFIGURE_DATE;
const char *lalCVSTag        = LAL_CVS_TAG;

/* <lalVerbatim file="LALVersionCP"> */
void
LALVersion( LALStatus *status, CHAR *message, UINT4 size, INT4 verbose )
{ /* </lalVerbatim> */
  INT4 nchar;
  INITSTATUS( status, "LALVersion", LALVERSIONC );

  ASSERT( message,  status, LALVERSIONH_ENULL, LALVERSIONH_MSGENULL );
  ASSERT( size > 0, status, LALVERSIONH_ESIZE, LALVERSIONH_MSGESIZE );

  nchar = verbose ?
    LALSnprintf( message, size,
        "LAL Version:         %s\n"
        "CVS Tag:             %s\n"
        "Build Date:          %s\n"
        "Configure Date:      %s\n"
        "Configure Arguments: %s\n"
        "(RCS %s)\n",
        lalVersion, LAL_CVS_TAG, lalBuildDate, lalConfigureDate,
        lalConfigureArgs, LALVERSIONC ) :
    LALSnprintf( message, size, "LAL Version: %s\n", lalVersion ) ;

  if ( nchar < 0 )
  {
    ABORT( status, LALVERSIONH_ESPRN, LALVERSIONH_MSGESPRN );
  }
  if ( nchar > (INT4) size )
  {
    ABORT( status, LALVERSIONH_ESHRT, LALVERSIONH_MSGESHRT );
  }

  RETURN( status );
}
