#if 0 /* autodoc block */
<lalVerbatim file="LALStdCV">
$Id$
</lalVerbatim>

<lalLaTeX>
\subsection{Module \texttt{LALStd.c}}

LAL replacement routines for standard C functions.  At present, this
just includes replacements for \texttt{snprintf}.

\subsection*{Prototypes}
\begin{verbatim}
int LALSnprintf( char *str, size_t size, const char *fmt, ... );
int LALVsnprintf( char *str, size_t size, const char *fmt, va_list ap );
\end{verbatim}
\idx{LALSnprintf()}
\idx{LALVsnprintf()}

\subsection*{Description}

The routines \verb+LALSnprintf()+ and \verb+LALVsnprintf()+ are
wrappers for \verb+snprintf()+ and \verb+vsnprintf()+, if these functions
are available, otherwise they are simply sprintf() or vsprintf().
It is strongly recommended that \verb+LALSnprintf()+ and \verb+LALVsnprintf()+
be used rather than \verb+sprintf()+ and \verb+vsprintf()+ as the latter are
prone to buffer-overflow problems.

\vfill{\footnotesize\input{LALStdCV}}

</lalLaTeX>
#endif /* autodoc block */

/* DO *NOT* INCLUDE STDIO... STDLIB has size_t */
#include <config.h>
#include <stdlib.h>
#include <stdarg.h>
#ifdef HAVE_VSNPRINTF
int vsnprintf( char *str, size_t size, const char *fmt, va_list ap );
#else /* dangerous: just use vsprintf!!! */
int vsprintf( char *str, const char *fmt, va_list ap );
#define vsnprintf( s, n, f, a ) vsprintf( s, f, a );
#endif

int LALSnprintf( char *str, size_t size, const char *fmt, ... );
int LALVsnprintf( char *str, size_t size, const char *fmt, va_list ap );
int LALSnprintf( char *str, size_t size, const char *fmt, ... )
{
  int n;
  va_list ap;
  va_start( ap, fmt );
  n = vsnprintf( str, size, fmt, ap );
  va_end( ap );
  return n;
}
int LALVsnprintf( char *str, size_t size, const char *fmt, va_list ap )
{
  return vsnprintf( str, size, fmt, ap );
}

/*
 * This needs to be in some source code somewhere,
 * so may as well put it here.
 */
#if defined(NDEBUG) || defined(LAL_NDEBUG)
const int lalNoDebug = 1;
#else
const int lalNoDebug = 0;
#endif
