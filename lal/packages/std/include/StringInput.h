/********************************** <lalVerbatim file="StringInputHV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\section{Header \texttt{StringInput.h}}
\label{s:StringInput.h}

Provides routines to parse \verb@CHARVector@s into other LAL
datatypes.

\subsection*{Synopsis}
\begin{verbatim}
#include "StringInput.h"
\end{verbatim}

\noindent This header provides prototypes for routines that construct
LAL data structures using the data from a character string.  As in
standard C, a \emph{string} is a block of non-null bytes of arbitrary
length, terminated by a null byte \verb@'\0'@, and referred to by a
value of type \verb@CHAR *@ pointing to the first byte in the string.
It is not to be confused with a \verb@CHARVector@, a LAL structure
referring to a block of data of a specified length, which may or may
not contain one or more instances of \verb@'\0'@.

In general, the routines under this header will have string inputs of
type \verb@const CHAR *@ (in order to allow, for instance, string
literals to be used as inputs), but will allocate \verb@CHARVector@
structures to store string outputs.  Unless otherwise specified, these
outputs are guaranteed to contain at least one \verb@'\0'@ character,
so their \verb@data@ fields are valid strings.  It is the
responsibility of the calling routine to ensure that the string input
contains a terminating \verb@'\0'@ within the memory segment pointed
to by the \verb@CHAR *@ input, in order to avoid segmentation
violation.

These routines are intended to work in conjunction with the functions
in \verb@StreamInput.h@ to add LAL robustness to otherwise ad-hoc data
input routines.  However, the functions in \verb@StringInput.h@ are
fully LAL-compliant and use only LAL types, so they are included in
\verb@liblal@ proper.

******************************************************* </lalLaTeX> */

#ifndef _STRINGINPUT_H
#define _STRINGINPUT_H

#include <lal/LALStdlib.h>

#ifdef __cplusplus
extern "C" {
#pragma }
#endif

NRCSID( STRINGINPUTH, "$Id$");

/********************************************************** <lalLaTeX>
\subsection*{Error conditions}
****************************************** </lalLaTeX><lalErrTable> */
#define STRINGINPUTH_ENUL 1
#define STRINGINPUTH_EOUT 2
#define STRINGINPUTH_EMEM 3

#define STRINGINPUTH_MSGENUL "Unexpected null pointer in arguments"
#define STRINGINPUTH_MSGEOUT "Output handle points to a non-null pointer"
#define STRINGINPUTH_MSGEMEM "Memory allocation error"
/******************************************** </lalErrTable><lalLaTeX>

\subsection*{Constants}
\idx[Constant]{LAL\_INT2\_FORMAT}
\idx[Constant]{LAL\_INT4\_FORMAT}
\idx[Constant]{LAL\_INT8\_FORMAT}
\idx[Constant]{LAL\_UINT2\_FORMAT}
\idx[Constant]{LAL\_UINT4\_FORMAT}
\idx[Constant]{LAL\_UINT8\_FORMAT}
\idx[Constant]{LAL\_REAL4\_FORMAT}
\idx[Constant]{LAL\_REAL8\_FORMAT}

The following constants are format strings that can be used by the
various C \verb@scanf()@ or \verb@printf()@ functions to parse or
write sequences of characters corresponding to base LAL datatypes.
Since the C datatypes (\verb@short@, \verb@int@, \verb@long@,
\verb@long long@, \verb@float@, \verb@double@, etc.) do not have fixed
mappings to LAL base datatypes (\verb@INT2@, \verb@INT4@, \verb@INT8@,
\verb@REAL4@, \verb@REAL8@, etc.), the appropriate format strings for
each LAL datatype must be determined at configuration time and set at
compile time.

These format strings give only the conversion character preceded by
any length modifier according to the type (\verb@short@, \verb@long@,
etc.).  In particular they do \emph{not} contain the initial
\verb@'%'@ character that initiates the conversion specification.
However, being \verb@#define@d string literals, they can be combined
with \verb@"%"@ string literals or more complicated format strings
through implicit concatenation.  Thus to scan \verb@string@ for a
\verb@UINT4@ number \verb@n@ one would write:
\begin{verbatim}
  sscanf( string, "%" LAL_UINT4_FORMAT, &n );
\end{verbatim}
Similarly, to print a \verb@REAL8@ number \verb@x@ with 12 digits
following the decimal place, one could use the following:
\begin{verbatim}
  printf( "%.12" LAL_REAL8_FORMAT, x );
\end{verbatim}
Of course, floating-point numbers are more commonly printed using the
\verb@"%e"@ conversion specifier, which does not generally require
type-dependent length modifiers.

\begin{center}
\begin{tabular}{|ll|}
\hline
Name & Usual value \\
\hline
\tt LAL\_INT2\_FORMAT  & {\tt "hd"}                \\
\tt LAL\_INT4\_FORMAT  & {\tt "d"}  or {\tt "ld"}  \\
\tt LAL\_INT8\_FORMAT  & {\tt "ld"} or {\tt "lld"} \\
\tt LAL\_UINT2\_FORMAT & {\tt "hu"}                \\
\tt LAL\_UINT4\_FORMAT & {\tt "u"}  or {\tt "lu"}  \\
\tt LAL\_UINT8\_FORMAT & {\tt "lu"} or {\tt "llu"} \\
\tt LAL\_REAL4\_FORMAT & {\tt "f"}                 \\
\tt LAL\_REAL8\_FORMAT & {\tt "lf"}                \\
\hline
\end{tabular}
\end{center}
******************************************************* </lalLaTeX> */
#if 2 == LAL_SIZEOF_SHORT
#define LAL_INT2_FORMAT "hd"
#define LAL_UINT2_FORMAT "hu"
#elif 2 == LAL_SIZEOF_INT
#define LAL_INT2_FORMAT "d"
#define LAL_UINT2_FORMAT "u"
#elif 2 == LAL_SIZEOF_LONG
#define LAL_INT2_FORMAT "ld"
#define LAL_UINT2_FORMAT "lu"
#elif 2 == LAL_SIZEOF_LONG_LONG
#define LAL_INT2_FORMAT "lld"
#define LAL_UINT2_FORMAT "llu"
#else
#define LAL_INT2_FORMAT "d"
#define LAL_UINT2_FORMAT "u"
#endif

#if 4 == LAL_SIZEOF_SHORT
#define LAL_INT4_FORMAT "hd"
#define LAL_UINT4_FORMAT "hu"
#elif 4 == LAL_SIZEOF_INT
#define LAL_INT4_FORMAT "d"
#define LAL_UINT4_FORMAT "u"
#elif 4 == LAL_SIZEOF_LONG
#define LAL_INT4_FORMAT "ld"
#define LAL_UINT4_FORMAT "lu"
#elif 4 == LAL_SIZEOF_LONG_LONG
#define LAL_INT4_FORMAT "lld"
#define LAL_UINT4_FORMAT "llu"
#else
#define LAL_INT4_FORMAT "d"
#define LAL_UINT4_FORMAT "u"
#endif

#if 8 == LAL_SIZEOF_SHORT
#define LAL_INT8_FORMAT "hd"
#define LAL_UINT8_FORMAT "hu"
#elif 8 == LAL_SIZEOF_INT
#define LAL_INT8_FORMAT "d"
#define LAL_UINT8_FORMAT "u"
#elif 8 == LAL_SIZEOF_LONG
#define LAL_INT8_FORMAT "ld"
#define LAL_UINT8_FORMAT "lu"
#elif 8 == LAL_SIZEOF_LONG_LONG
#define LAL_INT8_FORMAT "lld"
#define LAL_UINT8_FORMAT "llu"
#else
#define LAL_INT8_FORMAT "d"
#define LAL_UINT8_FORMAT "u"
#endif

#if 4 == LAL_SIZEOF_FLOAT
#define LAL_REAL4_FORMAT "f"
#elif 4 == LAL_SIZEOF_DOUBLE
#define LAL_REAL4_FORMAT "lf"
#else
#define LAL_REAL4_FORMAT "f"
#endif

#if 8 == LAL_SIZEOF_FLOAT
#define LAL_REAL8_FORMAT "f"
#elif 8 == LAL_SIZEOF_DOUBLE
#define LAL_REAL8_FORMAT "lf"
#else
#define LAL_REAL8_FORMAT "f"
#endif
/********************************************************** <lalLaTeX>

\subsection*{Types}

\subsubsection*{Structure \texttt{TokenList}}
\idx[Type]{TokenList}

This structure stores a number of null-terminated strings of arbitrary
length.  The entire list is stored flattened in a \verb@CHARVector@,
and individual tokens are pointed to by a \verb@CHAR *[]@ handle.  The
fields are:

\begin{description}
\item[\texttt{UINT4 nTokens}] The number of tokens in the list.

\item[\texttt{CHAR **tokens}] A list of pointers to the individual
tokens.  The elements \verb@tokens[0..nTokens-1]@ point to tokens, and
the element \verb@tokens[nTokens]@ is explicitly \verb@NULL@ (as is
the convention for an \verb@argv@ argument list).

\item[\texttt{CHARVector *list}] The flattened list of tokens,
separated by (and terminated with) \verb@'\0'@ characters.
\end{description}

******************************************************* </lalLaTeX> */

typedef struct tagTokenList {
  UINT4 nTokens;    /* number of tokens */
  CHAR **tokens;    /* list of pointers to tokens */
  CHARVector *list; /* flattened list of null-terminated tokens */
} TokenList;

/* <lalLaTeX>
\vfill{\footnotesize\input{StringInputHV}}
</lalLaTeX> */

/* Function prototypes. */

/* <lalLaTeX>
\newpage\input{StringTokenC}
</lalLaTeX> */
void
LALCreateTokenList( LALStatus  *status,
		    TokenList  **list,
		    const CHAR *string,
		    const CHAR *delimiters );

void
LALDestroyTokenList( LALStatus *status,
		     TokenList **list );

/* <lalLaTeX>
\newpage\input{StringConvertC}
</lalLaTeX> */
void
LALStringToU2( LALStatus *status, UINT2 *value, const CHAR *string, CHAR **endptr );

void
LALStringToU4( LALStatus *status, UINT4 *value, const CHAR *string, CHAR **endptr );

void
LALStringToU8( LALStatus *status, UINT8 *value, const CHAR *string, CHAR **endptr );

void
LALStringToI2( LALStatus *status, INT2 *value, const CHAR *string, CHAR **endptr );

void
LALStringToI4( LALStatus *status, INT4 *value, const CHAR *string, CHAR **endptr );

void
LALStringToI8( LALStatus *status, INT8 *value, const CHAR *string, CHAR **endptr );

void
LALStringToS( LALStatus *status, REAL4 *value, const CHAR *string, CHAR **endptr );

void
LALStringToD( LALStatus *status, REAL8 *value, const CHAR *string, CHAR **endptr );

void
LALStringToC( LALStatus *status, COMPLEX8 *value, const CHAR *string, CHAR **endptr );

void
LALStringToZ( LALStatus *status, COMPLEX16 *value, const CHAR *string, CHAR **endptr );

void
LALStringToGPS( LALStatus *status, LIGOTimeGPS *value, const CHAR *string, CHAR **endptr );

#ifdef __cplusplus
#pragma {
}
#endif

#endif /* _STRINGINPUT_H */
