/*
*  Copyright (C) 2007 Jolien Creighton
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

#if 0 /* autodoc block */

<lalVerbatim file="LALStdioHV">
$Id$
</lalVerbatim>

<lalLaTeX>

\section{Header \texttt{LALStdio.h}}
\label{s:LALStdio.h}

Provides LAL functions similar to the non-file functions in \verb+<stdio.h>+.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/LALStdio.h>
#include <lal/FileIO.h>
\end{verbatim}

\noindent This header provides the LALsnprintf function.

\vfill{\footnotesize\input{LALStdioHV}}
\newpage\input{LALStdC}
</lalLaTeX>
#endif /* autodoc block */

#ifndef _LALSTDIO_H
#define _LALSTDIO_H

#include <stdio.h>
#include <stdarg.h>
#include <lal/LALConfig.h>
#include <lal/LALRCSID.h>

#ifdef __cplusplus
extern "C" {
#endif

NRCSID( LALSTDIOH, "$Id$" );

#define LALFopen fopen
#define LALFclose fclose


#if LAL_SIZEOF_SHORT == 2
#define __LAL_INT2_PRI_PREFIX__ ""
#define __LAL_INT2_SCN_PREFIX__ "h"
#elif LAL_SIZEOF_INT == 2
#define __LAL_INT2_PRI_PREFIX__ ""
#define __LAL_INT2_SCN_PREFIX__ ""
#else
#error "ERROR: NO 2 BYTE INTEGER FOUND"
#endif

#if LAL_SIZEOF_INT == 4
#define __LAL_INT4_PRI_PREFIX__ ""
#define __LAL_INT4_SCN_PREFIX__ ""
#elif LAL_SIZEOF_LONG == 4
#define __LAL_INT4_PRI_PREFIX__ "l"
#define __LAL_INT4_SCN_PREFIX__ "l"
#else
#error "ERROR: NO 4 BYTE INTEGER FOUND"
#endif

#if LAL_SIZEOF_LONG == 8
#define __LAL_INT8_PRI_PREFIX__ "l"
#define __LAL_INT8_SCN_PREFIX__ "l"
#elif LAL_SIZEOF_LONG_LONG == 8
#define __LAL_INT8_PRI_PREFIX__ "ll"
#define __LAL_INT8_SCN_PREFIX__ "ll"
#else
#error "ERROR: NO 8 BYTE INTEGER FOUND"
#endif

#define LAL_INT2_PRId __LAL_INT2_PRI_PREFIX__ "d"
#define LAL_INT2_PRIi __LAL_INT2_PRI_PREFIX__ "i"
#define LAL_INT2_PRIo __LAL_INT2_PRI_PREFIX__ "o"
#define LAL_INT2_PRIu __LAL_INT2_PRI_PREFIX__ "u"
#define LAL_INT2_PRIx __LAL_INT2_PRI_PREFIX__ "x"
#define LAL_INT2_PRIX __LAL_INT2_PRI_PREFIX__ "X"

#define LAL_INT4_PRId __LAL_INT4_PRI_PREFIX__ "d"
#define LAL_INT4_PRIi __LAL_INT4_PRI_PREFIX__ "i"
#define LAL_INT4_PRIo __LAL_INT4_PRI_PREFIX__ "o"
#define LAL_INT4_PRIu __LAL_INT4_PRI_PREFIX__ "u"
#define LAL_INT4_PRIx __LAL_INT4_PRI_PREFIX__ "x"
#define LAL_INT4_PRIX __LAL_INT4_PRI_PREFIX__ "X"

#define LAL_INT8_PRId __LAL_INT8_PRI_PREFIX__ "d"
#define LAL_INT8_PRIi __LAL_INT8_PRI_PREFIX__ "i"
#define LAL_INT8_PRIo __LAL_INT8_PRI_PREFIX__ "o"
#define LAL_INT8_PRIu __LAL_INT8_PRI_PREFIX__ "u"
#define LAL_INT8_PRIx __LAL_INT8_PRI_PREFIX__ "x"
#define LAL_INT8_PRIX __LAL_INT8_PRI_PREFIX__ "X"

#define LAL_INT2_SCNd __LAL_INT2_SCN_PREFIX__ "d"
#define LAL_INT2_SCNi __LAL_INT2_SCN_PREFIX__ "i"
#define LAL_INT2_SCNo __LAL_INT2_SCN_PREFIX__ "o"
#define LAL_INT2_SCNu __LAL_INT2_SCN_PREFIX__ "u"
#define LAL_INT2_SCNx __LAL_INT2_SCN_PREFIX__ "x"

#define LAL_INT4_SCNd __LAL_INT4_SCN_PREFIX__ "d"
#define LAL_INT4_SCNi __LAL_INT4_SCN_PREFIX__ "i"
#define LAL_INT4_SCNo __LAL_INT4_SCN_PREFIX__ "o"
#define LAL_INT4_SCNu __LAL_INT4_SCN_PREFIX__ "u"
#define LAL_INT4_SCNx __LAL_INT4_SCN_PREFIX__ "x"

#define LAL_INT8_SCNd __LAL_INT8_SCN_PREFIX__ "d"
#define LAL_INT8_SCNi __LAL_INT8_SCN_PREFIX__ "i"
#define LAL_INT8_SCNo __LAL_INT8_SCN_PREFIX__ "o"
#define LAL_INT8_SCNu __LAL_INT8_SCN_PREFIX__ "u"
#define LAL_INT8_SCNx __LAL_INT8_SCN_PREFIX__ "x"

/* convenient versions of above that can be used in
 * either scanf or printf (decimal integers only) */
#define LAL_INT2_FORMAT  LAL_INT2_SCNd
#define LAL_INT4_FORMAT  LAL_INT4_SCNd
#define LAL_INT8_FORMAT  LAL_INT8_SCNd
#define LAL_UINT2_FORMAT LAL_INT2_SCNu
#define LAL_UINT4_FORMAT LAL_INT4_SCNu
#define LAL_UINT8_FORMAT LAL_INT8_SCNu
#define LAL_REAL4_FORMAT "g"
#define LAL_REAL8_FORMAT "lg"

int
LALSnprintf( char *, size_t, const char *, ... );

int
LALVsnprintf( char *, size_t, const char *, va_list );

#ifdef __cplusplus
}
#endif

#endif /* _LALSTDIO_H */
