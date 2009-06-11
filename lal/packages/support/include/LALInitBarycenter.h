/*
*  Copyright (C) 2007 Curt Cutler, Reinhard Prix, Teviet Creighton, Bernd Machenschalk
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

/********************************** <lalVerbatim file="LALInitBarycenterHV">
Author: Cutler, C.
$Id$
*********************************** </lalVerbatim> */

/* <lalLaTeX>

\section{Header \texttt{LALInitBarycenter.h}}
\label{s:LALInitBarycenter.h}

Provides a routine for reading Earth and Sun position information from
data files.


\subsection*{Synopsis}

\begin{verbatim}
#include "LALInitBarycenter.h"
\end{verbatim}

\noindent This header covers the routine \verb@LALInitBarycenter.c@.
Since it involves file I/O, it is placed in the \verb@support@
package, and included in the \verb@lalsupport@ library.

</lalLaTeX> */



#ifndef _LALINITBARYCENTER_H    /* Protect against double-inclusion */
#define _LALINITBARYCENTER_H

#include <lal/LALBarycenter.h>
#ifdef  __cplusplus
extern "C" {
#endif

NRCSID (LALINITBARYCENTERH,"$Id$");

/* <lalErrTable file="LALInitBarycenterHErrorTable"> */
#define LALINITBARYCENTERH_EOPEN    1
#define LALINITBARYCENTERH_EMEM     2
#define LALINITBARYCENTERH_EEPHFILE 32

#define LALINITBARYCENTERH_MSGEOPEN    "Could not open ephemeris file"
#define LALINITBARYCENTERH_MSGEMEM     "Out of memory"
#define LALINITBARYCENTERH_MSGEEPHFILE "Error in reading an ephemeris file"
/* </lalErrTable> */

/* <lalLaTeX>
\subsection*{Error conditions}
\vspace{0.1in}
\input{LALInitBarycenterHErrorTable}
</lalLaTeX> */

/* Function prototypes. */

void LALInitBarycenter(LALStatus *, EphemerisData *);


/********************************************************** <lalLaTeX>
\vfill{\footnotesize\input{LALInitBarycenterHV}}
\newpage\input{LALInitBarycenterC}
******************************************************* </lalLaTeX> */

#ifdef  __cplusplus
}
#endif      /* Close C++ protection */

#endif      /* Close double-include protection */























