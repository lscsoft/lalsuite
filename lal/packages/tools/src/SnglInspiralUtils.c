/*----------------------------------------------------------------------- 
 * 
 * File Name: SnglInspiralUtils.c
 *
 * Author: Brown, D. A.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="SnglInspiralUtilsCV">
Author: Brown, D. A.
$Id$
</lalVerbatim> 
#endif

#include <stdio.h>
#include <stdlib.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>

NRCSID( SNGLINSPIRALUTILSC, "$Id$" );

#if 0
<lalLaTeX>
\subsection{Module \texttt{SnglInspiralUtils.c}}

\noindent Blah.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{SnglInspiralUtilsCP}
\idx{LALSortSnglInspiral()}
\idx{LALCompareSnglInspiralByMass()}
\idx{LALCompareSnglInspiralByTime()}

\subsubsection*{Description}

\noindent Blah.

\subsubsection*{Algorithm}

\noindent None.

\subsubsection*{Uses}

\noindent None.

\subsubsection*{Notes}
%% Any relevant notes.
 
\vfill{\footnotesize\input{SnglInspiralUtilsCV}}

</lalLaTeX>
#endif

/* <lalVerbatim file="SnglInspiralUtilsCP"> */
void
LALSortSnglInspiral(
    LALStatus          *status,
    SnglInspiralTable **eventHead,
    int(*comparfunc)    (const void *, const void *)
    )
/* </lalVerbatim> */
{
}


/* <lalVerbatim file="SnglInspiralUtilsCP"> */
int
LALCompareSnglInspiralByMass(
    const void *a,
    const void *b
    )
/* </lalVerbatim> */
{
}


/* <lalVerbatim file="SnglInspiralUtilsCP"> */
int
LALCompareSnglInspiralByTime(
    const void *a,
    const void *b
    )
/* </lalVerbatim> */
{
}
