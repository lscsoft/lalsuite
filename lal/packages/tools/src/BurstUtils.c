/*----------------------------------------------------------------------- 
 * 
 * File Name: BurstUtils.c
 *
 * Author: Ray Majumder, S. K.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="BurstUtilsCV">
Author: Brown, D. A. Brady, P. R. Ray Majumder, S. K and Cannon, K. C.
$Id$
</lalVerbatim> 
#endif

#include <stdio.h>
#include <stdlib.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/Date.h>
#include <lal/LALInspiral.h>
#include <lal/BurstUtils.h>

NRCSID( SNGLBURSTUTILSC, "$Id$" );

#if 0
<lalLaTeX>
\subsection{Module \texttt{BurstUtils.c}}

\noindent Blah.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{BurstUtilsCP}
\idx{XLALMergerFrequency()}
\idx{XLALQnrFrequency()}
\idx{XLALMergerDuration()}
\idx{XLALMergerEnergyFraction()}
\idx{XLALMergerhrss()}

\subsubsection*{Description}

\noindent Blah.

\subsubsection*{Algorithm}

\noindent None.

\subsubsection*{Uses}

\noindent None.

\subsubsection*{Notes}
%% Any relevant notes.
 
\vfill{\footnotesize\input{BurstUtilsCV}}

</lalLaTeX>
#endif



/*
 * Sort a list of SnglBurstTable events into increasing order according to the
 * supplied comparison function.
 */

/* <lalVerbatim file="BurstUtilsCP"> */
REAL4
XLALMergerFrequency(
	InspiralTemplate *params
)
/* </lalVerbatim> */
{
	REAL4 fmerger;
    
	fmerger = 205 * (20/params->totalMass);

	return(fmerger);
}

/* <lalVerbatim file="BurstUtilsCP"> */
REAL4
XLALQnrFrequency(
	InspiralTemplate *params
)
/* </lalVerbatim> */
{
	REAL4 fqnr;
    
	fqnr = 1320 * (20/params->totalMass);

	return(fqnr);
}

REAL4
XLALMergerEnergyFraction(
	InspiralTemplate *params
)
/* </lalVerbatim> */
{
	REAL4 e_merger = 0.0;
    
	e_merger = 0.1 * 16 * params->mu * params->mu /params->totalMass;

	return(e_merger);
}

void 
XLALMergerDuration(
	InspiralTemplate *params,
	MergerDurationLimits *dur
)
/* </lalVerbatim> */
{
	dur->high = 50 * params->totalMass;
	dur->low = 10 * params->totalMass;
}
