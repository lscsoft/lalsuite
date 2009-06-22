/*
*  Copyright (C) 2008 Craig Robinson
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

#if 0
<lalVerbatim file="LALInspiralGetApproximantStringCV">
Author: Craig Robinson
$Id$
</lalVerbatim>

<lalLaTeX>
\subsection{Module \texttt{LALInspiralGetApproximantString.c}}

Function for creating the approximant string which gets written to output
files for a given approximant and PN order of the phasing.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{XLALInspiralGetApproximantStringCP}
\idx{XLALXLALInspiralGetApproximantString()}
\begin{itemize}
   \item \texttt{output,} Output, Pointer to the string in which to place the output
   \item \texttt{length,} Input, the length of the output string
   \item \texttt{approx,} Input, enumeration of the waveform approximant
   \item \texttt{order,} Input, post-Newtonian order of the phasing.
\end{itemize}

\subsubsection*{Description}

\subsubsection*{Algorithm}

\subsubsection*{Uses}

\begin{verbatim}
snprintf
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALInspiralRingdownWaveCV}}

</lalLaTeX>
#endif

#include <lal/LALStdlib.h>
#include <lal/LALError.h>

#include <lal/LALInspiral.h>

NRCSID (LALINSPIRALGETAPPROXIMANTSTRINGC, "$Id$");

int XLALInspiralGetApproximantString( CHAR        *output,
                                      UINT4       length,
                                      Approximant approx,
                                      LALPNOrder  order
                                    )
{

  static const char *func = "XLALInspiralGetApproximantString";

  CHAR approxString[LIGOMETA_SEARCH_MAX];
  CHAR orderString[LIGOMETA_SEARCH_MAX];

#ifndef LAL_NDEBUG
  if (!output)
    XLAL_ERROR( func, XLAL_EFAULT );

  if (length < 1)
    XLAL_ERROR( func, XLAL_EINVAL );
#endif

  /* Set the approximant string */

  switch ( approx )
  {
    case TaylorT1:
      snprintf( approxString, LIGOMETA_SEARCH_MAX * sizeof(CHAR),
          "TaylorT1" );
      break;

    case TaylorT2:
      snprintf( approxString, LIGOMETA_SEARCH_MAX * sizeof(CHAR),
          "TaylorT2" );
      break;

    case TaylorT3:
      snprintf( approxString, LIGOMETA_SEARCH_MAX * sizeof(CHAR),
          "TaylorT3" );
      break;

    case TaylorF1:
      snprintf( approxString, LIGOMETA_SEARCH_MAX * sizeof(CHAR),
          "TaylorF1" );
      break;

    case TaylorF2:
      snprintf( approxString, LIGOMETA_SEARCH_MAX * sizeof(CHAR),
          "TaylorF2" );
      break;

    case PadeT1:
      snprintf( approxString, LIGOMETA_SEARCH_MAX * sizeof(CHAR),
          "PadeT1" );
      break;

    case PadeF1:
      snprintf( approxString, LIGOMETA_SEARCH_MAX * sizeof(CHAR),
          "PadeF1" );
      break;

    case EOB:
      snprintf( approxString, LIGOMETA_SEARCH_MAX * sizeof(CHAR),
          "EOB" );
      break;

    case BCV:
      snprintf( approxString, LIGOMETA_SEARCH_MAX * sizeof(CHAR),
          "BCV" );
      break;

    case BCVSpin:
      snprintf( approxString, LIGOMETA_SEARCH_MAX * sizeof(CHAR),
          "BCVSpin" );
      break;

    case SpinTaylorT3:
      snprintf( approxString, LIGOMETA_SEARCH_MAX * sizeof(CHAR),
          "SpinTaylorT3" );
      break;

    case SpinTaylor:
      snprintf( approxString, LIGOMETA_SEARCH_MAX * sizeof(CHAR),
          "SpinTaylor" );
      break;

    case FindChirpSP:
      snprintf( approxString, LIGOMETA_SEARCH_MAX * sizeof(CHAR),
          "FindChirpSP" );
      break;

    case FindChirpPTF:
      snprintf( approxString, LIGOMETA_SEARCH_MAX * sizeof(CHAR),
          "FindChirpPTF" );
      break;

    case GeneratePPN:
      snprintf( approxString, LIGOMETA_SEARCH_MAX * sizeof(CHAR),
          "GeneratePPN" );
      break;

    case BCVC:
      snprintf( approxString, LIGOMETA_SEARCH_MAX * sizeof(CHAR),
          "BCVC" );
      break;

    case Eccentricity:
      snprintf( approxString, LIGOMETA_SEARCH_MAX * sizeof(CHAR),
          "Eccentricity" );
      break;

    case EOBNR:
      snprintf( approxString, LIGOMETA_SEARCH_MAX * sizeof(CHAR),
          "EOBNR" );
      break;
    case TaylorEt:
      snprintf( approxString, LIGOMETA_SEARCH_MAX * sizeof(CHAR),
          "TaylorEt" );
      break;
    case TaylorT4:
      snprintf( approxString, LIGOMETA_SEARCH_MAX * sizeof(CHAR),
          "TaylorT4" );
      break;

    case TaylorN:
      snprintf( approxString, LIGOMETA_SEARCH_MAX * sizeof(CHAR),
          "TaylorN" );
      break;

    default:
      XLALPrintError("Unknown or unsupported approximant.\n");
      XLAL_ERROR( func, XLAL_EINVAL );
      break;
  }


  /* Now set the order */
  switch ( order )
  {
    case LAL_PNORDER_NEWTONIAN:
      snprintf( orderString, LIGOMETA_SEARCH_MAX * sizeof(CHAR),
          "newtonian" );
      break;

    case LAL_PNORDER_HALF:
      snprintf( orderString, LIGOMETA_SEARCH_MAX * sizeof(CHAR),
          "oneHalfPN" );
      break;

    case LAL_PNORDER_ONE:
      snprintf( orderString, LIGOMETA_SEARCH_MAX * sizeof(CHAR),
          "onePN" );
      break;

    case LAL_PNORDER_ONE_POINT_FIVE:
      snprintf( orderString, LIGOMETA_SEARCH_MAX * sizeof(CHAR),
          "onePointFivePN" );
      break;

    case LAL_PNORDER_TWO:
      snprintf( orderString, LIGOMETA_SEARCH_MAX * sizeof(CHAR),
          "twoPN" );
      break;

    case LAL_PNORDER_TWO_POINT_FIVE:
      snprintf( orderString, LIGOMETA_SEARCH_MAX * sizeof(CHAR),
          "twoPointFivePN" );
      break;

    case LAL_PNORDER_THREE:
      snprintf( orderString, LIGOMETA_SEARCH_MAX * sizeof(CHAR),
          "threePN" );
      break;

    case LAL_PNORDER_THREE_POINT_FIVE:
      snprintf( orderString, LIGOMETA_SEARCH_MAX * sizeof(CHAR),
          "threePointFivePN" );
      break;

    case LAL_PNORDER_PSEUDO_FOUR:
      snprintf( orderString, LIGOMETA_SEARCH_MAX * sizeof(CHAR),
          "pseudoFourPN" );
      break;

    default:
      XLALPrintError( "Unknown or unsupported order.\n" );
      XLAL_ERROR( func, XLAL_EINVAL );
      break;
  }

  /* Now build the output and return */
  snprintf( output, length * sizeof(CHAR), "%s%s",
      approxString, orderString );

  return XLAL_SUCCESS;
}
