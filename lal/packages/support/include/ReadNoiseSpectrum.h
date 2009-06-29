/*
*  Copyright (C) 2007 Jolien Creighton, Patrick Brady
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

/*-----------------------------------------------------------------------
 *
 * File Name: ReadNoiseSpectrum.h
 *
 * Author: Brady, P. R.
 *
 * Revision: $Id$
 *
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="ReadNoiseSpectrumHV">
Author: Brady, P. R.
$Id$
</lalVerbatim>
<lalLaTeX>
\section{Header \texttt{ReadNoiseSpectrum.h}}
\label{s:ReadNoiseSpectrum.h}

Provides function to read in a file containing a possibly unequally sampled
noise amplitude spectrum ($\textrm{strain}/\sqrt(\textrm{Hz})$) and return as
a frequency series.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/ReadNoiseSpectrum.h>
\end{verbatim}

</lalLaTeX>
#endif

#ifndef _READNOISESPECTRUMH_H
#define _READNOISESPECTRUMH_H

#include <lal/LALDatatypes.h>
#include <lal/Date.h>
#include <lal/TwoDMesh.h>

#ifdef  __cplusplus
extern "C" {
#pragma }
#endif


NRCSID (READNOISESPECTRUMH, "$Id$");

/* <lalLaTeX>
\newpage\subsection*{Error codes}
</lalLaTeX> */
/* <lalErrTable> */
#define LALREADNOISESPECTRUMH_ENULL 1
#define LALREADNOISESPECTRUMH_ENNUL 2
#define LALREADNOISESPECTRUMH_EALOC 3
#define LALREADNOISESPECTRUMH_EOPEN 4
#define LALREADNOISESPECTRUMH_EFCLO 5
#define LALREADNOISESPECTRUMH_EPARS 8

#define LALREADNOISESPECTRUMH_MSGENULL "Null pointer"
#define LALREADNOISESPECTRUMH_MSGENNUL "Non-null pointer"
#define LALREADNOISESPECTRUMH_MSGEALOC "Memory allocation error"
#define LALREADNOISESPECTRUMH_MSGEOPEN "Error opening file"
#define LALREADNOISESPECTRUMH_MSGEFCLO "Error closing file"
#define LALREADNOISESPECTRUMH_MSGEPARS "Error parsing spectrum file"
/* </lalErrTable> */

#define LALREADNOISESPECTRUM_MAXLINELENGTH 2048

void LALReadNoiseSpectrum(
    LALStatus *status,
    REAL4FrequencySeries *spectrum,
    CHAR *fname
    );

#if 0
<lalLaTeX>
\vfill{\footnotesize\input{ReadNoiseSpectrumHV}}
\newpage\input{ReadNoiseSpectrumC}
</lalLaTeX> */
#endif

#ifdef  __cplusplus
#pragma {
}
#endif

#endif /* _READNOISESPECTRUMH_H */
