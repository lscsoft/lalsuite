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
