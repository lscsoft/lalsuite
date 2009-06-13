/*
*  Copyright (C) 2008 Santamaria L, Krishnan B, Whelan JT, Dias M, Parameswaran A
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
 * File Name: FindChirpPhenomCoeffs.h
 *
 * Author: Santamaria L, Krishnan B, Whelan JT, Dias M, Parameswaran A
 *
 * Revision: $Id$
 *
 *-----------------------------------------------------------------------
 */

/*
<lalVerbatim file="FindChirpPhenomCoeffsHV">
Author: Santamaria L, Krishnan B, Whelan JT, Dias M, Parameswaran A.
$Id$
</lalVerbatim>

<lalLaTeX>
\section{Header \texttt{FindChirpPhenomCoeffs.h}}
\label{s:FindChirpPhenomCoeffs.h}

Provides co\"{e}fficients for the phenomenological waveforms
introduced by Ajith et al. in arXiv:0710.2335 [gr-qc]

\subsection*{Synopsis}

\begin{verbatim}
#include <lal/FindChirpPhenomCoeffs.h>
\end{verbatim}

\input{FindChirpPhenomCoeffsDoc}
</lalLaTeX>
*/

#ifndef _BBHPHENOMCOEFFSH_H
#define _BBHPHENOMCOEFFSH_H

#ifdef  __cplusplus
extern "C" {
#pragma }
#endif

NRCSID (BBHPHENOMCOEFFSH, "$Id$");

/* This header contains the coeffs from the matching with the LONG */
/* Jena waveforms (those are not the ones published in the original paper */
/* but in the Amaldi 07 proceedings) */

#define BBHPHENOMCOEFFSH_FMERG_A   6.6389e-01
#define BBHPHENOMCOEFFSH_FMERG_B   -1.0321e-01
#define BBHPHENOMCOEFFSH_FMERG_C   1.0979e-01

#define BBHPHENOMCOEFFSH_FRING_A   1.3278e+00
#define BBHPHENOMCOEFFSH_FRING_B   -2.0642e-01
#define BBHPHENOMCOEFFSH_FRING_C   2.1957e-01

#define BBHPHENOMCOEFFSH_SIGMA_A   1.1383e+00
#define BBHPHENOMCOEFFSH_SIGMA_B   -1.7700e-01
#define BBHPHENOMCOEFFSH_SIGMA_C   4.6834e-02

#define BBHPHENOMCOEFFSH_FCUT_A   1.7086e+00
#define BBHPHENOMCOEFFSH_FCUT_B   -2.6592e-01
#define BBHPHENOMCOEFFSH_FCUT_C   2.8236e-01

#define BBHPHENOMCOEFFSH_PSI0_X   -1.5829e-01
#define BBHPHENOMCOEFFSH_PSI0_Y   8.7016e-02
#define BBHPHENOMCOEFFSH_PSI0_Z   -3.3382e-02

#define BBHPHENOMCOEFFSH_PSI2_X   3.2967e+01
#define BBHPHENOMCOEFFSH_PSI2_Y   -1.9000e+01
#define BBHPHENOMCOEFFSH_PSI2_Z   2.1345e+00

#define BBHPHENOMCOEFFSH_PSI3_X   -3.0849e+02
#define BBHPHENOMCOEFFSH_PSI3_Y   1.8211e+02
#define BBHPHENOMCOEFFSH_PSI3_Z   -2.1727e+01

#define BBHPHENOMCOEFFSH_PSI4_X   1.1525e+03
#define BBHPHENOMCOEFFSH_PSI4_Y   -7.1477e+02
#define BBHPHENOMCOEFFSH_PSI4_Z   9.9692e+01

#define BBHPHENOMCOEFFSH_PSI6_X   1.2057e+03
#define BBHPHENOMCOEFFSH_PSI6_Y   -8.4233e+02
#define BBHPHENOMCOEFFSH_PSI6_Z   1.8046e+02

#define BBHPHENOMCOEFFSH_PSI7_X   -0.0000e+00
#define BBHPHENOMCOEFFSH_PSI7_Y   0.0000e+00
#define BBHPHENOMCOEFFSH_PSI7_Z   0.0000e+00

#if 0
<lalLaTeX>
\vfill{\footnotesize\input{FindChirpPhenomCoeffsHV}}
</lalLaTeX>
#endif

#ifdef  __cplusplus
#pragma {
}
#endif

#endif /* _BBHPHENOMCOEFFSH_H */
