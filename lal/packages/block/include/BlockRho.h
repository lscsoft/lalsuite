/*
*  Copyright (C) 2007 Jolien Creighton, Matt Tibbits
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

/****************** <lalVerbatim file="BlockRhoHV">
Author: Tibbits, M. M.
$Id$
********************************* </lalVerbatim> */

/******************************* <lalLaTeX file="BlockRhoH">
\section{Header \texttt{BlockRho.h}}
\label{s:BlockRho.h}

Provides structures and definitions global to the Rho algorithms,
notably \texttt{BlockSearchParams}.

\subsection*{Error conditions}

\input{BlockRhoHErrTab}

\subsection*{Structures and Unions}

*************************************************** </lalLaTeX> */

#ifndef _BLOCK_RHO_H_
#define _BLOCK_RHO_H_

#include <math.h>
#include <lal/LALRCSID.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/LALConstants.h>
#include <lal/LALMoment.h>
#include <lal/Matrix.h>

/* C++ protection */
#ifdef  __cplusplus
extern "C" {
#endif

NRCSID( BLOCKRHOH, "$Id$");


/******** <lalErrTable file="BlockRhoHErrTab"> ********/
#define BLOCKRHOH_ENULL 1
#define BLOCKRHOH_ENNUL 2
#define BLOCKRHOH_EALOC 3
#define BLOCKRHOH_ENUMZ 4
#define	BLOCKRHOH_EARG  5
#define BLOCKRHOH_EDATA 6

#define BLOCKRHOH_MSGENULL "Null pointer"
#define BLOCKRHOH_MSGENNUL "Non-null pointer"
#define BLOCKRHOH_MSGEALOC "Memory allocation error"
#define BLOCKRHOH_MSGENUMZ "Data segment length is zero"
#define	BLOCKRHOH_MSGEARG  "Error parsing command-line arguments"
#define BLOCKRHOH_MSGEDATA "Too few input data points to define a Rho statistic"

/******** </lalErrTable> ********/

void LALBlockRho2 (
 LALStatus		*status,
 REAL8			*result,
 REAL8			*rpeak,
 INT4			*myindex,
 REAL8Sequence		*data,
 UINT4			*marginOfExclusion			
);

void LALBlockRho3 (
 LALStatus               *status,
 REAL8                   *result,
 REAL8                   *rpeak,
 INT4                    *myindex,
 REAL8Sequence           *data
);


/* C++ protection */
#ifdef  __cplusplus
}
#endif

/******************************* <lalLaTeX file="BlockRhoH">
  \vfill{\footnotesize\input{BlockRhoHV}}
  \newpage\input{BlockRho2C}
  %\newpage\input{BlockRho3C}
*************************************************** </lalLaTeX> */

#endif /* _BLOCK_RHO_H_ */
