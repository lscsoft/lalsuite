/*
*  Copyright (C) 2007 Badri Krishnan, Jolien Creighton
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
 * File Name: RngMedBias.h
 *
 * Authors:  Krishnan, B.
 *
 * Revision: $Id$
 *
 * History:   Created by Krishnan Feb 24, 2004
 *
 *
 *-----------------------------------------------------------------------
 */

/* *********************************** <lalVerbatim file="RngMedBiasHV">
Author: Krishnan, B., Itoh, Y.
$Id$
************************************* </lalVerbatim> */

/* <lalLaTeX>  *********************************************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\section{Header \texttt{RngMedBias.h}}
\label{s:RngMedBias.h}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection*{Synopsis}

\begin{verbatim}
#include <lal/RngMedBias.h>
\end{verbatim}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection*{Error conditions}
\vspace{0.1in}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\vfill{\footnotesize\input{RngMedBiasHV}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\newpage\input{RngMedBiasC}
%%%%%%%%%% Test program. %%
\newpage\input{RngMedBiasTestC}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

*************************************************</lalLaTeX> */

/*
 * 4.  Protection against double inclusion (include-loop protection)
 *     Note the naming convention!
 */

#ifndef _RNGMEDBIAS_H
#define _RNGMEDBIAS_H

/*
 * 5. Includes. This header may include others; if so, they go immediately
 *    after include-loop protection. Includes should appear in the following
 *    order:
 *    a. Standard library includes
 *    b. LDAS includes
 *    c. LAL includes
 */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <lal/LALStdlib.h>

/*
 *   Protection against C++ name mangling
 */

#ifdef  __cplusplus
extern "C" {
#endif

 /*
 * 6. Assignment of Id string using NRCSID()
 */

NRCSID (RNGMEDBIASH, "$Id$");

/*
 * 7. Error codes and messages. This must be auto-extracted for
 *    inclusion in the documentation.
 */

/* <lalErrTable file="RngMedBiasHErrorTable"> */

#define RNGMEDBIASH_ENULL 1
#define RNGMEDBIASH_EVAL 5

#define RNGMEDBIASH_MSGENULL "Null pointer"
#define RNGMEDBIASH_MSGEVAL  "Invalid value"

/* </lalErrTable>  */


/* ******************************************************
 * 8. Macros. But, note that macros are deprecated.
 *    They could be moved to the modules where are needed
 */


/* *******************************************************
 * 9. Constant Declarations. (discouraged)
 */



/* **************************************************************
 * 10. Structure, enum, union, etc., typdefs.
 */

/*
 * 11. Extern Global variables. (discouraged)
 */


/*
 * 12. Functions Declarations (i.e., prototypes).
 */


void LALRngMedBias (LALStatus   *status,
		 REAL8       *biasFactor,
		 INT4        blkSize
		 );

#ifdef  __cplusplus
}                /* Close C++ protection */
#endif

#endif     /* Close double-include protection _RNGMEDBIAS_H */








