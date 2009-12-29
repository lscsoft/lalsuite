/*
 *  Copyright (C) 2005 Badri Krishnan, Alicia Sintes
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

/**
 *  \file Velocity.h
 *  \author Badri Krishnan
 *  \brief Header file for velocity and position routines.
 *  \ingroup moduleHoughPulsar
 *
 *  $Id$
 *
 *  Created by Badri Krishnan on May 23, 2003
***************************************************** */


/************************************<lalVerbatim file="VelocityHV">
Author: Krishnan, B., Sintes, A.M.
$Id$
*************************************</lalVerbatim> */

/* <lalLaTeX> *************************************************

\section{Header \texttt{Velocity.h}}
\label{s:Velocity.h}
Computation of instant and averaged velocities for a given detector and the
like.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection*{Synopsis}

\begin{verbatim}
#include <lal/Velocity.h>
\end{verbatim}

To find the velocity of a given detetector at a given time, or the averaged
velocity  of a detector in a certain time interval.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection*{Error conditions}
\vspace{0.1in}
\input{VelocityHErrorTable}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection*{Structures}
\begin{verbatim}
struct VelocityPar
\end{verbatim}
\index{\texttt{VelocityPar}}
\noindent This structure stores the parameters required by LALBarycenter to calculate
Earth velocity at a given detector location.
\begin{description}
\item[\texttt{LALDetector    detector}]
\item[\texttt{EphemerisData  *edat}]  ephemeris data pointer from LALInitBarycenter
\item[\texttt{LIGOTimeGPS  startTime}]  start of time interval
\item[\texttt{REAL8          tBase}]  duration of interval
\item[\texttt{REAL8          vTol}]   fractional accuracy required for velocity
\end{description}


\begin{verbatim}
struct AvgVelPar
\end{verbatim}
\index{\texttt{AvgVelPar}}
\begin{description}
\item[\texttt{LALDetector    detector}]
\item[\texttt{EphemerisData  *edat}]
\end{description}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\vfill{\footnotesize\input{VelocityHV}}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\newpage\input{VelocityC}
%%%%%%%%%%Test program. %%
\newpage\input{TestVelocityC}

*************************************************** </lalLaTeX> */


/* double inclusion protection */
#ifndef _VELOCITY_H
#define _VELOCITY_H

/* *************
 *    Includes. This header may include others; if so, they go immediately
 *    after include-loop protection. Includes should appear in the following
 *    order:
 *    a. Standard library includes
 *    b. LDAS includes
 *    c. LAL includes
 */
#include<lal/Date.h>
#include<lal/LALDatatypes.h>
#include<lal/ComputeSky.h>
#include<lal/LALInitBarycenter.h>
#include<lal/LALBarycenter.h>
#include<lal/LALStdlib.h>
#include<lal/LALConstants.h>

/*  ****************************************************
 *   Protection against C++ name mangling
 */

#ifdef  __cplusplus
extern "C" {
#endif
/* ***************************************
 *   Assignment of Id string using NRCSID()
 */

NRCSID( VELOCITYH, "$Id$");

/* ***************************************
 *   Error codes and messages. This must be auto-extracted for
 *    inclusion in the documentation.
 */
/* <lalErrTable file="VelocityHErrorTable"> */
#define VELOCITYH_ENULL 1
#define VELOCITYH_EVAL 2
#define VELOCITYH_MSGENULL "Null Pointer"
#define VELOCITYH_MSGEVAL "Invalid Value"
/* </lalErrTable>  */

/* *****************************************************
 *   Structure, enum, union, etc., typdefs.
 */

/** parameters required by LALBarycenter to calculate Earth velocity at
   a given detector location */
typedef struct tagVelocityPar {
  LALDetector    detector; /**< the detector */
  EphemerisData  *edat;  /**< ephemeris data pointer from LALInitBarycenter */
  LIGOTimeGPS    startTime; /**< start of time interval */
  REAL8          tBase; /**< duration of interval */
  REAL8          vTol;  /**< fractional accuracy required for velocity (redundant for average velocity calculation) */
} VelocityPar;


/* ***************************************************
 *  Functions Declarations (i.e., prototypes).
 */
void LALAvgDetectorVel(LALStatus    *status,
		    REAL8        v[3], /* output vector representing average velocity */
		    VelocityPar  *in); /* parameters required to calculate V */

void LALAvgDetectorPos(LALStatus    *status,
		    REAL8        x[3], /* output vector representing average position */
		    VelocityPar  *in); /* parameters required to calculate position */

void LALDetectorVel(LALStatus   *status,
		 REAL8       v[3],  /* output velocity vector */
		 LIGOTimeGPS *time0, /* time at which velocity is calculated */
		 LALDetector  detector, /* detector */
		 EphemerisData *edat);

void LALDetectorPos(LALStatus   *status,
		 REAL8       x[3],  /* output position vector */
		 LIGOTimeGPS *time0, /* time at which position is calculated */
		 LALDetector  detector, /* detector*/
		 EphemerisData *edat);

/* ****************************************************** */

#ifdef  __cplusplus
}                /* Close C++ protection */
#endif

#endif  /* end of double inclusion protection */

