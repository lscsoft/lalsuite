/*
*  Copyright (C) 2007 Alexander Dietz, Robert Adam Mercer
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




#ifndef _SKYANNULUS_H  /* Double-include protection. */
#define _SKYANNULUA_H

#if 0
<lalLaTeX>
\subsection{Module \texttt{SkyAnnulus.c}}

\noindent
Provides methods to calculate the annulus given two or three detectors with a certain time differe
nce.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{SkyAnnulusCP}
\idx{LALComputeSingleAnnulus(
                   LALStatus          *status,
                   INT2Vector         *sky,
                   LALPlaceAndGPS     *placetime1,
                   LALPlaceAndGPS     *placetime2,
                   REAL4              dtErr,
                   INT2               nx,
                   INT2               ny)}
\idx{LALComputeTripleAnnuli( LALStatus          *status,
                       INT2Vector         **sky,
                       LALPlaceAndGPS     *placetime1,
                       LALPlaceAndGPS     *placetime2,
                       LALPlaceAndGPS     *placetime3,
                       REAL4              dtErr,
                       INT2               nx,
                       INT2               ny)}
\idx{LALExistTripleAnnuli( LALStatus
                           *status,
                     INT2               *exist,
                     LALPlaceAndGPS     *placetime1,
                     LALPlaceAndGPS     *placetime2,
                     LALPlaceAndGPS     *placetime3,
                     REAL4              dtErr,
                     INT2               nx,
                     INT2               ny)}

\subsubsection*{Description}

%\noident
\texttt{LALComputeSingleAnnulus} calculates the single annulus when two detectors and two times are given in \texttt{placetime1} and \texttt{placetime2} assuming an error in time resolution of \texttt{dtErr}. The result is a vector \texttt{sky} with dimensions \texttt{nx}x\texttt{ny} representing the whole sky. The index of the vector ($\mathrm{index} = i + j\cdot nx$) represents the coordinates ($\alpha=360^\circ \cdot i / nx$ and $\delta= 180^\circ \cdot i / ny -90$.
If the desired position lies within the annulus the element of the vector is 1, 0 else.
\texttt{LALComputeTripleAnnuli} calculates the annuli when {\em three} detectors and times are given.
The result is a vector of four vectors, with the first vector (index 0) containing the complete overlap of all annuli and the next three vectors (indixes 1,2,3) containing the annuli of two of the time intervals each.
\texttt{LALExistTripleAnnuli} just checks if there is at least one position in the sky overlapped by all three annuli created with \texttt{LALComputeTripleAnnuli}.
If that is the case, \texttt{exist} will contain the value 1, 0 else.


\subsubsection*{Algorithm}

%\noindent

\subsubsection*{Uses}

																																																		\noindent For each module it is not needed to allocate memory for the structure \texttt{sky}, this is done in the code. The only array needed for \texttt{LALComputeTripleAnnuli} is an array \texttt{sky} of pointers to \texttt{INT2Vector}  with 4 elements.

\subsubsection*{Notes}
%% Any relevant notes.

\vfill{\footnotesize\input{SkyAnnulusCV}}

</lalLaTeX>
#endif


#include <lal/DetectorSite.h>
#include <lal/SkyCoordinates.h>
#include <lal/Date.h>
#include <lal/TimeDelay.h>

NRCSID( COMPUTESKYPOSH, "$Id$" );



void
LALComputeSingleAnnulus (
		   LALStatus          *status,
		   INT2Vector         *sky,
		   LALPlaceAndGPS     *placetime1,
		   LALPlaceAndGPS     *placetime2,
		   REAL4              dtErr,
		   INT2               nx,
		   INT2               ny
		   );

void
LALComputeTripleAnnuli (
		       LALStatus          *status,
		       INT2Vector         **sky,
		       LALPlaceAndGPS     *placetime1,
		       LALPlaceAndGPS     *placetime2,
		       LALPlaceAndGPS     *placetime3,
		       REAL4              dtErr,
		       INT2               nx,
		       INT2               ny
		       );

void
LALExistTripleAnnuli(
		     LALStatus          *status,
		     INT2               *exist,
		     LALPlaceAndGPS     *placetime1,
		     LALPlaceAndGPS     *placetime2,
		     LALPlaceAndGPS     *placetime3,
		     REAL4              dtErr,
		     INT2               nx,
		     INT2               ny
		     );



#endif
