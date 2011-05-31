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

/**

\author Tibbits, M. M.
\file

\heading{Header \ref LALMoment.h}
\latexonly\label{s_LALMoment_h}\endlatexonly

\code
The LALDMoment() and LALSMoment() associated header file.
(S - single precision )
(D - double precision )
\endcode

\heading{Synopsis}
\code
#include <lal/LALMoment.h>
\endcode

This header provides the prototype for the LALDMoment() and LALSMoment() function.

*/

/* Double Include Protection */
#ifndef _LALMOMENT_H
#define _LALMOMENT_H

#include <lal/LALStdlib.h>
#include <lal/LALDatatypes.h>


/* C++ protection. */
#ifdef  __cplusplus
extern "C" {
#endif


NRCSID( LALMOMENTH, "$Id$");


/**

\heading{Error codes}



*/

/**\name Error Codes */ /*@{*/

#define LALMOMENTH_ENULL 1
#define	LALMOMENTH_ENNUL 2
#define LALMOMENTH_ELNTH 3
#define LALMOMENTH_ESEGZ 4
#define LALMOMENTH_ENUMZ 5
#define LALMOMENTH_EALOC 6

#define LALMOMENTH_MSGENULL "NULL pointer."
#define	LALMOMENTH_MSGENNUL "Non-NULL pointer."
#define LALMOMENTH_MSGELNTH "Must have more than one data point."
#define LALMOMENTH_MSGESEGZ "Invalid number of segments"
#define LALMOMENTH_MSGENUMZ "Invalid number of points in segment"
#define LALMOMENTH_MSGEALOC "Memory Allocation Error"

/*@}*/


/* Function prototypes */






void LALSMoment
(
	LALStatus		*status,
	REAL4			*result,
	REAL4Sequence		*data,
	INT4			whichMoment
);


void LALDMoment
(
	LALStatus		*status,
	REAL8			*result,
	REAL8Sequence		*data,
	INT4			whichMoment
);





#ifdef  __cplusplus
}
#endif  /* C++ protection. */
#endif  /* Double-include protection. _LALMOMENT_H_ */
