/*
*  Copyright (C) 2007 Jolien Creighton
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

/* <lalVerbatim> */
#ifndef _LALPRIMER_H  /* Double-include protection. */
#define _LALPRIMER_H

#include <lal/LALDatatypes.h>

#ifdef  __cplusplus   /* C++ protection. */
extern "C" {
#endif

#define LALPRIMERH_ENULL 1
#define LALPRIMERH_EDIV0 2
#define LALPRIMERH_MSGENULL "Null pointer"
#define LALPRIMERH_MSGEDIV0 "Division by zero"

void
REAL4Invert( LALStatus *status, REAL4 *output, REAL4 input );

void
REAL4Divide( LALStatus *status, REAL4 *output, REAL4 numer, REAL4 denom);


#ifdef  __cplusplus
}
#endif  /* C++ protection. */
#endif  /* Double-include protection. */
/* </lalVerbatim> */
