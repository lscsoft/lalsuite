/*
*  Copyright (C) 2007 Jolien Creighton, John Whelan
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
\author Allen, B.; generalized by J. T. Whelan <jtwhelan@loyno.edu>
\file

\heading{Header \ref PrintVector.h}
\latexonly\label{s_PrintVector_h}\endlatexonly

This is a simple utility to print vectors into a file.

\heading{Synopsis}
\code
#include <lal/PrintVector.h>
\endcode

Contains the prototypes for the
LAL\f$\langle\mbox{DT}\rangle\f$PrintVector functions





*/

#ifndef _PRINTVECTOR_H
#define _PRINTVECTOR_H

#include <lal/LALRCSID.h>

#ifndef _LALSTDLIB_H
#include <lal/LALStdlib.h>
#ifndef _LALSTDLIB_H
#define _LALSTDLIB_H
#endif
#endif

#ifdef  __cplusplus
extern "C" {
#endif


NRCSID( PRINTVECTORH, "$Id$" );

void LALCHARPrintVector( CHARVector *vector );
void LALI2PrintVector( INT2Vector *vector );
void LALI4PrintVector( INT4Vector *vector );
void LALI8PrintVector( INT8Vector *vector );
void LALU2PrintVector( UINT2Vector *vector );
void LALU4PrintVector( UINT4Vector *vector );
void LALU8PrintVector( UINT8Vector *vector );
void LALPrintVector( REAL4Vector *vector );
void LALSPrintVector( REAL4Vector *vector );
void LALDPrintVector( REAL8Vector *vector );
void LALCPrintVector( COMPLEX8Vector *vector );
void LALZPrintVector( COMPLEX16Vector *vector );

#ifdef  __cplusplus
}
#endif

#endif /* _PRINTVECTOR_H */
