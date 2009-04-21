/*
*  Copyright (C) 2007 Craig Robinson 
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
 * File Name: TrigScanEThincaCommon.h
 *
 * Author: Robinson, C. A. K.
 * 
 * Revision: $Id: LIGOMetadataUtils.h,v 1.156 2008/08/29 11:20:58 sfairhur Exp $
 * 
 *-----------------------------------------------------------------------
 */

#include <lal/CoincInspiralEllipsoid.h>

#if 0
<lalVerbatim file="TrigScanEThincaCommonHV">
Author: Robinson, C. A. K.
$Id: LIGOMetadataUtils.h,v 1.156 2008/08/29 11:20:58 sfairhur Exp $
</lalVerbatim> 
<lalLaTeX>
\section{Header \texttt{TrigScanEThincaCommon.h}}
\label{s:TrigScanEThincaCommon.h}

Provides helper functions common to TrigScan and E-thinca.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/TrigScanEThincaCommon.h>
\end{verbatim}

\noindent This header provides functions used for creating and destroying the
linked lists used in TrigScan and E-thinca.

</lalLaTeX>
#endif

#ifndef _TRIGSCANETHINCACOMMON_H
#define _TRIGSCANETHINCACOMMON_H

#ifdef  __cplusplus
extern "C" {
#endif

NRCSID( TRIGSCANETHINCACOMMONH, "$Id: CoincInspiralEllipsoid.h,v 1.7 2007/06/08 14:41:56 bema Exp $" );

TriggerErrorList * XLALCreateTriggerErrorList( SnglInspiralTable *tableHead,
                                               REAL8             scaleFactor,
                                               REAL8             *tcMax );

void XLALDestroyTriggerErrorList( TriggerErrorList *errorListHead );

#ifdef  __cplusplus
}
#endif

#endif /* _TRIGSCANETHINCACOMMON_H */
