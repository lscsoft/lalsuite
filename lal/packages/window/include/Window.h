/*
 * Copyright (C) 2007 Bruce Allen, Duncan Brown, Jolien Creighton, Kipp
 * Cannon, Teviet Creighton
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with with program; see the file COPYING. If not, write to the Free
 * Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 */

#ifndef _WINDOW_H
#define _WINDOW_H

#include <lal/LALDatatypes.h>

#ifdef  __cplusplus
extern "C" {
#endif


NRCSID (WINDOWH, "$Id$");


typedef struct tagREAL4Window {
	REAL4Sequence *data;
	REAL8          sumofsquares;
	REAL8          sum;
} REAL4Window;


typedef struct tagREAL8Window {
	REAL8Sequence *data;
	REAL8          sumofsquares;
	REAL8          sum;
} REAL8Window;


REAL4Window *XLALCreateREAL4WindowFromSequence(REAL4Sequence *sequence);
REAL8Window *XLALCreateREAL8WindowFromSequence(REAL8Sequence *sequence);

REAL4Window *XLALCreateRectangularREAL4Window(UINT4 length);
REAL4Window *XLALCreateHannREAL4Window(UINT4 length);
REAL4Window *XLALCreateWelchREAL4Window(UINT4 length);
REAL4Window *XLALCreateBartlettREAL4Window(UINT4 length);
REAL4Window *XLALCreateParzenREAL4Window(UINT4 length);
REAL4Window *XLALCreatePapoulisREAL4Window(UINT4 length);
REAL4Window *XLALCreateHammingREAL4Window(UINT4 length);
REAL4Window *XLALCreateKaiserREAL4Window(UINT4 length, REAL4 beta);
REAL4Window *XLALCreateCreightonREAL4Window(UINT4 length, REAL4 beta);
REAL4Window *XLALCreateTukeyREAL4Window(UINT4 length, REAL4 beta);
REAL4Window *XLALCreateGaussREAL4Window(UINT4 length, REAL4 beta);

REAL8Window *XLALCreateRectangularREAL8Window(UINT4 length);
REAL8Window *XLALCreateHannREAL8Window(UINT4 length);
REAL8Window *XLALCreateWelchREAL8Window(UINT4 length);
REAL8Window *XLALCreateBartlettREAL8Window(UINT4 length);
REAL8Window *XLALCreateParzenREAL8Window(UINT4 length);
REAL8Window *XLALCreatePapoulisREAL8Window(UINT4 length);
REAL8Window *XLALCreateHammingREAL8Window(UINT4 length);
REAL8Window *XLALCreateKaiserREAL8Window(UINT4 length, REAL8 beta);
REAL8Window *XLALCreateCreightonREAL8Window(UINT4 length, REAL8 beta);
REAL8Window *XLALCreateTukeyREAL8Window(UINT4 length, REAL8 beta);
REAL8Window *XLALCreateGaussREAL8Window(UINT4 length, REAL8 beta);

void XLALDestroyREAL4Window(REAL4Window *window);
void XLALDestroyREAL8Window(REAL8Window *window);

REAL4Sequence *XLALUnitaryWindowREAL4Sequence(REAL4Sequence *sequence, const REAL4Window *window);
COMPLEX8Sequence *XLALUnitaryWindowCOMPLEX8Sequence(COMPLEX8Sequence *sequence, const REAL4Window *window);
REAL8Sequence *XLALUnitaryWindowREAL8Sequence(REAL8Sequence *sequence, const REAL8Window *window);
COMPLEX16Sequence *XLALUnitaryWindowCOMPLEX16Sequence(COMPLEX16Sequence *sequence, const REAL8Window *window);

/*
 * LEGACY CODE:  DO NOT USE!!!
 *
 * FIXME: remove as soon as possible
 */


typedef enum {
	Rectangular,
	Hann,
	Welch,
	Bartlett,
	Parzen,
	Papoulis,
	Hamming,
	Kaiser,
	Creighton,
	Tukey
} WindowType;

typedef struct tagLALWindowParams {
	INT4        length;
	WindowType  type;
	REAL4       beta;
} LALWindowParams;

void LALWindow(LALStatus *, REAL4Vector *, LALWindowParams *);
void LALCreateREAL4Window(LALStatus *, REAL4Window **, LALWindowParams *);


#ifdef  __cplusplus
}
#endif

#endif /* _WINDOW_H */
