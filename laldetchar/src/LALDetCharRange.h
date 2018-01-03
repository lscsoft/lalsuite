/*
 *  Copyright (C) 2012 Duncan M. Macleod
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

#ifndef _LALDETCHARRANGE_H
#define _LALDETCHARRANGE_H

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

#include <lal/LALDatatypes.h>

REAL4
XLALREAL4InspiralRange(
    const REAL4FrequencySeries * spectrum,
    REAL4 snr,
    REAL4 mass_1,
    REAL4 mass_2,
    REAL4 f_min
);

REAL8
XLALREAL8InspiralRange(
    const REAL8FrequencySeries * spectrum,
    REAL8 snr,
    REAL8 mass_1,
    REAL8 mass_2,
    REAL8 f_min
);

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _LALDETCHARRANGE_H */
