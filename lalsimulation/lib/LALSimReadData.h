/*
*  Copyright (C) 2013 Jolien Creighton
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

#ifndef _LALSIMREADDATA_H
#define _LALSIMREADDATA_H

#include <stddef.h>
#include <lal/FileIO.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
}       /* so that editors will match preceding brace */
#endif

LALFILE *XLALSimReadDataFileOpen(const char *fname);
size_t XLALSimReadDataFile2Col(double **xdat, double **ydat, LALFILE * fp);
size_t XLALSimReadDataFileNCol(double **data, size_t *ncol, LALFILE * fp);

#if 0
{       /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _LALSIMREADDATA_H */
