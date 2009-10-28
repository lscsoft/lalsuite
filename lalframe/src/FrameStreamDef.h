/*
*  Copyright (C) 2007 Bernd Machenschalk, Jolien Creighton
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

#ifndef _FRAMESTREAMDEF_H
#define _FRAMESTREAMDEF_H

#include <lal/LALFrameL.h>
#include <lal/LALDatatypes.h>

#include <lal/LALRCSID.h>
NRCSID (FRAMESTREAMDEFH,"$Id$");

/* Useful macros */
#define SECNAN_TO_I8TIME( sec, nan ) \
  ((INT8)1000000000*(INT8)(sec)+(INT8)(nan))
/* Dangerous!!!: */
#define EPOCH_TO_I8TIME( epoch ) \
  SECNAN_TO_I8TIME( (epoch).gpsSeconds, (epoch).gpsNanoSeconds )
#define SET_EPOCH( pepoch, i8time ) \
  do { INT8 t=(i8time); LIGOTimeGPS *pe=(pepoch); \
    pe->gpsSeconds=t/(INT8)1000000000; pe->gpsNanoSeconds=t%(INT8)1000000000; \
  } while( 0 )

typedef struct
tagFrFileInfo
{
  INT4  ind;
  CHAR *url;
  INT4  t0;
  INT4  dt;
}
FrFileInfo;

/* Definition of FrStream */
struct
tagFrStream
{
  FrFileInfo     *filelist;
  UINT4           numfiles;
  UINT4           filenum;
  struct FrFile  *frfile;
  struct FrameH  *frame;
  LIGOTimeGPS     epoch;
  INT4            end;
  INT4            err;
  INT4            gap;
};
#endif /* _FRAMESTREAMDEF_H */
