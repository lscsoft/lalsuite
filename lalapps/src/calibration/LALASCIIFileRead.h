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

#include <stdio.h>
#include <lal/LALDatatypes.h>

typedef struct tagLALDataFileNameFields
{
  CHAR  site[FILENAME_MAX];
  CHAR  description[FILENAME_MAX];
  INT4  tstart;
  INT4  duration;
  CHAR  extension[FILENAME_MAX];
}
LALDataFileNameFields;

typedef struct tagLALCalRefFileNameDescriptionFields
{
  CHAR ifo[FILENAME_MAX];
  CHAR channelPostfix[FILENAME_MAX];
  CHAR run[FILENAME_MAX];
  INT4 version;
}
LALCalRefFileNameDescriptionFields;

typedef struct tagLALCalFacFileNameDescriptionFields
{
  CHAR ifo[FILENAME_MAX];
  CHAR run[FILENAME_MAX];
  INT4 version;
  INT4 deltaT;
}
LALCalFacFileNameDescriptionFields;

int XLALDataFileNameParse( LALDataFileNameFields *fields, const char *fname );

int XLALCalRefFileNameDescriptionParse( LALCalRefFileNameDescriptionFields *fields, const char *description );
int XLALCalFacFileNameDescriptionParse( LALCalFacFileNameDescriptionFields *fields, const char *description );

int XLALASCIIFileCountRows( const char *fname );
REAL8VectorSequence * XLALASCIIFileReadColumns( INT4 ncol, const char *fname );

REAL4 XLALASCIIFileReadCalFacHeader( const char *fname );
REAL4 XLALASCIIFileReadCalRefHeader( const char *fname );
int XLALASCIIFileReadCalFac( REAL4TimeSeries **alpha, REAL4TimeSeries **gamma, const char *fname );
int XLALASCIIFileReadCalRef( COMPLEX8FrequencySeries **series, REAL8 *duration, const char *fname );
