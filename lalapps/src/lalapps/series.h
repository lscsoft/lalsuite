/*
*  Copyright (C) 2007 Duncan Brown, Stephen Fairhurst
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

#ifndef SERIES_H_
#define SERIES_H_

#include <lal/LALDatatypes.h>

#ifdef  __cplusplus
extern "C" {
#endif

typedef enum { Time, Freq, Trans } domain;

#define IS_TIME( domain_ ) ( domain_ == Time )
#define IS_FREQ( domain_ ) ( ( domain_ == Freq ) || ( domain_ == Trans ) )
#define IS_TRANS( domain_ ) ( domain_ == Trans )

struct series
{
  char *name;
  LIGOTimeGPS tbeg;
  LIGOTimeGPS tend;
  domain      dom;
  int         type;
  double      step;
  float       f0;
  const char *unit;
  size_t      size;
  float      *data;
  double     *ddata;
};

double epoch_diff( const LIGOTimeGPS *t2, const LIGOTimeGPS *t1 );
void epoch_add( LIGOTimeGPS *t1, LIGOTimeGPS *t0, double dt );
int write_ilwd( const char *fname, const struct series *ser );
struct FrameH *fr_add_proc_data( struct FrameH *frame, const struct series *ser );

#ifdef  __cplusplus
}
#endif

#endif
