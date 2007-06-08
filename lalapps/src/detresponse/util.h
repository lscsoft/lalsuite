/*
*  Copyright (C) 2007 David Chin, Jolien Creighton
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

/*
 * Author: David Chin <dwchin@umich.edu> +1-734-709-9119
 * $Id$
 */
#ifndef _DETRESPONSE_UTIL_H
#define _DETRESPONSE_UTIL_H

/* macro for minimum of two arguments */
#define DETRESPONSE_MIN(a, b) (((a) < (b)) ? (a) : (b))

/* wrap fopen(3) and allocation functions a la GNU */
FILE *xfopen(const char *path, const char *mode);
int xfclose(FILE * stream);
void *xmalloc(size_t length);
void *xrealloc(void *p, size_t length);
void *xcalloc(size_t nmemb, size_t length);

/* print detector parameters */
void PrintLALDetector(LALDetector * const detector);

/* print source params */
void print_source(const LALSource * source);

/* print time info */
void print_time_info(const LALTimeIntervalAndNSample * time_info);

/* print detector response */
void print_response(const LALDetAMResponse *resp);

int mystrncasecmp(char *s1, char *s2, unsigned int n);

/* strlcpy is non-standard, so emulate it here */
size_t mystrlcpy(char *dst, const char *src, size_t size);
/* strlcat is non-standard, so emulate it here */
size_t mystrlcat(char *dst, const char *src, size_t size);

void square_timeseries(REAL4TimeSeries *ts);
void add_timeseries(REAL4TimeSeries * sum, REAL4TimeSeries * a,
                    REAL4TimeSeries * b);

void set_detector_params(LALStatus * status,
                         LALFrDetector * frdet, LALDetector * det,
                         const char * name,
                         REAL8 vertex_longitude,
                         REAL8 vertex_latitude,
                         REAL4 vertex_elevation,
                         REAL4 x_altitude,
                         REAL4 x_azimuth,
                         REAL4 y_altitude,
                         REAL4 y_azimuth);

REAL8 deg_to_rad(REAL8 degrees);

#endif
