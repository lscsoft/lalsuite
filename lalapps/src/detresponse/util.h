/*
 * Author: David Chin <dwchin@umich.edu> +1-734-709-9119
 * $Id$
 */
#ifndef _DETRESPONSE_UTIL_H
#define _DETRESPONSE_UTIL_H

/* macro for minimum of two arguments */
#define DETRESPONSE_MIN(a, b) (((a) < (b)) ? (a) : (b))

/* wrap fopen(3) and fclose(3) */
FILE *xfopen(const char *path, const char *mode);
int   xfclose(FILE *stream);

/* wrap strncpy(3) */
char *strlcpy(char *dst, const char *src, size_t len);

/* print detector parameters */
void PrintLALDetector(LALDetector * const detector);

/* print source params */
void print_source(const LALSource * source);

/* print time info */
void print_time_info(const LALTimeIntervalAndNSample * time_info);

REAL8 deg_to_rad(REAL8 degrees);

#endif
