/*
 * Author: David Chin <dwchin@umich.edu> +1-734-709-9119
 * $Id$
 */
/* error-handler-wrapped versions of std c library functions */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>

#include "detresponse.h"

FILE *
xfopen(const char *path, const char *mode)
{
  FILE * f = (FILE *) NULL;

  f = fopen(path, mode);

  if (f == (FILE *)NULL)
    {
      fprintf(stderr, "%s: ", path);
      perror("could not open file");
      exit(errno);
    }
  
  return f;
}

int
xfclose(FILE * stream)
{
  if (stream != (FILE *)NULL)
    return fclose(stream);
  else
    return 0;
}

void *
xcalloc(size_t nmemb, size_t size)
{
  void *p;
  p = (void *)calloc(nmemb, size);

  if (p == (void *)NULL)
    {
      perror("xcalloc: out of memory");
      exit(EXIT_FAILURE);
    }

  return p;
}

void *
xmalloc(size_t length)
{
  void *p;
  p = (void *)malloc(length);

  if (p == (void *)NULL)
    {
      perror("xmalloc: out of memory");
      exit(EXIT_FAILURE);
    }

  return p;
}

void *
xrealloc(void *p, size_t length)
{
  p = (void *)realloc(p, length);

  if (p == (void *)NULL)
    {
      perror("xrealloc: out of memory");
      exit(EXIT_FAILURE);
    }

  return p;
}
void
print_source(const LALSource * source)
{
  printf("Source name: %s\n", source->name);
  printf("RA:          % 14.20e rad\n",
         source->equatorialCoords.longitude);
  printf("Dec:         % 14.20e rad\n",
         source->equatorialCoords.latitude);
  printf("Orientation: % 14.20e rad\n", source->orientation);

}

void
PrintLALDetector(LALDetector * const detector)
{
  printf("Detector  = \n");
  switch(detector->type) {
    case LALDETECTORTYPE_IFODIFF:
      printf("IFO, differential mode\n");
      break;
    case LALDETECTORTYPE_IFOXARM:
      printf("IFO, one-arm mode (X-arm)\n");
      break;
    case LALDETECTORTYPE_IFOYARM:
      printf("IFO, one-arm mode (Y-arm\n");
      break;
    case LALDETECTORTYPE_IFOCOMM:
      printf("IFO, common mode\n");
      break;
    case LALDETECTORTYPE_CYLBAR:
      printf("Cylindrical bar\n");
      break;
    case LALDETECTORTYPE_ABSENT:
      printf("Undefined\n");
      break;
    default:
      printf("Undefined\n");
      break;
  } 
  printf("{\n");
  printf("  { % 22.15e, % 22.15e, % 22.15e },\n",
         detector->location[0],
         detector->location[1],
         detector->location[2]);
  printf("  {\n");
  printf("    { % 22.15e, % 22.15e, % 22.15e },\n",
         detector->response[0][0],
         detector->response[0][1],
         detector->response[0][2]);
  printf( "    { % 22.15e, % 22.15e, % 22.15e },\n",
          detector->response[1][0],
          detector->response[1][1],
          detector->response[1][2]);
  printf( "    { % 22.15e, % 22.15e, % 22.15e }\n  },\n",
          detector->response[2][0],
          detector->response[2][1],
          detector->response[2][2]);
  printf("  {\n");
  printf("    \"%s\",\n", detector->frDetector.name);
  printf("    vertex longitude = % 22.15e,\n",
         detector->frDetector.vertexLongitudeRadians);
  printf("    vertex latitude  = % 22.15e,\n",
         detector->frDetector.vertexLatitudeRadians);
  printf("    vertex elevation = % 22.15e,\n",
         detector->frDetector.vertexElevation);
  printf("    X-arm altitude   = % 22.15e,\n",
         detector->frDetector.xArmAltitudeRadians);
  printf("    X-arm azimuth    = % 22.15e,\n",
         detector->frDetector.xArmAzimuthRadians);
  printf("    Y-arm altitude   = % 22.15e,\n",
         detector->frDetector.yArmAltitudeRadians);
  printf("    Y-arm azimuth    = % 22.15e\n  }\n}\n",
         detector->frDetector.yArmAzimuthRadians);
  return;
}

void print_time_info(const LALTimeIntervalAndNSample * p_time_info)
{
  printf("Accuracy: ");
  switch (p_time_info->accuracy) {
    case LALLEAPSEC_STRICT:
      printf("strict\n");
      break;
    case LALLEAPSEC_LOOSE:
      printf("loose\n");
      break;
    default:
      printf("undefined\n");
      break;
  }
  printf("time_info.epoch.gpsSeconds     = % 14d\n",
         (*p_time_info).epoch.gpsSeconds);
  printf("time_info.epoch.gpsNanoSeconds = % 14d\n",
         (*p_time_info).epoch.gpsNanoSeconds);
  printf("time_info.deltaT               = % 14e\n",
         (*p_time_info).deltaT);
  printf("time_info.nSample              = % 14d\n",
         (*p_time_info).nSample);
}


void print_response(const LALDetAMResponse *r)
{
  printf("Plus:  % 20.14e\n", r->plus);
  printf("Cross: % 20.14e\n", r->cross);
}

/*
 * stolen from the doxygen source
 */
int 
mystrncasecmp(char *s1, char *s2, unsigned int n)
{
  if (n == 0)
    return 0;
  
  while ((n-- != 0)
         && (tolower(*(unsigned char *)s1) == tolower(*(unsigned char *)s2)))
    {
      if (n == 0 || *s1 == '\0' || *s2 == '\0')
        return 0;
      s1++;
      s2++;
    }
  
  return tolower(*(unsigned char *) s1) - tolower(*(unsigned char *) s2);
}



/*
 * strlcpy is non-standard; implementation is
 * copied from OpenBSD source.
 */
/*
 * Copyright (c) 1998 Todd C. Miller <Todd.Miller@courtesan.com>
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. The name of the author may not be used to endorse or promote products
 *    derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES,
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
 * AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL
 * THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
 * OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
 * OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
/*
 * Copy src to string dst of size siz.  At most siz-1 characters
 * will be copied.  Always NUL terminates (unless siz == 0).
 * Returns strlen(src); if retval >= siz, truncation occurred.
 */
size_t mystrlcpy(char *dst, const char *src, size_t siz)
{
        char *d = dst;
        const char *s = src;
        size_t n = siz;

        /* Copy as many bytes as will fit */
        if (n != 0 && --n != 0) {
                do {
                        if ((*d++ = *s++) == 0)
                                break;
                } while (--n != 0);
        }

        /* Not enough room in dst, add NUL and traverse rest of src */
        if (n == 0) {
                if (siz != 0)
                        *d = '\0';              /* NUL-terminate dst */
                while (*s++)
                        ;
        }

        return(s - src - 1);    /* count does not include NUL */
}



/*
 * Copyright (c) 1998 Todd C. Miller <Todd.Miller@courtesan.com>
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. The name of the author may not be used to endorse or promote products
 *    derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES,
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
 * AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL
 * THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
 * OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
 * OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/*
 * Appends src to string dst of size siz (unlike strncat, siz is the
 * full size of dst, not space left).  At most siz-1 characters
 * will be copied.  Always NUL terminates (unless siz <= strlen(dst)).
 * Returns strlen(src) + MIN(siz, strlen(initial dst)).
 * If retval >= siz, truncation occurred.
 */
size_t
mystrlcat(char *dst, const char *src, size_t siz)
{
        char *d = dst;
        const char *s = src;
        size_t n = siz;
        size_t dlen;

        /* Find the end of dst and adjust bytes left but don't go past end */
        while (n-- != 0 && *d != '\0')
                d++;
        dlen = d - dst;
        n = siz - dlen;

        if (n == 0)
                return(dlen + strlen(s));
        while (*s != '\0') {
                if (n != 1) {
                        *d++ = *s;
                        n--;
                }
                s++;
        }
        *d = '\0';

        return(dlen + (s - src));       /* count does not include NUL */
}



void
square_timeseries(REAL4TimeSeries * in)
{
  size_t i;

  for (i = 0; i < in->data->length; ++i)
    {
      in->data->data[i] *= in->data->data[i];
    }
}



void
add_timeseries(REAL4TimeSeries * sum, REAL4TimeSeries * a,
               REAL4TimeSeries * b)
{
  size_t i;

  /* assume that all the info about the timeseries a & b are the same */
  (void)mystrlcpy(sum->name, a->name, LALNameLength);
  sum->epoch = a->epoch;
  sum->deltaT = a->deltaT;
  sum->f0 = a->f0;
  sum->sampleUnits = a->sampleUnits;

  for (i = 0; i < sum->data->length; ++i)
    {
      sum->data->data[i] = a->data->data[i] + b->data->data[i];
    }
}



void set_detector_params(LALStatus * status,
                         LALFrDetector * frdet, LALDetector * det,
                         const char * name,
                         REAL8 vertex_longitude,
                         REAL8 vertex_latitude,
                         REAL4 vertex_elevation,
                         REAL4 x_altitude,
                         REAL4 x_azimuth,
                         REAL4 y_altitude,
                         REAL4 y_azimuth)
{
  (void)mystrlcpy(frdet->name, name, LALNameLength);
  frdet->vertexLongitudeRadians = vertex_longitude;
  frdet->vertexLatitudeRadians  = vertex_latitude;
  frdet->vertexElevation        = vertex_elevation;
  frdet->xArmAltitudeRadians    = x_altitude;
  frdet->xArmAzimuthRadians     = x_azimuth;
  frdet->yArmAltitudeRadians    = y_altitude;
  frdet->yArmAzimuthRadians     = y_azimuth;

  LALCreateDetector(status, det, frdet, LALDETECTORTYPE_IFODIFF);
}


  
REAL8
deg_to_rad(REAL8 degrees)
{
  return degrees * (REAL8)(LAL_PI / 180.);
}

