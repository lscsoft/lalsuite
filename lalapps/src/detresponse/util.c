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

char *
mystrlcpy(char       * dst,
          const char * src,
          size_t       len)
{
  char * retval = strncpy(dst, src, len);
  dst[DETRESPONSE_MIN(len, strlen(src)+1) - 1] = '\0';
  return retval;
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
  printf( "Detector  = \n");
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
  printf("time_info.epoch.gpsSeconds     = % 9d\n",
         (*p_time_info).epoch.gpsSeconds);
  printf("time_info.epoch.gpsNanoSeconds = % 9d\n",
         (*p_time_info).epoch.gpsNanoSeconds);
  printf("time_info.deltaT               = % 14.20e\n",
         (*p_time_info).deltaT);
  printf("time_info.nSample              = % 9d\n",
         (*p_time_info).nSample);
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


  
REAL8
deg_to_rad(REAL8 degrees)
{
  return degrees * (REAL8)(LAL_PI / 180.);
}

