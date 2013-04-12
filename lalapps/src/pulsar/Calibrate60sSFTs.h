/*
*  Copyright (C) 2007 Xavier Siemens
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

#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <glob.h>
#include <lal/LALDatatypes.h>
#include <lal/LALBarycenter.h>
#include <lal/LALInitBarycenter.h>
#include <lal/LALDemod.h>

#define MAXFILES 70000          /* Maximum # of files in a directory */
#define MAXFILENAMELENGTH 256   /* Maximum # of characters of a SFT filename*/

#define MAXLINESRS   500000     /* Maximum # of lines in a Response or Sensing file */
#define MAXLINESF    100000     /* Maximum # of lines in a Factors file */

struct CommandLineArgsTag 
{
  char *directory;
  char *caldirectory;
  char *run;
  char *IFO; 
  char *outputdirectory;
} CommandLineArgs;

struct headertag 
{
  REAL8 endian;
  INT4  gps_sec;
  INT4  gps_nsec;
  REAL8 tbase;
  INT4  firstfreqindex;
  INT4  nsamples;
} header;

typedef struct ResponseFunctionTag
{
  REAL4 Frequency[MAXLINESRS];
  REAL4 Magnitude[MAXLINESRS];
  REAL4 Phase[MAXLINESRS];
  REAL4 re[MAXLINESRS];
  REAL4 im[MAXLINESRS];
} Response;

typedef struct SensingFunctionTag
{
  REAL4 Frequency[MAXLINESRS];
  REAL4 Magnitude[MAXLINESRS];
  REAL4 Phase[MAXLINESRS];
  REAL4 re[MAXLINESRS];
  REAL4 im[MAXLINESRS];
} Sensing;

struct FactorsTag
{
  INT4  t[MAXLINESF];      /* in GPS seconds */
  REAL4 alpha[MAXLINESF];
  REAL4 alpha_beta[MAXLINESF];
} Factors;

int ReadCommandLine(int argc,char *argv[],struct CommandLineArgsTag *CLA);
int ReadSFTDirectory(struct CommandLineArgsTag CLA);
int ReadCalibrationFiles(struct CommandLineArgsTag CLA);
int ComputeInitialRSFunctions(void);
int CalibrateSfts(struct CommandLineArgsTag CLA);
int Freemem(void);


COMPLEX8 tmpa, tmpb, tmpc; 
REAL4 tmpx, tmpy;

#define cmul( a, b ) \
( tmpa = (a), tmpb = (b), \
  tmpc.realf_FIXME = crealf(tmpa) * crealf(tmpb) - cimagf(tmpa) * cimagf(tmpb), \
  tmpc.imagf_FIXME = crealf(tmpa) * cimagf(tmpb) + cimagf(tmpa) * crealf(tmpb), \
  tmpc )

#define cdiv( a, b ) \
( tmpa = (a), tmpb = (b), \
  fabs( crealf(tmpb) ) >= fabs( cimagf(tmpb) ) ? \
    ( tmpx = cimagf(tmpb) / crealf(tmpb), \
      tmpy = crealf(tmpb) + tmpx * cimagf(tmpb), \
      tmpc.realf_FIXME = ( crealf(tmpa) + tmpx * cimagf(tmpa) ) / tmpy, \
      tmpc.imagf_FIXME = ( cimagf(tmpa) - tmpx * crealf(tmpa) ) / tmpy, \
      tmpc ) : \
    ( tmpx = crealf(tmpb) / cimagf(tmpb), \
      tmpy = cimagf(tmpb) + tmpx * crealf(tmpb), \
      tmpc.realf_FIXME = ( crealf(tmpa) * tmpx + cimagf(tmpa) ) / tmpy, \
      tmpc.imagf_FIXME = ( cimagf(tmpa) * tmpx - crealf(tmpa) ) / tmpy, \
      tmpc ) )

