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

#define MAXFILES 60000         /* Maximum # of files in a directory */
#define MAXFILENAMELENGTH 256   /* Maximum # of characters of a SFT filename */

struct CommandLineArgsTag 
{
  char *directory;
  char *outputfile;
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

int ReadCommandLine(int argc,char *argv[],struct CommandLineArgsTag *CLA);
int ReadSFTDirectory(struct CommandLineArgsTag CLA);
int ComputePSD(struct CommandLineArgsTag CLA);


