/*
*  Copyright (C) 2008 Matt Pitkin, Curt Cutler, E. Myles Standish, David
*  Hoffman.
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


/* Matt Pitkin - 07/04/08 */

/**
 * \file
 * \ingroup pulsarApps
 * \author Matt Pitkin
 * \brief
 * This code will take in JPL binary ephemeris files containing (e.g.
 * DE200.1950.2050), and output the position and velocities of a given solar
 * system body over a given time range. The code is based on Curt Cutler's
 * FORTRAN code, which itself uses E. Miles. Standish's JPL functions. However,
 * large parts of the ephemeris file reading function and interpolation
 * functions have been taken from David Hoffman's ephem_read.c code found at
 * http://www.projectpluto.com/jpl_eph.htm - the code found there offers far
 * more functionality then my code.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <getopt.h>

#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/LALVCSInfo.h>
#include <LALAppsVCSInfo.h>

#define MERCURY 0
#define VENUS 1
#define EARTH 2                                    /* Earth-Moon barycenter */
#define MARS 3
#define JUPITER 4
#define SATURN 5
#define URANUS 6
#define NEPTUNE 7
#define PLUTO 8
#define MOON 9                                      /* Relative to geocenter */
#define SUN 10

#define CURT_AU 1.4959787066e8  /* 1 AU defined from Curt's code (the same as
                                 * the AU definition in the JPL DE200 file) */ 

#define USAGE \
"Usage: %s [options]\n\n"\
" --help              display this message\n"\
" --verbose           display all error messages\n"\
" --ephem-file        path to and name of JPL binary ephemeris file\n\
                     e.g. /home/matthew/tempo_ephem/DE405.1950.2050\n"\
" --output-file        path to and name of file to output the ephemeris\n"\
" --year              year over which to calculate the ephemeris e.g. 2008 \n"\
" --interval          time step between successive output ephemeris points\n\
                     (in integer hours)\n"\
" --num-years         number of years over which the ephemeris will be\n\
                     created (integer)\n"\
" --overlap           number of days overlap with previous and next year\n"\
" --target            the target solar system body (e.g. SUN, EARTH)\n"\
"\n"\
" --test              compare output with current ephemeris files\n"\
"\n"

#define TESTFAIL 11

/* struct for the binary header records (adpated from
http://www.projectpluto.com/jpl_eph.htm file ephem_types.h) */
typedef struct tagheaderData1{
  CHAR ttl[3][84];
  CHAR cnam[400][6];
  REAL8 ss[3];
  INT4 con;
  REAL8 au;
  REAL8 emrat;
  INT4 ipt[12][3];
  INT4 numde;
  INT4 libratPtr[3]; /* for libration */
} headerData1;

typedef struct tagheaderRecord1{
  headerData1 data;

  CHAR pad[10000]; /* padding for overflow data */
} headerRecord1;

typedef struct tagheaderData2{
  REAL8 cval[400];
} headerData2;

typedef struct tagheaderRecord2{
  headerData2 data;

  CHAR pad[10000];
} headerRecord2;

typedef struct taginputParams_t{
  CHAR ephemfile[256]; /* path and name of JPL binary ephemeris file */
  CHAR outputfile[256]; /* path and name of output ephemeris file */

  INT4 year;  /* year of ephemeris file to extract */
  REAL8 nhre;  /* number of hours between successive data point for output */
  INT4 noverlap; /* number of days overlap with previos years */
  INT4 nyears; /* number of years over which to create ephemeris */

  INT4 target; /* the solar system body */
  CHAR *targName;
} inputParams_t;

/* header and data records from ephemeris file */
headerRecord1 head1;
headerRecord2 head2;

REAL8 Tbeg=0., Tend=0., Tspan=0.;

/* global variables */
INT4 verbose=0;
INT4 test=0;

UINT4 ARRAY_SIZE=0;

/* function to convert a GPS time in JD to TDB (similar to my routine
LALTDBMJDtoGPS */
void convert(REAL8 *gps_JD, REAL8 *time);

/* function to read the record length of the binary file - this differs for
   different ephemeris files (equivalent to Standish's FSIZER2 */
INT4 fsizer(FILE *fp);

/* functions to swap endianness - taken from SFTfileIO.c file */
void endian_swap(CHAR *pdata, size_t dsize, size_t nelements);

/* function to read in coefficients from the ephemeris file for a give time */
void read_coeffs(REAL8 *coeffArray, REAL8 *time, FILE *fp);

void interpolate_state(REAL8 *coeffArray, REAL8 *time, INT4 target,
  REAL8 *state, FILE *fp);

void pleph(REAL8 *coeffArray, REAL8 *time, INT4 target, REAL8 *state, FILE *fp);

/* function to get the input arguments from the command line */
void get_input_args(inputParams_t *inputParams, INT4 argc, CHAR *argv[]);

int main(int argc, char **argv){
  FILE *fpe=NULL;
  FILE *fp=NULL; /* file pointer for ephemeris file */

  /* variable names taken from Curt's code */
  REAL8 jd_2000=2451544.5;
  INT4 gps=0, gps_2000=630720000, gps_yr=0;
  REAL8 fgps=0.; /* float version of gps */
  REAL8 time[2], R[6], A[3], gps_JD[2];
  REAL8 Vlast[3], Vnow[3], Vnext[3], Rnow[3];
  INT4 nyr=0, ndays=0;
  REAL8 nhr=0.;

  REAL8 finterval=0, halfinterval_jd=0;

  REAL8 *coeffArray=NULL;

  REAL8 day=86400., hour=3600.;
  REAL8 halfhour_jd=1800./86400.;
  REAL8 gps_JD_start1=0.;

  INT4 leapdays=0, extraleaps=0;

  INT4 nentries=0;

  INT4 i=0;

  CHAR *lalpath=NULL, earthEphem[256];
  CHAR *tempopath=NULL;
  REAL8 pos[3], vel[3], acc[3], title[3];

  inputParams_t inputs;

  /* allocate memory for headers structures */
  /* head1 = calloc(sizeof(headerRecord1), 1);
  head2 = calloc(sizeof(headerRecord2), 1); */

  /* get input parameters */
  get_input_args(&inputs, argc, argv);

  if(test){
    if(verbose)
      fprintf(stderr, "Performing test against the existing ephemeris.\n");

    /* get the path to LALPULSAR from the environment variables */
    if((lalpath = getenv("LALPULSAR_PREFIX")) == NULL){
      fprintf(stderr, "LALPULSAR_PREFIX environment variable not set. Cannot perform test!\n");
      return TESTFAIL;
    }

    /* set the name of one of the earth ephemeris files in LAL */
    sprintf(earthEphem, "%s/share/lalpulsar/earth98.dat", lalpath);

    /* get JPL ephemeris files - only try if TEMPO2 is installed and the TEMPO2
       environment variable is set */
    if((tempopath = getenv("TEMPO2")) == NULL){
      fprintf(stderr, "You must have the TEMPO2 environment variable set to\
 find the JPL ephemeris.\n");
      return TESTFAIL;
    }
    sprintf(inputs.ephemfile, "%s/ephemeris/DE405.1950.2050", tempopath);

    /* set other inputs */
    CHAR tmp[] = "EARTH";
    inputs.targName = tmp;
    inputs.target = EARTH;
    inputs.year = 1998;
    inputs.nhre = 4;
    inputs.nyears = 1;
    inputs.noverlap = 10;

    if((fpe = fopen(earthEphem, "r")) == NULL){
      fprintf(stderr, "Error... couldn't open file %s for reading!\n",
      earthEphem);
      return TESTFAIL;
    }
  }

  if(verbose){
    fprintf(stderr, "Input information:\n");
    fprintf(stderr, "  Reading from JPL ephmeris file:- %s\n",
inputs.ephemfile);
    fprintf(stderr, "  Targeting solar system body:- %s\n", inputs.targName);
    fprintf(stderr, "  Ephemeris year:- %d\n", inputs.year);
    fprintf(stderr, "  Number of years:- %d\n", inputs.nyears);
    fprintf(stderr, "  Output interval (hours):- %.2lf\n", inputs.nhre);
    fprintf(stderr, "  Outout overlap with adjacent years (days):- %d\n",
inputs.noverlap);
  }

  nyr = inputs.year;

  /* calculate the number of leap days between Jan on the given year and Jan in
     2000 */
  if(nyr >= 2001)
    leapdays = (INT4)(nyr-1997)/4;
  else
    leapdays = -(INT4)(2000-nyr)/4;

  /* work out how many leap days there will be over the period of the ephemeris
     - nyears */
  if(nyr + inputs.nyears >= 2001)
    extraleaps = (INT4)(nyr + inputs.nyears - 1997)/4;
  else
    extraleaps = -(INT4)(2000 - (nyr + inputs.nyears))/4;

  extraleaps -= leapdays;
  ndays = 365*inputs.nyears + extraleaps;

  gps_yr = gps_2000 + day*(365*(nyr-2000) + leapdays);
  // unused: INT4 gps_start = gps_yr - inputs.noverlap*day;
  gps_JD_start1 = jd_2000 + 365.*(REAL8)(nyr-2000.) + (REAL8)leapdays -
    (REAL8)inputs.noverlap;

  /* calculate the number of entries in each ephemeris */
  /* nentries = (INT4)((8760 + inputs.noverlap*48)/inputs.nhre) + 1; */
  nentries = (INT4)((ndays*24 + inputs.noverlap*48)/inputs.nhre) + 1;

  /* open ephemeris file for reading */
  if((fp = fopen(inputs.ephemfile, "r")) == NULL){
    fprintf(stderr, "Error... couldn't open file %s for reading!\n",
      inputs.ephemfile);
    exit(1);
  }

  /* read the header info from the ephemeris files an work out the array size */
  if( (ARRAY_SIZE = fsizer(fp)) == 0 ){
    fprintf(stderr, "Array size is zero!\n");
    return 1;
  }

  /* allocate memory for chebychev polynomial coefficients */
  if( (coeffArray = realloc(coeffArray, ARRAY_SIZE*sizeof(REAL8))) == NULL ){
    fprintf(stderr, "Error allocating memory.\n");
    return 1;
  }

  if(verbose){
    fprintf(stderr, "The array size for ephemeris file %s is %d.\n",
      inputs.ephemfile , ARRAY_SIZE);
  }

  /* get first set of coefficients from ephemeris file */
  /*time=0.; */
  time[0] = 0.;
  time[1] = 0.;
  read_coeffs(coeffArray, time, fp);

  if(test){
    if( fscanf(fpe, "%lf%lf%lf", &title[0], &title[1], &title[2]) != 3 ){
      fprintf(stderr, "Error reading in ephemeris header.\n");
      return TESTFAIL;
    }

    /* check that values are the same as in the existing ephemeris file */
    if(gps_yr - title[0] != 0. || inputs.nhre*hour - title[1] != 0. ||
       nentries - title[2] != 0.){
      fprintf(stderr, "Start time, or length of ephemeris is not the same as \
in the existing file!\n");
      return TESTFAIL;
    }
  }
  else{
    /* open the Earth and Sun files for writing */
    if((fpe = fopen(inputs.outputfile, "w")) == NULL){
      fprintf(stderr, "Error... can't open output ephemeris file for \
writing!\n");
      exit(1);
    }

    /* output header information on lines starting with a # comment */
    fprintf(fpe, "# Build information for %s\n", argv[0]);
    fprintf(fpe, "# Author: "LALAPPS_VCS_AUTHOR"\n");
    fprintf(fpe, "# LALApps Commit ID: "LALAPPS_VCS_ID"\n");
    fprintf(fpe, "# LALApps Commit Date: "LALAPPS_VCS_DATE"\n");
    fprintf(fpe, "#\n# Ephemeris creation command:-\n#\t");
    for( INT4 k=0; k<argc; k++ ) fprintf(fpe, "%s ", argv[k]);
    fprintf(fpe, "\n");
    
    CHAR *efile = strrchr(inputs.ephemfile, '/');
    
    fprintf(fpe, "#\n# JPL ephemeris file %s from TEMPO2 \
(http://www.atnf.csiro.au/research/pulsar/tempo2/)\n", efile+1);

    /* some information about the data */
    fprintf(fpe, "#\n# This file consists of a header line containing:\n");
    fprintf(fpe, "#\tGPS time of the start year, interval between entries \
(secs), no. of entries\n");
    fprintf(fpe, "# Each entry consists of:\n");
    fprintf(fpe, "#\tGPS time\t\tPos. x (lt sec)\t\tPos. y (lt sec)\n\
#\tPos. z (lt sec)\t\tVel. x (lt sec/sec)\tVel. y (lt sec/sec)\n\
#\tVel. z (lt sec/sec)\tAcc. x (lt sec/sec^2)\tAcc. y (lt sec/sec^2)\n\
#\tAcc. z (lt sec/sec^2)\n");
    fprintf(fpe, "# with entries calculated at an approximation to \
\"ephemeris time\" via the conversion:\n\
#    dT = (GPS(JD) + 51.148/86400) - 2451545\n\
#    Teph = GPS(JD) + (51.148 + 0.001658*sin(6.24008 + 0.017202*dT)\n\
#           + 1.4e-5*sin(2*(6.24008 + 0.017202*dT)))/86400.\n\
# where GPS(JD) is the GPS time converted into Julian Days, and with 1 A.U. \
defined as %.3f km\n", CURT_AU);
    

    fprintf(fpe, "\t%d\t%lf\t%d\n", gps_yr, inputs.nhre*hour, nentries);
  }

  nhr = inputs.nhre;

  finterval = hour*nhr;
  halfinterval_jd=halfhour_jd*nhr;

  gps_JD[0] = gps_JD_start1;
  gps_JD[1] = 0.;

  gps = gps_2000 + (gps_JD[0] - jd_2000)*day;

  convert(gps_JD, time);

  pleph(coeffArray, time, inputs.target, R, fp);

  Rnow[0] = R[0];
  Rnow[1] = R[1];
  Rnow[2] = R[2];
  Vnow[0] = R[3];
  Vnow[1] = R[4];
  Vnow[2] = R[5];

  gps_JD[0] = gps_JD_start1 - 1.;
  gps_JD[1] = 1. - halfinterval_jd;

  convert(gps_JD, time);

  /* get values of velocity around present time to calculate the accelerations
  */
  pleph(coeffArray, time, inputs.target, R, fp);
  Vlast[0] = R[3];
  Vlast[1] = R[4];
  Vlast[2] = R[5];

  gps_JD[0] = gps_JD_start1;
  gps_JD[1] = halfinterval_jd;

  convert(gps_JD, time);

  pleph(coeffArray, time, inputs.target, R, fp);
  Vnext[0] = R[3];
  Vnext[1] = R[4];
  Vnext[2] = R[5];

  /* calculate accelerations */
  A[0] = (Vnext[0] - Vlast[0])/finterval;
  A[1] = (Vnext[1] - Vlast[1])/finterval;
  A[2] = (Vnext[2] - Vlast[2])/finterval;

  fgps = (REAL8)gps;

  if(test){
    /* compare the ephemeris positions, velocites and accelerations */
    if( fscanf(fpe, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", &fgps, &pos[0], &pos[1],
      &pos[2], &vel[0], &vel[1], &vel[2], &acc[0], &acc[1], &acc[2]) != 10){
      fprintf(stderr, "Error reading in ephemeris values.\n");
      return TESTFAIL;
    }

    if(fabs(Rnow[0]-pos[0]) > 1e-12 || fabs(Rnow[1]-pos[1]) > 1e-12 ||
      fabs(Rnow[2]-pos[2]) > 1e-12 || fabs(Vnow[0]-vel[0]) > 1e-19 ||
      fabs(Vnow[1]-vel[1]) > 1e-19 || fabs(Vnow[2]-vel[2]) > 1e-19 ||
      fabs(A[0]-acc[0]) > 1e-23 || fabs(A[1]-acc[1]) > 1e-23 ||
      fabs(A[2]-acc[2]) > 1e-23){
      fprintf(stderr, "Position, velocity or acceleration is not the same as \
in the exiting ephmeris file.\n");
        return TESTFAIL;
    }
  }
  else{
    fprintf(fpe,
"\t%.7lf\t%.15lf\t%.15lf\n\t%.15lf\t%.16le\t%.16le\n\t%.16le\t%.16le\t%.16le\n\
\t%.16le\n", fgps, Rnow[0], Rnow[1], Rnow[2], Vnow[0], Vnow[1], Vnow[2],
A[0], A[1], A[2]);
  }

  /* reset Vnext to Vlast */
  Vlast[0] = Vnext[0];
  Vlast[1] = Vnext[1];
  Vlast[2] = Vnext[2];

  for(i=1;i< nentries;i++){
    gps_JD[1] += halfinterval_jd;

    if(gps_JD[1] >= 1.){
      gps_JD[1] -= 1.;
      gps_JD[0] += 1.;

      /* eliminate any residual round-off error */
      if(fabs(gps_JD[1]) < 1.e-10) gps_JD[1] = 0.;
    }

    gps = gps_2000 + (gps_JD[0] - jd_2000)*day + (INT4)(gps_JD[1]*day + 1.e-4);
/* the int call and 1.d-4 are just to make sure gps is the ``right'' integer*/
    convert(gps_JD, time);
    
    pleph(coeffArray, time, inputs.target, R, fp);
    Rnow[0] = R[0];
    Rnow[1] = R[1];
    Rnow[2] = R[2];
    Vnow[0] = R[3];
    Vnow[1] = R[4];
    Vnow[2] = R[5];
    
    gps_JD[1] += halfinterval_jd;
    if(gps_JD[1] >= 1.){
      gps_JD[1] -= 1.;
      gps_JD[0] += 1.;
    }

    convert(gps_JD, time);

    pleph(coeffArray, time, inputs.target, R, fp);
    Vnext[0] = R[3];
    Vnext[1] = R[4];
    Vnext[2] = R[5];
    A[0] = (Vnext[0] - Vlast[0])/finterval;
    A[1] = (Vnext[1] - Vlast[1])/finterval;
    A[2] = (Vnext[2] - Vlast[2])/finterval;

    fgps = (REAL8)gps;

    if(test){
      /* compare the ephemeris positions, velocites and accelerations */
      if( fscanf(fpe, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", &fgps, &pos[0], &pos[1],
        &pos[2], &vel[0], &vel[1], &vel[2], &acc[0], &acc[1], &acc[2]) != 10){
        fprintf(stderr, "Error reading in ephemeris values.\n");
        return TESTFAIL;
      }

      if(fabs(Rnow[0]-pos[0]) > 1e-12 || fabs(Rnow[1]-pos[1]) > 1e-12 ||
         fabs(Rnow[2]-pos[2]) > 1e-12 || fabs(Vnow[0]-vel[0]) > 1e-19 ||
         fabs(Vnow[1]-vel[1]) > 1e-19 || fabs(Vnow[2]-vel[2]) > 1e-19 ||
         fabs(A[0]-acc[0]) > 1e-23 || fabs(A[1]-acc[1]) > 1e-23 ||
         fabs(A[2]-acc[2]) > 1e-23){
        fprintf(stderr, "Position, velocity or acceleration is not the same as \
in the exiting ephmeris file.\n");
        return TESTFAIL;
      }
    }
    else{
      fprintf(fpe,
"\t%.7lf\t%.15lf\t%.15lf\n\t%.15lf\t%.16le\t%.16le\n\t%.16le\t%.16le\t%.16le\n\
\t%.16le\n", fgps, Rnow[0], Rnow[1], Rnow[2], Vnow[0], Vnow[1], Vnow[2], A[0],
A[1], A[2]);
    }

    /* reset Vnext to Vlast */
    Vlast[0] = Vnext[0];
    Vlast[1] = Vnext[1];
    Vlast[2] = Vnext[2];
  }

  /* free memory */
  fclose(fp);
  fclose(fpe);

  free(coeffArray);

  if(test && verbose)
    fprintf(stderr, "The code has passed the test. Hoorah!\n");

  return 0;
}

/* function to determine the record size and initialise the ephemeris */
INT4 fsizer(FILE *fp){
  INT4 kmx=0, khi=0;
  UINT4 i=0;
  INT4 nd=0;
  UINT4 size=0;

  /* read in header info */

  /* read in data structure one part at a time so 32vs64bit alignment
     issues don't become a problem */
  if( fread(&head1.data.ttl, sizeof(head1.data.ttl), 1, fp) != 1 ){
    fprintf(stderr, "Error reading in ephemeris ttl header info.\n");
    return 0;
  }
  if( fread(&head1.data.cnam, sizeof(head1.data.cnam), 1, fp) != 1 ){
    fprintf(stderr, "Error reading in ephemeris cnam header info.\n");
    return 0;
  }
  if( fread(&head1.data.ss, sizeof(head1.data.ss), 1, fp) != 1 ){
    fprintf(stderr, "Error reading in ephemeris ss header info.\n");
    return 0;
  }
  if( fread(&head1.data.con, sizeof(head1.data.con), 1, fp) != 1 ){
    fprintf(stderr, "Error reading in ephemeris con header info.\n");
    return 0;
  }
  if( fread(&head1.data.au, sizeof(head1.data.au), 1, fp) != 1 ){
    fprintf(stderr, "Error reading in ephemeris au header info.\n");
    return 0;
  }
  if( fread(&head1.data.emrat, sizeof(head1.data.emrat), 1, fp) != 1 ){
    fprintf(stderr, "Error reading in ephemeris emrat header info.\n");
    return 0;
  }
  if( fread(&head1.data.ipt, sizeof(head1.data.ipt), 1, fp) != 1 ){
    fprintf(stderr, "Error reading in ephemeris ipt header info.\n");
    return 0;
  }
  if( fread(&head1.data.numde, sizeof(head1.data.numde), 1, fp) != 1 ){
    fprintf(stderr, "Error reading in ephemeris numde header info.\n");
    return 0;
  }
  if( fread(&head1.data.libratPtr, sizeof(head1.data.libratPtr), 1, fp) != 1 ){
    fprintf(stderr, "Error reading in ephemeris libratPtr header info.\n");
    return 0;
  }

  /* flip bytes of ipt */
  endian_swap((CHAR*)&head1.data.ipt, sizeof(INT4), 36);
  endian_swap((CHAR*)&head1.data.libratPtr, sizeof(INT4), 3);
  
  /*** Calculate array size in the ephemeris */
  for (i=0;i<12;i++){
    if(head1.data.ipt[i][0] > kmx){
      kmx = head1.data.ipt[i][0];
      khi = i;
    }
  }

  if(head1.data.libratPtr[0] > kmx){
    kmx = head1.data.libratPtr[0];
    khi = 12;
  }

  nd = 3;
  if (khi == 11)
    nd = 2;

  if(khi < 12){
    size = (head1.data.ipt[khi][0] +
      nd*head1.data.ipt[khi][1]*head1.data.ipt[khi][2]-1);
  }
  else if(khi == 12){
    size = (head1.data.libratPtr[0] +
      nd*head1.data.libratPtr[1]*head1.data.libratPtr[2]-1);
  }
  /*****************************************/

  /* set file pointer to position of header data 2 */
  fseek(fp, size*sizeof(REAL8), SEEK_SET);

  /* read in values of head2 */
  if( fread(&head2.data.cval, sizeof(head2.data.cval), 1, fp) != 1 ){
    fprintf(stderr, "Error reading in ephemeris cval header info.\n");
    return 0;
  }

  /* Initialise to ephemeris to the point at which the coefficient values start
  */
  fseek(fp, 2*size*sizeof(REAL8), SEEK_SET);
  
  /* flip bytes of values */
  endian_swap((CHAR*)&head1.data.au, sizeof(REAL8), 1);
  endian_swap((CHAR*)&head1.data.emrat, sizeof(REAL8), 1);
  endian_swap((CHAR*)&head1.data.ss, sizeof(REAL8), 3);
  endian_swap((CHAR*)&head1.data.con, sizeof(INT4), 1);
  endian_swap((CHAR*)&head1.data.numde, sizeof(INT4), 1);
  endian_swap((CHAR*)&head2.data.cval, sizeof(REAL8), 400);

  if(verbose){
    fprintf(stderr, "Check value of AU is correct:\n");
    if(head1.data.au > 149597870. && head1.data.au < 149597871.)
      fprintf(stderr, "  Correct: 1 AU = %.16lf km\n", head1.data.au);
    else{
      fprintf(stderr, "Error: value of AU is wrong = %.4lf km ... abort!\n",
        head1.data.au);
      exit(1);
    }
  }

  /** DO NOT CLOSE FILE POINTER AS IT IS USED LATER ON */

  return size;
}

void convert(REAL8 *gps_JD, REAL8 *time){
  REAL8 g=0.;
  REAL8 tdt2=0.;
  REAL8 diff=0., diff1=0.;

  tdt2 = gps_JD[1] + 51.184/86400.;
  diff1 = gps_JD[0] - 2451545.0;

  diff = diff1+tdt2;
  g = 6.24008 + 0.017202*diff;

  time[0] = gps_JD[0];
  time[1] = tdt2 + (0.001658*sin(g) + 1.4e-5*sin(2.*g))/86400.;

  if(time[1] >= 1.){
    time[1] = time[1] - 1.;
    time[0] = time[0] + 1.;
  }
}

/* this is taken from David Hoffman's ephem_read functions from
   http://www.projectpluto.com/jpl_eph.htm */
void read_coeffs(REAL8 *coeffArray, REAL8 *time, FILE *fp){
  REAL8 Tdelta = 0.;
  INT4 offset = 0;

  /* for first time in this time should be set to zero to read in the first set
     of coefficients */
  if( time[0] == 0. ){
    if( fread(&coeffArray[0], sizeof(REAL8), ARRAY_SIZE, fp) != ARRAY_SIZE ){
      fprintf(stderr, "Error reading in first set of coefficients.\n");
      exit(1);
    }

    /* swap bits */
    endian_swap((CHAR*)&coeffArray[0], sizeof(REAL8), ARRAY_SIZE);

    /* set beginning and end time of these coefficients */
    Tbeg = coeffArray[0];
    Tend = coeffArray[1];
    Tspan = Tend - Tbeg;

    if(verbose){
      fprintf(stdout, "First record in ephemeris file:\n  start time %lf\n  end\
 time %lf\n", Tbeg, Tend);
    }
  }
  else{
    /* find the record in the file that contains the input time conpared to the
       previous record (which gave the values of Tbeg and Tend) */
    if( time[0]+time[1] < Tbeg ){ /* this will be backwards from last point */
      Tdelta = Tbeg - time[0] + time[1];
      offset = -(INT4)ceil(Tdelta/Tspan);
    }

    if( time[0]+time[1] > Tend ){
      Tdelta = time[0] - Tend + time[1];
      offset = (INT4)ceil(Tdelta/Tspan);
    }

    /* get ephemeris coefficients from that record */
    fseek(fp, (offset-1)*ARRAY_SIZE*sizeof(REAL8), SEEK_CUR);
    if( fread(&coeffArray[0], sizeof(REAL8), ARRAY_SIZE, fp) != ARRAY_SIZE ){
      fprintf(stderr, "Error reading in of coefficients for the record.\n");
      exit(1);
    }

    /* swap bits */
    endian_swap((CHAR*)&coeffArray[0], sizeof(REAL8), ARRAY_SIZE);

    /* set new values of Tbeg and Tend */
    Tbeg = coeffArray[0];
    Tend = coeffArray[1];
    Tspan = Tend - Tbeg;
  }
}

void pleph(REAL8 *coeffArray, REAL8 *time, INT4 target, REAL8 *state, FILE *fp){
  REAL8 stateTemp[6];

  INT4 i=0;

  REAL8 au=CURT_AU;
  /* Curt's AU value from his code - this is slightly different from the one
     in the JPL binary ephemeris files (of order a few hundred metres), so I
     need to multiply things be a factor of au_Curt/au_JPLephem */

  /* if we're getting the Earth data then correct for the Moon */
  if( target == EARTH ){
    /* get the Earth-Moon barycenter */
    interpolate_state(coeffArray, time, target, stateTemp, fp);

    for (i=0; i<6; i++)
      state[i] = stateTemp[i];

    /* get moon position */
    interpolate_state(coeffArray, time, MOON, stateTemp, fp);

    /* remove the effect of the Moon */
    for (i=0; i<6; i++)
      state[i] -= stateTemp[i]/(1. + head1.data.emrat);
  }
  else if( target == MOON ){
    /* get the Moon position */
    interpolate_state(coeffArray, time, target, stateTemp, fp);

    for (i=0; i<6; i++)
      state[i] = stateTemp[i];

    /* get Earth-Moon barycenter position */
    interpolate_state(coeffArray, time, EARTH, stateTemp, fp);

    /* remove the effect of the Earth */
    for (i=0; i<6; i++)
      state[i] = stateTemp[i] + state[i]*(1.-1./(1.+head1.data.emrat));
  }
  else
    interpolate_state(coeffArray, time, target, state, fp);

  /* put values in light seconds and light seconds per second, from km and km
     per day*/
  for (i=0;i<3;i++){
    state[i] *= (au/head1.data.au)*1.e3/(REAL8)LAL_C_SI;
    state[i+3] *= (au/head1.data.au)*1.e3/((REAL8)LAL_C_SI*86400.);
  }
}

/* this function will compute the position and velocity vector of a given
planetary body from the Chebyshev coefficients - it does not compute nutations
of librations */
void interpolate_state(REAL8 *coeffArray, REAL8 *time, INT4 target,
  REAL8 *state, FILE *fp){
  REAL8 A[50], Cp[50], Psum[3], Vsum[3], Up[50];
  // unused: REAL8 B[50];
  REAL8 Tbreak = 0., Tseg = 0., Tsub=0., Tc=0.;

  INT4 i=0, j=0;
  INT4 C=0, G=0, N=0, offset=0;

  /* initialise local arrays */
  for(i=0; i<50 ;i++){
    A[i] = 0.;
    // unused: B[i] = 0.;
  }

  /* determin if we need a new record to be read, or if the current one is OK */
  if( time[0]+time[1] < Tbeg || time[0]+time[1] > Tend)
    read_coeffs(coeffArray, time, fp);

  /* read coefficients from the header record */
  C = head1.data.ipt[target][0] - 1;        /* coeff array entry point */
  N = head1.data.ipt[target][1];            /* number of coeff's */
  G = head1.data.ipt[target][2];            /* granules in current record */

  /* compute the normalised time of the Chebyshev polynomials */
  if ( G == 1 ){
    /* Tc = 2.*(time - Tbeg)/Tspan - 1.;*/
    Tc = (2.*(time[0]-Tbeg)/Tspan) + (2.*time[1]/Tspan) - 1.;

    for (i=C ; i<(C+3*N) ; i++)  A[i-C] = coeffArray[i];
  }
  else if ( G > 1 ){
    Tsub = Tspan/(REAL8)G;          /* compute subgranule interval */

    for ( j=G ; j>0 ; j-- ){
      Tbreak = Tbeg + ((REAL8) j-1) * Tsub;

      if ( time[0]+time[1] > Tbreak ){
        Tseg  = Tbreak;
        offset = j-1;
        break;
      }
    }

    /* Tc = 2.*(time - Tseg)/Tsub - 1.; */
    Tc = (2.*(time[0] - Tseg)/Tsub) + (2.*time[1]/Tsub) - 1.;
    C  = C + 3 * offset * N;

    for (i=C ; i<(C+3*N) ; i++) A[i-C] = coeffArray[i];
  }
  else{     /* something has gone terribly wrong */
    fprintf(stderr, "Error... number of granules must be >= 1: check header \
data.\n");
    exit(1);
  }

  /* calculate the position and velocity */
  for ( i=0 ; i<3 ; i++ ){    /* compute interpolating polynomials */
    Cp[0] = 1.;
    Cp[1] = Tc;
    Cp[2] = 2.*Tc*Tc - 1.;

    Up[0] = 0.;
    Up[1] = 1.;
    Up[2] = 4.*Tc;

    for ( j=3 ; j<N ; j++ ){
      Cp[j] = 2.*Tc*Cp[j-1] - Cp[j-2];
      Up[j] = 2.*Tc*Up[j-1] + 2.*Cp[j-1] - Up[j-2];
    }

    Psum[i] = 0.;   /* compute interpolated position & velocity */
    Vsum[i] = 0.;

    for ( j=N-1 ; j>-1 ; j-- )  Psum[i] = Psum[i] + A[j+i*N] * Cp[j];
    for ( j=N-1 ; j>0  ; j-- )  Vsum[i] = Vsum[i] + A[j+i*N] * Up[j];

    state[i] = Psum[i];
    state[i+3] = Vsum[i]*2.*((REAL8) G)/(Tspan);
  }
}


/* a little endian-swapper */
void endian_swap(CHAR *pdata, size_t dsize, size_t nelements){
  UINT4 i, j, indx;
  CHAR tempbyte;

  if (dsize <= 1) return;

  for (i=0; i<nelements; i++){
    indx = dsize;
    for (j=0; j<dsize/2; j++){
      tempbyte = pdata[j];
      indx = indx - 1;
      pdata[j] = pdata[indx];
      pdata[indx] = tempbyte;
    }

    pdata = pdata + dsize;
  }

  return;

} /* endian swap */


void get_input_args(inputParams_t *inputParams, INT4 argc, CHAR *argv[]){
  struct option long_options[] =
  {
    { "help",          no_argument,       0, 'h' },
    { "verbose",       no_argument, &verbose, 1 },
    { "ephem-file",    required_argument, 0, 'e'},
    { "output-file",   required_argument, 0, 'E'},
    { "year",          required_argument, 0, 'y'},
    { "interval",      required_argument, 0, 'i'},
    { "overlap",       required_argument, 0, 'n'},
    { "target",        required_argument, 0, 't'},
    { "num-years",     required_argument, 0, 'N'},
    { "test",          no_argument, &test, 1},
    { 0, 0, 0, 0 }
  };

  CHAR args[] = "he:E:y:i:n:t:N:";
  CHAR *program = argv[0];

  if(argc == 1){
    fprintf(stderr, "Error... no input parameters given! Try again!!!\n");
    exit(0);
  }

  /* set defaults */
  inputParams->nhre = 4; /* default to 4 hours between points */

  inputParams->noverlap = 0; /* default to no overlap */
  inputParams->nyears = 1;   /* default to one year of data */

  /* parse input arguements */
  while ( 1 ){
    INT4 option_index = 0;
    INT4 c;

    c = getopt_long_only( argc, argv, args, long_options, &option_index );
    if ( c == -1 ) /* end of options */
      break;

    switch( c ){
      case 0:
        if( long_options[option_index].flag )
          break;
        else
          fprintf(stderr, "Error passing option %s with argument %s\n",
            long_options[option_index].name, optarg);
      case 'h': /* help message */
        fprintf(stderr, USAGE, program);
        exit(0);
      case 'e':
        sprintf(inputParams->ephemfile, "%s", optarg);
        break;
      case 'E':
        sprintf(inputParams->outputfile, "%s", optarg);
        break;
      case 'y':
        inputParams->year = atoi(optarg);
        break;
      case 'i':
        inputParams->nhre = atof(optarg);
        break;
      case 'n':
        inputParams->noverlap = atoi(optarg);
        break;
      case 't':
        inputParams->targName = NULL;

        if((inputParams->targName=(CHAR*)strstr("MERCURY",optarg))!=NULL)
          inputParams->target = MERCURY;
        else if((inputParams->targName=(CHAR*)strstr("VENUS",optarg))!=NULL)
          inputParams->target = VENUS;
        else if((inputParams->targName=(CHAR*)strstr("EARTH",optarg))!=NULL)
          inputParams->target = EARTH;
        else if((inputParams->targName=(CHAR*)strstr("MARS",optarg))!=NULL)
          inputParams->target = MARS;
        else if((inputParams->targName=(CHAR*)strstr("JUPITER",optarg)) !=
          NULL)
          inputParams->target = JUPITER;
        else if((inputParams->targName=(CHAR*)strstr("SATURN",optarg)) !=
          NULL)
          inputParams->target = SATURN;
        else if((inputParams->targName=(CHAR*)strstr("URANUS",optarg)) !=
          NULL)
          inputParams->target = URANUS;
        else if((inputParams->targName=(CHAR*)strstr("NEPTUNE",optarg)) !=
          NULL)
          inputParams->target = NEPTUNE;
        else if((inputParams->targName=(CHAR*)strstr("PLUTO",optarg))!=NULL)
          inputParams->target = PLUTO;
        else if((inputParams->targName=(CHAR*)strstr("MOON",optarg))!=NULL)
          inputParams->target = MOON;
        else if((inputParams->targName=(CHAR*)strstr("SUN",optarg))!=NULL)
          inputParams->target = SUN;
        else{
          fprintf(stderr, "You must enter a valid solar system body!\n");
          exit(0);
        }
        break;
      case 'N':
        inputParams->nyears = atoi(optarg);
        break;
      case '?':
        fprintf(stderr, "Unknown error while parsing options\n");
      default:
        fprintf(stderr, "Unknown error while parsing options\n");
    }
  }
}
