/*
*  Copyright (C) 2012 Matt Pitkin
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

/* Many of these functions are (with occasional modification) taken directly 
 * from the TEMPO2 software package http://www.atnf.csiro.au/research/pulsar/tempo2/
 * written by George Hobbs and Russell Edwards */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <getopt.h>
#include <time.h>
#include <lal/Date.h>

/* TEMPO style time delay file */
#define TT2TDB_FILE "/ephemeris/TDB.1950.2050"
/* TEMPO2 style time delay file */
#define IFTEPH_FILE "/ephemeris/TIMEEPH_short.te405" 

#define MAXFNAME 512

#define MJDEPOCH 44244.0
#define DAYSTOSEC 86400.0
#define GPS2TT 51.184

/* macros taken from TEMPO2's ifteph.h file */
/* constants determined by Irwin & Fukushima */
#define IFTE_TEPH0 -65.564518e-6

//This is the value used by if99 : #define IFTE_KM1  1.55051979154e-8 */
// However we should use the IAU value of L_B that follows from
// their definition of L_G: L_B = 1.55051976772e-8, K=1/(1-L_B)
#define IFTE_KM1 1.55051979176e-8
#define IFTE_K (((long double)1.0) + ((long double)IFTE_KM1))

#define USAGE \
"Usage: %s [options]\n\n"\
" --help              display this message\n"\
" --verbose           display all error messages\n"\
" --ephem-type        TEMPO/TDB or TEMPO2/TCB/Teph\n"\
" --output-path       path into which to output the ascii ephemeris (the\n\
                     file name will get constructed from the ephemeris type\n\
                     and the start and end years of the input times.)\n"\
" --start             a start GPS time\n"\
" --end               an end GPS time\n"\
" --interval          an interval (in seconds) between entries\n"\
"\n"

typedef enum {
  TT2TDB,
  TT2TCB
} etype;

typedef struct taginputParams_t{
  char *ephemtype; /* type of ephemeris */
  char ephemfile[MAXFNAME]; /* path and name of binary ephemeris file */
  char *outputpath; /* path to output ephemeris file */

  double startT;   /* a start GPS time */
  double endT;     /* an end GPS time */
  double interval; /* number of seconds between output entries */
  
  etype et;
} inputParams_t;

void get_input_args(inputParams_t *inputParams, int argc, char *argv[]);

/*** FUNCTIONS TAKEN FROM TEMPO2 *****/

/* Fortran-equivalent functions */
double fortran_mod(double a, double p);
int open_file(char fname[MAXFNAME]);
double read_double(void);
int read_int(void);
void close_file(void);

double FB_deltaT(long double mjd_tt, char fname[MAXFNAME]);

double IF_deltaT(long double mjd_tt);

/* functions and structures for reading in Irwin and Fukushima ephemeris file */

struct IFTE_interpolation_info {
   double pc[18],vc[18], twot;
   int np, nv;
};

typedef struct {
  char title[256];
  double startJD, endJD, stepJD;
  int ephver;
  double L_C;
  int swap_endian;
  int reclen;
  int irec;
  double buf[322];
  FILE *f;
  struct IFTE_interpolation_info iinfo;
  int ipt[2][3];
} IFTEphemeris;

static IFTEphemeris ifte;

void IFTE_init(const char fname[MAXFNAME]);
void IFTE_close_file(void);
static void IFTEinterp( struct IFTE_interpolation_info *iinfo,
                        const double coef[], const double t[2], const int ncf,
                        const int ncm, const int na, const int ifl,
                        double posvel[] );

void IFTE_get_Vals(double JDeph0, double JDeph1, int kind,
                   double *res);
void IFTE_get_DeltaT_DeltaTDot(double Teph0, double Teph1,
                               double *DeltaT, double *DeltaTDot);
double IFTE_DeltaT(double Teph0, double Teph1);

/* functions to perform endian swapping */
void IFTswap8(char *dword);
void IFTswapDouble(double *dbl);
void IFTswapInts(int *word, int n);
void IFTswapInt(int *word);
void IFTswap4(char *word);
void IFTswap8N(char *dwords, int n);
void IFTswapDoubles(double *dbl, int N);
