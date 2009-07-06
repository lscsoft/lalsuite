/*
*  Copyright (C) 2007 Chris Messenger
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

/*-----------------------------------------------------------------------
 *
 * File Name: TEMPOcomparsion.c
 * Authors:  C. Messenger
 *
 * Revision: $Id$
 *
 * History:   Created by Messenger
 *           
 *
 *-----------------------------------------------------------------------
 */
#include <glob.h> 
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h> 
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/PulsarDataTypes.h>
#include <lal/UserInput.h>
#include <lal/LALInitBarycenter.h>
#include <lal/BinaryPulsarTiming.h>
#include <lal/GeneratePulsarSignal.h>
#include <lal/Random.h>

#include <lal/LogPrintf.h>

#include <lalapps.h>

RCSID( "$Id$");

/* ---------- Error codes and messages ---------- */
#define TEMPOCOMPARISONC_ENORM 0
#define TEMPOCOMPARISONC_ESUB  1
#define TEMPOCOMPARISONC_EARG  2
#define TEMPOCOMPARISONC_EBAD  3
#define TEMPOCOMPARISONC_EFILE 4
#define TEMPOCOMPARISONC_ENULL 5
#define TEMPOCOMPARISONC_EMEM  6
#define TEMPOCOMPARISONC_EINPUT 7

#define TEMPOCOMPARISONC_MSGENORM "Normal exit"
#define TEMPOCOMPARISONC_MSGESUB  "Subroutine failed"
#define TEMPOCOMPARISONC_MSGEARG  "Error parsing arguments"
#define TEMPOCOMPARISONC_MSGEBAD  "Bad argument values"
#define TEMPOCOMPARISONC_MSGEFILE "Could not create output file"
#define TEMPOCOMPARISONC_MSGENULL "Null Pointer"
#define TEMPOCOMPARISONC_MSGEMEM  "Out of memory"
#define TEMPOCOMPARISONC_MSGEINPUT  "Bad input"

/*---------- local defines ---------- */

#define TRUE (1==1)
#define FALSE (1==0)

#define GBT_LOCATION_X 882589.65    
#define GBT_LOCATION_Y -4924872.32    
#define GBT_LOCATION_Z 3943729.348

#define NARRABRI_LOCATION_X -4752329.7000
#define NARRABRI_LOCATION_Y 2790505.9340 
#define NARRABRI_LOCATION_Z -3200483.7470

#define ARECIBO_LOCATION_X 2390490.0
#define ARECIBO_LOCATION_Y  -5564764.0
#define ARECIBO_LOCATION_Z 1994727.0

#define NANSHAN_LOCATION_X -228310.702
#define NANSHAN_LOCATION_Y  4631922.905
#define NANSHAN_LOCATION_Z 4367064.059

#define DSS_43_LOCATION_X -4460892.6
#define DSS_43_LOCATION_Y 2682358.9 
#define DSS_43_LOCATION_Z -3674756.0

#define PARKES_LOCATION_X  -4554231.5
#define PARKES_LOCATION_Y 2816759.1
#define PARKES_LOCATION_Z -3454036.3

#define JODRELL_LOCATION_X 3822252.643
#define JODRELL_LOCATION_Y -153995.683
#define JODRELL_LOCATION_Z 5086051.443

#define VLA_LOCATION_X  -1601192.
#define VLA_LOCATION_Y -5041981.4
#define VLA_LOCATION_Z 3554871.4

#define NANCAY_LOCATION_X  4324165.81
#define NANCAY_LOCATION_Y 165927.11
#define NANCAY_LOCATION_Z 4670132.83

#define EFFELSBERG_LOCATION_X 4033949.5
#define EFFELSBERG_LOCATION_Y 486989.4    
#define EFFELSBERG_LOCATION_Z 4900430.8

#define JODRELLM4_LOCATION_X 3822252.643
#define JODRELLM4_LOCATION_Y -153995.683
#define JODRELLM4_LOCATION_Z 5086051.443

#define GB300_LOCATION_X 881856.58
#define GB300_LOCATION_Y -4925311.86 
#define GB300_LOCATION_Z 3943459.70

#define GB140_LOCATION_X 882872.57
#define GB140_LOCATION_Y -4924552.73
#define GB140_LOCATION_Z 3944154.92

#define GB853_LOCATION_X 882315.33
#define GB853_LOCATION_Y -4925191.41
#define GB853_LOCATION_Z 3943414.05

#define LA_PALMA_LOCATION_X 5327021.651
#define LA_PALMA_LOCATION_Y -1719555.576
#define LA_PALMA_LOCATION_Z 3051967.932

#define Hobart_LOCATION_X -3950077.96
#define Hobart_LOCATION_Y 2522377.31
#define Hobart_LOCATION_Z -4311667.52

#define Hartebeesthoek_LOCATION_X 5085442.780
#define Hartebeesthoek_LOCATION_Y 2668263.483
#define Hartebeesthoek_LOCATION_Z -2768697.034

#define WSRT_LOCATION_X 3828445.659
#define WSRT_LOCATION_Y 445223.600000
#define WSRT_LOCATION_Z 5064921.5677

#define COE_LOCATION_X  0.0
#define COE_LOCATION_Y 1.0
#define COE_LOCATION_Z 0.0

/*    882589.65    -4924872.32     3943729.348      GBT                 gbt   */
/* -4752329.7000  2790505.9340     -3200483.7470    NARRABRI            atca */
/*   2390490.0     -5564764.0      1994727.0        ARECIBO             ao */
/*  -228310.702   4631922.905      4367064.059      NANSHAN             nanshan */
/*  -4460892.6     2682358.9       -3674756.0       DSS_43              tid43 */
/*  -4554231.5     2816759.1       -3454036.3       PARKES              pks */
/*  3822252.643   -153995.683      5086051.443      JODRELL             jb */
/*  -1601192.     -5041981.4       3554871.4        VLA                 vla */
/*  4324165.81     165927.11       4670132.83       NANCAY              ncy */
/*  4033949.5      486989.4        4900430.8        EFFELSBERG          eff */
/*  3822252.643   -153995.683      5086051.443      JODRELLM4           jbm4 */
/*  881856.58     -4925311.86      3943459.70       GB300               gb300 */
/*  882872.57     -4924552.73      3944154.92       GB140               gb140 */
/*  882315.33     -4925191.41      3943414.05       GB853               gb853 */
/*  5327021.651    -1719555.576    3051967.932      LA_PALMA            lap */
/*  -3950077.96    2522377.31     -4311667.52       Hobart              hob */
/*  5085442.780   2668263.483     -2768697.034      Hartebeesthoek      hart */
/*  3828445.659  445223.600000  5064921.5677        WSRT                wsrt */
/* # */
/* ####### From Jodrell obsys.dat file */
/* # */
/*  383395.727    -173759.585      5077751.313      MKIII               j   */
/*  3817176.557   -162921.170      5089462.046      TABLEY              k   */
/*  3828714.504   -169458.987      5080647.749      DARNHALL            l   */
/*  3859711.492   -201995.082      5056134.285      KNOCKIN             m   */
/*  3923069.135   -146804.404      5009320.570      DEFFORD             n   */
/*        0.0           1.0              0.0        COE                 coe   */
/* # */
/* ###### Telescope ID changed from the Jodrell obsys.dat file */
/* # */
/*  3822473.365   -153692.318      5085851.303      JB_MKII             jbmk2 */
/*  3822294.825   -153862.275      5085987.071      JB_42ft             jb42 */
/* # */
/* # New telescopes */
/* # 284543.5      -175324.040      2400.            LA_PALMA            p */
/*  1719555.576   5327021.651      3051967.932      LA_PALMA            p */

/* ---------- local types ---------- */

/** user input variables */
typedef struct 
{ 
  BOOLEAN help;

  CHAR *RAJ;
  CHAR *DECJ;
  REAL8 TstartUTCMJD;
  REAL8 TrefTDBMJD;
  REAL8 DeltaTMJD;
  REAL8 DurationMJD;
  CHAR *det;
  CHAR *ephemyear;
  CHAR *ephemdir;
  REAL8 f0;
  REAL8 fdot;
  CHAR *PSRJ;
  CHAR *Observatory;
  BOOLEAN randparams;
  BOOLEAN randsampling;
  INT4 seed;

} UserVariables_t;

/* a time structure to accurately store MJD times so as not to lose any precision */
typedef struct 
{ 
  INT4 days;       /* integer MJD days */
  REAL8 fracdays;  /* fractional day */
} MJDTime;

/*---------- empty structs for initializations ----------*/
UserVariables_t empty_UserVariables;
/* ---------- global variables ----------*/

extern int vrbflg;
extern INT4 lalDebugLevel;


/* ---------- local prototypes ---------- */
void initUserVars (LALStatus *status, int argc, char *argv[], UserVariables_t *uvar);
void TDBMJDtoGPS(LALStatus *status,LIGOTimeGPS *GPS,MJDTime MJD);
void AddIntervaltoMJD(LALStatus *status,LALTimeInterval interval,MJDTime *MJDout, MJDTime MJDin);
void REAL8toMJD(LALStatus *status,MJDTime *MJD,REAL8 x);
void UTCMJDtoGPS(LALStatus *status,LIGOTimeGPS *GPS,MJDTime MJD,INT4 leap);
void UTCGPStoMJD(LALStatus *status,MJDTime *MJD,LIGOTimeGPS *GPS,INT4 leap);
void TDBGPStoMJD(LALStatus *status,MJDTime *MJD,LIGOTimeGPS GPS,INT4 leap);
void MJDtoREAL8(LALStatus *status,REAL8 *x,MJDTime MJD);
void deltaMJD(LALStatus *status,MJDTime *deltaMJD,MJDTime *x,MJDTime *y);
void LALRadstoRA(LALStatus *status, CHAR *RA, REAL8 rads);
void LALRadstoDEC(LALStatus *status, CHAR *DEC, REAL8 rads);
void GPStoTDBMJD(LALStatus *status,MJDTime *TDBMJD,LIGOTimeGPS GPS);
void MultiplyInterval(LALStatus *status,LALTimeInterval *newinterval,LALTimeInterval interval,INT4 factor);

/*============================================================
 * FUNCTION definitions
 *============================================================*/

int 
main(int argc, char *argv[]){ 

  static LALStatus       status;  /* LALStatus pointer */ 
  UserVariables_t uvar = empty_UserVariables;
  FILE *fp = NULL;
  BarycenterInput baryinput = empty_BarycenterInput;
  REAL8 alpha, delta;
  EphemerisData *edat = NULL;
  CHAR EphemEarth[1024], EphemSun[1024];
  LALLeapSecFormatAndAcc formatAndAcc = {LALLEAPSEC_GPSUTC, LALLEAPSEC_LOOSE};
  INT4 leap0,leap;
  LIGOTimeGPS epoch;
  LIGOTimeGPS TstartSSB, TendSSB, TendGPS;
  INT4 n;
  LIGOTimeGPS *TSSB = NULL;
  MJDTime TstartUTCMJD;
  LIGOTimeGPS TDET;
  REAL8 temp;
  INT4 i;
  MJDTime tempTOA;
  REAL8 dt;
  LIGOTimeGPS TstartGPS;
  MJDTime *TOA = NULL;
  CHAR tempstr[15];
  CHAR *tempstr2;
  CHAR TrefMJDstr[20],TstartMJDstr[20],TfinishMJDstr[20],TOAstr[20];
  PulsarSignalParams pulsarparams;
  CHAR parfile[256];
  CHAR timfile[256];
  CHAR detcode[16];
  REAL8 TstartUTCMJDtest;
  REAL8 diff;
  MJDTime MJDdiff, MJDtest;
  CHAR RA[256],DEC[256];
  MJDTime TrefTDBMJD;
  LIGOTimeGPS TrefSSB_TDB_GPS;
  LALTimeInterval dtinterval;

  /* LALDebugLevel must be called before anything else */
  lalDebugLevel = 1;
  vrbflg = 1;	/* verbose error-messages */

  /* set LAL error-handler */
  lal_errhandler = LAL_ERR_EXIT;

  LAL_CALL( LALGetDebugLevel( &status, argc, argv, 'd'), &status);

  /* set log-level */
  LogSetLevel ( lalDebugLevel );

  LAL_CALL (initUserVars (&status, argc, argv, &uvar), &status);	  

  /* exit if help was required */
  if (uvar.help)
    exit(0); 

  /* define sky position and/or convert sky position to radians */
  /* first - if we want a random sky position */
  if (uvar.randparams) {
    REAL8 sindelta;
    REAL4Vector *vector = NULL;
    RandomParams *params = NULL;
    REAL8 alphatest, deltatest;

    vector = XLALCreateREAL4Vector(2);
    LALCreateRandomParams( &status, &params, uvar.seed );
    
    /* fill vector with uniform deviates */
    for (i=0;i<(INT4)vector->length;i++) LALUniformDeviate(&status,&(vector->data[i]),params);
    
    alpha = LAL_TWOPI*vector->data[0];
    sindelta = -1.0 + 2.0*vector->data[1];
    delta = asin(sindelta);
    if (lalDebugLevel) fprintf(stdout,"STATUS : Randomly generate a sky position - alpha = %6.12f delta = %6.12f (rads)\n",alpha,delta);

    /* convert back to HMS representation for output */
    LALRadstoRA(&status,RA,alpha);
    LALRadstoDEC(&status,DEC,delta);
    if (lalDebugLevel) fprintf(stdout,"STATUS : Randomly generate a sky position - alpha = %s delta = %s (hms)\n",RA,DEC);

    /* convert back to rads to test conversion */
    alphatest = LALDegsToRads(RA, "alpha");
    deltatest = LALDegsToRads(DEC, "delta");
    if (lalDebugLevel) fprintf(stdout,"STATUS : Converted back to rads for testing - alpha = %6.12f delta = %6.12f (rads)\n",alphatest,deltatest);
    
    /* free memory */
    XLALDestroyVector(vector);
    LALDestroyRandomParams(&status,&params);

  }
  /* if a user has defined a sky position */
  else if ((LALUserVarWasSet(uvar.RAJ))&&(LALUserVarWasSet(uvar.DECJ))) {
    if (lalDebugLevel) fprintf(stdout,"STATUS : User defined sky position - alpha = %s delta = %s (hms)\n",uvar.RAJ,uvar.DECJ);
    alpha = LALDegsToRads(uvar.RAJ, "alpha");
    delta = LALDegsToRads(uvar.DECJ, "delta");
    if (lalDebugLevel) fprintf(stdout,"STATUS : Converted user defined sky position to - alpha = %6.12f delta = %6.12f (rads)\n",alpha,delta);
  }
  else {
    fprintf(stderr,"ERROR : must either set random params or define a sky position. Exiting.\n");
    return(TEMPOCOMPARISONC_EINPUT);
  }
  
  /* define start time in an MJD structure */
  REAL8toMJD(&status,&TstartUTCMJD,uvar.TstartUTCMJD);
  if (lalDebugLevel) fprintf(stdout,"STATUS : TstartUTCMJD converted to MJD days = %d fracdays = %6.12f\n",TstartUTCMJD.days,TstartUTCMJD.fracdays);
  
  /* convert back to test conversions */
  MJDtoREAL8(&status,&TstartUTCMJDtest,TstartUTCMJD);
  diff = uvar.TstartUTCMJD - TstartUTCMJDtest;
  if ( fabs(diff) > 1e-9) {
    fprintf(stderr,"ERROR : Time conversion gives discrepancy of %e sec. Exiting.\n",diff);
    return(TEMPOCOMPARISONC_ESUB);
  }
  if (lalDebugLevel) fprintf(stdout,"STATUS : MJD conversion gives discrepancies of %e sec\n",diff);

  /* use start time to define an epoch for the leap seconds */
  /* Note that epochs are defined in TDB !!! but here we only need to be rough to get a leap second value */
  TDBMJDtoGPS(&status,&epoch,TstartUTCMJD);
  if (lalDebugLevel) fprintf(stdout,"STATUS : leap second epoch = %d %d\n",epoch.gpsSeconds,epoch.gpsNanoSeconds);

  /* deal with ephemeris files and compute leap seconds */
  edat = LALCalloc(1, sizeof(EphemerisData));
  LALSnprintf(EphemEarth, 1024, "%s/earth%s.dat", uvar.ephemdir, uvar.ephemyear);
  LALSnprintf(EphemSun, 1024, "%s/sun%s.dat", uvar.ephemdir, uvar.ephemyear);
  edat->ephiles.earthEphemeris = EphemEarth;
  edat->ephiles.sunEphemeris = EphemSun;
  LALLeapSecs (&status, &leap0, &epoch, &formatAndAcc);
  edat->leap = (INT2) leap0;
  if (lalDebugLevel) fprintf(stdout,"STATUS : leap seconds = %d\n",leap0);
  
  /* initialise the barycenter routines */
  LALInitBarycenter(&status, edat);
  if (lalDebugLevel) fprintf(stdout,"STATUS : Initiated Barycenter\n");
  
  /* select detector location */
  if (strcmp(uvar.Observatory,"GBT")==0) {
    baryinput.site.location[0] = GBT_LOCATION_X;
    baryinput.site.location[1] = GBT_LOCATION_Y;
    baryinput.site.location[2] = GBT_LOCATION_Z;
    sprintf(detcode,"gbt"); 
  }
  else if (strcmp(uvar.Observatory,"NARRABRI")==0) {
    baryinput.site.location[0] = NARRABRI_LOCATION_X;
    baryinput.site.location[1] = NARRABRI_LOCATION_Y;
    baryinput.site.location[2] = NARRABRI_LOCATION_Z;
    sprintf(detcode,"atca"); 
  }
   else if (strcmp(uvar.Observatory,"ARECIBO")==0) {
    baryinput.site.location[0] = ARECIBO_LOCATION_X;
    baryinput.site.location[1] = ARECIBO_LOCATION_Y;
    baryinput.site.location[2] = ARECIBO_LOCATION_Z;
    sprintf(detcode,"ao"); 
  }
   else if (strcmp(uvar.Observatory,"NANSHAN")==0) {
    baryinput.site.location[0] = NANSHAN_LOCATION_X;
    baryinput.site.location[1] = NANSHAN_LOCATION_Y;
    baryinput.site.location[2] = NANSHAN_LOCATION_Z;
    sprintf(detcode,"nanshan"); 
  }
   else if (strcmp(uvar.Observatory,"DSS_43")==0) {
    baryinput.site.location[0] = DSS_43_LOCATION_X;
    baryinput.site.location[1] = DSS_43_LOCATION_Y;
    baryinput.site.location[2] = DSS_43_LOCATION_Z;
    sprintf(detcode,"tid43"); 
  }
  else if (strcmp(uvar.Observatory,"PARKES")==0) {
    baryinput.site.location[0] = PARKES_LOCATION_X;
    baryinput.site.location[1] = PARKES_LOCATION_Y;
    baryinput.site.location[2] = PARKES_LOCATION_Z;
    sprintf(detcode,"pks"); 
  }
   else if (strcmp(uvar.Observatory,"JODRELL")==0) {
    baryinput.site.location[0] = JODRELL_LOCATION_X;
    baryinput.site.location[1] = JODRELL_LOCATION_Y;
    baryinput.site.location[2] = JODRELL_LOCATION_Z;
    sprintf(detcode,"jb"); 
  }
   else if (strcmp(uvar.Observatory,"VLA")==0) {
    baryinput.site.location[0] = VLA_LOCATION_X;
    baryinput.site.location[1] = VLA_LOCATION_Y;
    baryinput.site.location[2] = VLA_LOCATION_Z;
    sprintf(detcode,"vla"); 
  }
   else if (strcmp(uvar.Observatory,"NANCAY")==0) {
    baryinput.site.location[0] = NANCAY_LOCATION_X;
    baryinput.site.location[1] = NANCAY_LOCATION_Y;
    baryinput.site.location[2] = NANCAY_LOCATION_Z;
    sprintf(detcode,"ncy"); 
  }
  else if (strcmp(uvar.Observatory,"EFFELSBERG")==0) {
    baryinput.site.location[0] = EFFELSBERG_LOCATION_X;
    baryinput.site.location[1] = EFFELSBERG_LOCATION_Y;
    baryinput.site.location[2] = EFFELSBERG_LOCATION_Z;
    sprintf(detcode,"eff"); 
  }
  else if (strcmp(uvar.Observatory,"JODRELLM4")==0) {
    baryinput.site.location[0] = JODRELLM4_LOCATION_X;
    baryinput.site.location[1] = JODRELLM4_LOCATION_Y;
    baryinput.site.location[2] = JODRELLM4_LOCATION_Z;
    sprintf(detcode,"jbm4"); 
  }
  else if (strcmp(uvar.Observatory,"GB300")==0) {
    baryinput.site.location[0] = GB300_LOCATION_X;
    baryinput.site.location[1] = GB300_LOCATION_Y;
    baryinput.site.location[2] = GB300_LOCATION_Z;
    sprintf(detcode,"gb300"); 
  }
  else if (strcmp(uvar.Observatory,"GB140")==0) {
    baryinput.site.location[0] = GB140_LOCATION_X;
    baryinput.site.location[1] = GB140_LOCATION_Y;
    baryinput.site.location[2] = GB140_LOCATION_Z;
    sprintf(detcode,"gb140"); 
  }
  else if (strcmp(uvar.Observatory,"GB853")==0) {
    baryinput.site.location[0] = GB853_LOCATION_X;
    baryinput.site.location[1] = GB853_LOCATION_Y;
    baryinput.site.location[2] = GB853_LOCATION_Z;
    sprintf(detcode,"gb853"); 
  }
  else if (strcmp(uvar.Observatory,"LA_PALMA")==0) {
    baryinput.site.location[0] = LA_PALMA_LOCATION_X;
    baryinput.site.location[1] = LA_PALMA_LOCATION_Y;
    baryinput.site.location[2] = LA_PALMA_LOCATION_Z;
    sprintf(detcode,"lap"); 
  }
  else if (strcmp(uvar.Observatory,"Hobart")==0) {
    baryinput.site.location[0] = Hobart_LOCATION_X;
    baryinput.site.location[1] = Hobart_LOCATION_Y;
    baryinput.site.location[2] = Hobart_LOCATION_Z;
    sprintf(detcode,"hob"); 
  }
  else if (strcmp(uvar.Observatory,"Hartebeesthoek")==0) {
    baryinput.site.location[0] = Hartebeesthoek_LOCATION_X;
    baryinput.site.location[1] = Hartebeesthoek_LOCATION_Y;
    baryinput.site.location[2] = Hartebeesthoek_LOCATION_Z;
    sprintf(detcode,"hart"); 
  }
  else if (strcmp(uvar.Observatory,"WSRT")==0) {
    baryinput.site.location[0] = WSRT_LOCATION_X;
    baryinput.site.location[1] = WSRT_LOCATION_Y;
    baryinput.site.location[2] = WSRT_LOCATION_Z;
    sprintf(detcode,"wsrt"); 
  }
  else if (strcmp(uvar.Observatory,"COE")==0) {
    baryinput.site.location[0] = COE_LOCATION_X;
    baryinput.site.location[1] = COE_LOCATION_Y;
    baryinput.site.location[2] = COE_LOCATION_Z;
    sprintf(detcode,"coe"); 
  }
  else if (strcmp(uvar.Observatory,"SSB")!=0) {
    fprintf(stderr,"ERROR. Unknown Observatory %s. Exiting.\n",uvar.Observatory);
    return(TEMPOCOMPARISONC_EFILE);
  }
  if (lalDebugLevel) fprintf(stdout,"STATUS : selected observatory %s - observatoryt code = %s\n",uvar.Observatory,detcode);
  if (lalDebugLevel) fprintf(stdout,"STATUS : baryinput location = %6.12f %6.12f %6.12f\n",baryinput.site.location[0],baryinput.site.location[1],baryinput.site.location[2]);

  /* convert start time to UTC GPS */
  UTCMJDtoGPS(&status,&TstartGPS,TstartUTCMJD,leap0);
  if (lalDebugLevel) fprintf(stdout,"STATUS : TstartGPS = %d %d\n",TstartGPS.gpsSeconds,TstartGPS.gpsNanoSeconds);

  /* convert back to test conversion */
  UTCGPStoMJD(&status,&MJDtest,&TstartGPS,leap0);
  deltaMJD(&status,&MJDdiff,&MJDtest,&TstartUTCMJD);
  diff = (MJDdiff.days+MJDdiff.fracdays)*86400;
  if ( fabs(diff)  > 1e-9) {
    fprintf(stderr,"ERROR : Time conversion gives discrepancy of %e sec. Exiting.\n",diff);
    return(TEMPOCOMPARISONC_ESUB);
  }
  if (lalDebugLevel) fprintf(stdout,"STATUS : MJD conversion gives discrepancies of %e sec\n",diff);
  
  /* define reference time in an MJD structure */
  REAL8toMJD(&status,&TrefTDBMJD,uvar.TrefTDBMJD);
  if (lalDebugLevel) fprintf(stdout,"STATUS : TrefTDBMJD converted to MJD days = %d fracdays = %6.12f\n",TrefTDBMJD.days,TrefTDBMJD.fracdays);

  /* convert reference time to TDB GPS */
  TDBMJDtoGPS(&status,&TrefSSB_TDB_GPS,TrefTDBMJD);
  if (lalDebugLevel) fprintf(stdout,"STATUS : TrefSSB_TDB_GPS = %d %d\n",TrefSSB_TDB_GPS.gpsSeconds,TrefSSB_TDB_GPS.gpsNanoSeconds);

  /* convert back to test conversion */
  TDBGPStoMJD(&status,&MJDtest,TrefSSB_TDB_GPS,leap0);
  deltaMJD(&status,&MJDdiff,&MJDtest,&TrefTDBMJD);
  diff = (MJDdiff.days+MJDdiff.fracdays)*86400;
  if ( fabs(diff)  > 1e-9) {
    fprintf(stderr,"ERROR : Time conversion gives discrepancy of %e sec. Exiting.\n",diff);
    return(TEMPOCOMPARISONC_ESUB);
  }
  if (lalDebugLevel) fprintf(stdout,"STATUS : MJD conversion gives discrepancies of %e sec\n",diff);

  /* fill in required pulsar params structure for Barycentering */
  pulsarparams.site = NULL;
  pulsarparams.site = (LALDetector *)LALMalloc(sizeof(LALDetector));
  pulsarparams.site->location[0] = baryinput.site.location[0];
  pulsarparams.site->location[1] = baryinput.site.location[1];
  pulsarparams.site->location[2] = baryinput.site.location[2];
  pulsarparams.pulsar.position.longitude = alpha;
  pulsarparams.pulsar.position.latitude = delta;
  pulsarparams.pulsar.position.system = COORDINATESYSTEM_EQUATORIAL;
  pulsarparams.ephemerides = edat;
  
  /* generate SSB initial TOA in GPS */
  LALConvertGPS2SSB (&status,&TstartSSB,TstartGPS,&pulsarparams); 
  if (lalDebugLevel) fprintf(stdout,"STATUS : TstartSSB = %d %d\n",TstartSSB.gpsSeconds,TstartSSB.gpsNanoSeconds);
  
  /* force integer seconds */
  /* TstartSSB.gpsNanoSeconds = 0; */
/*   if (lalDebugLevel) fprintf(stdout,"STATUS : (after rounding down to integer seconds) TstartSSB = %d %d\n",TstartSSB.gpsSeconds,TstartSSB.gpsNanoSeconds); */

  /* define TOA end time in GPS */
  temp = uvar.DurationMJD*86400.0;
  TendGPS = TstartGPS;
  XLALGPSAdd(&TendGPS, temp);
  if (lalDebugLevel) fprintf(stdout,"STATUS : GPS end time of TOAs = %d %d\n",TendGPS.gpsSeconds,TendGPS.gpsNanoSeconds);

  /* generate SSB end time in GPS (force integer seconds) */
  LALConvertGPS2SSB (&status,&TendSSB,TendGPS,&pulsarparams); 
  if (lalDebugLevel) fprintf(stdout,"STATUS : TendSSB = %d %d\n",TendSSB.gpsSeconds,TendSSB.gpsNanoSeconds);
  
  /* force integer seconds */
  /* TendSSB.gpsNanoSeconds = 0; */
/*   if (lalDebugLevel) fprintf(stdout,"STATUS : (after rounding down to integer seconds) TendSSB = %d %d\n",TendSSB.gpsSeconds,TendSSB.gpsNanoSeconds); */

  /* define TOA seperation in the SSB */ 
  dt = uvar.DeltaTMJD*86400.0;
  LALFloatToInterval(&status,&dtinterval,&dt);
  n = (INT4)ceil(uvar.DurationMJD/uvar.DeltaTMJD);
  if (lalDebugLevel) fprintf(stdout,"STATUS : TOA seperation at SSB = %d sec %d nano\n",dtinterval.seconds,dtinterval.nanoSeconds);
  if (lalDebugLevel) fprintf(stdout,"STATUS : number of TOAs to generate = %d\n",n);
  
  /* allocate memory for artificial SSB TOAs */
  TSSB = (LIGOTimeGPS *)LALMalloc(n*sizeof(LIGOTimeGPS));
  TOA = (MJDTime *)LALMalloc(n*sizeof(MJDTime));
  
  /* generate artificial SSB TOAs given the phase model phi = 2*pi*(f0*(t-tref) + 0.5*fdot*(t-tref)^2) */
  for  (i=0;i<n;i++) {

    LALTimeInterval dtstart;
    REAL8 dtref,fnow,cyclefrac,dtcor;
    LIGOTimeGPS tnow;

    /* define current interval */
    MultiplyInterval(&status,&dtstart,dtinterval,i);
    if (lalDebugLevel) fprintf(stdout,"STATUS : current (t-tstart) = %d %d\n",dtstart.seconds,dtstart.nanoSeconds);

    /* define current t */
    LALIncrementGPS(&status,&tnow,&TstartSSB,&dtstart);
    if (lalDebugLevel) fprintf(stdout,"STATUS : current t = %d %d\n",tnow.gpsSeconds,tnow.gpsNanoSeconds);

    /* define current t-tref */
    dtref = XLALGPSDiff(&tnow,&TrefSSB_TDB_GPS);
    if (lalDebugLevel) fprintf(stdout,"STATUS : current (t - tref) = %9.12f\n",dtref);

    dtcor = 1;
    while (dtcor>1e-9) {
      
      /* define actual cycle fraction at requested time */
      cyclefrac = fmod(uvar.f0*dtref + 0.5*uvar.fdot*dtref*dtref,1.0);
      if (lalDebugLevel) fprintf(stdout,"STATUS : cyclefrac = %9.12f\n",cyclefrac);

      /* define instantaneous frequency */
      fnow = uvar.f0 + uvar.fdot*dtref;
      if (lalDebugLevel) fprintf(stdout,"STATUS : instananeous frequency = %9.12f\n",fnow);

      /* add correction to time */
      dtcor = cyclefrac/fnow;
      dtref -= dtcor;
      if (lalDebugLevel) fprintf(stdout,"STATUS : timing correction = %9.12f\n",dtcor);
      if (lalDebugLevel) fprintf(stdout,"STATUS : corrected dtref to = %9.12f\n",dtref);
    }

    /* define time of zero phase */
    TSSB[i] = TrefSSB_TDB_GPS;
    XLALGPSAdd(&TSSB[i], dtref);
    if (lalDebugLevel) fprintf(stdout,"STATUS : TSSB[%d] = %d %d\n",i,TSSB[i].gpsSeconds,TSSB[i].gpsNanoSeconds);
 
  }
  
  /* loop over SSB time of arrivals and compute detector time of arrivals */
  for (i=0;i<n;i++) {
    
      LIGOTimeGPS TSSBtest;
      LIGOTimeGPS GPStest;
      LALTimeInterval interval;

      /* convert SSB to Detector time */
      LALConvertSSB2GPS (&status,&TDET,TSSB[i],&pulsarparams); 
      if (lalDebugLevel) fprintf(stdout,"STATUS : converted SSB TOA %d %d -> Detector TOA %d %d\n",TSSB[i].gpsSeconds,TSSB[i].gpsNanoSeconds,TDET.gpsSeconds,TDET.gpsNanoSeconds);
      
      /* convert back for testing conversion */
      LALConvertGPS2SSB (&status,&TSSBtest,TDET,&pulsarparams); 
      LALDeltaGPS (&status,&interval,&TSSBtest,&TSSB[i]);
      diff = (REAL8)interval.seconds + 1e-9*(REAL8)interval.nanoSeconds;
      if ( fabs(diff)  > 1e-9) {
	fprintf(stderr,"ERROR : Time conversion gives discrepancy of %e sec. Exiting.\n",diff);
	return(TEMPOCOMPARISONC_ESUB);
      }
      if (lalDebugLevel) fprintf(stdout,"STATUS : SSB -> detector conversion gives discrepancies of %e sec\n",diff);
      
      /* recompute leap seconds incase they've changed */
      LALLeapSecs (&status, &leap, &TDET, &formatAndAcc); 
      
      /* must now convert to an MJD time for TEMPO */
      /* Using UTC conversion as used by Matt in his successful comparison */
      UTCGPStoMJD (&status,&tempTOA,&TDET,leap);
      if (lalDebugLevel) fprintf(stdout,"STATUS : output MJD time = %d %6.12f\n",tempTOA.days,tempTOA.fracdays);

      /* convert back to test conversion */
      UTCMJDtoGPS (&status,&GPStest,tempTOA,leap);
      LALDeltaGPS (&status,&interval,&TDET,&GPStest);
      diff = (REAL8)interval.seconds + 1e-9*(REAL8)interval.nanoSeconds;
      if ( fabs(diff)  > 1e-9) {
	fprintf(stderr,"ERROR. Time conversion gives discrepancy of %e sec. Exiting.\n",diff);
	return(TEMPOCOMPARISONC_ESUB);
      }
      if (lalDebugLevel) fprintf(stdout,"STATUS : MJD time conversion gives discrepancies of %e sec\n",diff);

      /* fill in results */
      TOA[i].days = tempTOA.days;
      TOA[i].fracdays = tempTOA.fracdays;
       
  }
     
  /* convert the start time of the TOAs to a TDB time for use as a reference time in the TEMPO par file */
  /* as we have a monochromatic signal this reference time is meaningless */
  /* GPStoTDBMJD(&status,&TrefMJD,TstartGPS); */

  /* convert MJDTime structures into strings for precision output */
  /* snprintf(tempstr,15,"%1.13f",TrefMJD.fracdays); */
/*   tempstr2 = tempstr+2; */
/*   snprintf(TrefMJDstr,19,"%d.%s",TrefMJD.days,tempstr2);  */
/*   if (lalDebugLevel) fprintf(stdout,"STATUS : Converted reference MJD %d %6.12f to the string %s\n",TrefMJD.days,TrefMJD.fracdays,TrefMJDstr); */

  snprintf(tempstr,15,"%1.13f",TOA[0].fracdays);
  tempstr2 = tempstr+2;
  snprintf(TstartMJDstr,19,"%d.%s",TOA[0].days,tempstr2); 
  if (lalDebugLevel) fprintf(stdout,"STATUS : Converted initial TOA MJD %d %6.12f to the string %s\n",TOA[0].days,TOA[0].fracdays,TstartMJDstr);

  snprintf(tempstr,15,"%1.13f",TOA[n-1].fracdays);
  tempstr2 = tempstr+2;
  snprintf(TfinishMJDstr,19,"%d.%s",TOA[n-1].days,tempstr2); 
  if (lalDebugLevel) fprintf(stdout,"*** Converted MJD to a string %s\n",TfinishMJDstr);
  if (lalDebugLevel) fprintf(stdout,"STATUS : Converted final TOA MJD %d %6.12f to the string %s\n",TOA[n-1].days,TOA[n-1].fracdays,TfinishMJDstr);

  /* define output file names */
  sprintf(parfile,"%s.par",uvar.PSRJ);
  sprintf(timfile,"%s.tim",uvar.PSRJ);

  /* output to par file in format required by TEMPO 2 */
  if ((fp = fopen(parfile,"w")) == NULL) {
    fprintf(stderr,"ERROR. Could not open file %s. Exiting.\n",parfile);
    return(TEMPOCOMPARISONC_EFILE);
  }
  fprintf(fp,"PSRJ\t%s\n",uvar.PSRJ);
  if (uvar.randparams) {
    fprintf(fp,"RAJ\t%s\t1\n",RA);
    fprintf(fp,"DECJ\t%s\t1\n",DEC);
  }
  else {
    fprintf(fp,"RAJ\t%s\t1\n",uvar.RAJ);
    fprintf(fp,"DECJ\t%s\t1\n",uvar.DECJ);
  }
  fprintf(fp,"PEPOCH\t%6.12f\n",uvar.TrefTDBMJD);
  fprintf(fp,"POSEPOCH\t%6.12f\n",uvar.TrefTDBMJD);
  fprintf(fp,"DMEPOCH\t%6.12f\n",uvar.TrefTDBMJD);
  fprintf(fp,"DM\t0.0\n");
  fprintf(fp,"F0\t%6.16f\t1\n",uvar.f0);
  fprintf(fp,"F1\t%6.16f\t0\n",uvar.fdot);
  fprintf(fp,"START\t%s\n",TstartMJDstr);
  fprintf(fp,"FINISH\t%s\n",TfinishMJDstr);
  fprintf(fp,"TZRSITE\t%s\n",detcode);
  fprintf(fp,"CLK\tUTC(NIST)\n");
  fprintf(fp,"EPHEM\tDE405\n");
  fprintf(fp,"UNITS\tTDB\n");
  fprintf(fp,"MODE\t0\n");
  
  /* close par file */
  fclose(fp);
  
  /* output to tim file in format required by TEMPO 2 */
  if ((fp = fopen(timfile,"w")) == NULL) {
    fprintf(stderr,"ERROR. Could not open file %s. Exiting.\n",timfile);
    return(TEMPOCOMPARISONC_EFILE);
  }

  fprintf(fp,"FORMAT 1\n");
  for (i=0;i<n;i++) {
    
    /* convert each TOA to a string for output */ 
    snprintf(tempstr,18,"%1.16f",TOA[i].fracdays);
    tempstr2 = tempstr+2;
    snprintf(TOAstr,22,"%d.%s",TOA[i].days,tempstr2);
    fprintf(fp,"blank.dat\t1000.0\t%s\t1.0\t%s\n",TOAstr,detcode);
    if (lalDebugLevel) fprintf(stdout,"STATUS : Converting MJD time %d %6.16f to string %s\n",TOA[i].days,TOA[i].fracdays,TOAstr);
    
  }
  
  /* close tim file */
  fclose(fp);

  /* free memory */
  LALFree(TSSB);
  LALFree(TOA);
  LALFree(pulsarparams.site);
  LALFree(edat->ephemS);
  LALFree(edat->ephemE);
  LALFree(edat);
  LALDestroyUserVars (&status);
  LALCheckMemoryLeaks(); 
  
  return TEMPOCOMPARISONC_ENORM;
  
} /* main() */


/** register all "user-variables" */
void
initUserVars (LALStatus *status, int argc, char *argv[], UserVariables_t *uvar)
{
  INITSTATUS( status, "initUserVars", rcsid );
  ATTATCHSTATUSPTR (status);

  /* set a few defaults */
  uvar->help = FALSE;

  uvar->RAJ = (CHAR*)LALMalloc(512);
  sprintf(uvar->RAJ,"00:00.00.0000");

  uvar->DECJ = (CHAR*)LALMalloc(512);
  sprintf(uvar->DECJ,"00:00.00.0000");

  uvar->TstartUTCMJD = 53400;
  uvar->TrefTDBMJD = 53400;
  uvar->DeltaTMJD = 1;
  uvar->DurationMJD = 1800;

  uvar->ephemyear = (CHAR*)LALMalloc(512);
  sprintf(uvar->ephemyear,"05-09");

  uvar->f0 = 1.0;
  uvar->fdot = 0.0;

  uvar->PSRJ = (CHAR*)LALMalloc(512);
  sprintf(uvar->PSRJ,"TEST");

  uvar->Observatory = (CHAR*)LALMalloc(512);
  sprintf(uvar->Observatory,"JODRELL");

  uvar->randparams = FALSE;
  uvar->seed = 0;

  /* register user input variables */
  LALregBOOLUserStruct ( status, 	help, 		'h', UVAR_HELP,    	"Print this message" );
  LALregSTRINGUserStruct ( status, 	RAJ, 	        'r', UVAR_OPTIONAL, 	"Right ascension hh:mm.ss.ssss [Default = 00:00.00.0000]");
  LALregSTRINGUserStruct ( status, 	DECJ, 	        'j', UVAR_OPTIONAL, 	"Declination dd:mm.ss.ssss [Default = 00:00.00.0000]");
  LALregSTRINGUserStruct ( status,      ephemdir,       'e', UVAR_REQUIRED, 	"Name of ephemeris directory");
  LALregSTRINGUserStruct ( status,      ephemyear,      'y', UVAR_OPTIONAL, 	"Ephemeris year [Default = 05-09]");
  LALregREALUserStruct ( status, 	f0,     	'f', UVAR_OPTIONAL, 	"The signal frequency in Hz at SSB at the reference time [Default = 1.0]");
  LALregREALUserStruct ( status, 	fdot,     	'p', UVAR_OPTIONAL, 	"The signal frequency derivitive in Hz at SSB at the reference time [Default = 0.0]");
  LALregREALUserStruct ( status, 	TrefTDBMJD, 	'R', UVAR_OPTIONAL, 	"Reference time at the SSB in TDB in MJD [Default = 53400 ~ Jan 2005]");
  LALregREALUserStruct ( status, 	TstartUTCMJD, 	'T', UVAR_OPTIONAL, 	"Start time of output TOAs in UTC [Default = 53400 ~ Jan 2005]");
  LALregREALUserStruct ( status, 	DeltaTMJD, 	't', UVAR_OPTIONAL, 	"Time inbetween TOAs (in days) [DEFAULT = 1]");
  LALregREALUserStruct ( status, 	DurationMJD, 	'D', UVAR_OPTIONAL, 	"Full duration of TOAs (in days) [Default = 1800 ~ 5 years]");
  LALregSTRINGUserStruct ( status,      PSRJ,           'n', UVAR_OPTIONAL, 	"Name of pulsar [Default = TEST]");
  LALregSTRINGUserStruct ( status,      Observatory,    'O', UVAR_OPTIONAL, 	"TEMPO observatory name (GBT,ARECIBO,NARRABRI,NANSHAN,DSS_43,PARKES,JODRELL,VLA,NANCAY,COE,SSB) [Default = JODRELL]");
  LALregBOOLUserStruct ( status,        randparams,     'r', UVAR_OPTIONAL, 	"Override sky position with random values [Default = FALSE]");
  LALregINTUserStruct ( status, 	seed,     	'o', UVAR_OPTIONAL, 	"The random seed (integer) [Default = 0 = clock]");

  /* read all command line variables */
  TRY( LALUserVarReadAllInput(status->statusPtr, argc, argv), status);

  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* initUserVars() */



/* ------------------------------------------------------------------------------- */
/* the following functions have been written to maintain a high level of precision */
/* in storing MJD times */ 
/* ------------------------------------------------------------------------------- */

/* this function takes a REAL8 input MJD time and returns the same time */
/* but stored in an MJFTime structure to avoid loss of precision.       */
/* A REAL8 can only store a present day MJD time to 12 decimal figure   */
/* corresponding to 1e-7 sec */
void REAL8toMJD(LALStatus *status,MJDTime *MJD,REAL8 x) {

  INITSTATUS( status, "REAL8toMJD", rcsid );
  ATTATCHSTATUSPTR (status);

  /* take integer part of input time */
  MJD->days = (INT4)floor(x);
  MJD->fracdays = fmod(x,1.0);

  DETATCHSTATUSPTR (status);
  RETURN (status);

}

void MJDtoREAL8(LALStatus *status,REAL8 *x,MJDTime MJD) {

  INITSTATUS( status, "MJDtoREAL8", rcsid );
  ATTATCHSTATUSPTR (status);

  *x = (REAL8)MJD.days + (REAL8)MJD.fracdays*86400.0;

  DETATCHSTATUSPTR (status);
  RETURN (status);

}

void deltaMJD(LALStatus *status,MJDTime *dMJD,MJDTime *x,MJDTime *y) {

  MJDTime z;
  
  INITSTATUS( status, "deltaMJD", rcsid );
  ATTATCHSTATUSPTR (status);

  /* remove the days part from the other time */
  z.days = x->days - y->days;
  z.fracdays = x->fracdays;

  /* remove frac days part */
  z.fracdays = z.fracdays - y->fracdays;
  if (z.fracdays<0) {
    z.fracdays = 1.0 + z.fracdays;
    z.days--;
  }

  dMJD->days = z.days;
  dMJD->fracdays = z.fracdays;

  DETATCHSTATUSPTR (status);
  RETURN (status);

}

/* this function adds a LALTimeInterval in sec.nanosec to an MJDTime structure */
void AddIntervaltoMJD(LALStatus *status,LALTimeInterval interval,MJDTime *MJDout,MJDTime MJDin) {

  REAL8 fracdays = 0.0;
  INT4 days = 0;
  REAL8 temp,sfd,nfd;

  INITSTATUS( status, "AddIntervaltoMJD", rcsid );
  ATTATCHSTATUSPTR (status);

  /* convert seconds part to fractional days */
  sfd = (REAL8)interval.seconds/86400.0;

  temp = MJDin.fracdays + sfd;
  fracdays = fmod(temp,1.0); 
  days = MJDin.days + (INT4)floor(temp);
 
  /* convert nanoseconds part to fractional days */
  nfd = (REAL8)interval.nanoSeconds/(1e9*86400.0);

  temp = fracdays + nfd;
  fracdays = fmod(temp,1.0);
  days += (INT4)floor(temp);

  MJDout->days = days;
  MJDout->fracdays = fracdays;

  DETATCHSTATUSPTR (status);
  RETURN (status);
  
}
void GPStoTDBMJD(LALStatus *status,MJDTime *TDBMJD,LIGOTimeGPS GPS) {

  REAL8 Tdiff, dtrel;
  REAL8 MJDtemp;
  MJDTime MJDtest;
  LALTimeInterval interval;
  LIGOTimeGPS GPStemp;

  INITSTATUS( status, "GPStoTDBMJD", rcsid );
  ATTATCHSTATUSPTR (status);
  
  /* straight forward conversion from GPS to MJD */
  /* we do not need the more accurate MJDtime structure */
  MJDtemp = ((REAL8)GPS.gpsSeconds + 1e-9*(REAL8)GPS.gpsNanoSeconds)/86400.0 + 44244.0;

  /* Check not before the start of GPS time (MJD 44222) */
  if(MJDtemp < 44244.){
    fprintf(stderr, "Input time is not in range.\n");
    exit(0);
  } 
  
  /* compute the relativistic effect in TDB */
  Tdiff = MJDtemp + (2400000.5-2451545.0);
  /* meanAnomaly = 357.53 + 0.98560028*Tdiff; */ /* mean anomaly in degrees */
  /* meanAnomaly *= LAL_PI_180; */ /* mean anomaly in rads */
  /* dtrel = 0.001658*sin(meanAnomaly) + 0.000014*sin(2.*meanAnomaly); */ /* time diff in seconds */
  dtrel = 0.0016568*sin((357.5 + 0.98560028*Tdiff) * LAL_PI_180) +
    0.0000224*sin((246.0 + 0.90251882*Tdiff) * LAL_PI_180) +
    0.0000138*sin((355.0 + 1.97121697*Tdiff) * LAL_PI_180) +
    0.0000048*sin((25.0 + 0.08309103*Tdiff) * LAL_PI_180) +
    0.0000047*sin((230.0 + 0.95215058*Tdiff) *LAL_PI_180);

  /* define interval that is the magical number factor of 32.184 + 19 leap seconds to the */
  /* start of GPS time plus the TDB-TT correction and add it to the GPS input MJD */
  /* max absolute value of TDB-TT correction is < 0.184 sec */
  interval.seconds = 51;
  interval.nanoSeconds = (INT4)(0.184*1e9) + (INT4)(dtrel*1e9);

  /* add this to the input GPS time */
  /* this time is now in TDB but represented as a GPS structure */
  LALIncrementGPS(status->statusPtr,&GPStemp,&GPS,&interval); 

  /* convert to an MJD structure */
  MJDtest.days = (INT4)floor((REAL8)GPStemp.gpsSeconds/86400.0) + 44244;
  MJDtest.fracdays = fmod((REAL8)GPStemp.gpsSeconds/86400.0,1.0);
  interval.seconds = 0;
  interval.nanoSeconds = GPStemp.gpsNanoSeconds;
  AddIntervaltoMJD(status->statusPtr,interval,TDBMJD,MJDtest);

  DETATCHSTATUSPTR (status);
  RETURN (status);

}

void TDBMJDtoGPS(LALStatus *status,LIGOTimeGPS *GPS,MJDTime MJD){
 
  REAL8 Tdiff, TDBtoTT;
  REAL8 MJDtemp;
  LALTimeInterval interval;
  LIGOTimeGPS GPStemp;
  INT4 tempsec, tempnano;

  INITSTATUS( status, "TDBMJDtoGPS", rcsid );
  ATTATCHSTATUSPTR (status);

  /* convert MJD to REAL8 for calculation of the TDB-TT factor which */
  /* is small enough to allow not using the more accurate MJDtime structure */
  MJDtemp = (REAL8)MJD.days + MJD.fracdays;
 /*  printf("\tMJDtemp = %6.12f\n",MJDtemp); */

  /* Check not before the start of GPS time (MJD 44222) */
  if(MJDtemp < 44244.){
    fprintf(stderr, "Input time is not in range.\n");
    exit(0);
  } 
  
  Tdiff = MJDtemp + (2400000.5-2451545.0);
  /* meanAnomaly = 357.53 + 0.98560028*Tdiff; */ /* mean anomaly in degrees */
 /*  printf("\tmeanAnomaly (deg) = %6.12f\n",meanAnomaly); */
  /* meanAnomaly *= LAL_PI_180; */ /* mean anomaly in rads */
  /* TDBtoTT = 0.001658*sin(meanAnomaly) + 0.000014*sin(2.*meanAnomaly); */ /* time diff in seconds */
  TDBtoTT = 0.0016568*sin((357.5 + 0.98560028*Tdiff) * LAL_PI_180) +
    0.0000224*sin((246.0 + 0.90251882*Tdiff) * LAL_PI_180) +
    0.0000138*sin((355.0 + 1.97121697*Tdiff) * LAL_PI_180) +
    0.0000048*sin((25.0 + 0.08309103*Tdiff) * LAL_PI_180) +
    0.0000047*sin((230.0 + 0.95215058*Tdiff) *LAL_PI_180);

/*  printf("\tTdiff = %6.12f meanAnomoly (rads) = %6.12f TDBtoTT = %6.12f\n",Tdiff,meanAnomaly,TDBtoTT); */

  /* convert MJDtime to GPS with no corrections (yet) */
  tempsec = (MJD.days-44244)*86400;
 /*  printf("\ttempsec = %d\n",tempsec); */
  tempsec += (INT4)floor(MJD.fracdays*86400.0);
 /*  printf("\ttempsec = %d\n",tempsec); */
  tempnano = (INT4)(1e9*(MJD.fracdays*86400.0 - (INT4)floor(MJD.fracdays*86400.0)));
 /*  printf("\ttempnano = %d\n",tempnano); */
  GPStemp.gpsSeconds = tempsec;
  GPStemp.gpsNanoSeconds = tempnano;

  /* define interval that is the magical number factor of 32.184 + 19 leap seconds to the */
  /* start of GPS time plus the TDB-TT correction and minus it from the input MJD */
  /* max absolute value of TDB-TT correction is < 0.184 sec */
  interval.seconds = 51;
  interval.nanoSeconds = 184000000 + (INT4)(TDBtoTT*1e9);
 /*  printf("\tinterval nanoseconds = %d\n",interval.nanoSeconds); */
  LALDecrementGPS(status->statusPtr,GPS,&GPStemp,&interval);

  DETATCHSTATUSPTR (status);
  RETURN (status);
  
}

void UTCMJDtoGPS(LALStatus *status,LIGOTimeGPS *GPS,MJDTime MJD,INT4 leap){
 
  REAL8 MJDtemp;
  INT4 tempsec, tempnano;

  INITSTATUS( status, "UTCMJDtoGPS", rcsid );
  ATTATCHSTATUSPTR (status);

  /* convert MJD to REAL8 for error checking */
  MJDtemp = (REAL8)MJD.days + MJD.fracdays;
 
  /* Check not before the start of GPS time (MJD 44222) */
  if(MJDtemp < 44244.){
    fprintf(stderr, "Input time is not in range.\n");
    exit(0);
  } 

  /* convert MJDtime to GPS with no corrections (yet) */
  tempsec = (MJD.days-44244)*86400 + leap;
/*   printf("\ttempsec = %d\n",tempsec); */
  tempsec += (INT4)floor(MJD.fracdays*86400.0);
 /*  printf("\ttempsec = %d\n",tempsec); */
  tempnano = (INT4)(1e9*(MJD.fracdays*86400.0 - (INT4)floor(MJD.fracdays*86400.0)));
/*   printf("\ttempnano = %d\n",tempnano); */
  GPS->gpsSeconds = tempsec;
  GPS->gpsNanoSeconds = tempnano;

  DETATCHSTATUSPTR (status);
  RETURN (status);
  
}

void UTCGPStoMJD(LALStatus *status,MJDTime *MJD,LIGOTimeGPS *GPS,INT4 leap){
 
  LIGOTimeGPS tempGPS;
  LALTimeInterval interval;

  INITSTATUS( status, "UTCGPStoMJD", rcsid );
  ATTATCHSTATUSPTR (status);

  /* compute integer MJD days */
  MJD->days = 44244 + (INT4)floor(((REAL8)GPS->gpsSeconds + 1e-9*(REAL8)GPS->gpsNanoSeconds - (REAL8)leap)/86400.0);
  
  /* compute corresponding exact GPS for this integer days MJD */
  tempGPS.gpsSeconds = (MJD->days-44244)*86400 + leap;
  tempGPS.gpsNanoSeconds = 0;

  /* compute the difference as an interval */
  LALDeltaGPS (status->statusPtr,&interval,GPS,&tempGPS);

  /* convert the interval to fractional MJD days */
  MJD->fracdays = ((REAL8)interval.seconds + 1e-9*(REAL8)interval.nanoSeconds)/86400.0;

  /* printf("inside MJD = %d %6.12f\n",MJD->days,MJD->fracdays); */

  DETATCHSTATUSPTR (status);
  RETURN (status);
  
}

void TDBGPStoMJD(LALStatus *status,MJDTime *MJD,LIGOTimeGPS GPS,INT4 leap) {

  MJDTime MJDrough;
  LALTimeInterval interval;
  LIGOTimeGPS GPSguess;
  
  INITSTATUS( status, "TDBGPStoMJD", rcsid );
  ATTATCHSTATUSPTR (status);
  
  /* do rough conversion to MJD */
  UTCGPStoMJD(status->statusPtr,&MJDrough,&GPS,leap);

  while (1) {
    
    /* convert rough MJD guess correctly back to GPS */
    TDBMJDtoGPS(status->statusPtr,&GPSguess,MJDrough);
    
    /* compute difference */
    LALDeltaGPS (status->statusPtr,&interval,&GPS,&GPSguess);
    
    /* add difference to MJD */
    AddIntervaltoMJD(status->statusPtr,interval,MJD,MJDrough);
   
    /* update current guess */
    MJDrough = *MJD;

    if ((interval.seconds==0)&&(interval.nanoSeconds<2)) { 
      DETATCHSTATUSPTR (status);
      RETURN (status);
    }

  }

}

void LALRadstoRA(LALStatus *status, CHAR *RA, REAL8 rads) {

  INT4 hours, mins, secs, fracsecs;
  REAL8 x;

  INITSTATUS( status, "LALRadstoRA", rcsid );
  ATTATCHSTATUSPTR (status);

  if ((rads<0.0)||(rads>=LAL_TWOPI)) {
    fprintf(stderr, "alpha not in range [0,2PI)\n");
    exit(0);
  }

  x = rads*24.0/LAL_TWOPI;
  hours = (INT4)floor(x);
  x = 60.0*(x - (REAL8)hours);
  mins = (INT4)floor(x);
  x = 60.0*(x - (REAL8)mins);
  secs = (INT4)floor(x);
  x = x - (REAL8)secs;
  printf("x = %6.12f\n",x);
  fracsecs = (INT4)floor(1e8*x + 0.5);

  snprintf(RA,18,"%d:%02d:%02d.%08d",hours,mins,secs,fracsecs);
  printf("RA = %s\n",RA);
 
  DETATCHSTATUSPTR (status);
  RETURN (status);

}

void LALRadstoDEC(LALStatus *status, CHAR *DEC, REAL8 rads) {

  INT4 degs, arcmins, arcsecs, fracarcsecs;
  REAL8 x;
  INT4 sign = 1;

  INITSTATUS( status, "LALRadstoDEC", rcsid );
  ATTATCHSTATUSPTR (status);

  if ((rads<-(LAL_PI/2.0))||(rads>=(LAL_PI/2.0))) {
    fprintf(stderr, "delta not in range [-PI/2,PI/2)\n");
    exit(0);
  }

  if (rads<0.0) sign = -1;

  x = fabs(rads)*360.0/LAL_TWOPI;
  degs = (INT4)floor(x);
  x = 60.0*(x - (REAL8)degs);
  arcmins = (INT4)floor(x);
  x = 60.0*(x - (REAL8)arcmins);
  arcsecs = (INT4)floor(x);
  x = x - (REAL8)arcsecs;
  printf("x = %6.12f\n",x);
  fracarcsecs = (INT4)floor(1e7*x + 0.5);

  

  snprintf(DEC,18,"%d:%02d:%02d.%07d",sign*degs,arcmins,arcsecs,fracarcsecs);
  printf("DEC = %s\n",DEC);

  DETATCHSTATUSPTR (status);
  RETURN (status);

}

void MultiplyInterval(LALStatus *status,LALTimeInterval *newinterval,LALTimeInterval interval,INT4 factor) {

  INT4 temp;

  INITSTATUS( status, "MultiplyInterval", rcsid );
  ATTATCHSTATUSPTR (status);

  newinterval->seconds = interval.seconds*factor; 
  temp = interval.nanoSeconds*factor;
  newinterval->seconds += (INT4)floor(temp*1e-9);
  newinterval->nanoSeconds = (INT4)fmod(temp,1e9);

  DETATCHSTATUSPTR (status);
  RETURN (status);

}
