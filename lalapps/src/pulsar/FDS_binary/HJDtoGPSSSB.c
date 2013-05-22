/*
 * Copyright (C) 2010 Chris Messenger
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

/** \author C.Messenger
 * \ingroup pulsarApps
 * \file
 * \brief
 * This code is converts an HJD (Heliocentric Julian Date) to a GPS time defined at the
 * SSB (solar system barycenter).  
 *
 * It simply takes the dot product of the vector defining the sun location in the SSB 
 * frame with the unit vector defining the sky position of the source.  The quantity 
 * (divided by c) is the time difference between a wave-front passing the heliocenter
 * and the SSB i.e the SBS time is equal to the heliocentric time *plus* this correction.
 * 
 * The accuracy of this code is only good to of order 0.01 seconds.
 *
 */

/***********************************************************************************************/
/* includes */
#include <lal/LALDatatypes.h>
#include <lal/UserInput.h>
#include <lal/LogPrintf.h>
#include <lalapps.h> 
#include <lal/LALInitBarycenter.h> 
#include <lal/FindRoot.h>
#include <lal/BinaryPulsarTiming.h>

/* gsl includes */ 
#include <gsl/gsl_spline.h>

/***********************************************************************************************/
/* some global constants */

#define STRINGLENGTH 256              /* the length of general string */
#define LONGSTRINGLENGTH 1024         /* the length of general string */
 
/***********************************************************************************************/
/* some useful macros */

#define MYMAX(x,y) ( (x) > (y) ? (x) : (y) )
#define MYMIN(x,y) ( (x) < (y) ? (x) : (y) )

/***********************************************************************************************/
/* define internal structures */

/** A structure that stores user input variables 
 */
typedef struct { 
  BOOLEAN help;		            /**< trigger output of help string */
  REAL8 HJD;                        /**< the heliocentric julian date */
  CHAR *HJDconv;                    /**< the timing convention of the input HJD time */
  CHAR *ra;                         /**< the right ascension in hh:mm:ss.s format */
  CHAR *dec;                        /**< the declination in dd:mm:ss.s format */
  REAL8 ra_rads;                    /**< the right ascension in radians */
  REAL8 dec_rads;                   /**< the declination in radians */
  CHAR *ephemdir;                   /**< the ephemeris directory */
  CHAR *ephemyear;                  /**< the ephemeris year */
  BOOLEAN version;	            /**< output version-info */
} UserInput_t;

/***********************************************************************************************/
/* global variables */
extern int vrbflg;	 	/**< defined in lalapps.c */

/* gsl stuff needs to be seen globally for the root finding */
gsl_interp_accel *acc_x = NULL;               /* gsl interpolation accelerator */
gsl_interp_accel *acc_y = NULL;               /* gsl interpolation accelerator */
gsl_interp_accel *acc_z = NULL;               /* gsl interpolation accelerator */
gsl_spline *spline_x = NULL;                  /* gsl spline variable */
gsl_spline *spline_y = NULL;                  /* gsl spline variable */
gsl_spline *spline_z = NULL;                  /* gsl spline variable */
REAL8 n[3];                                   /* the source position unit vector */
REAL8 *t = NULL;
REAL8 *x = NULL;
REAL8 *y = NULL;
REAL8 *z = NULL;

/***********************************************************************************************/
/* define functions */
int main(int argc,char *argv[]);
int XLALReadUserVars(int argc,char *argv[],UserInput_t *uvar,CHAR **clargs);
static void TimeOfEmission(LALStatus *status,REAL8 *tr,REAL8 lE,void *tr0);
REAL8 LALUTCMJDtoGPS(REAL8 MJDutc);

/***********************************************************************************************/
/* empty initializers */
UserInput_t empty_UserInput;

/** The main function of semicoherentbinary.c
 *
 */
int main( int argc, char *argv[] )
{
  LALStatus status = blank_status;              /* empty LAL status structure */
  UserInput_t uvar = empty_UserInput;           /* user input variables */
  CHAR *clargs = NULL;                          /* store the command line args */
  UINT4 i;                                      /* counter */
  EphemerisData *ephemeris = NULL;              /* the ephemeris data */
  CHAR earthfile[LONGSTRINGLENGTH];             /* name of earth ephemeris file */
  CHAR sunfile[LONGSTRINGLENGTH];               /* name of sun ephemeris file */ 
  DFindRootIn input;                            /* the input structure for the root finding procedure */
  REAL8 helio_gps;                              /* the gps version of the input HJD time */
  REAL8 te;                                     /* the GPS emission time at the earth */
  REAL8 alpha,delta;                            /* the sky position angles in radians */

  vrbflg = 1;	                        /* verbose error-messages */

  /* turn off default GSL error handler */
  gsl_set_error_handler_off();

  /* setup LAL debug level */
  LogSetLevel(lalDebugLevel);

  /* register and read all user-variables */
  if (XLALReadUserVars(argc,argv,&uvar,&clargs)) {
    LogPrintf(LOG_CRITICAL,"%s : XLALReadUserVars() failed with error = %d\n",__func__,xlalErrno);
    return 1;
  }
  LogPrintf(LOG_DEBUG,"%s : read in uservars\n",__func__);
 
  /* if coordinates input in hh:mm:ss.s format then convert to radians */
  if (XLALUserVarWasSet(&(uvar.ra))) alpha = XLALhmsToRads(uvar.ra);
  else alpha = uvar.ra_rads;
  if (XLALUserVarWasSet(&(uvar.dec))) delta = XLALdmsToRads(uvar.dec);
  else delta = uvar.dec_rads;

  {
    REAL8 sinTheta = sin(LAL_PI/2.0-delta);
    n[2]=cos(LAL_PI/2.0-delta);   /* n is vector that points towards source */
    n[1]=sinTheta*sin(alpha);     /* in Cartesian coords based on J2000 */
    n[0]=sinTheta*cos(alpha);     /* 0=x,1=y,2=z */
  }

  /* convert MJD to TT GPS time */
  if (strstr(uvar.HJDconv,"TT")) {
    helio_gps = XLALTTMJDtoGPS(uvar.HJD-2400000.5);
    LogPrintf(LOG_DEBUG,"%s : heliocentric MJD %6.12f -> GPS(TT) %6.12f\n",__func__,uvar.HJD,helio_gps);
  }
  else if (strstr(uvar.HJDconv,"TDB")) {
    helio_gps = XLALTDBMJDtoGPS(uvar.HJD-2400000.5);
    LogPrintf(LOG_DEBUG,"%s : heliocentric MJD %6.12f -> GPS(TDB) %6.12f\n",__func__,uvar.HJD,helio_gps);
  }
  else if (strstr(uvar.HJDconv,"UTC")) {
    helio_gps = LALUTCMJDtoGPS(uvar.HJD-2400000.5);
    LogPrintf(LOG_DEBUG,"%s : heliocentric MJD %6.12f -> GPS(UTC) %6.12f\n",__func__,uvar.HJD,helio_gps);
  }
  else {
    LogPrintf(LOG_CRITICAL,"%s : Can only convert HJD times in TT, TDB and UTC at present.  Exiting.\n",__func__);
    return 1;
  }

  /* define names of ephemeris files */
  snprintf(earthfile,LONGSTRINGLENGTH,"%s/earth%s.dat",uvar.ephemdir,uvar.ephemyear);
  snprintf(sunfile,LONGSTRINGLENGTH,"%s/sun%s.dat",uvar.ephemdir,uvar.ephemyear);

  /* initialise the barycentering routines */
  if ((ephemeris = XLALInitBarycenter(earthfile,sunfile)) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : XLALInitBaryCenter() failed with error = %d\n",__func__,xlalErrno);
    return 1;
  }
  LogPrintf(LOG_DEBUG,"%s : read in ephemeris files\n",__func__);

  /* allocate memory for gsl input vectors */
  if ((t = XLALCalloc(ephemeris->nentriesS,sizeof(REAL8))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : XLALCalloc() failed with error = %d\n",__func__,xlalErrno);
    return 1;
  }
  if ((x = XLALCalloc(ephemeris->nentriesS,sizeof(REAL8))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : XLALCalloc() failed with error = %d\n",__func__,xlalErrno);
    return 1;
  }
  if ((y = XLALCalloc(ephemeris->nentriesS,sizeof(REAL8))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : XLALCalloc() failed with error = %d\n",__func__,xlalErrno);
    return 1;
  }
  if ((z = XLALCalloc(ephemeris->nentriesS,sizeof(REAL8))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : XLALCalloc() failed with error = %d\n",__func__,xlalErrno);
    return 1;
  }

  /* setup gsl interpolation for each of the 3 earth position components */
  for (i=0;i<(UINT4)ephemeris->nentriesS;i++) {
    t[i] = ephemeris->ephemS[i].gps;
    x[i] = ephemeris->ephemS[i].pos[0];
    y[i] = ephemeris->ephemS[i].pos[1];
    z[i] = ephemeris->ephemS[i].pos[2];
  }  
  LogPrintf(LOG_DEBUG,"%s : time spans %6.12f -> %6.12f\n",__func__,t[0],t[ephemeris->nentriesS-1]);
  acc_x = gsl_interp_accel_alloc();
  acc_y = gsl_interp_accel_alloc();
  acc_z = gsl_interp_accel_alloc();
  spline_x = gsl_spline_alloc(gsl_interp_cspline,ephemeris->nentriesS);
  spline_y = gsl_spline_alloc(gsl_interp_cspline,ephemeris->nentriesS);
  spline_z = gsl_spline_alloc(gsl_interp_cspline,ephemeris->nentriesS);
  gsl_spline_init(spline_x,t,x,ephemeris->nentriesS);
  gsl_spline_init(spline_y,t,y,ephemeris->nentriesS);
  gsl_spline_init(spline_z,t,z,ephemeris->nentriesS);
     
  /* compute gps emission time using a root finding procedure */
  input.function = TimeOfEmission;     /* This is the name of the function we must solve to find the time of emission */
  input.xmin = helio_gps - 500;        /* We know that E will be found within the possible light travel time window */
  input.xmax = helio_gps + 500;
  input.xacc = 1e-6;                   /* The accuracy of the root finding procedure */
  
  /* expand domain until a root is bracketed */
  LALDBracketRoot(&status,&input,&helio_gps);
  
  /* bisect domain to find eccentric anomoly E corresponding to the SSB time of the midpoint of this SFT */
  LALDBisectionFindRoot(&status,&te,&input,&helio_gps);

  /* check if result is within original ephemeris boundaries */
  if ((te>=t[0]) && (te<=t[ephemeris->nentriesS-1])) {
    LogPrintf(LOG_NORMAL,"%s : solved to find SSB time = %6.2f +/- 0.01 (sec)\n",__func__,te); 
  }
  else {
    LogPrintf(LOG_CRITICAL,"%s : input ephemeris files do not span desired epoch.  Exptroplated result is unreliable.\n",__func__,te); 
    return 1;
  }

  /* free ephemeris info */
  XLALDestroyEphemerisData(ephemeris);
 
  /* free gsl interpolation */
  gsl_spline_free(spline_x);
  gsl_spline_free(spline_y);
  gsl_spline_free(spline_z);
  gsl_interp_accel_free(acc_x);
  gsl_interp_accel_free(acc_y);
  gsl_interp_accel_free(acc_z);
  XLALFree(t);
  XLALFree(x);
  XLALFree(y);
  XLALFree(z);
   
 /* Free config-Variables and userInput stuff */
  XLALDestroyUserVars();
  XLALFree(clargs);

  /* did we forget anything ? */
  LALCheckMemoryLeaks();
  LogPrintf(LOG_DEBUG,"%s : successfully checked memory leaks.\n",__func__);

  LogPrintf(LOG_DEBUG,"%s : successfully completed.\n",__func__);
  return 0;
  
} /* end of main */

/** For a given set of binary parameters we solve the following function for
 *  the eccentric anomoly E
 */
static void TimeOfEmission(LALStatus *status,
			   REAL8 *tr,
			   REAL8 tssb,
			   void *tr0
			   )
{
  INITSTATUS(status);
  ASSERT(tr0,status, 1, "Null pointer");

  /* interpolate to find the x,y,z coordinates of the heliocenter in the SSB frame (in seconds) */
  REAL8 xr = gsl_spline_eval(spline_x,tssb,acc_x);
  REAL8 yr = gsl_spline_eval(spline_y,tssb,acc_y);
  REAL8 zr = gsl_spline_eval(spline_z,tssb,acc_z);

  /* compute projection along the line of sight vector to the source */
  REAL8 helio_deltat = xr*n[0] + yr*n[1] + zr*n[2];

  /* this is the function relating the time of passage through the heliocenter to the passage through the SSB */
  *tr = *(REAL8 *)tr0*(-1.0) + tssb - helio_deltat;

  RETURN(status);

}

/** Read in input user arguments
 *
 */
int XLALReadUserVars(int argc,            /**< [in] the command line argument counter */ 
		     char *argv[],        /**< [in] the command line arguments */
		     UserInput_t *uvar,   /**< [out] the user input structure */
		     CHAR **clargs        /**< [out] the command line args string */
		     )
{
  CHAR *version_string;
  INT4 i;

 
  /* ---------- register all user-variables ---------- */
  XLALregBOOLUserStruct(help, 		        'h', UVAR_HELP,     "Print this message");
  XLALregREALUserStruct(HJD,                    'H', UVAR_REQUIRED, "The time in Heliocentric Julian Days");
  XLALregSTRINGUserStruct(HJDconv, 	        'c', UVAR_REQUIRED, "The timing convention for the HJD value (TT, TDB or UTC)");
  XLALregSTRINGUserStruct(ra, 	                'R', UVAR_OPTIONAL, "The source right ascension in hh:mm:ss.s format");
  XLALregSTRINGUserStruct(dec, 	                'D', UVAR_OPTIONAL, "The source declination in dd:mm:ss.s format");
  XLALregREALUserStruct(ra_rads,                'r', UVAR_OPTIONAL, "The source right ascension in radians");
  XLALregREALUserStruct(dec_rads,               'd', UVAR_OPTIONAL, "The source declination in radians");
  XLALregSTRINGUserStruct(ephemdir, 	        'E', UVAR_REQUIRED, "The ephemeris directory");
  XLALregSTRINGUserStruct(ephemyear, 	        'y', UVAR_REQUIRED, "The ephemeris year string"); 
  XLALregBOOLUserStruct(version,                'V', UVAR_SPECIAL,  "Output code version");

  /* do ALL cmdline and cfgfile handling */
  if (XLALUserVarReadAllInput(argc, argv)) {
    LogPrintf(LOG_CRITICAL,"%s : XLALUserVarReadAllInput() failed with error = %d\n",__func__,xlalErrno);
    XLAL_ERROR(XLAL_EINVAL);
  }
  
  /* if help was requested, we're done here */
  if (uvar->help) exit(0);

  if ((version_string = XLALGetVersionString(0)) == NULL) {
    XLALPrintError("XLALGetVersionString(0) failed.\n");
    exit(1);
  }
  
  if (uvar->version) {
    printf("%s\n",version_string);
    exit(0);
  }
  XLALFree(version_string);

  if ((!XLALUserVarWasSet(&(uvar->ra))) && (!XLALUserVarWasSet(&(uvar->ra_rads)))) {
    LogPrintf(LOG_CRITICAL,"%s : user must specify at least ra or ra_rads !  Exiting.\n",__func__,xlalErrno);
    XLAL_ERROR(XLAL_EINVAL);
  }
  if ((!XLALUserVarWasSet(&(uvar->dec))) && (!XLALUserVarWasSet(&(uvar->dec_rads)))) {
    LogPrintf(LOG_CRITICAL,"%s : user must specify at least dec or dec_rads !  Exiting.\n",__func__,xlalErrno);
    XLAL_ERROR(XLAL_EINVAL);
  }

  /* put clargs into string */
  *clargs = XLALCalloc(1,sizeof(CHAR));
  for (i=0;i<argc;i++) {
    INT4 len = 2 + strlen(argv[i]) + strlen(*clargs);
    *clargs = XLALRealloc(*clargs,len*sizeof(CHAR));
    strcat(*clargs,argv[i]);
    strcat(*clargs," ");
  }

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",__func__);
  return XLAL_SUCCESS;
  
}

/** Compute the GPS time from an input UTC time in MJD format 
 *
 * This is a copy of an octapps code that Reinhard prix originally
 * copied from a lalapps code that no longer exists.
 *
 */
REAL8 LALUTCMJDtoGPS(REAL8 MJDutc     /**< [in] the UTC time in MJD format */
		     )
{
 
  REAL8 REF_GPS_SECS = 793130413.0;
  REAL8 REF_MJD = 53423.75;
  
  REAL8 GPS = (-REF_MJD + MJDutc) * 86400.0 + REF_GPS_SECS;

  return GPS;
  
}

