/****************************** <lalVerbatim file="LALBarycenterTestCV">
Author: Cutler, C.
$Id$
******************************* </lalVerbatim> */

/* <lalLaTeX>

\subsection{Program \texttt{LALBarycenterTest.c}}
\label{ss:LALBarycenterTest.c}

Tests the routine \verb@LALBarycenter()@.  Exercises some of the error
conditions and makes sure that they work.

\subsubsection*{Usage}
\begin{verbatim}
LALBarycenterTest 
\end{verbatim}

\subsubsection*{Description}

This program demonstrates the use of \verb@LALBarycenter.c@.
The two ephemeris files specified in the \verb@EphemerisFilenames@
structure (e.g., for data taken in 1998, \verb@sun98.dat@ and \verb@earth98.dat@)
are assumed to be in the same directory as the program as
the test program. 

\subsubsection*{Exit codes}
\input{LALBarycenterTestCE}

\subsubsection*{Uses}
\begin{verbatim}
lalDebugLevel
LALMalloc()
LALFree()
LALBarycenterInit()
LALBarycenterEarth()
LALBarycenter()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALBarycenterTestCV}}

</lalLaTeX> */

#include <lal/LALBarycenter.h>
#include "LALInitBarycenter.h"

NRCSID(LALBARYCENTERTESTC,"$Id$");

/***************************** <lalErrTable file="LALBarycenterTestCE"> */
#define LALBARYCENTERTESTC_ENOM 0
#define LALBARYCENTERTESTC_EARG 1
#define LALBARYCENTERTESTC_ESUB 2

#define LALBARYCENTERTESTC_MSGENOM "Nominal exit"
#define LALBARYCENTERTESTC_MSGEARG "Error parsing command-line arguments"
#define LALBARYCENTERTESTC_MSGESUB "Subroutine returned error"
/***************************** </lalErrTable> */

/*
  int lalDebugLevel=0;
*/
  INT4 lalDebugLevel=7;

BarycenterInput baryinput;
EmissionTime  emit;  
LIGOTimeGPS tGPS;
EarthState earth;

INT4 t2000 = 630720013; /* gps time at Jan 1, 2000 00:00:00 UTC */
INT4 t1998 = 630720013-730*86400-1;/* gps at Jan 1,1998 00:00:00 UTC*/

int
main()
{
  static LALStatus stat;
  
  INT4 i,k; /*dummy indices*/
  EphemerisData *edat = NULL;
  
  
  REAL8 alpha,delta;  /* RA and DEC (radians) in 
			 ICRS realization of J2000 coords.*/
  
  edat = (EphemerisData *)LALMalloc(sizeof(EphemerisData));    
  
  (*edat).ephiles.earthEphemeris = "earth98.dat";
  (*edat).ephiles.sunEphemeris = "sun98.dat";
  
  
  LALInitBarycenter(&stat, edat);
  printf("stat.statusCode = %d\n",stat.statusCode); 
  REPORTSTATUS(&stat);

 /* The routines using LALBarycenter package, the code above, leading 
    up LALInitBarycenter call, should be near top of main. The idea is
    that ephemeris data is read into RAM once, at the beginning.
 */   

  
  { /*Now getting coords. for GEO detector. Could eventually be done
      with a call to LALCreateDetector */

  REAL8 flatFac = 0.00335281e0; /*flattening factor for ellipsoidal Earth,
                                 from Explan. Supp. to Astronom. Almanac */ 
  
  REAL8 longitude,latitude,elev,theta,rd;
  
  /*Note I'm assuming values below are GEO's long and lat in geocentric
    coords--BUT THAT may not be right!! E.g., in obsys.dat, tempo takes
    uses long and lat in geodetic coords (see p. 703 of Supp to Astron Almanac)
    Also I set elevation to zero; MUST find and right elevation!! 
  */
  
  latitude = (52.e0 + 15.e0/60)*LAL_PI/180; /*detector latitude for GEO 
					      (radians north of equator) */  
  
  longitude = (9.e0 + 48.e0/60 + 36.e0/3600)*LAL_PI/180; 
                 /*detector longitude for GEO (radians east of Greenwich) */
  
  elev =  0.e0;  /*detector elev above mean sea level (sec) */
  
  /*-------------------------------------------------------------------------
   *converting latitude from geodetic to geocentric (from p. 700 of Ex. Supp.)
   * if necessary:  
   
   latitude = latitude - (692.74e0*LAL_PI/(3600*180))*sin(2.e0*latitutde)
   + (1.16e0*LAL_PI/(3600*180))*sin(4.e0*latitude);
   *---------------------------------------------------------------------------
   */
  
  theta=LAL_PI/2.e0 - latitude;
  
  rd=elev + (LAL_REARTH_SI/LAL_C_SI)*(1.e0-flatFac)
    /sqrt(1.e0 - (sin(theta))*(sin(theta))*flatFac*(2.e0-flatFac) );
  
  baryinput.site.location[0]=rd*sin(theta)*cos(longitude);
  baryinput.site.location[1]=rd*sin(theta)*sin(longitude);
  baryinput.site.location[2]=rd*cos(theta);
  }  
  
  /* next: outer loop over pulse arrival times; LALBarycenterEarth
    called ONCE per arrival time */
  
  for (i=0;i < 10; i++){
    
    /*GPS time(sec) =  tGPS.gpsSeconds + 1.e-9*tGPS.gpsNanoSeconds  */ 
    
    tGPS.gpsSeconds = t1998;
    tGPS.gpsSeconds +=i*3600*50; 
    tGPS.gpsNanoSeconds = 0;
    
    LALBarycenterEarth(&stat, &earth, &tGPS, edat);
    REPORTSTATUS(&stat);
    
    /*next: inner loop over different sky positions, for each arrival time;
     LALBarycenter called ONCE per sky position (or ONCE per detector) */
    
    for (k=0;k<3;k++){
      
      alpha=(LAL_PI/12.0)*(14.e0 + 51.e0/60.e0 + 
			   +38.56024702e0/3.6e3) + LAL_PI*k/10.e0;
      delta=(LAL_PI/180.e0)*(12.e0+ 19.e0/60.e0
			     +59.1434800e0/3.6e3); 
      
      baryinput.alpha = alpha;
      baryinput.delta = delta;
      
      baryinput.tgps.gpsSeconds = tGPS.gpsSeconds;
      baryinput.tgps.gpsNanoSeconds = tGPS.gpsNanoSeconds;
      
      LALBarycenter(&stat, &emit, &baryinput, &earth);
      REPORTSTATUS(&stat);
      
      printf("%d %d %d %25.17e %25.17e\n", k, 
	     tGPS.gpsSeconds,  tGPS.gpsNanoSeconds,
	     (emit.deltaT + tGPS.gpsSeconds + tGPS.gpsNanoSeconds*1.e-9), 
             emit.tDot);

      printf("%d %d %25.17e\n",
	     emit.te.gpsSeconds, emit.te.gpsNanoSeconds, emit.deltaT);

      printf("%25.17e %25.17e %25.17e\n",
	     emit.rDetector[0],emit.rDetector[1],emit.rDetector[2]);

      printf("%25.17e %25.17e %25.17e\n",
	     emit.vDetector[0],emit.vDetector[1],emit.vDetector[2]);
    }    
  }
  LALFree(edat);
  return 0;
}
