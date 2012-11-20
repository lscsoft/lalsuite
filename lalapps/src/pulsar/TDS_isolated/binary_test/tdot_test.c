/* 
 * Matt Pitkin 31/10/2012
 * 
 * test that the time derivative from the new Barycentering function is 
 * still correct 
 * 
 * compile with:
 * gcc tdot_test.c -o tdot_test -L${LSCSOFT_LOCATION}/lib
-I${LSCSOFT_LOCATION}/include -lm -llalsupport -llal -llalpulsar -std=c99
 * 
 * Change ttype to TYPE_TDB for TDB test, TYPE_TCB for TCB test or
 * TYPE_ORIGINAL to test original XLALBarycenterEarth function 
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <lal/LALBarycenter.h>
#include <lal/LALInitBarycenter.h>
#include <lal/Date.h>
#include <lal/LALConstants.h>
#include <lal/LALString.h>

int main (int argv, char **argc){
  FILE *fpout = NULL;
  
  BarycenterInput baryinput;
  EarthState earth;
  EmissionTime emit;
  EphemerisData *edat=NULL;
  TimeCorrectionData *tdat=NULL;
  char *earthFile=NULL, *sunFile=NULL, *tcFile=NULL, *lalpath=NULL;

  TimeCorrectionType ttype = TYPE_TDB;
 
  double tstart = 900000001., dt = 600., dstep = 1.;
  int npoints = 50000, i = 0;
  
  /* location of the Parkes telescope */
  baryinput.site.location[0] = -4554231.5/LAL_C_SI;
  baryinput.site.location[1] = 2816759.1/LAL_C_SI;
  baryinput.site.location[2] = -3454036.3/LAL_C_SI;
  
  if((lalpath = getenv("LALPULSAR_PREFIX")) == NULL){
    fprintf(stderr, "LALPULSAR_PREFIX environment variable not set!\n");
    exit(1);
  }
  
  earthFile = XLALStringDuplicate(lalpath);
  sunFile = XLALStringDuplicate(lalpath);
  
  earthFile = XLALStringAppend( earthFile, 
                                "/share/lalpulsar/earth00-19-DE405.dat.gz");
  sunFile = XLALStringAppend( sunFile,
                              "/share/lalpulsar/sun00-19-DE405.dat.gz");

  edat = XLALInitBarycenter( earthFile, sunFile );

  fpout = fopen("tdot.txt", "w");

  tcFile = XLALStringDuplicate(lalpath);
  
  /* read in the time correction file */
  if( ttype == TYPE_TEMPO2 || ttype == TYPE_TCB ){
    tcFile = XLALStringAppend( tcFile,
                               "/share/lalpulsar/te405_2000-2019.dat.gz" );
  }
  else if ( ttype == TYPE_TDB ){
    tcFile = XLALStringAppend( tcFile,
                               "/share/lalpulsar/tdb_2000-2019.dat.gz" ); 
  }
  
  tdat = XLALInitTimeCorrections( tcFile );

  baryinput.delta = LAL_PI*(66-90)/180.; /* declination (on ecliptic) */
  baryinput.alpha = LAL_TWOPI*18./24.; /* right ascension */
  baryinput.dInv = 0.0;  /* no parallax */
  
  for( i=0; i<npoints; i++ ){
    double t = tstart + (double)i*dt;
    double tplus, tminus;
    
    /* get time derivative at time t */
    XLALGPSSetREAL8( &baryinput.tgps, t );
    XLALBarycenterEarthNew( &earth, &baryinput.tgps, edat, tdat, ttype );
    XLALBarycenter( &emit, &baryinput, &earth );
    
    fprintf(fpout, "%.5lf\t%.12lf\t", t, emit.tDot);
    
    /* get time derivative by finding the gradient of the time delay manually */
    XLALGPSSetREAL8( &baryinput.tgps, t+dstep/2. );
    XLALBarycenterEarthNew( &earth, &baryinput.tgps, edat, tdat, ttype );
    XLALBarycenter( &emit, &baryinput, &earth );
    tplus = emit.deltaT;
    
    XLALGPSSetREAL8( &baryinput.tgps, t-dstep/2. );
    XLALBarycenterEarthNew( &earth, &baryinput.tgps, edat, tdat, ttype );
    XLALBarycenter( &emit, &baryinput, &earth );
    tminus = emit.deltaT;
    
    fprintf(fpout, "%.12lf\n", 1. + (tplus-tminus)/dstep);
  }
  
  fclose(fpout);
  
  return 0;
}
