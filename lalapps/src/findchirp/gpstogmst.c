#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/Date.h>

double greenwich_mean_sidereal_time( int gpssec, int gpsnan, int taiutc )
{
  /* cf. S. Aoki et al., A&A 105, 359 (1982) eqs. 13 & 19 */
  /* also cf. http://aa.usno.navy.mil */
  /* Note: 00h UT 01 Jan 2000 has JD=2451544.5 and GPS=630720013 */
  const double JD_12h_01_Jan_2000     = 2451545.0;
  const double JD_00h_01_Jan_2000     = 2451544.5;
  const double GPS_00h_01_Jan_2000    = 630720013;
  const double TAIUTC_00h_01_Jan_2000 = 32; /* leap seconds: TAI - UTC */

  double t;
  double dpU;
  double TpU;
  double gmst;

  /* compute number of seconds since 00h UT 01 Jan 2000 */
  t  = gpssec - GPS_00h_01_Jan_2000;
  t += 1e-9 * gpsnan;
  t += taiutc - TAIUTC_00h_01_Jan_2000;

  /* compute number of days since 12h UT 01 Jan 2000 */
  dpU  = floor( t / ( 24.0 * 3600.0 ) ); /* full days since 0h UT 01 Jan 2000 */
  dpU += JD_00h_01_Jan_2000 - JD_12h_01_Jan_2000; /* i.e., -0.5 */

  /* compute number of centuries since 12h UT 31 Dec 1899 */
  TpU = dpU / 36525.0;

  /* compute the gmst at 0h of the current day */
  gmst = 24110.54841
    + TpU * ( 8640184.812866
        + TpU * ( 0.093104
          - TpU * 6.2e-6 ) ); /* seconds */

  /* add the sidereal time since the start of the day */
  t = fmod( t, 24.0 * 3600.0 ); /* seconds since start of day */
  gmst += t * 1.002737909350795; /* corrections omitted */

  /* convert to fractions of a day and to radians */
  gmst = fmod( gmst / ( 24.0 * 3600.0 ), 1.0 ); /* fraction of day */
  gmst *= 2.0 * LAL_PI; /* radians */
  return gmst;
}

int main ( int argc, char *argv[] )
{
  LALStatus status;
  LALMSTUnitsAndAcc units_and_acc;
  LIGOTimeGPS gps;
  double j_gmst, d_gmst;

  if ( argc != 3 )
  {
    fprintf( stderr, "Usage: %s sec ns\n", argv[0] );
    return 1;
  }
  
  memset( &status, 0, sizeof(LALStatus) );
  
  gps.gpsSeconds = atoi( argv[1] );
  gps.gpsNanoSeconds = atoi( argv[2] );
  fprintf( stdout, "GPS time = %d sec %d ns\n\n", 
      gps.gpsSeconds, gps.gpsNanoSeconds );

  j_gmst = 
    greenwich_mean_sidereal_time( gps.gpsSeconds, gps.gpsNanoSeconds, 32 );

  units_and_acc.units = MST_RAD;
  units_and_acc.accuracy = LALLEAPSEC_STRICT;
  LALGPStoGMST1(&status, &d_gmst, &gps, &units_and_acc);

  fprintf( stdout, "inspinj GMST = %22.16e rads\n", j_gmst );
  fprintf( stdout, "LAL GMST     = %22.16e rads\n\n", d_gmst );

  j_gmst = j_gmst * 12.0 / LAL_PI;

  units_and_acc.units = MST_HRS;
  LALGPStoGMST1(&status, &d_gmst, &gps, &units_and_acc);

  fprintf( stdout, "inspinj GMST = %22.16e h\n", j_gmst );
  fprintf( stdout, "LAL GMST     = %22.16e h\n", d_gmst );

  return 0;
}
