#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <lal/UserInput.h>
#include <lal/LALStdlib.h>
#include <lal/Date.h>
#include <lal/AVFactories.h>

#define TRUE (1==1)
#define FALSE (1==0)


INT4 lalDebugLevel;

NRCSID (UTCTOGPSC, "$Id$");

int main(int argc, char *argv[])
{
  static LALStatus     status;
  LIGOTimeGPS          gpsTime;
  LIGOTimeGPS          refGPS;
  LALDate              utcDate;
  LALLeapSecAccuracy   accuracy = LALLEAPSEC_LOOSE;
  CHARVector          *timestamp = NULL;
  time_t               sec;
  BOOLEAN uvar_help;
  INT4 uvar_year, uvar_month, uvar_day, uvar_hour, uvar_min, uvar_sec;

  /*  set up the default parameters  */
  lalDebugLevel = 0;
  /* LALDebugLevel must be called before anything else */
  LALGetDebugLevel( &status, argc, argv, 'd');

  uvar_help = FALSE;
  uvar_year = 1995;
  uvar_month = 1;
  uvar_day = 1;
  uvar_hour = 0;
  uvar_min = 0;
  uvar_sec = 0;

  LALRegisterBOOLUserVar(   &status, "help",            'h', UVAR_HELP,     "Print this message",  &uvar_help );  
  LALRegisterINTUserVar(    &status, "year",            'y', UVAR_OPTIONAL, "Year (XXXX)",         &uvar_year );
  LALRegisterINTUserVar(    &status, "month",           'm', UVAR_OPTIONAL, "Month (XX)",          &uvar_month );
  LALRegisterINTUserVar(    &status, "day",             'd', UVAR_OPTIONAL, "day (XX)",            &uvar_day ) ;
  LALRegisterINTUserVar(    &status, "month",           'H', UVAR_OPTIONAL, "Hour (XX)",           &uvar_hour );
  LALRegisterINTUserVar(    &status, "month",           'M', UVAR_OPTIONAL, "Min (XX)",            &uvar_min ) ;
  LALRegisterINTUserVar(    &status, "month",           's', UVAR_OPTIONAL, "Sec (XX)",            &uvar_sec ) ;
  

  /* read all command line variables */
  LALUserVarReadAllInput(&status, argc, argv);

  /* exit if help was required */
  if (uvar_help)
    exit(0); 


  /*LALCHARCreateVector(&status, &timestamp, (UINT4)64);*/

  if (lalDebugLevel > 0)
    REPORTSTATUS(&status);
  
  /* Set the date */
  utcDate.unixDate.tm_year = uvar_year - 1900;
  utcDate.unixDate.tm_mon  = uvar_month;
  utcDate.unixDate.tm_mday = uvar_day;
  utcDate.unixDate.tm_hour = uvar_hour;
  utcDate.unixDate.tm_min  = uvar_min;
  utcDate.unixDate.tm_sec  = uvar_sec;
  utcDate.residualNanoSeconds = 0;

  LALUTCtoGPS(&status, &gpsTime, &utcDate, &accuracy);

  fprintf(stdout, "%d\n", gpsTime.gpsSeconds);

  LALDestroyUserVars(&status);   

  /* normal exit */
  return 0;

}
