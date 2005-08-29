/********************************
Author: Shawhan, P. S.
$Id$
**************************************************** */

/* <lalVerbatim> */
#include <stdlib.h>
#include <lal/LALStdlib.h>
#include <lal/StringInput.h>

NRCSID( STRINGCONVERTTESTC, "$Id$" );

#define LAL_INT4_MAX 2147483647
#define LAL_INT4_ABSMIN LAL_UINT8_C(2147483648)

#define MAXGPSSTRINGS 256
#define SETGPSCASE( string, sec, ns, remain )                    \
    istring++;                                                   \
    if (istring >= MAXGPSSTRINGS) {                              \
        fprintf(stderr,"Too many GPS test cases; ABORTING\n");   \
        exit -1;                                                 \
    }                                                            \
    strcpy( gpsString[istring], string );                       \
    gpsOutSec[istring] = sec ;                                   \
    gpsOutNS[istring] = ns ;                                     \
    strcpy( gpsOutRemainder[istring], remain );

int lalDebugLevel = 0;

int
main( int argc, char **argv )
{
  static LALStatus stat;
  CHAR gpsString[MAXGPSSTRINGS][256];
  INT4 gpsOutSec[MAXGPSSTRINGS];
  INT4 gpsOutNS[MAXGPSSTRINGS];
  CHAR gpsOutRemainder[MAXGPSSTRINGS][256];

  LIGOTimeGPS gps;
  CHAR *endptr;
  INT4 istring, remlength;
  INT4 nfailures = 0;

  /*------ Parse input line. ------*/
  if ( argc == 2 )
    lalDebugLevel = atoi( argv[1] );
  else if ( argc != 1 )
    {
      fprintf( stderr, "Usage: %s [ lalDebugLevel ]\n", argv[0] );
      return 0; /* so that test script won't fail */
    }

  /*------ Initialize arrays of example cases and expected outputs ------*/
  for ( istring=0; istring<MAXGPSSTRINGS; istring++ ) {
    gpsString[istring][0] = '\0';
    gpsOutSec[istring] = 0;
    gpsOutNS[istring] = 0;
    gpsOutRemainder[istring][0] = '\0';
  }

  /*------ Fill arrays of example cases and expected outputs ------*/
  istring = -1;    /* Crucial for proper operation of SETGPSCASE macro */
  SETGPSCASE( "1234.5", 1234, 500000000, "" );
  SETGPSCASE( "712345678", 712345678, 0, "" );
  SETGPSCASE( "00000000712346678", 712346678, 0, "" );
  SETGPSCASE( "000000000000000000000000000000000712347678", 712347678, 0, "" );
  SETGPSCASE( "000000000000000000712348678.00000000000000", 712348678, 0, "" );
  SETGPSCASE( "000000000000000000712349678.00000000000001", 712349678, 0, "" );
  SETGPSCASE( "722345678.", 722345678, 0, "" );
  SETGPSCASE( "1722346678.", 1722346678, 0, "" );
  SETGPSCASE( "01722347678.", 1722347678, 0, "" );
  SETGPSCASE( "001722348678.", 1722348678, 0, "" );
  SETGPSCASE( "732345678.0", 732345678, 0, "" );
  SETGPSCASE( "742345678.7", 742345678, 700000000, "" );
  SETGPSCASE( "752345678.000861", 752345678, 861000, "" );
  SETGPSCASE( "762345678.000862547", 762345678, 862547, "" );
  SETGPSCASE( "772345678.0008635474", 772345678, 863547, "" );
  SETGPSCASE( "782345678.0008645475", 782345678, 864548, "" );
  SETGPSCASE( "792345678.000865547687287", 792345678, 865548, "" );
  SETGPSCASE( "702345678.9999999994", 702345678, 999999999, "" );
  SETGPSCASE( "712345678.9999999995", 712345679, 0, "" );
  SETGPSCASE( "722345678.9999999996", 722345679, 0, "" );
  SETGPSCASE( "2000000000", 2000000000, 0, "" );
  SETGPSCASE( "7323456785", LAL_INT4_MAX, 999999999, "" );
  SETGPSCASE( "7423456785234", LAL_INT4_MAX, 999999999, "" );
  SETGPSCASE( "752345678e0", 752345678, 0, "" );
  SETGPSCASE( "762345678e+0", 762345678, 0, "" );
  SETGPSCASE( "772345678e-0", 772345678, 0, "" );
  SETGPSCASE( "782345678e00", 782345678, 0, "" );
  SETGPSCASE( "792345678e+00", 792345678, 0, "" );
  SETGPSCASE( "702345678e-00", 702345678, 0, "" );
  SETGPSCASE( "712345678.e0", 712345678, 0, "" );
  SETGPSCASE( "722345678.e+0", 722345678, 0, "" );
  SETGPSCASE( "732345678.e-0", 732345678, 0, "" );
  SETGPSCASE( "742345678.00e0", 742345678, 0, "" );
  SETGPSCASE( "752345678.00e+0", 752345678, 0, "" );
  SETGPSCASE( "762345678.00e-0", 762345678, 0, "" );
  SETGPSCASE( "772345678.06e0", 772345678, 60000000, "" );
  SETGPSCASE( "782345678.06e+0", 782345678, 60000000, "" );
  SETGPSCASE( "792345678.06e-0", 792345678, 60000000, "" );
  SETGPSCASE( "7023.45678e5", 702345678, 0, "" );
  SETGPSCASE( "7123.457785255e+05", 712345778, 525500000, "" );
  SETGPSCASE( "7223458785255e-4", 722345878, 525500000, "" );
  SETGPSCASE( "43d", 43, 0, "d" );
  SETGPSCASE( "44.3873qr", 44, 387300000, "qr" );
  SETGPSCASE( "45.3973 qr", 45, 397300000, " qr" );
  SETGPSCASE( "46.3073 e2", 46, 307300000, " e2" );
  SETGPSCASE( "47.3173e2", 4731, 730000000, "" );
  SETGPSCASE( "6.85e7", 68500000, 0, "" );
  SETGPSCASE( "6.9512345678901e7", 69512345, 678901000, "" );
  SETGPSCASE( "6.05e7dkjf", 60500000, 0, "dkjf" );
  SETGPSCASE( "6.15ex0", 6, 150000000, "ex0" );
  SETGPSCASE( "6.25E7", 62500000, 0, "" );
  SETGPSCASE( "6.35E7dkjf", 63500000, 0, "dkjf" );
  SETGPSCASE( "6.45Ex0", 6, 450000000, "Ex0" );
  SETGPSCASE( "752345678.5433e258", LAL_INT4_MAX, 999999999, "" );
  SETGPSCASE( "762345678.5533e258r574", LAL_INT4_MAX, 999999999, "r574" );
  SETGPSCASE( "772345678.5633e.258", 772345678, 563300000, "e.258" );
  SETGPSCASE( "782345678.5733.258", 782345678, 573300000, ".258" );
  SETGPSCASE( "792345678.5833+258", 792345678, 583300000, "+258" );
  SETGPSCASE( "702345678.5933-258", 702345678, 593300000, "-258" );
  SETGPSCASE( "712345678.5033.258E02", 712345678, 503300000, ".258E02" );
  SETGPSCASE( "-722345678.5133", -722345679, 486700000, "" );
  SETGPSCASE( "-73234567800.5233", (INT4) (-LAL_INT4_ABSMIN), 0, "" );
  SETGPSCASE( "-742345678.000000625", -742345679, 999999375, "" );
  SETGPSCASE( "-743345678.9999999994", -743345679, 1, "" );
  SETGPSCASE( "-744345678.9999999995", -744345679, 0, "" );
  SETGPSCASE( "-752345678.9999999996", -752345679, 0, "" );
  SETGPSCASE( "5e-2", 0, 50000000, "" );
  SETGPSCASE( "7e-7", 0, 700, "" );
  SETGPSCASE( "6e-10", 0, 1, "" );
  SETGPSCASE( "8e-11", 0, 0, "" );
  SETGPSCASE( "-7e-12", 0, 0, "" );
  SETGPSCASE( "-4e-6", -1, 999996000, "" );
  SETGPSCASE( "-4.2e-2", -1, 958000000, "" );
  SETGPSCASE( ".5244", 0, 524400000, "" );
  SETGPSCASE( "-.5244", -1, 475600000, "" );
  SETGPSCASE( "0", 0, 0, "" );
  SETGPSCASE( "+", 0, 0, "+" );
  SETGPSCASE( "-", 0, 0, "-" );
  SETGPSCASE( "e", 0, 0, "e" );
  SETGPSCASE( "e3", 0, 0, "e3" );
  SETGPSCASE( "x", 0, 0, "x" );

  /*------ Loop over GPS test strings ------*/
  for ( istring=0; istring<MAXGPSSTRINGS; istring++ ) {
    if ( gpsString[istring][0] == '\0' ) continue;

    LALStringToGPS( &stat, &gps, gpsString[istring], &endptr );

    /* Check for a LAL error condition */
    if ( stat.statusCode ) {
      if ( lalDebugLevel > 0 ) {
	fprintf( stderr,
		 "Error[0] 1: program %s, file %s, line %i, %s\n"
		 "         Function LALStringToGPS() failed\n",
		 argv[0], __FILE__, __LINE__, STRINGCONVERTTESTC );
	REPORTSTATUS( &stat );
	return stat.statusCode;
      }
    }

    /* Print result */
    if ( !stat.statusCode ) {
      fprintf( stdout, "For string '%s':\n", gpsString[istring] );
      fprintf( stdout,
	       "  GPS time (sec,nsec) = (%11d,%10d), remainder '%s' ",
	       gps.gpsSeconds, gps.gpsNanoSeconds, endptr );
    }
    /* Pad with spaces to line up PASS/FAIL text */
    remlength = strlen(endptr);
    while ( remlength < 10 ) { fprintf(stdout," "); remlength++; }

    /* Check correctness of result */
    if ( gps.gpsSeconds == gpsOutSec[istring] &&
	 gps.gpsNanoSeconds == gpsOutNS[istring] &&
	 strcmp( endptr, gpsOutRemainder[istring] ) == 0 ) {
      fprintf( stdout, "Pass\n" );
    } else {
      fprintf( stdout, "*Fail*\n" );
      nfailures++;
    }


  }

  /*-- Report if there were any failures --*/
  if ( nfailures ) {
    fprintf( stdout, "Summary of GPS string conversion tests: %d FAILURES\n",
	     nfailures );
    return 9;
  } else {
    fprintf(stdout,"Summary of GPS string conversion tests: all succeeded\n");
  }

  return 0;
}
/* </lalVerbatim> */
