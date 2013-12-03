/*
*  Copyright (C) 2007 Peter Shawhan
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

#include <stdlib.h>
#include <lal/LALStatusMacros.h>
#include <lal/SegmentsIO.h>

/**
 * \file
 * \ingroup SegmentsIO_h
 * \author Peter Shawhan
 *
 * \brief Tests the routines in \ref SegmentsIO.c.
 *
 * ### Usage ###
 *
 * \code
 * SegmentsIOTest [ lalDebugLevel ]
 * \endcode
 *
 * The default value of \c lalDebugLevel is 4.
 *
 * If the \c lalDebugLevel argument is omitted, the test program sets it
 * to 4 to turn on info messages only.
 * Note that this default value does not enable LAL/XLAL error messages,
 * since many of the tests intentionally create error conditions and verify that
 * the proper error codes are generated.  If you want to turn on the LAL/XLAL
 * error and warning messages, specify a \c lalDebugLevel value of 7,
 * or 23 if you also want informational messages related to memory checking.
 * If you specify 0, then no messages will be printed under any circumstances.
 * However, in all cases, the return status of the program will be 0 if all
 * tests passed, 1 if one or more tests failed.
 *
 * \note This test program does not currently do an exhaustive test of
 * functionality and failure modes; it is more like a starting point for
 * spot-checking by modifying, recompiling and running this test program
 * and inspecting the output.
 *
 * ### Exit codes ###
 *
 * <ul>
 * <li>0 if all tests passed.</li>
 *
 * <li>1 if one or more tests failed.</li>
 * </ul>
 *
 */

/** \cond DONT_DOXYGEN */

/*-- Macros for this test program --*/

#define RETPASS(testname,status) XLALPrintInfo( \
    "Pass return check for %s: return=%d, xlalErrno=%d\n", \
    testname, status, xlalErrno );

#define RETFAIL(testname,status) XLALPrintInfo( \
    "*FAIL* return check for %s: return=%d, xlalErrno=%d\n", \
    testname, status, xlalErrno ); nfailures++;

#define FUNCPASS(testname) XLALPrintInfo( \
    "Pass functional check for %s\n", testname );

#define FUNCFAIL(testname,reason) XLALPrintInfo( \
    "*FAIL* functional check for %s: %s\n", testname, reason ); nfailures++;

/*-- Default debug level includes info messages (4), but not
     memory checking (16), error messages (1), or warning messages (2) --*/

/*===========================================================================*/


/*===========================================================================*/
int main(void)
{
  INT4 nfailures = 0;
  static LALStatus status;
  INT4 xstatus;
  LALSegList seglist1;
  LALSeg seg;
  LIGOTimeGPS segstart1 = {710000000, 123456789};
  LIGOTimeGPS segend1 =   {710000234, 555555555};

  /*-- Default debug level includes info messages (4), but not
     memory checking (16), error messages (1), or warning messages (2) --*/

  /*-------------------------------------------------------------------------*/
  XLALPrintInfo("\n========== Initial setup \n");
  /*-------------------------------------------------------------------------*/

  /*-------------------------------------------------------------------------*/
  /* Initialize the segment list */
  xstatus = XLALSegListInit( &seglist1 );
  if ( xstatus ) {
    RETFAIL( "Initializing segment list", xstatus );
    XLALClearErrno();
  }

  /* Add one segment to the segment list */
  xstatus = XLALSegSet( &seg, &segstart1, &segend1, 29 );
  if ( xstatus ) {
    RETFAIL( "Creating segment from GPS times", xstatus );
    XLALClearErrno();
  } else {
    xstatus = XLALSegListAppend( &seglist1, &seg );
    if ( xstatus ) {
      RETFAIL( "Appending initial segment to segment list", xstatus );
      XLALClearErrno();
    }
  }


  /*-------------------------------------------------------------------------*/
  XLALPrintInfo("\n========== LALSegListRead tests \n");
  /*-------------------------------------------------------------------------*/

  LALSegListRead( &status, &seglist1, TEST_DATA_DIR "SegmentsInput1.data", "" );
  if ( status.statusCode ) {
    RETFAIL( "LALSegListRead with standard segment list file",
	     status.statusCode );
    REPORTSTATUS( &status );
  }

#if 0
  /*-- Dump out list of segments that was read --*/
  printf( "Segments:\n" );
  for ( segp=seglist1.segs; segp<seglist1.segs+seglist1.length; segp++ ) {
    printf( "  %10d.%09d - %10d.%09d  %d\n",
	    segp->start.gpsSeconds, segp->start.gpsNanoSeconds,
	    segp->end.gpsSeconds, segp->end.gpsNanoSeconds,
	    segp->id );
  }
#endif

  /*-------------------------------------------------------------------------*/
  XLALPrintInfo("\n========== LALSegListWrite tests \n");
  /*-------------------------------------------------------------------------*/

  LALSegListWrite( &status, &seglist1, "SegmentsOutput1.data", "adi" );
  if ( status.statusCode ) {
    RETFAIL( "LALSegListWrite with standard segment list", status.statusCode );
    REPORTSTATUS( &status );
    nfailures++;
  } else {
    XLALPrintInfo( "Wrote segment list file SegmentsOutput1.data\n" );
  }

  /*-------------------------------------------------------------------------*/
  /* Clean up leftover seg lists */
  if ( seglist1.segs ) { XLALSegListClear( &seglist1 ); }

  /*-------------------------------------------------------------------------*/
  LALCheckMemoryLeaks();

  /*-------------------------------------------------------------------------*/
  if ( nfailures == 0 ) {
    XLALPrintInfo("\n= = = = All tests passed = = = =\n\n");
    return 0;
  } else {
    XLALPrintInfo("\n= = = = %d tests FAILED = = = =\n\n", nfailures);
    return 1;
  }
}
/** \endcond */
