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

/**
 * \author Peter Shawhan
 * \file
 * \ingroup Segments_h
 * \brief Tests the segment and segment list manipulation functions.
 *
 * ### Usage ###
 *
 * \code
 * SegmentsTest [ lalDebugLevel ]
 * \endcode
 *
 * The default value of \c lalDebugLevel is 4.
 *
 * ### Description ###
 *
 * This program tests the various XLAL functions which deal with segments
 * and segment lists.
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
 */

#include <stdlib.h>
#include <stdio.h>
#include <lal/LALStdlib.h>
#include <lal/Segments.h>
#include <lal/Date.h>

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



int main( int argc, char *argv[] )
{
  INT4 nfailures = 0;
  INT4 status;
  int compval;
  CHAR testname[256];
  LALSeg seg;
  LALSeg segc1, segc2;
  LALSeg *segptr;
  LALSegList seglist1, seglist2;
  INT4 seglist1ok;
  INT4 nsegs, ilast, itime, iseg, nsearchtests;
  INT4 nbigsize = 1000000, nfullbinary = 1000000;

  /*-- Test times --*/
  LIGOTimeGPS time1  = { 794285000, 602350000 };
  LIGOTimeGPS time2  = { 794285010, 902350000 };
  LIGOTimeGPS time3  = {1794286020, 702351111 };
  LIGOTimeGPS time3p = { 394285040, 502351234 };
  LIGOTimeGPS time4a = { 794285020, 400000002 };
  LIGOTimeGPS time4b = { 794285020, 555555555 };
  LIGOTimeGPS time4c = { 794285030, 702351111 };

  LIGOTimeGPS dtime[10] = { {794000000,         0},
			    {794000000, 100000000},
			    {794000000, 220000000},
			    {794000000, 333000000},
			    {794000000, 444400000},
			    {794000000, 555550000},
			    {794000000, 666666000},
			    {794000000, 777777700},
			    {794000000, 888888880},
			    {794000000, 999999999} };

  /*-- Test segments --*/
  LALSeg seg1   = {  {794285000,602350000},  {794285010,902350000}, 1 };
  LALSeg seg2   = {  {794285020,400000002},  {794285030,702351111}, 2 };
  /* seg3p is a point in time */
  LALSeg seg3p  = {  {394285040,502351234},  {394285040,502351234}, 3 };
  /* seg4a and seg4b touch */
  LALSeg seg4a  = {  {794285050,602350000},  {794285051,300200100}, 4 };
  LALSeg seg4b  = {  {794285051,300200100},  {794285055,902350000}, 44 };
  /* seg5a and seg5b have the same start time but different end times */
  LALSeg seg5a  = { {1794285060,604350000}, {1794285061,400300200}, 5 };
  LALSeg seg5b  = { {1794285060,604350000}, {1794285061,500400300}, 55 };
  /* seg6a and seb6c overlap with seg6b, but not with each other */
  LALSeg seg6a  = {  {794285062,444444444},  {794285064,333333333}, 6 };
  LALSeg seg6b  = {  {794285063,666666666},  {794285066,555555555}, 66 };
  LALSeg seg6c  = {  {794285065,888888888},  {794285067,777777777}, 666 };
  /* segbad is an invalid segment (end time earlier than start time) */
  LALSeg segbad = {  {794285020,602350000},  {794285010,902350000}, 999 };

  /* The following segments are used to build a segment list for testing
     XLALSegListSearch.  The segments are added one at a time, and the
     search function is tested after each one is added.  As the list is
     built, it starts out disjoint, then disjoint but with touching segments,
     then not disjoint but sorted, then not sorted. */
  LALSeg sseg[] = {  {{794285000,0},  {794285001,0}, 1},
		     {{794285002,0},  {794285003,0}, 2},
		     {{794285004,0},  {794285005,0}, 3},
		     {{794285005,0},  {794285006,0}, 4},
		     {{794285007,0},  {794285009,0}, 5},
		     {{794285008,0},  {794285010,0}, 6},
		     {{794285012,0},  {794285014,0}, 7},
		     {{794284998,0},  {794284999,0}, 8} };

  /* The following array is a list of GPS times for testing XLALSegListSearch.
     The 'infirst' field of the structure indicates when the GPS should first
     be found by the search, e.g. if the 'infirst' value is equal to 4, then
     that time should be found once the segment list has been built up to the
     point of having 4 or more segments. */

  typedef struct tagTestTime {
    LIGOTimeGPS gps;
    INT4 infirst;
  } TestTime;

  TestTime lalstime[] = {  {{794284997,0}, 0},
		        {{794284998,0}, 8},
		        {{794284998,5}, 8},
		        {{794284999,5}, 0},
		        {{794285000,0}, 1},
		        {{794285000,5}, 1},
		        {{794285000,999999999}, 1},
		        {{794285001,0}, 0},
		        {{794285001,5}, 0},
		        {{794285002,0}, 2},
		        {{794285002,5}, 2},
		        {{794285003,0}, 0},
		        {{794285003,5}, 0},
		        {{794285004,0}, 3},
		        {{794285004,5}, 3},
		        {{794285005,0}, 4},
		        {{794285005,5}, 4},
		        {{794285006,0}, 0},
		        {{794285006,5}, 0},
		        {{794285007,0}, 5},
		        {{794285007,5}, 5},
		        {{794285008,0}, 5},
		        {{794285008,5}, 5},
		        {{794285009,0}, 6},
		        {{794285009,5}, 6},
		        {{794285010,5}, 0},
		        {{794285011,0}, 0},
		        {{794285012,0}, 7},
		        {{794285012,5}, 7},
		        {{794285015,0}, 0} };

  /*-- Default debug level includes info messages (4), but not
     memory checking (16), error messages (1), or warning messages (2) --*/

  /*------ Parse input line. ------*/
  if ( argc != 1 )
    {
      fprintf( stderr, "Usage: %s [ lalDebugLevel ]\n", argv[0] );
      return 0; /* so that test script won't fail */
    }

  /*-------------------------------------------------------------------------*/
  XLALPrintInfo("\n========== XLALSegSet tests \n");
  /*-------------------------------------------------------------------------*/

  /*------------------------------*/
  strcpy( testname, "XLALSegSet normal operation" );
  status = XLALSegSet( &seg, &time1, &time2, 5 );
  if ( status ) {
    RETFAIL( testname, status );
  } else if ( XLALGPSCmp(&seg.start,&time1) == 0 &&
	      XLALGPSCmp(&seg.end,&time2) == 0 &&
	      seg.id == 5 ) {
    FUNCPASS( testname );
  } else {
    FUNCFAIL( testname, "Set wrong start time, end time, and/or id" );
  }

  /*------------------------------*/
  strcpy( testname, "XLALSegSet zero-length segment" );
  status = XLALSegSet( &seg, &time3, &time3, 6 );
  if ( status ) {
    RETFAIL( testname, status );
  } else if ( XLALGPSCmp(&seg.start,&time3) == 0 &&
	      XLALGPSCmp(&seg.end,&time3) == 0 &&
	      seg.id == 6 ) {
    FUNCPASS( testname );
  } else {
    FUNCFAIL( testname, "Set wrong start time, end time, and/or id" );
  }

  /*------------------------------*/
  strcpy( testname, "XLALSegSet invalid segment" );
  status = XLALSegSet( &seg, &time2, &time1, 0 );
  if ( status == XLAL_FAILURE && xlalErrno == XLAL_EDOM ) {
    RETPASS( testname, status );
  } else {
    RETFAIL( testname, status );
  }

  /*------------------------------*/
  strcpy( testname, "XLALSegSet null segment pointer" );
  status = XLALSegSet( NULL, &time1, &time2, 0 );
  if ( status == XLAL_FAILURE && xlalErrno == XLAL_EFAULT ) {
    RETPASS( testname, status );
  } else {
    RETFAIL( testname, status );
  }

  /*------------------------------*/
  strcpy( testname, "XLALSegSet null time pointer 1" );
  status = XLALSegSet( &seg, NULL, &time3, 0 );
  if ( status == XLAL_FAILURE && xlalErrno == XLAL_EFAULT ) {
    RETPASS( testname, status );
  } else {
    RETFAIL( testname, status );
  }

  /*------------------------------*/
  strcpy( testname, "XLALSegSet null time pointer 2" );
  status = XLALSegSet( &seg, &time2, NULL, 0 );
  if ( status == XLAL_FAILURE && xlalErrno == XLAL_EFAULT ) {
    RETPASS( testname, status );
  } else {
    RETFAIL( testname, status );
  }

  /*-------------------------------------------------------------------------*/
  XLALPrintInfo("\n========== XLALSegCreate tests \n");
  /*-------------------------------------------------------------------------*/

  /*------------------------------*/
  strcpy( testname, "XLALSegCreate normal operation" );
  segptr = XLALSegCreate( &time1, &time2, 5 );
  if ( segptr == NULL ) {
    RETFAIL( testname, segptr );
  } else if ( XLALGPSCmp(&(segptr->start),&time1) == 0 &&
	      XLALGPSCmp(&(segptr->end),&time2) == 0 &&
	      segptr->id == 5 ) {
    FUNCPASS( testname );
  } else {
    FUNCFAIL( testname, "Set wrong start time, end time, and/or id" );
  }
  if (segptr) LALFree( segptr );

  /*------------------------------*/
  strcpy( testname, "XLALSegCreate zero-length segment" );
  segptr = XLALSegCreate( &time3, &time3, 6 );
  if ( segptr == NULL ) {
    RETFAIL( testname, segptr );
  } else if ( XLALGPSCmp(&(segptr->start),&time3) == 0 &&
	      XLALGPSCmp(&(segptr->end),&time3) == 0 &&
	      segptr->id == 6 ) {
    FUNCPASS( testname );
  } else {
    FUNCFAIL( testname, "Set wrong start time, end time, and/or id" );
  }
  if (segptr) LALFree( segptr );

  /*------------------------------*/
  strcpy( testname, "XLALSegCreate invalid segment" );
  segptr = XLALSegCreate( &time2, &time1, 0 );
  if ( segptr == NULL && xlalErrno == XLAL_EDOM ) {
    RETPASS( testname, segptr );
  } else {
    RETFAIL( testname, segptr );
  }
  if (segptr) LALFree( segptr );

  /*------------------------------*/
  strcpy( testname, "XLALSegCreate null time pointer 1" );
  segptr = XLALSegCreate( NULL, &time3, 0 );
  if ( segptr == NULL && xlalErrno == XLAL_EFAULT ) {
    RETPASS( testname, segptr );
  } else {
    RETFAIL( testname, segptr );
  }
  if (segptr) LALFree( segptr );

  /*------------------------------*/
  strcpy( testname, "XLALSegCreate null time pointer 2" );
  segptr = XLALSegCreate( &time2, NULL, 0 );
  if ( segptr == NULL && xlalErrno == XLAL_EFAULT ) {
    RETPASS( testname, segptr );
  } else {
    RETFAIL( testname, segptr );
  }
  if (segptr) LALFree( segptr );


  /*-------------------------------------------------------------------------*/
  XLALPrintInfo("\n========== XLALGPSInSeg tests \n");
  /*-------------------------------------------------------------------------*/

  /*------------------------------*/
  strcpy( testname, "XLALGPSInSeg normal comparison 1" );
  compval = XLALGPSInSeg( &time1, &seg2 );
  if ( xlalErrno == 0 ) {
    if ( compval < 0 ) {
      FUNCPASS( testname );
    } else {
      FUNCFAIL( testname, "Returned wrong comparison value" );
    }
  } else {
    RETFAIL( testname, compval );
  }

  /*------------------------------*/
  strcpy( testname, "XLALGPSInSeg normal comparison 3" );
  compval = XLALGPSInSeg( &time3, &seg2 );
  if ( xlalErrno == 0 ) {
    if ( compval > 0 ) {
      FUNCPASS( testname );
    } else {
      FUNCFAIL( testname, "Returned wrong comparison value" );
    }
  } else {
    RETFAIL( testname, compval );
  }

  /*------------------------------*/
  strcpy( testname, "XLALGPSInSeg normal comparison 4a" );
  compval = XLALGPSInSeg( &time4a, &seg2 );
  if ( xlalErrno == 0 ) {
    if ( compval == 0 ) {
      FUNCPASS( testname );
    } else {
      FUNCFAIL( testname, "Returned wrong comparison value" );
    }
  } else {
    RETFAIL( testname, compval );
  }

  /*------------------------------*/
  strcpy( testname, "XLALGPSInSeg normal comparison 4b" );
  compval = XLALGPSInSeg( &time4b, &seg2 );
  if ( xlalErrno == 0 ) {
    if ( compval == 0 ) {
      FUNCPASS( testname );
    } else {
      FUNCFAIL( testname, "Returned wrong comparison value" );
    }
  } else {
    RETFAIL( testname, compval );
  }

  /*------------------------------*/
  strcpy( testname, "XLALGPSInSeg normal comparison 4c" );
  compval = XLALGPSInSeg( &time4c, &seg2 );
  if ( xlalErrno == 0 ) {
    if ( compval > 0 ) {
      FUNCPASS( testname );
    } else {
      FUNCFAIL( testname, "Returned wrong comparison value" );
    }
  } else {
    RETFAIL( testname, compval );
  }

  /*------------------------------*/
  strcpy( testname, "XLALGPSInSeg normal comparison 3p" );
  compval = XLALGPSInSeg( &time3p, &seg3p );
  if ( xlalErrno == 0 ) {
    if ( compval == 0 ) {
      FUNCPASS( testname );
    } else {
      FUNCFAIL( testname, "Returned wrong comparison value" );
    }
  } else {
    RETFAIL( testname, compval );
  }

  /*------------------------------*/
  strcpy( testname, "XLALGPSInSeg null GPS pointer" );
  compval = XLALGPSInSeg( NULL, &seg2 );
  if ( xlalErrno == XLAL_EFAULT ) {
    RETPASS( testname, compval );
  } else {
    RETFAIL( testname, compval );
  }

  /*------------------------------*/
  strcpy( testname, "XLALGPSInSeg null seg pointer" );
  compval = XLALGPSInSeg( &time1, NULL );
  if ( xlalErrno == XLAL_EFAULT ) {
    RETPASS( testname, compval );
  } else {
    RETFAIL( testname, compval );
  }


  /*-------------------------------------------------------------------------*/
  XLALPrintInfo("\n========== XLALSegCmp tests \n");
  /*-------------------------------------------------------------------------*/

  /*------------------------------*/
  strcpy( testname, "XLALSegCmp normal comparison (1,2)" );
  compval = XLALSegCmp( &seg1, &seg2 );
  if ( xlalErrno == 0 ) {
    if ( compval < 0 ) {
      FUNCPASS( testname );
    } else {
      FUNCFAIL( testname, "Returned wrong comparison value" );
    }
  } else {
    RETFAIL( testname, compval );
  }

  /*------------------------------*/
  strcpy( testname, "XLALSegCmp normal comparison (2,1)" );
  compval = XLALSegCmp( &seg2, &seg1 );
  if ( xlalErrno == 0 ) {
    if ( compval > 0 ) {
      FUNCPASS( testname );
    } else {
      FUNCFAIL( testname, "Returned wrong comparison value" );
    }
  } else {
    RETFAIL( testname, compval );
  }

  /*------------------------------*/
  strcpy( testname, "XLALSegCmp normal comparison (2,2)" );
  compval = XLALSegCmp( &seg2, &seg2 );
  if ( xlalErrno == 0 ) {
    if ( compval == 0 ) {
      FUNCPASS( testname );
    } else {
      FUNCFAIL( testname, "Returned wrong comparison value" );
    }
  } else {
    RETFAIL( testname, compval );
  }

  /*------------------------------*/
  strcpy( testname, "XLALSegCmp normal comparison (1,3p)" );
  compval = XLALSegCmp( &seg1, &seg3p );
  if ( xlalErrno == 0 ) {
    if ( compval > 0 ) {
      FUNCPASS( testname );
    } else {
      FUNCFAIL( testname, "Returned wrong comparison value" );
    }
  } else {
    RETFAIL( testname, compval );
  }

  /*------------------------------*/
  strcpy( testname, "XLALSegCmp normal comparison (3p,3p)" );
  compval = XLALSegCmp( &seg3p, &seg3p );
  if ( xlalErrno == 0 ) {
    if ( compval == 0 ) {
      FUNCPASS( testname );
    } else {
      FUNCFAIL( testname, "Returned wrong comparison value" );
    }
  } else {
    RETFAIL( testname, compval );
  }

  /*------------------------------*/
  strcpy( testname, "XLALSegCmp normal comparison (5a,5b)" );
  compval = XLALSegCmp( &seg5a, &seg5b );
  if ( xlalErrno == 0 ) {
    if ( compval < 0 ) {
      FUNCPASS( testname );
    } else {
      FUNCFAIL( testname, "Returned wrong comparison value" );
    }
  } else {
    RETFAIL( testname, compval );
  }

  /*------------------------------*/
  strcpy( testname, "XLALSegCmp normal comparison (5b,5a)" );
  compval = XLALSegCmp( &seg5b, &seg5a );
  if ( xlalErrno == 0 ) {
    if ( compval > 0 ) {
      FUNCPASS( testname );
    } else {
      FUNCFAIL( testname, "Returned wrong comparison value" );
    }
  } else {
    RETFAIL( testname, compval );
  }

  /*------------------------------*/
  strcpy( testname, "XLALSegCmp null seg pointer 1" );
  compval = XLALSegCmp( NULL, &seg5a );
  if ( xlalErrno == XLAL_EFAULT ) {
    RETPASS( testname, compval );
  } else {
    RETFAIL( testname, compval );
  }

  /*------------------------------*/
  strcpy( testname, "XLALSegCmp null seg pointer 2" );
  compval = XLALSegCmp( &seg5b, NULL );
  if ( xlalErrno == XLAL_EFAULT ) {
    RETPASS( testname, compval );
  } else {
    RETFAIL( testname, compval );
  }


  /*-------------------------------------------------------------------------*/
  XLALPrintInfo("\n========== XLALSegListInit tests \n");
  /*-------------------------------------------------------------------------*/

  /*------------------------------*/
  strcpy( testname, "XLALSegListInit normal operation" );
  seglist1.initMagic = -123;
  status = XLALSegListInit( &seglist1 );
  if ( status ) {
    RETFAIL( testname, status );
  } else if ( seglist1.initMagic == SEGMENTSH_INITMAGICVAL ) {
    FUNCPASS( testname );
  } else {
    FUNCFAIL( testname, "Failed to set initMagic correctly" );
  }

  /*------------------------------*/
  strcpy( testname, "XLALSegListInit null seglist pointer" );
  status = XLALSegListInit( NULL );
  if ( status == XLAL_FAILURE && xlalErrno == XLAL_EFAULT ) {
    RETPASS( testname, status );
  } else {
    RETFAIL( testname, status );
  }


  /*-------------------------------------------------------------------------*/
  XLALPrintInfo("\n========== XLALSegListClear tests \n");
  /*-------------------------------------------------------------------------*/

  /*------------------------------*/
  strcpy( testname, "XLALSegListClear normal operation (empty)" );
  seglist1.sorted = 0;
  status = XLALSegListClear( &seglist1 );
  if ( status ) {
    RETFAIL( testname, status );
  } else if ( seglist1.sorted == 1 ) {
    FUNCPASS( testname );
  } else {
    FUNCFAIL( testname, "Failed to reset sorted correctly" );
  }

  /*------------------------------*/
  strcpy( testname, "XLALSegListClear null seglist pointer" );
  status = XLALSegListClear( NULL );
  if ( status == XLAL_FAILURE && xlalErrno == XLAL_EFAULT ) {
    RETPASS( testname, status );
  } else {
    RETFAIL( testname, status );
  }

  /*------------------------------*/
  strcpy( testname, "XLALSegListClear uninitialized seglist" );
  seglist1.initMagic = -444;
  status = XLALSegListClear( &seglist1 );
  if ( status == XLAL_FAILURE && xlalErrno == XLAL_EINVAL ) {
    RETPASS( testname, status );
  } else {
    RETFAIL( testname, status );
  }

  /*------------------------------*/
  strcpy( testname, "XLALSegListClear normal operation (filled) - init" );
  status = XLALSegListInit( &seglist1 );
  if ( status ) {
    RETFAIL( testname, status );
  } else if ( seglist1.segs != NULL ) {
    FUNCFAIL( testname, "Failed to init segment list" );
  } else {

    strcpy( testname, "XLALSegListClear normal operation (filled) - fill" );
    status = XLALSegListAppend( &seglist1, &seg2 );
    if ( status ) {
      RETFAIL( testname, status );
    } else if ( seglist1.segs == NULL ) {
      FUNCFAIL( testname, "Failed to store segment 2" );
    } else {

      status = XLALSegListAppend( &seglist1, &seg1 );
      if ( status ) {
	RETFAIL( testname, status );
      } else if ( seglist1.sorted == 1 ) {
	FUNCFAIL( testname, "Failed to clear sorted flag" );
      } else {

	/*-- Finally, we're ready to test clearing --*/
	strcpy( testname, "XLALSegListClear normal operation (filled)" );
	status = XLALSegListClear( &seglist1 );
	if ( status ) {
	  RETFAIL( testname, status );
	} else if ( seglist1.sorted == 1 ) {
	  FUNCPASS( testname );
	} else {
	  FUNCFAIL( testname, "Failed to reset sorted correctly" );
	}

      }
    }
  }


  /*-------------------------------------------------------------------------*/
  XLALPrintInfo("\n========== XLALSegListAppend tests \n");
  /*-------------------------------------------------------------------------*/

  /*------------------------------*/
  strcpy( testname, "XLALSegListAppend null seglist pointer" );
  status = XLALSegListAppend( NULL, &seg1 );
  if ( status == XLAL_FAILURE && xlalErrno == XLAL_EFAULT ) {
    RETPASS( testname, status );
  } else {
    RETFAIL( testname, status );
  }

  /*------------------------------*/
  strcpy( testname, "XLALSegListAppend null seg pointer" );
  status = XLALSegListAppend( &seglist1, NULL );
  if ( status == XLAL_FAILURE && xlalErrno == XLAL_EFAULT ) {
    RETPASS( testname, status );
  } else {
    RETFAIL( testname, status );
  }

  /*------------------------------*/
  strcpy( testname, "XLALSegListAppend uninitialized seglist" );
  seglist1.initMagic = -444;
  status = XLALSegListAppend( &seglist1, &seg1 );
  if ( status == XLAL_FAILURE && xlalErrno == XLAL_EINVAL ) {
    RETPASS( testname, status );
  } else {
    RETFAIL( testname, status );
  }

  /*------------------------------*/
  strcpy( testname, "XLALSegListAppend invalid segment - init" );
  status = XLALSegListInit( &seglist1 );
  if ( status ) {
    RETFAIL( testname, status );
  } else if ( seglist1.segs != NULL ) {
    FUNCFAIL( testname, "Failed to init segment list" );
  } else {

    strcpy( testname, "XLALSegListAppend invalid segment" );
    status = XLALSegListAppend( &seglist1, &segbad );
    if ( status == XLAL_FAILURE && xlalErrno == XLAL_EDOM ) {
      RETPASS( testname, status );
    } else {
      RETFAIL( testname, status );
    }

    /* Clean up */
    XLALSegListClear( &seglist1 );

  }

  /*------------------------------*/
  strcpy( testname, "XLALSegListAppend decimal place record - init" );
  status = XLALSegListInit( &seglist1 );
  if ( status ) {
    RETFAIL( testname, status );
  } else if ( seglist1.segs != NULL ) {
    FUNCFAIL( testname, "Failed to init segment list" );
  } else if ( seglist1.dplaces != 0 ) {
    FUNCFAIL( testname, "Failed to init dplaces for segment list" );
  } else {
    /* Loop over segments */
    strcpy( testname, "XLALSegListAppend decimal place record" );
    for ( iseg=0; iseg<=9; iseg++ ) {
      XLALSegSet( &seg, dtime+0, dtime+iseg, 0 );
      XLALSegListAppend( &seglist1, &seg );
      if ( seglist1.dplaces != 3*(((UINT4)iseg+2)/3) ) {
	FUNCFAIL( testname, "Wrong dplaces value after appending segment" );
      }
    }
    for ( iseg=9; iseg>=0; iseg-- ) {
      XLALSegSet( &seg, dtime+0, dtime+iseg, 0 );
      XLALSegListAppend( &seglist1, &seg );
      if ( seglist1.dplaces != 9 ) {
	FUNCFAIL( testname, "Wrong dplaces value after re-appending segment" );
      }
    }

    /* Clean up */
    XLALSegListClear( &seglist1 );

  }

  /*------------------------------*/
  strcpy( testname, "XLALSegListAppend normal operation (empty) - init" );
  status = XLALSegListInit( &seglist1 );
  if ( status ) {
    RETFAIL( testname, status );
    seglist1ok = 0;
  } else if ( seglist1.segs != NULL ) {
    FUNCFAIL( testname, "Failed to init segment list" );
    seglist1ok = 0;
  } else {
    seglist1ok = 1;
  }

  strcpy( testname, "XLALSegListAppend first segment" );
  if ( seglist1ok ) {
    status = XLALSegListAppend( &seglist1, &seg1 );
    if ( status ) {
      RETFAIL( testname, status );
      seglist1ok = 0;
    } else if ( seglist1.segs && seglist1.arraySize == SEGMENTSH_ALLOCBLOCK &&
		seglist1.length == 1 ) {
      FUNCPASS( testname );
    } else {
      FUNCFAIL( testname, "Failed to store first segment" );
      seglist1ok = 0;
    }
  } else {
    XLALPrintInfo( "Skip test: %s\n", testname );
  }

  /*-- Continue with a second segment --*/
  strcpy( testname, "XLALSegListAppend sorted, disjoint segment" );
  if ( seglist1ok ) {
    status = XLALSegListAppend( &seglist1, &seg6a );
    if ( status ) {
      RETFAIL( testname, status );
      seglist1ok = 0;
    } else if ( seglist1.length == 2 && seglist1.sorted &&
		seglist1.disjoint ) {
      FUNCPASS( testname );
    } else {
      FUNCFAIL( testname, "Failed to store second segment" );
      seglist1ok = 0;
    }
  } else {
    XLALPrintInfo( "Skip test: %s\n", testname );
  }


  /*-- Continue with a third segment which is not disjoint --*/
  strcpy( testname, "XLALSegListAppend non-disjoint segment" );
  if ( seglist1ok ) {
    status = XLALSegListAppend( &seglist1, &seg6b );
    if ( status ) {
      RETFAIL( testname, status );
      seglist1ok = 0;
    } else if ( seglist1.length == 3 && seglist1.sorted &&
		! seglist1.disjoint ) {
      FUNCPASS( testname );
    } else {
      FUNCFAIL( testname, "Failed to store second segment" );
      seglist1ok = 0;
    }
  } else {
    XLALPrintInfo( "Skip test: %s\n", testname );
  }

  /*-- Continue with a fourth segment which is not sorted --*/
  strcpy( testname, "XLALSegListAppend non-sorted segment" );
  if ( seglist1ok ) {
    status = XLALSegListAppend( &seglist1, &seg4a );
    if ( status ) {
      RETFAIL( testname, status );
      seglist1ok = 0;
    } else if ( seglist1.length == 4 && ! seglist1.sorted &&
		! seglist1.disjoint ) {
      FUNCPASS( testname );
    } else {
      FUNCFAIL( testname, "Failed to store second segment" );
      seglist1ok = 0;
    }
  } else {
    XLALPrintInfo( "Skip test: %s\n", testname );
  }


  /*-------------------------------------------------------------------------*/
  XLALPrintInfo("\n========== XLALSegListSort tests \n");
  /*-------------------------------------------------------------------------*/

  /*------------------------------*/
  strcpy( testname, "XLALSegListSort null seglist pointer" );
  status = XLALSegListSort( NULL );
  if ( status == XLAL_FAILURE && xlalErrno == XLAL_EFAULT ) {
    RETPASS( testname, status );
  } else {
    RETFAIL( testname, status );
  }

  /*------------------------------*/
  strcpy( testname, "XLALSegListSort uninitialized seglist" );
  seglist2.initMagic = -444;
  status = XLALSegListSort( &seglist2 );
  if ( status == XLAL_FAILURE && xlalErrno == XLAL_EINVAL ) {
    RETPASS( testname, status );
  } else {
    RETFAIL( testname, status );
  }

  /*------------------------------*/
  strcpy( testname, "XLALSegListSort normal operation" );
  /*-- Sort the existing segment list --*/
  if ( seglist1ok ) {
    status = XLALSegListSort( &seglist1 );
    if ( status ) {
      RETFAIL( testname, status );
      seglist1ok = 0;
    } else if ( seglist1.length == 4 &&
		XLALSegCmp( &(seglist1.segs[0]), &seg1 ) == 0 &&
		XLALSegCmp( &(seglist1.segs[1]), &seg4a ) == 0 &&
		XLALSegCmp( &(seglist1.segs[2]), &seg6a ) == 0 &&
		XLALSegCmp( &(seglist1.segs[3]), &seg6b ) == 0 &&
		seglist1.sorted && ! seglist1.disjoint ) {
      FUNCPASS( testname );
    } else {
      FUNCFAIL( testname, "Did not sort correctly" );
      seglist1ok = 0;
    }
  } else {
    XLALPrintInfo( "Skip test: %s\n", testname );
  }


  /*-------------------------------------------------------------------------*/
  XLALPrintInfo("\n========== XLALSegListCoalesce tests \n");
  /*-------------------------------------------------------------------------*/

  /*------------------------------*/
  strcpy( testname, "XLALSegListCoalesce null seglist pointer" );
  status = XLALSegListCoalesce( NULL );
  if ( status == XLAL_FAILURE && xlalErrno == XLAL_EFAULT ) {
    RETPASS( testname, status );
  } else {
    RETFAIL( testname, status );
  }

  /*------------------------------*/
  strcpy( testname, "XLALSegListCoalesce uninitialized seglist" );
  seglist2.initMagic = -444;
  status = XLALSegListCoalesce( &seglist2 );
  if ( status == XLAL_FAILURE && xlalErrno == XLAL_EINVAL ) {
    RETPASS( testname, status );
  } else {
    RETFAIL( testname, status );
  }

  /*------------------------------*/
  strcpy( testname, "XLALSegListCoalesce normal operation 1" );
  /*-- Coalesce the existing segment list --*/
  if ( seglist1ok ) {
    status = XLALSegListCoalesce( &seglist1 );
    /*-- Construct segment for comparison purposes --*/
    XLALSegSet( &segc1, &(seg6a.start), &(seg6b.end), 8 );
    if ( status ) {
      RETFAIL( testname, status );
      seglist1ok = 0;
    } else if ( seglist1.length == 3 &&
		XLALSegCmp( &(seglist1.segs[0]), &seg1 ) == 0 &&
		XLALSegCmp( &(seglist1.segs[1]), &seg4a ) == 0 &&
		XLALSegCmp( &(seglist1.segs[2]), &segc1 ) == 0 &&
		seglist1.sorted && seglist1.disjoint ) {
      FUNCPASS( testname );
    } else {
      FUNCFAIL( testname, "Did not coalesce correctly" );
      seglist1ok = 0;
    }
  } else {
    XLALPrintInfo( "Skip test: %s\n", testname );
  }


  /*------------------------------*/
  /*-- Add some segments to the existing segment list --*/
  strcpy( testname, "XLALSegListCoalesce normal operation 2 - add 6c" );
  if ( seglist1ok ) {
    status = XLALSegListAppend( &seglist1, &seg6c );
    if ( status ) {
      RETFAIL( testname, status );
      seglist1ok = 0;
    } else if ( seglist1.length == 4 && seglist1.sorted &&
		! seglist1.disjoint ) {
      FUNCPASS( testname );
    } else {
      FUNCFAIL( testname, "Failed to store second segment" );
      seglist1ok = 0;
    }
  } else {
    XLALPrintInfo( "Skip test: %s\n", testname );
  }

  strcpy( testname, "XLALSegListCoalesce normal operation 2 - add 4b" );
  if ( seglist1ok ) {
    status = XLALSegListAppend( &seglist1, &seg4b );
    if ( status ) {
      RETFAIL( testname, status );
      seglist1ok = 0;
    } else if ( seglist1.length == 5 && ! seglist1.sorted &&
		! seglist1.disjoint ) {
      FUNCPASS( testname );
    } else {
      FUNCFAIL( testname, "Failed to store second segment" );
      seglist1ok = 0;
    }
  } else {
    XLALPrintInfo( "Skip test: %s\n", testname );
  }

  strcpy( testname, "XLALSegListCoalesce normal operation 2 - add 2" );
  if ( seglist1ok ) {
    status = XLALSegListAppend( &seglist1, &seg2 );
    if ( status ) {
      RETFAIL( testname, status );
      seglist1ok = 0;
    } else if ( seglist1.length == 6 && ! seglist1.sorted &&
		! seglist1.disjoint ) {
      FUNCPASS( testname );
    } else {
      FUNCFAIL( testname, "Failed to store second segment" );
      seglist1ok = 0;
    }
  } else {
    XLALPrintInfo( "Skip test: %s\n", testname );
  }

  strcpy( testname, "XLALSegListCoalesce normal operation 2 - add 2 again" );
  if ( seglist1ok ) {
    status = XLALSegListAppend( &seglist1, &seg2 );
    if ( status ) {
      RETFAIL( testname, status );
      seglist1ok = 0;
    } else if ( seglist1.length == 7 && ! seglist1.sorted &&
		! seglist1.disjoint ) {
      FUNCPASS( testname );
    } else {
      FUNCFAIL( testname, "Failed to store second segment" );
      seglist1ok = 0;
    }
  } else {
    XLALPrintInfo( "Skip test: %s\n", testname );
  }

  strcpy( testname, "XLALSegListCoalesce normal operation 2 - coalesce" );
  if ( seglist1ok ) {
    status = XLALSegListCoalesce( &seglist1 );
    /*-- Construct segments for comparison purposes --*/
    XLALSegSet( &segc1, &(seg4a.start), &(seg4b.end), 8 );
    XLALSegSet( &segc2, &(seg6a.start), &(seg6c.end), 8 );
    if ( status ) {
      RETFAIL( testname, status );
      seglist1ok = 0;
    } else if ( seglist1.length == 4 &&
		XLALSegCmp( &(seglist1.segs[0]), &seg1 ) == 0 &&
		XLALSegCmp( &(seglist1.segs[1]), &seg2 ) == 0 &&
		XLALSegCmp( &(seglist1.segs[2]), &segc1 ) == 0 &&
		XLALSegCmp( &(seglist1.segs[3]), &segc2 ) == 0 &&
		seglist1.sorted && seglist1.disjoint ) {
      FUNCPASS( testname );
    } else {
      FUNCFAIL( testname, "Did not coalesce correctly" );
      seglist1ok = 0;
    }
  } else {
    XLALPrintInfo( "Skip test: %s\n", testname );
  }


  /*-------------------------------------------------------------------------*/
  XLALPrintInfo("\n========== XLALSegListSearch tests \n");
  /*-------------------------------------------------------------------------*/

  /*------------------------------*/
  strcpy( testname, "XLALSegListSearch null seglist pointer" );
  segptr = XLALSegListSearch( NULL, &time1 );
  if ( segptr == NULL && xlalErrno == XLAL_EFAULT ) {
    RETPASS( testname, segptr );
  } else {
    RETFAIL( testname, segptr );
  }

  /*------------------------------*/
  strcpy( testname, "XLALSegListSearch null GPS pointer" );
  segptr = XLALSegListSearch( &seglist1, NULL );
  if ( segptr == NULL && xlalErrno == XLAL_EFAULT ) {
    RETPASS( testname, segptr );
  } else {
    RETFAIL( testname, segptr );
  }

  /*------------------------------*/
  strcpy( testname, "XLALSegListSearch uninitialized seglist" );
  seglist2.initMagic = -444;
  segptr = XLALSegListSearch( &seglist2, &time1 );
  if ( segptr == NULL && xlalErrno == XLAL_EINVAL ) {
    RETPASS( testname, segptr );
  } else {
    RETFAIL( testname, segptr );
  }

  /*------------------------------*/
  strcpy( testname, "XLALSegListSearch empty seglist - init" );
  status = XLALSegListInit( &seglist2 );
  if ( status ) {
    RETFAIL( testname, status );
  } else if ( seglist2.initMagic == SEGMENTSH_INITMAGICVAL ) {
    FUNCPASS( testname );
  } else {
    FUNCFAIL( testname, "Failed to set initMagic correctly" );
  }

  strcpy( testname, "XLALSegListSearch empty seglist - search" );
  segptr = XLALSegListSearch( &seglist2, &time1 );
  if ( segptr == NULL && xlalErrno == 0 ) {
    RETPASS( testname, segptr );
  } else {
    RETFAIL( testname, segptr );
  }

  nsearchtests = 0;

  /*------------------------------*/
  /* Loop over segment list length */

  for ( nsegs=1; nsegs<=(INT4)sizeof(sseg)/(INT4)sizeof(LALSeg); nsegs++ ) {

    /* Append a segment */
    sprintf( testname,
	     "XLALSegListSearch - appending to make list of %d segments",
	     nsegs );
    status = XLALSegListAppend( &seglist2, sseg+nsegs-1 );
    if ( status ) {
      RETFAIL( testname, status );
      break;
    }

    /* Loop over test times */
    for ( itime=0; itime < (INT4)sizeof(lalstime)/(INT4)sizeof(TestTime); itime++ ) {

      /* Loop over lastFound states */
      for ( ilast = -1; ilast < nsegs; ilast++ ) {
	if ( ilast == -1 ) {
	  seglist2.lastFound = NULL;
	} else {
	  seglist2.lastFound = &(seglist2.segs[ilast]);
	}

	nsearchtests++;

	sprintf( testname,
		 "XLALSegListSearch - nsegs=%d, itime=%d, ilast=%d",
		 nsegs, itime, ilast );
	/* Do the search */
	segptr = XLALSegListSearch( &seglist2, &(lalstime[itime].gps) );

	/* Check whether the function had an error */
	if ( xlalErrno != 0 ) {
	  RETFAIL( testname, segptr );
	  continue;
	}

	/* Check that the result is correct */
	if ( lalstime[itime].infirst > 0 && nsegs >= lalstime[itime].infirst ) {
	  /* The search should have found a match */
	  if ( ! segptr ) {
	    FUNCFAIL( testname, "Failed to find a match" );
	    continue;
	  }
	  /* The matched segment should contain the time */
	  if ( XLALGPSInSeg(&(lalstime[itime].gps),segptr) != 0 ) {
	    FUNCFAIL( testname, "Segment returned does not contain the time" );
	    continue;
	  }
	} else {
	  /* The search should NOT have found a match */
	  if ( segptr ) {
	    FUNCFAIL( testname, "Found inappropriate match" );
	    continue;
	  }
	}

      } /* End loop over lastFound states */

    } /* End loop over test times */

  } /* End loop over number of segments in list */

  XLALPrintInfo( "Tested %d cases for XLALSegListSearch for correctness\n",
		 nsearchtests );

  /* Construct a big segment list */
  XLALSegListClear( &seglist2 );
  XLALPrintInfo( "\nConstructing a segment list with %d segments...\n",
		 nbigsize );
  for ( iseg=0; iseg<nbigsize; iseg++ ) {
    sprintf( testname, "Appending segment %d to big list", iseg );
    time1.gpsSeconds = 800000000 + iseg;
    time1.gpsNanoSeconds = 0;
    time2.gpsSeconds = 800000000 + iseg;
    time2.gpsNanoSeconds = 500000000;
    status = XLALSegSet( &seg, &time1, &time2, 0 );
    if ( status != 0 ) {
      RETFAIL( testname, status );
      break;
    }
    status = XLALSegListAppend( &seglist2, &seg );
    if ( status != 0 ) {
      RETFAIL( testname, status );
      break;
    }
  }

  /* Do a bunch of binary searches in the big segment list */
  XLALPrintInfo( "\nDoing %d full binary searches in the segment list...\n",
		 nfullbinary );
  time1.gpsSeconds = 800000000 + nbigsize/3 ;
  time1.gpsNanoSeconds = 250000000;
  for ( itime=0; itime<nfullbinary; itime++ ) {
    /* Clear the lastFound pointer, to force a full binary search */
    seglist2.lastFound = NULL;
    segptr = XLALSegListSearch( &seglist2, &time1 );
  }


  /*-------------------------------------------------------------------------*/
  /* Clean up leftover seg lists */
  if ( seglist1.segs ) { XLALSegListClear( &seglist1 ); }
  if ( seglist2.segs ) { XLALSegListClear( &seglist2 ); }

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
