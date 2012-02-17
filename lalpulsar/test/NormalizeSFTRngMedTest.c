/*
*  Copyright (C) 2007 Badri Krishnan
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
 * File Name: SFTCleanTest.c
 * Authors:  Krishnan, B.
 *
 *
 * History:   Created by Krishnan August 2005
 *
 *
 *-----------------------------------------------------------------------
 */

#define LAL_USE_OLD_COMPLEX_STRUCTS
#include <lal/FrequencySeries.h>
#include <lal/NormalizeSFTRngMed.h>
#include <lal/Units.h>

#define INIT_MEM(x) memset(&(x), 0, sizeof((x)))
#define REL_ERR(x,y) ( fabs((x) - (y)) / fabs( (x) ) )

/* Default parameters. */

INT4 lalDebugLevel=3;
LALStatus empty_status;

REAL8 tol = 1e-16;

int main ( void )
{
  const char *fn = __func__;

  LALStatus status = empty_status;

  SFTtype *mySFT;
  LIGOTimeGPS epoch = { 731210229, 0 };
  REAL8 dFreq = 1.0 / 1800.0;
  REAL8 f0 = 150.0 - 2.0 * dFreq;

  /* init data array */
  COMPLEX8 vals[] = {
    { -1.249241e-21,   1.194085e-21 },
    {  2.207420e-21,   2.472366e-22 },
    {  1.497939e-21,   6.593609e-22 },
    {  3.544089e-20,  -9.365807e-21 },
    {  1.292773e-21,  -1.402466e-21 }
  };
  /* reference result for 3-bin block running-median computed in octave:
octave> sft = [ \
        -1.249241e-21 +  1.194085e-21i, \
         2.207420e-21 +  2.472366e-22i, \
         1.497939e-21 +  6.593609e-22i, \
         3.544089e-20 -  9.365807e-21i, \
         1.292773e-21 -  1.402466e-21i  \
         ];
octave> periodo = abs(sft).^2;
octave> m1 = median ( periodo(1:3) ); m2 = median ( periodo(2:4) ); m3 = median ( periodo (3:5 ) );
octave> rngmed = [ m1, m1, m2, m3, m3 ];
octave> printf ("rngmedREF3 = { %.16g, %.16g, %.16g, %.16g, %.16g };\n", rngmed );
printf ("rngmedREF3[] = { %.16g, %.16g, %.16g, %.16g, %.16g };\n", rngmed );
        rngmedREF3[] = { 2.986442063306e-42, 2.986442063306e-42, 4.933828992779561e-42, 3.638172910684999e-42, 3.638172910684999e-42 };
  */
  REAL8 rngmedREF3[] = { 2.986442063306e-42, 2.986442063306e-42, 4.933828992779561e-42, 3.638172910684999e-42, 3.638172910684999e-42 };

  UINT4 numBins = sizeof ( vals ) / sizeof(vals[0] );

  if ( (mySFT = XLALCreateSFT ( numBins )) == NULL ) {
    XLALPrintError ("%s: Failed to create test-SFT using XLALCreateSFT(), xlalErrno = %d\n", fn, xlalErrno );
    return XLAL_EFAILED;
  }
  /* init header */
  strcpy ( mySFT->name, "H1;testSFTRngmed" );
  mySFT->epoch = epoch;
  mySFT->f0 = f0;
  mySFT->deltaF = dFreq;

  /* we simply copy over these data-values into the SFT */
  UINT4 iBin;
  for ( iBin = 0; iBin < numBins; iBin ++ )
    mySFT->data->data[iBin] = vals[iBin];

  /* test running-median PSD estimation is simple blocksize cases */
  UINT4 blockSize = 3;

  /* get median->mean bias correction, needed for octave-reference results, to make
   * them comparable to the bias-corrected results from LALSFTtoRngmed()
   */
  REAL8 medianBias;
  LALRngMedBias ( &status, &medianBias, blockSize );
  if ( status.statusCode != 0 ) {
    XLALPrintError ("Something failed in LALRngMedBias(), statusCode = %d\n", status.statusCode );
    XLAL_ERROR ( XLAL_EFAILED );
  }

  /* get memory for running-median vector */
  REAL8FrequencySeries rngmed;
  INIT_MEM ( rngmed );
  if ( (rngmed.data = XLALCreateREAL8Vector ( numBins )) == NULL )
    XLAL_ERROR ( XLAL_EFAILED );
  /* compute running median */
  LALSFTtoRngmed ( &status, &rngmed, mySFT, blockSize );
  if ( status.statusCode != 0 ) {
    XLALPrintError ("Something failed in LALSFTtoRngmed(), statusCode = %d\n", status.statusCode );
    XLAL_ERROR ( XLAL_EFAILED );
  }

  BOOLEAN pass = 1;
  const CHAR *passStr;
  printf ("%4s %22s %22s %8s    <%g\n", "Bin", "rngmed(LAL)", "rngmed(Octave)", "relError", tol);
  for (iBin=0; iBin < numBins; iBin ++ )
    {
      REAL8 rngmedVAL = rngmed.data->data[iBin];
      REAL8 rngmedREF = rngmedREF3[iBin] / medianBias;	// apply median-bias correction
      REAL8 relErr = REL_ERR ( rngmedREF, rngmedVAL );
      if ( relErr > tol ) {
        pass = 0;
        passStr = "fail";
      } else {
        passStr = "OK.";
      }

      printf ("%4d %22.16g %22.16g %8.1g    %s\n", iBin, rngmedVAL, rngmedREF, relErr, passStr );

    } /* for iBin < numBins */

  /* free memory */
  XLALDestroyREAL8Vector ( rngmed.data );
  XLALDestroySFT ( mySFT );

  LALCheckMemoryLeaks();

  if ( !pass )
    {
      printf ("Test failed! Difference exceeded tolerance.\n");
      return XLAL_EFAILED;
    }
  else
    {
      printf ("Test passed.\n");
      return XLAL_SUCCESS;
    }

} /* main() */
