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

#include <lal/FrequencySeries.h>
#include <lal/NormalizeSFTRngMed.h>
#include <lal/Units.h>

/*********************************************************************/
/* Macros for printing errors & testing subroutines (from Creighton) */
/*********************************************************************/

#define ERROR( code, msg, statement )                                \
do {                                                                 \
  if ( lalDebugLevel & LALERROR )                                    \
    XLALPrintError( "Error[0] %d: program %s, file %s, line %d, %s\n" \
                   "        %s %s\n", (code), *argv, __FILE__,       \
              __LINE__, "$Id$", statement ? statement :  \
                   "", (msg) );                                      \
} while (0)

#define INFO( statement )                                            \
do {                                                                 \
  if ( lalDebugLevel & LALINFO )                                     \
    XLALPrintError( "Info[0]: program %s, file %s, line %d, %s\n"     \
                   "        %s\n", *argv, __FILE__, __LINE__,        \
              "$Id$", (statement) );                     \
} while (0)

#define SUB( func, statusptr )                                       \
do {                                                                 \
  if ( (func), (statusptr)->statusCode ) {                           \
    ERROR( NORMALIZESFTRNGMEDC_ESUB, NORMALIZESFTRNGMEDC_MSGESUB,      \
           "Function call \"" #func "\" failed:" );                  \
    return NORMALIZESFTRNGMEDC_ESUB;                                  \
  }                                                                  \
} while (0)
/******************************************************************/

/* Error codes and messages */
#define NORMALIZESFTRNGMEDC_ENORM 0
#define NORMALIZESFTRNGMEDC_ESUB  1
#define NORMALIZESFTRNGMEDC_EARG  2
#define NORMALIZESFTRNGMEDC_EBAD  3
#define NORMALIZESFTRNGMEDC_EFILE 4

#define NORMALIZESFTRNGMEDC_MSGENORM "Normal exit"
#define NORMALIZESFTRNGMEDC_MSGESUB  "Subroutine failed"
#define NORMALIZESFTRNGMEDC_MSGEARG  "Error parsing arguments"
#define NORMALIZESFTRNGMEDC_MSGEBAD  "Bad argument values"
#define NORMALIZESFTRNGMEDC_MSGEFILE "Could not create output file"

/* Default parameters. */

INT4 lalDebugLevel=3;
LALStatus empty_status;

#define NUM_BINS  5

int main ( void )
{
  const char *fn = __func__;

  LALStatus status = empty_status;

  SFTtype *mySFT;
  LIGOTimeGPS epoch = { 731210229, 0 };
  REAL8 dFreq = 1.0 / 1800.0;
  REAL8 f0 = 150.0 - 2.0 * dFreq;

  if ( (mySFT = XLALCreateSFT ( NUM_BINS )) == NULL ) {
    XLALPrintError ("%s: Failed to create test-SFT using XLALCreateSFT(), xlalErrno = %d\n", fn, xlalErrno );
    return NORMALIZESFTRNGMEDC_ESUB;
  }
  /* init header */
  strcpy ( mySFT->name, "H1;testSFTRngmed" );
  mySFT->epoch = epoch;
  mySFT->f0 = f0;
  mySFT->deltaF = dFreq;

  /* init data array */
  COMPLEX8 vals[NUM_BINS] = {
    { -1.249241e-21,   1.194085e-21 },
    {  2.207420e-21,   2.472366e-22 },
    {  1.497939e-21,   6.593609e-22 },
    {  3.544089e-20,  -9.365807e-21 },
    {  1.292773e-21,  -1.402466e-21 }
  };
  /* we simply copy over these data-values into the SFT */
  UINT4 iBin;
  for ( iBin = 0; iBin < NUM_BINS; iBin ++ )
    mySFT->data->data[iBin] = vals[iBin];

  /* test running-median PSD estimation is simple blocksize cases */
  UINT4 blockSize = 3;

  /* get median->mean bias correction */
  REAL8 medianBias;
  SUB ( LALRngMedBias( &status, &medianBias, blockSize ), &status);






  /* free memory */
  XLALDestroySFT ( mySFT );

  LALCheckMemoryLeaks();

  return NORMALIZESFTRNGMEDC_ENORM;

} /* main() */
