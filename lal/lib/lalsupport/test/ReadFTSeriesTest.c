/*
*  Copyright (C) 2007 Cristina Valeria Torres, Jolien Creighton, John Whelan
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
 * \file
 * \ingroup ReadFTSeries_h
 * \author Torres, C. W.
 *
 * \brief Tests the routines in \ref ReadTimeSeries_c and \ref ReadFrequencySeries_c.
 *
 * ### Usage ###
 *
 * \code
 * ReadFTSeriesTest
 * \endcode
 *
 * ### Description ###
 *
 * For each of the real and complex datatypes (single and double
 * precision), this program fills a Time- and FrequencySeries, prints
 * them to disk with the appropriate
 * \c LALPrint<datatype><tt>FrequencySeries()</tt>
 * and \c LALPrint<datatype><tt>TimeSeries()</tt>
 * routines, then reads them back in with the appropriate
 * \c LALRead<datatype><tt>FrequencySeries()</tt>
 * and \c LALRead<datatype><tt>TimeSeries()</tt>
 * routines and checks to make sure the resulting series agree, printing the results to standard error.
 */

/**\name Error Codes */ /*@{*/
#define READFTSERIESTESTC_ENOM 0        /**< Nominal exit */
#define READFTSERIESTESTC_ECHK 1        /**< Error checking failed to catch bad data */
#define READFTSERIESTESTC_EFUN 2        /**< Subroutine returned error for valid data */
#define READFTSERIESTESTC_EFLS 3        /**< Subroutine returned unexpected results */
/*@}*/

/** \cond DONT_DOXYGEN */
#define READFTSERIESTESTC_MSGENOM "Nominal exit"
#define READFTSERIESTESTC_MSGECHK "Error checking failed to catch bad data"
#define READFTSERIESTESTC_MSGEFUN "Subroutine returned error for valid data"
#define READFTSERIESTESTC_MSGEFLS "Subroutine returned unexpected results"


#include <lal/Units.h>
#include <lal/ReadFTSeries.h>
#include <lal/PrintFTSeries.h>
#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#define  READFTSERIESTEST_TOL 1e6
#define  READFTSERIESTEST_LEN 20



int main( void )
{
  static LALStatus status;

  REAL4Sequence       *sSequenceIn;
  REAL4Sequence       *sSequenceOut;
  REAL8Sequence       *dSequenceIn;
  REAL8Sequence       *dSequenceOut;

  COMPLEX8Sequence    *cSequenceIn;
  COMPLEX8Sequence    *cSequenceOut;
  COMPLEX16Sequence   *zSequenceIn;
  COMPLEX16Sequence   *zSequenceOut;


  REAL4FrequencySeries      sFrequencySeries;
  REAL8FrequencySeries      dFrequencySeries;
  COMPLEX8FrequencySeries   cFrequencySeries;
  COMPLEX16FrequencySeries  zFrequencySeries;

  REAL4FrequencySeries      sFrequencySeries2;
  REAL8FrequencySeries      dFrequencySeries2;
  COMPLEX8FrequencySeries   cFrequencySeries2;
  COMPLEX16FrequencySeries  zFrequencySeries2;

  REAL4TimeSeries           sTimeSeries;
  REAL8TimeSeries           dTimeSeries;
  COMPLEX8TimeSeries        cTimeSeries;
  COMPLEX16TimeSeries       zTimeSeries;

  REAL4TimeSeries           sTimeSeries2;
  REAL8TimeSeries           dTimeSeries2;
  COMPLEX8TimeSeries        cTimeSeries2;
  COMPLEX16TimeSeries       zTimeSeries2;

  REAL4            *sData;
  REAL8            *dData;
  COMPLEX8         *cData;
  COMPLEX16        *zData;

  /* Boolean Vars */
  BOOLEAN     unitComp;

  /* variables for units */
  RAT4        raise;
  LALUnitPair unitPair;
  LALUnitPair inputUnits;
  LALUnit     strainToMinus2;
  LALUnit     adcToMinus2;
  LALUnit     adcStrain;


  /* This routine should generate a file with data */
  /* to be read by ReadFTSeries.c*/
  LIGOTimeGPS  t;
  UINT4         i;


  /* Data Test Variable */
  UINT4   j;


  fprintf(stderr,"Testing value of LALUnitTextSize ... ");
  if ( (int)LALSupportUnitTextSize != (int)LALUnitTextSize )
  {
    fprintf(stderr,"UnitTextSize mismatch: [ReadFTSeriesTest:%s]\n",
	    READFTSERIESTESTC_MSGEFLS);
    return READFTSERIESTESTC_EFLS;
  }
  fprintf(stderr,"PASS\n");

  t.gpsSeconds = 0;
  t.gpsNanoSeconds = 0;

  fprintf(stderr,"Testing Print/Read COMPLEX8FrequencySeries ... ");

  cSequenceIn = NULL;
  LALCCreateVector(&status, &cSequenceIn, READFTSERIESTEST_LEN);
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }
  for ( i=1, cData=cSequenceIn->data; i<=READFTSERIESTEST_LEN ; i++, cData++ )
  {
    *(cData) = crectf( i, i+1 );
    }

  strncpy(cFrequencySeries.name,"Complex frequency series",LALNameLength);
  cFrequencySeries.sampleUnits = lalHertzUnit;
  cFrequencySeries.deltaF = 1;
  cFrequencySeries.epoch = t;
  cFrequencySeries.f0 = 5;

  cFrequencySeries.data = cSequenceIn;

  LALCPrintFrequencySeries(&cFrequencySeries, "cFSInput.dat");
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }

  cSequenceOut = NULL;
  LALCCreateVector( &status, &cSequenceOut, READFTSERIESTEST_LEN);
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }

  cFrequencySeries2.data = cSequenceOut;

  LALCReadFrequencySeries(&status, &cFrequencySeries2, "./cFSInput.dat");
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }

  if (abs(cFrequencySeries.deltaF - cFrequencySeries.deltaF)/
      cFrequencySeries.deltaF > READFTSERIESTEST_TOL)
  {
    fprintf(stderr,"DeltaF MisMatch [ReadFTSeriesTest:%s]\n",
	    READFTSERIESTESTC_MSGEFLS);
    return READFTSERIESTESTC_EFLS;
  }
  if (strcmp(cFrequencySeries.name,cFrequencySeries2.name) != 0)
  {
    fprintf(stderr,"Name Mismatch [ReadFTSeriesTest:%s]\n",
	    READFTSERIESTESTC_MSGEFLS);
    return READFTSERIESTESTC_EFLS;
  }
  if ((cFrequencySeries.epoch.gpsSeconds) !=
      (cFrequencySeries2.epoch.gpsSeconds))
  {
    fprintf(stderr,"Epoch Second Mismatch [ReadFTSeriesTest:%s]\n",
	    READFTSERIESTESTC_MSGEFLS);
    return READFTSERIESTESTC_EFLS;
  }
  if ((cFrequencySeries.epoch.gpsNanoSeconds) !=
      (cFrequencySeries2.epoch.gpsNanoSeconds))
  {
    fprintf(stderr,"Epoch NanosecondMismatch [ReadFTSeriesTest:%s]\n",
	    READFTSERIESTESTC_MSGEFLS);
    return READFTSERIESTESTC_EFLS;
  }
  if ((cFrequencySeries.f0) ?
      (fabs(cFrequencySeries.f0 - cFrequencySeries2.f0)/cFrequencySeries.f0) :
      (abs(cFrequencySeries.f0 - cFrequencySeries2.f0)  >
      READFTSERIESTEST_TOL))
  {
    fprintf(stderr,"f0 Mismatch [ReadFTSeriesTest:%s]\n",
	    READFTSERIESTESTC_MSGEFLS);
    return READFTSERIESTESTC_EFLS;
  }
  inputUnits.unitOne = &cFrequencySeries.sampleUnits;
  inputUnits.unitTwo = &cFrequencySeries2.sampleUnits;
  LALUnitCompare(&status,&unitComp,&inputUnits);
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }

  if (!unitComp)
  {
    fprintf(stderr,"Units Mismatch [ReadFTSeriesTest:%s]\n",
	    READFTSERIESTESTC_MSGEFLS);
    return READFTSERIESTESTC_EFLS;
  }

  for (j = 0; j < cSequenceIn->length;j++)
  {
    if ((crealf(cSequenceIn->data[j]) ?
	 fabs((crealf(cSequenceIn->data[j]) - crealf(cSequenceOut->data[j]))
	      /crealf(cSequenceIn->data[j]))
	 :fabs(crealf(cSequenceIn->data[j]) - crealf(cSequenceOut->data[j]))) >
	 READFTSERIESTEST_TOL)
    {
      fprintf(stderr,"Data Tolerance Exceeded [ReadFTSeriesTest:%s]\n",
	      READFTSERIESTESTC_MSGEFLS);
      return READFTSERIESTESTC_EFLS;
    }
    if ((cimagf(cSequenceIn->data[j]) ?
	 fabs((cimagf(cSequenceIn->data[j]) - cimagf(cSequenceOut->data[j]))
	      /cimagf(cSequenceIn->data[j]))
	 :fabs(cimagf(cSequenceIn->data[j]) - cimagf(cSequenceOut->data[j]))) >
	 READFTSERIESTEST_TOL)
    {
      fprintf(stderr,"Data Tolerance Exceeded [ReadFTSeriesTest:%s]\n",
	      READFTSERIESTESTC_MSGEFLS);
      return READFTSERIESTESTC_EFLS;
    }
  }

  fprintf(stderr,"PASS\n");

  fprintf(stderr,"Testing Print/Read COMPLEX16FrequencySeries ... ");

  /* Test case 2 */

  t.gpsSeconds = 45678;
  t.gpsNanoSeconds = 89065834;

  zSequenceIn = NULL;
  LALZCreateVector( &status, &zSequenceIn, READFTSERIESTEST_LEN);
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }


  for ( i=1, zData=zSequenceIn->data; i<=READFTSERIESTEST_LEN ; i++, zData++ )
  {
    *(zData) = crect( i/4.0, (i+1)/5.0 );
  }
  zFrequencySeries.sampleUnits = lalDimensionlessUnit;
  strncpy(zFrequencySeries.name,"Complex frequency series",LALNameLength);
  zFrequencySeries.deltaF = 1.3;
  zFrequencySeries.epoch = t;
  zFrequencySeries.f0 = 0;
  zFrequencySeries.data = zSequenceIn;

  LALZPrintFrequencySeries(&zFrequencySeries, "zFSInput.dat");
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }

  zSequenceOut = NULL;
  LALZCreateVector( &status, &zSequenceOut, READFTSERIESTEST_LEN);
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }

  zFrequencySeries2.data = zSequenceOut;

  LALZReadFrequencySeries(&status, &zFrequencySeries2, "./zFSInput.dat");
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }

  if ((zFrequencySeries.deltaF) != (zFrequencySeries2.deltaF))
  {
    fprintf(stderr,"DeltaF Mismatch [ReadFTSeriesTest:%s]\n",
	    READFTSERIESTESTC_MSGEFLS);
    return READFTSERIESTESTC_EFLS;
  }
  if (strcmp(zFrequencySeries.name,zFrequencySeries2.name) != 0)
  {
    fprintf(stderr,"Name Mismatch [ReadFTSeriesTest:%s]\n",
	    READFTSERIESTESTC_MSGEFLS);
    return READFTSERIESTESTC_EFLS;
  }
  if ((zFrequencySeries.epoch.gpsSeconds) !=
      (zFrequencySeries2.epoch.gpsSeconds))
  {
    fprintf(stderr,"Epoch Seconds Mismatch [ReadFTSeriesTest:%s]\n",
	    READFTSERIESTESTC_MSGEFLS);
    return READFTSERIESTESTC_EFLS;
  }
  if ((zFrequencySeries.epoch.gpsNanoSeconds) !=
      (zFrequencySeries2.epoch.gpsNanoSeconds))
  {
    fprintf(stderr,"Epoch NanoSeconds Mismatch [ReadFTSeriesTest:%s]\n",
	    READFTSERIESTESTC_MSGEFLS);
    return READFTSERIESTESTC_EFLS;
  }
  if (zFrequencySeries.f0 ?
       (fabs(zFrequencySeries.f0 - zFrequencySeries2.f0)/zFrequencySeries.f0)
       : (fabs(zFrequencySeries.f0 - zFrequencySeries2.f0)) >
       READFTSERIESTEST_TOL)
  {
    fprintf(stderr,"f0 Mismatch [ReadFTSeriesTest:%s]\n",
	    READFTSERIESTESTC_MSGEFLS);
    return READFTSERIESTESTC_EFLS;
  }
  inputUnits.unitOne = &zFrequencySeries.sampleUnits;
  inputUnits.unitTwo = &zFrequencySeries2.sampleUnits;
  LALUnitCompare(&status,&unitComp,&inputUnits);
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }

  if (!unitComp)
  {
    fprintf(stderr,"Units Mismatch [ReadFTSeriesTest:%s]\n",
	    READFTSERIESTESTC_MSGEFLS);
    return READFTSERIESTESTC_EFLS;
  }

  for (j = 0; j < zSequenceIn->length;j++)
  {
    if ((creal(zSequenceIn->data[j]) ?
	fabs((creal(zSequenceIn->data[j]) - creal(zSequenceOut->data[j]))
	     /creal(zSequenceIn->data[j])) :
	fabs(creal(zSequenceIn->data[j]) - creal(zSequenceOut->data[j]))) >
	READFTSERIESTEST_TOL)
    {
      fprintf(stderr,"Data Tolerance Exceeded [ReadFTSeriesTest:%s]\n",
	      READFTSERIESTESTC_MSGEFLS);
      return READFTSERIESTESTC_EFLS;
    }
   if ((cimag(zSequenceIn->data[j]) ?
	fabs((cimag(zSequenceIn->data[j]) - cimag(zSequenceOut->data[j]))
	     /cimag(zSequenceIn->data[j])) :
	fabs(cimag(zSequenceIn->data[j]) - cimag(zSequenceOut->data[j]))) >
	READFTSERIESTEST_TOL)
    {
      fprintf(stderr,"Data Tolerance Exceeded [ReadFTSeriesTest:%s]\n",
	      READFTSERIESTESTC_MSGEFLS);
      return READFTSERIESTESTC_EFLS;
    }
  }

  fprintf(stderr,"PASS\n");

  fprintf(stderr,"Testing Print/Read REAL8FrequencySeries ... ");

  /* Test case 3 */
  t.gpsSeconds = 45678;
  t.gpsNanoSeconds = 89065834;

  dSequenceIn = NULL;
  LALDCreateVector( &status, &dSequenceIn, READFTSERIESTEST_LEN);
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }

  for ( i=1, dData=dSequenceIn->data; i<=READFTSERIESTEST_LEN ; i++, dData++ )
  {
    *(dData) = 0.005;
  }

  strncpy(dFrequencySeries.name,"Complex frequency series",LALNameLength);
  /* set units */
  raise.numerator = -2;
  raise.denominatorMinusOne = 0;
  LALUnitRaise(&status, &strainToMinus2, &lalStrainUnit, &raise);
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }

  LALUnitRaise(&status, &adcToMinus2, &lalADCCountUnit, &raise);
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }

  unitPair.unitOne = &strainToMinus2;
  unitPair.unitTwo = &adcToMinus2;
  LALUnitMultiply(&status, &(adcStrain), &unitPair);
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }

  unitPair.unitOne = &adcStrain;
  unitPair.unitTwo = &lalHertzUnit;
  LALUnitMultiply(&status, &dFrequencySeries.sampleUnits, &unitPair);
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }


  dFrequencySeries.deltaF = 1.3;
  dFrequencySeries.epoch = t;
  dFrequencySeries.f0 = 0;
  dFrequencySeries.data = dSequenceIn;
  LALDPrintFrequencySeries(&dFrequencySeries, "dFSInput.dat");
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }


  dSequenceOut = NULL;
  LALDCreateVector( &status, &dSequenceOut, READFTSERIESTEST_LEN);
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }


  dFrequencySeries2.data = dSequenceOut;

  LALDReadFrequencySeries(&status, &dFrequencySeries2, "./dFSInput.dat");
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }

  if ((dFrequencySeries.deltaF) != (dFrequencySeries.deltaF))
  {
    fprintf(stderr,"DeltaF Mismatch [ReadFTSeriesTest:%s]\n",
	    READFTSERIESTESTC_MSGEFLS);
    return READFTSERIESTESTC_EFLS;
  }
  if (strcmp(dFrequencySeries.name,dFrequencySeries2.name) != 0)
  {
    fprintf(stderr,"Name Mismatch [ReadFTSeriesTest:%s]\n",
	    READFTSERIESTESTC_MSGEFLS);
    return READFTSERIESTESTC_EFLS;
  }
  if ((dFrequencySeries.epoch.gpsSeconds)
      != (dFrequencySeries2.epoch.gpsSeconds))
  {
    fprintf(stderr,"Epoch Seconds Mismatch [ReadFTSeriesTest:%s]\n",
	    READFTSERIESTESTC_MSGEFLS);
    return READFTSERIESTESTC_EFLS;
  }
  if ((dFrequencySeries.epoch.gpsSeconds)
      != (dFrequencySeries2.epoch.gpsSeconds))
  {
    fprintf(stderr,"Epoch NanoSeconds Mismatch [ReadFTSeriesTest:%s]\n",
	    READFTSERIESTESTC_MSGEFLS);
    return READFTSERIESTESTC_EFLS;
  }
  if (dFrequencySeries.f0 ?
       (fabs(dFrequencySeries.f0 - dFrequencySeries2.f0)/dFrequencySeries.f0)
       : (fabs(dFrequencySeries.f0 - dFrequencySeries2.f0)) >
       READFTSERIESTEST_TOL)
  {
    fprintf(stderr,"f0 Mismatch [ReadFTSeriesTest:%s]\n",
	    READFTSERIESTESTC_MSGEFLS);
    return READFTSERIESTESTC_EFLS;
  }

  inputUnits.unitOne = &dFrequencySeries.sampleUnits;
  inputUnits.unitTwo = &dFrequencySeries2.sampleUnits;
  LALUnitCompare(&status,&unitComp,&inputUnits);
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }


  if (!unitComp)
  {
    fprintf(stderr,"Unit Mismatch [ReadFTSeriesTest:%s]\n",
	    READFTSERIESTESTC_MSGEFLS);
    return READFTSERIESTESTC_EFLS;
  }


  for (j = 0; j < dSequenceIn->length;j++)
  {
    if ((dSequenceIn->data[j] ?
	 fabs((dSequenceIn->data[j] - dSequenceOut->data[j])
	      /dSequenceIn->data[j])
	 :fabs(dSequenceIn->data[j] - dSequenceOut->data[j])) >
	 READFTSERIESTEST_TOL)
    {
      fprintf(stderr,"Data Tolerance Exceeded [ReadFTSeriesTest:%s]\n",
	      READFTSERIESTESTC_MSGEFLS);
      return READFTSERIESTESTC_EFLS;
    }

  }


  fprintf(stderr,"PASS\n");

  fprintf(stderr,"Testing Print/Read REAL4FrequencySeries ... ");


 /* Test case 4 */
  t.gpsSeconds = 45678;
  t.gpsNanoSeconds = 89065834;

  sSequenceIn = NULL;
  LALSCreateVector( &status, &sSequenceIn, READFTSERIESTEST_LEN);
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }

  for ( i=1, sData=sSequenceIn->data; i<=READFTSERIESTEST_LEN ; i++, sData++ )
  {
    *(sData) = 0.005;
  }

  strncpy(sFrequencySeries.name,"Complex frequency series",LALNameLength);
  /* set units */
  raise.numerator = -2;
  raise.denominatorMinusOne = 0;
  LALUnitRaise(&status, &strainToMinus2, &lalStrainUnit, &raise);
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }

  LALUnitRaise(&status, &adcToMinus2, &lalADCCountUnit, &raise);
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }

  unitPair.unitOne = &strainToMinus2;
  unitPair.unitTwo = &adcToMinus2;
  LALUnitMultiply(&status, &(adcStrain), &unitPair);
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }

  unitPair.unitOne = &adcStrain;
  unitPair.unitTwo = &lalHertzUnit;
  LALUnitMultiply(&status, &sFrequencySeries.sampleUnits, &unitPair);
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }


  sFrequencySeries.deltaF = 1.3;
  sFrequencySeries.epoch = t;
  sFrequencySeries.f0 = 5;
  sFrequencySeries.data = sSequenceIn;

  sSequenceOut = NULL;
  LALSCreateVector( &status, &sSequenceOut, READFTSERIESTEST_LEN);
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }

  sFrequencySeries2.data = sSequenceOut;
  LALSPrintFrequencySeries(&sFrequencySeries, "sFSInput.dat");
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }

  LALSReadFrequencySeries(&status, &sFrequencySeries2, "./sFSInput.dat");
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }

  if (fabs(sFrequencySeries.deltaF - sFrequencySeries2.deltaF)
      /sFrequencySeries.deltaF > READFTSERIESTEST_TOL)
  {
    fprintf(stderr,"Deltaf Mismatch [ReadFTSeriesTest:%s]\n",
	    READFTSERIESTESTC_MSGEFLS);
    return READFTSERIESTESTC_EFLS;
  }

  if (strcmp(sFrequencySeries.name,sFrequencySeries2.name)!=0)
  {
    fprintf(stderr,"Name Mismatch [ReadFTSeriesTest:%s]\n",
	    READFTSERIESTESTC_MSGEFLS);
    return READFTSERIESTESTC_EFLS;
  }

  if ((sFrequencySeries.epoch.gpsSeconds) !=
      (sFrequencySeries2.epoch.gpsSeconds))
  {
    fprintf(stderr,"Epoch Seconds Mismatch [ReadFTSeriesTest:%s]\n",
	    READFTSERIESTESTC_MSGEFLS);
    return READFTSERIESTESTC_EFLS;
  }
  if ((sFrequencySeries.epoch.gpsNanoSeconds) !=
      (sFrequencySeries2.epoch.gpsNanoSeconds))
  {
    fprintf(stderr,"Epoch NanoSeconds Mismatch [ReadFTSeriesTest:%s]\n",
	    READFTSERIESTESTC_MSGEFLS);
    return READFTSERIESTESTC_EFLS;
  }
  if (sFrequencySeries.f0 ?
       (fabs(sFrequencySeries.f0 - sFrequencySeries2.f0)/sFrequencySeries.f0)
       : (fabs(sFrequencySeries.f0 - sFrequencySeries2.f0)) >
       READFTSERIESTEST_TOL)
  {
    fprintf(stderr,"f0 Mismatch [ReadFTSeriesTest:%s]\n",
	    READFTSERIESTESTC_MSGEFLS);
    return READFTSERIESTESTC_EFLS;
  }

  inputUnits.unitOne = &sFrequencySeries.sampleUnits;
  inputUnits.unitTwo = &sFrequencySeries2.sampleUnits;
  LALUnitCompare(&status,&unitComp,&inputUnits);
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }

  if (!unitComp)
  {
    fprintf(stderr,"Unit Mismatch [ReadFTSeriesTest:%s]\n",
	    READFTSERIESTESTC_MSGEFLS);
    return READFTSERIESTESTC_EFLS;
  }

  for (j = 0; j < sSequenceIn->length;j++)
  {
    if ((sSequenceIn->data[j] ?
	 fabs((sSequenceIn->data[j] - sSequenceOut->data[j])
	      /sSequenceIn->data[j])
	 :fabs(sSequenceIn->data[j] - sSequenceOut->data[j])) >
	 READFTSERIESTEST_TOL)
    {
      fprintf(stderr,"Data Tolerance Exceeded [ReadFTSeriesTest:%s]\n",
	      READFTSERIESTESTC_MSGEFLS);
      return READFTSERIESTESTC_EFLS;
    }
  }


  LALCDestroyVector(&status, &cSequenceIn);
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }

  LALCDestroyVector(&status, &cSequenceOut);
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }

  LALZDestroyVector(&status, &zSequenceIn);
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }
  LALZDestroyVector(&status, &zSequenceOut);
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }

  LALDDestroyVector(&status, &dSequenceIn);
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }

  LALDDestroyVector(&status, &dSequenceOut);
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }

  LALSDestroyVector(&status, &sSequenceIn);
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }

  LALSDestroyVector(&status, &sSequenceOut);
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }

  fprintf(stderr,"PASS\n");


  fprintf(stderr,"Testing Print/Read REAL4TimeSeries ... ");

  /* Here is where testing for ReadTimeSeries is done */
  /* S Case ReadTimeSeries */

  raise.numerator = -2;
  raise.denominatorMinusOne = 0;
  LALUnitRaise(&status, &strainToMinus2, &lalStrainUnit, &raise);
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }

  LALUnitRaise(&status, &adcToMinus2, &lalADCCountUnit, &raise);
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }

  unitPair.unitOne = &strainToMinus2;
  unitPair.unitTwo = &adcToMinus2;
  LALUnitMultiply(&status, &(adcStrain), &unitPair);
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }

  unitPair.unitOne = &adcStrain;
  unitPair.unitTwo = &lalHertzUnit;
  LALUnitMultiply(&status, &sTimeSeries.sampleUnits, &unitPair);
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }


  t.gpsSeconds = 45678;
  t.gpsNanoSeconds = 89065834;

  sSequenceIn = NULL;
  LALSCreateVector( &status, &sSequenceIn, READFTSERIESTEST_LEN);
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }


  for ( i=1, sData=sSequenceIn->data; i<=READFTSERIESTEST_LEN ; i++, sData++ )
  {
    *(sData) = 0.005;
  }
  strncpy(sTimeSeries.name,"Real4 Time series",LALNameLength);
  sTimeSeries.deltaT = 1.3;
  sTimeSeries.epoch = t;
  sTimeSeries.data = sSequenceIn;
  sTimeSeries.f0 = 5;
  LALSPrintTimeSeries(&sTimeSeries, "sTSInput.dat");
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }

  sSequenceOut = NULL;
  LALSCreateVector( &status, &sSequenceOut, READFTSERIESTEST_LEN);
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }

  sTimeSeries2.data = sSequenceOut;

  LALSReadTimeSeries(&status, &sTimeSeries2, "./sTSInput.dat");
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }

  if (fabs(sTimeSeries.deltaT-sTimeSeries2.deltaT) /
      sTimeSeries.deltaT > READFTSERIESTEST_TOL)
  {
    fprintf(stderr,"DeltaT Mismatch [ReadFTSeriesTest:%s]\n",
	    READFTSERIESTESTC_MSGEFLS);
    return READFTSERIESTESTC_EFLS;
  }

  if (strcmp(sFrequencySeries.name,sFrequencySeries2.name) != 0)
  {
    fprintf(stderr,"Name Mismatch [ReadFTSeriesTest:%s]\n",
	    READFTSERIESTESTC_MSGEFLS);
    return READFTSERIESTESTC_EFLS;
  }

  if ((sTimeSeries.epoch.gpsSeconds) != (sTimeSeries2.epoch.gpsSeconds))
  {
    fprintf(stderr,"Epoch Seconds Mismatch [ReadFTSeriesTest:%s]\n",
	    READFTSERIESTESTC_MSGEFLS);
    return READFTSERIESTESTC_EFLS;
  }
  if ((sTimeSeries.epoch.gpsNanoSeconds) != (sTimeSeries2.epoch.gpsNanoSeconds))
  {
    fprintf(stderr,"Epoch NanoSeconds Mismatch [ReadFTSeriesTest:%s]\n",
	    READFTSERIESTESTC_MSGEFLS);
    return READFTSERIESTESTC_EFLS;
  }
  /*  printf("%f ... %f  f0 value\n",sTimeSeries.f0,sTimeSeries2.f0);*/
  if (sTimeSeries.f0 ?
       (fabs(sTimeSeries.f0 - sTimeSeries2.f0)/sTimeSeries.f0)
       : (fabs(sTimeSeries.f0 - sTimeSeries2.f0)) >
       READFTSERIESTEST_TOL)
  {
    fprintf(stderr,"f0 Mismatch [ReadFTSeriesTest:%s]\n",
	    READFTSERIESTESTC_MSGEFLS);
    return READFTSERIESTESTC_EFLS;
  }

  inputUnits.unitOne = &sTimeSeries.sampleUnits;
  inputUnits.unitTwo = &sTimeSeries2.sampleUnits;
  LALUnitCompare(&status,&unitComp,&inputUnits);
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }

  if (!unitComp)
  {
    fprintf(stderr,"Units Mismatch [ReadFTSeriesTest:%s]\n",
	    READFTSERIESTESTC_MSGEFLS);
    return READFTSERIESTESTC_EFLS;
  }
  for (j = 0; j < sSequenceIn->length;j++)
  {
    if ((sSequenceIn->data[j] ?
	 fabs((sSequenceIn->data[j] - sSequenceOut->data[j])
	      /sSequenceIn->data[j])
	 :fabs(sSequenceIn->data[j] - sSequenceOut->data[j])) >
	 READFTSERIESTEST_TOL)
    {
      fprintf(stderr,"Data Tolerance Exceeded [ReadFTSeriesTest:%s]\n",
	      READFTSERIESTESTC_MSGEFLS);
      return READFTSERIESTESTC_EFLS;
    }
  }

  fprintf(stderr,"PASS\n");

  fprintf(stderr,"Testing Print/Read COMPLEX16TimeSeries ... ");

  /* Z case ReadTimeSeries*/

  raise.numerator = -2;
  raise.denominatorMinusOne = 0;
  LALUnitRaise(&status, &strainToMinus2, &lalStrainUnit, &raise);
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }

  LALUnitRaise(&status, &adcToMinus2, &lalADCCountUnit, &raise);
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }

  unitPair.unitOne = &strainToMinus2;
  unitPair.unitTwo = &adcToMinus2;
  LALUnitMultiply(&status, &(adcStrain), &unitPair);
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }

  unitPair.unitOne = &adcStrain;
  unitPair.unitTwo = &lalHertzUnit;
  LALUnitMultiply(&status, &zTimeSeries.sampleUnits, &unitPair);
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }

  t.gpsSeconds = 45678;
  t.gpsNanoSeconds = 89065834;

  zSequenceIn = NULL;
  LALZCreateVector( &status, &zSequenceIn, READFTSERIESTEST_LEN);
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }


  for ( i=1, zData=zSequenceIn->data; i<=READFTSERIESTEST_LEN ; i++, zData++ )
  {
    *(zData) = crect( 0.005, 1 );
  }
  strncpy(zTimeSeries.name,"Complex16 Time series",LALNameLength);
  zTimeSeries.deltaT = 1.3;
  zTimeSeries.epoch = t;
  zTimeSeries.data = zSequenceIn;
  zTimeSeries.f0 = 0;
  LALZPrintTimeSeries(&zTimeSeries, "zTSInput.dat");
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }

  zSequenceOut = NULL;
  LALZCreateVector( &status, &zSequenceOut, READFTSERIESTEST_LEN);
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }

  zTimeSeries2.data = zSequenceOut;

  LALZReadTimeSeries(&status, &zTimeSeries2, "./zTSInput.dat");
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }

  if (fabs(zTimeSeries.deltaT) - (zTimeSeries2.deltaT)/zTimeSeries.deltaT >
      READFTSERIESTEST_TOL)
  {
    fprintf(stderr,"Mismatch DeltaT [ReadFTSeriesTest:%s]\n",
	    READFTSERIESTESTC_MSGEFLS);
    return READFTSERIESTESTC_EFLS;
  }


  if (strcmp(zTimeSeries.name,zTimeSeries2.name) != 0)
  {
    fprintf(stderr,"Name Mismatch [ReadFTSeriesTest:%s]\n",
	    READFTSERIESTESTC_MSGEFLS);
    return READFTSERIESTESTC_EFLS;
  }

  if ((zTimeSeries.epoch.gpsSeconds) != (zTimeSeries2.epoch.gpsSeconds))
  {
    fprintf(stderr,"Epoch Second Mismatch [ReadFTSeriesTest:%s]\n",
	    READFTSERIESTESTC_MSGEFLS);
    return READFTSERIESTESTC_EFLS;
  }
  if ( (zTimeSeries.epoch.gpsNanoSeconds)
       != (zTimeSeries2.epoch.gpsNanoSeconds) )
  {
    fprintf(stderr,"Epoch Nanosecond Mismatch [ReadFTSeriesTest:%s]\n",
	    READFTSERIESTESTC_MSGEFLS);
    return READFTSERIESTESTC_EFLS;
  }

  if (zTimeSeries.f0 ?
       (fabs(zTimeSeries.f0 - zTimeSeries2.f0)/zTimeSeries.f0)
       : (fabs(zTimeSeries.f0 - zTimeSeries2.f0)) >
       READFTSERIESTEST_TOL)
  {
    fprintf(stderr,"f0 Mismatch [ReadFTSeriesTest:%s]\n",
	    READFTSERIESTESTC_MSGEFLS);
    return READFTSERIESTESTC_EFLS;
  }

  inputUnits.unitOne = &zTimeSeries.sampleUnits;
  inputUnits.unitTwo = &zTimeSeries2.sampleUnits;
  LALUnitCompare(&status,&unitComp,&inputUnits);
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }

  if (!unitComp)
  {
    fprintf(stderr,"Unit Mismatch [ReadFTSeriesTest:%s]\n",
	    READFTSERIESTESTC_MSGEFLS);
    return READFTSERIESTESTC_EFUN;
  }

  for (j = 0; j < zSequenceIn->length;j++)
  {
    if ((creal(zSequenceIn->data[j]) ?
	 fabs((creal(zSequenceIn->data[j]) - creal(zSequenceOut->data[j]))
	      /creal(zSequenceIn->data[j]))
	 :fabs(creal(zSequenceIn->data[j]) - creal(zSequenceOut->data[j]))) >
	 READFTSERIESTEST_TOL)
    {
      fprintf(stderr,"Data Tolerance Exceeded [ReadFTSeriesTest:%s]\n",
	      READFTSERIESTESTC_MSGEFLS);
      return READFTSERIESTESTC_EFLS;
    }
    if ((cimag(zSequenceIn->data[j]) ?
	 fabs((cimag(zSequenceIn->data[j]) - cimag(zSequenceOut->data[j]))
	      /cimag(zSequenceIn->data[j]))
	 :fabs(cimag(zSequenceIn->data[j]) - cimag(zSequenceOut->data[j]))) >
	 READFTSERIESTEST_TOL)
    {
      fprintf(stderr,"Data Tolerance Exceeded [ReadFTSeriesTest:%s]\n",
	      READFTSERIESTESTC_MSGEFLS);
      return READFTSERIESTESTC_EFLS;
    }
  }

  fprintf(stderr,"PASS\n");

  fprintf(stderr,"Testing Print/Read REAL8TimeSeries ... ");

  /* D case  ReadTimeSeries*/

  raise.numerator = -2;
  raise.denominatorMinusOne = 0;
  LALUnitRaise(&status, &strainToMinus2, &lalStrainUnit, &raise);
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }
  LALUnitRaise(&status, &adcToMinus2, &lalADCCountUnit, &raise);
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }

  unitPair.unitOne = &strainToMinus2;
  unitPair.unitTwo = &adcToMinus2;
  LALUnitMultiply(&status, &(adcStrain), &unitPair);
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }

  unitPair.unitOne = &adcStrain;
  unitPair.unitTwo = &lalHertzUnit;
  LALUnitMultiply(&status, &sTimeSeries.sampleUnits, &unitPair);
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }

  strncpy(dTimeSeries.name,"REAL8 Time series",LALNameLength);
  dTimeSeries.sampleUnits = lalHertzUnit;
  t.gpsSeconds = 4578;
  t.gpsNanoSeconds = 890634;

  dSequenceIn = NULL;
  LALDCreateVector( &status, &dSequenceIn, READFTSERIESTEST_LEN);
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }


  for ( i=1, dData=dSequenceIn->data; i<=READFTSERIESTEST_LEN ; i++, dData++ )
  {
    *(dData) = 0.005;
  }

  dTimeSeries.deltaT = 1.3;
  dTimeSeries.epoch = t;
  dTimeSeries.data = dSequenceIn;
  dTimeSeries.f0 = 0;
  LALDPrintTimeSeries(&dTimeSeries, "dTSInput.dat");
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }

  dSequenceOut = NULL;
  LALDCreateVector( &status, &dSequenceOut, READFTSERIESTEST_LEN);
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }

  dTimeSeries2.data = dSequenceOut;
  LALDReadTimeSeries(&status, &dTimeSeries2, "./dTSInput.dat");
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }

  if (fabs(dTimeSeries.deltaT) - (dTimeSeries2.deltaT)/dTimeSeries.deltaT
      > READFTSERIESTEST_TOL)
  {
    fprintf(stderr,"DeltaT Mismatch [ReadFTSeriesTest:%d,%s]\n",status.statusCode,
	    status.statusDescription );
    return READFTSERIESTESTC_EFLS;
  }
  if (strcmp(dTimeSeries.name,dTimeSeries2.name) != 0)
  {
    fprintf(stderr,"Name Mismatch [ReadFTSeriesTest:%d,%s]\n",status.statusCode,
	    status.statusDescription );
    return READFTSERIESTESTC_EFLS;
  }

  if ((dTimeSeries.epoch.gpsSeconds) != (dTimeSeries2.epoch.gpsSeconds))
  {
    fprintf(stderr,"Epoch Seconds Mismatch [ReadFTSeriesTest:%s]\n",
	    READFTSERIESTESTC_MSGEFLS);
    return READFTSERIESTESTC_EFLS;
  }

  if ((dTimeSeries.epoch.gpsNanoSeconds)
      != (dTimeSeries2.epoch.gpsNanoSeconds))
  {
    fprintf(stderr,"Epoch Nanoseconds Mismatch [ReadFTSeriesTest:%s]\n",
	    READFTSERIESTESTC_MSGEFLS);
    return READFTSERIESTESTC_EFLS;
  }

  if (dTimeSeries.f0 ?
       (fabs(dTimeSeries.f0 - dTimeSeries2.f0)/dTimeSeries.f0)
       : (fabs(dTimeSeries.f0 - dTimeSeries2.f0)) >
       READFTSERIESTEST_TOL)
  {
    fprintf(stderr,"f0 Mismatch [ReadFTSeriesTest:%s]\n",
	    READFTSERIESTESTC_MSGEFLS);
    return READFTSERIESTESTC_EFLS;
  }

  inputUnits.unitOne = &dTimeSeries.sampleUnits;
  inputUnits.unitTwo = &dTimeSeries2.sampleUnits;
  LALUnitCompare(&status,&unitComp,&inputUnits);
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }

  if (!unitComp)
  {
    fprintf(stderr,"Unit Mismatch [ReadFTSeriesTest:%s]\n",
	    READFTSERIESTESTC_MSGEFLS);
    return READFTSERIESTESTC_EFLS;
  }

  for (j = 0; j < dSequenceIn->length;j++)
  {
    if ((dSequenceIn->data[j] ?
	 fabs((dSequenceIn->data[j] - dSequenceOut->data[j])
	      /dSequenceIn->data[j])
	 :fabs(dSequenceIn->data[j] - dSequenceOut->data[j])) >
	 READFTSERIESTEST_TOL)
    {
      fprintf(stderr,"Data Tolerance Exceeded [ReadFTSeriesTest:%s]\n",
	      READFTSERIESTESTC_MSGEFLS);
      return READFTSERIESTESTC_EFLS;
    }
  }

  fprintf(stderr,"PASS\n");

  fprintf(stderr,"Testing Print/Read COMPLEX8TimeSeries ... ");

  /* C case ReadTimeSeries*/

  raise.numerator = -2;
  raise.denominatorMinusOne = 0;
  LALUnitRaise(&status, &strainToMinus2, &lalStrainUnit, &raise);
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }

  LALUnitRaise(&status, &adcToMinus2, &lalADCCountUnit, &raise);
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }

  unitPair.unitOne = &strainToMinus2;
  unitPair.unitTwo = &adcToMinus2;
  LALUnitMultiply(&status, &(adcStrain), &unitPair);
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }

  unitPair.unitOne = &adcStrain;
  unitPair.unitTwo = &lalHertzUnit;
  LALUnitMultiply(&status, &cTimeSeries.sampleUnits, &unitPair);
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }

  t.gpsSeconds = 45678;
  t.gpsNanoSeconds = 89065834;

  cSequenceIn = NULL;
  LALCCreateVector( &status, &cSequenceIn, READFTSERIESTEST_LEN);
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }

  for ( i=1, cData=cSequenceIn->data; i<=READFTSERIESTEST_LEN ; i++, cData++ )
  {
    *(cData) = crectf( 0.005, 1 );
  }
  strncpy(cTimeSeries.name,"Complex8 Time series",LALNameLength);
  cTimeSeries.deltaT = 1.3;
  cTimeSeries.epoch = t;
  cTimeSeries.data = cSequenceIn;
  cTimeSeries.f0 = 0;
  cSequenceOut = NULL;
  LALCPrintTimeSeries(&cTimeSeries, "cTSInput.dat");
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }
  LALCCreateVector( &status, &cSequenceOut, READFTSERIESTEST_LEN);
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }

  cTimeSeries2.data = cSequenceOut;

  LALCReadTimeSeries(&status, &cTimeSeries2, "./cTSInput.dat");
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }

  if (fabs(cTimeSeries.deltaT - cTimeSeries2.deltaT)/cTimeSeries.deltaT
      > READFTSERIESTEST_TOL)
  {
    fprintf(stderr,"DeltaT Mismatch [ReadFTSeriesTest:%s]\n",
	    READFTSERIESTESTC_MSGEFLS);
    return READFTSERIESTESTC_EFLS;
  }

  if (strcmp(cTimeSeries.name,cTimeSeries2.name) != 0)
  {
    fprintf(stderr,"Name Mismatch [ReadFTSeriesTest:%s]\n",
	    READFTSERIESTESTC_MSGEFLS);
    return READFTSERIESTESTC_EFLS;
  }
  if ((cTimeSeries.epoch.gpsSeconds) != (cTimeSeries2.epoch.gpsSeconds))
  {
    fprintf(stderr,"Epoch Seconds Mismatch [ReadFTSeriesTest:%s]\n",
	    READFTSERIESTESTC_MSGEFLS);
    return READFTSERIESTESTC_EFLS;
  }
  if ((cTimeSeries.epoch.gpsNanoSeconds)!=(cTimeSeries2.epoch.gpsNanoSeconds))
  {
    fprintf(stderr,"Epoch Nanoseconds Mismatch [ReadFTSeriesTest:%s]\n",
	    READFTSERIESTESTC_MSGEFLS);
    return READFTSERIESTESTC_EFLS;
  }
  if (cTimeSeries.f0 ?
       (fabs(cTimeSeries.f0 - cTimeSeries2.f0)/cTimeSeries.f0)
       : (fabs(cTimeSeries.f0 - cTimeSeries2.f0)) >
       READFTSERIESTEST_TOL)
  {
    fprintf(stderr,"f0 Mismatch [ReadFTSeriesTest:%s]\n",
	    READFTSERIESTESTC_MSGEFLS);
    return READFTSERIESTESTC_EFLS;
  }

  inputUnits.unitOne = &cTimeSeries.sampleUnits;
  inputUnits.unitTwo = &cTimeSeries2.sampleUnits;
  LALUnitCompare(&status,&unitComp,&inputUnits);
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }

  if (!unitComp)
  {
    fprintf(stderr,"Units Mismatch [ReadFTSeriesTest:%s]\n",
	    READFTSERIESTESTC_MSGEFLS);
    return READFTSERIESTESTC_EFLS;
  }

  for (j = 0; j < cSequenceIn->length;j++)
  {
    if ((crealf(cSequenceIn->data[j]) ?
	 fabs((crealf(cSequenceIn->data[j]) - crealf(cSequenceOut->data[j]))
	      /crealf(cSequenceIn->data[j]))
	 :fabs(crealf(cSequenceIn->data[j]) - crealf(cSequenceOut->data[j]))) >
	 READFTSERIESTEST_TOL)
    {
      fprintf(stderr,"Data Tolerance Exceeded [ReadFTSeriesTest:%s]\n",
	      READFTSERIESTESTC_MSGEFLS);
      return READFTSERIESTESTC_EFLS;
    }
    if ((cimagf(cSequenceIn->data[j]) ?
	 fabs((cimagf(cSequenceIn->data[j]) - cimagf(cSequenceOut->data[j]))
	      /cimagf(cSequenceIn->data[j]))
	 :fabs(cimagf(cSequenceIn->data[j]) - cimagf(cSequenceOut->data[j]))) >
	 READFTSERIESTEST_TOL)
    {
      fprintf(stderr,"Data Tolerance Exceeded [ReadFTSeriesTest:%s]\n",
	      READFTSERIESTESTC_MSGEFLS);
      return READFTSERIESTESTC_EFLS;
    }
  }

  fprintf(stderr,"PASS\n");

  /* *******************Deallocate all memory****************** */
  LALCDestroyVector(&status, &cSequenceIn);
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }

  LALCDestroyVector(&status, &cSequenceOut);
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }

  LALZDestroyVector(&status, &zSequenceIn);
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }

  LALZDestroyVector(&status, &zSequenceOut);
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }

  LALDDestroyVector(&status, &dSequenceIn);
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }

  LALDDestroyVector(&status, &dSequenceOut);
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }

  LALSDestroyVector(&status, &sSequenceIn);
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }

  LALSDestroyVector(&status, &sSequenceOut);
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [ReadFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, READFTSERIESTESTC_MSGEFUN);
    return READFTSERIESTESTC_EFUN;
  }

  LALCheckMemoryLeaks();

  fprintf(stderr,"ReadFTSeries passed all tests.\n");

  return READFTSERIESTESTC_ENOM;

}
/** \endcond */
