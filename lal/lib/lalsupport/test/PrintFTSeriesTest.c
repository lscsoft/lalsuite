/*
*  Copyright (C) 2007 Jolien Creighton, John Whelan
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
   \file
   \ingroup PrintFTSeries_h
   \author J. T. Whelan <jtwhelan@loyno.edu>

   \brief Tests the routines in \ref PrintTimeSeries_c and \ref PrintFrequencySeries_c.

   \heading{Usage}
   \code
PrintFTSeriesTest
   \endcode

   \heading{Description}

This program generates and prints a sequence of frequency and time
series; the program only detects errors coming from other LAL
functions, so more in-depth testing requires  examinaton of
   the output files.  (The program \c ReadFTSeriesTest also tests
   the routines in \ref PrintFrequencySeries_c and
   \ref ReadFrequencySeries_c.)

   \heading{Notes}

The program as written generates and prints single and double
precision real and complex time and frequency series.  The routines
for integers are not tested.
*/

/**\name Error Codes */ /*@{*/
#define PRINTFTSERIESTESTC_ENOM 0       /**< Nominal exit */
#define PRINTFTSERIESTESTC_EFUN 1       /**< Error from LAL function */
/*@}*/

/** \cond DONT_DOXYGEN */

#define PRINTFTSERIESTESTC_MSGENOM "Nominal exit"
#define PRINTFTSERIESTESTC_MSGEFUN "Error from LAL function"


#include <complex.h>
#include <lal/Units.h>
#include <lal/PrintFTSeries.h>
#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>


int main( void )
{
  static LALStatus status;

  COMPLEX16Sequence   *zSequence;
  COMPLEX8Sequence    *cSequence;
  REAL8Sequence       *dSequence;
  REAL4Sequence       *sSequence;

  COMPLEX16TimeSeries   zTimeSeries;
  COMPLEX8TimeSeries    cTimeSeries;
  REAL8TimeSeries       dTimeSeries;
  REAL4TimeSeries       sTimeSeries;

  COMPLEX16FrequencySeries   zFrequencySeries;
  COMPLEX8FrequencySeries    cFrequencySeries;
  REAL8FrequencySeries       dFrequencySeries;
  REAL4FrequencySeries       sFrequencySeries;

  COMPLEX16        *z;
  COMPLEX8         *c;
  REAL8            *d;
  REAL4            *s;

  INT2             n;
  LIGOTimeGPS      t80;
  LIGOTimeGPS      t00;
  LIGOTimeGPS      t10;
  UINT4            i;

  t80.gpsSeconds = 0;
  t80.gpsNanoSeconds = 0;

  t00.gpsSeconds = 3600 * 24 * (15 * 365 + 5 * 366);
  t00.gpsNanoSeconds = 0;

  t10.gpsSeconds = 3600 * 24 * (22 * 365 + 8 * 366);
  t10.gpsNanoSeconds = 0;

  fprintf(stderr,"Printing COMPLEX16TimeSeries to zTS.dat\n");

  zSequence = NULL;

  LALZCreateVector( &status, &zSequence, 8 );
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [PrintFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, PRINTFTSERIESTESTC_MSGEFUN);
    return PRINTFTSERIESTESTC_EFUN;
  }
  for ( n=zSequence->length, z=zSequence->data; n > 0 ; --n, ++z ) {
    *z = sinh(90.0*(4-n));
    *z += - I / (1e-300 + creal(*z));
  }
  strncpy(zTimeSeries.name,"Complex time series",LALNameLength);
  zTimeSeries.sampleUnits = lalDimensionlessUnit;
  zTimeSeries.deltaT = 1e-3;
  zTimeSeries.f0 = 0;
  zTimeSeries.epoch = t80;
  zTimeSeries.data = zSequence;

  LALZPrintTimeSeries(&zTimeSeries, "zTS.dat");
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [PrintFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, PRINTFTSERIESTESTC_MSGEFUN);
    return PRINTFTSERIESTESTC_EFUN;
  }

  fprintf(stderr,"Printing COMPLEX16FrequencySeries to zFS.dat\n");

  strncpy(zFrequencySeries.name,"Complex frequency series",LALNameLength);
  zFrequencySeries.sampleUnits = lalDimensionlessUnit;
  zFrequencySeries.deltaF = 1;
  zFrequencySeries.epoch = t80;
  zFrequencySeries.f0 = 0;
  zFrequencySeries.data = zSequence;

  LALZPrintFrequencySeries(&zFrequencySeries, "zFS.dat");
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [PrintFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, PRINTFTSERIESTESTC_MSGEFUN);
    return PRINTFTSERIESTESTC_EFUN;
  }

  LALZDestroyVector( &status, &zSequence );
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [PrintFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, PRINTFTSERIESTESTC_MSGEFUN);
    return PRINTFTSERIESTESTC_EFUN;
  }

  cSequence = NULL;

  fprintf(stderr,"Printing COMPLEX8TimeSeries to cTS.dat\n");

  LALCCreateVector( &status, &cSequence, 8 );
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [PrintFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, PRINTFTSERIESTESTC_MSGEFUN);
    return PRINTFTSERIESTESTC_EFUN;
  }
  for ( n=cSequence->length, c=cSequence->data; n > 0 ; --n, ++c ) {
    *c = sinh(9.0*(4-n));
    *c += - I / (1e-30 + creal(*c));
  }
  strncpy(cTimeSeries.name,"Complex time series",LALNameLength);
  cTimeSeries.sampleUnits = lalDimensionlessUnit;
  cTimeSeries.deltaT = 1.0/1024.0;
  cTimeSeries.f0 = 0;
  cTimeSeries.epoch = t00;
  cTimeSeries.data = cSequence;

  LALCPrintTimeSeries(&cTimeSeries, "cTS.dat");
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [PrintFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, PRINTFTSERIESTESTC_MSGEFUN);
    return PRINTFTSERIESTESTC_EFUN;
  }

  fprintf(stderr,"Printing COMPLEX8FrequencySeries to cFS.dat\n");

  strncpy(cFrequencySeries.name,"Complex frequency series",LALNameLength);
  cFrequencySeries.sampleUnits = lalDimensionlessUnit;
  cFrequencySeries.deltaF = 1;
  cFrequencySeries.epoch = t80;
  cFrequencySeries.f0 = 0;
  cFrequencySeries.data = cSequence;

  LALCPrintFrequencySeries(&cFrequencySeries, "cFS.dat");
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [PrintFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, PRINTFTSERIESTESTC_MSGEFUN);
    return PRINTFTSERIESTESTC_EFUN;
  }

  fprintf(stderr,"Printing REAL8TimeSeries to dTS.dat\n");

  dSequence = NULL;

  LALDCreateVector( &status, &dSequence, 8 );
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [PrintFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, PRINTFTSERIESTESTC_MSGEFUN);
    return PRINTFTSERIESTESTC_EFUN;
  }
  for ( n=dSequence->length, d=dSequence->data; n > 0 ; --n, ++d ) {
    *d = sinh(90.0*(4-n));
  }
  strncpy(dTimeSeries.name,"Real time series",LALNameLength);
  dTimeSeries.sampleUnits = lalDimensionlessUnit;
  dTimeSeries.f0 = 0;
  dTimeSeries.deltaT = 1.0/1024.0;
  dTimeSeries.epoch = t00;
  dTimeSeries.data = dSequence;

  LALDPrintTimeSeries(&dTimeSeries, "dTS.dat");
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [PrintFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, PRINTFTSERIESTESTC_MSGEFUN);
    return PRINTFTSERIESTESTC_EFUN;
  }

  fprintf(stderr,"Printing REAL8FrequencySeries to dFS.dat\n");

/*    for ( n=dSequence->length, d=dSequence->data; n > 0 ; --n, ++d ) { */
/*      *d = 1 / (1e-300 + *d); */
/*    } */
  for ( n=dSequence->length, d=dSequence->data; n > 0 ; --n, ++d ) {
    *d = 1 / (1e-300 + *d);
  }
  strncpy(dFrequencySeries.name,"Real frequency series",LALNameLength);
  dFrequencySeries.sampleUnits = lalDimensionlessUnit;
  dFrequencySeries.f0 = 0     ;
  dFrequencySeries.deltaF = 128;
  dFrequencySeries.epoch = t00;
  dFrequencySeries.data = dSequence;
  LALDPrintFrequencySeries(&dFrequencySeries, "dFS.dat");
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [PrintFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, PRINTFTSERIESTESTC_MSGEFUN);
    return PRINTFTSERIESTESTC_EFUN;
  }

  LALDDestroyVector( &status, &dSequence );
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [PrintFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, PRINTFTSERIESTESTC_MSGEFUN);
    return PRINTFTSERIESTESTC_EFUN;
  }

  fprintf(stderr,"Printing REAL4TimeSeries to sFS.dat\n");

  sSequence = NULL;

  LALSCreateVector( &status, &sSequence, 8 );
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [PrintFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, PRINTFTSERIESTESTC_MSGEFUN);
    return PRINTFTSERIESTESTC_EFUN;
  }
  for ( n=sSequence->length, s=sSequence->data; n > 0 ; --n, ++s ) {
    *s = sinh(9.0*(4-n));
  }
  strncpy(sFrequencySeries.name,"Real time series",LALNameLength);
  sTimeSeries.sampleUnits = lalDimensionlessUnit;
  sTimeSeries.deltaT = 1.0/1024.0;
  sTimeSeries.f0 = 0;
  sTimeSeries.epoch = t10;
  sTimeSeries.data = sSequence;
  LALSPrintTimeSeries(&sTimeSeries, "sTS.dat");
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [PrintFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, PRINTFTSERIESTESTC_MSGEFUN);
    return PRINTFTSERIESTESTC_EFUN;
  }

  fprintf(stderr,"Printing REAL4FrequencySeries to sFS.dat\n");

  for ( n=sSequence->length, s=sSequence->data; n > 0 ; --n, ++s ) {
    *s = 1 / (1e-30 + *s);
  }
  strncpy(sFrequencySeries.name,"Real frequency series",LALNameLength);
  sFrequencySeries.sampleUnits = lalDimensionlessUnit;
  sFrequencySeries.f0 = 0;
  sFrequencySeries.deltaF = 128;
  sFrequencySeries.epoch = t10;
  sFrequencySeries.data = sSequence;
  LALSPrintFrequencySeries(&sFrequencySeries, "sFS.dat");
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [PrintFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, PRINTFTSERIESTESTC_MSGEFUN);
    return PRINTFTSERIESTESTC_EFUN;
  }

  LALSDestroyVector( &status, &sSequence );
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [PrintFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, PRINTFTSERIESTESTC_MSGEFUN);
    return PRINTFTSERIESTESTC_EFUN;
  }

  fprintf(stderr,"Printing heterodyned REAL8FrequencySeries to hFSe.dat\n");

  strncpy(cFrequencySeries.name,"Heterodyned frequency series",LALNameLength);
  cFrequencySeries.sampleUnits = lalDimensionlessUnit;
  cFrequencySeries.f0 = 500.0;
  cFrequencySeries.deltaF = 50.0;
  cFrequencySeries.epoch = t00;
  cFrequencySeries.data = cSequence;

  for (i=0; i<cSequence->length; ++i) {
    cSequence->data[i] = 1.0*i;
    cSequence->data[i] += I * (cFrequencySeries.f0 + cFrequencySeries.deltaF * i);
  }

  LALCPrintFrequencySeries(&cFrequencySeries, "hFSe.dat");
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [PrintFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, PRINTFTSERIESTESTC_MSGEFUN);
    return PRINTFTSERIESTESTC_EFUN;
  }

  LALCDestroyVector( &status, &cSequence );
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [PrintFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, PRINTFTSERIESTESTC_MSGEFUN);
    return PRINTFTSERIESTESTC_EFUN;
  }

  fprintf(stderr,"Printing heterodyned REAL8FrequencySeries to hFSo.dat\n");

  LALCCreateVector( &status, &cSequence, 9 );
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [PrintFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, PRINTFTSERIESTESTC_MSGEFUN);
    return PRINTFTSERIESTESTC_EFUN;
  }

  cFrequencySeries.sampleUnits = lalDimensionlessUnit;
  cFrequencySeries.f0 = 500.0;
  cFrequencySeries.deltaF = 50.0;
  cFrequencySeries.epoch = t00;
  cFrequencySeries.data = cSequence;

  for (i=0; i<cSequence->length; ++i) {
    cSequence->data[i] = 1.0*i;
    cSequence->data[i] += I * (cFrequencySeries.f0 + cFrequencySeries.deltaF*i);
  }

  LALCPrintFrequencySeries(&cFrequencySeries, "hFSo.dat");
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [PrintFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, PRINTFTSERIESTESTC_MSGEFUN);
    return PRINTFTSERIESTESTC_EFUN;
  }
  LALCDestroyVector( &status, &cSequence );
  if (status.statusCode != 0)
  {
    fprintf(stderr,"[%i]: %s [PrintFTSeriesTest:%s]\n",status.statusCode,
	    status.statusDescription, PRINTFTSERIESTESTC_MSGEFUN);
    return PRINTFTSERIESTESTC_EFUN;
  }

  LALCheckMemoryLeaks();

  return 0;
}
/** \endcond */
