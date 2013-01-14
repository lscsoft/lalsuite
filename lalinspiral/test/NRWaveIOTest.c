/*
 * Copyright (C) 2006 S.Fairhurst, B.Krishnan, L.Santamaria
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


/** \file NRWaveIOTest.c
 * \ingroup NRWaveIO_h
 * \author S.Fairhurst, B.Krishnan, L.Santamaria
 *
 * \brief Test-code for NRWaveIO
 *
 *
 */

/*---------- INCLUDES ----------*/
#include <config.h>
#include <lal/NRWaveIO.h>

/*---------- DEFINES ----------*/

/** \name Error codes */
/*@{*/
#define NRWAVEIOTESTC_ENORM 	0
#define NRWAVEIOTESTC_ESUB  	1
#define NRWAVEIOTESTC_EARG  	2
#define NRWAVEIOTESTC_EBAD  	3
#define NRWAVEIOTESTC_EFILE 	4
#define NRWAVEIOTESTC_ESFTDIFF 5

#define NRWAVEIOTESTC_MSGENORM "Normal exit"
#define NRWAVEIOTESTC_MSGESUB  "Subroutine failed"
#define NRWAVEIOTESTC_MSGEARG  "Error parsing arguments"
#define NRWAVEIOTESTC_MSGEBAD  "Bad argument values"
#define NRWAVEIOTESTC_MSGEFILE "Could not create output file"
#define NRWAVEIOTESTC_MSGESFTDIFF "initial and final SFTs differ"
/*@}*/


#ifndef SRCDIR
#define SRCDIR "."
#endif

#define TESTDIR SRCDIR "/"


/* Default parameters. */

/*********************************************************************/
/* Macros for printing errors & testing subroutines (from Creighton) */
/*********************************************************************/

#define ERROR( code, msg, statement )                                \
do {                                                                 \
  if ( lalDebugLevel & LALERROR )                                    \
    LALPrintError( "Error[0] %d: program %s, file %s, line %d, %s\n" \
                   "        %s %s\n", (code), *argv, __FILE__,       \
              __LINE__, "$Id$", statement ? statement :  \
                   "", (msg) );                                      \
} while (0)

#define INFO( statement )                                            \
do {                                                                 \
  if ( lalDebugLevel & LALINFO )                                     \
    LALPrintError( "Info[0]: program %s, file %s, line %d, %s\n"     \
                   "        %s\n", *argv, __FILE__, __LINE__,        \
              "$Id$", (statement) );                     \
} while (0)

#define SUB( func, statusptr )                                       \
do {                                                                 \
  if ( (func), (statusptr)->statusCode ) {                           \
    ERROR( NRWAVEIOTESTC_ESUB, NRWAVEIOTESTC_MSGESUB,      \
           "Function call \"" #func "\" failed:" );                  \
    return NRWAVEIOTESTC_ESUB;                                  \
  }                                                                  \
} while (0)


#define SHOULD_FAIL( func, statusptr )							\
do { 											\
  if ( func, ! (statusptr)->statusCode ) {						\
    ERROR( NRWAVEIOTESTC_ESUB, NRWAVEIOTESTC_MSGESUB,      				\
          "Function call '" #func "' should have failed for this SFT but didn't!\n");	\
    return NRWAVEIOTESTC_ESUB;   			                               	\
   }											\
} while(0)

#define SHOULD_FAIL_WITH_CODE( func, statusptr, code )					\
do { 											\
  if ( func, (statusptr)->statusCode != code) {						\
    LALPrintError( "Function call '" #func "' should have failed with code " #code ", but returned %d instead.\n",	\
		   (statusptr)->statusCode );						\
    return NRWAVEIOTESTC_ESUB;   			                               	\
   }											\
} while(0)


#define SHOULD_WORK( func, statusptr )							\
do { 											\
  if ( func, (statusptr)->statusCode ) {						\
    ERROR( NRWAVEIOTESTC_ESUB, NRWAVEIOTESTC_MSGESUB,      				\
          "Function call '" #func "' failed but should have worked for this SFT!");	\
    return NRWAVEIOTESTC_ESUB;   			                               	\
   }											\
} while(0)




/*---------- empty initializers ---------- */
LALStatus empty_status;

/*---------- Global variables ----------*/

INT4 lalDebugLevel = 3;

/* ----------------------------------------------------------------------*/
int main(int argc, char *argv[])
{
  LALStatus status = empty_status;

  CHAR filename[128];
  NRWaveCatalog nrcatalog;

  UINT4 k, length;
  REAL4TimeVectorSeries *nrdata=NULL;

  argc=0;

  sprintf(filename, "NRWaveIOTest.data");
  SHOULD_WORK (LALReadNRWave(&status, &nrdata, 10.0, filename), &status);

  length = nrdata->data->vectorLength;

  for (k = 0; k < length; k++) {
    fprintf(stdout, "%e  %e  %e\n", k*nrdata->deltaT, nrdata->data->data[k],
	    nrdata->data->data[length+k]);
  }

  fprintf(stdout, "%%filename=%s, deltaT=%e sec, Heterodyne Freq.=%e, length=%d \n",
	  nrdata->name, nrdata->deltaT, nrdata->f0, nrdata->data->vectorLength);

  XLALDestroyREAL4VectorSequence ( nrdata->data );
  LALFree(nrdata);


  /* test config file reading */
  sprintf(filename, "example.cfg");
  SHOULD_WORK (LALNRDataFind( &status, &nrcatalog, ".", filename), &status);
  fprintf(stdout, "config file %s contains %d waves\n", filename, nrcatalog.length);

  for ( k = 0; k < nrcatalog.length; k++) {
    fprintf(stdout,"\n\n Wave # %d:\n", k);
    fprintf(stdout, "Mass ratio = %f, spin1=(%f,%f,%f), spin2=(%f,%f, %f), l=%d, m=%d, datafile=%s\n",
	    nrcatalog.data[k].massRatio, nrcatalog.data[k].spin1[0], nrcatalog.data[k].spin1[1],
	    nrcatalog.data[k].spin1[2], nrcatalog.data[k].spin2[0], nrcatalog.data[k].spin1[1],
	    nrcatalog.data[k].spin1[2], nrcatalog.data[k].mode[0], nrcatalog.data[k].mode[1],
	    nrcatalog.data[k].filename);
  }

  SHOULD_WORK (LALReadNRWave(&status, &nrdata, 10.0, nrcatalog.data[0].filename), &status);

  for (k = 0; k < length; k++) {
    fprintf(stdout, "%e  %e  %e\n", k*nrdata->deltaT, nrdata->data->data[k],
	    nrdata->data->data[length+k]);
  }

  fprintf(stdout, "%%filename=%s, deltaT=%e sec, Heterodyne Freq.=%e, length=%d \n",
	  nrdata->name, nrdata->deltaT, nrdata->f0, nrdata->data->vectorLength);

  XLALDestroyREAL4VectorSequence ( nrdata->data );
  LALFree(nrdata);


  LALFree(nrcatalog.data);

  LALCheckMemoryLeaks();

  INFO( NRWAVEIOTESTC_MSGENORM );
  return NRWAVEIOTESTC_ENORM;
}
