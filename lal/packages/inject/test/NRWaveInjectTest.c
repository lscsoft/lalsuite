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
 * \ingroup NRWaveIO
 * \author S.Fairhurst, B.Krishnan, L.Santamaria
 *
 * \brief Test-code for NRWaveIO
 *
 * $Id$
 *
 */

/*---------- INCLUDES ----------*/
#include <config.h>
#include <lal/NRWaveIO.h>
#include <lal/NRWaveInject.h>

NRCSID (NRWAVEIOTESTC, "$Id$");

/*---------- DEFINES ----------*/

/** \name Error codes */
/*@{*/
#define NRWAVEINJECTTESTC_ENORM 	0
#define NRWAVEINJECTTESTC_ESUB  	1
#define NRWAVEINJECTTESTC_EARG  	2
#define NRWAVEINJECTTESTC_EBAD  	3
#define NRWAVEINJECTTESTC_EFILE 	4
#define NRWAVEINJECTTESTC_ESFTDIFF 5

#define NRWAVEINJECTTESTC_MSGENORM "Normal exit"
#define NRWAVEINJECTTESTC_MSGESUB  "Subroutine failed"
#define NRWAVEINJECTTESTC_MSGEARG  "Error parsing arguments"
#define NRWAVEINJECTTESTC_MSGEBAD  "Bad argument values"
#define NRWAVEINJECTTESTC_MSGEFILE "Could not create output file"
#define NRWAVEINJECTTESTC_MSGESFTDIFF "initial and final SFTs differ"
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
              __LINE__, NRWAVEINJECTTESTC, statement ? statement :  \
                   "", (msg) );                                      \
} while (0)

#define INFO( statement )                                            \
do {                                                                 \
  if ( lalDebugLevel & LALINFO )                                     \
    LALPrintError( "Info[0]: program %s, file %s, line %d, %s\n"     \
                   "        %s\n", *argv, __FILE__, __LINE__,        \
              NRWAVEINJECTTESTC, (statement) );                     \
} while (0)

#define SUB( func, statusptr )                                       \
do {                                                                 \
  if ( (func), (statusptr)->statusCode ) {                           \
    ERROR( NRWAVEINJECTTESTC_ESUB, NRWAVEINJECTTESTC_MSGESUB,      \
           "Function call \"" #func "\" failed:" );                  \
    return NRWAVEINJECTTESTC_ESUB;                                  \
  }                                                                  \
} while (0)


#define SHOULD_FAIL( func, statusptr )							\
do { 											\
  if ( func, ! (statusptr)->statusCode ) {						\
    ERROR( NRWAVEINJECTTESTC_ESUB, NRWAVEINJECTTESTC_MSGESUB,      				\
          "Function call '" #func "' should have failed for this SFT but didn't!\n");	\
    return NRWAVEINJECTTESTC_ESUB;   			                               	\
   }											\
} while(0)

#define SHOULD_FAIL_WITH_CODE( func, statusptr, code )					\
do { 											\
  if ( func, (statusptr)->statusCode != code) {						\
    LALPrintError( "Function call '" #func "' should have failed with code " #code ", but returned %d instead.\n",	\
		   (statusptr)->statusCode );						\
    return NRWAVEINJECTTESTC_ESUB;   			                               	\
   }											\
} while(0)


#define SHOULD_WORK( func, statusptr )							\
do { 											\
  if ( func, (statusptr)->statusCode ) {						\
    ERROR( NRWAVEINJECTTESTC_ESUB, NRWAVEINJECTTESTC_MSGESUB,      				\
          "Function call '" #func "' failed but should have worked for this SFT!");	\
    return NRWAVEINJECTTESTC_ESUB;   			                               	\
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
  CHAR *str=NULL;
  NRWaveCatalog nrcatalog;

  SimInspiralTable  inj;

  /* test config file reading */
  sprintf(filename, "NRWaveIOTest.cfg");
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


  /* choose some mass values */
  inj.mass1 = 10;
  inj.mass2 = 20;

  str = XLALFindNRFile( &nrCatalog, &inj, 2, 2);

  fprintf(stdout, "Closest waveform in mass ratio is %s\n", str);

  LALFree(nrcatalog.data);

  LALCheckMemoryLeaks();

  INFO( NRWAVEINJECTTESTC_MSGENORM );
  return NRWAVEINJECTTESTC_ENORM;
}
