/*
*  Copyright (C) 2007 Jolien Creighton, Thomas Cokelaer
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
 * \author Cokelaer, T.
 * \file
 *
 * \brief Interface to generate any kind of gravitational waves signal.
 *
 * \heading{Usage}
 * \code
 * \c InterfaceTest
 * \endcode
 *
 * \heading{Description}
 * Right now this test file read an input xml (injection.xml) file and for
 * each valid line it computes the approriate waveform. Those waveforms are
 * either produce within the inject package (PPN waveform) or within the
 * inspiral package (TaylorT[1,2,3], EOB, SpinTaylor, PadeT1).
 *
 * Then the function LALGeneralInspiral create the amplitude, freq and phi
 * vectors needed for further injection including the data itself (noise, h(t))
 *
 * Finally, the injections are stored in a vector which is saved in "injection.dat"
 * file.
 *
 * \heading{Uses}
 * \code
 *
 * \endcode
 * \heading{Notes}
 */

/**\name Error Codes */ /*@{*/
#define INJECTIONINTERFACETESTC_ENORM 	0
#define INJECTIONINTERFACETESTC_ESUB  	1
#define INJECTIONINTERFACETESTC_EARG  	2
#define INJECTIONINTERFACETESTC_EVAL  	3
#define INJECTIONINTERFACETESTC_EFILE 	4
#define INJECTIONINTERFACETESTC_EMEM  	5
#define INJECTIONINTERFACETESTC_EINJECT  6

#define INJECTIONINTERFACETESTC_MSGENORM 	"Normal exit"
#define INJECTIONINTERFACETESTC_MSGESUB  	"Subroutine failed"
#define INJECTIONINTERFACETESTC_MSGEARG  	"Error parsing arguments"
#define INJECTIONINTERFACETESTC_MSGEVAL  	"Input argument out of valid range"
#define INJECTIONINTERFACETESTC_MSGEFILE 	"Could not open file"
#define INJECTIONINTERFACETESTC_MSGEMEM  	"Out of memory"
#define INJECTIONINTERFACETESTC_MSGEINJECT  	"No valid injection to do ... ? "
/*@}*/

/* --- the names of the files to be used --- */
#define INJECTIONINTERFACETEST_INJECTIONXMLFILE    "injection.xml"
#define INJECTIONINTERFACETEST_INJECTIONOUTPUTFILE "injection.dat"

/* --- include files --- */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include <lal/GenerateInspiral.h>
#include <lal/VectorOps.h>
#include <lal/SeqFactories.h>

#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/LALDatatypes.h>
#include <lal/AVFactories.h>

#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOLwXMLInspiralRead.h>

#include <lal/Date.h>
#include <lal/Units.h>
#include <lal/FindChirp.h>
#include <lal/PrintFTSeries.h>

#define ERROR( code, msg, statement )                                \
do                                                                   \
if ( lalDebugLevel & LALERROR )                                      \
{                                                                    \
  LALPrintError( "Error[0] %d: program %s, file %s, line %d, %s\n"   \
		 "        %s %s\n", (code), program, __FILE__,       \
		 __LINE__, "$Id$", statement ? statement :   \
                 "", (msg) );                                        \
}                                                                    \
while (0)

#define SUB( func, statusptr )                                       \
do                                                                   \
if ( (func), (statusptr)->statusCode )                               \
{                                                                    \
  ERROR( INJECTIONINTERFACETESTC_ESUB, INJECTIONINTERFACETESTC_MSGESUB,                \
         "Function call \"" #func "\" failed:" );                    \
  exit( INJECTIONINTERFACETESTC_ESUB );                                       \
}                                                                    \
while (0)

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

char *program;



/* --- Main program start here --- */
int main(int UNUSED argc, char **argv)
{
  /* some variables to create the data */
  UINT4 	startTime;
  UINT4		endTime;
  UINT4 	k;
  UINT4 	numPoints 	= 524288  ; /* arbitrary length of the data*/
  UINT4         numInjections   = 0; /* by default no injection. */
  REAL8 	sampling	= 2048.;
  static  LALStatus 	status;
  /* injection structure */
  SimInspiralTable    *injections = NULL;
  SimInspiralTable    *thisInj    = NULL;
  /* the data */
  REAL4TimeSeries               ts ; /* A time series to store the injection */
  COMPLEX8FrequencySeries       fs ; /* A freq series to store the psd 	*/
  FILE		*output;             /* output result*/
  CHAR *filename;
  UINT4 lenfile;

  /* --- for debugging --- */
  program       = *argv;

  /* --- Start Main part here --- */
  /* First, we test if the injection xml file exist */
  if ((output	= fopen(INJECTIONINTERFACETEST_INJECTIONXMLFILE,"r")) == NULL){
	ERROR(INJECTIONINTERFACETESTC_EFILE,INJECTIONINTERFACETESTC_MSGEFILE, 0);
	exit(0);
  }
  else {
	fclose(output); /*if it exist let's close it */
  }

  /* then let's start to open the output file */
  if ((output 	= fopen(INJECTIONINTERFACETEST_INJECTIONOUTPUTFILE,"w")) == NULL){
      ERROR(INJECTIONINTERFACETESTC_EFILE,INJECTIONINTERFACETESTC_MSGEFILE, 0);
      exit(0);
    }

  /* let's allocate memory for to create some dummy data to
     inject the waveforms */
  memset(&ts, 0, sizeof(REAL4TimeSeries));
  LALSCreateVector( &status, &(ts.data), numPoints);

  memset( &fs, 0, sizeof(COMPLEX8FrequencySeries));
  LALCCreateVector( &status, &(fs.data), numPoints / 2 + 1 );


  /* --- Let's fix some variables we have to be in agreement with the xml file data --- */
  ts.epoch.gpsSeconds 	= 729273600;		       	/* gps time of the time series		*/
  startTime 		= 729273600;		       	/* start time and end time of ..	*/
  endTime   		= startTime + 200;	       	/* ..injection; should be in agreement..*/
  						       	/* ..with the xml file			*/
  ts.sampleUnits 	= lalADCCountUnit;	       	/*  UNITY ?? 				*/
  ts.deltaT 		= 1./sampling;		       	/* sampling				*/
  ts.name[0]='H'; 					/* L, H, G, V or T for the detector
							   otherwise it is optimally oriented	*/


  fs.deltaF 		= sampling / numPoints;	       /* idem for fs				*/
  fs.sampleUnits 	= lalADCCountUnit;

  /* --- the psd is flat for simplicity --- */
  for( k = 0 ; k< numPoints/2; k++){
    fs.data->data[k] = 1;
  }
  /* --- and time series is zero --- */
  for( k = 0 ; k< numPoints; k++){
    ts.data->data[k] = 0.;
  }


  lenfile = strlen(INJECTIONINTERFACETEST_INJECTIONXMLFILE) + 1;
  filename = (CHAR *) calloc (lenfile, sizeof(CHAR *));
  memcpy( filename, INJECTIONINTERFACETEST_INJECTIONXMLFILE, lenfile);


  /* --- read injection  here --- */
  SUB(numInjections = SimInspiralTableFromLIGOLw( &injections,
						  filename,
						  startTime,
						  endTime), &status);

  /* any injection to do ? */
  if ( numInjections <= 0 )
    {
      ERROR(INJECTIONINTERFACETESTC_EINJECT, INJECTIONINTERFACETESTC_MSGEINJECT, 0);
      exit( 1 );
    }
  /* --- inject here --- */
  SUB( LALFindChirpInjectSignals( &status,
				  &ts,
				  injections,
				  &fs ), &status);


  /* --- now we save the results --- */
  for (k = 0; k < numPoints; k++){
    fprintf(output,"%15.12f %e\n",  (float)k/sampling, ts.data->data[k]);
  }
  fclose(output);

  SUB( LALSDestroyVector( &status, &(ts.data) ), &status );
  SUB( LALCDestroyVector( &status, &(fs.data) ), &status );

  /* --- and finally free memory --- */
  while ( injections )
    {
      thisInj 		= injections;
      injections 	= injections->next;
      LALFree(thisInj);
    }

  LALCheckMemoryLeaks();
  return 0;
}
