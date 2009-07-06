/*
*  Copyright (C) 2007 Badri Krishnan, Jolien Creighton, Alicia Sintes Olives
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
 * File Name: TestPeak2PHMD.c
 *
 * Authors: Sintes, A.M., Krishnan, B.
 *
 * Revision: $Id$
 *
 * History:   Created by Sintes June 7, 2001
 *            Modified...
 *
 *-----------------------------------------------------------------------
 */

/*
 * 1.  An author and Id block
 */

/************************************ <lalVerbatim file="TestPeak2PHMDCV">
Author: Sintes, A. M., Krishnan, B.
$Id$
************************************* </lalVerbatim> */

/*
 * 2. Commented block with the documetation of this module
 */


/* ************************************************ <lalLaTeX>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection{Program \ \texttt{TestPeak2PHMD.c}}
\label{s:TestPeak2PHMD.c}
Tests the construction of Partial-Hough-Map-Derivatives ({\sc phmd})

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Usage}
\begin{verbatim}
TestPeak2PHMD [-d debuglevel] [-o outfile] [-f f0] [-p alpha delta]
\end{verbatim}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Description}

%TO BE CHANGED

This program generates a patch grid, calculates the parameters needed for
building a {\sc lut}, builds the {\sc lut}, constructs a {\sc phmd} at a
certain frequency (shifted from the frequency at which the {\sc lut} was built),
 and outputs the {\sc phmd} into a file. The sky patch is set at the south pole,
no spin-down parameters are assumed for the demodulation and
 every third  peak in the spectrum is selected. The peak-gram frequency interval
 is large enough to ensure compatibility with the {\sc lut} and the frequency of
 the {\sc phmd}. \\

 By default, running this program with no arguments simply tests the subroutines,
producing an output file called \verb@OutHough.asc@.  All default parameters are set from
\verb@#define@d constants.\\

The \verb@-d@ option sets the debug level to the specified value
\verb@debuglevel@.  The \verb@-o@ flag tells the program to print the partial Hough map
derivative  to the specified data file \verb@outfile@.  The
\verb@-f@ option sets the intrinsic frequency \verb@f0@ at which build the {\sc
lut}.   The \verb@-p@ option sets the velocity orientation of the detector
\verb@alpha@, \verb@delta@ (in radians).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Exit codes}
\vspace{0.1in}
\input{TESTPEAK2PHMDCErrorTable}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Uses}
\begin{verbatim}
LALHOUGHPatchGrid()
LALHOUGHParamPLUT()
LALHOUGHConstructPLUT()
LALHOUGHPeak2PHMD()
LALPrintError()
LALMalloc()
LALFree()
LALCheckMemoryLeaks()
\end{verbatim}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\subsubsection*{Notes}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\vfill{\footnotesize\input{TestPeak2PHMDCV}}

********************************************   </lalLaTeX> */



#include <lal/PHMD.h>



NRCSID (TESTPEAK2PHMDC, "$Id$");


/* Error codes and messages */

/************** <lalErrTable file="TESTPEAK2PHMDCErrorTable"> */
#define TESTPEAK2PHMDC_ENORM 0
#define TESTPEAK2PHMDC_ESUB  1
#define TESTPEAK2PHMDC_EARG  2
#define TESTPEAK2PHMDC_EBAD  3
#define TESTPEAK2PHMDC_EFILE 4

#define TESTPEAK2PHMDC_MSGENORM "Normal exit"
#define TESTPEAK2PHMDC_MSGESUB  "Subroutine failed"
#define TESTPEAK2PHMDC_MSGEARG  "Error parsing arguments"
#define TESTPEAK2PHMDC_MSGEBAD  "Bad argument values"
#define TESTPEAK2PHMDC_MSGEFILE "Could not create output file"
/******************************************** </lalErrTable> */


/* Default parameters. */

INT4 lalDebugLevel=0;

#define F0 500.0          /*  frequency to build the LUT. */
#define TCOH 100000.0     /*  time baseline of coherent integration. */
#define DF    (1./TCOH)   /*  frequency  resolution. */
#define ALPHA 0.0
#define DELTA 0.0
#define MWR 1             /*.minWidthRatio */
#define FILEOUT "OutHough.asc"      /* file output */
#define PIXELFACTOR 2
/* Usage format string. */

#define USAGE "Usage: %s [-d debuglevel] [-o outfile] [-f f0] [-p alpha delta] [-s patchSizeX patchSizeY]\n"

/*********************************************************************/
/* Macros for printing errors & testing subroutines (from Creighton) */
/*********************************************************************/

#define ERROR( code, msg, statement )                                \
do {                                                                 \
  if ( lalDebugLevel & LALERROR )                                    \
    LALPrintError( "Error[0] %d: program %s, file %s, line %d, %s\n" \
                   "        %s %s\n", (code), *argv, __FILE__,       \
              __LINE__, TESTPEAK2PHMDC, statement ? statement :  \
                   "", (msg) );                                      \
} while (0)

#define INFO( statement )                                            \
do {                                                                 \
  if ( lalDebugLevel & LALINFO )                                     \
    LALPrintError( "Info[0]: program %s, file %s, line %d, %s\n"     \
                   "        %s\n", *argv, __FILE__, __LINE__,        \
              TESTPEAK2PHMDC, (statement) );                     \
} while (0)

#define SUB( func, statusptr )                                       \
do {                                                                 \
  if ( (func), (statusptr)->statusCode ) {                           \
    ERROR( TESTPEAK2PHMDC_ESUB, TESTPEAK2PHMDC_MSGESUB,      \
           "Function call \"" #func "\" failed:" );                  \
    return TESTPEAK2PHMDC_ESUB;                                  \
  }                                                                  \
} while (0)
/******************************************************************/

/* A global pointer for debugging. */
#ifndef NDEBUG
char *lalWatch;
#endif

/* the Hough Map derivative pixel type */
/*typedef CHAR HoughDT; */


/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvv------------------------------------ */
int main(int argc, char *argv[]){

  static LALStatus       status;  /* LALStatus pointer */
  static HOUGHptfLUT     lut;     /* the Look Up Table */
  static HOUGHPatchGrid  patch;   /* Patch description */
  static HOUGHParamPLUT  parLut;  /* parameters needed to build lut  */
  static HOUGHResolutionPar parRes;
  static HOUGHSizePar    parSize;
  static HOUGHDemodPar   parDem;  /* demodulation parameters */
  static HOUGHPeakGram   pg;
  static HOUGHphmd       phmd; /* the partial Hough map derivative */
  /* ------------------------------------------------------- */

  INT8   f0Bin;           /* freq. bin to construct LUT */
  UINT2  xSide, ySide;
  UINT2  maxNBins, maxNBorders;

  /* the Hough derivative map. The patch containing at most
     SIDEX*SIDEY pixels */
  /* HoughDT PHMD[SIDEY][SIDEX+1]; */

  HoughDT *PHMD;

  HoughDT *pointer;

  CHAR *fname = NULL;               /* The output filename */
  FILE *fp=NULL;                    /* Output file */

  INT4 arg;                         /* Argument counter */
  INT4 i,j,k;                       /* Index counter, etc */

  REAL8 f0, alpha, delta, veloMod;
  REAL8 patchSizeX, patchSizeY;

  /************************************************************/
  /* Set up the default parameters. */
  /************************************************************/

  f0 =  F0;
  f0Bin = F0*TCOH;

  parRes.f0Bin =  f0Bin;
  parRes.deltaF = DF;
  parRes.patchSkySizeX  = patchSizeX = 1.0/(TCOH*F0*VEPI);
  parRes.patchSkySizeY  = patchSizeY = 1.0/(TCOH*F0*VEPI);
  parRes.pixelFactor = PIXELFACTOR;
  parRes.pixErr = PIXERR;
  parRes.linErr = LINERR;
  parRes.vTotC = VTOT;


  parDem.deltaF = DF;
  parDem.skyPatch.alpha = 0.0;
  parDem.skyPatch.delta = -LAL_PI_2;

  alpha = ALPHA;
  delta = DELTA;
  veloMod = VTOT;


  /********************************************************/
  /* Parse argument list.  i stores the current position. */
  /********************************************************/
  arg = 1;
  while ( arg < argc ) {
    /* Parse debuglevel option. */
    if ( !strcmp( argv[arg], "-d" ) ) {
      if ( argc > arg + 1 ) {
        arg++;
        lalDebugLevel = atoi( argv[arg++] );
      } else {
        ERROR( TESTPEAK2PHMDC_EARG, TESTPEAK2PHMDC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return TESTPEAK2PHMDC_EARG;
      }
    }
    /* Parse output file option. */
    else if ( !strcmp( argv[arg], "-o" ) ) {
      if ( argc > arg + 1 ) {
        arg++;
        fname = argv[arg++];
      } else {
        ERROR( TESTPEAK2PHMDC_EARG, TESTPEAK2PHMDC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return TESTPEAK2PHMDC_EARG;
      }
    }
    /* Parse frequency option. */
    else if ( !strcmp( argv[arg], "-f" ) ) {
      if ( argc > arg + 1 ) {
        arg++;
	f0 = atof(argv[arg++]);
	f0Bin = f0*TCOH;
	parRes.f0Bin =  f0Bin;
      } else {
        ERROR( TESTPEAK2PHMDC_EARG, TESTPEAK2PHMDC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return TESTPEAK2PHMDC_EARG;
      }
    }
    /* Parse velocity position options. */
    else if ( !strcmp( argv[arg], "-p" ) ) {
      if ( argc > arg + 2 ) {
        arg++;
	alpha = atof(argv[arg++]);
	delta = atof(argv[arg++]);
      } else {
        ERROR( TESTPEAK2PHMDC_EARG, TESTPEAK2PHMDC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return TESTPEAK2PHMDC_EARG;
      }
    }
     /* Parse patch size option. */
    else if ( !strcmp( argv[arg], "-s" ) ) {
      if ( argc > arg + 2 ) {
        arg++;
	parRes.patchSkySizeX = patchSizeX = atof(argv[arg++]);
        parRes.patchSkySizeY = patchSizeY = atof(argv[arg++]);
      } else {
        ERROR( TESTPEAK2PHMDC_EARG, TESTPEAK2PHMDC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return TESTPEAK2PHMDC_EARG;
      }
    }
    /* Unrecognized option. */
    else {
      ERROR( TESTPEAK2PHMDC_EARG, TESTPEAK2PHMDC_MSGEARG, 0 );
      LALPrintError( USAGE, *argv );
      return TESTPEAK2PHMDC_EARG;
    }
  } /* End of argument parsing loop. */
  /******************************************************************/

  if ( f0 < 0 ) {
    ERROR( TESTPEAK2PHMDC_EBAD, TESTPEAK2PHMDC_MSGEBAD, "freq<0:" );
    LALPrintError( USAGE, *argv  );
    return TESTPEAK2PHMDC_EBAD;
  }

  /******************************************************************/
  /* create patch grid */
  /******************************************************************/

  SUB( LALHOUGHComputeSizePar( &status, &parSize, &parRes ),  &status );

  xSide = parSize.xSide;
  ySide = parSize.ySide;
  maxNBins = parSize.maxNBins;
  maxNBorders = parSize.maxNBorders;

  /* allocate memory based on xSide and ySide */
  patch.xSide = xSide;
  patch.ySide = ySide;

  patch.xCoor = NULL;
  patch.yCoor = NULL;
  patch.xCoor = (REAL8 *)LALMalloc(xSide*sizeof(REAL8));
  patch.yCoor = (REAL8 *)LALMalloc(ySide*sizeof(REAL8));

  SUB( LALHOUGHFillPatchGrid( &status, &patch, &parSize ), &status );

  /******************************************************************/
  /* memory allocation and settings */
  /******************************************************************/

  lut.maxNBins = maxNBins;
  lut.maxNBorders = maxNBorders;
  lut.border =
         (HOUGHBorder *)LALMalloc(maxNBorders*sizeof(HOUGHBorder));
  lut.bin =
         (HOUGHBin2Border *)LALMalloc(maxNBins*sizeof(HOUGHBin2Border));

  phmd.maxNBorders = maxNBorders;
  phmd.leftBorderP =
       (HOUGHBorder **)LALMalloc(maxNBorders*sizeof(HOUGHBorder *));
  phmd.rightBorderP =
       (HOUGHBorder **)LALMalloc(maxNBorders*sizeof(HOUGHBorder *));

  phmd.ySide = ySide;
  phmd.firstColumn = NULL;
  phmd.firstColumn = (UCHAR *)LALMalloc(ySide*sizeof(UCHAR));

  PHMD = (HoughDT *)LALMalloc((xSide+1)*ySide*sizeof(HoughDT));

  for (i=0; i<maxNBorders; ++i){
    lut.border[i].ySide = ySide;
    lut.border[i].xPixel = (COORType *)LALMalloc(ySide*sizeof(COORType));
  }

  /******************************************************************/
  /* Case: no spins, patch at south pole */
  /******************************************************************/

  parDem.veloC.x = veloMod*cos(delta)*cos(alpha);
  parDem.veloC.y = veloMod*cos(delta)*sin(alpha);
  parDem.veloC.z = veloMod*sin(delta);

  parDem.positC.x = 0.0;
  parDem.positC.y = 0.0;
  parDem.positC.z = 0.0;
  parDem.timeDiff = 0.0;
  parDem.spin.length = 0;
  parDem.spin.data = NULL;

  /******************************************************************/
  /* Frequency-bin  of the Partial Hough Map*/
  /******************************************************************/

  phmd.fBin = f0Bin + 21; /* a bit shifted from the LUT */

  /******************************************************************/
  /* A Peakgram for testing                                         */
  /******************************************************************/
  pg.deltaF = DF;
  pg.fBinIni = (phmd.fBin) - maxNBins ;
  pg.fBinFin = (phmd.fBin) + 5*maxNBins;
  pg.length = 1; /* could be much smaller */
  pg.peak = NULL;
  pg.peak = (INT4 *)LALMalloc( (pg.length) * sizeof(INT4));

  /* for (ii=0; ii< pg.length; ++ii){  pg.peak[ii] = 8*ii;  } */
  /* this peakgram is for testing */
  pg.peak[0] = maxNBins;


  /******************************************************************/
  /* calculate parameters needed for buiding the LUT */
  /******************************************************************/
  SUB( LALHOUGHParamPLUT( &status, &parLut, &parSize, &parDem ),  &status );

  /******************************************************************/
  /* build the LUT */
  /******************************************************************/
  SUB( LALHOUGHConstructPLUT( &status, &lut, &patch, &parLut ), &status );

  /******************************************************************/
  /* build a PHMD from a peakgram and LUT  */
  /******************************************************************/

  SUB( LALHOUGHPeak2PHMD( &status, &phmd, &lut, &pg ), &status );

  /******************************************************************/
  /* construct  PHMD[i][j] accordingly  */
  /*******************************************************/

  /* initializing output  space */
  pointer = &( PHMD[0]);
  for ( k=0; k< (xSide+1)*ySide; ++k ){
    *pointer = 0;
    ++pointer;
  }

  /* first column correction */
  for ( k=0; k< ySide; ++k ){
    PHMD[k*(xSide+1) +0] = phmd.firstColumn[k];
  }

  /* left borders =>  +1 */
  for (k=0; k< phmd.lengthLeft; ++k){
    INT2 xindex, yLower,yUpper;
    COORType    *xPixel;

    yLower = (*(phmd.leftBorderP[k])).yLower;
    yUpper = (*(phmd.leftBorderP[k])).yUpper;
    xPixel =  &( (*(phmd.leftBorderP[k])).xPixel[0] );

    for(j=yLower; j<=yUpper;++j){
      xindex = xPixel[j];
      PHMD[j*(xSide+1) + xindex] += 1;
    }
  }

  /* right borders =>  -1 */
  for (k=0; k< phmd.lengthRight; ++k){
    INT2 xindex, yLower,yUpper;
    COORType    *xPixel;

    yLower = (*(phmd.rightBorderP[k])).yLower;
    yUpper = (*(phmd.rightBorderP[k])).yUpper;
    xPixel =  &( (*(phmd.rightBorderP[k])).xPixel[0] );

    for(j=yLower; j<=yUpper;++j){
      xindex = xPixel[j];
      PHMD[j*(xSide+1) + xindex] -= 1;
    }
  }


  /******************************************************************/
  /* printing the results into a particular file                    */
  /* if the -o option was given, or into  FILEOUT                   */
  /******************************************************************/

  if ( fname ) {
    fp = fopen( fname, "w" );
  } else {
    fp = fopen( FILEOUT , "w" );
  }

  if ( !fp ){
    ERROR( TESTPEAK2PHMDC_EFILE, TESTPEAK2PHMDC_MSGEFILE, 0 );
    return TESTPEAK2PHMDC_EFILE;
  }


  for(j=ySide-1; j>=0; --j){
    for(i=0;i<xSide;++i){
      fprintf( fp ," %f",  PHMD[j*(xSide+1) + i]);
      fflush( fp );
    }
    fprintf( fp ," \n");
    fflush( fp );
  }

  fclose( fp );


  /******************************************************************/
  /* Free memory and exit */
  /******************************************************************/

  LALFree(pg.peak);

  for (i=0; i<maxNBorders; ++i){
    LALFree( lut.border[i].xPixel);
  }

  LALFree( lut.border);
  LALFree( lut.bin);

  LALFree( phmd.leftBorderP);
  LALFree( phmd.rightBorderP);
  LALFree( phmd.firstColumn);
  LALFree( PHMD);

  LALFree( patch.xCoor);
  LALFree( patch.yCoor);

  LALCheckMemoryLeaks();

  INFO( TESTPEAK2PHMDC_MSGENORM );
  return TESTPEAK2PHMDC_ENORM;
}

/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
