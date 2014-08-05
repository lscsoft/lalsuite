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
 * File Name: TestConstructPLUT.c
 *
 * Authors: Sintes, A.M.,
 *
 *
 * History:   Created by Sintes June 7, 2001
 *            Modified by Badri Krishnan Feb 2003
 *
 *-----------------------------------------------------------------------
 */

/*
 * 1.  An author and Id block
 */

/**
 * \author Sintes, A. M.
 * \file
 * \ingroup LUT_h
 * \brief Tests the construction of the Look up Table (\c LUT)
 *
 * ### Usage ###
 *
 * \code
 * TestConstructPLUT [-d debuglevel] [-o outfile] [-f f0] [-p alpha delta] [-s patchSizeX patchSizeY]
 * \endcode
 *
 * ### Description ###
 *
 * \%TO BE CHANGED
 *
 * This program generates a patch grid, calculates the parameters needed for
 * building the \c LUT, builds the \c LUT and outputs a partial Hough map
 * derivative into a file. The sky patch is set at the south pole,
 * no spin-down parameters are assumed for the demodulation and
 * every third  peak in the spectrum is selected.
 *
 * By default, running this program with no arguments simply tests the subroutines,
 * producing an output file called <tt>OutHough.asc</tt>.  All default parameters are set from
 * <tt>\#define</tt>d constants.
 *
 * The <b>-d</b> option sets the debug level to the specified value
 * \c debuglevel.  The <b>-o</b> flag tells the program to print the partial Hough map
 * derivative  to the specified data file \c outfile.  The
 * <b>-f</b> option sets the intrinsic frequency \c f0 at which build the <tt>LUT</tt>.
 * The <b>-p</b> option sets the velocity orientation of the detector
 * \c alpha, \c delta (in radians).
 *
 * ### Uses ###
 *
 * \code
 * LALHOUGHCalcParamPLUT()
 * LALHOUGHConstructPLUT()
 * LALPrintError()
 * LALMalloc()
 * LALFree()
 * LALCheckMemoryLeaks()
 * \endcode
 *
 */


#include <lal/LUT.h>

/* Error codes and messages */

/**\name Error Codes */ /*@{*/
#define TESTCONSTRUCTPLUTC_ENORM 0
#define TESTCONSTRUCTPLUTC_ESUB  1
#define TESTCONSTRUCTPLUTC_EARG  2
#define TESTCONSTRUCTPLUTC_EBAD  3
#define TESTCONSTRUCTPLUTC_EFILE 4

#define TESTCONSTRUCTPLUTC_MSGENORM "Normal exit"
#define TESTCONSTRUCTPLUTC_MSGESUB  "Subroutine failed"
#define TESTCONSTRUCTPLUTC_MSGEARG  "Error parsing arguments"
#define TESTCONSTRUCTPLUTC_MSGEBAD  "Bad argument values"
#define TESTCONSTRUCTPLUTC_MSGEFILE "Could not create output file"
/*@}*/


/** \cond DONT_DOXYGEN */

/* Default parameters. */


#define F0 500.0          /*  frequency to build the LUT. */
#define TCOH 100000.0     /*  time baseline of coherent integration. */
#define DF    (1./TCOH)   /*  frequency  resolution. */
#define ALPHA 0.0
#define DELTA 0.0
#define MWR 1             /*.minWidthRatio */
#define FILEOUT "OutHough.asc"      /* file output */

/* Usage format string. */

#define USAGE "Usage: %s [-d debuglevel] [-o outfile] [-f f0] [-p alpha delta][-s patchSizeX patchSizeY]\n"

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
    ERROR( TESTCONSTRUCTPLUTC_ESUB, TESTCONSTRUCTPLUTC_MSGESUB,      \
           "Function call \"" #func "\" failed:" );                  \
    return TESTCONSTRUCTPLUTC_ESUB;                                  \
  }                                                                  \
} while (0)
/******************************************************************/

/* A global pointer for debugging. */
#ifndef NDEBUG
char *lalWatch;
#endif

#define PIXELFACTOR 2

/* the Hough Map derivative pixel type */
typedef CHAR HoughDT;


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
  INT4 i,j,k,binPoint;              /* Index counter, etc */
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
      } else {
        ERROR( TESTCONSTRUCTPLUTC_EARG, TESTCONSTRUCTPLUTC_MSGEARG, 0 );
        XLALPrintError( USAGE, *argv );
        return TESTCONSTRUCTPLUTC_EARG;
      }
    }
    /* Parse output file option. */
    else if ( !strcmp( argv[arg], "-o" ) ) {
      if ( argc > arg + 1 ) {
        arg++;
        fname = argv[arg++];
      } else {
        ERROR( TESTCONSTRUCTPLUTC_EARG, TESTCONSTRUCTPLUTC_MSGEARG, 0 );
        XLALPrintError( USAGE, *argv );
        return TESTCONSTRUCTPLUTC_EARG;
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
        ERROR( TESTCONSTRUCTPLUTC_EARG, TESTCONSTRUCTPLUTC_MSGEARG, 0 );
        XLALPrintError( USAGE, *argv );
        return TESTCONSTRUCTPLUTC_EARG;
      }
    }
    /* Parse velocity position options. */
    else if ( !strcmp( argv[arg], "-p" ) ) {
      if ( argc > arg + 2 ) {
        arg++;
	alpha = atof(argv[arg++]);
	delta = atof(argv[arg++]);
      } else {
        ERROR( TESTCONSTRUCTPLUTC_EARG, TESTCONSTRUCTPLUTC_MSGEARG, 0 );
        XLALPrintError( USAGE, *argv );
        return TESTCONSTRUCTPLUTC_EARG;
      }
    }
     /* Parse patch size option. */
    else if ( !strcmp( argv[arg], "-s" ) ) {
      if ( argc > arg + 2 ) {
        arg++;
	parRes.patchSkySizeX = patchSizeX = atof(argv[arg++]);
        parRes.patchSkySizeY = patchSizeY = atof(argv[arg++]);
      } else {
        ERROR( TESTCONSTRUCTPLUTC_EARG, TESTCONSTRUCTPLUTC_MSGEARG, 0 );
        XLALPrintError( USAGE, *argv );
        return TESTCONSTRUCTPLUTC_EARG;
      }
    }
    /* Unrecognized option. */
    else {
      ERROR( TESTCONSTRUCTPLUTC_EARG, TESTCONSTRUCTPLUTC_MSGEARG, 0 );
      XLALPrintError( USAGE, *argv );
      return TESTCONSTRUCTPLUTC_EARG;
    }
  } /* End of argument parsing loop. */
  /******************************************************************/

  if ( f0 < 0 ) {
    ERROR( TESTCONSTRUCTPLUTC_EBAD, TESTCONSTRUCTPLUTC_MSGEBAD, "freq<0:" );
    XLALPrintError( USAGE, *argv  );
    return TESTCONSTRUCTPLUTC_EBAD;
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

  PHMD = (HoughDT *)LALMalloc((xSide+1)*ySide*sizeof(HoughDT));

  for (i=0; i < maxNBorders; ++i){
    lut.border[i].ySide = ySide;
    lut.border[i].xPixel = (COORType *)LALMalloc(ySide*sizeof(COORType));
  }


  /******************************************************************/
  /* Case: no spins, patch at south pole   */
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
  /* calculate parameters needed for buiding the LUT */
  /******************************************************************/
  SUB( LALHOUGHCalcParamPLUT( &status, &parLut, &parSize, &parDem ),  &status );

  /******************************************************************/
  /* build the LUT */
  /******************************************************************/
  SUB( LALHOUGHConstructPLUT( &status, &lut, &patch, &parLut ), &status );

  /******************************************************************/
  /* construct  PHMD[i][j] accordingly  */
  /*******************************************************/

  /* initializing output  space */
  pointer = &( PHMD[0]);
  for ( k=0; k< (xSide+1)*ySide; ++k ){
    *pointer = 0;
    ++pointer;
  }

  binPoint = lut.iniBin + lut.nBin;

   /* just as a test examples */

  for( k= lut.iniBin; k < binPoint ; k+=2 ){
  /* this should be for plotting each two bins!
     so one can see all border with +1 or -1 */

  /* for( k= -2; k < 3 ; k+=2 ){ */
        /* now just 3 peaks selected -2,0,+2 */


    INT2 lb1,rb1,lb2,rb2; /* The border index. If zero means that */
        /* it does not intersect the patch, or nothing to clip */
    INT2 max1,min1,max2,min2;
    INT2 xindex;

    /* conversion of "peak index" (separation to the f0 frequency)
       into the bin in the LUT */

    i = k;
    if( k < 0) i = binPoint -1-k;

    /*reading the bin information */

    lb1 = lut.bin[i].leftB1;
    rb1 = lut.bin[i].rightB1;
    lb2 = lut.bin[i].leftB2;
    rb2 = lut.bin[i].rightB2;

    max1 = lut.bin[i].piece1max;
    min1 = lut.bin[i].piece1min;
    max2 = lut.bin[i].piece2max;
    min2 = lut.bin[i].piece2min;

    /* drawing the annuli borders */
    if(lb1){
      for(j = lut.border[lb1].yLower;
	  j <= lut.border[lb1].yUpper; j++ ){
	xindex =  lut.border[lb1].xPixel[j];
	PHMD[j*(xSide+1) + xindex] += 1;
      }
    }
    if(lb2){
      for(j = lut.border[lb2].yLower;
	  j <= lut.border[lb2].yUpper; j++ ){
	xindex =  lut.border[lb2].xPixel[j];
	PHMD[j*(xSide+1) + xindex] += 1;
      }
    }
    if(rb1){
      for(j = lut.border[rb1].yLower;
	  j <= lut.border[rb1].yUpper; j++ ){
	xindex =  lut.border[rb1].xPixel[j];
	PHMD[j*(xSide+1) + xindex] -= 1;
      }
    }
    if(rb2){
      for(j = lut.border[rb2].yLower;
	  j <= lut.border[rb2].yUpper; j++ ){
	xindex =  lut.border[rb2].xPixel[j];
	PHMD[j*(xSide+1) + xindex] -= 1;
      }
    }

    /* correcting border effects */

    /* note: if max1<min1, nothing should be done! */
    for(j=min1; j<=max1; ++j){ PHMD[j*(xSide+1)] += 1; }

    for(j=min2; j<=max2; ++j){ PHMD[j*(xSide+1)] += 1; }

  }


  /*******************************************************/
  /* printing the results into a particular file         */
  /* if the -o option was given, or into  FILEOUT        */
  /*******************************************************/

  if ( fname ) {
    fp = fopen( fname, "w" );
  } else {
    fp = fopen( FILEOUT , "w" );
  }

  if ( !fp ){
    ERROR( TESTCONSTRUCTPLUTC_EFILE, TESTCONSTRUCTPLUTC_MSGEFILE, 0 );
    return TESTCONSTRUCTPLUTC_EFILE;
  }


  for(j=ySide-1; j>=0; --j){
    for(i=0;i<xSide;++i){
      fprintf( fp ," %d",  PHMD[j*(xSide+1) + i] );
      fflush( fp );
    }
    fprintf( fp ," \n");
    fflush( fp );
  }

  fclose( fp );

  /******************************************************************/
  /* Free memory and exit */
  /******************************************************************/

  for (i=0; i<maxNBorders; ++i){
    LALFree( lut.border[i].xPixel);
  }

  LALFree( lut.border);
  LALFree( lut.bin);

  LALFree( PHMD);

  LALFree( patch.xCoor);
  LALFree( patch.yCoor);

  LALCheckMemoryLeaks();

  INFO( TESTCONSTRUCTPLUTC_MSGENORM );
  return TESTCONSTRUCTPLUTC_ENORM;
}

/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */

/** \endcond */
