/*
*  Copyright (C) 2007 Badri Krishnan, Alicia Sintes Olives
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
 * File Name: TestNDHoughMap.c
 *
 * Authors: Sintes, A.M., Krishnan, B.
 *
 *
 * History:   Created by Sintes June 25, 2001
 *            Modified     August 6, 2001
 *            Modified by Badri Krishnan Feb 2003
 *-----------------------------------------------------------------------
 */

/*
 * 1.  An author and Id block
 */

/**
 * \author Sintes, A. M., Krishnan, B.
 * \file
 * \ingroup HoughMap_h
 * \brief Tests the construction of Hough maps.
 *
 * \heading{Program \ref TestNDHoughMap.c}
 *
 * \heading{Usage}
 * \code
 * TestNDHoughMap [-d debuglevel] [-o outfile] [-f f0] [-p alpha delta] [-s patchSizeX patchSizeY]
 * \endcode
 *
 * \heading{Description}
 *
 * \%TO BE CHANGED
 *
 * Similar to the previous ones, this program generates a patch grid, calculates
 * the parameters needed for
 * building a \c lut, and  builds the \c lut. Then, given a peak-gram
 * constructs a \c phmd at a
 * certain frequency (shifted from the frequency at which the \c lut was built).
 * The sky patch is set at the south pole,
 * no spin-down parameters are assumed for the demodulation and
 * every third  peak in the spectrum is selected. The peak-gram frequency interval
 * is large enough to ensure compatibility with the \c lut and the frequency of
 * the \c phmd.
 *
 * Moreover, this program initializes a Hough map \c ht
 * and the Hough
 * map derivative space \c hd, adds one \c phmd into the Hough map
 * derivative \c hd,
 * constructs the total Hough map \c ht by integrating the \c hd,
 * and outputs the \c ht into a file.
 *
 * By default, running this program with no arguments simply tests the subroutines,
 * producing an output file called <tt>OutHough.asc</tt>.  All default parameters are set from
 * <tt>\#define</tt>d constants.
 *
 * The <b>-d</b> option sets the debug level to the specified value
 * \c debuglevel.  The <b>-o</b> flag tells the program to print the partial Hough map
 * derivative  to the specified data file \c outfile.  The
 * <b>-f</b> option sets the intrinsic frequency \c f0 at which build the <tt>lut</tt>.
 * The <b>-p</b> option sets the velocity orientation of the detector
 * \c alpha, \c delta (in radians).
 *
 * \heading{Uses}
 * \code
 * LALNDHOUGHParamPLUT()
 * LALHOUGHConstructPLUT()
 * LALHOUGHPeak2PHMD()
 * LALHOUGHInitializeHT()
 * LALHOUGHInitializeHD()
 * LALHOUGHAddPHMD2HD()
 * LALHOUGHIntegrHD2HT()
 * LALPrintError()
 * LALMalloc()
 * LALFree()
 * LALCheckMemoryLeaks()
 * \endcode
 *
 */


#include <lal/HoughMap.h>

/* Error codes and messages */

/**\name Error Codes */ /*@{*/
#define TESTNDHOUGHMAPC_ENORM 0
#define TESTNDHOUGHMAPC_ESUB  1
#define TESTNDHOUGHMAPC_EARG  2
#define TESTNDHOUGHMAPC_EBAD  3
#define TESTNDHOUGHMAPC_EFILE 4

#define TESTNDHOUGHMAPC_MSGENORM "Normal exit"
#define TESTNDHOUGHMAPC_MSGESUB  "Subroutine failed"
#define TESTNDHOUGHMAPC_MSGEARG  "Error parsing arguments"
#define TESTNDHOUGHMAPC_MSGEBAD  "Bad argument values"
#define TESTNDHOUGHMAPC_MSGEFILE "Could not create output file"
/*@}*/

/** \cond DONT_DOXYGEN */

/* Default parameters. */


#define F0 500.0          /*  frequency to build the LUT. */
#define TCOH 1800.0     /*  time baseline of coherent integration. */
#define DF    (1./TCOH)   /*  frequency  resolution. */
#define ALPHA 0.0
#define DELTA 0.0
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
    ERROR( TESTNDHOUGHMAPC_ESUB, TESTNDHOUGHMAPC_MSGESUB,      \
           "Function call \"" #func "\" failed:" );                  \
    return TESTNDHOUGHMAPC_ESUB;                                  \
  }                                                                  \
} while (0)
/******************************************************************/

/* A global pointer for debugging. */
#ifndef NDEBUG
char *lalWatch;
#endif


/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvv------------------------------------ */
int main(int argc, char *argv[]){

  static LALStatus       status;  /* LALStatus pointer */
  static HOUGHptfLUT     lut;     /* the Look Up Table */
  static HOUGHPatchGrid  patch;   /* Patch description */
  static HOUGHParamPLUT  parLut;  /* parameters needed to build lut  */
  static HOUGHResolutionPar parRes;
  static HOUGHDemodPar   parDem;  /* demodulation parameters */
  static HOUGHSizePar    parSize;
  static HOUGHPeakGram   pg;
  static HOUGHphmd       phmd; /* the partial Hough map derivative */
  static HOUGHMapDeriv   hd;   /* the Hough map derivative */
  static HOUGHMapTotal   ht;   /* the total Hough map */
  /* ------------------------------------------------------- */

  UINT2  maxNBins, maxNBorders;

  INT8   f0Bin;           /* freq. bin to construct LUT */
  UINT2  xSide, ySide;

  CHAR *fname = NULL;               /* The output filename */
  FILE *fp=NULL;                    /* Output file */

  INT4 arg;                         /* Argument counter */
  INT4 i,j;                       /* Index counter, etc */
  UINT4 k;

  REAL8 f0, alpha, delta, veloMod;
  REAL8 patchSizeX, patchSizeY;

  /************************************************************/
  /* Set up the default parameters. */
  /************************************************************/

  ht.map = NULL;
  hd.map = NULL;

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
        ERROR( TESTNDHOUGHMAPC_EARG, TESTNDHOUGHMAPC_MSGEARG, 0 );
        XLALPrintError( USAGE, *argv );
        return TESTNDHOUGHMAPC_EARG;
      }
    }
    /* Parse output file option. */
    else if ( !strcmp( argv[arg], "-o" ) ) {
      if ( argc > arg + 1 ) {
        arg++;
        fname = argv[arg++];
      } else {
        ERROR( TESTNDHOUGHMAPC_EARG, TESTNDHOUGHMAPC_MSGEARG, 0 );
        XLALPrintError( USAGE, *argv );
        return TESTNDHOUGHMAPC_EARG;
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
        ERROR( TESTNDHOUGHMAPC_EARG, TESTNDHOUGHMAPC_MSGEARG, 0 );
        XLALPrintError( USAGE, *argv );
        return TESTNDHOUGHMAPC_EARG;
      }
    }
    /* Parse velocity position options. */
    else if ( !strcmp( argv[arg], "-p" ) ) {
      if ( argc > arg + 2 ) {
        arg++;
	alpha = atof(argv[arg++]);
	delta = atof(argv[arg++]);
      } else {
        ERROR( TESTNDHOUGHMAPC_EARG, TESTNDHOUGHMAPC_MSGEARG, 0 );
        XLALPrintError( USAGE, *argv );
        return TESTNDHOUGHMAPC_EARG;
      }
    }
     /* Parse patch size option. */
    else if ( !strcmp( argv[arg], "-s" ) ) {
      if ( argc > arg + 2 ) {
        arg++;
	parRes.patchSkySizeX = patchSizeX = atof(argv[arg++]);
        parRes.patchSkySizeY = patchSizeY = atof(argv[arg++]);
      } else {
        ERROR( TESTNDHOUGHMAPC_EARG, TESTNDHOUGHMAPC_MSGEARG, 0 );
        XLALPrintError( USAGE, *argv );
        return TESTNDHOUGHMAPC_EARG;
      }
    }
    /* Unrecognized option. */
    else {
      ERROR( TESTNDHOUGHMAPC_EARG, TESTNDHOUGHMAPC_MSGEARG, 0 );
      XLALPrintError( USAGE, *argv );
      return TESTNDHOUGHMAPC_EARG;
    }
  } /* End of argument parsing loop. */
  /******************************************************************/

  if ( f0 < 0 ) {
    ERROR( TESTNDHOUGHMAPC_EBAD, TESTNDHOUGHMAPC_MSGEBAD, "freq<0:" );
    XLALPrintError( USAGE, *argv  );
    return TESTNDHOUGHMAPC_EBAD;
  }


  /******************************************************************/
  /* create patch grid */
  /******************************************************************/

  SUB( LALHOUGHComputeNDSizePar( &status, &parSize, &parRes ),  &status );

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

  ht.xSide = xSide;
  ht.ySide = ySide;
  ht.map   = (HoughTT *)LALMalloc(xSide*ySide*sizeof(HoughTT));

  hd.xSide = xSide;
  hd.ySide = ySide;
  hd.map   = (HoughDT *)LALMalloc((xSide+1)*ySide*sizeof(HoughDT));

  phmd.ySide = ySide;
  phmd.firstColumn = NULL;
  phmd.firstColumn = (UCHAR *)LALMalloc(ySide*sizeof(UCHAR));


  for (i=0; i<maxNBorders; ++i){
    lut.border[i].ySide = ySide;
    lut.border[i].xPixel = (COORType *)LALMalloc(ySide*sizeof(COORType));
  }


  /******************************************************************/
  /* Case: no spins, patch at south pole */
  /************************************************************/

  parDem.veloC.x = veloMod*cos(delta)*cos(alpha);
  parDem.veloC.y = veloMod*cos(delta)*sin(alpha);
  parDem.veloC.z = veloMod*sin(delta);

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
  pg.fBinFin = (phmd.fBin)+ 5*maxNBins;
  pg.length = maxNBins; /* could be much smaller */
  pg.peak = NULL;
  pg.peak = (INT4 *)LALMalloc( (pg.length) * sizeof(INT4));

  for (k=0; k< pg.length; ++k){  pg.peak[k] = 3*k;  } /* a test */


  /******************************************************************/
  /* calculate parameters needed for buiding the LUT */
  /******************************************************************/
  SUB( LALNDHOUGHParamPLUT( &status, &parLut, &parSize, &parDem ),  &status );

  /******************************************************************/
  /* build the LUT */
  /******************************************************************/
  SUB( LALHOUGHConstructPLUT( &status, &lut, &patch, &parLut ), &status );

  /******************************************************************/
  /* build a PHMD from a peakgram and LUT  */
  /******************************************************************/

  SUB( LALHOUGHPeak2PHMD( &status, &phmd, &lut, &pg ), &status );

  /******************************************************************/
  /* initializing the Hough map space */
  /******************************************************************/

  SUB( LALHOUGHInitializeHT( &status, &ht, &patch ), &status );

  /******************************************************************/
  /* initializing the Hough map derivative space */
  /******************************************************************/

  SUB( LALHOUGHInitializeHD( &status, &hd), &status );

  /******************************************************************/
  /* sum a partial-HMD into a HD */
  /******************************************************************/
  SUB( LALHOUGHAddPHMD2HD( &status, &hd, &phmd ), &status );

  /******************************************************************/
  /* construction of a total Hough map: integration of a HM-deriv.  */
  /******************************************************************/

  SUB( LALHOUGHIntegrHD2HT( &status, &ht, &hd ), &status );

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
    ERROR( TESTNDHOUGHMAPC_EFILE, TESTNDHOUGHMAPC_MSGEFILE, 0 );
    return TESTNDHOUGHMAPC_EFILE;
  }


  for(j=ySide-1; j>=0; --j){
    for(i=0;i<xSide;++i){
      fprintf( fp ," %f", ht.map[j*xSide +i]);
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

  LALFree( ht.map);
  LALFree( hd.map);

  LALFree( patch.xCoor);
  LALFree( patch.yCoor);

  LALCheckMemoryLeaks();

  INFO( TESTNDHOUGHMAPC_MSGENORM );
  return TESTNDHOUGHMAPC_ENORM;
}

/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */

/** \endcond */
