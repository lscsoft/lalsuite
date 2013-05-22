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
 * File Name: TestDriveNDHough.c
 *
 * Authors: Sintes, A.M., Krishnan, B.
 *
 *
 * History:   Created by Sintes August 7, 2001
 *            Modified by Badri Krishnan Feb 2003
 *
 *-----------------------------------------------------------------------
 */

/*
 * 1.  An author and Id block
 */

/**
\author Sintes, A. M., Krishnan, B.
\file
\ingroup LALHough_h
\brief Tests the construction

\heading{Program \ref TestDriveNDHough.c}

\heading{Usage}
\code
TestDriveNDHough [-d debuglevel] [-o outfile] [-f f0] [-p alpha delta] [-s patchSizeX patchSizeY]
\endcode


\heading{Description}

\%TO BE CHANGED

This program generates  a patch grid, a vector of  \c luts by changing the
alpha component of the velocity orientation of the detector by a fixed amount
in each of them, and a vector of
peak-grams (all of them containing the same information). Similar to the previous
test codes the patch is set at the south pole.

Then the program builds the set
of \c phmd, updates the cylinder and computes a Hough map at a given
frequency using only one horizontal line set of \c phmd, and outputs the
result into a file.


By default, running this program with no arguments simply tests the subroutines,
producing an output file called <tt>OutHough.asc</tt>.  All default parameters are set from
<tt>\#define</tt>d constants.

The <b>-d</b> option sets the debug level to the specified value
\c debuglevel.  The <b>-o</b> flag tells the program to print the partial Hough map
derivative  to the specified data file \c outfile.  The
<b>-f</b> option sets the intrinsic frequency \c f0 at which build the <tt>lut</tt>.
The <b>-p</b> option sets the velocity orientation of the detector
\c alpha, \c delta (in radians) for the first \c lut (time-stamp).

\heading{Uses}
\code
LALNDHOUGHParamPLUT()
LALHOUGHConstructPLUT()
LALHOUGHConstructSpacePHMD()
LALHOUGHupdateSpacePHMDup()
LALHOUGHInitializeHT()
LALHOUGHConstructHMT()
LALPrintError()
LALMalloc()
LALFree()
LALCheckMemoryLeaks()
\endcode

*/

#include <lal/LALHough.h>

/* Error codes and messages */

/**\name Error Codes */ /*@{*/
#define TESTDRIVENDHOUGHC_ENORM 0
#define TESTDRIVENDHOUGHC_ESUB  1
#define TESTDRIVENDHOUGHC_EARG  2
#define TESTDRIVENDHOUGHC_EBAD  3
#define TESTDRIVENDHOUGHC_EFILE 4

#define TESTDRIVENDHOUGHC_MSGENORM "Normal exit"
#define TESTDRIVENDHOUGHC_MSGESUB  "Subroutine failed"
#define TESTDRIVENDHOUGHC_MSGEARG  "Error parsing arguments"
#define TESTDRIVENDHOUGHC_MSGEBAD  "Bad argument values"
#define TESTDRIVENDHOUGHC_MSGEFILE "Could not create output file"
/*@}*/


/** \cond DONT_DOXYGEN */

/* Default parameters. */


#define F0 500.0          /*  frequency to build the LUT. */
#define TCOH 1800.0     /*  time baseline of coherent integration. */
#define DF    (1./TCOH)   /*  frequency  resolution. */
#define ALPHA 0.0
#define DELTA 0.0
#define FILEOUT "OutHough.asc"      /* file output */
#define MOBSCOH 10
#define NFSIZE  5
#define STEPALPHA 0.005
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
    ERROR( TESTDRIVENDHOUGHC_ESUB, TESTDRIVENDHOUGHC_MSGESUB,      \
           "Function call \"" #func "\" failed:" );                  \
    return TESTDRIVENDHOUGHC_ESUB;                                  \
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

  static LALStatus           status;  /* LALStatus pointer */
  static HOUGHptfLUTVector   lutV; /* the Look Up Table vector*/
  static HOUGHPeakGramVector pgV;
  static PHMDVectorSequence  phmdVS;  /* the partial Hough map derivatives */
  static UINT8FrequencyIndexVector freqInd;

  static HOUGHResolutionPar parRes;
  static HOUGHPatchGrid  patch;   /* Patch description */

  static HOUGHParamPLUT  parLut;  /* parameters needed to build lut  */
  static HOUGHDemodPar   parDem;  /* demodulation parameters */
  static HOUGHSizePar    parSize;

  static HOUGHMapTotal   ht;   /* the total Hough map */
  /* ------------------------------------------------------- */

  UINT2  maxNBins, maxNBorders;

  INT8   f0Bin;           /* freq. bin to construct LUT */
  INT8   fBin;

  UINT2 xSide, ySide;

  CHAR *fname = NULL;               /* The output filename */
  FILE *fp=NULL;                    /* Output file */

  INT4 arg;                         /* Argument counter */
  UINT4 i,j;                       /* Index counter, etc */
  INT4 k;
  REAL8 f0, alpha, delta, veloMod;
  REAL8 patchSizeX, patchSizeY;

  /************************************************************/
  /* Set up the default parameters. */
  /* **********************************************************/

  lutV.length    = MOBSCOH;
  pgV.length     = MOBSCOH;
  phmdVS.length  = MOBSCOH;
  freqInd.length = MOBSCOH;
  phmdVS.nfSize  = NFSIZE;

  freqInd.deltaF = DF;
  phmdVS.deltaF  = DF;

  lutV.lut = NULL;
  pgV.pg = NULL;
  phmdVS.phmd = NULL;
  freqInd.data = NULL;
  ht.map = NULL;

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
        ERROR( TESTDRIVENDHOUGHC_EARG, TESTDRIVENDHOUGHC_MSGEARG, 0 );
        XLALPrintError( USAGE, *argv );
        return TESTDRIVENDHOUGHC_EARG;
      }
    }
    /* Parse output file option. */
    else if ( !strcmp( argv[arg], "-o" ) ) {
      if ( argc > arg + 1 ) {
        arg++;
        fname = argv[arg++];
      } else {
        ERROR( TESTDRIVENDHOUGHC_EARG, TESTDRIVENDHOUGHC_MSGEARG, 0 );
        XLALPrintError( USAGE, *argv );
        return TESTDRIVENDHOUGHC_EARG;
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
        ERROR( TESTDRIVENDHOUGHC_EARG, TESTDRIVENDHOUGHC_MSGEARG, 0 );
        XLALPrintError( USAGE, *argv );
        return TESTDRIVENDHOUGHC_EARG;
      }
    }
    /* Parse velocity position options. */
    else if ( !strcmp( argv[arg], "-p" ) ) {
      if ( argc > arg + 2 ) {
        arg++;
	alpha = atof(argv[arg++]);
	delta = atof(argv[arg++]);
      } else {
        ERROR( TESTDRIVENDHOUGHC_EARG, TESTDRIVENDHOUGHC_MSGEARG, 0 );
        XLALPrintError( USAGE, *argv );
        return TESTDRIVENDHOUGHC_EARG;
      }
    }
    /* Parse patch size option. */
    else if ( !strcmp( argv[arg], "-s" ) ) {
      if ( argc > arg + 2 ) {
        arg++;
	parRes.patchSkySizeX = patchSizeX = atof(argv[arg++]);
        parRes.patchSkySizeY = patchSizeY = atof(argv[arg++]);
      } else {
        ERROR( TESTDRIVENDHOUGHC_EARG, TESTDRIVENDHOUGHC_MSGEARG, 0 );
        XLALPrintError( USAGE, *argv );
        return TESTDRIVENDHOUGHC_EARG;
      }
    }
    /* Unrecognized option. */
    else {
      ERROR( TESTDRIVENDHOUGHC_EARG, TESTDRIVENDHOUGHC_MSGEARG, 0 );
      XLALPrintError( USAGE, *argv );
      return TESTDRIVENDHOUGHC_EARG;
    }
  } /* End of argument parsing loop. */
  /******************************************************************/

  if ( f0 < 0 ) {
    ERROR( TESTDRIVENDHOUGHC_EBAD, TESTDRIVENDHOUGHC_MSGEBAD, "freq<0:" );
    XLALPrintError( USAGE, *argv  );
    return TESTDRIVENDHOUGHC_EBAD;
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

  /* allocate memory based on xSide and ySide */
  patch.xCoor = NULL;
  patch.yCoor = NULL;
  patch.xCoor = (REAL8 *)LALMalloc(xSide*sizeof(REAL8));
  patch.yCoor = (REAL8 *)LALMalloc(ySide*sizeof(REAL8));

  SUB( LALHOUGHFillPatchGrid( &status, &patch, &parSize ), &status );

  /******************************************************************/
  /* memory allocation and settings */
  /******************************************************************/

  lutV.lut = (HOUGHptfLUT *)LALMalloc(MOBSCOH*sizeof(HOUGHptfLUT));
  pgV.pg = (HOUGHPeakGram *)LALMalloc(MOBSCOH*sizeof(HOUGHPeakGram));
  phmdVS.phmd =(HOUGHphmd *)LALMalloc(MOBSCOH*NFSIZE*sizeof(HOUGHphmd));
  freqInd.data =  ( UINT8 *)LALMalloc(MOBSCOH*sizeof(UINT8));

  for(j=0; j<lutV.length; ++j){
    lutV.lut[j].maxNBins = maxNBins;
    lutV.lut[j].maxNBorders = maxNBorders;
    lutV.lut[j].border =
         (HOUGHBorder *)LALMalloc(maxNBorders*sizeof(HOUGHBorder));
    lutV.lut[j].bin =
         (HOUGHBin2Border *)LALMalloc(maxNBins*sizeof(HOUGHBin2Border));
  }

  for(j=0; j<phmdVS.length * phmdVS.nfSize; ++j){
    phmdVS.phmd[j].maxNBorders = maxNBorders;
    phmdVS.phmd[j].leftBorderP =
       (HOUGHBorder **)LALMalloc(maxNBorders*sizeof(HOUGHBorder *));
    phmdVS.phmd[j].rightBorderP =
       (HOUGHBorder **)LALMalloc(maxNBorders*sizeof(HOUGHBorder *));
  }


  ht.xSide = xSide;
  ht.ySide = ySide;
  ht.map   = NULL;
  ht.map   = (HoughTT *)LALMalloc(xSide*ySide*sizeof(HoughTT));

  for(j=0; j<phmdVS.length * phmdVS.nfSize; ++j){
    phmdVS.phmd[j].ySide = ySide;
    phmdVS.phmd[j].firstColumn = NULL;
    phmdVS.phmd[j].firstColumn = (UCHAR *)LALMalloc(ySide*sizeof(UCHAR));
  }

  for (j=0; j<lutV.length ; ++j){
    for (i=0; i<maxNBorders; ++i){
      lutV.lut[j].border[i].ySide = ySide;
      lutV.lut[j].border[i].xPixel =
                            (COORType *)LALMalloc(ySide*sizeof(COORType));
    }
  }


  /******************************************************************/
  /* Case: no spins, patch at south pole */
  /************************************************************/
  parDem.deltaF = DF;
  parDem.skyPatch.alpha = 0.0;
  parDem.skyPatch.delta = -LAL_PI_2;


  parDem.timeDiff = 0.0;
  parDem.spin.length = 0;
  parDem.spin.data = NULL;

           /*************************************************/
  for (j=0;j< MOBSCOH;++j){  /* create all the LUTs */
    parDem.veloC.x = veloMod*cos(delta)*cos(alpha);
    parDem.veloC.y = veloMod*cos(delta)*sin(alpha);
    parDem.veloC.z = veloMod*sin(delta);

    alpha +=  STEPALPHA; /* shift alpha several degrees */

    /* calculate parameters needed for buiding the LUT */
    SUB( LALNDHOUGHParamPLUT( &status, &parLut, &parSize, &parDem ),  &status );

    /* build the LUT */
    SUB( LALHOUGHConstructPLUT( &status, &(lutV.lut[j]), &patch, &parLut ),
	 &status );
  }


  /******************************************************************/
  /* create Peakgrams for testing                                         */
  /******************************************************************/

  fBin = f0Bin + 21;  /* a Frequency-bin  shifted from the LUT */

  for (j=0;j< MOBSCOH;++j){  /* create all the peakgrams */
   pgV.pg[j].deltaF = DF;
   pgV.pg[j].fBinIni = (fBin) - maxNBins;
   pgV.pg[j].fBinFin = (fBin)+ 5*maxNBins;
   pgV.pg[j].length = maxNBins; /* could be much smaller */
   pgV.pg[j].peak = NULL;
   pgV.pg[j].peak = (INT4 *)LALMalloc( ( pgV.pg[j].length) * sizeof(INT4));

   for (i=0; i< pgV.pg[j].length; ++i){ pgV.pg[j].peak[i] = 3*i; } /* test */
  }


  /******************************************************************/
  /* build the set of  PHMD  */
  /******************************************************************/

  phmdVS.fBinMin = fBin;
  SUB( LALHOUGHConstructSpacePHMD(&status, &phmdVS, &pgV, &lutV), &status );

  /* shift the structure one frequency bin */
  SUB( LALHOUGHupdateSpacePHMDup(&status, &phmdVS, &pgV, &lutV), &status );


  /******************************************************************/
  /* initializing the Hough map space */
  /******************************************************************/

  SUB( LALHOUGHInitializeHT( &status, &ht, &patch ), &status );


  /******************************************************************/
  /* construction of a total Hough map  */
  /******************************************************************/

  for (j=0;j< MOBSCOH;++j){
    freqInd.data[j]= fBin+2;
  }

  SUB( LALHOUGHConstructHMT( &status, &ht, &freqInd, &phmdVS ), &status );

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
    ERROR( TESTDRIVENDHOUGHC_EFILE, TESTDRIVENDHOUGHC_MSGEFILE, 0 );
    return TESTDRIVENDHOUGHC_EFILE;
  }


  for(k=ySide-1; k>=0; --k){
    for(i=0;i<xSide;++i){
      fprintf( fp ," %f", ht.map[k*xSide +i]);
      fflush( fp );
    }
    fprintf( fp ," \n");
    fflush( fp );
  }

  fclose( fp );


  /******************************************************************/
  /* Free memory and exit */
  /******************************************************************/
  for (j=0;j< MOBSCOH;++j){
    LALFree( pgV.pg[j].peak);  /* All of them */
  }


  for (j=0; j<lutV.length ; ++j){
    for (i=0; i<maxNBorders; ++i){
      LALFree( lutV.lut[j].border[i].xPixel);
    }
    LALFree( lutV.lut[j].border);
    LALFree( lutV.lut[j].bin);
  }

  for(j=0; j<phmdVS.length * phmdVS.nfSize; ++j){
    LALFree( phmdVS.phmd[j].leftBorderP);
    LALFree( phmdVS.phmd[j].rightBorderP);
    LALFree( phmdVS.phmd[j].firstColumn);
  }

  LALFree(lutV.lut);
  LALFree(pgV.pg);
  LALFree(phmdVS.phmd);
  LALFree(freqInd.data);

  LALFree(ht.map);

  LALFree(patch.xCoor);
  LALFree(patch.yCoor);

  LALCheckMemoryLeaks();

  INFO( TESTDRIVENDHOUGHC_MSGENORM );
  return TESTDRIVENDHOUGHC_ENORM;
}

/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */



/** \endcond */
