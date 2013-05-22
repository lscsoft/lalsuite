/*
*  Copyright (C) 2007 Alicia Sintes Olives
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
 * \ingroup pulsarApps
 * \author Badri Krishnan, Alicia Sintes 
 *
 * History:   Created by Sintes June 16, 2003
 *    to test part of the Hough-Driver code.
 *    Case: Non demodulated. Search with no spin-down, source frequency as the 
 *          first search frequency.
 *    No input from SFT data yet implemented here.
 */


#include "./Validation1.h"

/* ************************************************************
 * Usage format string. 
 */

#define USAGE "Usage: %s [-d debuglevel]  [-i ifo (1,2,3)][-o outfile-basename] [-f firstSearchFrequency (in Hz)] [-p alpha delta (in radians)] [-s patchSizeX patchSizeY (in radians)]\n"


/* ***************************************************************
 * Constant Declarations.  Default parameters.
 */

#define EARTHEPHEMERIS "/afs/aeiw/grawave/Linux/lal/lal/packages/pulsar/test/earth03.dat"
#define SUNEPHEMERIS "/afs/aeiw/grawave/Linux/lal/lal/packages/pulsar/test/sun03.dat"

#define IFO 1         /*  detector, 1:GEO, 2:LLO, 3:LHO */

#define F0 500.0          /*  frequency to build the LUT and start search */
/* #define FBAND 1.0           search frequency band  (in Hz) */
#define TCOH 1800.0     /*  time baseline of coherent integration (in seconds) */
#define DF    (1./TCOH)   /*  frequency  resolution. */
#define ALPHA 0.0		/* center of the sky patch (in radians) */
#define DELTA  (-LAL_PI_2)
#define PATCHSIZEX (LAL_PI*0.99) /* patch size */
#define PATCHSIZEY (LAL_PI*0.99)

#define FILEOUT "OutHough.asc"      /* file output */
#define MOBSCOH 60
#define NFSIZE  5

#define T0SEC 730000044
/* #define T0SEC 714153733*/
#define T0NSEC 0
#define JUMPTIME 84600 /* a day  86400-1800 */
/* #define JUMPTIME 604800 a week */

#define ACCURACY 0.01

/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvv------------------------------------ */
int main(int argc, char *argv[]){

  static LALStatus           status;  /* LALStatus pointer */
  
  static LALDetector         detector;
  static LIGOTimeGPSVector   timeV;
  static REAL8Cart3CoorVector velV;
  
  static HOUGHptfLUTVector   lutV; /* the Look Up Table vector*/
  static HOUGHPeakGramVector pgV;
  static PHMDVectorSequence  phmdVS;  /* the partial Hough map derivatives */
  static UINT8FrequencyIndexVector freqInd;

  static HOUGHResolutionPar parRes;
  static HOUGHPatchGrid  patch;   /* Patch description */
  static HOUGHParamPLUT  parLut;  /* parameters needed to build lut  */
  static HOUGHDemodPar   parDem;  /* demodulation parameters or  */
  static HOUGHSizePar    parSize;
  static VelocityPar     velPar;

  static HOUGHMapTotal   ht;   /* the total Hough map */
  /* ------------------------------------------------------- */

  CHAR  *earthEphemeris = NULL;
  CHAR  *sunEphemeris = NULL;
  INT4   ifo;
  REAL8  vel[3];
  INT4   mObsCoh;
  UINT2  maxNBins, maxNBorders;

  INT8   f0Bin;           /* freq. bin to construct LUT */
  INT8   fBin;

  UINT2  xSide, ySide;

  CHAR  *fname = NULL;               /* The output filename */
  FILE  *fp=NULL;                    /* Output file */
  INT4  arg;                         /* Argument counter */
  UINT4 i,j;                       /* Index counter, etc */
  INT4  k;
  REAL8 f0, alpha, delta;
  REAL8 patchSizeX, patchSizeY;
  REAL8 Xx,Xy,Xz;

  /******************************************************************/
  /*    Set up the default parameters.      */
  /* ****************************************************************/

  detector = lalCachedDetectors[LALDetectorIndexGEO600DIFF]; /* default */
  ifo = IFO;
  
  if (ifo ==1) detector=lalCachedDetectors[LALDetectorIndexGEO600DIFF];
  if (ifo ==2) detector=lalCachedDetectors[LALDetectorIndexLLODIFF];
  if (ifo ==3) detector=lalCachedDetectors[LALDetectorIndexLHODIFF];

  earthEphemeris = EARTHEPHEMERIS;
  sunEphemeris = SUNEPHEMERIS;

  mObsCoh = MOBSCOH;
  
  timeV.length   = mObsCoh;
  velV.length    = mObsCoh;
  lutV.length    = mObsCoh;
  pgV.length     = mObsCoh;
  phmdVS.length  = mObsCoh;
  freqInd.length = mObsCoh;
  phmdVS.nfSize  = NFSIZE;

  freqInd.deltaF = DF;
  phmdVS.deltaF  = DF;

  timeV.time = NULL;
  velV.data = NULL;
  lutV.lut = NULL;
  pgV.pg = NULL;
  phmdVS.phmd = NULL;
  freqInd.data = NULL;
  ht.map = NULL;

  f0 =  F0;
  f0Bin = F0*TCOH;

  parRes.f0Bin =  f0Bin;
  parRes.deltaF = DF;
/*
 *   parRes.patchSkySizeX  = patchSizeX = 1.0/(TCOH*F0*VEPI);
 *   parRes.patchSkySizeY  = patchSizeY = 1.0/(TCOH*F0*VEPI);
 */
  parRes.patchSkySizeX  = patchSizeX = PATCHSIZEX;
  parRes.patchSkySizeY  = patchSizeY = PATCHSIZEY;
  parRes.pixelFactor = PIXELFACTOR;
  parRes.pixErr = PIXERR;
  parRes.linErr = LINERR;
  parRes.vTotC = VTOT;

  /* Case: no spins & Non demodulation */
  parDem.deltaF = DF;
  parDem.skyPatch.alpha = ALPHA;
  parDem.skyPatch.delta = DELTA;
  parDem.timeDiff = 0.0;
  parDem.spin.length = 0;
  parDem.spin.data = NULL;
  parDem.positC.x = 0.0;
  parDem.positC.y = 0.0;
  parDem.positC.z = 0.0;

  velPar.detector = detector;
  velPar.tBase = TCOH;
  velPar.vTol = ACCURACY;
  
  alpha = ALPHA;
  delta = DELTA;
  
  /*****************************************************************/
  /*****************************************************************/
  /*    Parse argument list.  i stores the current position.       */
  /*****************************************************************/
  arg = 1;
  while ( arg < argc ) {
    /* Parse debuglevel option. */
    if ( !strcmp( argv[arg], "-d" ) ) {
      if ( argc > arg + 1 ) {
        arg++;
      } else {
        ERROR( VALIDATION1_EARG, VALIDATION1_MSGEARG, 0 );
        XLALPrintError( USAGE, *argv );
        return VALIDATION1_EARG;
      }
    }
    /* Parse interferometer option. */
    else if ( !strcmp( argv[arg], "-i" ) ) {
      if ( argc > arg + 1 ) {
        arg++;
        ifo = atoi( argv[arg++] );
	if (ifo ==1) detector=lalCachedDetectors[LALDetectorIndexGEO600DIFF];
        if (ifo ==2) detector=lalCachedDetectors[LALDetectorIndexLLODIFF];
        if (ifo ==3) detector=lalCachedDetectors[LALDetectorIndexLHODIFF];
        velPar.detector = detector;
      } else {
        ERROR( VALIDATION1_EARG, VALIDATION1_MSGEARG, 0 );
        XLALPrintError( USAGE, *argv );
        return VALIDATION1_EARG;
      }
    }
      /* Parse filename of earth  ephemeris data option. */
      else if ( !strcmp( argv[arg], "-E" ) ) {
        if ( argc > arg + 1 ) {
          arg++;
          earthEphemeris = argv[arg++];
        } else {
          ERROR( VALIDATION1_EARG, VALIDATION1_MSGEARG, 0 );
          XLALPrintError( USAGE, *argv );
          return VALIDATION1_EARG;
        }
      }
      /* Parse filename of sun ephemeris data option. */
      else if ( !strcmp( argv[arg], "-S" ) ) {
        if ( argc > arg + 1 ) {
          arg++;
          sunEphemeris = argv[arg++];
        } else {
          ERROR( VALIDATION1_EARG, VALIDATION1_MSGEARG, 0 );
          XLALPrintError( USAGE, *argv );
          return VALIDATION1_EARG;
        }
      }
    /* Parse output file option. */
    else if ( !strcmp( argv[arg], "-o" ) ) {
      if ( argc > arg + 1 ) {
        arg++;
        fname = argv[arg++];
      } else {
        ERROR( VALIDATION1_EARG, VALIDATION1_MSGEARG, 0 );
        XLALPrintError( USAGE, *argv );
        return VALIDATION1_EARG;
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
        ERROR( VALIDATION1_EARG, VALIDATION1_MSGEARG, 0 );
        XLALPrintError( USAGE, *argv );
        return VALIDATION1_EARG;
      }
    }
    /* Parse sky position options. */
    else if ( !strcmp( argv[arg], "-p" ) ) {
      if ( argc > arg + 2 ) {
        arg++;
        alpha = atof(argv[arg++]);
        delta = atof(argv[arg++]);
	parDem.skyPatch.alpha = alpha;
        parDem.skyPatch.delta = delta;
      } else {
        ERROR( VALIDATION1_EARG, VALIDATION1_MSGEARG, 0 );
        XLALPrintError( USAGE, *argv );
        return VALIDATION1_EARG;
      }
    }
    /* Parse patch size option. */
    else if ( !strcmp( argv[arg], "-s" ) ) {
      if ( argc > arg + 2 ) {
        arg++;
        parRes.patchSkySizeX = patchSizeX = atof(argv[arg++]);
        parRes.patchSkySizeY = patchSizeY = atof(argv[arg++]);
      } else {
        ERROR( VALIDATION1_EARG, VALIDATION1_MSGEARG, 0 );
        XLALPrintError( USAGE, *argv );
        return VALIDATION1_EARG;
      }
    }
    /* Unrecognized option. */
    else {
      ERROR( VALIDATION1_EARG, VALIDATION1_MSGEARG, 0 );
      XLALPrintError( USAGE, *argv );
      return VALIDATION1_EARG;
    }
  } /* End of argument parsing loop. */
  /******************************************************************/

  if ( f0 < 0 ) {
    ERROR( VALIDATION1_EBAD, VALIDATION1_MSGEBAD, "freq<0:" );
    XLALPrintError( USAGE, *argv  );
    return VALIDATION1_EBAD;
  }

  /******************************************************************/
  /******************************************************************/
  /* create time stamps (for a test) */
  /******************************************************************/
  timeV.time = (LIGOTimeGPS *)LALMalloc(mObsCoh*sizeof(LIGOTimeGPS));
  timeV.time[0].gpsSeconds = T0SEC;
  timeV.time[0].gpsNanoSeconds = T0NSEC;
  for(j=1; j<timeV.length; ++j){
    timeV.time[j].gpsSeconds = timeV.time[j-1].gpsSeconds + TCOH + JUMPTIME;
    timeV.time[j].gpsNanoSeconds = T0NSEC;
  }

  /******************************************************************/
  /* compute detector velocity for those time stamps (for a test) */
  /******************************************************************/

  velV.data = (REAL8Cart3Coor *)LALMalloc(mObsCoh*sizeof(REAL8Cart3Coor));
  velPar.edat = NULL; 
  {
    EphemerisData    *edat=NULL;
   
    /*  ephemeris info */
    edat = (EphemerisData *)LALMalloc(sizeof(EphemerisData));
   (*edat).ephiles.earthEphemeris = earthEphemeris;
   (*edat).ephiles.sunEphemeris = sunEphemeris;

    /* read in ephemeris data */
    SUB( LALInitBarycenter( &status, edat), &status);
    velPar.edat = edat;
    
    for(j=0; j<velV.length; ++j){
      velPar.startTime.gpsSeconds     = timeV.time[j].gpsSeconds;
      velPar.startTime.gpsNanoSeconds = timeV.time[j].gpsNanoSeconds;
    
      SUB( LALAvgDetectorVel ( &status, vel, &velPar), &status );
      velV.data[j].x= vel[0];
      velV.data[j].y= vel[1];
      velV.data[j].z= vel[2];   
    }
    LALFree(edat->ephemE);
    LALFree(edat->ephemS);
    LALFree(edat);
  }

  /******************************************************************/  
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
  /******************************************************************/
  /* memory allocation and settings */
  /******************************************************************/

  lutV.lut = (HOUGHptfLUT *)LALMalloc(mObsCoh*sizeof(HOUGHptfLUT));
  pgV.pg = (HOUGHPeakGram *)LALMalloc(mObsCoh*sizeof(HOUGHPeakGram));
  phmdVS.phmd =(HOUGHphmd *)LALMalloc(mObsCoh*NFSIZE*sizeof(HOUGHphmd));
  freqInd.data =  ( UINT8 *)LALMalloc(mObsCoh*sizeof(UINT8));

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
  /******************************************************************/
  /* create Peakgrams for testing                                   */
  /******************************************************************/

  /* let us fix a source at a given position 
     (not at the center of the patch, but a pixel center) */
  { 
    UINT2              xPos, yPos;
    REAL8Cart2Coor     sourceProjected;
    REAL8UnitPolarCoor sourceRotated;
    REAL8UnitPolarCoor skyPatchCenter;
    REAL8UnitPolarCoor sourceLocation;
   
    xPos = xSide/2;
    yPos = ySide/2;

     sourceProjected.x = patch.xCoor[xPos];
    sourceProjected.y = patch.yCoor[yPos];
    skyPatchCenter.alpha = parDem.skyPatch.alpha;
    skyPatchCenter.delta = parDem.skyPatch.delta;
    
    /* invert the stereographic projection for a point on the projected plane */
    SUB( LALStereoInvProjectCart( &status, &sourceRotated, &sourceProjected ), &status );
 
    /* undo roation in case the patch is not centered at the south pole */
    SUB( LALInvRotatePolarU( &status, &sourceLocation, &sourceRotated, &skyPatchCenter ), &status );

    Xx= cos(sourceLocation.delta)* cos(sourceLocation.alpha);
    Xy= cos(sourceLocation.delta)* sin(sourceLocation.alpha);
    Xz= sin(sourceLocation.delta); 
  }
  
  for (j=0;j< mObsCoh; ++j) {  /* create all the peakgrams */
    REAL8 veldotX;
    
    pgV.pg[j].deltaF = DF;
/*
 *     pgV.pg[j].fBinIni = f0Bin/2; 
 *     pgV.pg[j].fBinFin = f0Bin*2;
 */
    pgV.pg[j].fBinIni = f0Bin-maxNBins; 
    pgV.pg[j].fBinFin = f0Bin+3*maxNBins;
    /* pgV.pg[j].length = 2; */
    pgV.pg[j].length = 1; 
    pgV.pg[j].peak = NULL;
    pgV.pg[j].peak = (INT4 *)LALMalloc( ( pgV.pg[j].length) * sizeof(INT4));
    
    veldotX = Xx*velV.data[j].x + Xy*velV.data[j].y + Xz*velV.data[j].z;
    /*  fBin = [ f0Bin * ( 1 + veldotX ) ] */
    pgV.pg[j].peak[0]= floor( f0Bin*(1.0+veldotX) +0.5) - pgV.pg[j].fBinIni;
     /* pgV.pg[j].peak[1]= pgV.pg[j].peak[0]+2; */
   
  }
 

 
  /******************************************************************/
  /******************************************************************/
  /* create            all the LUTs                           */
  /******************************************************************/
         
  for (j=0;j< mObsCoh;++j){  /* create all the LUTs */
    parDem.veloC.x = velV.data[j].x;
    parDem.veloC.y = velV.data[j].y;
    parDem.veloC.z = velV.data[j].z;
    
    /* calculate parameters needed for buiding the LUT */
    SUB( LALNDHOUGHParamPLUT( &status, &parLut, &parSize, &parDem ), &status );

    /* build the LUT */
    SUB( LALHOUGHConstructPLUT( &status, &(lutV.lut[j]), &patch, &parLut ),
         &status );
  }

  /******************************************************************/
  /* starting the search  (A very simple case for only 1 frequency) */
  /******************************************************************/
  /* build the set of  PHMD  starting at f0Bin*/
  /******************************************************************/
  fBin = f0Bin;
  phmdVS.fBinMin = fBin;
  SUB( LALHOUGHConstructSpacePHMD(&status, &phmdVS, &pgV, &lutV), &status );

  /* shift the structure one frequency bin */
  /*  SUB( LALHOUGHupdateSpacePHMDup(&status, &phmdVS, &pgV, &lutV), &status );*/


  /******************************************************************/
  /* initializing the Hough map space */
  /******************************************************************/

  SUB( LALHOUGHInitializeHT( &status, &ht, &patch ), &status );

  /******************************************************************/
  /* construction of a total Hough map  */
  /******************************************************************/

  for (j=0;j< mObsCoh;++j){
    freqInd.data[j]= fBin;
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
    ERROR( VALIDATION1_EFILE, VALIDATION1_MSGEFILE, 0 );
    return VALIDATION1_EFILE;
  }


  for(k=ySide-1; k>=0; --k){
    for(i=0;i<xSide;++i){
      fprintf( fp ," %d", ht.map[k*xSide +i]);
      fflush( fp );
    }
    fprintf( fp ," \n");
    fflush( fp );
  }

  fclose( fp );


  /******************************************************************/
  /* Free memory and exit */
  /******************************************************************/
  for (j=0;j< mObsCoh;++j){
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

  LALFree(timeV.time);
  LALFree(velV.data);

  LALFree(lutV.lut);
  LALFree(pgV.pg);
  LALFree(phmdVS.phmd);
  LALFree(freqInd.data);

  LALFree(ht.map);

  LALFree(patch.xCoor);
  LALFree(patch.yCoor);

  LALCheckMemoryLeaks();

  INFO( VALIDATION1_MSGENORM );
  return VALIDATION1_ENORM;
}

/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */




