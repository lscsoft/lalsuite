/*
*  Copyright (C) 2007 Chris Messenger
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

#ifndef _SIDEBAND_H
#define _SIDEBAND_H

#include <lal/LALDatatypes.h>
#include <lal/DetectorSite.h>
#include <lal/LALBarycenter.h>
#include <lal/LALDemod.h>
#include <lal/Date.h>

/* C++ protection. */
#ifdef  __cplusplus
extern "C" {
#endif
 
#define NPARAMS 9
#define NPARAMS_IL 3

/**
 * A structure designed to store all the information required to describe a binary orbit
 */
  typedef struct {
    REAL8 OrbitalSemiMajorAxis; /**< */
    REAL8 OrbitalPeriod;
    REAL8 OrbitalEccentricity;
    REAL8 ArgumentofPeriapse;
    LIGOTimeGPS TimeofSSBPeriapsePassage;
    REAL8 alpha;
    REAL8 delta;
    REAL8 f0;
    REAL8 phi0;
    REAL8 psi;
    REAL8 cosi;
    REAL8 h0;
  } BinarySourceParams;
  
  /** A structure for storing ABCcoeffients */
  typedef struct {
    REAL8 omega0;
    REAL8 aco[3];
    REAL8 bco[3];
    REAL8 uco[5];
    REAL8 vco[5];
    REAL8 wco[5];
    REAL8 apco[3];
    REAL8 bpco[3];
    REAL8 upco[5];
    REAL8 vpco[5];
    REAL8 wpco[5];
  } ABCcoefficients;
  
  /** A structure designed to store general input template parameters */
  typedef struct {
    INT4 Mmax;
    INT4 Mmin;
    LIGOTimeGPS tstartSSB;
    LIGOTimeGPS tstart;
    LIGOTimeGPS reftime;
    ABCcoefficients *ABCco;
    REAL8 T;
    REAL8 Tobs;
    REAL8 sqrtSh;
    REAL8 alpha;
    REAL8 delta;
    REAL8 w0;
    REAL8Vector *freqsamples;
    BOOLEAN local;
    INT4 pmax;
    COMPLEX16Vector *wa;
    COMPLEX16Vector *wb;
    REAL8 dfwindow;
    INT4 windowrange;
    INT4Vector *timestamps;
    INT4Vector *gapvectorstart;
    INT4Vector *gapvectorend;
    INT4 mmin;
    INT4 mmax;
    INT4 nsft;
  } SideBandTemplateParams;
  
  /** A structure for storing a frequency domain signal template */
  typedef struct {
    REAL8 minfreq;
    REAL8 deltaf;
    UINT4 length;
    LIGOTimeGPS epoch;
    COMPLEX16Vector *fourier;
  } SideBandTemplate;

 
  /** A structure for storing ABCcoeffient parameters */
  typedef struct {
    LIGOTimeGPS tstart;
    REAL8 alpha;
    REAL8 delta;
  } ABCcoParams;

  /** A structure for storing an MCMC signal parameter vector */
  typedef struct {
    REAL8 f0;
    REAL8 period;
    REAL8 a;
    LIGOTimeGPS tp;
    REAL8 argp;
    REAL8 e;
    REAL8 psi;
    REAL8 cosi;
    REAL8 phi0;
    REAL8 h0;
    REAL8 logL;
    REAL8 logWt;
    BOOLEAN accepted;
  } SideBandMCMCVector;

  /** A structure for storing a fourier dataset  */
  typedef struct {
    COMPLEX16Vector *fourier;
    REAL8Vector *freq;
  } SideBandDataSet;

  /** A structure for storing the MCMCM parameter ranges  */
  typedef struct {
    REAL8 f0min;
    REAL8 f0max;
    REAL8 periodmin;
    REAL8 periodmax;
    REAL8 amin;
    REAL8 amax;
    LIGOTimeGPS tpmin;
    LIGOTimeGPS tpmax;
    REAL8 argpmin;
    REAL8 argpmax;
    REAL8 emin;
    REAL8 emax;
    REAL8 h0min;
    REAL8 h0max;
    REAL8 psimin;
    REAL8 psimax;
    REAL8 cosimin;
    REAL8 cosimax;
    REAL8 phi0min;
    REAL8 phi0max;
  } SideBandMCMCRanges;

  /** A structure for storing the MCMCM parameter jump sizes  */
  typedef struct {
    REAL8 jump[3];
    REAL8 prob[3];
  } SideBandMCMCJumpProb;

  /** A structure for storing the MCMCM parameter jump sizes  */
  typedef struct {
    SideBandMCMCJumpProb f0;
    SideBandMCMCJumpProb period;
    SideBandMCMCJumpProb a;
    SideBandMCMCJumpProb x;
    SideBandMCMCJumpProb y;
    SideBandMCMCJumpProb e;
    SideBandMCMCJumpProb h0;
    SideBandMCMCJumpProb psi;
    SideBandMCMCJumpProb cosi;
    SideBandMCMCJumpProb phi0;
  } SideBandMCMCJumpProbs;
  


  /* A structure to store parameters required for selecting frequencies for use in the MCMC */
  typedef struct {
    SideBandMCMCRanges ranges;
    INT4 mmax;
    INT4 mmin;
    REAL8Vector *minf;
    REAL8Vector *maxf;
    BOOLEAN am;
    REAL8 df;
  } SelectSideBandFrequencyParams;

  /* A structure to store parameters required for reading fourier data */
  typedef struct {
    CHAR file[256];
    REAL8 minf;
    REAL8 maxf;
    REAL8 Tobs;
    REAL8 df;
    INT4 nsft;
  } ReadSideBandDataParams;

  /* A structure to store parameters required for reading fourier data */
  typedef struct {
    REAL8 minfreq;
    REAL8 maxfreq;
    REAL8Vector *minf;
    REAL8Vector *maxf;
    INT4 Nthresh;
    REAL8 sqrtSh;
    REAL8 safedf;
  } EstimateSideBandNoiseParams;

  /* LALDemod functions now put into CFSLALDemod.c */
  void GenerateSideBandTemplate(LALStatus *,BinarySourceParams *,SideBandTemplateParams *,SideBandTemplate **);
  void ComputeABCcoefficients (LALStatus *,ABCcoParams *,LALDetector *,ABCcoefficients **);
  void ReadTimeStamps(LALStatus *,CHAR *,INT4, SideBandTemplateParams **);
  void ComputeSideBandWindow(LALStatus *,ABCcoefficients *, CHAR *,SideBandTemplateParams **);
  void ComputeSideBandLikelihood(LALStatus *,SideBandMCMCVector *,SideBandDataSet *,SideBandTemplate **,SideBandTemplateParams *);
  void InitEphemeris(LALStatus *,EphemerisData *,const CHAR *,const CHAR *);
  void ReadSideBandPriors(LALStatus *status,CHAR *,SideBandMCMCRanges *,SideBandMCMCJumpProbs *);
  void SelectSideBandFrequencies (LALStatus *,SideBandDataSet **,SelectSideBandFrequencyParams *,SideBandDataSet **);
  void ReadSideBandData (LALStatus *,ReadSideBandDataParams *,SideBandDataSet **);
  void EstimateSideBandNoise (LALStatus *, SideBandDataSet **,EstimateSideBandNoiseParams *);
  REAL4 bessj(INT4, REAL4);
  REAL4 bessj1(REAL4);
  REAL4 bessj0(REAL4);
      

#ifdef  __cplusplus
}
#endif  
/* C++ protection. */

#endif  /* Double-include protection. */
