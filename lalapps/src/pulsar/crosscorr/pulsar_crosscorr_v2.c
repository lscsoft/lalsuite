/*
 *  Copyright (C) 2013 Badri Krishnan, Shane Larson, John Whelan
 *  Copyright (C) 2013, 2014 Badri Krishnan, John Whelan, Yuanhao Zhang
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
 * \author B.Krishnan, S.Larson, J.T.Whelan, Y.Zhang
 * \date 2013, 2014
 * \file pulsar_crosscorr_v2.c
 * \ingroup pulsarApps
 * \brief Perform CW cross-correlation search - version 2
 *
 */


/*lalapps includes */
#include <lalapps.h>
#include <lal/UserInput.h>
#include <lal/SFTfileIO.h>
#include <lal/SFTutils.h>
#include <lal/LogPrintf.h>
#include <lal/DopplerScan.h>
#include <lal/DetectorStates.h>
#include <lal/ExtrapolatePulsarSpins.h>
#include <lal/LALInitBarycenter.h>
#include <lal/NormalizeSFTRngMed.h>
#include <lal/LALString.h>
#include <lal/PulsarCrossCorr_v2.h>
#include "CrossCorrToplist.h"

/* introduce mismatch in f and all 5 binary parameters */

/* user input variables */
typedef struct{
  BOOLEAN help; /**< if the user wants a help message */
  INT4    startTime;          /**< desired start GPS time of search */
  INT4    endTime;            /**< desired end GPS time */
  REAL8   maxLag;             /**< maximum lag time in seconds between SFTs in correlation */
  BOOLEAN inclAutoCorr;       /**< include auto-correlation terms (an SFT with itself) */
  REAL8   fStart;             /**< start frequency in Hz */
  REAL8   fBand;              /**< frequency band to search over in Hz */
  /* REAL8   fdotStart;          /\**< starting value for first spindown in Hz/s*\/ */
  /* REAL8   fdotBand;           /\**< range of first spindown to search over in Hz/s *\/ */
  REAL8   alphaRad;           /**< right ascension in radians */
  REAL8   deltaRad;           /**< declination in radians */
  REAL8   refTime;            /**< reference time for pulsar phase definition */
  REAL8   orbitAsiniSec;      /**< start projected semimajor axis in seconds */
  REAL8   orbitAsiniSecBand;  /**< band for projected semimajor axis in seconds */
  REAL8   orbitPSec;          /**< binary orbital period in seconds */
  REAL8   orbitTimeAsc;       /**< start time of ascension for binary orbit */
  REAL8   orbitTimeAscBand;   /**< band for time of ascension for binary orbit */
  CHAR    *sftLocation;       /**< location of SFT data */
  CHAR    *ephemEarth;	      /**< Earth ephemeris file to use */
  CHAR    *ephemSun;	      /**< Sun ephemeris file to use */
  INT4    rngMedBlock;        /**< running median block size */
  INT4    numBins;            /**< number of frequency bins to include in sum */
  REAL8   mismatchF;          /**< mismatch for frequency spacing */
  REAL8   mismatchA;          /**< mismatch for spacing in semi-major axis */
  REAL8   mismatchT;          /**< mismatch for spacing in time of periapse passage */
  REAL8   mismatchP;          /**< mismatch for spacing in period */
  INT4    numCand;            /**< number of candidates to keep in output toplist */
  BOOLEAN readlist;           /**Flag indicate whether read in SFT-pair list, if FALSE, output SFT-pair list*/
  CHAR    *toplistFilename;   /**< output filename containing candidates in toplist */
} UserInput_t;

/* struct to store useful variables */
typedef struct{
  SFTCatalog *catalog; /**< catalog of SFTs */
  EphemerisData *edat; /**< ephemeris data */
} ConfigVariables;

#define TRUE (1==1)
#define FALSE (1==0)
#define MAXFILENAMELENGTH 512

/* local function prototypes */
int XLALInitUserVars ( UserInput_t *uvar );
int XLALInitializeConfigVars (ConfigVariables *config, const UserInput_t *uvar);
int XLALDestroyConfigVars (ConfigVariables *config);
int GetNextCrossCorrTemplate(BOOLEAN *binaryParamsFlag, PulsarDopplerParams *dopplerpos, PulsarDopplerParams *binaryTemplateSpacings, PulsarDopplerParams *minBinaryTemplate, PulsarDopplerParams *maxBinaryTemplate);


int main(int argc, char *argv[]){

  UserInput_t XLAL_INIT_DECL(uvar);
  static ConfigVariables config;

  /* sft related variables */
  MultiSFTVector *inputSFTs = NULL;
  MultiPSDVector *multiPSDs = NULL;
  MultiNoiseWeights *multiWeights = NULL;
  MultiLIGOTimeGPSVector *multiTimes = NULL;
  MultiLALDetector multiDetectors;
  MultiDetectorStateSeries *multiStates = NULL;
  MultiAMCoeffs *multiCoeffs = NULL;
  SFTIndexList *sftIndices = NULL;
  SFTPairIndexList *sftPairs = NULL;
  REAL8Vector *shiftedFreqs = NULL;
  UINT4Vector *lowestBins = NULL;
  REAL8Vector *kappaValues = NULL;
  REAL8Vector *signalPhases = NULL;

  PulsarDopplerParams XLAL_INIT_DECL(dopplerpos);
  PulsarDopplerParams thisBinaryTemplate, binaryTemplateSpacings;
  PulsarDopplerParams minBinaryTemplate, maxBinaryTemplate;
  SkyPosition XLAL_INIT_DECL(skyPos);
  MultiSSBtimes *multiSSBTimes = NULL;
  MultiSSBtimes *multiBinaryTimes = NULL;

  INT4  k;
  UINT4 j;
  REAL8 fMin, fMax; /* min and max frequencies read from SFTs */
  REAL8 deltaF; /* frequency resolution associated with time baseline of SFTs */

  REAL8Vector *curlyGUnshifted = NULL;
  REAL8 diagff = 0; /*diagonal metric components*/
  REAL8 diagaa = 0;
  REAL8 diagTT = 0;
  REAL8 diagpp = 1;
  REAL8 ccStat = 0;
  REAL8 evSquared=0;
  REAL8 estSens=0; /*estimated sensitivity(4.13)*/
  BOOLEAN dopplerShiftFlag = FALSE;
  toplist_t *ccToplist=NULL;
  CrossCorrBinaryOutputEntry thisCandidate;
  UINT4 checksum;

  LogPrintf (LOG_CRITICAL, "starting time\n", 0);/*for debug convenience to record calculating time*/
  /* initialize and register user variables */
  if ( XLALInitUserVars( &uvar ) != XLAL_SUCCESS ) {
    LogPrintf ( LOG_CRITICAL, "%s: XLALInitUserVars() failed with errno=%d\n", __func__, xlalErrno );
    XLAL_ERROR( XLAL_EFUNC );
  }

  /* read user input from the command line or config file */
  if ( XLALUserVarReadAllInput ( argc, argv ) != XLAL_SUCCESS ) {
    LogPrintf ( LOG_CRITICAL, "%s: XLALUserVarReadAllInput() failed with errno=%d\n", __func__, xlalErrno );
    XLAL_ERROR( XLAL_EFUNC );
  }

  if (uvar.help)	/* if help was requested, then exit */
    return 0;

  /* configure useful variables based on user input */
  if ( XLALInitializeConfigVars ( &config, &uvar) != XLAL_SUCCESS ) {
    LogPrintf ( LOG_CRITICAL, "%s: XLALInitUserVars() failed with errno=%d\n", __func__, xlalErrno );
    XLAL_ERROR( XLAL_EFUNC );
  }

  deltaF = config.catalog->data[0].header.deltaF;
  REAL8 Tsft = 1.0 / deltaF;

  /* create the toplist */
  create_crossCorrBinary_toplist( &ccToplist, uvar.numCand);
  /* now read the data */

  /* /\* get SFT parameters so that we can initialise search frequency resolutions *\/ */
  /* /\* calculate deltaF_SFT *\/ */
  /* deltaF_SFT = catalog->data[0].header.deltaF;  /\* frequency resolution *\/ */
  /* timeBase= 1.0/deltaF_SFT; /\* sft baseline *\/ */

  /* /\* catalog is ordered in time so we can get start, end time and tObs *\/ */
  /* firstTimeStamp = catalog->data[0].header.epoch; */
  /* lastTimeStamp = catalog->data[catalog->length - 1].header.epoch; */
  /* tObs = XLALGPSDiff( &lastTimeStamp, &firstTimeStamp ) + timeBase; */

  /* /\*set pulsar reference time *\/ */
  /* if (LALUserVarWasSet ( &uvar_refTime )) { */
  /*   XLALGPSSetREAL8(&refTime, uvar_refTime); */
  /* }  */
  /* else {	/\*if refTime is not set, set it to midpoint of sfts*\/ */
  /*   XLALGPSSetREAL8(&refTime, (0.5*tObs) + XLALGPSGetREAL8(&firstTimeStamp));  */
  /* } */

  /* /\* set frequency resolution defaults if not set by user *\/ */
  /* if (!(LALUserVarWasSet (&uvar_fResolution))) { */
  /*   uvar_fResolution = 1/tObs; */
  /* } */

  /* { */
  /*   /\* block for calculating frequency range to read from SFTs *\/ */
  /*   /\* user specifies freq and fdot range at reftime */
  /*      we translate this range of fdots to start and endtime and find */
  /*      the largest frequency band required to cover the  */
  /*      frequency evolution  *\/ */
  /*   PulsarSpinRange spinRange_startTime; /\**< freq and fdot range at start-time of observation *\/ */
  /*   PulsarSpinRange spinRange_endTime;   /\**< freq and fdot range at end-time of observation *\/ */
  /*   PulsarSpinRange spinRange_refTime;   /\**< freq and fdot range at the reference time *\/ */

  /*   REAL8 startTime_freqLo, startTime_freqHi, endTime_freqLo, endTime_freqHi, freqLo, freqHi; */

  /*   REAL8Vector *fdotsMin=NULL; */
  /*   REAL8Vector *fdotsMax=NULL; */

  /*   UINT4 k; */

  /*   fdotsMin = (REAL8Vector *)LALCalloc(1, sizeof(REAL8Vector)); */
  /*   fdotsMin->length = N_SPINDOWN_DERIVS; */
  /*   fdotsMin->data = (REAL8 *)LALCalloc(fdotsMin->length, sizeof(REAL8)); */

  /*   fdotsMax = (REAL8Vector *)LALCalloc(1, sizeof(REAL8Vector)); */
  /*   fdotsMax->length = N_SPINDOWN_DERIVS; */
  /*   fdotsMax->data = (REAL8 *)LALCalloc(fdotsMax->length, sizeof(REAL8)); */

  /*   XLAL_INIT_MEM(spinRange_startTime); */
  /*   XLAL_INIT_MEM(spinRange_endTime); */
  /*   XLAL_INIT_MEM(spinRange_refTime); */

  /*   spinRange_refTime.refTime = refTime; */
  /*   spinRange_refTime.fkdot[0] = uvar_f0; */
  /*   spinRange_refTime.fkdotBand[0] = uvar_fBand; */
  /* } */

  /* FIXME: need to correct fMin and fMax for Doppler shift, rngmedian bins and spindown range */
  /* this is essentially just a place holder for now */
  /* FIXME: this running median buffer is overkill, since the running median block need not be centered on the search frequency */
  REAL8 vMax = LAL_TWOPI * uvar.orbitAsiniSec /uvar.orbitPSec + LAL_TWOPI * LAL_REARTH_SI / (LAL_DAYSID_SI * LAL_C_SI) + LAL_TWOPI * LAL_AU_SI/(LAL_YRSID_SI * LAL_C_SI); /*calculate the maximum relative velocity in speed of light*/
  fMin = uvar.fStart * (1 - vMax) - 0.5 * uvar.rngMedBlock * deltaF;
  fMax = (uvar.fStart + uvar.fBand) * (1 + vMax) + 0.5 * uvar.rngMedBlock * deltaF;

  /* read the SFTs*/
  if ((inputSFTs = XLALLoadMultiSFTs ( config.catalog, fMin, fMax)) == NULL){
    LogPrintf ( LOG_CRITICAL, "%s: XLALLoadMultiSFTs() failed with errno=%d\n", __func__, xlalErrno );
    XLAL_ERROR( XLAL_EFUNC );
  }

  /* calculate the psd and normalize the SFTs */
  if (( multiPSDs =  XLALNormalizeMultiSFTVect ( inputSFTs, uvar.rngMedBlock, NULL )) == NULL){
    LogPrintf ( LOG_CRITICAL, "%s: XLALNormalizeMultiSFTVect() failed with errno=%d\n", __func__, xlalErrno );
    XLAL_ERROR( XLAL_EFUNC );
  }

  /* compute the noise weights for the AM coefficients */
  if (( multiWeights = XLALComputeMultiNoiseWeights ( multiPSDs, uvar.rngMedBlock, 0 )) == NULL){
    LogPrintf ( LOG_CRITICAL, "%s: XLALComputeMultiNoiseWeights() failed with errno=%d\n", __func__, xlalErrno );
    XLAL_ERROR( XLAL_EFUNC );
  }

  /* read the timestamps from the SFTs */
  if ((multiTimes = XLALExtractMultiTimestampsFromSFTs ( inputSFTs )) == NULL){
    LogPrintf ( LOG_CRITICAL, "%s: XLALExtractMultiTimestampsFromSFTs() failed with errno=%d\n", __func__, xlalErrno );
    XLAL_ERROR( XLAL_EFUNC );
  }

  /* read the detector information from the SFTs */
  if ( XLALMultiLALDetectorFromMultiSFTs ( &multiDetectors, inputSFTs ) != XLAL_SUCCESS){
    LogPrintf ( LOG_CRITICAL, "%s: XLALMultiLALDetectorFromMultiSFTs() failed with errno=%d\n", __func__, xlalErrno );
    XLAL_ERROR( XLAL_EFUNC );
  }

  /* Find the detector state for each SFT */
  if ((multiStates = XLALGetMultiDetectorStates ( multiTimes, &multiDetectors, config.edat, 0.0 )) == NULL){
    LogPrintf ( LOG_CRITICAL, "%s: XLALGetMultiDetectorStates() failed with errno=%d\n", __func__, xlalErrno );
    XLAL_ERROR( XLAL_EFUNC );
  }

  /* Note this is specialized to a single sky position */
  /* This might need to be moved into the config variables */
  skyPos.system = COORDINATESYSTEM_EQUATORIAL;
  skyPos.longitude = uvar.alphaRad;
  skyPos.latitude  = uvar.deltaRad;

  /* Calculate the AM coefficients (a,b) for each SFT */
  if ((multiCoeffs = XLALComputeMultiAMCoeffs ( multiStates, multiWeights, skyPos )) == NULL){
    LogPrintf ( LOG_CRITICAL, "%s: XLALComputeMultiAMCoeffs() failed with errno=%d\n", __func__, xlalErrno );
    XLAL_ERROR( XLAL_EFUNC );
  }

  /* Construct the flat list of SFTs (this sort of replicates the
     catalog, but there's not an obvious way to get the information
     back) */

  if ( ( XLALCreateSFTIndexListFromMultiSFTVect( &sftIndices, inputSFTs ) != XLAL_SUCCESS ) ) {
    LogPrintf ( LOG_CRITICAL, "%s: XLALCreateSFTIndexListFromMultiSFTVect() failed with errno=%d\n", __func__, xlalErrno );
    XLAL_ERROR( XLAL_EFUNC );
  }

  /* Construct the list of SFT pairs */
#define SFTPAIRFILENAME "SFTpairlist.dat"
#define SFTFILENAME "SFTlist.dat"
#define PCC_SFTPAIR_HEADER "#The length of SFT-pair list is %u\n"
#define PCC_SFT_HEADER "#The length of SFT list is %u\n"
  FILE *fp = NULL;

  if (uvar.readlist == FALSE) /*Probably no need to introduce a flag, while it can just test whether the file exists or not, but it works well now/Harry*/
    {
      if ( ( XLALCreateSFTPairIndexList( &sftPairs, sftIndices, inputSFTs, uvar.maxLag, uvar.inclAutoCorr ) != XLAL_SUCCESS ) ) {
	LogPrintf ( LOG_CRITICAL, "%s: XLALCreateSFTPairIndexList() failed with errno=%d\n", __func__, xlalErrno );
	XLAL_ERROR( XLAL_EFUNC );
      }

      fp = fopen(SFTPAIRFILENAME,"w");
      fprintf(fp,PCC_SFTPAIR_HEADER, sftPairs->length ); /*output the length of SFT-pair list to the header*/
      for(j = 0; j < sftPairs->length; j++){
	fprintf(fp, "%u %u\n", sftPairs->data[j].sftNum[0], sftPairs->data[j].sftNum[1]);
      }
      fclose(fp);

      fp = fopen(SFTFILENAME,"w");
      fprintf(fp,PCC_SFTPAIR_HEADER, sftPairs->length ); /*output the length of SFT list to the header*/
      for(j = 0; j < sftIndices->length; j++){ /*output the SFT list */
	fprintf(fp,"%u %u\n", sftIndices->data[j].detInd, sftIndices->data[j].sftInd);
      }
      fclose(fp);
    }
  else
    {
      if( (fp = fopen(SFTPAIRFILENAME,"r")) == NULL){
	LogPrintf ( LOG_CRITICAL, "didn't find SFT-pari list while readin = TRUE\n", 0);
	XLAL_ERROR( XLAL_EFUNC );
      }
      fscanf(fp,PCC_SFTPAIR_HEADER, &sftPairs->length);
      for(j = 0; j < sftPairs->length; j++){ /*read in  the SFT-pair list */
	fscanf(fp, "%u %u\n", &sftPairs->data[j].sftNum[0], &sftPairs->data[j].sftNum[1]);
      }
      fclose(fp);
    }

  /* Get weighting factors for calculation of metric */
  /* note that the sigma-squared is now absorbed into the curly G
     because the AM coefficients are noise-weighted. */
  if ( ( XLALCalculateAveCurlyGAmpUnshifted( &curlyGUnshifted, sftPairs, sftIndices, multiCoeffs)  != XLAL_SUCCESS ) ) {
    LogPrintf ( LOG_CRITICAL, "%s: XLALCalculateAveCurlyGUnshifted() failed with errno=%d\n", __func__, xlalErrno );
    XLAL_ERROR( XLAL_EFUNC );
  }
  /*initialize binary parameters structure*/
  XLAL_INIT_MEM(minBinaryTemplate);
  XLAL_INIT_MEM(maxBinaryTemplate);
  XLAL_INIT_MEM(thisBinaryTemplate);
  XLAL_INIT_MEM(binaryTemplateSpacings);
  /*fill in minbinaryOrbitParams*/
  XLALGPSSetREAL8( &minBinaryTemplate.tp, uvar.orbitTimeAsc);
  minBinaryTemplate.argp = 0.0;
  minBinaryTemplate.asini = uvar.orbitAsiniSec;
  minBinaryTemplate.ecc = 0.0;
  minBinaryTemplate.period = uvar.orbitPSec;
  minBinaryTemplate.fkdot[0] = uvar.fStart;
  /*fill in maxBinaryParams*/
  XLALGPSSetREAL8( &maxBinaryTemplate.tp, uvar.orbitTimeAsc + uvar.orbitTimeAscBand);
  maxBinaryTemplate.argp = 0.0;
  maxBinaryTemplate.asini = uvar.orbitAsiniSec + uvar.orbitAsiniSecBand;
  maxBinaryTemplate.ecc = 0.0;
  maxBinaryTemplate.period = uvar.orbitPSec;
  maxBinaryTemplate.fkdot[0] = uvar.fStart + uvar.fBand;
  /*fill in thisBinaryTemplate*/
  XLALGPSSetREAL8( &thisBinaryTemplate.tp, uvar.orbitTimeAsc);
  thisBinaryTemplate.argp = 0.0;
  thisBinaryTemplate.asini = 0.5*(minBinaryTemplate.asini + maxBinaryTemplate.asini);
  thisBinaryTemplate.ecc = 0.0;
  thisBinaryTemplate.period =0.5*(minBinaryTemplate.period + maxBinaryTemplate.period);
  thisBinaryTemplate.fkdot[0]=0.5*(minBinaryTemplate.fkdot[0] + maxBinaryTemplate.fkdot[0]);

  /*Get metric diagonal components, also estimate sensitivity i.e. E[rho]/(h0)^2 (4.13)*/
  if ( (XLALFindLMXBCrossCorrDiagMetric(&estSens, &diagff, &diagaa, &diagTT, thisBinaryTemplate, curlyGUnshifted, sftPairs, sftIndices, inputSFTs/*, kappaValues*/)  != XLAL_SUCCESS ) ) {
    LogPrintf ( LOG_CRITICAL, "%s: XLALFindLMXBCrossCorrDiagMetric() failed with errno=%d\n", __func__, xlalErrno );
    XLAL_ERROR( XLAL_EFUNC );
  }

  /* initialize the doppler scan struct which stores the current template information */
  XLALGPSSetREAL8(&dopplerpos.refTime, uvar.refTime);
  dopplerpos.Alpha = uvar.alphaRad;
  dopplerpos.Delta = uvar.deltaRad;
  dopplerpos.fkdot[0] = uvar.fStart - sqrt(uvar.mismatchF / diagff);/*minus a  frequency spacing for testing convenience to search at true parameter, may need to modify*/
  /* set all spindowns to zero */
  for (k=1; k < PULSAR_MAX_SPINS; k++)
    dopplerpos.fkdot[k] = 0.0;

  /* now set the initial values of binary parameters */
  thisBinaryTemplate.asini = uvar.orbitAsiniSec;
  thisBinaryTemplate.period = uvar.orbitPSec;
  XLALGPSSetREAL8( &thisBinaryTemplate.tp, uvar.orbitTimeAsc);
  thisBinaryTemplate.ecc = 0.0;
  thisBinaryTemplate.argp = 0.0;
  /* copy to dopplerpos */
  dopplerpos.asini = thisBinaryTemplate.asini;
  dopplerpos.period = thisBinaryTemplate.period;
  dopplerpos.tp = thisBinaryTemplate.tp;
  dopplerpos.ecc = thisBinaryTemplate.ecc;
  dopplerpos.argp = thisBinaryTemplate.argp;

  /* spacing in frequency from diagff */ /* set spacings in new dopplerparams struct */
  binaryTemplateSpacings.fkdot[0] = sqrt(uvar.mismatchF / diagff);
  binaryTemplateSpacings.asini = sqrt(uvar.mismatchA / diagaa);
  binaryTemplateSpacings.period = sqrt(uvar.mismatchP / diagpp);
  /* this is annoying: tp is a GPS time while we want a difference
     in time which should be just REAL8 */
  XLALGPSSetREAL8( &binaryTemplateSpacings.tp, sqrt(uvar.mismatchT / diagTT));
  /* metric elements for eccentric case not considered? */

  /* Calculate SSB times (can do this once since search is currently only for one sky position, and binary doppler shift is added later) */
  if ((multiSSBTimes = XLALGetMultiSSBtimes ( multiStates, skyPos, dopplerpos.refTime, SSBPREC_RELATIVISTICOPT )) == NULL){
    LogPrintf ( LOG_CRITICAL, "%s: XLALGetMultiSSBtimes() failed with errno=%d\n", __func__, xlalErrno );
    XLAL_ERROR( XLAL_EFUNC );
  }

  /* Allocate structure for binary doppler-shifting information */
  if ((multiBinaryTimes = XLALDuplicateMultiSSBtimes ( multiSSBTimes )) == NULL){
    LogPrintf ( LOG_CRITICAL, "%s: XLALDuplicateMultiSSBtimes() failed with errno=%d\n", __func__, xlalErrno );
    XLAL_ERROR( XLAL_EFUNC );
  }

  UINT8 numSFTs = sftIndices->length;
  if ((shiftedFreqs = XLALCreateREAL8Vector ( numSFTs ) ) == NULL){
    LogPrintf ( LOG_CRITICAL, "%s: XLALCreateREAL8Vector() failed with errno=%d\n", __func__, xlalErrno );
    XLAL_ERROR( XLAL_EFUNC );
  }
  if ((lowestBins = XLALCreateUINT4Vector ( numSFTs ) ) == NULL){
    LogPrintf ( LOG_CRITICAL, "%s: XLALCreateUINT4Vector() failed with errno=%d\n", __func__, xlalErrno );
    XLAL_ERROR( XLAL_EFUNC );
  }
  if ((kappaValues = XLALCreateREAL8Vector ( numSFTs ) ) == NULL){
    LogPrintf ( LOG_CRITICAL, "%s: XLALCreateREAL8Vector() failed with errno=%d\n", __func__, xlalErrno );
    XLAL_ERROR( XLAL_EFUNC );
  }
  if ((signalPhases = XLALCreateREAL8Vector ( numSFTs ) ) == NULL){
    LogPrintf ( LOG_CRITICAL, "%s: XLALCreateREAL8Vector() failed with errno=%d\n", __func__, xlalErrno );
    XLAL_ERROR( XLAL_EFUNC );
  }

  /* args should be : spacings, min and max doppler params */
  while ( (GetNextCrossCorrTemplate(&dopplerShiftFlag, &dopplerpos, &binaryTemplateSpacings, &minBinaryTemplate, &maxBinaryTemplate) == 0) )
    {
      /* do useful stuff here*/

      /* Apply additional Doppler shifting using current binary orbital parameters */
      /* Might want to be clever about checking whether we've changed the orbital parameters or only the frequency */
      if (dopplerShiftFlag == TRUE)
	{
	  if ( (XLALAddMultiBinaryTimes( &multiBinaryTimes, multiSSBTimes, &dopplerpos )  != XLAL_SUCCESS ) ) {
	    LogPrintf ( LOG_CRITICAL, "%s: XLALAddMultiBinaryTimes() failed with errno=%d\n", __func__, xlalErrno );
	    XLAL_ERROR( XLAL_EFUNC );
	  }
	}


      if ( (XLALGetDopplerShiftedFrequencyInfo( shiftedFreqs, lowestBins, kappaValues, signalPhases, uvar.numBins, &dopplerpos, sftIndices, inputSFTs, multiBinaryTimes, Tsft )  != XLAL_SUCCESS ) ) {
	LogPrintf ( LOG_CRITICAL, "%s: XLALGetDopplerShiftedFrequencyInfo() failed with errno=%d\n", __func__, xlalErrno );
	XLAL_ERROR( XLAL_EFUNC );
      }

      if ( (XLALCalculatePulsarCrossCorrStatistic( &ccStat, &evSquared, curlyGUnshifted, signalPhases, lowestBins, kappaValues, uvar.numBins, sftPairs, sftIndices, inputSFTs, multiWeights)  != XLAL_SUCCESS ) ) {
	LogPrintf ( LOG_CRITICAL, "%s: XLALCalculateAveCrossCorrStatistic() failed with errno=%d\n", __func__, xlalErrno );
	XLAL_ERROR( XLAL_EFUNC );
      }

      /* fill candidate struct and insert into toplist if necessary */
      thisCandidate.freq = dopplerpos.fkdot[0];
      thisCandidate.tp = XLALGPSGetREAL8( &dopplerpos.tp );
      thisCandidate.argp = dopplerpos.argp;
      thisCandidate.asini = dopplerpos.asini;
      thisCandidate.ecc = dopplerpos.ecc;
      thisCandidate.period = dopplerpos.period;
      thisCandidate.rho = ccStat;
      thisCandidate.evSquared = evSquared;
      thisCandidate.estSens = estSens * estSens;

      insert_into_crossCorrBinary_toplist(ccToplist, thisCandidate);

    } /* end while loop over templates */


  /* write candidates to file */
  sort_crossCorrBinary_toplist( ccToplist );
  /* add error checking */

  final_write_crossCorrBinary_toplist_to_file( ccToplist, uvar.toplistFilename, &checksum);


  XLALDestroyREAL8Vector ( signalPhases );
  XLALDestroyREAL8Vector ( kappaValues );
  XLALDestroyUINT4Vector ( lowestBins );
  XLALDestroyREAL8Vector ( shiftedFreqs );
  XLALDestroyMultiSSBtimes ( multiBinaryTimes );
  XLALDestroyMultiSSBtimes ( multiSSBTimes );
  XLALDestroyREAL8Vector ( curlyGUnshifted );
  XLALDestroySFTPairIndexList( sftPairs );
  XLALDestroySFTIndexList( sftIndices );
  XLALDestroyMultiAMCoeffs ( multiCoeffs );
  XLALDestroyMultiDetectorStateSeries ( multiStates );
  XLALDestroyMultiTimestamps ( multiTimes );
  XLALDestroyMultiNoiseWeights ( multiWeights );
  XLALDestroyMultiPSDVector ( multiPSDs );
  XLALDestroyMultiSFTVector ( inputSFTs );

  /* de-allocate memory for configuration variables */
  XLALDestroyConfigVars ( &config );

  /* de-allocate memory for user input variables */
  XLALDestroyUserVars();

  /* free toplist memory */
  free_crossCorr_toplist(&ccToplist);

  /* check memory leaks if we forgot to de-allocate anything */
  LALCheckMemoryLeaks();

  LogPrintf (LOG_CRITICAL, "The metric element g_ff=%.9f\n", diagff);
  LogPrintf (LOG_CRITICAL, "The metric element g_aa=%.9f\n", diagaa);
  LogPrintf (LOG_CRITICAL, "The metric element g_TT=%.9f\n", diagTT);
  LogPrintf (LOG_CRITICAL, "End time\n", 0);/*for debug convenience to record calculating time*/
  return 0;

} /* main */


/* initialize and register user variables */
int XLALInitUserVars (UserInput_t *uvar)
{

  /* initialize with some defaults */
  uvar->help = FALSE;
  uvar->startTime = 814838413;	/* 1 Nov 2005, ~ start of S5 */
  uvar->endTime = uvar->startTime + (INT4) round ( LAL_YRSID_SI ) ;	/* 1 year of data */
  uvar->maxLag = 0.0;
  uvar->inclAutoCorr = FALSE;
  uvar->fStart = 100.0;
  uvar->fBand = 0.1;
  /* uvar->fdotStart = 0.0; */
  /* uvar->fdotBand = 0.0; */
  uvar->alphaRad = 0.0;
  uvar->deltaRad = 0.0;
  uvar->rngMedBlock = 50;
  uvar->numBins = 1;
  uvar->readlist = FALSE;
  /* default for reftime is in the middle */
  uvar->refTime = 0.5*(uvar->startTime + uvar->endTime);

  /* zero binary orbital parameters means not a binary */
  uvar->orbitAsiniSec = 0.0;
  uvar->orbitAsiniSecBand = 0.0;
  uvar->orbitPSec = 0.0;
  uvar->orbitTimeAsc = 0;
  uvar->orbitTimeAscBand = 0;

  /*default mismatch values */
  /* set to 0.1 by default -- for no real reason */
  /* make 0.1 a macro? */
  uvar->mismatchF = 0.1;
  uvar->mismatchA = 0.1;
  uvar->mismatchT = 0.1;
  uvar->mismatchP = 0.1;

  uvar->ephemEarth = XLALStringDuplicate("earth00-19-DE405.dat.gz");
  uvar->ephemSun = XLALStringDuplicate("sun00-19-DE405.dat.gz");

  uvar->sftLocation = XLALCalloc(1, MAXFILENAMELENGTH+1);

  /* initialize number of candidates in toplist -- default is just to return the single best candidate */
  uvar->numCand = 1;
  uvar->toplistFilename = XLALStringDuplicate("toplist_crosscorr.dat");


  /* register  user-variables */
  XLALregBOOLUserStruct  ( help, 	   'h',  UVAR_HELP, "Print this message");
  XLALregINTUserStruct   ( startTime,       0,  UVAR_OPTIONAL, "Desired start time of analysis in GPS seconds");
  XLALregINTUserStruct   ( endTime,         0,  UVAR_OPTIONAL, "Desired end time of analysis in GPS seconds");
  XLALregREALUserStruct  ( maxLag,          0,  UVAR_OPTIONAL, "Maximum lag time in seconds between SFTs in correlation");
  XLALregBOOLUserStruct  ( inclAutoCorr,    0,  UVAR_OPTIONAL, "Include auto-correlation terms (an SFT with itself)");
  XLALregREALUserStruct  ( fStart,          0,  UVAR_OPTIONAL, "Start frequency in Hz");
  XLALregREALUserStruct  ( fBand,           0,  UVAR_OPTIONAL, "Frequency band to search over in Hz ");
  /* XLALregREALUserStruct  ( fdotStart,     0,  UVAR_OPTIONAL, "Start value of spindown in Hz/s"); */
  /* XLALregREALUserStruct  ( fdotBand,      0,  UVAR_OPTIONAL, "Band for spindown values in Hz/s"); */
  XLALregREALUserStruct  ( alphaRad,        0,  UVAR_OPTIONAL, "Right ascension for directed search (radians)");
  XLALregREALUserStruct  ( deltaRad,        0,  UVAR_OPTIONAL, "Declination for directed search (radians)");
  XLALregREALUserStruct  ( refTime,         0,  UVAR_OPTIONAL, "SSB reference time for pulsar-parameters [Default: midPoint]");
  XLALregREALUserStruct  ( orbitAsiniSec,   0,  UVAR_OPTIONAL, "Start of search band for projected semimajor axis (seconds) [0 means not a binary]");
  XLALregREALUserStruct  ( orbitAsiniSecBand, 0,  UVAR_OPTIONAL, "Width of search band for projected semimajor axis (seconds)");
  XLALregREALUserStruct  ( orbitPSec,       0,  UVAR_OPTIONAL, "Binary orbital period (seconds) [0 means not a binary]");
  XLALregREALUserStruct  ( orbitTimeAsc,    0,  UVAR_OPTIONAL, "Start of orbital time-of-ascension band in GPS seconds");
  XLALregREALUserStruct  ( orbitTimeAscBand, 0,  UVAR_OPTIONAL, "Width of orbital time-of-ascension band (seconds)");
  XLALregSTRINGUserStruct( ephemEarth,      0,  UVAR_OPTIONAL, "Earth ephemeris file to use");
  XLALregSTRINGUserStruct( ephemSun,        0,  UVAR_OPTIONAL, "Sun ephemeris file to use");
  XLALregSTRINGUserStruct( sftLocation,     0,  UVAR_REQUIRED, "Filename pattern for locating SFT data");
  XLALregINTUserStruct   ( rngMedBlock,     0,  UVAR_OPTIONAL, "Running median block size for PSD estimation");
  XLALregINTUserStruct   ( numBins,         0,  UVAR_OPTIONAL, "Number of frequency bins to include in calculation");
  XLALregREALUserStruct  ( mismatchF,       0,  UVAR_OPTIONAL, "Desired mismatch for frequency spacing");
  XLALregREALUserStruct  ( mismatchA,       0,  UVAR_OPTIONAL, "Desired mismatch for asini spacing");
  XLALregREALUserStruct  ( mismatchT,       0,  UVAR_OPTIONAL, "Desired mismatch for periapse passage time spacing");
  XLALregREALUserStruct  ( mismatchP,       0,  UVAR_OPTIONAL, "Desired mismatch for period spacing");
  XLALregINTUserStruct   ( numCand,         0,  UVAR_OPTIONAL, "Number of candidates to keep in toplist");
  XLALregSTRINGUserStruct( toplistFilename, 0,  UVAR_OPTIONAL, "Output filename containing candidates in toplist");
  XLALregBOOLUserStruct  ( readlist,        0,  UVAR_OPTIONAL, "Read in the SFT-pair list(FALSE) or output(TRUE)");

  if ( xlalErrno ) {
    XLALPrintError ("%s: user variable initialization failed with errno = %d.\n", __func__, xlalErrno );
    XLAL_ERROR ( XLAL_EFUNC );
  }

  return XLAL_SUCCESS;
}



/* initialize and register user variables */
int XLALInitializeConfigVars (ConfigVariables *config, const UserInput_t *uvar)
{

  static SFTConstraints constraints;
  LIGOTimeGPS startTime, endTime;

  /* set sft catalog constraints */
  constraints.detector = NULL;
  constraints.timestamps = NULL;
  constraints.minStartTime = &startTime;
  constraints.maxStartTime = &endTime;
  XLALGPSSet( constraints.minStartTime, uvar->startTime, 0);
  XLALGPSSet( constraints.maxStartTime, uvar->endTime,0);

  /* This check doesn't seem to work, since XLALGPSSet doesn't set its
     first argument.

     if ( (constraints.minStartTime == NULL)&& (constraints.maxStartTime == NULL) ) {
     LogPrintf ( LOG_CRITICAL, "%s: XLALGPSSet() failed with errno=%d\n", __func__, xlalErrno );
     XLAL_ERROR( XLAL_EFUNC );
     }

  */

  /* get catalog of SFTs */
  if ((config->catalog = XLALSFTdataFind (uvar->sftLocation, &constraints)) == NULL){
    LogPrintf ( LOG_CRITICAL, "%s: XLALSFTdataFind() failed with errno=%d\n", __func__, xlalErrno );
    XLAL_ERROR( XLAL_EFUNC );
  }

  /* initialize ephemeris data*/
  XLAL_CHECK ( (config->edat = XLALInitBarycenter ( uvar->ephemEarth, uvar->ephemSun )) != NULL, XLAL_EFUNC );

  return XLAL_SUCCESS;

}
/* XLALInitializeConfigVars() */

/* deallocate memory associated with config variables */
int XLALDestroyConfigVars (ConfigVariables *config)
{
  XLALDestroySFTCatalog(config->catalog);
  XLALDestroyEphemerisData(config->edat);
  return XLAL_SUCCESS;
}
/* XLALDestroyConfigVars() */

/* getting the next template */
/** FIXME: spacings and min, max values of binary parameters are not used yet */


int GetNextCrossCorrTemplate(BOOLEAN *binaryParamsFlag, PulsarDopplerParams *dopplerpos, PulsarDopplerParams *binaryTemplateSpacings, PulsarDopplerParams *minBinaryTemplate, PulsarDopplerParams *maxBinaryTemplate)
{

  REAL8 new_freq = dopplerpos->fkdot[0];
  REAL8 new_asini = dopplerpos->asini;
  REAL8 new_tp = XLALGPSGetREAL8(&(dopplerpos->tp));
  REAL8 tp_hi = XLALGPSGetREAL8(&(maxBinaryTemplate->tp));

  /* basic sanity checks */
  if (binaryTemplateSpacings == NULL)
    return -1;

  if (minBinaryTemplate == NULL)
    return -1;

  if (maxBinaryTemplate == NULL)
    return -1;

  /* check spacings not negative */

  if (new_freq < maxBinaryTemplate->fkdot[0])    /*loop over f at first*/
    {
      new_freq = dopplerpos->fkdot[0] + binaryTemplateSpacings->fkdot[0];
      dopplerpos->fkdot[0] = new_freq;
      *binaryParamsFlag = FALSE;
      return 0;
    }
  else
    {
      if (new_asini < maxBinaryTemplate->asini)  /*after looping all f, initialize f and loop over a_p*/
	{
	  new_asini = dopplerpos->asini + binaryTemplateSpacings->asini;
	  dopplerpos->asini = new_asini;
	  new_freq = minBinaryTemplate->fkdot[0];
	  dopplerpos->fkdot[0] = new_freq;
	  *binaryParamsFlag = TRUE;
	  return 0;
	}
      else
	{
	  if (new_tp < tp_hi)                    /*after looping the plane of f and a_p, initialize f, a_p and loop over T*/
	    {
	      new_tp = XLALGPSGetREAL8(XLALGPSAddGPS(&(dopplerpos->tp), &(binaryTemplateSpacings->tp)));
	      XLALGPSSetREAL8(&(dopplerpos->tp),new_tp);
	      new_asini = minBinaryTemplate->asini;
	      dopplerpos->asini = new_asini;
	      new_freq = minBinaryTemplate->fkdot[0];
	      dopplerpos->fkdot[0] = new_freq;
	      *binaryParamsFlag = TRUE;
	      return 0;
	    }
	  else
	    {
	      return 1;
	    }
	}
    }
}
