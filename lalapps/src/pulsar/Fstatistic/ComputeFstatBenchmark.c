/*
*  Copyright (C) 2015 Reinhard Prix
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

#include <lalapps.h>
#include <lal/XLALError.h>
#include <lal/LALBarycenter.h>
#include <lal/LALInitBarycenter.h>
#include <lal/LogPrintf.h>
#include <lal/CWMakeFakeData.h>
#include <lal/LALConstants.h>
#include <lal/ExtrapolatePulsarSpins.h>
#include <lal/ComputeFstat.h>
#include <lal/DetectorStates.h>
#include <lal/LFTandTSutils.h>
#include <lal/LALString.h>
#include <lal/UserInput.h>
#include <lal/LALPulsarVCSInfo.h>

// benchmark ComputeFstat() functions for performance and memory usage
REAL8 XLALGetCurrentHeapUsageMB ( void );

typedef struct
{
  int FstatMethod;		//!< select which method/algorithm to use to compute the F-statistic
  REAL8Range Alpha;
  REAL8Range Delta;
  REAL8Range Freq;
  REAL8Range f1dot;
  REAL8Range f2dot;
  REAL8Range FreqResolution;
  INT4Range numFreqBins;
  INT4Range Tseg;
  INT4 numSegments;
  LALStringVector *IFOs;
  CHAR *outputInfo;
  INT4 numTrials;
  LIGOTimeGPS startTime;

  // ----- developer options
  CHAR *ephemEarth;		/**< Earth ephemeris file to use */
  CHAR *ephemSun;		/**< Sun ephemeris file to use */

  REAL8 Tsft;
  BOOLEAN sharedWorkspace;   	// useful for checking workspace sharing for Resampling
  BOOLEAN perSegmentSFTs;     	// Weave vs GCT convention: GCT loads SFT frequency ranges globally, Weave loads them per segment (more efficient)
  BOOLEAN resampFFTPowerOf2;
  INT4 Dterms;
  INT4 randSeed;

  BOOLEAN version;	// output code version
} UserInput_t;

// ---------- main ----------
int
main ( int argc, char *argv[] )
{

  CHAR *VCSInfoString;
  XLAL_CHECK_MAIN ( (VCSInfoString = XLALGetVersionString(0)) != NULL, XLAL_EFUNC );

  // ---------- handle user input ----------
  UserInput_t XLAL_INIT_DECL(uvar_s);
  UserInput_t *uvar = &uvar_s;

  uvar->FstatMethod = FMETHOD_RESAMP_BEST;
  uvar->Alpha[0] = 0;
  uvar->Alpha[1] = LAL_TWOPI;
  uvar->Delta[0] = -LAL_PI/2;
  uvar->Delta[1] =  LAL_PI/2;
  uvar->Freq[0] = 100;
  uvar->Freq[1] = 1000;
  uvar->f1dot[0] = -3e-9;
  uvar->f1dot[1] = 0;
  uvar->f2dot[0] = 0;
  uvar->f2dot[1] = 1e-16;
  uvar->FreqResolution[0] = 0.1;
  uvar->FreqResolution[1] = 1;
  uvar->numFreqBins[0] = 1000;
  uvar->numFreqBins[1] = 100000;
  uvar->Tseg[0] = 10 * 3600;
  uvar->Tseg[1] = 250 * 3600;

  uvar->numSegments = 90;
  uvar->numTrials = 1;
  uvar->startTime.gpsSeconds = 711595934;
  uvar->Tsft = 1800;
  uvar->sharedWorkspace = 1;
  uvar->resampFFTPowerOf2 = 1;
  uvar->perSegmentSFTs = 1;

  uvar->Dterms = 8;

  uvar->ephemEarth = XLALStringDuplicate("earth00-19-DE405.dat.gz");
  uvar->ephemSun = XLALStringDuplicate("sun00-19-DE405.dat.gz");
  uvar->randSeed = 1;

  XLAL_CHECK_MAIN ( (uvar->IFOs = XLALCreateStringVector ( "H1", NULL )) != NULL, XLAL_EFUNC );
  uvar->outputInfo = NULL;

  XLAL_CHECK ( XLALRegisterUvarAuxDataMember ( FstatMethod, UserEnum, XLALFstatMethodChoices(), 0, OPTIONAL, "F-statistic method to use" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN ( XLALRegisterUvarMember ( Alpha,          RAJRange,       0, OPTIONAL,  "Skyposition [drawn isotropically]: Range in 'Alpha' = right ascension)" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN ( XLALRegisterUvarMember ( Delta,          DECJRange,      0, OPTIONAL,  "Skyposition [drawn isotropically]: Range in 'Delta' = declination" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN ( XLALRegisterUvarMember ( Freq,           REAL8Range,     0, OPTIONAL,  "Search frequency in Hz [range to draw from]" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN ( XLALRegisterUvarMember ( f1dot,          REAL8Range,     0, OPTIONAL,  "Search 1st spindown f1dot in Hz/s [range to draw from]" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN ( XLALRegisterUvarMember ( f2dot,          REAL8Range,     0, OPTIONAL,  "Search 2nd spindown f2dot in Hz/s^2 [range to draw from]" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN ( XLALRegisterUvarMember ( FreqResolution, REAL8Range,     0, OPTIONAL,  "Frequency resolution 'R' in natural units 1/Tseg such that: dFreq = R/Tseg) [range to draw from]" ) == XLAL_SUCCESS, XLAL_EFUNC );

  XLAL_CHECK_MAIN ( XLALRegisterUvarMember ( startTime,      EPOCH,          0, OPTIONAL,  "Start time of first segment" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN ( XLALRegisterUvarMember ( Tseg,           INT4Range,      0, OPTIONAL,  "Coherent segment length in seconds [range to draw from]" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN ( XLALRegisterUvarMember ( numSegments,    INT4,           0, OPTIONAL,  "Number of semi-coherent segments" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN ( XLALRegisterUvarMember ( numFreqBins,    INT4Range,      0, OPTIONAL,  "Number of frequency bins to search [range to draw from]" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN ( XLALRegisterUvarMember ( IFOs,           STRINGVector,   0, OPTIONAL,  "IFOs to use [list]" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN ( XLALRegisterUvarMember ( numTrials,      INT4,           0, OPTIONAL,  "Number of repeated trials to run (with randomized parameters drawn from ranges)" ) == XLAL_SUCCESS, XLAL_EFUNC );

  XLAL_CHECK_MAIN ( XLALRegisterUvarMember ( sharedWorkspace,BOOLEAN,        0, OPTIONAL,  "Use workspace sharing across segments (only used in Resampling)" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN ( XLALRegisterUvarMember ( perSegmentSFTs, BOOLEAN,        0, OPTIONAL,  "Weave vs GCT: GCT determines and loads SFT frequency ranges globally, Weave does that per segment (more efficient)" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN ( XLALRegisterUvarMember ( resampFFTPowerOf2, BOOLEAN,     0, OPTIONAL,  "For Resampling methods: enforce FFT length to be a power of two (by rounding up)" ) == XLAL_SUCCESS, XLAL_EFUNC );

  XLAL_CHECK_MAIN ( XLALRegisterUvarMember ( Dterms,         INT4,           0, OPTIONAL,  "Number of kernel terms (single-sided) in\na) Dirichlet kernel if FstatMethod=Demod*\nb) sinc-interpolation if FstatMethod=Resamp*" ) == XLAL_SUCCESS, XLAL_EFUNC );

  XLAL_CHECK_MAIN ( XLALRegisterUvarMember ( outputInfo,     STRING,         0, OPTIONAL,  "Append Resampling internal info into this file") == XLAL_SUCCESS, XLAL_EFUNC );

  XLAL_CHECK_MAIN ( XLALRegisterUvarMember ( Tsft,           REAL8,          0, DEVELOPER, "SFT length" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN ( XLALRegisterUvarMember ( ephemEarth,     STRING,         0, DEVELOPER, "Earth ephemeris file to use") == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN ( XLALRegisterUvarMember ( ephemSun,       STRING,         0, DEVELOPER, "Sun ephemeris file to use") == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN ( XLALRegisterUvarMember ( randSeed,       INT4,           0, DEVELOPER, "Random seed to use for rand() to draw randomized parameters.") == XLAL_SUCCESS, XLAL_EFUNC );

  BOOLEAN should_exit = 0;
  XLAL_CHECK_MAIN ( XLALUserVarReadAllInput( &should_exit, argc, argv, lalPulsarVCSInfoList ) == XLAL_SUCCESS, XLAL_EFUNC );
  if ( should_exit ) {
    return EXIT_FAILURE;
  }

  // produce log-string (for output-file headers)
  CHAR *logstring = NULL;
  XLAL_CHECK_MAIN ( (logstring = XLALStringAppend ( logstring, "%% cmdline: " )) != NULL, XLAL_EFUNC );
  CHAR *cmdline;
  XLAL_CHECK_MAIN ( (cmdline = XLALUserVarGetLog ( UVAR_LOGFMT_CMDLINE )) != NULL, XLAL_EFUNC );
  XLAL_CHECK_MAIN ( (logstring = XLALStringAppend ( logstring, cmdline )) != NULL, XLAL_EFUNC );
  XLALFree ( cmdline );
  XLAL_CHECK_MAIN ( (logstring = XLALStringAppend ( logstring, "\n" )) != NULL, XLAL_EFUNC );
  XLAL_CHECK_MAIN ( (logstring = XLALStringAppend ( logstring, VCSInfoString )) != NULL, XLAL_EFUNC );

  // check user input
  XLAL_CHECK_MAIN ( uvar->numSegments >= 1, XLAL_EINVAL );
  XLAL_CHECK_MAIN ( uvar->Tsft > 1, XLAL_EINVAL );
  XLAL_CHECK_MAIN ( uvar->numTrials >= 1, XLAL_EINVAL );
  // ---------- end: handle user input ----------
  srand( uvar->randSeed );	// set random seed

  // common setup over repeated trials
  EphemerisData *ephem;
  XLAL_CHECK_MAIN ( (ephem = XLALInitBarycenter ( uvar->ephemEarth, uvar->ephemSun )) != NULL, XLAL_EFUNC );
  REAL8 memBase = XLALGetCurrentHeapUsageMB();

  UINT4 numDetectors = uvar->IFOs->length;
  // ----- setup injection and data parameters
  LIGOTimeGPSVector *startTime_l, *endTime_l;
  XLAL_CHECK_MAIN ( (startTime_l = XLALCreateTimestampVector ( uvar->numSegments )) != NULL, XLAL_EFUNC );
  XLAL_CHECK_MAIN ( (endTime_l = XLALCreateTimestampVector ( uvar->numSegments )) != NULL, XLAL_EFUNC );

  // ----- setup optional Fstat arguments
  FstatOptionalArgs optionalArgs = FstatOptionalArgsDefaults;
  MultiNoiseFloor XLAL_INIT_DECL(injectSqrtSX);
  injectSqrtSX.length = numDetectors;
  for ( UINT4 X=0; X < numDetectors; X ++ ) {
    injectSqrtSX.sqrtSn[X] = 1;
  }
  optionalArgs.injectSqrtSX = &injectSqrtSX;
  optionalArgs.FstatMethod = uvar->FstatMethod;
  optionalArgs.collectTiming = 1;
  optionalArgs.resampFFTPowerOf2 = uvar->resampFFTPowerOf2;
  optionalArgs.Dterms = uvar->Dterms;

  FILE *timingLogFILE = NULL;
  FILE *timingParFILE = NULL;
  if ( uvar->outputInfo != NULL )
    {
      char *parFname = NULL;
      XLAL_CHECK_MAIN ( (parFname = XLALStringAppend ( parFname, uvar->outputInfo )) != NULL, XLAL_EFUNC );
      XLAL_CHECK_MAIN ( (parFname = XLALStringAppend ( parFname, ".pars" )) != NULL, XLAL_EFUNC );

      XLAL_CHECK_MAIN ( (timingLogFILE = fopen (uvar->outputInfo, "ab")) != NULL, XLAL_ESYS, "Failed to open '%s' for appending\n", uvar->outputInfo );
      XLAL_CHECK_MAIN ( (timingParFILE = fopen (parFname, "ab"))         != NULL, XLAL_ESYS, "Failed to open '%s' for appending\n", parFname );
      XLALFree ( parFname );
      fprintf ( timingLogFILE, "%s\n", logstring );
      fprintf ( timingParFILE, "%s\n", logstring );
      fprintf ( timingParFILE, "%%%%%8s %10s %12s %12s %12s %12s %12s %12s %12s %12s\n",
                "Nseg", "Tseg", "Freq", "FreqBand", "dFreq", "f1dot", "f2dot", "Alpha", "Delta", "memUsageMB" );
    }
  FstatInputVector *inputs;
  FstatQuantities whatToCompute = (FSTATQ_2F | FSTATQ_2F_PER_DET);
  FstatResults *results = NULL;

#define drawFromREAL8Range(range) (range[0] + (range[1] - range[0]) * rand() / RAND_MAX )
#define drawFromINT4Range(range)  (range[0] + (INT4)round(1.0*(range[1] - range[0]) * rand() / RAND_MAX) )
  // ---------- main loop over repeated trials: randomize uniformly over input ranges  ----------
  for ( INT4 i = 0; i < uvar->numTrials; i ++ )
    {
      UINT4 Tseg_i        = drawFromINT4Range ( uvar->Tseg );

      SFTCatalog **catalogs;
      XLAL_CHECK_MAIN ( (catalogs = XLALCalloc ( uvar->numSegments, sizeof( catalogs[0] ))) != NULL, XLAL_ENOMEM );
      for ( INT4 l = 0; l < uvar->numSegments; l ++ )
        {
          startTime_l->data[l] = (l==0)? uvar->startTime : endTime_l->data[l-1];
          endTime_l->data[l]   = startTime_l->data[l];
          endTime_l->data[l].gpsSeconds += Tseg_i;

          MultiLIGOTimeGPSVector *multiTimestamps;
          XLAL_CHECK_MAIN ( (multiTimestamps = XLALMakeMultiTimestamps ( startTime_l->data[l], Tseg_i, uvar->Tsft, 0, numDetectors )) != NULL, XLAL_EFUNC );
          XLAL_CHECK_MAIN ( (catalogs[l] = XLALMultiAddToFakeSFTCatalog ( NULL, uvar->IFOs, multiTimestamps )) != NULL, XLAL_EFUNC );
          XLALDestroyMultiTimestamps ( multiTimestamps );
          startTime_l->data[l] = endTime_l->data[l];
        } // for l < numSegments

      PulsarSpinRange XLAL_INIT_DECL(spinRange_i);

      LIGOTimeGPS refTime = { uvar->startTime.gpsSeconds + 0.5 * uvar->numSegments * Tseg_i, 0 };
      spinRange_i.refTime = refTime;
      spinRange_i.fkdot[0] = drawFromREAL8Range ( uvar->Freq );
      spinRange_i.fkdot[1] = drawFromREAL8Range ( uvar->f1dot );
      spinRange_i.fkdot[2] = drawFromREAL8Range ( uvar->f2dot );
      REAL8 asini = 0, period = 0, ecc = 0, argp = 0;
      LIGOTimeGPS XLAL_INIT_DECL(tp);

      PulsarDopplerParams XLAL_INIT_DECL(Doppler_i);
      Doppler_i.refTime = refTime;
      Doppler_i.Alpha = drawFromREAL8Range ( uvar->Alpha );
      REAL8Range sDeltaRange;
      sDeltaRange[0] = sin ( uvar->Delta[0] );
      sDeltaRange[1] = sin ( uvar->Delta[1] );
      Doppler_i.Delta = asin ( drawFromREAL8Range ( sDeltaRange ) );
      memcpy ( &Doppler_i.fkdot, &spinRange_i.fkdot, sizeof(Doppler_i.fkdot) );;
      // not allowing to randomize or control binary-orbital parameters yet
      Doppler_i.period = period;
      Doppler_i.ecc = ecc;
      Doppler_i.asini = asini;
      Doppler_i.tp = tp;
      Doppler_i.argp = argp;

      UINT4 numFreqBins_i    = drawFromINT4Range ( uvar->numFreqBins );
      REAL8 FreqResolution_i = drawFromREAL8Range ( uvar->FreqResolution );

      REAL8 dFreq_i          = FreqResolution_i / Tseg_i;
      REAL8 FreqBand_i       = numFreqBins_i * dFreq_i;

      XLAL_CHECK_MAIN ( (inputs = XLALCreateFstatInputVector ( uvar->numSegments )) != NULL, XLAL_EFUNC );

      fprintf ( stderr, "trial %d/%d: Tseg = %.1f d, numSegments = %d, Alpha = %.2f rad, Delta = %.2f rad, Freq = %.6f Hz, f1dot = %.1e Hz/s, f2dot = %.1e Hz/s^2, R = %.2f, numFreqBins = %d [dFreq = %.2e Hz, FreqBand = %.2e Hz]\n",
                i+1, uvar->numTrials, Tseg_i / 86400.0, uvar->numSegments, Doppler_i.Alpha, Doppler_i.Delta, Doppler_i.fkdot[0], Doppler_i.fkdot[1], Doppler_i.fkdot[2], FreqResolution_i, numFreqBins_i, dFreq_i, FreqBand_i );

      spinRange_i.fkdotBand[0] = FreqBand_i;
      REAL8 minCoverFreq_il, maxCoverFreq_il;
      // GCT convention: determine global SFT frequency band for all segments
      if ( ! uvar->perSegmentSFTs ) {
        XLAL_CHECK_MAIN ( XLALCWSignalCoveringBand ( &minCoverFreq_il, &maxCoverFreq_il, &startTime_l->data[0], &endTime_l->data[uvar->numSegments-1], &spinRange_i, Doppler_i.asini, Doppler_i.period, Doppler_i.ecc ) == XLAL_SUCCESS, XLAL_EFUNC );
      }

      // create per-segment input structs
      for ( INT4 l = 0; l < uvar->numSegments; l ++ )
        {
          if ( uvar->sharedWorkspace && l > 0 ) {
            optionalArgs.prevInput = inputs->data[0];
          } else {
            optionalArgs.prevInput = NULL;
          }
          // Weave convention: determine per-segment SFT frequency band
          if ( uvar->perSegmentSFTs ) {
            XLAL_CHECK_MAIN ( XLALCWSignalCoveringBand ( &minCoverFreq_il, &maxCoverFreq_il, &startTime_l->data[l], &endTime_l->data[l], &spinRange_i, Doppler_i.asini, Doppler_i.period, Doppler_i.ecc ) == XLAL_SUCCESS, XLAL_EFUNC );
          }
          XLAL_CHECK_MAIN ( (inputs->data[l] = XLALCreateFstatInput ( catalogs[l], minCoverFreq_il, maxCoverFreq_il, dFreq_i, ephem, &optionalArgs )) != NULL, XLAL_EFUNC );
        }
      for ( INT4 l = 0; l < uvar->numSegments; l ++ ) {
        XLALDestroySFTCatalog ( catalogs[l] );
      }
      XLALFree ( catalogs );

      // ----- compute Fstatistics over segments
      for ( INT4 l = 0; l < uvar->numSegments; l ++ )
        {
          XLAL_CHECK_MAIN ( XLALComputeFstat ( &results, inputs->data[l], &Doppler_i, numFreqBins_i, whatToCompute ) == XLAL_SUCCESS, XLAL_EFUNC );

          // ----- output timing details to file if requested
          if ( timingLogFILE != NULL ) {
            XLAL_CHECK_MAIN ( XLALAppendFstatTiming2File ( inputs->data[l], timingLogFILE, (l == 0) && (i==0)) == XLAL_SUCCESS, XLAL_EFUNC );
          }
        } // for l < numSegments

      REAL8 memEnd = XLALGetCurrentHeapUsageMB();
      REAL8 memUsage = memEnd - memBase;
      const char *FmethodName = XLALGetFstatInputMethodName ( inputs->data[0] );
      fprintf (stderr, "%-15s: memoryUsage = %6.1f MB\n", FmethodName, memUsage );

      if ( timingParFILE != NULL )
        {
          fprintf ( timingParFILE, "%10d %10d %12g %12g %12g %12g %12g %12g %12g %12g\n",
                    uvar->numSegments, Tseg_i, Doppler_i.fkdot[0], FreqBand_i, dFreq_i, Doppler_i.fkdot[1], Doppler_i.fkdot[2], Doppler_i.Alpha, Doppler_i.Delta, memUsage
                    );
        }

      fprintf (stderr, "%-15s: memoryUsage = %6.1f MB\n", FmethodName, memEnd - memBase );
      XLALDestroyFstatInputVector ( inputs );
    } // for i < numTrials

  // ----- free memory ----------
  if ( timingLogFILE != NULL ) {
    fclose ( timingLogFILE );
  }
  if ( timingParFILE != NULL ) {
    fclose ( timingParFILE );
  }

  XLALDestroyFstatResults ( results );
  XLALDestroyUserVars();
  XLALDestroyEphemerisData ( ephem );
  XLALFree ( VCSInfoString );
  XLALFree ( logstring );
  XLALDestroyTimestampVector ( startTime_l );
  XLALDestroyTimestampVector ( endTime_l );

  LALCheckMemoryLeaks();

  return XLAL_SUCCESS;

} // main()


// --------------------------------------------------------------------------------
// code to read current process RSS memory usage from /proc, taken from
// https://stackoverflow.com/questions/63166/how-to-determine-cpu-and-memory-consumption-from-inside-a-process
static int parseLine(char* line);

REAL8
XLALGetCurrentHeapUsageMB ( void )
{
  FILE* file = fopen("/proc/self/status", "r");
  int result = -1;
  char line[128];
  if ( file == NULL ) {
    return result;
  }
  while ( fgets ( line, sizeof(line), file) != NULL )
    {
      if (strncmp(line, "VmRSS:", 6) == 0){
        result = parseLine(line);
        break;
      }
    }
  fclose(file);
  return (result / 1024.0);
} // XLALGetCurrentHeapUsageMB()

static int parseLine(char* line)
{
  // This assumes that a digit will be found and the line ends in " Kb".
  int i = strlen(line);
  const char* p = line;
  while (*p <'0' || *p > '9') p++;
  line[i-3] = '\0';
  i = atoi(p);
  return i;
}
// --------------------------------------------------------------------------------
