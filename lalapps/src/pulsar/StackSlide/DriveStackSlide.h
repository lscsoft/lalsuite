/*********************************************************************************/
/*                                                                               */
/* File: DriveStackSlide.h                                                       */
/* Purpose: Contains driver functions for StackSlide search                      */
/* Author: Mendell, G. and Landry, M.                                            */
/* Started: 03 December 2003                                                     */
/* $Id$               */
/*                                                                               */
/*********************************************************************************/
/* Notes: */
/* BLK = one Block of frequency domain data */
/* STK = one Stack of frequency domain data; Block are coherenty combined to make Stacks. */
/* SUM = Summed Stacks after Sliding */

/* REVISIONS: */
/* 01/06/04 gam; remove extraneous code that is now in StackSlide.c */
/* 01/14/04 gam; Change threshold1 and threshold2 from type REAL8 to type REAL4 */
/* 01/14/04 gam; Add code that searches SUMs for peaks */
/* 01/20/04 gam; Change findStackSlidePeaks to LALFindStackSlidePeaks; put params into struct */
/* 01/20/04 gam; Change unusedFlag2 to threshold3. */
/* 01/21/04 gam; Added BOOLEAN overlap_event to SnglStackSlidePeriodicTable struct; */
/* 01/27/04 gam; Change unusedFlag1, unusedParam1 and unusedParam2 to threshold4, threshold5, and threshold6. */
/* 01/31/04 gam; If (params->parameterSpaceFlag == 0) then set up sky with deltaRA(DEC) = deltaRA(0)/cos(DEC), where deltaRA(0) = params->deltaRA */
/* 02/02/04 gam; Make code clear that only currently supporting the case params->parameterSpaceFlag == 0. */
/* 02/04/04 gam; Add code to store peaks in the SnglStackSlidePeriodicTable struct and write these out in xml format. */
/* 02/09/04 gam; Clean up SnglStackSlidePeriodicTable and xml storage to remove unused or repetative columns */
/* 02/09/04 gam; Add xmlStream to params and output other tables into xml file */
/* 02/11/04 gam; Comment out ComputePowerStats, SUMStats, and FreqSeriesPowerStats */
/* 02/11/04 gam; Remove obsolete code associate with INCLUDE_BUILDSUMFRAME_CODE (BuildSUMFrameOutput) and INCLUDE_BUILDOUTPUT_CODE (BuildOutput). */
/* 02/11/04 gam; LALInitBarycenter lal support package rather than internal version */
/* 02/11/04 gam; Use power rather than amplitudes when finding peaks */
/* 02/12/04 gam; Add num_sums = numSUMsTotal, freq_index =f0SUM*params->tEffSUM, and num_bins = nBinsPerSUM to StackSlideSFTsSearchSummaryTable */
/*               Add freq_index to SnglStackSlidePeriodicTable.  Each of these are integers that can be stored exactly */
/* 02/17/04 gam; Implement new rules for thresholdFlag. */
/* 02/17/04 gam; Change outputFlag to outputSUMFlag; add outputEventFlag and keepThisNumber  */
/*               (Remove windowFilterParam2 and windowFilterParam3)                          */
/*               Implement new rules for outputSUMFlag, outputEventFlag. and keepThisNumber. */
/*               Add loudestPeaksArray:                                                      */
/*               Keep the loudest event for every eventBandWidth = params->bandSUM/keepThisNumber, which is the same as */
/*               keep the loudest event for every nBinsPerOutputEvent = params->nBinsPerSUM/keepThisNumber. */
/* 02/19/04 gam; make maxWidth type REAL8. */
/* 02/20/04 gam; In InitializePeaksArray, initialize frequency to frequency of first bin for each event */
/* 02/20/04 gam; change threshold6 to INT4 maxWidthBin; remove maxWidth; replace fMax and fMin with ifMax and ifMin */
/*               use threshold 3 as fracPwDrop; threshold 4 and 5 are currently unused */
/* 02/20/04 gam; Replace thresholdFlag with vetoWidePeaks and vetoOverlapPeaks in struct LALFindStackSlidePeakParams */
/* 03/02/04 gam; if (params->thresholdFlag & 8 > 0) update pw_mean_thissum, pw_stddev_thissum ignoring peak bins; update pwr_snr. */
/* 03/02/04 gam; add returnPeaks to LALFindStackSlidePeakParams; if false return info about SUM and peaks without returning peaks themselves. */
/* 03/02/04 gam; Add LALFindStackSlidePeakOutputs to LALFindStackSlidePeaks to output pwMeanWithoutPeaks, pwStdDevWithoutPeaks, binsWithoutPeaks, acceptedEventCount, rejectedEventCount */
/* 03/03/04 gam; change windowFilterParam1 to normalizationThreshold */
/* 04/14/04 gam; if (params->normalizationFlag & 4) > 0 normalize STKs using running median */
/* 04/15/04 gam; Add INT2 params->debugOptionFlag */
/* 04/15/04 gam; Add debugOptionFlag to struct StackSlideSkyPatchParams */
/* 04/26/04 gam; Change LALStackSlide to StackSlide and LALSTACKSLIDE to STACKSLIDE for initial entry to LALapps. */
/* 05/05/04 gam; Change params->normalizationThreshold to params->normalizationParameter.  If normalizing with running median use this to correct bias in median to get mean. */
/* 05/26/04 gam; Change finishedSUMs to finishedSUMs; add startSUMs; defaults are TRUE; use to control I/O during Monte Carlo */
/* 05/26/04 gam; Add whichMCSUM = which Monte Carlo SUM; default is -1. */
/* 07/09/04 gam; If using running median, use LALRngMedBias to set params->normalizationParameter to correct bias in the median. */
/* 08/30/04 gam; if (outputEventFlag & 4) > 0 set returnOneEventPerSUM to TRUE; only the loudest event from each SUM is then returned. */
/* 10/28/04 gam; if (params->weightFlag & 1) > 0 then use PowerFlux weights from running median. Must have (params->normalizationFlag & 4) > 0 */
/* 10/28/04 gam; change unused params->windowFilterFlag to REAL8 params->orientationAngle used to find F_+ and F_x with weightFlag or MC with fixed polarization angle */
/* 11/01/04 gam; if (params->weightFlag & 8) > 0 rescale STKs with threshold5 to prevent dynamic range issues. */
/* 12/06/04 gam; get params->sampleRate, = effective sample rate, from the SFTs; calculate params->deltaT after reading SFTs. */
/* 12/06/04 gam; add params->gpsEpochStartTimeNan; get gpsEpochStartTime, gpsEpochStartTimeNan, and gpsStartTime from command line; */
/* 12/06/04 gam; change calibrationFlag to cosInclinationAngle */
  
#ifndef _DRIVESTACKSLIDE_H
#define _DRIVESTACKSLIDE_H

/*********************************************/
/*                                           */
/* START SECTION: include header files       */
/*                                           */
/*********************************************/
#include <stdio.h>
#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/LALConstants.h>
#include <lal/LALDemod.h>
#include <lal/Units.h>
#include <lal/LALInitBarycenter.h>
/* 02/04/04 gam; next two are for xml I/O */
#include <lal/LIGOLwXML.h>
#include <lal/LIGOLwXMLHeaders.h>
/* 04/14/04 gam; next includes LAL running median code */ /* 07/09/04 gam; include RngMedBias.h */
#include <lal/LALRunningMedian.h>
#include <lal/RngMedBias.h>
/* 02/09/04 gam; next is needed for tables defined in LAL */
#include <lal/LIGOMetadataTables.h>
#include "StackSlide.h"
/* #include <lal/LALStackSlide.h> Will need to switch to this version when StackSlide is in LAL. */
/*********************************************/
/*                                           */
/* END SECTION: include header files         */
/*                                           */
/*********************************************/

NRCSID( DRIVESTACKSLIDEH, "$Id$");

/*********************************************/
/*                                           */
/* START SECTION: define constants           */
/*                                           */
/*********************************************/
#define DRIVESTACKSLIDEH_ENULL 1
#define DRIVESTACKSLIDEH_EGPSTINT 2
#define DRIVESTACKSLIDEH_EDELTAT 3
#define DRIVESTACKSLIDEH_ETBLK 5
#define DRIVESTACKSLIDEH_ERA 6
#define DRIVESTACKSLIDEH_EDEC 7
#define DRIVESTACKSLIDEH_EFREQ 8
#define DRIVESTACKSLIDEH_EFREQDERIV 9
#define DRIVESTACKSLIDEH_EALOC 10
#define DRIVESTACKSLIDEH_ENODATA 11
#define DRIVESTACKSLIDEH_ENDATA 12
#define DRIVESTACKSLIDEH_ETIMESTEP 13
#define DRIVESTACKSLIDEH_ESTARTTIME 14
#define DRIVESTACKSLIDEH_ESTOPTIME 15
#define DRIVESTACKSLIDEH_ENTBLK 17
#define DRIVESTACKSLIDEH_EMCOHBLKPERCALL 18
#define DRIVESTACKSLIDEH_EIFONICKNAME 19
#define DRIVESTACKSLIDEH_EIFO 20
#define DRIVESTACKSLIDEH_ETARGETNAME 21
#define DRIVESTACKSLIDEH_ECHANNEL 22
#define DRIVESTACKSLIDEH_EISTARTTIME 23
#define DRIVESTACKSLIDEH_EMISSINGBLKDATA 24
#define DRIVESTACKSLIDEH_ESTARTFREQ 25
#define DRIVESTACKSLIDEH_EFREQSTEPSIZE 26
#define DRIVESTACKSLIDEH_EPARAMSPACEFLAG 27
#define DRIVESTACKSLIDEH_EKEEPTHISNEVENTS 28
#define DRIVESTACKSLIDEH_EOUTPUTREQUEST 29
#define DRIVESTACKSLIDEH_ENORMPARAM 30

#define DRIVESTACKSLIDEH_MSGENULL           "Null pointer"
#define DRIVESTACKSLIDEH_MSGEGPSTINT        "Unexpected GPS time interval"
#define DRIVESTACKSLIDEH_MSGEDELTAT         "Invalid deltaT"
#define DRIVESTACKSLIDEH_MSGETBLK           "Invalid tBLK"
#define DRIVESTACKSLIDEH_MSGERA             "Invalid right ascension"
#define DRIVESTACKSLIDEH_MSGEDEC            "Invalid declination"
#define DRIVESTACKSLIDEH_MSGEFREQ           "Invalid frequency +/- 0.5*band (could be negative, outside LIGO band, or too high for sample rate)"
#define DRIVESTACKSLIDEH_MSGEFREQDERIV      "One of the frequency derivatives is possibly too large; frequency will evolve outside allowed band"
#define DRIVESTACKSLIDEH_MSGEALOC           "Memory allocation error"
#define DRIVESTACKSLIDEH_MSGENODATA         "No input data was found"
#define DRIVESTACKSLIDEH_MSGENDATA          "Invalid number of input data points"
#define DRIVESTACKSLIDEH_MSGETIMESTEP       "Incorrect input data time step"
#define DRIVESTACKSLIDEH_MSGESTARTTIME      "Incorrect input data start time"
#define DRIVESTACKSLIDEH_MSGESTOPTIME       "Incorrect input data stop time"
#define DRIVESTACKSLIDEH_MSGENTBLK          "nSamplesPerBLK and tBLK are inconsistent with deltaT"
#define DRIVESTACKSLIDEH_MSGEMCOHBLKPERCALL "Invalid numSUMsPerCall"
#define DRIVESTACKSLIDEH_MSGEIFONICKNAME    "Invalid or null ifoNickName"
#define DRIVESTACKSLIDEH_MSGEIFO            "Invalid or null IFO"
#define DRIVESTACKSLIDEH_MSGETARGETNAME     "Invalid or null Target Name"
#define DRIVESTACKSLIDEH_MSGECHANNEL        "Invalid or null tDomainChannel"
#define DRIVESTACKSLIDEH_MSGEISTARTTIME     "Requested GPS start time resulted in invalid index to input data."
#define DRIVESTACKSLIDEH_MSGEMISSINGBLKDATA "Some requested input BLK data is missing"
#define DRIVESTACKSLIDEH_MSGESTARTFREQ      "Input BLK start frequency does not agree with that requested"
#define DRIVESTACKSLIDEH_MSGEFREQSTEPSIZE   "Input BLK frequency step size does not agree with that expected"
#define DRIVESTACKSLIDEH_MSGEPARAMSPACEFLAG  "Value for parameterSpaceFlag is invalid or currently unsupported"
#define DRIVESTACKSLIDEH_MSGEKEEPTHISNEVENTS "2nd bit in outputEventFlag set to keep loudest but keepThisNumber was < 1!"
#define DRIVESTACKSLIDEH_MSGEOUTPUTREQUEST   "Cannot set thresholdFlag < 1 and outputEventFlag to output everything!"
#define DRIVESTACKSLIDEH_MSGENORMPARAM       "Cannot have normalizationParameter less than ln(2) or greater than 1 when using running median."

/* Limit on size of arrays holding channel names */
/* #define dbNameLimit 256; */ /* Should be defined in LAL? */

/*********************************************/
/*                                           */
/* END SECTION: define constants             */
/*                                           */
/*********************************************/

/*********************************************/
/*                                           */
/* START SECTION: define macros              */
/*                                           */
/*********************************************/

/* 02/04/04 gam; MACROS used for xml I/O */

/* 02/09/04 gam; add for StackSlideSFTsSearchSummaryTable */
/* 02/12/04 gam; Add num_sums = numSUMsTotal, freq_index = f0SUM*params->tEffSUM, and num_bins = nBinsPerSUM to StackSlideSFTsSearchSummaryTable */
#define LIGOLW_XML_LOCAL_SEARCHSUMMARY_STACKSLIDESFTS \
"   <Table Name=\"searchsummary_stackslidesftsgroup:searchsummary_stackslidesfts:table\">\n" \
"      <Column Name=\"searchsummary_stackslidesftsgroup:searchsummary_stackslidesfts:process_id\" Type=\"ilwd:char\"/>\n" \
"      <Column Name=\"searchsummary_stackslidesftsgroup:searchsummary_stackslidesfts:ifo\" Type=\"lstring\"/>\n" \
"      <Column Name=\"searchsummary_stackslidesftsgroup:searchsummary_stackslidesfts:data_directory\" Type=\"lstring\"/>\n" \
"      <Column Name=\"searchsummary_stackslidesftsgroup:searchsummary_stackslidesfts:comment\" Type=\"lstring\"/>\n" \
"      <Column Name=\"searchsummary_stackslidesftsgroup:searchsummary_stackslidesfts:start_time\" Type=\"int_4s\"/>\n" \
"      <Column Name=\"searchsummary_stackslidesftsgroup:searchsummary_stackslidesfts:start_time_ns\" Type=\"int_4s\"/>\n" \
"      <Column Name=\"searchsummary_stackslidesftsgroup:searchsummary_stackslidesfts:duration\" Type=\"real_8\"/>\n" \
"      <Column Name=\"searchsummary_stackslidesftsgroup:searchsummary_stackslidesfts:sft_baseline\" Type=\"real_8\"/>\n" \
"      <Column Name=\"searchsummary_stackslidesftsgroup:searchsummary_stackslidesfts:num_sfts\" Type=\"int_4s\"/>\n" \
"      <Column Name=\"searchsummary_stackslidesftsgroup:searchsummary_stackslidesfts:start_freq\" Type=\"real_8\"/>\n" \
"      <Column Name=\"searchsummary_stackslidesftsgroup:searchsummary_stackslidesfts:band\" Type=\"real_8\"/>\n" \
"      <Column Name=\"searchsummary_stackslidesftsgroup:searchsummary_stackslidesfts:sum_baseline\" Type=\"real_8\"/>\n" \
"      <Column Name=\"searchsummary_stackslidesftsgroup:searchsummary_stackslidesfts:freq_index\" Type=\"int_4s\"/>\n" \
"      <Column Name=\"searchsummary_stackslidesftsgroup:searchsummary_stackslidesfts:num_bins\" Type=\"int_4s\"/>\n" \
"      <Column Name=\"searchsummary_stackslidesftsgroup:searchsummary_stackslidesfts:num_sums\" Type=\"int_4s\"/>\n" \
"      <Stream Name=\"searchsummary_stackslidesftsgroup:searchsummary_stackslidesfts:table\" Type=\"Local\" Delimiter=\",\">\n"

#define LOCAL_SEARCHSUMMARY_STACKSLIDESFTS_ROW \
"         \"process:process_id:0\",\"%s\",\"%s\",\"%s\",%d,%d,%22.16e,%22.16e,%d,%22.16e,%22.16e,%22.16e,%d,%d,%d"

/* 02/09/04 gam; remove fderiv_2-5; remove pw_max_thissum and freq_max_thissum; change false_alarm_prob_upperlimit to false_alarm_prob */
/* 02/11/04 gam; Use power rather than amplitudes when finding peaks; change amplitude to power and snr to pwr_snr */
/* 02/12/04 gam; Add freq_index to SnglStackSlidePeriodicTable */
#define LIGOLW_XML_SNGL_LOCAL_STACKSLIDEPERIODIC \
"   <Table Name=\"sngl_stackslideperiodicgroup:sngl_stackslideperiodic:table\">\n" \
"      <Column Name=\"sngl_stackslideperiodicgroup:sngl_stackslideperiodic:process_id\" Type=\"ilwd:char\"/>\n" \
"      <Column Name=\"sngl_stackslideperiodicgroup:sngl_stackslideperiodic:sum_no\" Type=\"int_4s\"/>\n" \
"      <Column Name=\"sngl_stackslideperiodicgroup:sngl_stackslideperiodic:event_no_thissum\" Type=\"int_4s\"/>\n" \
"      <Column Name=\"sngl_stackslideperiodicgroup:sngl_stackslideperiodic:overlap_event\" Type=\"int_4s\"/>\n" \
"      <Column Name=\"sngl_stackslideperiodicgroup:sngl_stackslideperiodic:frequency\" Type=\"real_8\"/>\n" \
"      <Column Name=\"sngl_stackslideperiodicgroup:sngl_stackslideperiodic:freq_index\" Type=\"int_4s\"/>\n" \
"      <Column Name=\"sngl_stackslideperiodicgroup:sngl_stackslideperiodic:power\" Type=\"real_4\"/>\n" \
"      <Column Name=\"sngl_stackslideperiodicgroup:sngl_stackslideperiodic:width_bins\" Type=\"int_4s\"/>\n" \
"      <Column Name=\"sngl_stackslideperiodicgroup:sngl_stackslideperiodic:num_subpeaks\" Type=\"int_4s\"/>\n" \
"      <Column Name=\"sngl_stackslideperiodicgroup:sngl_stackslideperiodic:pwr_snr\" Type=\"real_4\"/>\n" \
"      <Column Name=\"sngl_stackslideperiodicgroup:sngl_stackslideperiodic:false_alarm_prob\" Type=\"real_4\"/>\n" \
"      <Column Name=\"sngl_stackslideperiodicgroup:sngl_stackslideperiodic:goodness_of_fit\" Type=\"real_4\"/>\n" \
"      <Column Name=\"sngl_stackslideperiodicgroup:sngl_stackslideperiodic:sky_ra\" Type=\"real_8\"/>\n" \
"      <Column Name=\"sngl_stackslideperiodicgroup:sngl_stackslideperiodic:sky_dec\" Type=\"real_8\"/>\n" \
"      <Column Name=\"sngl_stackslideperiodicgroup:sngl_stackslideperiodic:fderiv_1\" Type=\"real_8\"/>\n" \
"      <Column Name=\"sngl_stackslideperiodicgroup:sngl_stackslideperiodic:pw_mean_thissum\" Type=\"real_4\"/>\n" \
"      <Column Name=\"sngl_stackslideperiodicgroup:sngl_stackslideperiodic:pw_stddev_thissum\" Type=\"real_4\"/>\n" \
"      <Stream Name=\"sngl_stackslideperiodicgroup:sngl_stackslideperiodic:table\" Type=\"Local\" Delimiter=\",\">\n"

/* 02/09/04 gam; remove fderiv_2-5; remove pw_max_thissum and freq_max_thissum; change false_alarm_prob_upperlimit to false_alarm_prob */
#define SNGL_LOCAL_STACKSLIDEPERIODIC_ROW \
"         \"process:process_id:0\",%d,%d,%d,%22.16e,%d,%e,%d,%d,%e,%e,%e,%22.16e,%23.16e,%23.16e,%e,%e"

/*********************************************/
/*                                           */
/* END SECTION: define macros                */
/*                                           */
/*********************************************/

/*********************************************/
/*                                           */
/* START SECTION: define structs             */
/*                                           */
/*********************************************/
/* Define structure to hold statistics about a power series */ /* 02/11/04 gam; comment out */
/* typedef struct tagFreqSeriesPowerStats {
	REAL8 pwMax;
	REAL8 freqPWMax;
	REAL8 pwMean;
	REAL8 pwStdDev;
}
FreqSeriesPowerStats; */

/* 02/04/04 gam; Move items here that need to recorded only once per job */ /* 02/09/04 gam; summarize for stack slide of SFTs */
/* 02/12/04 gam; Add num_sums = numSUMsTotal, freq_index =f0SUM*params->tEffSUM, and num_bins = nBinsPerSUM to StackSlideSFTsSearchSummaryTable */
typedef struct
tagStackSlideSFTsSearchSummaryTable
{
  struct tagStackSlideSFTsSearchSummaryTable *next;
  CHAR          ifo[3];
  CHAR          data_directory[256];
  CHAR          comment[256];
  UINT4         start_time;
  UINT4         start_time_ns;
  REAL8         duration;
  REAL8         sft_baseline;
  INT4          num_sfts;
  REAL8         start_freq;
  REAL8         band;
  REAL8         sum_baseline;
  INT4          freq_index;
  INT4          num_bins;
  INT4          num_sums;
}
StackSlideSFTsSearchSummaryTable;

/* 01/14/04 gam; Structure that holds data for significant peak found in SUMs */ /* 02/04/04 gam Reduce this table to essentials */
/* 02/09/04 gam; remove fderiv_2-5; remove pw_max_thissum and freq_max_thissum; change false_alarm_prob_upperlimit to false_alarm_prob */
/* 02/11/04 gam; Use power rather than amplitudes when finding peaks; change amplitude to power and snr to pwr_snr */
/* 02/12/04 gam; Add freq_index to SnglStackSlidePeriodicTable */
typedef struct
tagSnglStackSlidePeriodicTable
{
  struct tagSnglStackSlidePeriodicTable *next;
  INT4          sum_no;
  INT4          event_no_thissum;
  INT4          overlap_event;
  REAL8         frequency;
  INT4          freq_index;
  REAL4         power;
  INT4          width_bins;
  INT4          num_subpeaks;
  REAL4         pwr_snr;
  REAL4         false_alarm_prob;
  REAL4         goodness_of_fit;
  REAL8         sky_ra;
  REAL8         sky_dec;
  REAL8         fderiv_1;
  REAL4         pw_mean_thissum;
  REAL4         pw_stddev_thissum;
}
SnglStackSlidePeriodicTable;

/* 02/04/04 gam; comment out and just use SnglStackSlidePeriodicTable for peaks */
/* typedef struct
tagStackSlidePeak
{
  struct tagStackSlidePeak *next;
  REAL8         frequency;
  REAL4         amplitude;
  REAL4         width;
  REAL4         snr;
  REAL4         false_alarm_prob;
  REAL4         goodness_of_fit;
}
StackSlidePeak; */

/* 03/02/04 gam; if (params->thresholdFlag & 8 > 0) set returnMeanStdWithoutPeaks */
typedef struct
tagLALFindStackSlidePeakOutputs
{
  SnglStackSlidePeriodicTable **pntrToPeaksPntr;
  REAL4   *pwMeanWithoutPeaks;
  REAL4   *pwStdDevWithoutPeaks;
  INT4    *binsWithoutPeaks;
  INT4    *acceptedEventCount;
  INT4    *rejectedEventCount;
}
LALFindStackSlidePeakOutputs;

/* 01/20/04 gam; Change findStackSlidePeaks to LALFindStackSlidePeaks; put params into struct */
/* 02/20/04 gam; change threshold6 to INT4 maxWidthBin; remove maxWidth; replace fMax and fMin with ifMax and ifMin */
/* 03/02/04 gam; if (params->thresholdFlag & 8 > 0) set updateMeanStdDev true */
/* 03/02/04 gam; if returnPeaks return peaks; if false returns just pwMeanWithoutPeaks and pwStdDevWithoutPeaks. */
typedef struct
tagLALFindStackSlidePeakParams
{
  BOOLEAN returnPeaks;        /* 03/02/04 gam */
  BOOLEAN updateMeanStdDev;   /* 03/02/04 gam */  
  BOOLEAN vetoWidePeaks;      /* 02/20/04 gam */
  BOOLEAN vetoOverlapPeaks;   /* 02/20/04 gam */
  BOOLEAN returnOneEventPerSUM; /* 08/30/04 gam */
  REAL4   threshold1;
  REAL4   threshold2;
  INT4    maxWidthBins;   /* 02/17/04 gam; maximum width in bins */
  INT4    ifMin;          /* 02/20/04 gam; index of minimum freq not in overlap region */
  INT4    ifMax;          /* 02/20/04 gam; index of maximum freq not in overlap region */
  REAL4   fracPwDrop;     /* 01/27/04 gam; input as threshold3 */
  INT4    iSUM;           /* 02/04/04 gam; index of the Sum used for bookkeeping */
  INT4    iSky;           /* 02/04/04 gam; index into skyPosData to find sky position */
  INT4    iDeriv;         /* 02/04/04 gam; index into freqDerivData to find freq derivatives */
  REAL8 **skyPosData;     /* 02/04/04 gam; skyPosData[iSky][0] = RA, skyPosData[iSky][1] = DEC */
  REAL8 **freqDerivData;  /* 02/04/04 gam; freqDerivData[iDeriv][i] = freq deriv i */
  INT4    numSpinDown;    /* 02/04/04 gam; number of spindown */
}
LALFindStackSlidePeakParams;

/* 02/20/04 gam; change threshold6 to INT4 maxWidthBin; remove maxWidth; replace fMax and fMin with ifMax and ifMin */
typedef struct
tagLALUpdateLoudestStackSlideParams
{
  INT4    arraySize;           /* 02/17/04 gam; size of the loudestPeaksArray */
  INT4    nBinsPerOutputEvent; /* 02/17/04 gam; each loudest event corresponds to this number of event */
  INT4    ifMin;               /* 02/20/04 gam; index of minimum freq not in overlap region */
  INT4    ifMax;               /* 02/20/04 gam; index of maximum freq not in overlap region */
  INT4    iSUM;                /* 02/17/04 gam; index of the Sum used for bookkeeping */
  INT4    iSky;                /* 02/17/04 gam; index into skyPosData to find sky position */
  INT4    iDeriv;              /* 02/17/04 gam; index into freqDerivData to find freq derivatives */
  REAL8 **skyPosData;          /* 02/17/04 gam; skyPosData[iSky][0] = RA, skyPosData[iSky][1] = DEC */
  REAL8 **freqDerivData;       /* 02/17/04 gam; freqDerivData[iDeriv][i] = freq deriv i */
  INT4    numSpinDown;         /* 02/17/04 gam; number of spindown */
}
LALUpdateLoudestStackSlideParams;

/* 01/31/04 gam; struct that holds info on patch of the sky code is working on. */
typedef struct
tagStackSlideSkyPatchParams
{
  REAL8   startRA;                   /* Starting Right Ascension in radians */
  REAL8   stopRA;                    /* Ending Right Ascension in radians */
  REAL8   deltaRA;                   /* Right Ascension step size in radians */
  INT4    numRA;                     /* Number of RA to compute SUMs for */

  REAL8   startDec;                  /* Starting Declination in radians */
  REAL8   stopDec;                   /* Ending Declination in radians  */
  REAL8   deltaDec;                  /* Declination step size in radians */
  INT4    numDec;                    /* Number of Dec to compute SUMs for */
  
  INT2    debugOptionFlag;           /* 04/15/04 gam; for debugging */
}
StackSlideSkyPatchParams;

typedef struct tagStackSlideSearchParams {

  /******************************************/
  /*                                        */
  /* START SECTION: parameters passed as    */
  /* arguments.  Those indented are         */
  /* computed from these.                   */
  /*                                        */
  /******************************************/
  
  /* 12/06/04 gam; add params->gpsEpochStartTimeNan; get gpsEpochStartTime, gpsEpochStartTimeNan, and gpsStartTime from command line; */
  UINT4   gpsEpochStartTimeSec;      /* GPS start time of data requested seconds */
  UINT4   gpsEpochStartTimeNan;      /* GPS start time of data requested nanoseconds */
  UINT4   gpsStartTimeSec;           /* GPS start time of data requested seconds */
    UINT4   gpsStartTimeNan;         /* GPS start time of data requested nanoseconds; currently fixed as zero. */
  REAL8   duration;                  /* Total time being analyzed  */

  INT4    numBLKs;                   /* Number of input BLKs.  Not duration/tBLK if gaps are present */
  REAL8   tBLK;                      /* duration in seconds of BLKs */
  REAL8   tEffBLK;                   /* Effective duration of BLK such that dfBLK = 1.0/tEffBLK  */
    REAL8   dfBLK;                     /* Freq resolution of BLK = 1.0/tEffBLK  Could be different from 1/tBLK due to oversampling  */
  REAL8   f0BLK;                     /* Start frequency of the input BLKs */
  REAL8   bandBLK;                   /* Band width of the input BLKs. */
  INT4    nBinsPerBLK;               /* Number of data points in the input BLK in the input frequency band */

  INT4    numBLKsPerSTK;             /* Number of BLKs to use to make one STK */
    INT4    numSTKs;                   /* Number of actual STKs. Not duration/tSTK if gaps are present */
    REAL8   tSTK;                      /* duration in seconds of STKs = tBLK*numBLKsPerSTK */
  REAL8   tEffSTK;                   /* Effective duration of STK such that dfSTK = 1.0/tEffSTK */
    REAL8   dfSTK;                     /* Freq resolution of STK = 1.0/tEffSTK. Could be different from 1/tSTK due to oversampling  */
  REAL8   f0STK;                     /* Start frequency of STKs */
  REAL8   bandSTK;                   /* Band width of the STKs. */
  INT4    nBinsPerSTK;               /* Number of data points in the STKs in the STK frequency band */

  INT4    numSTKsPerSUM;             /* Number of STKs to use to make one SUM (Usually duration/tSTK) */
    INT4    numSUMsPerParamSpacePt;  /* Number of output SUMs per parameter space point = params->duration/params->tSUM. (Usually will = 1) */
    REAL8   tSUM;                      /* duration in seconds of output SUMs = tSTK*numSTKsPerSUM. (Usually = duration) */
  REAL8   tEffSUM;                   /* Effective duration of SUM such that dfSUM = 1.0/tEffSUM */
    REAL8   dfSUM;                     /* Freq resolution of output SUMs.  dfSUM = 1.0/tEffSUM. Could be different from dfSTK in some algorithms?  */
  REAL8   f0SUM;                     /* Start frequency of output SUMs */
  REAL8   bandSUM;                   /* Band width of the output SUMs. */
  INT4    nBinsPerSUM;               /* Number of data points in the output SUMs in the input frequency band */

  CHAR    *ifoNickName;              /* 2 character ifoNickName = H1, H2, or L1; note IFO is the site = LHO, LLO, GEO, etc.... */
  CHAR    *IFO;                      /* Identifies the interferometer site = LHO, LLO, GEO, etc.... */
  CHAR    *patchName;                /* Name of the patch in parameter space of the search (e.g, Galactic Center) */
  INT4    patchNumber;               /* Index of the patch in parameter space of the search (will be -1 until indexing implemented) */

  /* CHAR    tDomainChannel[dbNameLimit]; */ /* Time domain channel name used to produce blocks */
  /* CHAR    fDomainChannel[dbNameLimit]; */ /* Frequency domain channel name that blocks are stored as */
  
   CHAR    tDomainChannel[256]; /* Time domain channel name used to produce blocks */
   CHAR    fDomainChannel[256]; /* Frequency domain channel name that blocks are stored as */
  
  INT2    inputDataTypeFlag;         /* 0 = SFTs, 1 = PSDs, 2 = F-stat */

  INT2    parameterSpaceFlag;        /* 0 = use input delta for each param, 1 = params are input as vectors of data, 2 = use input param metric, 3 = create param metric */
    INT4    numParamSpacePts;          /* Total number of points in the parameter space to cover */
    INT4    numSUMsTotal;              /* Total Number of Sums output = numSUMsPerParamPt*numParamSpacePts */

  INT2    stackTypeFlag;             /* 0 means stacks are PSDs from SFTs, 1 means stacks are F-statistic from SFTs */
  INT4    Dterms;                    /* = terms used in the computation of the dirichlet kernel in LALDemod (Default value is 32) */
  
  /* 01/14/04 gam; Change threshold1 and threshold2 from type REAL8 to type REAL4 */  
  /* 02/20/04 gam; reorganize */
  INT2    thresholdFlag;             /* How to apply the thresholds */
  REAL4   threshold1;                /* if thresholdFlag > 0 then peak found if power above this */
  REAL4   threshold2;                /* if thresholdFlag > 0 then peak ends if power drops below this */
  REAL4   threshold3;                /* if thresholdFlag > 0 then this is fracPwDrop that indicates a new peak. */
  REAL4   threshold4;                /* unused */
  REAL4   threshold5;                /* unused */
  INT4    maxWidthBins;              /* maximum width in bins */

  /* INT2    calibrationFlag; */     /* 12/06/04 gam */ /* Flag that specifies what calibration to do; -1 means Blks are already calibrated, 0 means leave uncalibrated; 1 calibrate  */

  INT2    weightFlag;                /* Flag that specifies whether to weight BLKs or STKs with a(t) or b(t).  */

  REAL8   orientationAngle;          /* 10/28/04 gam; change unused params->windowFilterFlag to REAL8 params->orientationAngle used to find F_+ and F_x with weightFlag or MC with fixed polarization angle */
  REAL8   cosInclinationAngle;       /* 12/06/04 gam */
          
  /* INT2    windowFilterFlag;       */   /* 10/28/04 gam */ /* Flag that specifies whether any filtering or windowing was done or should be done. (If < 0 then specifies what was done in the time domain) */
  /* REAL8   windowFilterParam1;     */   /* 03/03/04 gam */ /* 1st parameter to use in windowing or filtering */
  /* REAL8   windowFilterParam2;     */   /* 02/17/04 gam */ /* 2nd paramter to use in windowing or filtering */
  /* REAL8   windowFilterParam3;     */   /* 02/17/04 gam */ /* 3rd paramter to use in windowing or filtering */

  INT2    normalizationFlag;         /* Flag that specifies what normalization to do.  If < 0 then specifies what normalization was already done. */
  REAL8   f0NRM;                     /* Start frequency to normalize over to create NRMs */
  REAL8   bandNRM;                   /* Band width to normalize over to create NRMs. */
  INT4    nBinsPerNRM;               /* Number of data points to normalize over to create to creat NRMs */
  REAL4   normalizationParameter;    /* 03/03/04 gam; change windowFilterParam1 to normalizationThreshold */ /* 05/05/04 gam; Change normalizationThreshold to normalizationParameter. */

  INT2    testFlag;                  /* Specifies any tests or debugging to do; 0 means no test */

  INT4   numSUMsPerCall;            /* If > 0 then = number of SUMs (parameter space points) to compute during each call to ApplySearch */
  
  /* INT2    outputFlag; */          /* Flag that specifies what to output; if > 0 then will output Sums into frame file */
  INT2   outputSUMFlag;             /* 02/17/04 gam; Flag that determines whether to output SUMs e.g., in ascii. */
  INT2   outputEventFlag;           /* 02/17/04 gam; Flag that deterines xml output */
  INT4   keepThisNumber;            /* 02/17/04 gam; How many event events to keep; keep loudest in each bandSUM/keepThisNumber band */

  /* 01/27/04 gam; Change unusedFlag1, unusedParam1 and unusedParam2 to threshold4, threshold5, and threshold6. */
  /*  INT2    unusedFlag1;             */  /* Place holder for future flag */
  /*  REAL8   unusedParam1;            */  /* Place holder for future parameter */
  /* 01/20/04 gam; Change unusedFlag2 to threshold3. */
  /* INT2    unusedFlag2;            */  /* Place holder for future flag */
  /*  REAL8   unusedParam2;            */ /* Place holder for future parameter */

  StackSlideSkyPatchParams *stksldSkyPatchData;  /* 01/31/04 gam; struct that holds info on patch of the sky code is working on. */
  /* REAL8   startRA;                 */  /* Starting Right Ascension in radians */
  /* REAL8   stopRA;                  */  /* Ending Right Ascension in radians */
  /* REAL8   deltaRA;                 */  /* Right Ascension step size in radians */
  /* INT4    numRA;                   */  /* Number of RA to compute SUMs for */

  /* REAL8   startDec;                */  /* Starting Declination in radians */
  /* REAL8   stopDec;                 */  /* Ending Declination in radians  */
  /* REAL8   deltaDec;                */  /* Declination step size in radians */
  /* INT4    numDec;                  */  /* Number of Dec to compute SUMs for */
  
  INT4    numSpinDown;               /* Number of nonzero spindown parameters */

  REAL8   startFDeriv1;              /* Starting 1st deriv of freq in Hz/s */
  REAL8   stopFDeriv1;               /* Ending 1st deriv of freq in Hz/s */
  REAL8   deltaFDeriv1;              /* 1st deriv of freq step size in Hz/s */
  INT4    numFDeriv1;                /* Number of 1st derivs to compute SUMs for */

  REAL8   startFDeriv2;              /* Starting 2nd deriv of freq in Hz/s^2 */
  REAL8   stopFDeriv2;               /* Ending 2nd deriv of freq in Hz/s^2 */
  REAL8   deltaFDeriv2;              /* 2nd deriv of freq step size in Hz/s^2 */
  INT4    numFDeriv2;                /* Number of 2nd derivs to compute SUMs for */

  REAL8   startFDeriv3;              /* Starting 3rd deriv of freq in Hz/s^3 */
  REAL8   stopFDeriv3;               /* Ending 3rd deriv of freq in Hz/s^3 */
  REAL8   deltaFDeriv3;              /* 3rd deriv of freq step size in Hz/s^3 */
  INT4    numFDeriv3;                /* Number of 3rd derivs to compute SUMs for */

  REAL8   startFDeriv4;              /* Starting 4th deriv of freq in Hz/s^4 */
  REAL8   stopFDeriv4;               /* Ending 4th deriv of freq in Hz/s^4 */
  REAL8   deltaFDeriv4;              /* 4th deriv of freq step size in Hz/s^4 */
  INT4    numFDeriv4;                /* Number of 4th derivs to compute SUMs for */

  REAL8   startFDeriv5;              /* Starting 5th deriv of freq in Hz/s^5 */
  REAL8   stopFDeriv5;               /* Ending 5th deriv of freq in Hz/s^5 */
  REAL8   deltaFDeriv5;              /* 5th deriv of freq step size in Hz/s^5 */
  INT4    numFDeriv5;                /* Number of 5th derivs to compute SUMs for */
  
  CHAR    *sunEdatFile;              /* File with sun ephemeris data */  
  CHAR    *earthEdatFile;            /* File with earth ephemeris data */
  CHAR    *sftDirectory;             /* Directory with SFTs */    
  CHAR    *outputFile;               /* Name of file to output results in */  

  INT2    debugOptionFlag;           /* 04/15/04 gam; Add INT2 params->debugOptionFlag */
  /******************************************/
  /*                                        */
  /* END SECTION: parameters passed as      */
  /* arguments.  Those indented are         */
  /* computed from these.                   */
  /*                                        */
  /******************************************/

  /******************************************/
  /*                                        */
  /* START SECTION: other parameters        */
  /*                                        */
  /******************************************/
  
  /* 12/06/04 gam; get params->sampleRate, = effective sample rate, from the SFTs; calculate params->deltaT after reading SFTs. */
  REAL8   sampleRate;                /* Sample rate of the time-domain data used to make the frequency-domain input data Blocks */
  REAL8   deltaT;                    /* Time step size in seconds = 1.0/sampleRate */
  
  CHAR  *dsoName;         /* Name of this DSO */ /* 11/05/01 gam */

  /* Basic beowulf node descriptors */
  BOOLEAN searchMaster;           /* TRUE on the search master */
  UINT4   rank;                   /* Rank of this slave */
  UINT4   numSlaves;              /* Total number of slaves */
  UINT4   curSlaves;              /* Current number of slaves still working */
  UINT4   numNodes;               /*  07/12/02 gam  Should equal number of slaves + 1 for search master */

  /* UINT4  numDBOutput; */ /* 02/11/04 gam */  /* Used to keep track of the number of rows written to the database */

  /* gpsTimeInterval times;       */ /* The GPS start and end time of the actual data */
  
  /* 02/09/04 gam; use to hold times for search summary table */
  LIGOTimeGPS  actualStartTime; 
  LIGOTimeGPS  actualEndTime;
      
  LALUnit     unitBLK;            /* Unit BLKs are stored in (e.g., strain per root Hz) */

  FFT **BLKData;                   /* Container for BLKs */
  REAL4FrequencySeries **STKData;  /* Container for STKs */
  REAL4FrequencySeries **SUMData;  /* Container for SUMs */

  REAL4Vector **savSTKData;              /* 10/28/04 gam; save STKs for reuse with powerFlux style weighting of STKs for each sky position */
  REAL4Vector *detResponseTStampMidPts;  /* 10/28/04 gam; container for detector response F_+ or F_x for one sky position, one polarization angle, for midpoints of a timeStamps */
  REAL4Vector **inverseSquareMedians;    /* 10/28/04 gam; container with inverse square medians for each STK for each frequency bin; for use with powerFlux style weighting of STKs */
  REAL4Vector *sumInverseSquareMedians;  /* 10/28/04 gam; container with sum of inverse square medians for each frequency bin; for use with powerFlux style weighting of STKs. */
  
  INT4 iMinBLK;      /* Index of minimum frequency in BLK band */
  INT4 iMaxBLK;      /* Index of maximum frequency in BLK band */
  INT4 iMinNRM;      /* Index of mimimum frequency to include when normalizing BLKs */
  INT4 iMaxNRM;      /* Index of maximum frequency to include when normalizing BLKs */

  LIGOTimeGPS *timeStamps;          /* Container for time stamps giving initial time for which each BLK was computed. */
  /* gpsTimeInterval *timeIntervals; */  /* Array of time intervals; needed for building frame output */

  EphemerisData *edat; /* 07/10/02 gam Add EphemerisData *edat = pointer to ephemeris data to StackSlideSearchParams */

  BOOLEAN finishedBLKs;   /* Set equal to true when all BLKS for this job have been found in input data */
  BOOLEAN finishedSTKs;   /* Set equal to true when all STKS for this job have been created */
  /* BOOLEAN finishedSUMs; */ /* 05/26/04 gam */
  BOOLEAN startSUMs;    /* 05/26/04 gam; use to control I/O during Monte Carlo  */
  BOOLEAN finishSUMs;   /* 05/26/04 gam; use to control I/O during Monte Carlo  */

  INT4 whichMCSUM;      /* 05/26/04 gam; which SUM the Monte Carlo Simulation is running on. */

  INT4 whichSTK;          /* which STK does BLK go with  */
  INT4 lastWhichSTK;      /* Last value of whichSTK does BLK go with  */

  REAL8 **skyPosData;     /* Container for Parameter Space Data */
  REAL8 **freqDerivData;  /* Container for Frequency Derivative Data */
  INT4 numSkyPosTotal;    /* Total number of Sky positions to cover */
  INT4 numFreqDerivTotal; /* Total number of Frequency evolution models to cover */

  /* FreqSeriesPowerStats *SUMStats; */ /* 02/11/04 gam */ /* Container for statistics about each SUM */
  
  /* SnglStackSlidePeriodicTable *significantPeaks; 02/04/04 gam */ /* 01/14/04 gam; Structure that holds data for significant peak found in SUMs */
  /* StackSlidePeak *peaks; 02/04/04 gam */ /* 01/14/04 gam; Structure that holds data for peaks found in SUMs */
  SnglStackSlidePeriodicTable  *peaks; /* 02/04/04 gam; Use this structure to hold link list of data for peaks found in SUMs */
  
  LIGOLwXMLStream *xmlStream;  /* 02/09/04 gam; for xml output */
  
  /******************************************/
  /*                                        */
  /* END SECTION: other parameters          */
  /*                                        */
  /******************************************/

}
StackSlideSearchParams;

/*********************************************/
/*                                           */
/* END SECTION: define structs               */
/*                                           */
/*********************************************/

/******************************************/
/*                                        */
/* START SECTION: prototype declarations  */
/*                                        */
/******************************************/
void StackSlideInitSearch(
    LALStatus              *status,
    StackSlideSearchParams *params,
    int                     argc,
    char                   *argv[]    
    );     
void StackSlideConditionData(
    LALStatus              *status,
    StackSlideSearchParams *params
    );
void StackSlideApplySearch(
    LALStatus              *status,
    StackSlideSearchParams *params
    );
void StackSlideFinalizeSearch(
    LALStatus              *status,
    StackSlideSearchParams *params
    );  
/* void ComputePowerStats( const REAL4Vector *fsVec, FreqSeriesPowerStats *stats, REAL8 fStart, REAL8 df ); */ /* 02/11/04 gam */
/* void InitializePeaksArray( SnglStackSlidePeriodicTable *loudestPeaksArray, INT4 arraySize ); */ /* 02/20/04 gam */ /* 02/17/04 gam; this and next two lines */
void InitializePeaksArray( SnglStackSlidePeriodicTable *loudestPeaksArray, INT4 arraySize, REAL8 f0, REAL8 df, INT4 nBinsPerOutputEvent );
void LALUpdateLoudestFromPeaks( SnglStackSlidePeriodicTable *loudestPeaksArray, const SnglStackSlidePeriodicTable *peaks, INT4 nBinsPerOutputEvent );
void LALUpdateLoudestFromSUMs( SnglStackSlidePeriodicTable *loudestPeaksArray, const REAL4FrequencySeries *oneSUM, const LALUpdateLoudestStackSlideParams *params );
/* 01/14/04 gam; Function for finding peaks in StackSlide SUMs */
/* void findStackSlidePeaks( const REAL4FrequencySeries *oneSUM, StackSlidePeak *peaks, REAL4 threshold1, REAL4 threshold2 ) */
/* 01/20/04 gam; Change findStackSlidePeaks to LALFindStackSlidePeaks; put params into struct */
/* void LALFindStackSlidePeaks( StackSlidePeak *peaks, const REAL4FrequencySeries *oneSUM, LALFindStackSlidePeakParams *params ); */ /* 02/04/04 gam */
/* void LALFindStackSlidePeaks( SnglStackSlidePeriodicTable *peaks, const REAL4FrequencySeries *oneSUM, LALFindStackSlidePeakParams *params ); */ /* 02/04/04 gam */
/* void LALFindStackSlidePeaks( SnglStackSlidePeriodicTable **pntrToPeaksPntr, const REAL4FrequencySeries *oneSUM, LALFindStackSlidePeakParams *params ); */ /* 03/02/04 gam */
void LALFindStackSlidePeaks( LALFindStackSlidePeakOutputs *outputs, const REAL4FrequencySeries *oneSUM, LALFindStackSlidePeakParams *params );
void CountOrAssignSkyPosData(REAL8 **skyPosData, INT4 *numSkyPosTotal, BOOLEAN returnData, StackSlideSkyPatchParams *params); /* 01/31/04 gam*/
/* 10/28/04 gam; function that saves inverse square medians for weighting STKs */                 
void SaveInverseSquareMedians(REAL4Vector **inverseSquareMedians, REAL4Vector *sumInverseSquareMedians, INT2 weightFlag, REAL4Vector *medians, INT4 k, INT4 mediansOffset1, INT4 mediansOffset2, INT4 mediansLengthm1, INT4 nBinsPerSTK);
/* 10/28/04 gam; savSTKDATA for reuse with powerFlux style weighting of STKs for each sky position */
void SaveSTKData(REAL4Vector ** savSTKData, REAL4FrequencySeries **STKData, INT4 numSTKs, INT4 nBinsPerSTK);
/* 10/28/04 gam; apply powerFlux style weights */
void WeightSTKsWithoutBeamPattern(REAL4FrequencySeries **STKData, REAL4Vector **inverseSquareMedians, REAL4Vector *sumInverseSquareMedians, INT4 numSTKs, INT4 nBinsPerSTK);
/* 10/28/04 gam; get squared detector response F_+^2 or F_x^2 for one sky position, one polarization angle, for midpoints of a timeStamps */
void GetDetResponseTStampMidPts(LALStatus *stat, REAL4Vector *detResponseTStampMidPts, LIGOTimeGPS *timeStamps, INT4 numSTKs, REAL8 tSTK,
     LALDetector *cachedDetector, REAL8 *skyPosData, REAL8 orientationAngle, CoordinateSystem coordSystem, INT2 plusOrCross);
/* 10/28/04 gam; apply powerFlux style weights including detector beam pattern response */
void WeightSTKsIncludingBeamPattern(REAL4FrequencySeries **STKData,
           REAL4Vector ** savSTKData,
           REAL4Vector **inverseSquareMedians,
           REAL4Vector *sumInverseSquareMedians,
           REAL4Vector *detResponseTStampMidPts,
           INT4 numSTKs, INT4 nBinsPerSTK, REAL8 tSTK);
/* 11/01/04 gam; if (params->weightFlag & 8) > 0 rescale STKs with threshold5 to prevent dynamic range issues. */
void RescaleSTKData(REAL4FrequencySeries **STKData, INT4 numSTKs, INT4 nBinsPerSTK,REAL4 RescaleFactor);
/******************************************/
/*                                        */
/* END SECTION: prototype declarations    */
/*                                        */
/******************************************/

#endif /* _DRIVESTACKSLIDE_H */
