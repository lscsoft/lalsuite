/*
*  Copyright (C) 2007 Alexander Dietz, Drew Keppel, Duncan Brown, Eirini Messaritaki, Gareth Jones, Patrick Brady, Stephen Fairhurst, Craig Robinson , Thomas Cokelaer
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
 * File Name: thinca.c
 *
 * Author: Fairhurst, S.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <sys/types.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/Date.h>
#include <lal/TimeDelay.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOLwXMLRead.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/CoincInspiralEllipsoid.h>
#include <lal/Segments.h>
#include <lal/SegmentsIO.h>
#include <lalapps.h>
#include <lal/lalGitID.h>
#include <lalappsGitID.h>
#include <processtable.h>

RCSID("$Id$");

#define CVS_ID_STRING "$Id$"
#define CVS_NAME_STRING "$Name$"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define PROGRAM_NAME "thinca"

#define INCA_EARG   1
#define INCA_EROW   2
#define INCA_EFILE  3

#define INCA_MSGEARG   "Error parsing arguments"
#define INCA_MSGROW    "Error reading row from XML table"
#define INCA_MSGEFILE  "Could not open file"

#define MAXIFO 4

#define KAPPA 1000
#define EPSILON 2

#define ADD_PROCESS_PARAM( pptype, format, ppvalue ) \
  this_proc_param = this_proc_param->next = (ProcessParamsTable *) \
calloc( 1, sizeof(ProcessParamsTable) ); \
snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", \
    PROGRAM_NAME ); \
snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "--%s", \
    long_options[option_index].name ); \
snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "%s", pptype ); \
snprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, format, ppvalue );

int haveTrig[LAL_NUM_IFO];
int checkTimes = 0;
int multiIfoCoinc = 0;
int distCut = 0;
int iotaCut = 0;
int doPsi0Psi3Cut = 0;
int doAlphaFCut = 0;
int doBCV2H1H2Veto = 0;
int doBCVC = 0;
int h1h2Consistency = 0;
int doVeto = 0;
int completeCoincs = 0;

INT4 numExtTriggers = 0;
ExtTriggerTable   *exttrigHead = NULL;

/*
 * 
 * USAGE
 *
 */
#define USAGE(a, msg) \
fprintf( a, "Usage:  %s [options] [LIGOLW XML input files]\n" , msg);\
fprintf( a, "The following options are recognized.  Options not surrounded in [] are\n" );\
fprintf( a, "required.\n" );\
fprintf( a, "  [--help]                      display this message\n");\
fprintf( a, "  [--verbose]                   print progress information\n");\
fprintf( a, "  [--version]                   print version information and exit\n");\
fprintf( a, "  [--debug-level]   level       set the LAL debug level to LEVEL\n");\
fprintf( a, "  [--user-tag]      usertag     set the process_params usertag\n");\
fprintf( a, "  [--ifo-tag]       ifotag      set the ifo-tag - for file naming\n");\
fprintf( a, "  [--comment]       string      set the process table comment to STRING\n");\
fprintf( a, "  [--write-compress]            write a compressed xml file\n");\
fprintf( a, "\n");\
fprintf( a, "   --gps-start-time start_time  GPS second of data start time\n");\
fprintf( a, "   --gps-end-time   end_time    GPS second of data end time\n");\
fprintf( a, "  [--check-times]               Check that all times were analyzed\n");\
fprintf( a, "  [--multi-ifo-coinc]           Look for triple/quadruple ifo coincidence\n");\
fprintf( a, "  [--maximization-interval] max_dt set length of maximization interval in ms\n");\
fprintf( a, "\n");\
fprintf( a, "  [--g1-slide]      g1_slide    Slide G1 data by multiples of g1_slide\n");\
fprintf( a, "  [--h1-slide]      h1_slide    Slide H1 data by multiples of h1_slide\n");\
fprintf( a, "  [--h2-slide]      h2_slide    Slide H2 data by multiples of h2_slide\n");\
fprintf( a, "  [--l1-slide]      l1_slide    Slide L1 data by multiples of l1_slide\n");\
fprintf( a, "  [--t1-slide]      t1_slide    Slide T1 data by multiples of t1_slide\n");\
fprintf( a, "  [--v1-slide]      v1_slide    Slide V1 data by multiples of v1_slide\n");\
fprintf( a, "  [--num-slides]    num_slides  The number of time slides to perform\n");\
fprintf( a, "\n");\
fprintf( a, "  [--g1-triggers]               input triggers from G1\n");\
fprintf( a, "  [--h1-triggers]               input triggers from H1\n");\
fprintf( a, "  [--h2-triggers]               input triggers from H2\n");\
fprintf( a, "  [--l1-triggers]               input triggers from L1\n");\
fprintf( a, "  [--t1-triggers]               input triggers from T1\n");\
fprintf( a, "  [--v1-triggers]               input triggers from V1\n");\
fprintf( a, "\n");\
fprintf( a, "  [--do-veto]                    perform a veto on single IFO triggers\n");\
fprintf( a, "                                 at the times specified in the veto files below\n");\
fprintf( a, "  [--g1-veto-file]               veto file for G1\n");\
fprintf( a, "  [--h1-veto-file]               veto file for H1\n");\
fprintf( a, "  [--h2-veto-file]               veto file for H2\n");\
fprintf( a, "  [--l1-veto-file]               veto file for L1\n");\
fprintf( a, "  [--t1-veto-file]               veto file for T1\n");\
fprintf( a, "  [--v1-veto-file]               veto file for V1\n");\
fprintf( a, "\n");\
fprintf( a, "   --parameter-test     test    set parameters with which to test coincidence:\n");\
fprintf( a, "                                (m1_and_m2|mchirp_and_eta|mchirp_and_eta_ext|psi0_and_psi3|ellipsoid)\n");\
fprintf( a, "  [--g1-time-accuracy]  g1_dt   specify the timing accuracy of G1 in ms\n");\
fprintf( a, "  [--h1-time-accuracy]  h1_dt   specify the timing accuracy of H1 in ms\n");\
fprintf( a, "  [--h2-time-accuracy]  h2_dt   specify the timing accuracy of H2 in ms\n");\
fprintf( a, "  [--l1-time-accuracy]  l1_dt   specify the timing accuracy of L1 in ms\n");\
fprintf( a, "  [--t1-time-accuracy]  t1_dt   specify the timing accuracy of T1 in ms\n");\
fprintf( a, "  [--v1-time-accuracy]  v1_dt   specify the timing accuracy of V1 in ms\n");\
fprintf( a, "\n");\
fprintf( a, "  [--g1-mass-accuracy]  g1_dm   specify the mass accuracy of G1\n");\
fprintf( a, "  [--h1-mass-accuracy]  h1_dm   specify the mass accuracy of H1\n");\
fprintf( a, "  [--h2-mass-accuracy]  h2_dm   specify the mass accuracy of H2\n");\
fprintf( a, "  [--l1-mass-accuracy]  l1_dm   specify the mass accuracy of L1\n");\
fprintf( a, "  [--t1-mass-accuracy]  t1_dm   specify the mass accuracy of T1\n");\
fprintf( a, "  [--v1-mass-accuracy]  v1_dm   specify the mass accuracy of V1\n");\
fprintf( a, "\n");\
fprintf( a, "  [--g1-mchirp-accuracy] g1_dmchirp  specify the mchirp accuracy of G1\n");\
fprintf( a, "  [--h1-mchirp-accuracy] h1_dmchirp  specify the mchirp accuracy of H1\n");\
fprintf( a, "  [--h2-mchirp-accuracy] h2_dmchirp  specify the mchirp accuracy of H2\n");\
fprintf( a, "  [--l1-mchirp-accuracy] l1_dmchirp  specify the mchirp accuracy of L1\n");\
fprintf( a, "  [--t1-mchirp-accuracy] t1_dmchirp  specify the mchirp accuracy of T1\n");\
fprintf( a, "  [--v1-mchirp-accuracy] v1_dmchirp  specify the mchirp accuracy of V1\n");\
fprintf( a, "\n");\
fprintf( a, "  [--g1-eta-accuracy] g1_deta   specify the eta accuracy of G1\n");\
fprintf( a, "  [--h1-eta-accuracy] h1_deta   specify the eta accuracy of H1\n");\
fprintf( a, "  [--h2-eta-accuracy] h2_deta   specify the eta accuracy of H2\n");\
fprintf( a, "  [--l1-eta-accuracy] l1_deta   specify the eta accuracy of L1\n");\
fprintf( a, "  [--t1-eta-accuracy] t1_deta   specify the eta accuracy of T1\n");\
fprintf( a, "  [--v1-eta-accuracy] v1_deta   specify the eta accuracy of V1\n");\
fprintf( a, "\n");\
fprintf( a, "  [--g1-psi0-accuracy]  g1_dpsi0   specify the psi0 accuracy of G1\n");\
fprintf( a, "  [--h1-psi0-accuracy]  h1_dpsi0   specify the psi0 accuracy of H1\n");\
fprintf( a, "  [--h2-psi0-accuracy]  h2_dpsi0   specify the psi0 accuracy of H2\n");\
fprintf( a, "  [--l1-psi0-accuracy]  l1_dpsi0   specify the psi0 accuracy of L1\n");\
fprintf( a, "  [--t1-psi0-accuracy]  t1_dpsi0   specify the psi0 accuracy of T1\n");\
fprintf( a, "  [--v1-psi0-accuracy]  v1_dpsi0   specify the psi0 accuracy of V1\n");\
fprintf( a, "\n");\
fprintf( a, "  [--g1-psi3-accuracy]  g1_dpsi3   specify the psi3 accuracy of G1\n");\
fprintf( a, "  [--h1-psi3-accuracy]  h1_dpsi3   specify the psi3 accuracy of H1\n");\
fprintf( a, "  [--h2-psi3-accuracy]  h2_dpsi3   specify the psi3 accuracy of H2\n");\
fprintf( a, "  [--l1-psi3-accuracy]  l1_dpsi3   specify the psi3 accuracy of L1\n");\
fprintf( a, "  [--t1-psi3-accuracy]  t1_dpsi3   specify the psi3 accuracy of T1\n");\
fprintf( a, "  [--v1-psi3-accuracy]  v1_dpsi3   specify the psi3 accuracy of V1\n");\
fprintf( a, "\n");\
fprintf( a, "  [--e-thinca-parameter]  match    specify the e-thinca parameter\n");\
fprintf( a, "\n");\
fprintf( a, "  [--h1-h2-distance-cut]           perform H1-H2 distance cut\n");\
fprintf( a, "  [--h1-kappa]          h1_kappa   specify H1 kappa for eff dist test\n");\
fprintf( a, "  [--h2-kappa]          h2_kappa   specify H2 kappa for eff dist test\n");\
fprintf( a, "  [--h1-epsilon]        h1_epsilon specify H1 epsilon for eff dist test\n");\
fprintf( a, "  [--h2-epsilon]        h2_epsilon specify H1 epsilon for eff dist test\n");\
fprintf( a, "\n");\
fprintf( a, "  [--h1-h2-consistency]            perform H1-H2 consistency cut\n");\
fprintf( a, "  [--snr-cut]           snr        reject triggers below this snr\n");\
fprintf( a, "                                   needed when --h1-h2-consistency is given\n");\
fprintf( a, "\n");\
fprintf( a, "  [--dmchirp-high]      dmchirp    different chirp mass window for high masses\n");\
fprintf( a, "  [--high-mass]         mass       total mass for trigger above which\n");\
fprintf( a, "                                   the high mass mchirp is used\n");\
fprintf( a, "\n");\
fprintf( a, "                                   (both bounds must be specified)\n");\
fprintf( a, "  [--bcvc]                         perform cut based on alphaF\n");\
fprintf( a, "                                   with individual bounds given below.\n");\
fprintf( a, "  [--do-alphaf-cut]                perform cut based on alphaF\n");\
fprintf( a, "                                   with bounds given below\n");\
fprintf( a, "                                   (both bounds must be specified)\n");\
fprintf( a, "  [--h1-alphaf-hi]      alphaFhi   reject BCV triggers outside the specified \n");\
fprintf( a, "  [--h1-alphaf-lo]      alphaFlo   alphaF area for H1 (if \n");\
fprintf( a, "  [--h2-alphaf-hi]      alphaFhi   reject BCV triggers outside the specified \n");\
fprintf( a, "  [--h2-alphaf-lo]      alphaFlo   alphaF area for H2\n");\
fprintf( a, "  [--l1-alphaf-hi]      alphaFhi   reject BCV triggers outside the specified \n");\
fprintf( a, "  [--l1-alphaf-lo]      alphaFlo   alphaF area for L1\n");\
fprintf( a, "\n");\
fprintf( a, "  [--iota-cut-h1h2]    iotaCutH1H2 reject H1H2 triggers with iota above the specified value \n");\
fprintf( a, "  [--iota-cut-h1l1]    iotaCutH1L1 reject H1L1 triggers with iota above the specified value \n");\
fprintf( a, "\n");\
fprintf( a, "  [--do-bcvspin-h1h2-veto]         reject coincidences with H1 snr < H2 snr \n");\
fprintf( a, "\n");\
fprintf( a, "   --data-type          data_type  specify the data type, must be one of\n");\
fprintf( a, "                                   (playground_only|exclude_play|all_data)\n");\
fprintf( a, "   --complete-coincs               write out triggers from all non-vetoed ifos\n");\
fprintf( a, "                                   if not seen, snr is equal to zero\n");\
fprintf( a, "  [--exttrig]           source     enables the External-Trigger mode \n");\
fprintf( a, "                                         (using actual time delays) for a source\n");\
fprintf( a, "                                   specified in the source file\n");\
fprintf( a, "\n");\
fprintf( a, "[LIGOLW XML input files] list of the input trigger files.\n");\



/*
 * 
 * MAIN
 *
 */
  
int main( int argc, char *argv[] )
{
  static LALStatus      status;

  extern int vrbflg;

  LALPlaygroundDataMask dataType = unspecified_data_type;
  INT4  startCoincidence = -1;
  LIGOTimeGPS startCoinc = {0,0};
  INT4  endCoincidence = -1;
  LIGOTimeGPS endCoinc = {0,0};

  REAL8        slideStep[LAL_NUM_IFO];
  LIGOTimeGPS  slideTimes[LAL_NUM_IFO];
  LIGOTimeGPS  slideReset[LAL_NUM_IFO];
  INT4         numSlides = 0;
  INT4         slideNum  = 0;

  CHAR  ifoName[MAXIFO][LIGOMETA_IFO_MAX];
  CHAR  ifos[LIGOMETA_IFOS_MAX];
  CHAR  ifoA[LIGOMETA_IFO_MAX];
  CHAR  ifoB[LIGOMETA_IFO_MAX];
  CHAR  ifo[LIGOMETA_IFO_MAX];
  
  CHAR  comment[LIGOMETA_COMMENT_MAX];
  CHAR *userTag = NULL;
  CHAR *ifoTag = NULL;

  CHAR  fileName[FILENAME_MAX];
  CHAR  fileSlide[FILENAME_MAX];
  CHAR *vetoFileName[LAL_NUM_IFO] = {NULL, NULL, NULL, NULL, NULL, NULL};

  LALSegList vetoSegs[LAL_NUM_IFO];

  UINT4  numIFO = 0;
  UINT4  numTrigIFO = 0;
  UINT4  numTriggers = 0;
  UINT4  numCoinc = 0;
  UINT4  numDoubles = 0;
  UINT4  numTriples = 0;
  UINT4  numQuadruples = 0;
  UINT4  numTrigs[LAL_NUM_IFO];
  UINT4  N = 0;
  INT4  outCompress = 0;
  UINT4  slideH1H2Together = 0;

  LALDetector          aDet;

  SnglInspiralTable    *inspiralEventList = NULL;
  SnglInspiralTable    *thisInspiralTrigger = NULL;
  SnglInspiralTable    *snglOutput = NULL;
  SnglInspiralTable    *completionSngls = NULL;

  CoincInspiralTable   *coincInspiralList = NULL;
  CoincInspiralTable   *thisCoinc = NULL;

  EventIDColumn        *eventId;

  InspiralAccuracyList  accuracyParams;

  SearchSummvarsTable  *inputFiles = NULL;
  SearchSummvarsTable  *thisInputFile = NULL;

  SearchSummaryTable   *searchSummList = NULL;
  SearchSummaryTable   *thisSearchSumm = NULL;

  SummValueTable       *summValueList = NULL;
  SummValueTable       *thisSummValue = NULL;

  MetadataTable         proctable;
  MetadataTable         processParamsTable;
  MetadataTable         searchsumm;
  MetadataTable         summValueTable;
  MetadataTable         searchSummvarsTable;
  MetadataTable         inspiralTable;
  ProcessParamsTable   *this_proc_param = NULL;
  LIGOLwXMLStream       xmlStream;

  InterferometerNumber  ifoNumber = LAL_UNKNOWN_IFO;
  InterferometerNumber  ifoTwo    = LAL_UNKNOWN_IFO;
  INT4                  i;
  INT4                  loopVar;
  INT4                  maximizationInterval = 0;

  SnglInspiralBCVCalphafCut  alphafParams;

  /* by default we do not remove any triggers in the SNR Cut*/  
  
  REAL4                 snrCut = 0;   
  char*                 sourceFile=NULL;
 
  const CHAR                  *ifoArg[LAL_NUM_IFO] = 
                                   {"g1-triggers", "h1-triggers", 
                                    "h2-triggers", "l1-triggers", 
                                    "t1-triggers", "v1-triggers"};


  /* getopt arguments */
  struct option long_options[] =
  {
    {"verbose",             no_argument,   &vrbflg,                   1 },
    {"write-compress",      no_argument,   &outCompress,              1 },
    {"g1-triggers",         no_argument,   &(haveTrig[LAL_IFO_G1]),   1 },
    {"h1-triggers",         no_argument,   &(haveTrig[LAL_IFO_H1]),   1 },
    {"h2-triggers",         no_argument,   &(haveTrig[LAL_IFO_H2]),   1 },
    {"l1-triggers",         no_argument,   &(haveTrig[LAL_IFO_L1]),   1 },
    {"t1-triggers",         no_argument,   &(haveTrig[LAL_IFO_T1]),   1 },
    {"v1-triggers",         no_argument,   &(haveTrig[LAL_IFO_V1]),   1 },
    {"check-times",         no_argument,   &checkTimes,               1 },
    {"multi-ifo-coinc",     no_argument,   &multiIfoCoinc,            1 },
    {"h1-h2-distance-cut",  no_argument,   &distCut,                  1 },
    {"h1-h2-consistency",   no_argument,   &h1h2Consistency,          1 },
    {"do-bcvspin-h1h2-veto",no_argument,   &doBCV2H1H2Veto,           1 },
    {"do-alphaf-cut",       no_argument,   &doAlphaFCut,              1 },
    {"psi0-psi3-cut",       no_argument,   &doPsi0Psi3Cut,            1 },
    {"bcvc",                no_argument,   &doBCVC,                   1 },
    {"do-veto",             no_argument,   &doVeto,                   1 },
    {"complete-coincs",     no_argument,   &completeCoincs,           1 },
    {"g1-slide",            required_argument, 0,                    'b'},
    {"h1-slide",            required_argument, 0,                    'c'},
    {"h2-slide",            required_argument, 0,                    'd'},
    {"l1-slide",            required_argument, 0,                    'e'},
    {"t1-slide",            required_argument, 0,                    'f'},
    {"v1-slide",            required_argument, 0,                    'g'},
    {"num-slides",          required_argument, 0,                    'T'},
    {"g1-time-accuracy",    required_argument, 0,                    'A'},
    {"h1-time-accuracy",    required_argument, 0,                    'B'}, 
    {"h2-time-accuracy",    required_argument, 0,                    'C'}, 
    {"l1-time-accuracy",    required_argument, 0,                    'D'},
    {"t1-time-accuracy",    required_argument, 0,                    'E'}, 
    {"v1-time-accuracy",    required_argument, 0,                    'F'}, 
    {"g1-mass-accuracy",    required_argument, 0,                    'G'},
    {"h1-mass-accuracy",    required_argument, 0,                    'H'},
    {"h2-mass-accuracy",    required_argument, 0,                    'I'},
    {"l1-mass-accuracy",    required_argument, 0,                    'J'},
    {"t1-mass-accuracy",    required_argument, 0,                    'K'},
    {"v1-mass-accuracy",    required_argument, 0,                    'L'},
    {"g1-mchirp-accuracy",  required_argument, 0,                    'M'},
    {"h1-mchirp-accuracy",  required_argument, 0,                    'N'},
    {"h2-mchirp-accuracy",  required_argument, 0,                    'O'},
    {"l1-mchirp-accuracy",  required_argument, 0,                    'P'},
    {"t1-mchirp-accuracy",  required_argument, 0,                    'Q'},
    {"v1-mchirp-accuracy",  required_argument, 0,                    'R'},
    {"g1-eta-accuracy",     required_argument, 0,                    'm'},
    {"h1-eta-accuracy",     required_argument, 0,                    'n'},
    {"h2-eta-accuracy",     required_argument, 0,                    'o'},
    {"l1-eta-accuracy",     required_argument, 0,                    'p'},
    {"t1-eta-accuracy",     required_argument, 0,                    'q'},
    {"v1-eta-accuracy",     required_argument, 0,                    'r'},
    {"g1-psi0-accuracy",    required_argument, 0,                    '2'},
    {"h1-psi0-accuracy",    required_argument, 0,                    '3'},
    {"h2-psi0-accuracy",    required_argument, 0,                    '4'},
    {"l1-psi0-accuracy",    required_argument, 0,                    '5'},
    {"t1-psi0-accuracy",    required_argument, 0,                    '6'},
    {"v1-psi0-accuracy",    required_argument, 0,                    '7'},
    {"g1-psi3-accuracy",    required_argument, 0,                    '8'},
    {"h1-psi3-accuracy",    required_argument, 0,                    '9'},
    {"h2-psi3-accuracy",    required_argument, 0,                    '!'},
    {"l1-psi3-accuracy",    required_argument, 0,                    '-'},
    {"t1-psi3-accuracy",    required_argument, 0,                    '+'},
    {"v1-psi3-accuracy",    required_argument, 0,                    '='},
    {"e-thinca-parameter",  required_argument, 0,                    '`'},
    {"h1-kappa",            required_argument, 0,                    'W'},
    {"h2-kappa",            required_argument, 0,                    'Y'},
    {"h1-epsilon",          required_argument, 0,                    'w'},
    {"h2-epsilon",          required_argument, 0,                    'y'},
    {"parameter-test",      required_argument, 0,                    'a'},
    {"gps-start-time",      required_argument, 0,                    's'},
    {"gps-end-time",        required_argument, 0,                    't'},
    {"maximization-interval",required_argument, 0,                   '@'},    
    {"data-type",           required_argument, 0,                    'k'},
    {"comment",             required_argument, 0,                    'x'},
    {"user-tag",            required_argument, 0,                    'Z'},
    {"userTag",             required_argument, 0,                    'Z'},
    {"ifo-tag",             required_argument, 0,                    'i'},
    {"help",                no_argument,       0,                    'h'}, 
    {"debug-level",         required_argument, 0,                    'z'},
    {"version",             no_argument,       0,                    'V'},
    {"dmchirp-high",        required_argument, 0,                    '^'},
    {"high-mass",           required_argument, 0,                    '&'},
    {"h1-alphaf-lo",        required_argument, 0,                    'l'},
    {"h1-alphaf-hi",        required_argument, 0,                    'j'},
    {"h2-alphaf-lo",        required_argument, 0,                    'u'},
    {"h2-alphaf-hi",        required_argument, 0,                    'S'},
    {"l1-alphaf-lo",        required_argument, 0,                    'U'},
    {"l1-alphaf-hi",        required_argument, 0,                    'X'},
    {"iota-cut-h1h2",       required_argument, 0,                    '#'},
    {"iota-cut-h1l1",       required_argument, 0,                    '%'},
    {"snr-cut",             required_argument, 0,                    '*'},
    {"h1-veto-file",        required_argument, 0,                    '('},
    {"h2-veto-file",        required_argument, 0,                    ')'},
    {"g1-veto-file",        required_argument, 0,                    '{'},
    {"l1-veto-file",        required_argument, 0,                    '}'},
    {"t1-veto-file",        required_argument, 0,                    '['},
    {"v1-veto-file",        required_argument, 0,                    ']'},
    {"exttrig",             required_argument, 0,                    '_'},
    {0, 0, 0, 0}
  };
  int c;
  /* INFO: Remaining characters: * ( ) _ { } [ ]  (or even more...)*/

  alphafParams.h1_lo=-1e10;
  alphafParams.h1_hi=+1e10;
  alphafParams.h2_lo=-1e10;
  alphafParams.h2_hi=+1e10;
  alphafParams.l1_lo=-1e10;
  alphafParams.l1_hi=+1e10;
  alphafParams.psi0cut = 0; 
 
  /* set default values for those values.
     with those values ALL triggers will survive (i.e. no cut). */
  accuracyParams.exttrig=0;


  /*
   * 
   * initialize things
   *
   */

  lal_errhandler = LAL_ERR_EXIT;
  set_debug_level( "33" );
  setvbuf( stdout, NULL, _IONBF, 0 );

  /* create the process and process params tables */
  proctable.processTable = (ProcessTable *) calloc( 1, sizeof(ProcessTable) );
  XLALGPSTimeNow(&(proctable.processTable->start_time));
  if (strcmp(CVS_REVISION, "$Revi" "sion$"))
  {
    XLALPopulateProcessTable(proctable.processTable, PROGRAM_NAME,
        CVS_REVISION, CVS_SOURCE, CVS_DATE, 0);
  }
  else
  {
    XLALPopulateProcessTable(proctable.processTable, PROGRAM_NAME,
        lalappsGitCommitID, lalappsGitGitStatus, lalappsGitCommitDate, 0);
  }
  this_proc_param = processParamsTable.processParamsTable = 
    (ProcessParamsTable *) calloc( 1, sizeof(ProcessParamsTable) );
  memset( comment, 0, LIGOMETA_COMMENT_MAX * sizeof(CHAR) );

  /* create the search summary and zero out the summvars table */
  searchsumm.searchSummaryTable = (SearchSummaryTable *)
    calloc( 1, sizeof(SearchSummaryTable) );

  memset( &accuracyParams, 0, sizeof(InspiralAccuracyList) );
  accuracyParams.iotaCutH1H2=-1.0;
  /* by default, iotacutH1H2 is unphysical so that it must be provided 
     if --h1-h2-consistency-check option is given by the user (checked 
     in the parsing function)
  */

  memset( &aDet, 0, sizeof(LALDetector) );

  /* set the time slide data to zero */
  memset( &slideStep, 0, LAL_NUM_IFO * sizeof(REAL8) );
  memset( &slideTimes, 0, LAL_NUM_IFO * sizeof(LIGOTimeGPS) );
  memset( &slideReset, 0, LAL_NUM_IFO * sizeof(LIGOTimeGPS) );
  memset( &haveTrig, 0, LAL_NUM_IFO * sizeof(int) );

  /* initialize array first */
  for ( ifoNumber = 0; ifoNumber < LAL_NUM_IFO ; ++ifoNumber )
    {      
      numTrigs[ifoNumber]=0;
    }

  /* parse the arguments */
  while ( 1 )
  {
    /* getopt_long stores long option here */
    int option_index = 0;
    long int gpstime;
    size_t optarg_len;

    c = getopt_long_only( argc, argv, 
        "A:B:C:D:E:F:G:H:I:J:K:L:M:N:O:P:Q:R:S:T:U:VW:Y:Z:"
        "a:b:c:d:e:f:g:hi:j:k:l:m:n:o:p:q:r:s:t:u:v:w:x:y:z:"
        "2:3:4:5:6:7:8:9:`:!:-:+:=:@:^:&:*:(:):{:}:[:]:~", 
        long_options, &option_index );

    /* detect the end of the options */
    if ( c == -1 )
    {
      break;
    }

    switch ( c )
    {
      case 0:
        /* if this option set a flag, do nothing else now */
        if ( long_options[option_index].flag != 0 )
        {
          break;
        }
        else
        {
          fprintf( stderr, "Error parsing option %s with argument %s\n",
              long_options[option_index].name, optarg );
          exit( 1 );
        }
        break;

      case 'a':
        /* set the parameter test */
        if ( ! strcmp( "m1_and_m2", optarg ) )
        {
          accuracyParams.test = m1_and_m2;
        }
        else if ( ! strcmp( "psi0_and_psi3", optarg ) )
        {
          accuracyParams.test = psi0_and_psi3;
        }
        else if ( ! strcmp( "mchirp_and_eta", optarg ) )
        {
          accuracyParams.test = mchirp_and_eta;
        }
        else if ( ! strcmp( "ellipsoid", optarg ) )
        {
          accuracyParams.test = ellipsoid;
        }
        else
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "unknown test specified: "
              "%s (must be m1_and_m2, psi0_and_psi3, mchirp_and_eta or mchirp_and_eta_ext)\n",
              long_options[option_index].name, optarg );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'A':
        /* time accuracy G1, argument is in milliseconds */
        accuracyParams.ifoAccuracy[LAL_IFO_G1].dt = atof(optarg) * 1000000LL;
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;
      
      case 'B':
        /* time accuracy H1, argument is in milliseconds */
        accuracyParams.ifoAccuracy[LAL_IFO_H1].dt = atof(optarg) * 1000000LL;
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;
        
      case 'C':
        /* time accuracy H2, argument is in milliseconds */
        accuracyParams.ifoAccuracy[LAL_IFO_H2].dt = atof(optarg) * 1000000LL;
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;
        
      case 'D':
        /* time accuracy L1, argument is in milliseconds */
        accuracyParams.ifoAccuracy[LAL_IFO_L1].dt = atof(optarg) * 1000000LL;
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case 'E':
        /* time accuracy T1, argument is in milliseconds */
        accuracyParams.ifoAccuracy[LAL_IFO_T1].dt = atof(optarg) * 1000000LL;
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;
        
      case 'F':
        /* time accuracy V1, argument is in milliseconds */
        accuracyParams.ifoAccuracy[LAL_IFO_V1].dt = atof(optarg) * 1000000LL;
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case 'G':
        /* mass accuracy G1, argument is in solar masses */
        accuracyParams.ifoAccuracy[LAL_IFO_G1].dm = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;
      
      case 'H':
        /* mass accuracy H1, argument is in solar masses */
        accuracyParams.ifoAccuracy[LAL_IFO_H1].dm = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;
        
      case 'I':
        /* mass accuracy H2, argument is in solar masses */
        accuracyParams.ifoAccuracy[LAL_IFO_H2].dm = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;
        
      case 'J':
        /* mass accuracy L1, argument is in solar masses */
        accuracyParams.ifoAccuracy[LAL_IFO_L1].dm = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case 'K':
        /* mass accuracy T1, argument is in solar masses */
        accuracyParams.ifoAccuracy[LAL_IFO_T1].dm = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;
        
      case 'L':
        /* mass accuracy V1, argument is in solar masses */
        accuracyParams.ifoAccuracy[LAL_IFO_V1].dm = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case 'M':
        /* chirp mass accuracy G1, argument is in solar masses */
        accuracyParams.ifoAccuracy[LAL_IFO_G1].dmchirp = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;
      
      case 'N':
        /* chirp mass accuracy H1, argument is in solar masses */
        accuracyParams.ifoAccuracy[LAL_IFO_H1].dmchirp = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;
        
      case 'O':
        /* chirp mass accuracy H2, argument is in solar masses */
        accuracyParams.ifoAccuracy[LAL_IFO_H2].dmchirp = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;
        
      case 'P':
        /* chirp mass accuracy L1, argument is in solar masses */
        accuracyParams.ifoAccuracy[LAL_IFO_L1].dmchirp = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case 'Q':
        /* chirp mass accuracy T1, argument is in solar masses */
        accuracyParams.ifoAccuracy[LAL_IFO_T1].dmchirp = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;
        
      case 'R':
        /* chirp mass accuracy V1, argument is in solar masses */
        accuracyParams.ifoAccuracy[LAL_IFO_V1].dmchirp = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;
       
      case 'b':
        /* slide time for G1 */
        slideStep[LAL_IFO_G1] = atof( optarg );
        if ( slideStep[LAL_IFO_G1] < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "The slideStep must be positive\n"
              "(%f specified)\n",
              long_options[option_index].name, slideStep[LAL_IFO_G1] );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "double", "%f", slideStep[LAL_IFO_G1] );
        break;

      case 'c':
        /* slide time for H1 */
        slideStep[LAL_IFO_H1] = atof( optarg );
        if ( slideStep[LAL_IFO_H1] < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "The slideStep must be positive\n"
              "(%f specified)\n",
              long_options[option_index].name, slideStep[LAL_IFO_H1] );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "double", "%f", slideStep[LAL_IFO_H1] );
        break;
        
      case 'd':
        /* slide time for H2 */
        slideStep[LAL_IFO_H2] = atof( optarg );
        if ( slideStep[LAL_IFO_H2] < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "The slideStep must be positive\n"
              "(%f specified)\n",
              long_options[option_index].name, slideStep[LAL_IFO_H2] );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "double", "%f", slideStep[LAL_IFO_H2] );
        break;
        
      case 'e':
        /* slide time for L1 */
        slideStep[LAL_IFO_L1] = atof( optarg );
        if ( slideStep[LAL_IFO_L1] < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "The slideStep must be positive\n"
              "(%f specified)\n",
              long_options[option_index].name, slideStep[LAL_IFO_L1] );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "double", "%f", slideStep[LAL_IFO_L1] );
        break;
        
      case 'f':
        /* slide time for T1 */
        slideStep[LAL_IFO_T1] = atof( optarg );
        if ( slideStep[LAL_IFO_T1] < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "The slideStep must be positive\n"
              "(%f specified)\n",
              long_options[option_index].name, slideStep[LAL_IFO_T1] );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "double", "%f", slideStep[LAL_IFO_T1] );
        break;
        
      case 'g':
        /* slide time for V1 */
        slideStep[LAL_IFO_V1] = atof( optarg );
        if ( slideStep[LAL_IFO_V1] < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "The slideStep must be positive\n"
              "(%f specified)\n",
              long_options[option_index].name, slideStep[LAL_IFO_V1] );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "double", "%f", slideStep[LAL_IFO_V1] );
        break;

      case 'T':
        /* num slides*/
        numSlides = atoi( optarg );
        if ( numSlides < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "The number of time slides must be positive\n"
              "(%d specified)\n",
              long_options[option_index].name, numSlides );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%d", numSlides );
        break;

      case 'm':
        /* eta accuracy G1, argument is dimensionless */
        accuracyParams.ifoAccuracy[LAL_IFO_G1].deta = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;
      
      case 'n':
        /* eta accuracy H1, argument is dimensionless */
        accuracyParams.ifoAccuracy[LAL_IFO_H1].deta = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;
        
      case 'o':
        /* eta accuracy H2, argument is dimensionless */
        accuracyParams.ifoAccuracy[LAL_IFO_H2].deta = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;
        
      case 'p':
        /* eta accuracy L1, argument is dimensionless */
        accuracyParams.ifoAccuracy[LAL_IFO_L1].deta = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case 'q':
        /* eta accuracy T1, argument is dimensionless */
        accuracyParams.ifoAccuracy[LAL_IFO_T1].deta = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;
        
      case 'r':
        /* eta accuracy V1, argument is dimensionless */
        accuracyParams.ifoAccuracy[LAL_IFO_V1].deta = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case 's':
        /* start time coincidence window */
        gpstime = atol( optarg );
        if ( gpstime < 441417609 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "GPS start time is prior to " 
              "Jan 01, 1994  00:00:00 UTC:\n"
              "(%ld specified)\n",
              long_options[option_index].name, gpstime );
          exit( 1 );
        }
        if ( gpstime > 999999999 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "GPS start time is after " 
              "Sep 14, 2011  01:46:26 UTC:\n"
              "(%ld specified)\n", 
              long_options[option_index].name, gpstime );
          exit( 1 );
        }
        startCoincidence = (INT4) gpstime;
        ADD_PROCESS_PARAM( "int", "%ld", startCoincidence );
        break;

      case 't':
        /* end time coincidence window */
        gpstime = atol( optarg );
        if ( gpstime < 441417609 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "GPS start time is prior to " 
              "Jan 01, 1994  00:00:00 UTC:\n"
              "(%ld specified)\n",
              long_options[option_index].name, gpstime );
          exit( 1 );
        }
        if ( gpstime > 999999999 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "GPS start time is after " 
              "Sep 14, 2011  01:46:26 UTC:\n"
              "(%ld specified)\n", 
              long_options[option_index].name, gpstime );
          exit( 1 );
        }
        endCoincidence = (INT4) gpstime;
        ADD_PROCESS_PARAM( "int", "%ld", endCoincidence );
        break;

      case 'x':
        /* comment */
        if ( strlen( optarg ) > LIGOMETA_COMMENT_MAX - 1 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "comment must be less than %d characters\n",
              long_options[option_index].name, LIGOMETA_COMMENT_MAX );
          exit( 1 );
        }
        else
        {
          snprintf( comment, LIGOMETA_COMMENT_MAX, "%s", optarg);
        }
        break;
      

     case 'l':
        /* lower H1 alphaF cutoff */
        alphafParams.h1_lo = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case 'j':
        /* upper H1 alphaF cutoff */
        alphafParams.h1_hi = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case 'u':
        /* lower H2 alphaF cutoff */
        alphafParams.h2_lo = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case 'S':
        /* upper H2 alphaF cutoff */
        alphafParams.h2_hi = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case 'U':
        /* lower L1 alphaF cutoff */
        alphafParams.l1_lo = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case 'X':
        /* upper L1 alphaF cutoff */
        alphafParams.l1_hi = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;


     case '#':
        /* iota cut value for the H1H2 case (default=2) */
        accuracyParams.iotaCutH1H2 = atof(optarg);
        iotaCut  = 1;
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case '%':
        /* iota cut value for the H1L1 case  (default=2)*/
        accuracyParams.iotaCutH1L1 = atof(optarg);
        iotaCut = 1;
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case 'k':
        /* type of data to analyze */
        if ( ! strcmp( "playground_only", optarg ) )
        {
          dataType = playground_only;
        }
        else if ( ! strcmp( "exclude_play", optarg ) )
        {
          dataType = exclude_play;
        }
        else if ( ! strcmp( "all_data", optarg ) )
        {
          dataType = all_data;
        }
        else
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "unknown data type, %s, specified: "
              "(must be playground_only, exclude_play or all_data)\n",
              long_options[option_index].name, optarg );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;


      case 'h':
        /* help message */
        USAGE( stdout , "");
        exit( 1 );
        break;

      case 'z':
        set_debug_level( optarg );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'Z':
        /* create storage for the usertag */
        optarg_len = strlen(optarg) + 1;
        userTag = (CHAR *) calloc( optarg_len, sizeof(CHAR) );
        memcpy( userTag, optarg, optarg_len );

        this_proc_param = this_proc_param->next = (ProcessParamsTable *)
          calloc( 1, sizeof(ProcessParamsTable) );
        snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", 
            PROGRAM_NAME );
        snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "--userTag" );
        snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
        snprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, "%s",
            optarg );
        break;
        
      case 'i':
        /* create storage for the ifotag */
        optarg_len = strlen(optarg) + 1;
        ifoTag = (CHAR *) calloc( optarg_len, sizeof(CHAR) );
        memcpy( ifoTag, optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;


      case 'V':
        /* print version information and exit */
        fprintf( stdout, "The Hierarchical INspiral Coincidence Analysis\n" 
            "Steve Fairhurst\n"
            "CVS Version: " CVS_ID_STRING "\n"
            "CVS Tag: " CVS_NAME_STRING "\n" );
        fprintf( stdout, lalappsGitID );
        exit( 0 );
        break;

      case '?':
        USAGE( stderr , "");
        exit( 1 );
        break;

      case '2':
     /* psi0 mass accuracy G1  */
        accuracyParams.ifoAccuracy[LAL_IFO_G1].dpsi0 = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case '3':
     /* psi0 mass accuracy H1  */
        accuracyParams.ifoAccuracy[LAL_IFO_H1].dpsi0 = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case '4':
     /* psi0 mass accuracy H2  */
        accuracyParams.ifoAccuracy[LAL_IFO_H2].dpsi0 = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;
     
      case '5':
     /* psi0 mass accuracy L1  */
        accuracyParams.ifoAccuracy[LAL_IFO_L1].dpsi0 = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case '6':
     /* psi0 mass accuracy T1  */
        accuracyParams.ifoAccuracy[LAL_IFO_T1].dpsi0 = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case '7':
     /* psi0 mass accuracy V1  */
        accuracyParams.ifoAccuracy[LAL_IFO_V1].dpsi0 = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case '8':
     /* psi3 mass accuracy G1  */
        accuracyParams.ifoAccuracy[LAL_IFO_G1].dpsi3 = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case '9':
     /* psi3 mass accuracy H1  */
        accuracyParams.ifoAccuracy[LAL_IFO_H1].dpsi3 = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case '!':
     /* psi3 mass accuracy H2  */
        accuracyParams.ifoAccuracy[LAL_IFO_H2].dpsi3 = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case '-':
     /* psi3 mass accuracy L1  */
        accuracyParams.ifoAccuracy[LAL_IFO_L1].dpsi3 = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case '+':
     /* psi3 mass accuracy T1  */
        accuracyParams.ifoAccuracy[LAL_IFO_T1].dpsi3 = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case '=':
     /* psi3 mass accuracy V1  */
        accuracyParams.ifoAccuracy[LAL_IFO_V1].dpsi3 = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case '`':
      /* Ellipsoid scaling factor */
        accuracyParams.eMatch = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case '@':
        /* set the maximization window */
        maximizationInterval = atoi( optarg );
        if ( maximizationInterval < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "maximization interval must be positive:\n "
              "(%d ms specified)\n",
              long_options[option_index].name, maximizationInterval );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%ld",  maximizationInterval );
        break;
  
      case '^':
        for (loopVar=0; loopVar < LAL_NUM_IFO; loopVar++)
          accuracyParams.ifoAccuracy[loopVar].dmchirpHi = atof( optarg );
        break;  

      case '&':
        for (loopVar=0; loopVar < LAL_NUM_IFO; loopVar++) 
          accuracyParams.ifoAccuracy[loopVar].highMass = atof( optarg );
        break; 

      case 'W':
        /* kappa accuracy H1 */
        accuracyParams.ifoAccuracy[LAL_IFO_H1].kappa = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;
        
      case 'Y':
        /* kappa accuracy H2 */
        accuracyParams.ifoAccuracy[LAL_IFO_H2].kappa = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;
        
      case 'w':
        /* epsilon accuracy H1 */
        accuracyParams.ifoAccuracy[LAL_IFO_H1].epsilon = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;
        
      case 'y':
        /* epsilon accuracy H2 */
        accuracyParams.ifoAccuracy[LAL_IFO_H2].epsilon = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;
      
      case '*':
        /* snr cut */
        snrCut = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;
      
      case '(':
        /* veto filename */
        optarg_len = strlen( optarg ) + 1;
        vetoFileName[LAL_IFO_H1] = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( vetoFileName[LAL_IFO_H1], optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case ')':
        /* veto filename */
        optarg_len = strlen( optarg ) + 1;
        vetoFileName[LAL_IFO_H2] = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( vetoFileName[LAL_IFO_H2], optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;
 
      case '}':
        /* veto filename */
        optarg_len = strlen( optarg ) + 1;
        vetoFileName[LAL_IFO_L1] = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( vetoFileName[LAL_IFO_L1], optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;        
        
      case '{':
        /* veto filename */
        optarg_len = strlen( optarg ) + 1;
        vetoFileName[LAL_IFO_G1] = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( vetoFileName[LAL_IFO_G1], optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;        
        
      case '[':
        /* veto filename */
        optarg_len = strlen( optarg ) + 1;
        vetoFileName[LAL_IFO_T1] = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( vetoFileName[LAL_IFO_T1], optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;        
        
      case ']':
        /* veto filename */
        optarg_len = strlen( optarg ) + 1;
        vetoFileName[LAL_IFO_V1] = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( vetoFileName[LAL_IFO_V1], optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;       

      case '_':
        /* specifying GRB source file */
        optarg_len = strlen(optarg) + 1;
        sourceFile = (CHAR *) calloc( optarg_len, sizeof(CHAR) );
        memcpy( sourceFile, optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        accuracyParams.exttrig=1;
        break;
       
      default:
        fprintf( stderr, "Error: Unknown error while parsing options \n" );
        USAGE( stderr , argv[0]);
        exit( 1 );

    }
  }

  /*
   *
   * check the values of the arguments
   *
   */


  /* Start and End times  */

  if ( startCoincidence < 0 )
  {
    fprintf( stderr, "Error: --gps-start-time must be specified\n" );
    exit( 1 );
  }

  if ( endCoincidence < 0 )
  {
    fprintf( stderr, "Error: --gps-end-time must be specified\n" );
    exit( 1 );
  }

  /* set the gps times startCoinc and endCoinc */
  startCoinc.gpsSeconds = startCoincidence;
  endCoinc.gpsSeconds = endCoincidence;


  /* Parameter Test */
  if ( accuracyParams.test == no_test )
  {
    fprintf( stderr, "Error: --parameter-test must be specified\n" );
    exit( 1 );
  }

  /* Make sure the scaling parameter is set if the parameter test is 
   * set to ellipsoid.
   */
  if ( accuracyParams.test == ellipsoid &&
        accuracyParams.eMatch <= 0.0  )
  {
     fprintf( stderr, "Error: Invalid e-thinca parameter\n" );
     fprintf( stderr, "--e-thinca-parameter must be specified, "\
               "and must be larger than 0.\n" );
     exit(1);
  }    

  /* Data Type */
  if ( dataType == unspecified_data_type )
  {
    fprintf( stderr, "Error: --data-type must be specified\n");
    exit(1);
  }


  /* Store the IFOs we expect triggers from */
  for( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
  {
    if ( haveTrig[ifoNumber] )
    {
      /* write ifo name in ifoName list */
      XLALReturnIFO(ifo,ifoNumber);
      snprintf( ifoName[numIFO], LIGOMETA_IFO_MAX, ifo );
      numIFO++;

      /* store the argument in the process_params table */
      this_proc_param = this_proc_param->next = (ProcessParamsTable *)
        calloc( 1, sizeof(ProcessParamsTable) );
      snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, 
          "%s", PROGRAM_NAME );
      snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "--%s", 
          ifoArg[ifoNumber]);
      snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
      snprintf( this_proc_param->value, LIGOMETA_TYPE_MAX, " " );

      /* check that a non-zero timing accuracy was specified */
      if ( accuracyParams.test != ellipsoid &&
           ! accuracyParams.ifoAccuracy[ifoNumber].dt )
      {
        fprintf( stderr, "Error: --dt must be specified for %s\n", 
            ifoName[ifoNumber]);
        exit(1);
      }
    }
  }

  /* Check that if using different chirp mass accuracy at high mass, both 
     chirp mass accuracy and the high mass threshold are given */
  if ( (accuracyParams.ifoAccuracy[0].dmchirpHi || 
        accuracyParams.ifoAccuracy[0].highMass)  && 
       !(accuracyParams.ifoAccuracy[0].dmchirpHi 
         && accuracyParams.ifoAccuracy[0].highMass))
  {
    fprintf( stderr, 
        "Error: both --dmchirp-high and --high-mass must be specified\n");
    exit(1);
  }
    

  /* check that we have at least two IFOs specified, or can't do coincidence */
  if ( numIFO < 2 )
  {
    fprintf( stderr, "Must specify at least two IFOs to do coincidence\n"
        "%d specified\n", numIFO );
    exit ( 1 );
  }

  if ( numIFO > 2 && vrbflg )
  {
    if ( !multiIfoCoinc )
    {
      fprintf( stdout, 
          "Finding all double coincidences in %d IFO time.\n"
          "If you want triples/quadruples please specify --multi-ifo-coinc.\n",
          numIFO);
    }
    else
    {
      fprintf( stdout, 
          "Finding all double/triple/quadruple coincidences in %d IFO time.\n",
          numIFO);
                }
        }

  /* set ifos to be the alphabetical list of the ifos with triggers */
  if( numIFO == 2 )
  {
    snprintf( ifos, LIGOMETA_IFOS_MAX, "%s%s", ifoName[0], ifoName[1] );
  }
  else if ( numIFO == 3 )
  {
    snprintf( ifos, LIGOMETA_IFOS_MAX, "%s%s%s", ifoName[0], ifoName[1],
        ifoName[2] );
  }
  else if ( numIFO == 4 )
  {
    snprintf( ifos, LIGOMETA_IFOS_MAX, "%s%s%s%s", ifoName[0], ifoName[1],
        ifoName[2], ifoName[3]);
  }

  /* if numSlides is set, check that the slide times are different for
   * different ifos (for which we have triggers */
  if( numSlides )
  {
    for( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++ )
    {
      XLALReturnIFO( ifoA, ifoNumber );
     
      if( vrbflg && haveTrig[ifoNumber] ) fprintf( stdout, 
          "Performing a slide of multiples of %f seconds on %s\n", 
          slideStep[ifoNumber], ifoA);

      for( ifoTwo = ifoNumber + 1; ifoTwo < LAL_NUM_IFO; ifoTwo++ )
      {
        if ( haveTrig[ifoTwo] && haveTrig[ifoNumber] &&
            slideStep[ifoTwo] == slideStep[ifoNumber] )
        {
          XLALReturnIFO( ifoB, ifoTwo );

          if ( ( ! strcmp(ifoA,"H1") && ! strcmp(ifoB,"H2") ) ||
               ( ! strcmp(ifoA,"H2") && ! strcmp(ifoB,"H1") ) )
          {
            if( vrbflg ) fprintf( stdout,
              "The time slide specified for ifo %s is %f\n"
              "The time slide specified for ifo %s is also %f\n",
              ifoA, slideStep[ifoNumber], ifoB, slideStep[ifoTwo]);
            slideH1H2Together = 1;
          }
          else
          {
            fprintf( stderr,
              "The time slide specified for ifo %s is %f\n"
              "The time slide specified for ifo %s is also %f\n"
              "Must specify unique time slides for all instruments\n",
              ifoA, slideStep[ifoNumber], ifoB, slideStep[ifoTwo]);

            exit( 1 );
          }
        }
      }
    }
  }

  /* check if cut values has been applied  */ 
  if ( iotaCut )
  {
    if (accuracyParams.iotaCutH1H2<0 || accuracyParams.iotaCutH1H2>2 
        || accuracyParams.iotaCutH1L1<0 || accuracyParams.iotaCutH1L1>2) {
      fprintf( stderr,"The iota-cut values for H1H2 and/or H1L1 lie outside \n"
               "a meaning full range of [0;2]\n");
      exit( 1 );
    }

    
  }
  
  /* fill the comment, if a user has specified one, or leave it blank */
  if ( ! *comment )
  {
    snprintf( proctable.processTable->comment, LIGOMETA_COMMENT_MAX, " " );
    snprintf( searchsumm.searchSummaryTable->comment, LIGOMETA_COMMENT_MAX, 
        " " );
  } 
  else 
  {
    snprintf( proctable.processTable->comment, LIGOMETA_COMMENT_MAX,
        "%s", comment );
    snprintf( searchsumm.searchSummaryTable->comment, LIGOMETA_COMMENT_MAX,
        "%s", comment );
  }


  if ( outCompress )
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *)
      calloc( 1, sizeof(ProcessParamsTable) );
    snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, 
        "%s", PROGRAM_NAME );
    snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "--write-compress" );
    snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
    snprintf( this_proc_param->value, LIGOMETA_TYPE_MAX, " " );
  }
  /* store the check-times in the process_params table */
  if ( checkTimes )
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *)
      calloc( 1, sizeof(ProcessParamsTable) );
    snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, 
        "%s", PROGRAM_NAME );
    snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "--check-times" );
    snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
    snprintf( this_proc_param->value, LIGOMETA_TYPE_MAX, " " );
  }

  /* store the multi-ifo-coinc in the process_params table */
  if ( multiIfoCoinc )
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *)
      calloc( 1, sizeof(ProcessParamsTable) );
    snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, 
        "%s", PROGRAM_NAME );
    snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, 
        "--multi-ifo-coinc" );
    snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
    snprintf( this_proc_param->value, LIGOMETA_TYPE_MAX, " " );
  }

  /* store the distance cut in the process_params table */
  if ( distCut )
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *)
      calloc( 1, sizeof(ProcessParamsTable) );
    snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, 
        "%s", PROGRAM_NAME );
    snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, 
        "--h1-h2-distance-cut" );
    snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
    snprintf( this_proc_param->value, LIGOMETA_TYPE_MAX, " " );
  }

  /* store the h1h2 consistency option */ 
  if ( h1h2Consistency )
  {
    if ( !( snrCut > 0 ) )
    {
      fprintf( stderr,"The --snr-cut option must be specified and \n"
               "greater than zero when --h1-h2-consistency is given\n" );
      exit( 1 );
    }
    
    this_proc_param = this_proc_param->next = (ProcessParamsTable *)
      calloc( 1, sizeof(ProcessParamsTable) );
    snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, 
        "%s", PROGRAM_NAME );
    snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, 
        "--h1-h2-consistency" );
    snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
    snprintf( this_proc_param->value, LIGOMETA_TYPE_MAX, " " );
  }
  
  /* store the veto option */ 
  if ( doVeto )
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *)
      calloc( 1, sizeof(ProcessParamsTable) );
    snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, 
        "%s", PROGRAM_NAME );
    snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, 
        "--do-veto" );
    snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
    snprintf( this_proc_param->value, LIGOMETA_TYPE_MAX, " " );    
  }
  
  /* store the H1H2 snr cut for BCVSpin */
  if (doBCV2H1H2Veto)
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *)
      calloc( 1, sizeof(ProcessParamsTable) );
    snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, 
        "%s", PROGRAM_NAME );
    snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, 
        "--do-bcvspin-h1h2-veto");
    snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
    snprintf( this_proc_param->value, LIGOMETA_TYPE_MAX, " " );    
  }
  

  /* store the complete-coincs option */
  if (completeCoincs)
  {
    if (!multiIfoCoinc)
    {
      fprintf( stderr,
          "--multi-ifo-coinc must be specified when --complete-coincs is.\n" );
      exit( 1 );
    }

    this_proc_param = this_proc_param->next = (ProcessParamsTable *)
      calloc( 1, sizeof(ProcessParamsTable) );
    snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, 
        "%s", PROGRAM_NAME );
    snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, 
        "--complete-coincs");
    snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
    snprintf( this_proc_param->value, LIGOMETA_TYPE_MAX, " " );    
  }


  /* delete the first, empty process_params entry */
  this_proc_param = processParamsTable.processParamsTable;
  processParamsTable.processParamsTable = 
    processParamsTable.processParamsTable->next;
  free( this_proc_param );

  /*
   *
   * read in the input data from the rest of the arguments
   *
   */


  if ( optind < argc )
  {
    for( i = optind; i < argc; ++i )
    {
      INT4 numFileTriggers = 0;
      SnglInspiralTable   *inspiralFileList = NULL;
      SnglInspiralTable   *thisFileTrigger  = NULL;

      numFileTriggers = XLALReadInspiralTriggerFile( &inspiralFileList,
          &thisFileTrigger, &searchSummList, &inputFiles, argv[i] );

      if (numFileTriggers < 0)
      {
        fprintf(stderr, "Error reading triggers from file %s",
            argv[i]);
        exit( 1 );
      }
  
      /* read the summ value table as well. 
       * We do not need to do any sanity check.
       * It has been done in the previous call */
      XLALReadSummValueFile(&summValueList, argv[i]);
        
       
      /* maximize over a given interval */
      if ( maximizationInterval )
      {
        if (vrbflg)
        {
          fprintf( stdout, "Clustering triggers for over %d ms window\n",
              maximizationInterval);
        }
        XLALMaxSnglInspiralOverIntervals( &inspiralFileList, 
            (1.0e6 * maximizationInterval) );
        numFileTriggers = XLALCountSnglInspiral( inspiralFileList );
      }

      numTriggers += numFileTriggers;
      /* If there are any remaining triggers ... */
      if ( inspiralFileList )
      {
        /* add inspirals to list */
        if ( thisInspiralTrigger )
        {
          thisInspiralTrigger->next = inspiralFileList;
        }
        else
        {
          inspiralEventList = thisInspiralTrigger = inspiralFileList;
        }
        for( ; thisInspiralTrigger->next; 
            thisInspiralTrigger = thisInspiralTrigger->next);
      }

    }
    
  }
  else
  {
    fprintf( stderr, "Error: No trigger files specified.\n" );
    exit( 1 );
  }

  if ( vrbflg ) fprintf( stdout, "Read in a total of %d triggers.\n",
      numTriggers );

  /* We apply the SNR cut if requested */
  if ( snrCut > 0.0 )
  {    
    if ( vrbflg )
        fprintf(stdout, "Removing triggers with SNR lower than %f : ", snrCut); 
    LAL_CALL( LALSNRCutSingleInspiral( &status, &(inspiralEventList), 
          snrCut),  &status );
    numTriggers = XLALCountSnglInspiral(inspiralEventList);
    if ( vrbflg ) fprintf( stdout, "Have %d remaining triggers.\n",
          numTriggers );
  }
  
  /*    we initialise the veto segment list needed either by the h1h2 
        consistency check or the veto option itself. */
  if (h1h2Consistency || doVeto)
  {
    for ( ifoNumber = 0; ifoNumber< LAL_NUM_IFO; ifoNumber++)
    {
      /* to perform the veto, we  need the filename. */
      if ( vetoFileName[ifoNumber] && haveTrig[ifoNumber])
      {
        XLALSegListInit( &(vetoSegs[ifoNumber]) );
        LAL_CALL( LALSegListRead( &status, &(vetoSegs[ifoNumber]), 
              vetoFileName[ifoNumber], NULL),&status);
        XLALSegListCoalesce( &(vetoSegs[ifoNumber]) );

        /* keep only the segments that lie within the data-segment part */
        XLALSegListKeep(  &(vetoSegs[ifoNumber]), &startCoinc, &endCoinc );

        /* if the veto option is set, we remove single inspiral triggers 
           inside the list provided but we need to loop over the different 
           ifo name. */
        if (doVeto)
        {
          XLALReturnIFO(ifo,ifoNumber);
          if ( vrbflg ) fprintf( stdout, 
              "Applying veto segment (%s) list on ifo  %s \n ",
              vetoFileName[ifoNumber], ifo );
          inspiralEventList = XLALVetoSingleInspiral( inspiralEventList, 
              &(vetoSegs[ifoNumber]), ifo );
         /* count the triggers  */
          numTriggers = XLALCountSnglInspiral( inspiralEventList );
          if ( vrbflg ) fprintf( stdout, 
              " --> %d remaining triggers after veto segment list applied.\n",
              numTriggers );
        }
      }
    }
  }


  /*
   * for the case of BCV unconstrained-max, discard the triggers that
   * have alphaF greater than alphaFhi or lower than alphaFlo,
   * as those are specified in the command line.
   * If at least alphaFhi is not specified, do not discard any triggers.
   */

 
  /* perform the alphaf-cut */
  if (doBCVC & doAlphaFCut)
    {
      snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, 
         "--do-alphaf-cut" );
    if ( vrbflg ) fprintf( stdout,
       "Discarding triggers in H1 with alphaF > %f OR alphaF < %f (BCVC case) \n", 
        alphafParams.h1_lo, alphafParams.h1_hi );
    if ( vrbflg ) fprintf( stdout,
       "Discarding triggers in H2 with alphaF > %f OR alphaF < %f (BCVC case) \n", 
        alphafParams.h2_lo, alphafParams.h2_hi );
    if ( vrbflg ) fprintf( stdout,
       "Discarding triggers in L1 with alphaF > %f OR alphaF < %f (BCVC case) \n", 
        alphafParams.l1_lo, alphafParams.l1_hi );

     LAL_CALL( LALBCVCVetoSingleInspiral( &status, &(inspiralEventList),
        alphafParams ), &status );

     numTriggers = XLALCountSnglInspiral( inspiralEventList );
     if ( vrbflg ) fprintf( stdout, "%d remaining triggers after alphaF cut.\n",
         numTriggers );
    }
  
  /* check that we have read in data for all the requested times
     in all the requested instruments */
  if ( checkTimes )
  {
    if ( vrbflg ) fprintf( stdout, 
        "Checking that we have data for all times from all IFOs\n");
    for ( ifoNumber = 0; ifoNumber < (int)numIFO; ++ifoNumber )
    {
      LAL_CALL( LALCheckOutTimeFromSearchSummary ( &status, searchSummList, 
            ifoName[ifoNumber], &startCoinc, &endCoinc ), &status);
    }
  }

  if ( ! inspiralEventList )
  {
    /* no triggers, so no coincidences can be found */
    if ( vrbflg ) fprintf( stdout,
        "No triggers read in so no coincidences can be found\n" );

    goto cleanexit;
  }


  /* time sort the triggers */
  if ( vrbflg ) fprintf( stdout, "Sorting triggers\n" );
  LAL_CALL( LALSortSnglInspiral( &status, &(inspiralEventList),
        LALCompareSnglInspiralByTime ), &status );

  /* keep only triggers within the requested interval */
  if ( vrbflg ) fprintf( stdout, 
      "Discarding triggers outside requested interval\n" );
  LAL_CALL( LALTimeCutSingleInspiral( &status, &inspiralEventList,
        &startCoinc, &endCoinc), &status );


  /* keep play/non-play/all triggers */
  if ( dataType == playground_only && vrbflg ) fprintf( stdout, 
      "Keeping only playground triggers\n" );
  else if ( dataType == exclude_play && vrbflg ) fprintf( stdout, 
      "Keeping only non-playground triggers\n" );
  else if ( dataType == all_data && vrbflg ) fprintf( stdout, 
      "Keeping all triggers\n" );
  LAL_CALL( LALPlayTestSingleInspiral( &status, &inspiralEventList,
        &dataType ), &status );

  /* scroll to the end of the linked list of triggers, counting triggers and
   * freeing event_id's that were read in */
  thisInspiralTrigger = inspiralEventList;
  for (numTriggers = 0 ; thisInspiralTrigger; ++numTriggers,
      thisInspiralTrigger = thisInspiralTrigger->next )
  {  
    while ( thisInspiralTrigger->event_id )
    {
      eventId = (thisInspiralTrigger)->event_id;
      (thisInspiralTrigger)->event_id = (thisInspiralTrigger)->event_id->next;
      LALFree( eventId );
    }
  }

  if ( vrbflg ) fprintf( stdout, 
      "%d remaining triggers after time and data type cut.\n", numTriggers );


  if ( ! inspiralEventList )
  {
    /* no triggers remaining, so no coincidences can be found */
    if ( vrbflg ) fprintf( stdout,
        "No triggers remain after time/playground cuts\n"
        "No coincidences can be found\n" );

    goto cleanexit;
  }

  /* count the number of triggers for each IFO */
  for ( ifoNumber = 0; ifoNumber < LAL_NUM_IFO ; ++ifoNumber )
  {
    LALIfoCountSingleInspiral(&status, &numTrigs[ifoNumber], 
        inspiralEventList, ifoNumber);
    XLALReturnIFO(ifo,ifoNumber);
    if ( vrbflg ) fprintf( stdout, 
        "Have %d triggers from %s.\n", numTrigs[ifoNumber], 
        ifo );
    if ( numTrigs[ifoNumber] && !haveTrig[ifoNumber] )
    {
      fprintf( stderr, "Read in triggers from %s, none expected.\n",
          ifo);
      exit( 1 );
    }
    if ( haveTrig[ifoNumber] && numTrigs[ifoNumber] )
    {
      ++numTrigIFO;
    }
  }

  if ( !numTrigIFO )
  {
    if ( vrbflg ) fprintf( stdout, "Have no triggers from any IFOs\n"
        "Cannot be coincidences, so exiting without looking.\n");
    goto cleanexit;
  }
  else if ( numTrigIFO==1 )
  {
    if ( vrbflg ) fprintf( stdout, "Have triggers from only one IFO\n"
        "Cannot be coincidences, so exiting without looking.\n");
    goto cleanexit;
  }

  /* Populate the lightTravel matrix */
  if ( accuracyParams.exttrig )
  {
    LIGOTimeGPS timeTrigger;

    /* read the extTriggersTable from a file */
    numExtTriggers=LALExtTriggerTableFromLIGOLw( &exttrigHead, sourceFile,
                                                 0, 1);
    printf("Number of triggers read from the external trigger file: %d\n",
           numExtTriggers);
    
    if (numExtTriggers>1)
    {
      printf("WARNING: Only 1 external trigger expected in the file '%s'",
             sourceFile );
    }
    if (numExtTriggers==0)
    {
      printf("ERROR: No external trigger found in file '%s'",sourceFile );
      exit(1);
    } 

    /* extract the exttrig-time */
    timeTrigger.gpsSeconds     = exttrigHead->start_time;
    timeTrigger.gpsNanoSeconds = exttrigHead->start_time_ns;
    
    /* populate the accuracy params table */
    XLALPopulateAccuracyParamsExt( &accuracyParams, 
                                   &timeTrigger, exttrigHead->event_ra, 
                                   exttrigHead->event_dec );
  }
  else 
  {
    XLALPopulateAccuracyParams( &accuracyParams);
 
  }
 
  
  /* 
   *  
   * check for two IFO coincidence
   *
   */

 
  /* perform the time slides */
  for( slideNum = -numSlides; slideNum <= numSlides; slideNum++ )
  {
    SnglInspiralTable    *slideOutput = NULL;
    INT4                  numCoincInSlide = 0;

    coincInspiralList = NULL; 
      
    if ( numSlides ) 
    {
      for( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++ )
      {
        if( slideNum == -numSlides )
        {
          /* Initialize the slide-times with the initial value, 
             which is the negative extreme value */
          REAL8 tmpInt;
          INT4 tmpSlide, tmpSlideNS;
          tmpSlideNS = (INT4) (-1000000000*
              modf(numSlides * slideStep[ifoNumber], &tmpInt) );
          tmpSlide = (INT4) (-tmpInt);
          slideTimes[ifoNumber].gpsSeconds = tmpSlide;
          slideTimes[ifoNumber].gpsNanoSeconds = tmpSlideNS;
          slideReset[ifoNumber].gpsSeconds = (-tmpSlide);
          slideReset[ifoNumber].gpsNanoSeconds = (-tmpSlideNS);
        }
        else
        {
          /* Set the slide-times to the constant slide-step, 
             since this times refers to 'inspiralEventList', 
             which triggers are shifted in each slide by a constant amount. 
             The reset-time however refers to the coincidence list, so it must
             be decreased for each slide. */
          REAL8 tmpInt;
          INT4 tmpSlide, tmpSlideNS;
          tmpSlideNS = (INT4) (-1000000000*
              modf(slideStep[ifoNumber], &tmpInt) );
          tmpSlide = (INT4) tmpInt;

          slideTimes[ifoNumber].gpsSeconds = tmpSlide;
          slideTimes[ifoNumber].gpsNanoSeconds = tmpSlideNS;
          slideReset[ifoNumber].gpsSeconds -= tmpSlide;
          slideReset[ifoNumber].gpsNanoSeconds -= tmpSlideNS;
        }
      }
    }
    
    if ( vrbflg ) fprintf(stdout, "Performing time slide %d\n", slideNum );
    
    /* slide the data for any case */
    XLALTimeSlideSingleInspiral( inspiralEventList, &startCoinc, &endCoinc,
                                 slideTimes);
    LAL_CALL( LALSortSnglInspiral( &status, &(inspiralEventList),
              LALCompareSnglInspiralByTime ), &status );
    
    for ( ifoNumber = 0; ifoNumber< LAL_NUM_IFO; ifoNumber++) 
    {
      /* FIXME:  this code is executed even if there are no veto segments
       * (!?), so this function call fails.  It failed in the old code,
       * too, but there was no error checking so the code would keep
       * running.  Someone else will have to figure out how to untangle
       * this, for now I turn off error checking here too and clear the
       * error flag (this reproduces the old behaviour, but clearly
       * something else should be done) */
      if ( XLALTimeSlideSegList( &vetoSegs[ifoNumber], &startCoinc, &endCoinc,
                                 &slideTimes[ifoNumber] ) < 0 )
        /*exit(1)*/ XLALClearErrno();
    }

    /* don't analyze zero-lag if numSlides>0 */
    if ( numSlides && !slideNum ) continue;

    /* look for coincidences */    
    if ( accuracyParams.test == ellipsoid )
    {
      LAL_CALL( LALCreateTwoIFOCoincListEllipsoid(&status, &coincInspiralList,
              inspiralEventList, &accuracyParams ), &status);
    }
    else
    {
      LAL_CALL( LALCreateTwoIFOCoincList(&status, &coincInspiralList,
              inspiralEventList, &accuracyParams ), &status);
    }


    /* count the zero-lag coincidences */
    if ( !numSlides) 
    {
      if( coincInspiralList )
      {
        for (numCoinc = 0, thisCoinc = coincInspiralList;
               thisCoinc; ++numCoinc, thisCoinc = thisCoinc->next )
        {
        }
        if ( vrbflg ) fprintf( stdout, 
            "%d coincident triggers found before coincidence cuts..\n",
            numCoinc);
      }
    }
        
    /* BNS and Machos case */
    if( distCut && !doBCVC && haveTrig[LAL_IFO_H1] && haveTrig[LAL_IFO_H2] )
    /*S3*/
    {
      if ( vrbflg ) fprintf( stdout, 
         "Discarding triggers h1-h2-distance-cut using kappa=%f and epsilon = %f\n",                       
         accuracyParams.ifoAccuracy[LAL_IFO_H1].kappa,
         accuracyParams.ifoAccuracy[LAL_IFO_H1].epsilon);
       
      XLALInspiralDistanceCut( &coincInspiralList, &accuracyParams);
    }
    
    /* BCV case */  
    if (doBCVC)
    {     
      snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, 
         "--bcvc" );
      /* perform the iota cut */
      if( iotaCut  ) 
      {
        if ( vrbflg ) 
          fprintf( stdout, 
              "Applying iota Cut with iotaCutH1H2=%f and iotaCutH1L1=%f \n",
              accuracyParams.iotaCutH1H2,accuracyParams.iotaCutH1L1 ); 
            
        XLALInspiralIotaCutBCVC( &coincInspiralList, &accuracyParams );
        if ( vrbflg )  fprintf( stdout, "%d remaining coincident triggers.\n", 
            XLALCountCoincInspiral(coincInspiralList));
      }
    }
    
    if ( doPsi0Psi3Cut ) 
    {
      snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, 
         "--psi0-psi3-cut" );
      if ( vrbflg ) fprintf( stdout,
         "Discarding triggers using Dpsi0Dpsi3 cut \n");
      XLALInspiralPsi0Psi3CutBCVC( &coincInspiralList );
      if ( vrbflg ) fprintf( stdout, "%d remaining coincident triggers .\n", 
          XLALCountCoincInspiral(coincInspiralList));
    }
    
    
    /* BCVSpin case */
    if (doBCV2H1H2Veto)
    {
      if (vrbflg) fprintf (stdout, 
          "Discarding triggers with H2 snr > H1 snr \n" );
      XLALInspiralSNRCutBCV2( &coincInspiralList);
    }
       
    
    if ( multiIfoCoinc )
    {
      for( N = 3; N <= numIFO; N++)
      {
        LAL_CALL( LALCreateNIFOCoincList( &status, &coincInspiralList, 
              &accuracyParams, N ), &status );
      }
      LAL_CALL( LALRemoveRepeatedCoincs( &status, &coincInspiralList ), 
          &status );
    }
  
    /* perform the h1h2-consistency check */
    if ( h1h2Consistency && haveTrig[LAL_IFO_H1] && haveTrig[LAL_IFO_H2] )
    {
      if(vrbflg) 
      {
        if (vetoFileName[LAL_IFO_H1] && vetoFileName[LAL_IFO_H2])
        {
          fprintf(stdout, 
              "Using h1-h2-consistency with veto segment list %s and %s\n", 
              vetoFileName[LAL_IFO_H1], vetoFileName[LAL_IFO_H2]);
        }
        else
        { 
          fprintf(stdout, 
              "Using h1-h2-consistency without veto segment list. NOT RECOMMENDED\n");
        }
      }
      LAL_CALL( LALInspiralDistanceCutCleaning(&status,  &coincInspiralList,
           &accuracyParams, snrCut, &summValueList, &vetoSegs[LAL_IFO_H1], 
           &vetoSegs[LAL_IFO_H2]), &status);
      if ( vrbflg ) fprintf( stdout, 
          "%d remaining coincident triggers after h1-h2-consisteny .\n", 
          XLALCountCoincInspiral(coincInspiralList));
    }

    /* no time-slide */   
    if ( !slideNum ) 
    {
      /* count the coincs */
      if( coincInspiralList )
      {
        for (numCoinc = 0, thisCoinc = coincInspiralList;
            thisCoinc; ++numCoinc, thisCoinc = thisCoinc->next )
        {
          if ( thisCoinc->numIfos == 2 )
          {
            ++numDoubles;
          }
          else if ( thisCoinc->numIfos == 3 )
          {
            ++numTriples;
          }
          else if ( thisCoinc->numIfos == 4 )
          {
            ++numQuadruples;
          }
        }    
      }

      if ( vrbflg ) 
      {
        fprintf( stdout, "%d coincident triggers found.\n", numCoinc );
        fprintf( stdout, "%d double coincident triggers\n"
            "%d triple coincident triggers\n"
            "%d quadruple coincident triggers\n",
            numDoubles, numTriples, numQuadruples );
      }
    }          


    if ( coincInspiralList ) 
    {
      /* count the coincs, scroll to end of list */
      if( slideNum )
      {  

        /* keep only the requested coincs */
        if( slideH1H2Together )
        {
          if ( vrbflg ) fprintf( stdout,
              "Throwing out slide coincs found only as H1H2 doubles.\n" );
          numCoincInSlide = XLALCoincInspiralIfosDiscard(
              &coincInspiralList, "H1H2" );
          if ( vrbflg ) fprintf( stdout,
              "Kept %d non-H1H2 coincs in slide.\n", numCoincInSlide );
        }

        numCoincInSlide = 0;
        thisCoinc = coincInspiralList;
        while( thisCoinc )
        {
          ++numCoincInSlide;
          thisCoinc = thisCoinc->next;
        }
        
        if ( vrbflg ) fprintf( stdout,
            "%d coincident triggers found in slide.\n", numCoincInSlide );
        numCoinc += numCoincInSlide;          
      }

      /* complete the coincs */
      if ( completeCoincs )
      {
        completionSngls = XLALCompleteCoincInspiral ( coincInspiralList, 
            haveTrig);

        /* do a veto on the new sngls */
        for ( ifoNumber = 0; ifoNumber< LAL_NUM_IFO; ifoNumber++)
        {
        
          if (doVeto && vetoFileName[ifoNumber] && haveTrig[ifoNumber])
          {
            XLALReturnIFO(ifo,ifoNumber);
            if ( vrbflg ) fprintf( stdout, 
                "Applying veto list %s on completion sngls for ifo %s \n",
                vetoFileName[ifoNumber], ifo );
            completionSngls = XLALVetoSingleInspiral( completionSngls,
                &(vetoSegs[ifoNumber]), ifo );
          }
        }
      }

      /* for all: write out all coincs as singles with event IDs */
      LAL_CALL( LALExtractSnglInspiralFromCoinc( &status, &slideOutput, 
            coincInspiralList, &startCoinc, slideNum), &status );

      if ( numSlides && coincInspiralList ) 
      {
 
        /* the output triggers should be slid back to original time */
        XLALTimeSlideSingleInspiral( slideOutput, &startCoinc, &endCoinc,
                                     slideReset);
        LAL_CALL( LALSortSnglInspiral( &status, &(slideOutput),
              LALCompareSnglInspiralByTime ), &status );
      }
      
      while ( coincInspiralList )
      {
        thisCoinc = coincInspiralList;
        coincInspiralList = coincInspiralList->next;
        XLALFreeCoincInspiral( &thisCoinc );
      }

      while ( completionSngls )
      {
        SnglInspiralTable *thisSngl = NULL;
        thisSngl = completionSngls;
        completionSngls = completionSngls->next;
        XLALFreeSnglInspiral( &thisSngl );
      }
    }

  
    if ( snglOutput )
    {
      thisInspiralTrigger->next = slideOutput;
    }
    else
    {
      snglOutput = slideOutput;
    }
    
    if ( numSlides )
    {
      /* scroll to the end of the list */
      if ( slideOutput )
      {
        for( thisInspiralTrigger = slideOutput; thisInspiralTrigger->next; 
            thisInspiralTrigger = thisInspiralTrigger->next);
      }
    }
  }


  /*
   *
   * write the output xml file
   *
   */


  /* since we don't yet write coinc inspiral tables, we must make a list of
   * sngl_inspiral tables with the eventId's appropriately poplulated */
  
   
cleanexit:

  searchsumm.searchSummaryTable->in_start_time = startCoinc;
  searchsumm.searchSummaryTable->in_end_time = endCoinc;
  searchsumm.searchSummaryTable->out_start_time = startCoinc;
  searchsumm.searchSummaryTable->out_end_time = endCoinc;
  searchsumm.searchSummaryTable->nnodes = 1;

  if ( vrbflg ) fprintf( stdout, "writing output file... " );

  if ( userTag && ifoTag && !outCompress )
  {
    snprintf( fileName, FILENAME_MAX, "%s-THINCA_%s_%s-%d-%d.xml", 
        ifos, ifoTag, userTag, startCoincidence, 
        endCoincidence - startCoincidence );
    snprintf( fileSlide, FILENAME_MAX, "%s-THINCA_SLIDE_%s_%s-%d-%d.xml", 
        ifos, ifoTag, userTag, startCoincidence, 
        endCoincidence - startCoincidence );
  }
  else if ( !userTag && ifoTag && !outCompress )
  {
    snprintf( fileName, FILENAME_MAX, "%s-THINCA_%s-%d-%d.xml", ifos,
        ifoTag, startCoincidence, endCoincidence - startCoincidence );
    snprintf( fileSlide, FILENAME_MAX, "%s-THINCA_SLIDE_%s-%d-%d.xml", ifos,
        ifoTag, startCoincidence, endCoincidence - startCoincidence );
  }
  else if ( userTag && !ifoTag && !outCompress )
  {
    snprintf( fileName, FILENAME_MAX, "%s-THINCA_%s-%d-%d.xml", 
        ifos, userTag, startCoincidence, endCoincidence - startCoincidence );
    snprintf( fileSlide, FILENAME_MAX, "%s-THINCA_SLIDE_%s-%d-%d.xml", 
        ifos, userTag, startCoincidence, endCoincidence - startCoincidence );
  }
  else if ( userTag && ifoTag && outCompress )
  {     snprintf( fileName, FILENAME_MAX, "%s-THINCA_%s_%s-%d-%d.xml.gz",
        ifos, ifoTag, userTag, startCoincidence,
        endCoincidence - startCoincidence );
    snprintf( fileSlide, FILENAME_MAX, "%s-THINCA_SLIDE_%s_%s-%d-%d.xml.gz",
        ifos, ifoTag, userTag, startCoincidence,
        endCoincidence - startCoincidence );
  }
  else if ( !userTag && ifoTag && outCompress )
  {
    snprintf( fileName, FILENAME_MAX, "%s-THINCA_%s-%d-%d.xml.gz", ifos,
        ifoTag, startCoincidence, endCoincidence - startCoincidence );
    snprintf( fileSlide, FILENAME_MAX, "%s-THINCA_SLIDE_%s-%d-%d.xml.gz",
        ifos, ifoTag, startCoincidence, endCoincidence - startCoincidence );
  }
  else if ( userTag && !ifoTag && outCompress )
  {
    snprintf( fileName, FILENAME_MAX, "%s-THINCA_%s-%d-%d.xml.gz",
        ifos, userTag, startCoincidence, endCoincidence - startCoincidence );
    snprintf( fileSlide, FILENAME_MAX, "%s-THINCA_SLIDE_%s-%d-%d.xml.gz",
        ifos, userTag, startCoincidence, endCoincidence - startCoincidence );
  }
  else if ( !userTag && !ifoTag && outCompress )
  {
    snprintf( fileName, FILENAME_MAX, "%s-THINCA-%d-%d.xml.gz", ifos,
        startCoincidence, endCoincidence - startCoincidence );
    snprintf( fileSlide, FILENAME_MAX, "%s-THINCA_SLIDE-%d-%d.xml.gz", ifos,
        startCoincidence, endCoincidence - startCoincidence );
  }
  else
  {
    snprintf( fileName, FILENAME_MAX, "%s-THINCA-%d-%d.xml", ifos,
        startCoincidence, endCoincidence - startCoincidence );
    snprintf( fileSlide, FILENAME_MAX, "%s-THINCA_SLIDE-%d-%d.xml", ifos,
        startCoincidence, endCoincidence - startCoincidence );
  }
  searchsumm.searchSummaryTable->nevents = numCoinc;

  memset( &xmlStream, 0, sizeof(LIGOLwXMLStream) );

  if ( !numSlides )
  {
    LAL_CALL( LALOpenLIGOLwXMLFile( &status , &xmlStream, fileName ), 
        &status );
  }
  else
  {
    LAL_CALL( LALOpenLIGOLwXMLFile( &status , &xmlStream, fileSlide ), 
        &status );
  }
  /* write process table */

  snprintf( proctable.processTable->ifos, LIGOMETA_IFOS_MAX, ifos );

  XLALGPSTimeNow(&(proctable.processTable->end_time));
  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlStream, process_table ), 
      &status );
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, proctable, 
        process_table ), &status );
  LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlStream ), &status );

  /* write process_params table */
  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlStream, 
        process_params_table ), &status );
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, processParamsTable, 
        process_params_table ), &status );
  LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlStream ), &status );

  /* write search_summary table */
  snprintf( searchsumm.searchSummaryTable->ifos, LIGOMETA_IFOS_MAX, ifos );

  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlStream, 
        search_summary_table ), &status );
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, searchsumm, 
        search_summary_table ), &status );
  LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlStream ), &status );

  /* write the search_summvars table */
  LAL_CALL( LALBeginLIGOLwXMLTable( &status ,&xmlStream, 
        search_summvars_table), &status );
  searchSummvarsTable.searchSummvarsTable = inputFiles;
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, searchSummvarsTable,
        search_summvars_table), &status );
  LAL_CALL( LALEndLIGOLwXMLTable( &status, &xmlStream), &status );

  
  /* write the summ_value table */
  if( summValueList)
  {
    LAL_CALL( LALBeginLIGOLwXMLTable( &status ,&xmlStream, 
          summ_value_table), &status );
    summValueTable.summValueTable = summValueList;
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, summValueTable,
          summ_value_table), &status );
    LAL_CALL( LALEndLIGOLwXMLTable( &status, &xmlStream), &status );
  }

  /* write the sngl_inspiral table */
  if( snglOutput )
  {
    LAL_CALL( LALBeginLIGOLwXMLTable( &status ,&xmlStream, 
          sngl_inspiral_table), &status );
    inspiralTable.snglInspiralTable = snglOutput;
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, inspiralTable,
          sngl_inspiral_table), &status );
    LAL_CALL( LALEndLIGOLwXMLTable( &status, &xmlStream), &status );
  }

  LAL_CALL( LALCloseLIGOLwXMLFile( &status, &xmlStream), &status );

  if ( vrbflg ) fprintf( stdout, "done\n" );


  /*
   *
   * clean up the memory that has been allocated 
   *
   */


  if ( vrbflg ) fprintf( stdout, "freeing memory... " );

  free( proctable.processTable );
  free( searchsumm.searchSummaryTable );

  while ( processParamsTable.processParamsTable )
  {
    this_proc_param = processParamsTable.processParamsTable;
    processParamsTable.processParamsTable = this_proc_param->next;
    free( this_proc_param );
  }

  while ( inputFiles )
  {
    thisInputFile = inputFiles;
    inputFiles = thisInputFile->next;
    LALFree( thisInputFile );
  }

  while ( searchSummList )
  {
    thisSearchSumm = searchSummList;
    searchSummList = searchSummList->next;
    LALFree( thisSearchSumm );
  }

  while ( summValueList )
  {
    thisSummValue = summValueList;
    summValueList = summValueList->next;
    LALFree( thisSummValue );
  } 

  /* free the veto segment list. */
  for (ifoNumber=0; ifoNumber<LAL_NUM_IFO; ifoNumber++)
  {
    if ( vetoFileName[ifoNumber]  && haveTrig[ifoNumber])
    {
      free( vetoFileName[ifoNumber] );
    }

   if (vetoSegs[ifoNumber].initMagic == SEGMENTSH_INITMAGICVAL )
   {
        XLALSegListClear( &vetoSegs[ifoNumber] );
    }
  }
  
/* free the snglInspirals */
  while ( inspiralEventList )
  {
    thisInspiralTrigger = inspiralEventList;
    inspiralEventList = inspiralEventList->next;
    LAL_CALL( LALFreeSnglInspiral( &status, &thisInspiralTrigger ), &status );
  }

  while ( snglOutput )
  {
    thisInspiralTrigger = snglOutput;
    snglOutput = snglOutput->next;
    LAL_CALL( LALFreeSnglInspiral( &status, &thisInspiralTrigger ), &status );
  }
 
  while ( coincInspiralList )
  {
    thisCoinc = coincInspiralList;
    coincInspiralList = coincInspiralList->next;
    XLALFreeCoincInspiral( &thisCoinc );
  }



  if ( userTag ) free( userTag );
  if ( ifoTag ) free( ifoTag );

  if ( vrbflg ) fprintf( stdout, "done\n" );

  LALCheckMemoryLeaks();

  exit( 0 );
}
