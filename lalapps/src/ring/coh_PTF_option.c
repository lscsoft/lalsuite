/*
*  Copyright (C) 2007 Duncan Brown, Jolien Creighton, Lisa M. Goggin, Matt Pitkin
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

#include "coh_PTF.h"

/* parse command line arguments using getopt_long to get ring params */
int coh_PTF_parse_options(struct coh_PTF_params *params,int argc,char **argv )
{

  CHAR                         ifo[LIGOMETA_IFO_MAX];
  UINT4                        ifoNumber;
  static struct coh_PTF_params localparams;
  memset( &localparams.haveTrig, 0, LAL_NUM_IFO * sizeof(int) );
  struct option                long_options[] =
  {
    { "verbose",            no_argument, &vrbflg, 1 },
    { "strain-data",        no_argument, &localparams.strainData, 1 },
    { "zero-data",          no_argument, &localparams.zeroData, 1 },
    { "theoretical-spectrum",     no_argument, &localparams.whiteSpectrum, 1 },
    { "write-raw-data",     no_argument, &localparams.writeRawData, 1 },
    { "write-data",         no_argument, &localparams.writeProcessedData, 1 },
    { "write-inv-spectrum", no_argument, &localparams.writeInvSpectrum, 1 },
    { "write-segment",      no_argument, &localparams.writeSegment, 1 },
    { "write-filter-output",no_argument, &localparams.writeFilterOutput, 1 },
    { "analyze-inj-segs-only",no_argument, &localparams.analyzeInjSegsOnly, 1 },
    { "do-null-stream",     no_argument, &localparams.doNullStream, 1 },
    { "do-trace-snr",       no_argument, &localparams.doTraceSNR, 1 },
    { "do-bank-veto",       no_argument, &localparams.doBankVeto, 1 },
    { "do-auto-veto",       no_argument, &localparams.doAutoVeto, 1 },
    { "do-chi-square",      no_argument, &localparams.doChiSquare, 1 },
    { "do-sngl-chi-tests",  no_argument, &localparams.doSnglChiSquared, 1},
    { "do-clustering",      no_argument, &localparams.clusterFlag, 1},
/*    {"g1-data",             no_argument, &(haveTrig[LAL_IFO_G1]), 1 },*/
    {"h1-data",             no_argument, &(localparams.haveTrig[LAL_IFO_H1]),1},
    {"h2-data",             no_argument, &(localparams.haveTrig[LAL_IFO_H2]),1},
    {"l1-data",             no_argument, &(localparams.haveTrig[LAL_IFO_L1]),1},
/*    {"t1-data",             no_argument, &(haveTrig[LAL_IFO_T1]), 1 },*/
    {"v1-data",             no_argument, &(localparams.haveTrig[LAL_IFO_V1]),1},
    {"face-on-analysis",    no_argument, &(localparams.faceOnAnalysis),1},
    {"face-away-analysis",    no_argument, &(localparams.faceAwayAnalysis),1},
    {"dynamic-template-length",no_argument, &(localparams.dynTempLength),1},
    {"store-amplitude-params",no_argument, &(localparams.storeAmpParams),1},
    {"analyse-segment-end", no_argument, &(localparams.analSegmentEnd),1},
    {"do-short-slides", no_argument, &(localparams.doShortSlides),1},
    { "write-sngl-inspiral-table", no_argument, &(localparams.writeSnglInspiralTable),1},
    { "help",               no_argument, 0, 'h' },
    { "version",            no_argument, 0, 'V' },
    { "simulated-data",          required_argument, 0, '6' },
    { "gps-start-time",          required_argument, 0, 'a' },
    { "gps-start-time-ns",       required_argument, 0, 'A' },
    { "gps-end-time",            required_argument, 0, 'b' },
    { "gps-end-time-ns",         required_argument, 0, 'B' },
    { "trigger-time",            required_argument, 0, '<' },
    { "trigger-time-ns",         required_argument, 0, '>' },
    { "h1-channel-name",         required_argument, 0, 'c' },
    { "h1-frame-cache",          required_argument, 0, 'D' },
    { "h2-channel-name",         required_argument, 0, 'x' },
    { "h2-frame-cache",          required_argument, 0, 'X' },
    { "l1-channel-name",         required_argument, 0, 'y' },
    { "l1-frame-cache",          required_argument, 0, 'Y' },
    { "v1-channel-name",         required_argument, 0, 'z' },
    { "v1-frame-cache",          required_argument, 0, 'Z' },
    { "low-template-freq",       required_argument, 0, 'e' },
    { "low-filter-freq",         required_argument, 0, 'H' },
    { "high-filter-freq",        required_argument, 0, 'I' },
    { "highpass-frequency",      required_argument, 0, 'E' },
    { "injection-file",          required_argument, 0, 'i' },
    { "snr-threshold",           required_argument, 0, 'j' },
    { "spin-snr-threshold",      required_argument, 0, '2' },
    { "sngl-snr-threshold",      required_argument, 0, '1' },
    { "trig-time-window",        required_argument, 0, 'J' },
    { "user-tag",                required_argument, 0, 'k' },
    { "ifo-tag",                 required_argument, 0, 'K' },
    { "non-spin-snr2-threshold", required_argument, 0, 'l' },
    { "spin-snr2-threshold",     required_argument, 0, 'L' },
    { "spin-bank",               required_argument, 0, 'm' },
    { "non-spin-bank",           required_argument, 0, 'M' },
    { "only-segment-numbers",    required_argument, 0, 'n' },
    { "only-template-numbers",   required_argument, 0, 'N' },
    { "output-file",             required_argument, 0, 'o' },
    { "bank-file",               required_argument, 0, 'O' },
    { "num-auto-chisq-points",   required_argument, 0, 'p' },
    { "auto-veto-time-step",     required_argument, 0, 'P' },
    { "num-chi-square-bins",     required_argument, 0, 'q' },
    { "chi-square-threshold",    required_argument, 0, 'Q' },
    { "random-seed",             required_argument, 0, 'r' },
    { "dynamic-range-factor",    required_argument, 0, 'R' },
    { "sample-rate",             required_argument, 0, 's' },
    { "segment-duration",        required_argument, 0, 'S' },
    { "psd-segment-duration",        required_argument, 0, '9' },
    { "bank-veto-templates",     required_argument, 0, 't' },
    { "inverse-spec-length",     required_argument, 0, 'T' },
    { "trig-start-time",         required_argument, 0, 'u' },
    { "trig-end-time",           required_argument, 0, 'U' },
    { "block-duration",          required_argument, 0, 'w' },
    { "pad-data",                required_argument, 0, 'W' },
    { "right-ascension",         required_argument, 0, 'f' },
    { "declination",             required_argument, 0, 'F' },
    { "sky-error",               required_argument, 0, 'g' },
    { "timing-accuracy",         required_argument, 0, 'G' },
    { "approximant",             required_argument, 0, 'C' },
    { "order",                   required_argument, 0, 'v' },
    { "h1-slide-segment",        required_argument, 0, '!' }, 
    { "h2-slide-segment",        required_argument, 0, '&' },
    { "l1-slide-segment",        required_argument, 0, '(' },
    { "v1-slide-segment",        required_argument, 0, ')' },
    { "sky-positions-file",      required_argument, 0, '#' },
    { "fft-level",               required_argument, 0, '|' },
    { "cluster-window",          required_argument, 0, '4' },
    { "inj-search-window",       required_argument, 0, '3' },
    { "inj-mchirp-window",       required_argument, 0, '5' },
    { "ligo-calibrated-data",    required_argument, 0, '7' }, 
    { "virgo-calibrated-data",   required_argument, 0, '8' }, 
    { "short-slide-offset",      required_argument, 0, '@' },
    { 0, 0, 0, 0 }
  };
  char args[] = "a:A:b:B:c:C:D:e:E:f:F:g:G:h:H:i:I:j:J:k:K:l:L:m:M:n:N:o:O:p:P:q:Q:r:R:s:S:t:T:u:U:v:V:w:W:x:X:y:Y:z:Z:1:2:3:4:5:6:7:8:9:<:>:!:&:(:):#:|:@";
  char *program = argv[0];

  /* set default values for parameters before parsing arguments */
  coh_PTF_default_params( &localparams );

  while ( 1 )
  {
    int option_index = 0;
    int c;

    c = getopt_long_only( argc, argv, args, long_options, &option_index );
    if ( c == -1 ) /* end of options */
      break;

    switch ( c )
    {
      case 0: /* if option set a flag, nothing else to do */
        if ( long_options[option_index].flag )
          break;
        else
          error( "error parsing option %s with argument %s\n",
              long_options[option_index].name, optarg );
      case 'a': /* gps-start-time */
        localparams.startTime.gpsSeconds = atol( optarg );
        break;
      case 'A': /* gps-start-time-ns */
        localparams.startTime.gpsNanoSeconds = atol( optarg );
        break;
      case 'b': /* gps-end-time */
        localparams.endTime.gpsSeconds = atol( optarg );
        break;
      case 'B': /* gps-end-time-ns */
        localparams.endTime.gpsNanoSeconds = atol( optarg );
        break;
      case '<': /* trigger-time */
        localparams.trigTime.gpsSeconds = atol( optarg );
        break;
      case '>': /* trigger-time-ns */ 
        localparams.trigTime.gpsNanoSeconds = atol( optarg );
        break;
      case 'c': /* h1 channel-name */
        localparams.channel[LAL_IFO_H1] = optarg;
        break;
      case 'D': /* h1 frame-cache */
        localparams.dataCache[LAL_IFO_H1] = optarg;
        break;
      case 'y': /* l1 channel-name */
        localparams.channel[LAL_IFO_L1] = optarg;
        break;
      case 'Y': /* l1 frame-cache */
        localparams.dataCache[LAL_IFO_L1] = optarg;
        break;
      case 'z': /* v1 channel-name */
        localparams.channel[LAL_IFO_V1] = optarg;
        break;
      case 'Z': /* v1 frame-cache */
        localparams.dataCache[LAL_IFO_V1] = optarg;
        break;
      case 'x': /* h2 channel-name */
        localparams.channel[LAL_IFO_H2] = optarg;
        break;
      case 'X': /* h2 frame-cache */
        localparams.dataCache[LAL_IFO_H2] = optarg;
        break;
      case 'e': /* start frequency of template generation */
        localparams.lowTemplateFrequency = atof( optarg );
        break;
      case 'H': /* start frequency of matched filter */
        localparams.lowFilterFrequency = atof( optarg );
        break;
      case 'I': /* End frequency of matched filter */
        localparams.highFilterFrequency = atof( optarg );
        break;
      case 'E': /* highpass-frequency */
        localparams.highpassFrequency = atof( optarg );
        break;
      case 'C': /* waveform approximant */
        if ( ! strcmp( "FindChirpSP", optarg ) )
        {
          localparams.approximant = FindChirpSP;
        }
        else if ( ! strcmp( "FindChirpPTF", optarg ) )
        {
          localparams.approximant = FindChirpPTF;
        }
        else if ( ! strcmp( "TaylorT1", optarg) )
        {
          localparams.approximant = TaylorT1;
        }
        else if ( ! strcmp( "TaylorT2", optarg) )
        {
          localparams.approximant = TaylorT2;
        }
        else if ( ! strcmp( "TaylorT3", optarg) )
        {
          localparams.approximant = TaylorT3;
        }
        else if ( ! strcmp( "TaylorT4", optarg) )
        {
          localparams.approximant = TaylorT4;
        }
        else if ( ! strcmp( "GeneratePPN", optarg) )
        {
          localparams.approximant = GeneratePPN;
        }
        else if ( ! strcmp( "PadeT1", optarg) )
        {
          localparams.approximant = PadeT1;
        }
        else if ( ! strcmp( "EOB", optarg) )
        {
          localparams.approximant = EOB;
        }
        else if ( ! strcmp( "EOBNR", optarg) )
        {
          localparams.approximant = EOBNR;
        }
        else if ( ! strcmp( "IMRPhenomB", optarg) )
        {
          localparams.approximant = IMRPhenomB;
        }
        else
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "unknown order specified: "
              "%s (must be either FindChirpSP, FindChirpPTF or TaylorT4)\n",
              long_options[option_index].name, optarg );
          exit( 1 );
        }
        break;
      case '6': /* Simulated data option */
        localparams.simData = 1;
        if ( ! strcmp( "WhiteNoise",optarg))
        {
          localparams.simDataType = WHITE_PSD;
        }
        else if ( ! strcmp( "ILIGONoise",optarg))
        {
          localparams.simDataType = ILIGO_PSD;
        }
        else if ( ! strcmp( "ALIGONoise",optarg))
        {
          localparams.simDataType = ALIGO_PSD;
        }
        else
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "unknown data type specified:"
              "%s valid options are WhiteNoise, ILIGONoise or ALIGONoise",
              long_options[option_index].name, optarg );
          exit(1);
        }
        break;
      case 'v': /* PN order of waveform */        
        if ( ! strcmp( "twoPN", optarg ) )
        {
          localparams.order = LAL_PNORDER_TWO;
        }
        else if ( ! strcmp( "twoPointFivePN", optarg ) )
        {
          localparams.order = LAL_PNORDER_TWO_POINT_FIVE;
        }
        else if ( ! strcmp( "threePN", optarg ) )
        {
          localparams.order = LAL_PNORDER_THREE;
        }
        else if ( ! strcmp( "threePointFivePN", optarg ) )
        {
          localparams.order = LAL_PNORDER_THREE_POINT_FIVE;
        }
        else if ( ! strcmp( "pseudoFourPN", optarg ) )
        {
          localparams.order = LAL_PNORDER_PSEUDO_FOUR;
        }
        else
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "unknown order specified: "
              "%s (must be one of twoPN, twoPointFivePN, threePN, threePointFivePN, pseudoFourPN)\n",
              long_options[option_index].name, optarg );
          exit( 1 );
        }
        break;
      case 'f': /* right-ascension */
        localparams.rightAscension = atof( optarg ) * LAL_PI_180;
        break;
      case 'F': /* Declination */
        localparams.declination = atof( optarg ) * LAL_PI_180;
        break;
      case 'g': /* Error in declination */
        localparams.skyError = atof( optarg ) * LAL_PI_180;
        break;
      case 'G': /* timing accuracy of network */
        localparams.timingAccuracy = atof( optarg );
        break;
      case 'h': /* help */
        coh_PTF_usage( program );
        exit( 0 );
      case 'i': /* injection-file */
        localparams.injectFile = optarg;
        break;
      case 'j':
        localparams.threshold = atof(optarg); 
        break;
      case '2':
        localparams.spinThreshold = atof(optarg);
        break;
      case '1':
        localparams.snglSNRThreshold = atof(optarg);
        break;
      case 'J':
        localparams.timeWindow = atof(optarg);
        break;
      case 'k': /* user-tag */
        strncpy( localparams.userTag, optarg, sizeof( localparams.userTag ) - 1 );
        break;
      case 'K': /* ifo-tag */
        strncpy( localparams.ifoTag, optarg, sizeof( localparams.ifoTag ) - 1 );
        break;
      case 'l':
        localparams.nonspinSNR2threshold = atof(optarg);
        break;
      case 'L':
        localparams.spinSNR2threshold = atof(optarg);
        break;
      case 'm': /* spin bank */
        localparams.spinBank = 1;
        strncpy( localparams.spinBankName, optarg, sizeof( localparams.spinBankName ) - 1 );
        break;
      case 'M': /* non spin bank */
        localparams.noSpinBank = 1;
        strncpy( localparams.noSpinBankName, optarg, sizeof( localparams.noSpinBankName ) - 1 );
        break;
      case 'n': /* only-segment-numbers */
        localparams.segmentsToDoList = optarg;
        break;
      case 'N': /* only-template-number */
        localparams.templatesToDoList = optarg;
        break;
      case 'o': /* output-file */
        strncpy( localparams.outputFile, optarg, sizeof( localparams.outputFile ) - 1 );
        break;
      case 'O': /* bank-file */
        localparams.bankFile = optarg;
        break;
      case 'p': /* num auto chisq points */
        localparams.numAutoPoints = atoi( optarg );
        break;
      case 'P': /* Auto veto time step */
        localparams.autoVetoTimeStep = atof( optarg );
        break;
      case 'q': /* num chi square bins */
        localparams.numChiSquareBins = atoi( optarg );
        break;
      case 'Q': 
        localparams.chiSquareCalcThreshold = atof( optarg );
        break;
      case 'r': /* random seed */
        localparams.randomSeed = atoi( optarg );
        break;
      case 'R': /* dynamic range factor */
        localparams.dynRangeFac = atof( optarg );
        break;
      case 's': /* sample rate */
        localparams.sampleRate = atof( optarg );
        break;
      case 'S': /* segment-duration */
        localparams.segmentDuration = atof( optarg );
        break;
      case '9': /* PSD segment-duration */
        localparams.psdSegmentDuration = atof( optarg );
        break;
      case 't': /* bank veto template bank */
        localparams.bankVetoBankName = optarg;
        break;
      case 'T': /* inverse-spec-length */
        localparams.truncateDuration = atof( optarg );
        break;
      case 'u': /* trig-start-time */
        localparams.trigStartTimeNS = (INT8) atol( optarg ) * LAL_INT8_C(1000000000);
        break;
      case 'U': /* trig-end-time */
        localparams.trigEndTimeNS = (INT8) atol( optarg ) * LAL_INT8_C(1000000000);
        break;
      case 'w': /* block-duration */
        localparams.duration = atof( optarg );
        break;
      case 'W': /* pad-data */
        localparams.padData = atof( optarg );
        break;
      case '!': /* h1-slide-segment */
        localparams.slideSegments[LAL_IFO_H1] = atoi( optarg );
        break;
      case '&': /* h2-slide-segments */
        localparams.slideSegments[LAL_IFO_H2] = atoi( optarg );
        break;
      case '(': /* l1-slide-segments */
        localparams.slideSegments[LAL_IFO_L1] = atoi( optarg );
        break;
      case ')': /* v1-slide-segments */
        localparams.slideSegments[LAL_IFO_V1] = atoi( optarg );
        break;
      case '@': /* Short slide offset time */
        localparams.shortSlideOffset = atoi( optarg );
        break;
      case 'V': /* version */
        XLALOutputVersionString(stderr, 0);
        exit( 0 );
     case '#': /* sky grid file */
        localparams.skyPositionsFile = optarg;
        break;
     case '|': /* FFT-level for plans */
        localparams.fftLevel = atoi( optarg );
        break;
      case '4': /* Cluster window */
        localparams.clusterWindow = atof(optarg);
        break;
      case '3': /* Injection search window */
        localparams.injSearchWindow = atof( optarg );
        break;
      case '5': /* Injection search window */
        localparams.injMchirpWindow = atof( optarg );
        break;
      case '7':
        if (!strcmp("real_4", optarg))
        {
          localparams.ligoDoubleData = 0;
        }
        else if (!strcmp("real_8", optarg))
        {
          localparams.ligoDoubleData = 1;
        }
        else
        {
          fprintf(stderr, "invalid argument to --%s:\n"
                  "unknown data type specified;\n"
                  "%s (must be one of: real_4, real_8)\n",
                  long_options[option_index].name, optarg);
        }
        break;
      case '8':
        if (!strcmp("real_4", optarg))
        {
          localparams.virgoDoubleData = 0;
        }
        else if (!strcmp("real_8", optarg))
        {
          localparams.virgoDoubleData = 1;
        }
        else
        {
          fprintf(stderr, "invalid argument to --%s:\n"
                  "unknown data type specified;\n"
                  "%s (must be one of: real_4, real_8)\n",
                  long_options[option_index].name, optarg);
        }
        break;
      case '?':
        error( "unknown error while parsing options\n" );
      default:
        error( "unknown error while parsing options\n" );
    }
  }

  if ( optind < argc )
  {
    fprintf( stderr, "extraneous command line arguments:\n" );
    while ( optind < argc )
      fprintf( stderr, "%s\n", argv[optind++] );
    exit( 1 );
  }

  /* set number of ifos */
  localparams.numIFO = 0;

  for ( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++ )
  {
    if ( localparams.haveTrig[ifoNumber] )
    {
      XLALReturnIFO(ifo,ifoNumber);
      snprintf( localparams.ifoName[localparams.numIFO], LIGOMETA_IFO_MAX,\
                "%s", ifo );
      localparams.numIFO++;
    }
  }
  
  /* check for H1H2 */
  if (localparams.numIFO == 2)
  {
    if (! strcmp(localparams.ifoName[0],"H1"))
    {
      if (! strcmp(localparams.ifoName[1],"H2"))
      {
        localparams.singlePolFlag = 1;
      }
    }
  }
  /* Set the faceOn-faceAway flag */
  /* Otherwise it takes default value of 0 */
  if (localparams.faceOnAnalysis)
  {
    localparams.faceOnStatistic = 1;
  }
  else if (localparams.faceAwayAnalysis)
  {
    localparams.faceOnStatistic = 2;
  }

  /* Set the number of points in the time arrays */
  localparams.numTimePoints = floor(\
          localparams.segmentDuration * localparams.sampleRate + 0.5);
  /* Set the number of points in the frequency arrays */
  localparams.numFreqPoints = localparams.numTimePoints / 2 + 1;

  /* For now we stick to only analysing half of each segment */
  localparams.strideDuration = 0.5 * localparams.segmentDuration;

  /* FIXME: Hardcoded to 1s */
  localparams.numBufferPoints = floor(localparams.sampleRate + 0.5);

  /* Choose the start and end point of each segment for analysis */
  if (localparams.analSegmentEnd)
  { /* Want to analyse from the end of the segment */
    /* Start from the end */
    localparams.analEndPoint = localparams.numTimePoints;
    /* Remove the spectrum truncation */
    localparams.analEndPoint -= floor(\
        0.5 * localparams.truncateDuration * localparams.sampleRate + 0.5); 
    /* Remove the buffer points */
    localparams.analEndPoint -= localparams.numBufferPoints;
    /* And set the start point */
    localparams.analStartPoint = localparams.analEndPoint - \
        0.5*localparams.numTimePoints;
  }
  else
  { /* DEFAULT: Analyse the middle of the segment */
    localparams.analStartPoint = 1*localparams.numTimePoints/4;
    localparams.analEndPoint = (3*localparams.numTimePoints)/4;
  }

  localparams.analStartTime = localparams.analStartPoint / \
                                localparams.sampleRate;
  localparams.analStartPointBuf = localparams.analStartPoint\
                                  - localparams.numBufferPoints;
  localparams.analEndTime = localparams.analEndPoint / \
                                localparams.sampleRate;
  localparams.analEndPointBuf = localparams.analEndPoint\
                               + localparams.numBufferPoints;
  localparams.numAnalPoints = localparams.analEndPoint\
                             - localparams.analStartPoint;
  localparams.numAnalPointsBuf = localparams.analEndPointBuf\
                                - localparams.analStartPointBuf;
  /* Max template length is start length minus PSD truncation */
  localparams.maxTempLength = localparams.analStartTime;
  localparams.maxTempLength -= localparams.truncateDuration/2.;

  /* Determine the number of short slides */
  if (localparams.doShortSlides)
  {
    localparams.numShortSlides = 1 + localparams.numIFO * (int) floor( \
        localparams.strideDuration / \
        (localparams.shortSlideOffset * (localparams.numIFO-1)) );
  }
  else
  {
    localparams.numShortSlides = 1;
  }

  *params = localparams;

  return 0;
}

/* sets default values for parameters */
int coh_PTF_default_params( struct coh_PTF_params *params )
{
  /* overall, default values are zero */
  memset( params, 0, sizeof( *params ) );

  /* set start time */
  XLALGPSTimeNow(&params->jobStartTime);

  /* FFT plan defaults to 1 */
  params->fftLevel = 1;

  /* No injections unless supplied */
  params->injectList = NULL;

  /* set default sky location params */
  params->rightAscension = -1000.;
  params->declination = -1000.;
  params->skyError = 0.;
  params->timingAccuracy = 0.0005;
  params->skyPositionsFile = NULL;
  params->skyLooping = ALL_SKY;

  /* Default injection search window to 1s */
  params->injSearchWindow = 1.;

  /* dynamic range factor must be greater than zero */
  params->dynRangeFac = 1.0;

  /* Various frequencies must be set */
  params->highpassFrequency     = -1.0; 
  params->lowTemplateFrequency = -1.0;
  params->lowFilterFrequency = -1.0;
  params->highFilterFrequency = -1.0;

  /* segments and templates to do: all of them */
  params->segmentsToDoList  = "^-$";
  params->templatesToDoList = "^-$";

  /* flags specifying what to do: default is to do everything */
  params->getBank     = 1;
  params->getData     = 1;
  params->getSpectrum = 1;
  params->doFilter    = 1;

  /* FIXME: Some thresholding stuff, this should be command line options */
  params->bankVeton = 3.;
  params->bankVetoq = 4.;
  params->autoVeton = 3.;
  params->autoVetoq = 4.;
  params->nullStatThreshold = 5.25;
  params->nullStatGradOn = 30.;
  params->nullStatGradient = 50./700.;

  params->approximant = NumApproximants;
  params->order = LAL_PNORDER_NUM_ORDER;

  /* numeric type for data, set to S6 data options */
  params->ligoDoubleData = 1;
  params->virgoDoubleData = 0;

  return 0;
}

/* check sanity of parameters and sets appropriate values of unset parameters */
int coh_PTF_params_sanity_check( struct coh_PTF_params *params )
{
  UINT4 recordLength   = 0;
  UINT4 segmentLength  = 0;
  UINT4 segmentStride  = 0;
  UINT4 truncateLength = 0;
  UINT4 ifoNumber;
  INT8  startTime;
  INT8  endTime;
//  UINT4 slideSegments; Currently unused FIXME

  if ( params->getSpectrum ) /* need data and response if not strain data */
    sanity_check( params->getData && (params->strainData) );

  /* parameters required to get data */
  if ( params->getData )
  {
    /* checks on data duration */
    startTime = epoch_to_ns( &params->startTime );
    endTime   = epoch_to_ns( &params->endTime );
    sanity_check( startTime > 0 );
    sanity_check( endTime > startTime );
    sanity_check( params->duration > 0 );
    sanity_check( 1e9*params->duration == ((endTime - startTime)) );

    /* checks on size of data record */
    sanity_check( params->sampleRate > 0 );
    recordLength = params->duration * params->sampleRate;
    sanity_check( recordLength > 0 );
    for( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
    {
      if ( params->haveTrig[ifoNumber] )
      {
        sanity_check(params->channel[ifoNumber]);
        sanity_check(params->dataCache[ifoNumber]);
      }
    }
  }

  /* parameters required to get spectrum */
  if ( params->getSpectrum )
  {
    /* checks on size of data segments and stride */
    sanity_check( params->psdSegmentDuration > 0 );
    segmentLength = floor(params->psdSegmentDuration*params->sampleRate + 0.5);
    sanity_check( recordLength / segmentLength > 0 );
    params->psdStrideDuration = 0.5 * params->psdSegmentDuration;
    segmentStride = floor(params->psdStrideDuration * params->sampleRate + 0.5);
    sanity_check( segmentStride > 0 );
    sanity_check( params->truncateDuration > 0);
    truncateLength = floor(params->truncateDuration * params->sampleRate + 0.5);
    sanity_check( truncateLength > 0 );

    /* checks on data input information */
    /*sanity_check( params->channel );*/
    sanity_check( params->dynRangeFac > 0.0 );
  }

  /* Sanity checks on data stuff */
  sanity_check( params->segmentDuration > 0 );
  segmentLength = floor(params->segmentDuration * params->sampleRate + 0.5);
  sanity_check( recordLength / segmentLength > 0 );
  segmentStride = floor(params->strideDuration * params->sampleRate + 0.5);
  sanity_check( segmentStride > 0 );
  sanity_check( !( (recordLength - segmentLength) % segmentStride ) );
  params->numOverlapSegments = 1 + (recordLength - segmentLength)/segmentStride;
  sanity_check( params->numTimePoints > 0);
  sanity_check( params->numFreqPoints > 0);
  sanity_check( (params->analStartPoint < segmentLength) );
  sanity_check( (params->analEndPoint < segmentLength) );
  sanity_check( (params->analEndPoint - params->analStartPoint) \
                 == segmentStride);

  /* sky localisation params */
  if ( params->rightAscension != -1000. && params->declination != -1000. )
  {

    sanity_check( params->rightAscension >= 0.\
                  && params->rightAscension <= 2.*LAL_PI);
    sanity_check( params->declination >= -LAL_PI/2.\
                  && params->declination <= LAL_PI/2.);

    if ( params->skyError>0. )
    {

      sanity_check( params->skyError >= -LAL_PI/2.\
                    && params->skyError <=LAL_PI/2. );
      sanity_check( params->timingAccuracy > 0. );

      if ( params->numIFO == 2 )
        params->skyLooping = TWO_DET_SKY_PATCH;
      else
        params->skyLooping = SKY_PATCH;
    }
    else
    {
      params->skyLooping = SINGLE_SKY_POINT;
    }
  }
  else if (params->skyPositionsFile)
  {
    params->skyLooping = SKY_PATCH;
  }

  else
  {
    if ( params->numIFO == 2 )
      params->skyLooping = TWO_DET_ALL_SKY;
    else
      params->skyLooping = ALL_SKY;
  }

  if ( params->skyLooping != SINGLE_SKY_POINT && params->doNullStream )
  {
    error( "--do-null-stream and --sky-error are incompatbile. Null stream on more than one sky point is a bad idea.\n" );
  }

  /* Check that filter frequencies have been given */
  sanity_check( params->highpassFrequency > 0);
  sanity_check( params->lowTemplateFrequency > 0 || params->dynTempLength != 0);
  sanity_check( params->lowFilterFrequency > 0 && params->lowFilterFrequency >= params->lowTemplateFrequency);
  sanity_check( params->highFilterFrequency > params->lowFilterFrequency);

  sanity_check( params->approximant != NumApproximants);
  sanity_check( params->order != LAL_PNORDER_NUM_ORDER);
  sanity_check( params->dynTempLength == 0 || params->approximant == FindChirpSP);

// This needs fixing. Need a check on whether segmentsToDoList and 
// analyzeInjSegsOnly have been given.
//  sanity_check( ! ((params->segmentsToDoList  != "^-$") && (params->analyzeInjSegsOnly)));

  return 0;
}

/* Sanity check for coh_PTF_inspiral specific */
int coh_PTF_params_inspiral_sanity_check( struct coh_PTF_params *params )
{
  INT8 trigTime;
  trigTime = epoch_to_ns( &params->trigTime );
  sanity_check( trigTime > 0 );
  sanity_check( params->threshold );
  if ( params->spinBank )
  {
    sanity_check( params->spinThreshold );
  }
  sanity_check( params->timeWindow );
// This sanity check needs fixing!
//  sanity_check( params->outputFile );
  if ( params->bankFile )
  {
    fprintf(stderr,"Please use --spin-bank and/or --non-spin-bank with this ");
    fprintf(stderr,"code and not --bank-file.\n");
    sanity_check(! params->bankFile );
  }
  if ( params->doBankVeto && (! params->bankVetoBankName) )
  {
    fprintf(stderr, "When using --do-bank-veto you must also supply ");
    fprintf(stderr, "--bank-veto-templates. \n" );
    sanity_check(!( params->doBankVeto && (! params->bankVetoBankName)));
  }
  if ( params->bankVetoBankName && (! params->doBankVeto) )
  {
    fprintf(stderr, "Supplying --bank-veto-templates will do nothing if ");
    fprintf(stderr, "--do-bank-veto is not given. \n" );
  }
  if ( params->doAutoVeto && (! (params->autoVetoTimeStep && params->numAutoPoints)))
  {
    fprintf(stderr, "When using --do-auto-veto you must also supply ");
    fprintf(stderr, "--num-auto-chisq-points and --auto-veto-time-step\n");
    sanity_check(params->doAutoVeto && params->autoVetoTimeStep && params->numAutoPoints);
  }
  if ( params->doChiSquare && (! params->numChiSquareBins))
  {
    fprintf(stderr, "When using --do-chi-square you must also supply ");
    fprintf(stderr, "--num-chi-square-bins \n");
    sanity_check(params->doChiSquare && params->numChiSquareBins);
  }
  sanity_check(params->spinBank || params->noSpinBank);
  if ( params->clusterFlag)
    sanity_check( params->clusterWindow);

  if (params->numIFO == 0)
  {
    fprintf(stderr, "You have not specified any detectors to analyse");
    return 1;
  }
  else if (params->numIFO == 1)
  {
    fprintf(stdout, "You have only specified one detector, "
                    "why are you using the coherent code? \n");
  }

  sanity_check( params->numShortSlides > 0);
  sanity_check( ! (params->doShortSlides && params->shortSlideOffset == 0));

  sanity_check( ! (params->writeSnglInspiralTable && (params->numIFO != 1)));

  return 0;
}

/* Sanity check for coh_PTF_spin_checker specific */
int coh_PTF_params_spin_checker_sanity_check( struct coh_PTF_params *params )
{
  sanity_check( params->spinSNR2threshold > 0 );
  sanity_check( params->nonspinSNR2threshold > 0 );
  sanity_check( params->spinBank );
  sanity_check( params->noSpinBank);

  return 0;
}


/* prints a help message */
int coh_PTF_usage( const char *program )
{
  fprintf( stderr, "usage: %s options\n", program );
  fprintf( stderr, "\ngeneral options:\n" );
  fprintf( stderr, "--help                     print this message\n" );
  fprintf( stderr, "--version                  print the version of the code\n" );
  fprintf( stderr, "--verbose                  print verbose messages while running\n" );

  fprintf( stderr, "\ndata reading options:\n" );
  fprintf( stderr, "--h1-data                  Analyze h1 data \n" );
  fprintf( stderr, "--h2-data                  Analyze h2 data \n" );
  fprintf( stderr, "--l1-data                  Analyze l1 data \n" );
  fprintf( stderr, "--v1-data                  Analyze v1 data \n" );
  fprintf( stderr, "--h1-frame-cache=cachefile    name of the frame cache file\n" );
  fprintf( stderr, "--h2-frame-cache=cachefile    name of the frame cache file\n" );
  fprintf( stderr, "--l1-frame-cache=cachefile    name of the frame cache file\n" );
  fprintf( stderr, "--v1-frame-cache=cachefile    name of the frame cache file\n" );
  fprintf( stderr, "--h1-channel-name             data channel to analyze\n" );
  fprintf( stderr, "--h2-channel-name             data channel to analyze\n" );
  fprintf( stderr, "--l1-channel-name             data channel to analyze\n" );
  fprintf( stderr, "--v1-channel-name             data channel to analyze\n" );
  fprintf( stderr, "--gps-start-time=tstart    GPS start time of data to analyze (sec)\n" );
  fprintf( stderr, "--gps-start-time-ns=tstartns  nanosecond residual of start time\n" );
  fprintf( stderr, "--gps-end-time=tstop       GPS stop time of data to analyze (sec)\n" );
  fprintf( stderr, "--gps-end-time-ns=tstopns  nanosecond residual of stop time\n" );
  fprintf( stderr, "\nsimulated data options:\n" );
  fprintf( stderr, "--simulated-data=dataType  create simulated Gaussian noise. Can be WhiteNoise,ILIGONoise or ALIGONoise. \n" );
  fprintf( stderr, "--random-seed=seed         random number seed for simulated data\n" );
  fprintf( stderr, "--sample-rate=srate        sampling rate of simulated data (Hz)\n" );
  fprintf( stderr, "--zero-data                create a time series of zeros\n" );

  fprintf( stderr, "\ndata conditioning options:\n" );
  fprintf( stderr, "--highpass-frequency=fhi   high-pass filter data at frequency fhi (Hz)\n" );
  fprintf( stderr, "--sample-rate=srate        decimate data to be at sample rate srate (Hz)\n" );

  fprintf( stderr, "\ncalibration options:\n" );
  fprintf( stderr, "--strain-data              data is strain (already calibrated)\n" );
  fprintf( stderr, "--ligo-calibrated-data=TYPE   LIGO calibrated data of TYPE real_4 or real_8\n");
  fprintf( stderr, "--virgo-calibrated-data=TYPE   Virgo calibrated data of TYPE real_4 or real_8\n");
  fprintf( stderr, "--strain-data              data is strain (already calibrated)\n" );
  fprintf( stderr, "--dynamic-range-factor=dynfac  scale calibration by factor dynfac\n" );
  fprintf( stderr, "--fft-level=PLAN Set the fft plan to use level=PLAN\n" );

  fprintf( stderr, "\ndata segmentation options:\n" );
  fprintf( stderr, "--segment-duration=duration  duration of a data segment for filtering (sec)\n" );
  fprintf( stderr, "--psd-segment-duration=duration  duration of a data segment for PSD generation (sec)\n" );
  fprintf( stderr, "--block-duration=duration    duration of an analysis block (sec)\n" );
  fprintf( stderr, "--pad-data=duration          input data padding (sec)\n" );
  fprintf( stderr, "--h1-slide-segment=amount    amount to be slid H1\n" );
  fprintf( stderr, "--h2-slide-segment=amount    amount to be slid H2\n" );
  fprintf( stderr, "--l1-slide-segment=amount    amount to be slid L1\n" );
  fprintf( stderr, "--v1-slide-segment=amount    amount to be slid V1\n" );
  fprintf( stderr, "--do-short-slides  Enabling sliding within the analysis segments. \n");
  fprintf( stderr, "--short-slide-offset Sets the slide amount between ifos when doing short slides.\n");

  fprintf( stderr, "\npower spectrum options:\n" );
  fprintf( stderr, "--theoretical-spectrum      take the PSD as the PSD used to generate the simulated data\n" );
  fprintf( stderr, "--low-template-freq=fmin    low frequency cutoff for generation of templates (Hz)\n" );
  fprintf( stderr, "--low-filter-freq=f_low    low frequency cutoff for matched filtering (Hz)\n" );
  fprintf( stderr, "--high-filter-freq=f_max    high frequency cutoff for matched filtering (Hz)\n" );
  fprintf( stderr, "--inverse-spec-length=t    set length of inverse spectrum to t seconds\n" );
  fprintf( stderr, "\nbank generation options:\n" );
  fprintf( stderr, "--bank-file=name           Location of tmpltbank xml file\n" );
  fprintf( stderr, "--spin-bank=name   Location of output spin bank for spin checker or input spin bank for cohPTF_inspiral \n");
  fprintf( stderr, "--non-spin-bank=name   Location of output non spin bank for spin checker or input non spin bank for cohPTF_inspiral \n");
  fprintf( stderr, "\nfiltering options:\n" );
  fprintf( stderr, "--only-segment-numbers=seglist  list of segment numbers to compute\n" );
  fprintf( stderr, "--analyze-inj-segs-only  Only analyze times when injections have been made\n" );
  fprintf( stderr, "--only-template-numbers=tmpltlist  list of filter templates to use\n" );
  fprintf( stderr, "--face-on-analysis  Run with templates demanding inclination=0\n" );
  fprintf( stderr, "--face-away-analysis  Run with templates demanding inclination=pi/2\n" );
  fprintf( stderr, "--dynamic-template-length Run with templates whose length is dynamically set to be close to the maximum possible.\n");
  fprintf( stderr, "--analyse-segment-end Rather than analyse the middle half of analysis segments (1/4 to 3/4 of length), which is the default behaviour. Analyse from (1/2 - truncDuration to 1 - truncDuration of length) where truncDuration is the spectrum truncation which cannot be analysed. This is the closest to the end of the segment it is possible to analyse \n");

  fprintf( stderr, "\nsky location options:\n" );
  fprintf( stderr, "--right-ascension=ra       right ascension of external trigger in degrees\n" );
  fprintf( stderr, "--declination=dec          declination of external trigger in degrees\n" );
  fprintf( stderr, "--sky-error=err            1-sigma error radius in sky location of external trigger in degrees\n" );
  fprintf( stderr, "--timing-accuracy=t_acc    Accuracy ( in seconds ) of timing information\n" );
  fprintf( stderr, "--sky-positions-file=name  Location of sky locations file for IPN\n" );

  fprintf( stderr, "\ninjection options:\n" );
  fprintf( stderr, "--injection-file=file list of software injections to make into the data. If this option is not given injections are not made\n");
  fprintf( stderr, "--inj-search-window=arg    output injection triggers only within arg of the injections\n");
  fprintf( stderr, "--inj-mchirp-window=arg    search for injections only with templates with fractional mchirp distance within arg\n");

  fprintf( stderr, "\nTrigger extraction options:\n" );
  fprintf( stderr, "--snr-threshold=threshold Only keep triggers with a snr above threshold\n" );
  fprintf( stderr, "--spin-snr-threshold=threshold Only keep spinning triggers with a snr above threshold\n" );
  fprintf( stderr, "--sngl-snr-threshold Only keep triggers if at least one ifo has snr above value\n" );
  fprintf( stderr, "--non-spin-snr2-threshold=value SNR squared value over which a non spin trigger is considered found for spin checker program\n" );
  fprintf( stderr, "--spin-snr2-threshold=value SNR squared value over which a spin trigger is considered found for spin checker program\n" );
  fprintf( stderr, "--trig-time-window=window Keep loudest trigger within window seconds\n" );
  fprintf( stderr, "--do-null-stream Calculate Null SNR for potential triggers\n");
  fprintf( stderr, "--do-trace-snr Calculate Trace SNR for potential triggers \n");
  fprintf( stderr, "--do-bank-veto Calculate Bank Veto for potential triggers \n");
  fprintf( stderr, "--bank-veto-templates File containing templates to use for bank veto \n");
  fprintf( stderr, "--do-auto-veto Calculate Auto Veto for potential triggers \n");
  fprintf( stderr, "--do-chi-square Calculate the chi squared value for potential triggers \n");
  fprintf( stderr, "--num-auto-chisq-points Number of points to use in calculating auto veto \n");
  fprintf( stderr, "--auto-veto-time-step Seperation between points for auto veto \n");
  fprintf( stderr, "--num-chi-square-bins Number of bins to use to calculate chi square \n");
  fprintf (stderr, "--chi-square-threshold Only calculate chi square if detection statistic is above this threshold \n");
  fprintf (stderr, "--store-amplitude-params Calculate and store the amplitude params in the multi_inspiral table \n");

  fprintf( stderr, "\ntrigger output options:\n" );
  fprintf( stderr, "--output-file=outfile      output triggers to file outfile\n" );
  fprintf( stderr, "--trig-start-time=sec      output only triggers after GPS time sec. CURRENTLY NONFUNCTIONAL\n" );
  fprintf( stderr, "--trig-end-time=sec        output only triggers before GPS time sec. CURRENTLY NONFUNCTIONAL\n" );
  fprintf( stderr, "--ifo-tag=string           set ifotag to string for file naming\n" );
  fprintf( stderr, "--user-tag=string          set the process_params usertag to string\n" );
  fprintf( stderr, "--do-clustering            turn on clustering\n");
  fprintf( stderr, "--cluster-window=arg       cluster window length\n");

  fprintf( stderr, "\nintermediate data output options:\n" );
  fprintf( stderr, "--write-raw-data           write raw data before injection or conditioning\n" );
  fprintf( stderr, "--write-data               write data after injection and conditioning\n" );
  fprintf( stderr, "--write-inv-spectrum       write inverse power spectrum\n" );
  fprintf( stderr, "--write-segment            write overwhitened data segments\n" );
  fprintf( stderr, "--write-filter-output      write filtered data segments\n" );

  return 0;
}
