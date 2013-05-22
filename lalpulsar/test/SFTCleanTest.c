/*
*  Copyright (C) 2007 Badri Krishnan, Gregory Mendell
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
#include <lal/SFTClean.h>
#include <glob.h>

/* REVISIONS: */
/* 09/09/05 gam; make RandomParams *randPar a parameter for LALCleanCOMPLEX8SFT. Thus only need to */
/*               initialze RandomParams *randPar once and avoid repeatly opening /dev/urandom.  */

/**
\author Krishnan, B.
\file
\ingroup SFTClean_h
*/

/**\name Error Codes */
/*@{*/
#define SFTCLEANTESTC_ENORM     0
#define SFTCLEANTESTC_ESUB      1
#define SFTCLEANTESTC_EARG      2
#define SFTCLEANTESTC_EBAD      3
#define SFTCLEANTESTC_EFILE     4
#define SFTCLEANTESTC_ERANDFILE 5
#define SFTCLEANTESTC_ERANDSEED 6

#define SFTCLEANTESTC_MSGENORM     "Normal exit"
#define SFTCLEANTESTC_MSGESUB      "Subroutine failed"
#define SFTCLEANTESTC_MSGEARG      "Error parsing arguments"
#define SFTCLEANTESTC_MSGEBAD      "Bad argument values"
#define SFTCLEANTESTC_MSGEFILE     "Could not create output file"
#define SFTCLEANTESTC_MSGERANDFILE "Could not open /dev/urandom"
#define SFTCLEANTESTC_MSGERANDSEED "Error initializing random seed"
/*@}*/

/** \cond DONT_DOXYGEN */

/* Default parameters. */



#define MAXFILENAMELENGTH 256
/* defaults chosen for L1 */
#define HARMONICFILE "./harmonics.txt"
/*#define INPUTSFTDIR "/nfs/morbo/geo600/hannover/sft/S2-LIGO/S2_L1_Funky-v3Calv5DQ30MinSFTs"*/
#define INPUTSFTDIR "./S3H1data/"
/*#define OUTPUTSFTDIR "/nfs/morbo/geo600/hannover/sft/S2-LIGO-clean/S2_L1_Funky-v3Calv5DQ30MinSFTs-clean"*/
#define OUTPUTSFTDIR "./test/"
#define STARTFREQ 950.0
#define BANDFREQ 20.0
#define MAXFILES 3000 /* maximum number of files to read in a directory */
#define WINDOWSIZE 100
#define MAXBINS 20



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
    ERROR( SFTCLEANTESTC_ESUB, SFTCLEANTESTC_MSGESUB,      \
           "Function call \"" #func "\" failed:" );                  \
    return SFTCLEANTESTC_ESUB;                                  \
  }                                                                  \
} while (0)
/******************************************************************/

/* A global pointer for debugging. */
#ifndef NDEBUG
char *lalWatch;
#endif

#define TRUE (1==1)
#define FALSE (1==0)

/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvv------------------------------------ */
int main(int argc, char *argv[]){
  static LALStatus       status;  /* LALStatus pointer */
  static SFTtype         *sft;
  static LineNoiseInfo   lines, lines2;
  static LineHarmonicsInfo harmonics;

  INT4 j;
  INT4 nLines=0, count1, nHarmonicSets;
  INT4 mObsCoh;
  CHAR filelist[MAXFILES][MAXFILENAMELENGTH];
  CHAR tempstr1[256], tempstr2[256];

  /* log file and strings */
  FILE   *fpLog=NULL;
  CHAR   *fnameLog=NULL;
  CHAR   *logstr=NULL;

  /* user input variables */
  BOOLEAN uvar_help;
  CHAR *uvar_harmonicfname;        /* file with harmonics info */
  CHAR *uvar_inputSFTDir;    /* directory for unclean sfts */
  CHAR *uvar_outputSFTDir;   /* directory for cleaned sfts */
  REAL8 uvar_fStart, uvar_fBand;
  INT4  uvar_window, uvar_maxBins;

  /* 09/09/05 gam; randPar now a parameter for LALCleanCOMPLEX8SFT */
  FILE *fp=NULL;
  INT4 seed, ranCount;
  RandomParams *randPar=NULL;

  /* set defaults */

  /* LALDebugLevel must be called before anything else */
  SUB( LALGetDebugLevel( &status, argc, argv, 'd'), &status);

  uvar_help = FALSE;

  uvar_harmonicfname = (CHAR *)LALMalloc(256 * sizeof(CHAR));
  strcpy(uvar_harmonicfname,HARMONICFILE);

  uvar_inputSFTDir = (CHAR *)LALMalloc(256 * sizeof(CHAR));
  strcpy(uvar_inputSFTDir, INPUTSFTDIR);

  uvar_outputSFTDir = (CHAR *)LALMalloc(256 * sizeof(CHAR));
  strcpy(uvar_outputSFTDir, OUTPUTSFTDIR);

  uvar_fStart = STARTFREQ;
  uvar_fBand = BANDFREQ;
  uvar_window = WINDOWSIZE;
  uvar_maxBins = MAXBINS;

  /* register user input variables */
  SUB( LALRegisterBOOLUserVar(   &status, "help",            'h', UVAR_HELP,     "Print this message",                          &uvar_help),            &status);
  SUB( LALRegisterSTRINGUserVar( &status, "inputSFTDir",     'i', UVAR_OPTIONAL, "Input SFT Directory",                         &uvar_inputSFTDir),     &status);
  SUB( LALRegisterSTRINGUserVar( &status, "outputSFTDir",    'o', UVAR_OPTIONAL, "Output SFT Directory",                        &uvar_outputSFTDir),    &status);
  SUB( LALRegisterSTRINGUserVar( &status, "harmonicfname",   'H', UVAR_OPTIONAL, "File with list of lines",                     &uvar_harmonicfname),   &status);
  SUB( LALRegisterREALUserVar(   &status, "fStart",          'f', UVAR_OPTIONAL, "Frequency to start cleaning",                 &uvar_fStart),          &status);
  SUB( LALRegisterREALUserVar(   &status, "fBand",           'b', UVAR_OPTIONAL, "Frequency Band",                              &uvar_fBand),           &status);
  SUB( LALRegisterINTUserVar(    &status, "window",          'w', UVAR_OPTIONAL, "No. of bins for generating random numbers",   &uvar_window),          &status);
  SUB( LALRegisterINTUserVar(    &status, "maxBins",         'm', UVAR_OPTIONAL, "Max. No. of bins to clean",                   &uvar_maxBins),         &status);

  /* read all command line variables */
  SUB( LALUserVarReadAllInput(&status, argc, argv), &status);

  /* exit if help was required */
  if (uvar_help)
    exit(0);

  /* 09/09/05 gam; randPar now a parameter for LALCleanCOMPLEX8SFT */
  fp=fopen("/dev/urandom", "r");
  if (!fp) {
     fprintf(stderr,"Error in SFTCleanTest. Error code %i; desc: %s \n",SFTCLEANTESTC_ERANDFILE,SFTCLEANTESTC_MSGERANDFILE);
     exit(1);
  }
  ranCount = fread(&seed, sizeof(seed), 1, fp);
  if (!(ranCount==1)) {
     fprintf(stderr,"Error in SFTCleanTest. Error code %i; desc: %s \n",SFTCLEANTESTC_ERANDSEED,SFTCLEANTESTC_MSGERANDSEED);
     exit(1);
  }
  fclose(fp);
  SUB ( LALCreateRandomParams (&status, &randPar, seed), &status );

  /********logging the user input variables*********************/
  /* open log file for writing */
  fnameLog = (CHAR *)LALMalloc( 512*sizeof(CHAR));
  strcpy(fnameLog,uvar_outputSFTDir);
  strcat(fnameLog, "/CLEAN_log");
  if ((fpLog = fopen(fnameLog, "w")) == NULL) {
    fprintf(stderr, "Unable to open file %s for writing\n", fnameLog);
    LALFree(fnameLog);
    exit(1);
  }

  /* get the log string */
  SUB( LALUserVarGetLog(&status, &logstr, UVAR_LOGFMT_CFGFILE), &status);

  fprintf( fpLog, "## LOG FILE FOR SFT Cleaning\n\n");
  fprintf( fpLog, "# User Input:\n");
  fprintf( fpLog, "#-------------------------------------------\n");
  fprintf( fpLog, logstr);
  LALFree(logstr);

  /* copy contents of harmonics file into logfile */
  fprintf(fpLog, "\n\n# Contents of harmonics file:\n");
  fclose(fpLog);
  {
    CHAR command[1024] = "";
    sprintf(command, "cat %s >> %s", uvar_harmonicfname, fnameLog);
    system(command);
  }

  /* append an ident-string defining the exact CVS-version of the code used */
  fpLog = fopen(fnameLog, "a");
  {
    CHAR command[1024] = "";
    fprintf (fpLog, "\n\n# CVS-versions of executable:\n");
    fprintf (fpLog, "# -----------------------------------------\n");
    fclose (fpLog);

    sprintf (command, "ident %s | sort -u >> %s", argv[0], fnameLog);
    system (command);	/* we don't check this. If it fails, we assume that */
    			/* one of the system-commands was not available, and */
    			/* therefore the CVS-versions will not be logged */

    LALFree(fnameLog);
  }
  /*end of logging*********************************************************/

  SUB( LALFindNumberHarmonics (&status, &harmonics, uvar_harmonicfname), &status);
  nHarmonicSets = harmonics.nHarmonicSets;

  if (nHarmonicSets > 0)
    {
      harmonics.startFreq = (REAL8 *)LALMalloc(harmonics.nHarmonicSets * sizeof(REAL8));
      harmonics.gapFreq = (REAL8 *)LALMalloc(harmonics.nHarmonicSets * sizeof(REAL8));
      harmonics.numHarmonics = (INT4 *)LALMalloc(harmonics.nHarmonicSets * sizeof(INT4));
      harmonics.leftWing = (REAL8 *)LALMalloc(harmonics.nHarmonicSets * sizeof(REAL8));
      harmonics.rightWing = (REAL8 *)LALMalloc(harmonics.nHarmonicSets * sizeof(REAL8));


      SUB( LALReadHarmonicsInfo( &status, &harmonics, uvar_harmonicfname ), &status);

      nLines = 0;
      for (count1=0; count1 < nHarmonicSets; count1++)
	{
	  nLines += *(harmonics.numHarmonics + count1);
	}

      lines.nLines = nLines;
      lines.lineFreq = (REAL8 *)LALMalloc(nLines * sizeof(REAL8));
      lines.leftWing = (REAL8 *)LALMalloc(nLines * sizeof(REAL8));
      lines.rightWing = (REAL8 *)LALMalloc(nLines * sizeof(REAL8));

      SUB( LALHarmonics2Lines( &status, &lines, &harmonics), &status);


      lines2.nLines = nLines;
      lines2.lineFreq = (REAL8 *)LALMalloc(nLines * sizeof(REAL8));
      lines2.leftWing = (REAL8 *)LALMalloc(nLines * sizeof(REAL8));
      lines2.rightWing = (REAL8 *)LALMalloc(nLines * sizeof(REAL8));

      SUB( LALChooseLines( &status, &lines2, &lines, uvar_fStart, uvar_fStart + uvar_fBand), &status);
      nLines = lines2.nLines;
    }


{
    CHAR     command[256];
    glob_t   globbuf;
    INT4    jj;

    strcpy(command, uvar_inputSFTDir);
    strcat(command, "/*SFT*.*");

    globbuf.gl_offs = 1;
    glob(command, GLOB_ERR, NULL, &globbuf);

    if(globbuf.gl_pathc==0)
      {
	fprintf(stderr,"No SFTs in directory %s ... Exiting.\n", uvar_inputSFTDir);
	return 1;  /* stop the program */
      }

    /* we will read up to a certain number of SFT files, but not all
       if there are too many ! */
    mObsCoh = (MAXFILES < globbuf.gl_pathc) ? MAXFILES : globbuf.gl_pathc ;

    /* Remember to do the following:
       globfree(&globbuf); after reading the file names. The file names are
       globbuf.gl_pathv[fileno]   that one can copy into whatever as:
       strcpy(filelist[fileno],globbuf.gl_pathv[fileno]);  */

    for (jj=0; jj < mObsCoh; jj++){
      strcpy(filelist[jj],globbuf.gl_pathv[jj]);
    }
    globfree(&globbuf);
  }


  for (j=0; j < mObsCoh; j++)
    {
      sft=NULL;
      SUB (LALReadSFTfile (&status, &sft, uvar_fStart, uvar_fStart + uvar_fBand, filelist[j]), &status);

      /* clean the sft */
      if (nLines > 0)
/*      SUB( LALCleanCOMPLEX8SFT( &status, sft, uvar_maxBins, uvar_window, &lines2), &status); */ /* 09/09/05 gam */
	SUB( LALCleanCOMPLEX8SFT( &status, sft, uvar_maxBins, uvar_window, &lines2, randPar), &status);

      /* make the output sft filename */
      sprintf(tempstr1, "%d", sft->epoch.gpsSeconds);
      strcpy(tempstr2,uvar_outputSFTDir);
      strcat(tempstr2, "/CLEAN_SFT.");
      strcat(tempstr2, tempstr1);

      /* write the sft */
      SUB( LALWriteSFTfile( &status, sft, tempstr2),  &status );

      LALDestroySFTtype (&status, &sft);
    }


  /* Free memory */
  if (nLines > 0)
    {
      LALFree(lines.lineFreq);
      LALFree(lines.leftWing);
      LALFree(lines.rightWing);

      LALFree(lines2.lineFreq);
      LALFree(lines2.leftWing);
      LALFree(lines2.rightWing);

    }

  if (nHarmonicSets > 0)
    {
      LALFree(harmonics.startFreq);
      LALFree(harmonics.gapFreq);
      LALFree(harmonics.numHarmonics);
      LALFree(harmonics.leftWing);
      LALFree(harmonics.rightWing);
    }

  /* 09/09/05 gam; randPar now a parameter for LALCleanCOMPLEX8SFT */
  SUB ( LALDestroyRandomParams (&status, &randPar), &status);

  SUB (LALDestroyUserVars(&status), &status);

  LALCheckMemoryLeaks();

  INFO( SFTCLEANTESTC_MSGENORM );
  return SFTCLEANTESTC_ENORM;
}

/** \endcond */
