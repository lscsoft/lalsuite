/*
*  Copyright (C) 2007 Badri Krishnan
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
 * \author Badri Krishnan
 */

#include <lal/SFTClean.h>
#include <glob.h> 

/* Error codes and messages */

/**\name Error Codes */ /*@{*/
#define SFTCLEANC_ENORM 0
#define SFTCLEANC_ESUB  1
#define SFTCLEANC_EARG  2
#define SFTCLEANC_EBAD  3
#define SFTCLEANC_EFILE 4

#define SFTCLEANC_MSGENORM "Normal exit"
#define SFTCLEANC_MSGESUB  "Subroutine failed"
#define SFTCLEANC_MSGEARG  "Error parsing arguments"
#define SFTCLEANC_MSGEBAD  "Bad argument values"
#define SFTCLEANC_MSGEFILE "Could not create output file"
/*@}*/


/* Default parameters. */



#define MAXFILENAMELENGTH 256
/* defaults chosen for L1 */

#define HARMONICSFILE "./harmonicsS2LHO4K_200_400.txt" 
#define NSTARFILE "/home/badkri/S2results/ALLSKYMAX_H1"
#define CLEANNSTARFILE "./ALLSKYMAXCLEAN"
#define WINDOWSIZE 100
#define TRUE (1==1)
#define FALSE (1==0)


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
    ERROR( SFTCLEANC_ESUB, SFTCLEANC_MSGESUB,      \
           "Function call \"" #func "\" failed:" );                  \
    return SFTCLEANC_ESUB;                                  \
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
  static LineNoiseInfo   lines, lines2;
  static LineHarmonicsInfo harmonics; 
  
  INT4 nLines=0, count1, nHarmonicSets;
  INT4 j, r;
  FILE *fpNstar=NULL;
  FILE *fpOut=NULL;
  INT4 *nstarVec=NULL;
  REAL8 *freqVec=NULL;
  INT4 nstarNum=0, minBin, maxBin;
  REAL8 timeBase=1800.0, minFreq, maxFreq;
  INT4 flag= 0;
 
  /* user input variables */
  BOOLEAN uvar_help;
  CHAR *uvar_nstarfile=NULL;
  CHAR *uvar_harmonicsfile=NULL;
  CHAR *uvar_cleannstar=NULL;

  /* set defaults */

  uvar_help = FALSE;  

  uvar_nstarfile = (CHAR *)LALMalloc(512*sizeof(CHAR));
  strcpy(uvar_nstarfile, NSTARFILE);

  uvar_harmonicsfile = (CHAR *)LALMalloc(512*sizeof(CHAR));
  strcpy(uvar_harmonicsfile, HARMONICSFILE);

  uvar_cleannstar = (CHAR *)LALMalloc(512*sizeof(CHAR));
  strcpy(uvar_cleannstar, CLEANNSTARFILE);

  /* register user input variables */
  SUB( LALRegisterBOOLUserVar(   &status, "help",            'h', UVAR_HELP,     "Print this message",            &uvar_help),            &status);   
  SUB( LALRegisterSTRINGUserVar( &status, "nstarfile",       'n', UVAR_OPTIONAL, "File with max number counts",   &uvar_nstarfile),       &status);
  SUB( LALRegisterSTRINGUserVar( &status, "harmonicsfile",   'l', UVAR_OPTIONAL, "File with line information",    &uvar_harmonicsfile),   &status);
  SUB( LALRegisterSTRINGUserVar( &status, "cleannstarfile",  'o', UVAR_OPTIONAL, "Output file",                   &uvar_cleannstar),      &status);


  /* read all command line variables */
  SUB( LALUserVarReadAllInput(&status, argc, argv), &status);

  /* exit if help was required */
  if (uvar_help)
    exit(0); 
  

  /* find number of harmonics */
  SUB( LALFindNumberHarmonics (&status, &harmonics, uvar_harmonicsfile), &status); 
  nHarmonicSets = harmonics.nHarmonicSets; 

  /* convert harmonics to explicit lines */
  if (nHarmonicSets > 0)
    {
      harmonics.startFreq = (REAL8 *)LALMalloc(harmonics.nHarmonicSets * sizeof(REAL8));
      harmonics.gapFreq = (REAL8 *)LALMalloc(harmonics.nHarmonicSets * sizeof(REAL8));
      harmonics.numHarmonics = (INT4 *)LALMalloc(harmonics.nHarmonicSets * sizeof(INT4));
      harmonics.leftWing = (REAL8 *)LALMalloc(harmonics.nHarmonicSets * sizeof(REAL8));
      harmonics.rightWing = (REAL8 *)LALMalloc(harmonics.nHarmonicSets * sizeof(REAL8));
    

      SUB( LALReadHarmonicsInfo( &status, &harmonics, uvar_harmonicsfile ), &status);
      
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
    }


  /* read file with nstar information */
  fpNstar = fopen( uvar_nstarfile, "r");
  if (fpNstar == NULL)
    {
      fprintf(stderr, "nstar file %s not found\n", uvar_nstarfile);
      exit(0);
    }
  /* first find number of nstar values */ 
  nstarNum = 0;
  do {
    REAL8 temp1, temp2, temp3, temp4, temp5;
    r=fscanf(fpNstar,"%lf%lf%lf%lf%lf\n", &temp1, &temp2, &temp3, &temp4, &temp5); 
    /* make sure the line has the right number of entries or is EOF */
    if (r==5) nstarNum++;
  } while ( r != EOF);
  rewind(fpNstar);

  /* allocate memory for vector of nstar and freqstar */
  nstarVec = (INT4 *)LALMalloc(nstarNum*sizeof(INT4));
  freqVec  = (REAL8 *)LALMalloc(nstarNum*sizeof(REAL8));

  /* read nstarfile again */
  for (j=0; j<nstarNum; j++)
    {
      REAL8 temp1, temp2, temp3, temp4, temp5;
      r=fscanf(fpNstar,"%lf%lf%lf%lf%lf\n", &temp1, &temp2, &temp3, &temp4, &temp5); 
      if (r==5)
	{
	  nstarVec[j] = temp1;
	  freqVec[j] = temp2;
	}
    }

  /* now we are done with nstarfile */
  fclose(fpNstar);

  /* calculate min and max frequency bins in nstarfile */
  /* this assumes frequencies are ordered */
  minBin = floor(freqVec[0]*timeBase + 0.5);
  maxBin = floor(freqVec[nstarNum-1]*timeBase + 0.5);
  minFreq = freqVec[0];
  maxFreq = freqVec[nstarNum-1];

  /* choose lines affecting frequency range */
  lines2.nLines = nLines;
  lines2.lineFreq = (REAL8 *)LALMalloc(nLines * sizeof(REAL8));
  lines2.leftWing = (REAL8 *)LALMalloc(nLines * sizeof(REAL8));
  lines2.rightWing = (REAL8 *)LALMalloc(nLines * sizeof(REAL8));

  SUB( LALChooseLines (&status, &lines2, &lines, minFreq, maxFreq), &status); 
  SUB( LALCheckLines (&status, &flag, &lines2, 202.0), &status); 


  /* clean the nstarVec */
  if ( nstarNum > 0 )
    {
      /* loop over lines and clean nstarVec*/
      for (j=0; j<nLines; j++)
	{
	  REAL8 lineFreq, leftWing, rightWing;
	  INT4 lineBin, leftWingBins, rightWingBins;
	  INT4 leftCount, rightCount;

          lineFreq = lines.lineFreq[j];
	  leftWing = lines.leftWing[j] + lineFreq*VTOT;
	  rightWing = lines.rightWing[j] + lineFreq*VTOT;

	  lineBin = floor( lineFreq*timeBase + 0.5);
	  leftWingBins = floor( leftWing*timeBase + 0.5);
	  rightWingBins = floor( rightWing*timeBase + 0.5);

	  if ( (lineBin >= minBin) && (lineBin <= maxBin) )
	    {
	      /* clean lineBin */
	      nstarVec[lineBin-minBin] = 0;

	      /* now go left */	      
	      for (leftCount=0; leftCount < leftWingBins; leftCount++)
		{
		  if (lineBin-minBin-leftCount > 0)
		    {
		      nstarVec[lineBin-minBin-leftCount-1] = 0;
		    } 
		} /* end of loop over leftwing */

	      /* now go right */
	      for (rightCount = 0; rightCount < rightWingBins; rightCount++)
		{
		  if (maxBin - lineBin - rightCount > 0)
		    {
		      nstarVec[lineBin-minBin+rightCount+1]=0;
		    }
		} /* end of loop over rightwing */

	    } /* end of if statement */
	} /* end of loop over lines */
    } /* end of cleaning */
  

  /* now write the result */
  fpOut = fopen(uvar_cleannstar, "w");
  if (fpOut == NULL)
    {
      fprintf(stderr, "unable to open file %s for writing \n", uvar_cleannstar);
      exit(0);
    }
  for (j=0; j<nstarNum; j++)
    {
      fprintf(fpOut, "%d  %f \n", nstarVec[j], freqVec[j]);
    }
  fclose(fpOut);

  /* Free memory */
  LALFree(nstarVec);
  LALFree(freqVec);


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

  SUB (LALDestroyUserVars(&status), &status);
  LALCheckMemoryLeaks(); 

  INFO( SFTCLEANC_MSGENORM );
  return SFTCLEANC_ENORM;
}

/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */













