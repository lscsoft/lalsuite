/*-----------------------------------------------------------------------
 *
 * File Name: SFTClean.c
 *
 * Authors: Krishnan, B.
 *
 * Revision: $Id$
 *
 * History:   Created by B. Krishnan on Feb 22, 2004
 *            Moved from lalapps to lal on July 31, 2005
 *
 *-----------------------------------------------------------------------
 */

/************************************ <lalVerbatim file="SFTCleanCV">
Author: Sintes, A.M., Krishnan, B.
$Id$
************************************* </lalVerbatim> */

/* <lalLaTeX>  *******************************************************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\section{Module \texttt{SFTClean.c}}
\label{ss:SFTClean.c}
Routines for reading SFT binary files

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Prototypes}
\vspace{0.1in}
\input{SFTCleanD}
\index{\verb&ReadSFTCleanHeader1()&}
\index{\verb&ReadCOMPLEX8SFTCleanData1()&}
\index{\verb&ReadCOMPLEX16SFTCleanData1()&}
\index{\verb&COMPLEX8SFT2Periodogram1()&}
\index{\verb&COMPLEX16SFT2Periodogram1()&}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Description}

SFT cleaning routines 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Uses}
\begin{verbatim}
LALHO()
\end{verbatim}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Notes}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\vfill{\footnotesize\input{SFTCleanCV}}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

*********************************************** </lalLaTeX> */


#include <lal/SFTClean.h>



NRCSID (SFTCLEANC, "$Id$");

/*
 * The functions that make up the guts of this module
 */


/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */


/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
/********************************* <lalVerbatim file="SFTCleanD"> */
void LALFindNumberHarmonics (LALStatus    *status,
			  LineHarmonicsInfo   *harmonicInfo,
			  CHAR         *fname)
{/*   *********************************************  </lalVerbatim> */
 /* this function finds the number of harmonic sets in file "fname" and
    checks the file format */

  FILE *fp = NULL;
  CHAR  dump[128];
  INT4   harmonicCount, r, tempint; 
  REAL8 temp1, temp2, temp3, temp4;   


  INITSTATUS (status, "FindNumberHarmonics", SFTCLEANC);
  ATTATCHSTATUSPTR (status);

  /* make sure arguments are not null */
  ASSERT (harmonicInfo, status, SFTCLEANH_ENULL, SFTCLEANH_MSGENULL ); 
  ASSERT (harmonicInfo->nHarmonicSets > 0, status, SFTCLEANH_EVAL, SFTCLEANH_MSGEVAL);
  ASSERT (harmonicInfo->startFreq, status, SFTCLEANH_ENULL, SFTCLEANH_MSGENULL);
  ASSERT (harmonicInfo->gapFreq, status, SFTCLEANH_ENULL, SFTCLEANH_MSGENULL);
  ASSERT (harmonicInfo->numHarmonics, status, SFTCLEANH_ENULL, SFTCLEANH_MSGENULL);
  ASSERT (harmonicInfo->leftWing, status, SFTCLEANH_ENULL, SFTCLEANH_MSGENULL); 
  ASSERT (harmonicInfo->rightWing, status, SFTCLEANH_ENULL, SFTCLEANH_MSGENULL); 
  ASSERT (fname, status, SFTCLEANH_ENULL, SFTCLEANH_MSGENULL ); 

  /* open harmonics file for reading */
  fp = fopen( fname, "r");
  ASSERT (fname, status, SFTCLEANH_EFILE, SFTCLEANH_MSGEFILE);

  harmonicCount = 0;

  do {
    r=fscanf(fp,"%lf%lf%d%lf%lf%s\n", &temp1, &temp2, 
	     &tempint, &temp3, &temp4, dump);
    /* make sure the line has the right number of entries or is EOF */
    ASSERT( (r==6)||(r==EOF), status, SFTCLEANH_EHEADER, SFTCLEANH_MSGEVAL);
    if (r==6) harmonicCount++;
  } while ( r != EOF);

  harmonicInfo->nHarmonicSets = harmonicCount;

  fclose(fp);

  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);
}


/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/* *******************************  <lalVerbatim file="SFTCleanD"> */
void  LALReadHarmonicsInfo (LALStatus          *status,
			 LineHarmonicsInfo  *harmonicsInfo,
			 CHAR               *fname)
{/*   *********************************************  </lalVerbatim> */
  /* this reads the information about the lines: central frequency, left wing and 
     right wing */
  FILE    *fp = NULL;
  INT4    r, count, nHarmonicSets;
  REAL8   *startFreq=NULL;
  REAL8   *gapFreq=NULL;
  INT4    *numHarmonics=NULL;
  REAL8   *leftWing=NULL;
  REAL8   *rightWing=NULL;
  CHAR    dump[128];

  INITSTATUS (status, "ReadHarmonicsInfo", SFTCLEANC);
  ATTATCHSTATUSPTR (status);  

  /* make sure arguments are not null */
  ASSERT (harmonicsInfo, status, SFTCLEANH_ENULL, SFTCLEANH_MSGENULL); 
  ASSERT (harmonicsInfo->nHarmonicSets > 0, status, SFTCLEANH_EVAL, SFTCLEANH_MSGEVAL);
  ASSERT (harmonicsInfo->startFreq, status, SFTCLEANH_ENULL, SFTCLEANH_MSGENULL);
  ASSERT (harmonicsInfo->gapFreq, status, SFTCLEANH_ENULL, SFTCLEANH_MSGENULL);
  ASSERT (harmonicsInfo->numHarmonics, status, SFTCLEANH_ENULL, SFTCLEANH_MSGENULL);
  ASSERT (harmonicsInfo->leftWing, status, SFTCLEANH_ENULL, SFTCLEANH_MSGENULL); 
  ASSERT (harmonicsInfo->rightWing, status, SFTCLEANH_ENULL, SFTCLEANH_MSGENULL); 
  ASSERT (fname, status, SFTCLEANH_ENULL, SFTCLEANH_MSGENULL);
 
  /* open line noise file for reading */
  fp = fopen( fname, "r");
  ASSERT (fp, status, SFTCLEANH_EFILE,  SFTCLEANH_MSGEFILE);

  nHarmonicSets = harmonicsInfo->nHarmonicSets;
  startFreq = harmonicsInfo->startFreq;
  gapFreq = harmonicsInfo->gapFreq;
  numHarmonics = harmonicsInfo->numHarmonics; 
  leftWing = harmonicsInfo->leftWing;
  rightWing = harmonicsInfo->rightWing;

  /* read line information from file */
  for (count = 0; count < nHarmonicSets; count++){
    r=fscanf(fp,"%lf%lf%d%lf%lf%s\n", startFreq+count, gapFreq+count, numHarmonics+count, 
	     leftWing+count, rightWing+count, dump);
    ASSERT(r==6, status, SFTCLEANH_EHEADER, SFTCLEANH_MSGEVAL);
  }

  fclose(fp);

  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);

}

/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/* *******************************  <lalVerbatim file="SFTCleanD"> */
void  LALHarmonics2Lines (LALStatus          *status,
		       LineNoiseInfo      *lineInfo,
		       LineHarmonicsInfo  *harmonicsInfo)
{/*   *********************************************  </lalVerbatim> */
  /* this reads the information about the lines: central frequency, left wing and 
     right wing */
  
  INT4    count1, count2, nHarmonicSets, maxCount, position;
  REAL8   *startFreq;
  REAL8   *gapFreq;
  INT4    *numHarmonics;
  REAL8   *leftWing;
  REAL8   *rightWing;
  REAL8   f0, deltaf, leftDeltaf, rightDeltaf;


  INITSTATUS (status, "Harmonics2Lines", SFTCLEANC);
  ATTATCHSTATUSPTR (status);  

  /* make sure arguments are not null */
  ASSERT (harmonicsInfo, status, SFTCLEANH_ENULL, SFTCLEANH_MSGENULL); 
  ASSERT (harmonicsInfo->nHarmonicSets > 0, status, SFTCLEANH_EVAL, SFTCLEANH_MSGEVAL);
  ASSERT (harmonicsInfo->startFreq, status, SFTCLEANH_ENULL, SFTCLEANH_MSGENULL);
  ASSERT (harmonicsInfo->gapFreq, status, SFTCLEANH_ENULL, SFTCLEANH_MSGENULL);
  ASSERT (harmonicsInfo->numHarmonics, status, SFTCLEANH_ENULL, SFTCLEANH_MSGENULL);
  ASSERT (harmonicsInfo->leftWing, status, SFTCLEANH_ENULL, SFTCLEANH_MSGENULL); 
  ASSERT (harmonicsInfo->rightWing, status, SFTCLEANH_ENULL, SFTCLEANH_MSGENULL); 

  ASSERT (lineInfo, status, SFTCLEANH_ENULL, SFTCLEANH_MSGENULL); 
  ASSERT (lineInfo->nLines > 0, status, SFTCLEANH_EVAL, SFTCLEANH_MSGEVAL); 
  ASSERT (lineInfo->lineFreq, status, SFTCLEANH_ENULL, SFTCLEANH_MSGENULL); 
  ASSERT (lineInfo->leftWing, status, SFTCLEANH_ENULL, SFTCLEANH_MSGENULL); 
  ASSERT (lineInfo->rightWing, status, SFTCLEANH_ENULL, SFTCLEANH_MSGENULL); 
 
  nHarmonicSets = harmonicsInfo->nHarmonicSets;
  startFreq = harmonicsInfo->startFreq;
  gapFreq = harmonicsInfo->gapFreq;
  numHarmonics = harmonicsInfo->numHarmonics; 
  leftWing = harmonicsInfo->leftWing;
  rightWing = harmonicsInfo->rightWing;

  position = 0;
  for (count1=0; count1 < nHarmonicSets; count1++)
    {
      maxCount = *(numHarmonics + count1);
      f0 = *(startFreq + count1);
      deltaf = *(gapFreq + count1);
      leftDeltaf = *(leftWing + count1);
      rightDeltaf = *(rightWing + count1);
      for (count2 = 0; count2 < maxCount; count2++)
	{
	  *(lineInfo->lineFreq + count2 + position) = f0 + count2 * deltaf;
	  *(lineInfo->leftWing + count2 + position) = leftDeltaf;
	  *(lineInfo->rightWing + count2 + position) = rightDeltaf;
	}
      position += maxCount;
    }


  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);

}




/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
/* *******************************  <lalVerbatim file="SFTCleanD"> */
void LALFindNumberLines (LALStatus          *status,
		      LineNoiseInfo      *lineInfo,            
		      CHAR               *fname)
{/*   *********************************************  </lalVerbatim> */
  /* this function counts the number of lines present in the file "fname" and  
     checks that the format of the lines is correct */

  FILE *fp = NULL;
  REAL8 temp1,temp2,temp3;
  INT4  lineCount, r;
  CHAR  dump[128]; 

  INITSTATUS (status, "FindNumberLines", SFTCLEANC);
  ATTATCHSTATUSPTR (status);  

  /* make sure arguments are not null */
  ASSERT (lineInfo, status, SFTCLEANH_ENULL, SFTCLEANH_MSGENULL); 
  ASSERT (fname, status, SFTCLEANH_ENULL, SFTCLEANH_MSGENULL);

  /* open line noise file for reading */
  fp = fopen( fname, "r");
  ASSERT (fp, status, SFTCLEANH_EFILE,  SFTCLEANH_MSGEFILE);

  lineCount = 0;
  do {
    r=fscanf(fp,"%lf%lf%lf%s\n", &temp1, &temp2, &temp3, dump);
    /* make sure the line has the right number of entries or is EOF */
    ASSERT( (r==4)||(r==EOF), status, SFTCLEANH_EHEADER, SFTCLEANH_MSGEVAL);
    if (r==4) lineCount++;
  } while ( r != EOF);

  lineInfo->nLines = lineCount;

  fclose(fp);

  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);
}

/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/* *******************************  <lalVerbatim file="SFTCleanD"> */
void  LALReadLineInfo (LALStatus     *status,
		    LineNoiseInfo *lineInfo,
		    CHAR          *fname)
{/*   *********************************************  </lalVerbatim> */
  /* this reads the information about the lines: central frequency, left wing and 
     right wing */
  FILE    *fp = NULL;
  INT4    r, count, nLines;
  REAL8   *lineFreq=NULL;
  REAL8   *leftWing=NULL;
  REAL8   *rightWing=NULL;
  CHAR    dump[128];

  INITSTATUS (status, "ReadLineInfo", SFTCLEANC);
  ATTATCHSTATUSPTR (status);  

  /* make sure arguments are not null */
  ASSERT (lineInfo, status, SFTCLEANH_ENULL, SFTCLEANH_MSGENULL); 
  ASSERT (lineInfo->nLines > 0, status, SFTCLEANH_EVAL, SFTCLEANH_MSGEVAL);
  ASSERT (lineInfo->lineFreq, status, SFTCLEANH_ENULL, SFTCLEANH_MSGENULL); 
  ASSERT (lineInfo->leftWing, status, SFTCLEANH_ENULL, SFTCLEANH_MSGENULL); 
  ASSERT (lineInfo->rightWing, status, SFTCLEANH_ENULL, SFTCLEANH_MSGENULL); 
  ASSERT (fname, status, SFTCLEANH_ENULL, SFTCLEANH_MSGENULL);
 
  /* open line noise file for reading */
  fp = fopen( fname, "r");
  ASSERT (fp, status, SFTCLEANH_EFILE,  SFTCLEANH_MSGEFILE);

  nLines = lineInfo->nLines;
  lineFreq = lineInfo->lineFreq;
  leftWing = lineInfo->leftWing;
  rightWing = lineInfo->rightWing;
  /* read line information from file */
  for (count = 0; count < nLines; count++){
    r=fscanf(fp,"%lf%lf%lf%s\n", lineFreq+count, leftWing+count, rightWing+count, dump);
    ASSERT(r==4, status, SFTCLEANH_EHEADER, SFTCLEANH_MSGEVAL);
  }

  fclose(fp);

  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);

}


/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
void LALChooseLines (LALStatus        *status,
		  LineNoiseInfo    *outLine,
		  LineNoiseInfo    *inLine,
		  REAL8            fMin,
		  REAL8            fMax
		  )
{
  /* this function takes a set of lines and a frequency range as input */
  /* it returns a those lines which lie within the specified frequency range */
  INT4 nLinesIn, nLinesOut, j;
  REAL8 lineFreq, leftWing, rightWing;

  INITSTATUS (status, "ChooseLines", SFTCLEANC);
  ATTATCHSTATUSPTR (status);   

  /* some sanity checks */
  ASSERT (outLine, status, SFTCLEANH_ENULL, SFTCLEANH_MSGENULL);
  ASSERT (inLine, status, SFTCLEANH_ENULL, SFTCLEANH_MSGENULL);
  ASSERT (inLine->nLines > 0, status, SFTCLEANH_EVAL, SFTCLEANH_MSGEVAL);
  ASSERT (outLine->nLines > 0, status, SFTCLEANH_EVAL, SFTCLEANH_MSGEVAL);
  ASSERT (inLine->nLines - outLine->nLines == 0, status, SFTCLEANH_EVAL, SFTCLEANH_MSGEVAL);
  ASSERT (fMin > 0, status, SFTCLEANH_EVAL, SFTCLEANH_MSGEVAL);
  ASSERT (fMax > fMin, status, SFTCLEANH_EVAL, SFTCLEANH_MSGEVAL);

  nLinesIn = inLine->nLines;
  nLinesOut = 0;  
  /* loop over lines in inLine structure and see if they are within bounds */
  for(j=0; j<nLinesIn; j++)
    {
      lineFreq = inLine->lineFreq[j];
      leftWing = inLine->leftWing[j];
      rightWing = inLine->rightWing[j];
      if ( (lineFreq >= fMin) && (lineFreq <= fMax))
	{
	  outLine->lineFreq[nLinesOut] = lineFreq;
	  outLine->leftWing[nLinesOut] = leftWing;
	  outLine->rightWing[nLinesOut] = rightWing;
	  nLinesOut++;
	} 
      else if ( (lineFreq < fMin) && (lineFreq + rightWing > fMin) )
	{
	  outLine->lineFreq[nLinesOut] = lineFreq;
	  outLine->leftWing[nLinesOut] = leftWing;
	  outLine->rightWing[nLinesOut] = rightWing;
	  nLinesOut++;
	}
      else if ( (lineFreq > fMax) && (lineFreq - leftWing < fMax) )
	{
	  outLine->lineFreq[nLinesOut] = lineFreq;
	  outLine->leftWing[nLinesOut] = leftWing;
	  outLine->rightWing[nLinesOut] = rightWing;
	  nLinesOut++;
	}
    }

  /* if there are no lines inband then free memory */
  if (nLinesOut == 0)
    {
      LALFree(outLine->lineFreq);
      LALFree(outLine->leftWing);
      LALFree(outLine->rightWing);
    }
  else  /* else reallocate memory for outLine */
    {
      outLine->nLines = nLinesOut;
      outLine->lineFreq = (REAL8 *)LALRealloc(outLine->lineFreq, nLinesOut*sizeof(REAL8));
      outLine->leftWing = (REAL8 *)LALRealloc(outLine->leftWing, nLinesOut*sizeof(REAL8));
      outLine->rightWing = (REAL8 *)LALRealloc(outLine->rightWing, nLinesOut*sizeof(REAL8));
    }

  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);
}

#define TRUE (1==1)
#define FALSE (1==0)

/* function to count how many lines affect a given frequency */
/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
void LALCheckLines ( LALStatus           *status,
		  INT4                *countL,
		  LineNoiseInfo       *lines,
		  REAL8               freq)
{
  INT4 nLines, j;
  REAL8 lineFreq, leftWing, rightWing, doppler;

  INITSTATUS (status, "CheckLines", SFTCLEANC);
  ATTATCHSTATUSPTR (status);   

  /* sanity checks */
  ASSERT (lines, status, SFTCLEANH_ENULL, SFTCLEANH_MSGENULL);
  ASSERT (countL, status, SFTCLEANH_ENULL, SFTCLEANH_MSGENULL);
  ASSERT (lines->nLines > 0, status, SFTCLEANH_EVAL, SFTCLEANH_MSGEVAL);

  /* loop over lines and check if freq is affected by the lines */
  nLines = lines->nLines;
  *countL = 0;
  for (j=0; j<nLines; j++)
    {
      lineFreq = lines->lineFreq[j];
      leftWing = lines->leftWing[j];
      rightWing = lines->rightWing[j];
      /* add doppler band to wings */
      doppler = VTOT * (lineFreq + rightWing);
      leftWing += doppler;
      rightWing += doppler;
      /* now chech if freq lies in range */
      if ( (freq <= lineFreq + rightWing) && (freq >= lineFreq - leftWing))
	*countL += 1;  
    }
      
  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);
}

/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/* *******************************  <lalVerbatim file="SFTCleanD"> */
void LALCleanCOMPLEX8SFT (LALStatus          *status,
		       SFTtype            *sft,
		       INT4               width,       
		       INT4               window,
		       LineNoiseInfo      *lineInfo)
{ /*   *********************************************  </lalVerbatim> */
  /* function to clean the SFT based on the line information read earlier */
  
  INT4     nLines, count, leftCount, rightCount, lineBin, minBin, maxBin, k, tempk;
  INT4     leftWingBins, rightWingBins, length, sumBins;
  REAL8    deltaF, f0, tBase, bias;
  REAL8    stdPow, medianPow;
  REAL8    *tempDataPow=NULL;
  REAL8    *lineFreq=NULL;
  REAL8    *leftWing=NULL;
  REAL8    *rightWing=NULL;
  COMPLEX8 *inData;
  FILE *fp=NULL;   
  INT4 seed, ranCount;  
  RandomParams *randPar=NULL;
  static REAL4Vector *ranVector=NULL; 
  REAL4 *randVal;
  /* --------------------------------------------- */
  INITSTATUS (status, "CleanCOMPLEX8SFT", SFTCLEANC);
  ATTATCHSTATUSPTR (status);   
  
  /*   Make sure the arguments are not NULL: */ 
  ASSERT (sft,   status, SFTCLEANH_ENULL, SFTCLEANH_MSGENULL);
  inData = sft->data->data;
  ASSERT (inData, status, SFTCLEANH_ENULL, SFTCLEANH_MSGENULL);
  ASSERT (lineInfo, status, SFTCLEANH_ENULL, SFTCLEANH_MSGENULL);
  ASSERT (lineInfo->nLines, status, SFTCLEANH_EVAL, SFTCLEANH_MSGEVAL);
  ASSERT (lineInfo->lineFreq, status, SFTCLEANH_ENULL, SFTCLEANH_MSGENULL); 
  ASSERT (lineInfo->leftWing, status, SFTCLEANH_ENULL, SFTCLEANH_MSGENULL); 
  ASSERT (lineInfo->rightWing, status, SFTCLEANH_ENULL, SFTCLEANH_MSGENULL); 
  ASSERT (window > 0, status, SFTCLEANH_EVAL, SFTCLEANH_MSGEVAL);
  ASSERT (width >= 0, status, SFTCLEANH_EVAL, SFTCLEANH_MSGEVAL);
  length = sft->data->length;
  ASSERT (length > 0, status, SFTCLEANH_EHEADER, SFTCLEANH_MSGEVAL);

  /* get the value of RngMedBias from the window size */
  TRY( LALRngMedBias( status->statusPtr, &bias, 2*window ), status ); 
  
  /* copy pointers from input */
  nLines = lineInfo->nLines;
  lineFreq = lineInfo->lineFreq;
  leftWing = lineInfo->leftWing;
  rightWing = lineInfo->rightWing;
  
  deltaF = sft->deltaF;
  tBase = 1.0/deltaF;
  f0 = sft->f0;
  minBin = floor(f0/deltaF + 0.5);
  maxBin = minBin + length - 1;

  /* allocate memory for storing sft power */
  tempDataPow = LALMalloc(2*window*sizeof(REAL8));
 
  fp=fopen("/dev/urandom", "r");
  ASSERT(fp, status, SFTCLEANH_EFILE, SFTCLEANH_MSGEFILE);  

  ranCount = fread(&seed, sizeof(seed), 1, fp);
  ASSERT(ranCount==1, status, SFTCLEANH_EREAD, SFTCLEANH_MSGEREAD);  

  fclose(fp);

  /* calculate total number of bins to see how many random numbers must be generated */
  sumBins = 0;
  for (count = 0; count < nLines; count++)
    {
      { 
	INT4 tempSumBins;
	tempSumBins = floor(tBase*leftWing[count]) + floor(tBase*rightWing[count]);  
	sumBins += tempSumBins < 2*width ? tempSumBins : 2*width;
      } 
    }

  TRY ( LALCreateRandomParams (status->statusPtr, &randPar, seed), status);
  TRY ( LALCreateVector (status->statusPtr, &ranVector, 2*(sumBins + nLines)), status);
  TRY ( LALNormalDeviates (status->statusPtr, ranVector, randPar), status);
  TRY ( LALDestroyRandomParams (status->statusPtr, &randPar), status);  

  tempk = 0;  
  /* loop over the lines */
  for (count = 0; count < nLines; count++){
    
    /* find frequency bins for line frequency and wings */
    lineBin = floor(tBase * lineFreq[count] + 0.5);
    leftWingBins = floor(tBase * leftWing[count]);
    rightWingBins = floor(tBase * rightWing[count]);
    
    /* check that frequency is within band of sft and line is not too wide*/
    if ((lineBin >= minBin) && (lineBin <= maxBin) && (leftWingBins <= width) && (rightWingBins <= width)){

      /* estimate the sft power in "window" # of bins each side */
      for (k = 0; k < window ; k++){
	if (maxBin - lineBin - rightWingBins - k > 0)
	  inData = sft->data->data + lineBin - minBin + rightWingBins + k + 1;
	else
	  inData = sft->data->data + length - 1;

	tempDataPow[k] = (inData->re)*(inData->re) + (inData->im)*(inData->im);
	
	if (lineBin - minBin -leftWingBins - k > 0)
	  inData = sft->data->data + lineBin - minBin - leftWingBins - k - 1;
	else 
	  inData = sft->data->data;

	tempDataPow[k+window] = (inData->re)*(inData->re) + (inData->im)*(inData->im);   
      }
    
      gsl_sort( tempDataPow, 1, 2*window);
      medianPow = gsl_stats_median_from_sorted_data(tempDataPow, 1, 2*window);
      stdPow = sqrt(medianPow/(2 * bias));
      
      /* set sft value at central frequency to noise */
      inData = sft->data->data + lineBin - minBin;

      randVal = ranVector->data + tempk;  
      inData->re = stdPow * (*randVal); 
      tempk++;

      randVal = ranVector->data + tempk;  
      inData->im = stdPow * (*randVal); 
      tempk++;
    
      /* now go left and set the left wing to noise */
      /* make sure that we are always within the sft band */
      /* and set bins to zero only if Wing width is smaller than "width" */ 
      if ((leftWingBins <= width)){
	for (leftCount = 0; leftCount < leftWingBins; leftCount++){
	  if ( (lineBin - minBin - leftCount > 0)){
	    inData = sft->data->data + lineBin - minBin - leftCount - 1;

	    randVal = ranVector->data + tempk;  
	    inData->re = stdPow * (*randVal); 
	    tempk++;

	    randVal = ranVector->data + tempk;  
	    inData->im = stdPow * (*randVal); 
	    tempk++;
	  }
	}
      }
      
      /* now go right making sure again to stay within the sft band */
      if ((rightWingBins <= width )){
	for (rightCount = 0; rightCount < rightWingBins; rightCount++){
	  if ( (maxBin - lineBin - rightCount > 0)){
	    inData = sft->data->data + lineBin - minBin + rightCount + 1;

	    randVal = ranVector->data + tempk;  
	    inData->re = stdPow * (*randVal); 
            tempk++;

	    randVal = ranVector->data + tempk;  
	    inData->im = stdPow * (*randVal);
	    tempk++;
	  }
	}
      }
    }
  } /* end loop over lines */

  /* free memory */
  LALFree(tempDataPow);
  TRY (LALDestroyVector (status->statusPtr, &ranVector), status);

  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);
}









