/*
 *  Copyright (C) 2005 Badri Krishnan, Alicia Sintes, Greg Mendell  
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
 * \author Badri Krishnan, Alicia Sintes, Greg Mendell
 * \file SFTClean.c
 * \brief Module containing routines for dealing with spectral disturbances in SFTs

   \par Description 
   
   This module contains routines for dealing with lists of known spectral disturbances 
   in the frequency domain, and using them to clean SFTs.  

   The basic input is a text file containing a list of known spectral lines.  An example 
   is the following

   \verbatim
   0.0      0.25     4000     0.0        0.0   0.25Hzlines  
   0.0      60.0     20       0.5        0.5   60Hzlines     
   0.0      16.0     100      0.0        0.0   16Hzlines   
   166.7    0.0      1        0.0        0.0   Calibrationline 
   345.0    0.0      1        3.0        3.0   violinmodes   
   \endverbatim

   The file consists of rows with 6 columns each.  Each row has information about
   a set of spectral lines of the form \f$ f_n = f_0 + n\Delta f \f$.  The first column 
   is the start frequency \f$ f_0 \f$, the second column is the spacing \f$ \Delta f \f$,
   the third column is the total number of lines, the fourth column is the 
   left-width of each line (in Hz), the fifth column is the width on the right, and
   the last column is a brief comment string (no spaces).  If this file is meant to
   be used for cleaning SFTs, then certain features which the user must be aware of
   are explained in the documentation of the function LALCleanCOMPLEX8SFT().  

   \par Uses
   \code
   void LALFindNumberHarmonics (LALStatus, LineHarmonicsInfo, CHAR)
   void LALReadHarmonicsInfo (LALStatus, LineHarmonicsInfo, CHAR)
   void LALHarmonics2Lines (LALStatus, LineNoiseInfo, LineHarmonicsInfo)
   void LALChooseLines (LALStatus, LineNoiseInfo, LineNoiseInfo, REAL8, REAL8)
   void LALCheckLines ( LALStatus, INT4, LineNoiseInfo, REAL8)
   void LALFindNumberLines (LALStatus, LineNoiseInfo, CHAR)
   void LALReadLineInfo (LALStatus, LineNoiseInfo, CHAR)
   void LALCleanCOMPLEX8SFT (LALStatus, SFTtype, INT4, INT4, LineNoiseInfo, RandomParams)
   \endcode


*/

/************************************ <lalVerbatim file="SFTCleanCV">
Author: Sintes, A.M., Krishnan, B.
$Id$
************************************* </lalVerbatim> */

/* REVISIONS: */
/* 09/09/05 gam; if (nLinesOut == 0) still need outLine->nLines = nLinesOut; calling function needs to know this */
/* 09/09/05 gam; make RandomParams *randPar a parameter for CleanCOMPLEX8SFT. Thus only need to */
/*               initialze RandomParams *randPar once and avoid repeatly opening /dev/urandom.  */
/* 09/09/05 gam; prefix function names with LAL in init status macros */
/* 09/09/05 gam; only assert harmonicInfo and fname in LALFindNumberHarmonic and fix assert of fp. */
/*               Other pointers can be unititialized until nHarmonicSets is determined.            */

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



/**
 * Looks into the input file containing list of lines, does some checks on the 
 * file format, and calculates the number of harmonic sets in this file.  
 */
/*<lalVerbatim file="SFTCleanD"> */
void LALFindNumberHarmonics (LALStatus    *status,
			     LineHarmonicsInfo   *harmonicInfo, /**< list of harmonics */
			     CHAR         *fname /**< input filename */)
{/*   *********************************************  </lalVerbatim> */

  FILE *fp = NULL;
  CHAR  dump[128];
  INT4   harmonicCount, r, tempint; 
  REAL8 temp1, temp2, temp3, temp4;   


  INITSTATUS (status, "LALFindNumberHarmonics", SFTCLEANC);
  ATTATCHSTATUSPTR (status);

  /* make sure arguments are not null */
  ASSERT (harmonicInfo, status, SFTCLEANH_ENULL, SFTCLEANH_MSGENULL ); 
  /* ASSERT (harmonicInfo->nHarmonicSets > 0, status, SFTCLEANH_EVAL, SFTCLEANH_MSGEVAL); 
  ASSERT (harmonicInfo->startFreq, status, SFTCLEANH_ENULL, SFTCLEANH_MSGENULL);
  ASSERT (harmonicInfo->gapFreq, status, SFTCLEANH_ENULL, SFTCLEANH_MSGENULL);
  ASSERT (harmonicInfo->numHarmonics, status, SFTCLEANH_ENULL, SFTCLEANH_MSGENULL);
  ASSERT (harmonicInfo->leftWing, status, SFTCLEANH_ENULL, SFTCLEANH_MSGENULL); 
  ASSERT (harmonicInfo->rightWing, status, SFTCLEANH_ENULL, SFTCLEANH_MSGENULL); */ /* 09/09/05 gam */
  ASSERT (fname, status, SFTCLEANH_ENULL, SFTCLEANH_MSGENULL ); 

  /* open harmonics file for reading */
  fp = fopen( fname, "r");
  /* ASSERT (fname, status, SFTCLEANH_EFILE, SFTCLEANH_MSGEFILE); */ /* 09/09/05 gam */
  ASSERT (fp, status, SFTCLEANH_EFILE, SFTCLEANH_MSGEFILE);

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


/** Reads in the contents of the input line-info file and fills up
 *  the LineHarmonicsInfo structure.  Appropriate memory must be allocated for
 *  this structure before this function is called.  
 */
/*  <lalVerbatim file="SFTCleanD"> */
void  LALReadHarmonicsInfo (LALStatus          *status,
			    LineHarmonicsInfo  *harmonicsInfo, /**< list of harmonics */
			    CHAR               *fname /**< input file */)  
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

  INITSTATUS (status, "LALReadHarmonicsInfo", SFTCLEANC);
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

/** Converts the list of harmonic sets into an explicit list of spectral 
 *  lines.
 */
/*  <lalVerbatim file="SFTCleanD"> */
void  LALHarmonics2Lines (LALStatus          *status,
			  LineNoiseInfo      *lineInfo, /**< output list of explicit lines */
			  LineHarmonicsInfo  *harmonicsInfo) /**< input list of harmonics */
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


  INITSTATUS (status, "LALHarmonics2Lines", SFTCLEANC);
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




/** 
 * Finds total number of spectral-lines contained in case the input file is 
 * a list of explicit spectral lines -- obsolete.  
 * Use instead LALFindNumberHarmonics(). 
 */
/*  <lalVerbatim file="SFTCleanD"> */
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

  INITSTATUS (status, "LALFindNumberLines", SFTCLEANC);
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

/** Reads information from file containing list of explicit lines - obsolete. 
 *  Use instead LALReadHarmonicsInfo()  
 */
/*  <lalVerbatim file="SFTCleanD"> */
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

  INITSTATUS (status, "LALReadLineInfo", SFTCLEANC);
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


/** Takes a set of spectral lines and a frequency range as input and 
 * returns those lines which lie within the specified frequency range.  The output
 * is a reduced list of spectral lines which either lie completely within the 
 * frequency range or whose wings overlap with the frequency range.  This is useful
 * for discarding unnecessary lines to save computational cost and memory. 
 */
/*  <lalVerbatim file="SFTCleanD"> */
void LALChooseLines (LALStatus        *status,
		     LineNoiseInfo    *outLine,  /**< reduced list of lines */
		     LineNoiseInfo    *inLine,   /**< input list of lines */
		     REAL8            fMin,      /**< start of frequency band */
		     REAL8            fMax       /**< end of frequency band */
		  )
{/*   *********************************************  </lalVerbatim> */

  INT4 nLinesIn, nLinesOut, j;
  REAL8 lineFreq, leftWing, rightWing;

  INITSTATUS (status, "LALChooseLines", SFTCLEANC);
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
      outLine->nLines = nLinesOut; /* 09/09/05 gam; calling function needs to know this. */
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

/** Function to count how many lines affect a given frequency.  Input is a 
 * list of lines and a frequency.  The output is an integer which is equal to the 
 * number of lines which intersect this frequency.  An output of zero indicates
 * that the frequencty is not influenced by the lines.  Note that the doppler
 * broadening of the lines is taken into account while deciding whether the 
 * frequency is affected or not. 
 */
/*  <lalVerbatim file="SFTCleanD"> */
void LALCheckLines ( LALStatus           *status,
		     INT4                *countL, /**< number of lines affecting frequency */
		     LineNoiseInfo       *lines, /**< list of lines */
		     REAL8               freq)   /**< frequency to be checked */
{/*   *********************************************  </lalVerbatim> */

  INT4 nLines, j;
  REAL8 lineFreq, leftWing, rightWing, doppler;

  INITSTATUS (status, "LALCheckLines", SFTCLEANC);
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


/**
 * Function for cleaning a SFT given a set of known spectral disturbances.  
 * The algorithm is the following.  For each 
 * spectral line, first the central frequency of the line is converted to a 
 * bin index by floor(tBase * freq + 0.5).  If the wings are set to zero, then only this 
 * central bin is cleaned.  Note that if the frequency lies between two exactly 
 * resolved frequencies, then only one of these bins is cleaned.  The user must 
 * know about the SFT timebase and make sure that the central frequency is an
 * exactly resolved one.  The wingsize is calculated in bins according to 
 * floor(tBase * wingsize).  This is done separately for the left and right wings.  
 * Note the use of the floor function.  If the wingsize corresponds to say 2.5 bins, then
 * only 2 bins will be cleaned in addition to the central frequency.  
 *
 * Having decided which bins are to be cleaned, the next step is to produce random noise 
 * to replace the data in these bins. The fake random noise must mimic the 
 * behavior of the true noise in the vicinity of the spectral disturbance, and we must 
 * therefore estimate the noise floor in the vicinity of the disturbance.  
 * The user specifies a "window size" for cleaning,  and this determines how many data 
 * points from the SFT are to be used for estimating this noise floor.  
 * Consider the number of frequency bins on each side of the spectral disturbance  
 * (excluding the wings) given by the user specified window size.  The function calculates 
 * the SFT power in these bins and calculates their median.  The median is then converted
 * to a standard deviation for the real and imaginary parts of the SFT data, 
 * taking into account the median bias.  There is also a width parameter.  This is 
 * an upper limit on the number of bins that will be cleaned on either wing. Thus, if
 * width is specified to say 3, then only a maximum 3 bins will be cleaned on either side
 * irrespective of the actual wing size specified.  This parameter is present 
 * for historical reasons, and should probably be removed.  Currently, it is recommended
 * to set this sufficiently large (larger that any wing size in bins) so that 
 * it has no effect.  

 \par Some "features" that must be kept in mind
 
 The following is part of a email sent to pulgroup on 4/10/05.  Some of these points
 have already been mentioned above.  

 i) If the width of the line to be cleaned is specified to be zero and if the central frequency is midway between two exactly resolved bins, then the code
 will only clean one frequency bin and not both as it perhaps should.
 
 ii) The wings size is converted from Hz to frequency bins using the floor function. So some bins are being dropped at the edges due to rounding off.
 
 iii) The cleaning uses data from neighboring bins to generate random noise.  This is a problem if there are other spectral disturbances in the
 neighborhood or if the frequency bin is at the edge of the SFT.  Shouldn't one make sure that the data used to generate the random number are
 un-contaminated?
 
                                                                                                                                                             
 Regarding the first two points, the user is supposed to know the properties of the SFTs which are being cleaned and the list of lines is meant to be
 tailored for a particular set of SFTs. Therefore, the user should know the timebaseline of the SFT to make sure that the central frequency value being
 specified is a resolved frequency bin and
 the wings should be specified correctly.  In many cases, the size of the wings, when it is non-zero, is somewhat uncertain, and this uncertainty is often
 larger that this rounding off error.
 
 Perhaps one should specify the frequency bin index of the central frequency value and the wing size also in bins, but this makes the code more cumbersome
 and non-intuitive to use.
 
 Both of these points can be handled by documenting the code better, so that the user knows exactly what is being done and chooses the input
 appropriately .  I will do this as soon as possible.
 
 The third point is more difficult.  Currently, the code calculates the *median* of the power in the neighboring bins and uses that to generate a random
 number with the appropriate standard deviation.
 
 Let us first consider the cases when there are other spectral disturbances in the neighborhood so that the median estimate might be corrupted.  When
 there are only a small number of narrow (only few bins) disturbances in
 the neighborhood, the use of the median is meant to eliminate the effects of these lines.  However, this might not work for the cases when there is a
 broad spectral disturbance or when there are a large number of narrow disturbances.
 
 If there are a large number of narrow spectral disturbances very close to each other, it might make more sense to group them together and consider them
 to be one single broad disturbance; we probably won't trust a detection near these lines anyway.
 
 In the case of a broad spectral feature, it will certainly corrupt the median value and the random number generated won't be correct.  The alternative is
 to skip the broad feature and take only undisturbed bins further away. This might be ok, but recall that the purpose of the cleaning is to replace a
 spectral disturbance with a random number whose distribution is similar to the noise in the *neighborhood* of the disturbance.  For colored noise, if we
 have to go very far away in frequency to get a undisturbed frequency bin, then the random number distribution won't be correct either.
 Edge effects are taken care by making sure that the SFT being read in is large enough. There
 must be extra wings to the SFT equal to the window size used for cleaning.  If the edge effects are important, then the code using the cleaning routines
 must read in the extra bins.  
   
 */
/* *******************************  <lalVerbatim file="SFTCleanD"> */
void LALCleanCOMPLEX8SFT (LALStatus          *status,
			  SFTtype            *sft,  /**< SFT to be cleaned */
			  INT4               width, /**< maximum width to be cleaned -- set sufficiently large if all bins in each line are to be cleaned*/      
			  INT4               window,/**< window size for noise floor estimation in vicinity of a line */
			  LineNoiseInfo      *lineInfo, /**< list of lines */ 
			  RandomParams       *randPar /**< parameters for generating random noise */)
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
  /* FILE *fp=NULL;   
  INT4 seed, ranCount;  
  RandomParams *randPar=NULL; */ /* 09/09/05 gam; randPar now a function argument */
  static REAL4Vector *ranVector=NULL; 
  REAL4 *randVal;
  /* --------------------------------------------- */
  INITSTATUS (status, "LALCleanCOMPLEX8SFT", SFTCLEANC);
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
  
  /* 09/09/05 gam; randPar now a function argument */
  /* fp=fopen("/dev/urandom", "r");
  ASSERT(fp, status, SFTCLEANH_EFILE, SFTCLEANH_MSGEFILE);  

  ranCount = fread(&seed, sizeof(seed), 1, fp);
  ASSERT(ranCount==1, status, SFTCLEANH_EREAD, SFTCLEANH_MSGEREAD);  

  fclose(fp); */

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

  /* TRY ( LALCreateRandomParams (status->statusPtr, &randPar, seed), status); */ /* 09/09/05 gam; randPar now a function argument */
  TRY ( LALCreateVector (status->statusPtr, &ranVector, 2*(sumBins + nLines)), status);
  TRY ( LALNormalDeviates (status->statusPtr, ranVector, randPar), status);
  /* TRY ( LALDestroyRandomParams (status->statusPtr, &randPar), status);  */ /* 09/09/05 gam; randPar now a function argument */

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




/** Function to clean a sft vector -- calls LALCleanCOMPLEX8SFT repeatedly for each 
    sft in vector */
void LALCleanSFTVector (LALStatus       *status,
			SFTVector       *sftVect,  /**< SFTVector to be cleaned */
			INT4            width,     /**< maximum width to be cleaned -- set sufficiently large if all bins in each line are to be cleaned*/      
			INT4            window,    /**< window size for noise floor estimation in vicinity of a line */
			LineNoiseInfo   *lineInfo, /**< list of lines */ 
			RandomParams    *randPar   /**< parameters for generating random noise */)
{

  UINT4 k;

  INITSTATUS (status, "LALCleanSFTVector", SFTCLEANC);
  ATTATCHSTATUSPTR (status);   
  
  ASSERT (sftVect, status, SFTCLEANH_ENULL, SFTCLEANH_MSGENULL);  
  ASSERT (randPar, status, SFTCLEANH_ENULL, SFTCLEANH_MSGENULL);  
  ASSERT (sftVect->data, status, SFTCLEANH_ENULL, SFTCLEANH_MSGENULL);  
  ASSERT (sftVect->length > 0, status, SFTCLEANH_EVAL, SFTCLEANH_MSGEVAL);  
  ASSERT (lineInfo, status, SFTCLEANH_ENULL, SFTCLEANH_MSGENULL);
  ASSERT (lineInfo->nLines, status, SFTCLEANH_EVAL, SFTCLEANH_MSGEVAL);
  ASSERT (lineInfo->lineFreq, status, SFTCLEANH_ENULL, SFTCLEANH_MSGENULL); 
  ASSERT (lineInfo->leftWing, status, SFTCLEANH_ENULL, SFTCLEANH_MSGENULL); 
  ASSERT (lineInfo->rightWing, status, SFTCLEANH_ENULL, SFTCLEANH_MSGENULL); 
  ASSERT (window > 0, status, SFTCLEANH_EVAL, SFTCLEANH_MSGEVAL);
  ASSERT (width >= 0, status, SFTCLEANH_EVAL, SFTCLEANH_MSGEVAL);

  for ( k = 0; k < sftVect->length; k++) {
    TRY (LALCleanCOMPLEX8SFT (status->statusPtr, sftVect->data + k, width, window, lineInfo, randPar), status);    
  }

  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);
}


/** Function to clean a sft vector -- calls LALCleanCOMPLEX8SFT repeatedly for each 
    sft in vector */
void LALCleanMultiSFTVect (LALStatus       *status,
			   MultiSFTVector  *multVect, /**< SFTVector to be cleaned */
			   INT4            width,     /**< maximum width to be cleaned -- set sufficiently large if all bins in each line are to be cleaned*/      
			   INT4            window,    /**< window size for noise floor estimation in vicinity of a line */
			   LineNoiseInfo   *lineInfo, /**< list of lines */ 
			   RandomParams    *randPar   /**< parameters for generating random noise */)
{

  UINT4 k;

  INITSTATUS (status, "LALCleanMultiSFTVector", SFTCLEANC);
  ATTATCHSTATUSPTR (status);   
  
  ASSERT (multVect, status, SFTCLEANH_ENULL, SFTCLEANH_MSGENULL);  
  ASSERT (multVect->data, status, SFTCLEANH_ENULL, SFTCLEANH_MSGENULL);  
  ASSERT (multVect->length > 0, status, SFTCLEANH_EVAL, SFTCLEANH_MSGEVAL);  
  ASSERT (lineInfo, status, SFTCLEANH_ENULL, SFTCLEANH_MSGENULL);
  ASSERT (lineInfo->nLines, status, SFTCLEANH_EVAL, SFTCLEANH_MSGEVAL);
  ASSERT (lineInfo->lineFreq, status, SFTCLEANH_ENULL, SFTCLEANH_MSGENULL); 
  ASSERT (lineInfo->leftWing, status, SFTCLEANH_ENULL, SFTCLEANH_MSGENULL); 
  ASSERT (lineInfo->rightWing, status, SFTCLEANH_ENULL, SFTCLEANH_MSGENULL); 
  ASSERT (window > 0, status, SFTCLEANH_EVAL, SFTCLEANH_MSGEVAL);
  ASSERT (width >= 0, status, SFTCLEANH_EVAL, SFTCLEANH_MSGEVAL);


  for ( k = 0; k < multVect->length; k++) {
    TRY (LALCleanSFTVector (status->statusPtr, multVect->data[k], width, window, lineInfo, randPar), status);    
  }

  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);
}


/** function to remove lines from a sft vector given a file 
    containing list of lines */
void LALRemoveKnownLinesInSFTVect (LALStatus   *status,
				   SFTVector   *sftVect,  /**< SFTVector to be cleaned */
				   INT4        width,     /**< maximum width to be cleaned -- set sufficiently large if all bins in each line are to be cleaned*/      
				   INT4        window,    /**< window size for noise floor estimation in vicinity of a line */
				   CHAR        *linefile, /**< file with list of lines */
				   RandomParams *randPar) /**< for creating random numbers */ 
{


  static LineNoiseInfo   lines, lines2;
  static LineHarmonicsInfo harmonics; 
  INT4 nLines=0, count1, nHarmonicSets;

  REAL8 f_min, f_max, deltaF;
  UINT4 numBins;

  INITSTATUS (status, "LALRemoveKnownLinesInSFTVector", SFTCLEANC);
  ATTATCHSTATUSPTR (status);   


  ASSERT (sftVect, status, SFTCLEANH_ENULL, SFTCLEANH_MSGENULL);  
  ASSERT (sftVect->data, status, SFTCLEANH_ENULL, SFTCLEANH_MSGENULL);  
  ASSERT (sftVect->length > 0, status, SFTCLEANH_EVAL, SFTCLEANH_MSGEVAL);  
  ASSERT (linefile, status, SFTCLEANH_ENULL, SFTCLEANH_MSGENULL);  
  ASSERT (width > 0, status, SFTCLEANH_EVAL, SFTCLEANH_MSGEVAL);  
  ASSERT (window > 0, status, SFTCLEANH_EVAL, SFTCLEANH_MSGEVAL);  

  f_min = sftVect->data[0].f0;
  deltaF = sftVect->data->deltaF;
  numBins = sftVect->data->data->length;
  f_max = f_min + deltaF * numBins;

  TRY( LALFindNumberHarmonics (status->statusPtr, &harmonics, linefile), status); 
  nHarmonicSets = harmonics.nHarmonicSets; 

  if (nHarmonicSets > 0) /* nothing to do otherwise */
    {
      harmonics.startFreq = (REAL8 *)LALMalloc(harmonics.nHarmonicSets * sizeof(REAL8));
      harmonics.gapFreq = (REAL8 *)LALMalloc(harmonics.nHarmonicSets * sizeof(REAL8));
      harmonics.numHarmonics = (INT4 *)LALMalloc(harmonics.nHarmonicSets * sizeof(INT4));
      harmonics.leftWing = (REAL8 *)LALMalloc(harmonics.nHarmonicSets * sizeof(REAL8));
      harmonics.rightWing = (REAL8 *)LALMalloc(harmonics.nHarmonicSets * sizeof(REAL8));
    

      TRY( LALReadHarmonicsInfo( status->statusPtr, &harmonics, linefile ), status);
      
      nLines = 0;
      for (count1=0; count1 < nHarmonicSets; count1++)
	{
	  nLines += *(harmonics.numHarmonics + count1);
	}
      
      lines.nLines = nLines;
      lines.lineFreq = (REAL8 *)LALMalloc(nLines * sizeof(REAL8));
      lines.leftWing = (REAL8 *)LALMalloc(nLines * sizeof(REAL8));
      lines.rightWing = (REAL8 *)LALMalloc(nLines * sizeof(REAL8));

      TRY( LALHarmonics2Lines( status->statusPtr, &lines, &harmonics), status);


      lines2.nLines = nLines;
      lines2.lineFreq = (REAL8 *)LALMalloc(nLines * sizeof(REAL8));
      lines2.leftWing = (REAL8 *)LALMalloc(nLines * sizeof(REAL8));
      lines2.rightWing = (REAL8 *)LALMalloc(nLines * sizeof(REAL8));

      TRY( LALChooseLines( status->statusPtr, &lines2, &lines, f_min, f_max), status);
      nLines = lines2.nLines;
      
      /* clean the sft vector if there were any lines between f_min and f_max*/
      if ( nLines > 0 ) {
	TRY ( LALCleanSFTVector( status->statusPtr, sftVect, width, window, &lines2, randPar), status);

	/* if nLines2 == 0 then it is freed inside LALChooseLines -- change this? */
	LALFree(lines2.lineFreq);
	LALFree(lines2.leftWing);
	LALFree(lines2.rightWing);

      }
      
      /* free memory */
      LALFree(lines.lineFreq);
      LALFree(lines.leftWing);
      LALFree(lines.rightWing);

      LALFree(harmonics.startFreq);
      LALFree(harmonics.gapFreq);
      LALFree(harmonics.numHarmonics);
      LALFree(harmonics.leftWing);
      LALFree(harmonics.rightWing);

    } /* if nHarmonicSets > 0 */


  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);
}

/** top level function to remove lines from a multi sft vector given a list of files 
    containing list of known spectral lines */
void LALRemoveKnownLinesInMultiSFTVector (LALStatus        *status,
					  MultiSFTVector   *MultiSFTVect,  /**< SFTVector to be cleaned */
					  INT4             width,          /**< maximum width to be cleaned */      
					  INT4             window,         /**< window size for noise floor estimation in vicinity of a line */
					  LALStringVector *linefiles,      /**< file with list of lines */
					  RandomParams     *randPar)       /**< for creating random numbers */ 
{

  UINT4 k, j, numifos;
  CHAR *ifo;

  INITSTATUS (status, "LALRemoveKnownLinesInMultiSFTVector", SFTCLEANC);
  ATTATCHSTATUSPTR (status);   


  ASSERT (MultiSFTVect, status, SFTCLEANH_ENULL, SFTCLEANH_MSGENULL);  
  ASSERT (MultiSFTVect->data, status, SFTCLEANH_ENULL, SFTCLEANH_MSGENULL);  
  ASSERT (MultiSFTVect->length > 0, status, SFTCLEANH_EVAL, SFTCLEANH_MSGEVAL);  
  ASSERT (width > 0, status, SFTCLEANH_EVAL, SFTCLEANH_MSGEVAL);  
  ASSERT (window > 0, status, SFTCLEANH_EVAL, SFTCLEANH_MSGEVAL);  

  numifos = MultiSFTVect->length;

  if ( linefiles != NULL ) {

  ASSERT (linefiles->length > 0, status, SFTCLEANH_EVAL, SFTCLEANH_MSGEVAL);      
  ASSERT (linefiles->data, status, SFTCLEANH_ENULL, SFTCLEANH_MSGENULL);      

  /* loop over linefiles and clean the relevant SFTs */
  for ( k = 0; k < linefiles->length; k++) 
    { 
      ifo = NULL;
      /* try to get the ifo name from the linefile name */
      if ( (ifo = XLALGetChannelPrefix ( linefiles->data[k])) == NULL) { 
        ABORT ( status, SFTCLEANH_ELINENAME, SFTCLEANH_MSGELINENAME);  
      }

      /* loop over ifos and see if any matches */
      for ( j = 0; j < numifos; j++)
	{
	  if ( strncmp( ifo, MultiSFTVect->data[j]->data->name, 3) == 0) {
	    /* clean the sftvector which has matched */
	    TRY ( LALRemoveKnownLinesInSFTVect ( status->statusPtr, MultiSFTVect->data[j], 
						 width, window, linefiles->data[k], randPar), status);  
	  } 

	} /* loop over ifos */

      LALFree( ifo );

    } /* loop over linefiles */
  
  } /* if linefiles != NULL */


  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);
}

