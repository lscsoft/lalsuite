/*-----------------------------------------------------------------------
 *
 * File Name: SFTfileIO.c
 *
 * Authors: Sintes, A.M.,  Krishnan, B. & inspired from Siemens, X.
 *
 * Revision: $Id$
 *
 * History:   Created by Sintes May 21, 2003
 *            Modified by Krishnan on Feb 22, 2004
 *
 *-----------------------------------------------------------------------
 */

/************************************ <lalVerbatim file="SFTfileIOCV">
Author: Sintes, A.M., Krishnan, B., Prix, R.
$Id$
************************************* </lalVerbatim> */

/* <lalLaTeX>  *******************************************************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\section{MOdule \texttt{SFTfileIO.c}}
\label{ss:SFTfileIO.c}
Routines for reading SFT binary files

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Prototypes}
\idx{LALReadSFTheader()}
\idx{LALReadSFTtype()}
\idx{LALReadSFTfile()}
\idx{LALWriteSFTtoFile()}
\vspace{0.1in}
\input{SFTfileIOD}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Description}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Uses}
\begin{verbatim}
LALHO()
\end{verbatim}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Notes}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\vfill{\footnotesize\input{SFTfileIOCV}}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

*********************************************** </lalLaTeX> */

#include <lal/SFTfileIO.h>

NRCSID (SFTFILEIOC, "$Id$");


static void endian_swap(CHAR * pdata, size_t dsize, size_t nelements);
static void ByteSwapSFTHeader (LALStatus *stat, SFTHeader *header);

/*
 * The functions that make up the guts of this module
 */


/*----------------------------------------------------------------------*/

/* *******************************  <lalVerbatim file="SFTfileIOD"> */
void LALReadSFTheader (LALStatus  *status,
                       SFTHeader   *header,
		       const CHAR  *fname)
{ /*   *********************************************  </lalVerbatim> */
  
  FILE     *fp = NULL;
  SFTHeader  header1;
  size_t    errorcode;
  
  /* --------------------------------------------- */
  INITSTATUS (status, "LALReadSFTHeader", SFTFILEIOC);
  ATTATCHSTATUSPTR (status); 
  
  /*   Make sure the arguments are not NULL: */ 
  ASSERT (header, status, SFTFILEIOH_ENULL,  SFTFILEIOH_MSGENULL);
  ASSERT (fname,  status, SFTFILEIOH_ENULL,  SFTFILEIOH_MSGENULL);
  
  /* opening the SFT binary file */
  fp = fopen( fname, "r");
  ASSERT (fp, status, SFTFILEIOH_EFILE,  SFTFILEIOH_MSGEFILE);
  
  /* read the header */
  errorcode = fread( (void *)&header1, sizeof(SFTHeader), 1, fp);
  ASSERT (errorcode ==1, status, SFTFILEIOH_EHEADER,  SFTFILEIOH_MSGEHEADER);

  if (header1.version != 1) {
    ABORT (status, SFTFILEIOH_EVERSION, SFTFILEIOH_MSGEVERSION);
  }

  header->version     = header1.version;
  header->gpsSeconds = header1.gpsSeconds;
  header->gpsNanoSeconds = header1.gpsNanoSeconds;
  header->timeBase       = header1.timeBase;
  header->fminBinIndex   = header1.fminBinIndex ;
  header->length     = header1.length ;
  
  fclose(fp);
  
  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);

} /* LALReadSFTheader() */

/* ------------------------------------------------------------*/
/* *******************************  <lalVerbatim file="SFTfileIOD"> */
void LALReadSFTtype(LALStatus  *status,
		    SFTtype    *sft,    /* asumed  memory is allocated  */
		    const CHAR       *fname,
		    INT4       fminBinIndex)
{ /*   *********************************************  </lalVerbatim> */

  FILE        *fp = NULL;
  SFTHeader  header1;
  size_t     errorcode;
  INT4       offset;
  INT4       length;
  REAL8      deltaF,timeBase;
  
  /* --------------------------------------------- */
  INITSTATUS (status, "LALReadSFTtype", SFTFILEIOC);
  ATTATCHSTATUSPTR (status); 
  
  /*   Make sure the arguments are not NULL: */ 
  ASSERT (sft,   status, SFTFILEIOH_ENULL, SFTFILEIOH_MSGENULL);
  ASSERT (sft->data, status, SFTFILEIOH_ENULL, SFTFILEIOH_MSGENULL);
  ASSERT (fname, status, SFTFILEIOH_ENULL, SFTFILEIOH_MSGENULL);
  
  /* opening the SFT binary file */
  fp = fopen( fname, "r");
  ASSERT (fp, status, SFTFILEIOH_EFILE,  SFTFILEIOH_MSGEFILE);
  
  /* read the header */
  errorcode = fread( (void *)&header1, sizeof(SFTHeader), 1, fp);
  ASSERT (errorcode == 1,      status, SFTFILEIOH_EHEADER, SFTFILEIOH_MSGEHEADER);
  if (header1.version != 1) {
    ABORT (status, SFTFILEIOH_EVERSION, SFTFILEIOH_MSGEVERSION);
  }
  
  /* check that the data we want is in the file and it is correct */
  ASSERT (header1.timeBase > 0.0, status, SFTFILEIOH_EVAL, SFTFILEIOH_MSGEVAL);

  timeBase = header1.timeBase;
  deltaF  = 1.0/timeBase;

  sft->deltaF  = deltaF;
  sft->f0      = fminBinIndex*deltaF;
  sft->epoch.gpsSeconds     = header1.gpsSeconds;
  sft->epoch.gpsNanoSeconds = header1.gpsNanoSeconds;
  
  offset = fminBinIndex - header1.fminBinIndex;
  ASSERT (offset >= 0, status, SFTFILEIOH_EVAL, SFTFILEIOH_MSGEVAL);  

  length=sft->data->length;
    
  if (length > 0){
    COMPLEX8  *dataIn1 = NULL;
    COMPLEX8  *dataIn;
    COMPLEX8  *dataOut;
    INT4      n;
    REAL4     factor;
    
    ASSERT (sft->data->data, status, SFTFILEIOH_ENULL,  SFTFILEIOH_MSGENULL);
    ASSERT (header1.length >= offset+length, status, SFTFILEIOH_EVAL, SFTFILEIOH_MSGEVAL);
    dataIn1 = (COMPLEX8 *)LALMalloc(length*sizeof(COMPLEX8) );
    /* skip offset data points and read the required amount of data */
    ASSERT (0 == fseek(fp, offset * sizeof(COMPLEX8), SEEK_CUR), status, SFTFILEIOH_EVAL, SFTFILEIOH_MSGEVAL); 
    errorcode = fread( (void *)dataIn1, length*sizeof(COMPLEX8), 1, fp);

    dataIn = dataIn1;
    dataOut = sft->data->data;
    n= length;
    
    /* is this correct  ? */
    /* or if the data is normalized properly, there will be no need to consider. */
    factor = ((REAL4)length)/ ((REAL4) header1.length);
    
    /* let's normalize */
    while (n-- >0){
      dataOut->re = dataIn->re*factor;
      dataOut->im = dataIn->im*factor;
      ++dataOut;
      ++dataIn;
    }    
    LALFree(dataIn1);    
  }
  
  fclose(fp);
  
  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);
}

/* ----------------------------------------------------------------------*/
/* *******************************  <lalVerbatim file="SFTfileIOD"> */
void
LALWriteSFTtoFile (LALStatus  *status,
		   const SFTtype    *sft,
		   const CHAR       *outfname)
{/*   *********************************************  </lalVerbatim> */
  /* function to write an entire SFT to a file which has been read previously*/
  
  FILE  *fp = NULL;
  COMPLEX8  *inData;
  SFTHeader  header;
  INT4  i;
  REAL4  *rawdata;

  /* --------------------------------------------- */
  INITSTATUS (status, "WriteSFTtoFile", SFTFILEIOC);
  ATTATCHSTATUSPTR (status);   
 
  /*   Make sure the arguments are not NULL and perform basic checks*/ 
  ASSERT (sft,   status, SFTFILEIOH_ENULL, SFTFILEIOH_MSGENULL);
  ASSERT (sft->data,  status, SFTFILEIOH_EVAL, SFTFILEIOH_MSGEVAL);
  ASSERT (sft->data->length > 0,  status, SFTFILEIOH_EVAL, SFTFILEIOH_MSGEVAL);
  ASSERT (sft->data->data,   status, SFTFILEIOH_ENULL, SFTFILEIOH_MSGENULL);
  ASSERT (sft->deltaF, status, SFTFILEIOH_EVAL, SFTFILEIOH_MSGEVAL);
  ASSERT (outfname, status, SFTFILEIOH_ENULL, SFTFILEIOH_MSGENULL); 

  /* fill in the header information */
  header.version = 1.0;
  header.gpsSeconds = sft->epoch.gpsSeconds;
  header.gpsNanoSeconds = sft->epoch.gpsNanoSeconds;
  header.timeBase = 1.0 / sft->deltaF;
  header.fminBinIndex = (INT4)floor(sft->f0 / sft->deltaF + 0.5);
  header.length = sft->data->length; 

  /* open the file for writing */
  fp = fopen(outfname, "wb");
  if (fp == NULL) {
    ABORT (status, SFTFILEIOH_EFILE,  SFTFILEIOH_MSGEFILE);
  }

  /* write the header*/
  if( fwrite((void* )&header, sizeof(header), 1, fp) != 1) {
    ABORT (status, SFTFILEIOH_EWRITE, SFTFILEIOH_MSGEWRITE);
  }

  /* write data into a contiguous REAL4-array */
  rawdata = LALCalloc (1, 2* header.length * sizeof(REAL4));
  if (rawdata == NULL) {
    fclose (fp);
    ABORT (status, SFTFILEIOH_EMEM, SFTFILEIOH_MSGEMEM);    
  }

  inData = sft->data->data;
  for ( i = 0; i < header.length; i++)
    {
      rawdata[2 * i] = inData[i].re;
      rawdata[2 * i + 1] = inData[i].im;
    }
  
  if (fwrite( rawdata, 2 * header.length * sizeof(REAL4), 1, fp) != 1) {
    ABORT (status, SFTFILEIOH_EWRITE, SFTFILEIOH_MSGEWRITE);
  }
  LALFree (rawdata);
  fclose(fp);

  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);

} /* WriteSFTtoFile() */

/*----------------------------------------------------------------------
 * Reinhards draft versions of SFT IO functions, 
 * which combines Alicia's ReadSFTbinHeader1() and ReadSFTtype()
 *----------------------------------------------------------------------*/

/* <lalVerbatim file="SFTfileIOD"> */
void
LALReadSFTfile (LALStatus *status, 
		SFTtype **sft, 		/* output SFT */
		REAL8 fmin, 		/* lower frequency-limit */
		REAL8 fmax,		/* upper frequency-limit */
		const CHAR *fname)	/* filename */
{ /* </lalVerbatim> */
  FILE     *fp = NULL;
  SFTHeader  header;		/* SFT file-header version1 */
  REAL8 deltaF, f0, Band;
  UINT4 offset0, offset1, SFTlen;
  SFTtype *outputSFT;
  UINT4 i;
  REAL4 renorm;
  BOOLEAN swapEndian = 0;
  REAL8 version;
  REAL4 *rawdata = NULL;

  INITSTATUS (status, "ReadSFTfile", SFTFILEIOC);
  ATTATCHSTATUSPTR (status); 
  
  ASSERT (sft, status, SFTFILEIOH_ENULL,  SFTFILEIOH_MSGENULL);
  ASSERT (*sft == NULL, status, SFTFILEIOH_ENONULL, SFTFILEIOH_MSGENONULL);
  ASSERT (fname,  status, SFTFILEIOH_ENULL,  SFTFILEIOH_MSGENULL);
  
  /* opening the SFT binary file */
  if ( (fp = fopen( fname, "rb")) == NULL) {
    ABORT (status, SFTFILEIOH_EFILE, SFTFILEIOH_MSGEFILE);
  }
  
  /* read the header */
  if (fread( &header, sizeof(header), 1, fp) != 1) {
    ABORT (status, SFTFILEIOH_EHEADER,  SFTFILEIOH_MSGEHEADER);
  }
  version = header.version;
  /* check endian-ness and SFT-version */
  if (version > 1000000) {
    endian_swap ((CHAR*)&version, sizeof(version), 1);
    swapEndian = 1;
  }

  /* this is currently only supporting SFT version 1.0 ! */
  if (version > 1.0) {
    ABORT (status, SFTFILEIOH_EVERSION, SFTFILEIOH_MSGEVERSION);
  }
  
  /* do byte-swapping for header if needed */
  if (swapEndian) {
    TRY ( ByteSwapSFTHeader (status->statusPtr, &header), status);
  }
   
  deltaF = 1.0 / header.timeBase;
  f0 = header.fminBinIndex * deltaF;
  Band = header.length * deltaF;

  /* check that the required frequency-interval is part of the SFT */
  if ( (fmin < f0) || (fmax > f0 + Band) ) {
    ABORT (status, SFTFILEIOH_EFREQBAND, SFTFILEIOH_MSGEFREQBAND);
  }

  /* find the right bin offsets to read data from */
  /* the rounding here is chosen such that the required 
   * frequency-interval is _guaranteed_ to lie within the 
   * returned range 
   */
  offset0 = (UINT4)((fmin - f0) / deltaF);	/* round this down */
  offset1 = (UINT4) ceil((fmax - f0) / deltaF);	/* round this up! */
  SFTlen = offset1 - offset0;

  /* skip offset data points and read the required amount of data */
  if (fseek(fp, offset0 * sizeof(COMPLEX8), SEEK_CUR) != 0) {
    ABORT (status, SFTFILEIOH_ESEEK, SFTFILEIOH_MSGESEEK);
  }
  
  /* now allocate the SFT to be returned */
  outputSFT = LALCalloc (1, sizeof(SFTtype) );
  if (outputSFT == NULL) {
    fclose(fp);
    ABORT (status, SFTFILEIOH_EMEM, SFTFILEIOH_MSGEMEM);
  }
  LALCCreateVector (status->statusPtr, &(outputSFT->data), SFTlen);
  BEGINFAIL (status) {
    fclose(fp);
    LALFree (outputSFT);
  } ENDFAIL (status);

  rawdata = LALCalloc (1, 2*SFTlen*sizeof(REAL4));
  if (rawdata == NULL)
    {
      LALCDestroyVector (status->statusPtr, &(outputSFT->data));
      LALFree (outputSFT);
      ABORT (status, SFTFILEIOH_EMEM, SFTFILEIOH_MSGEMEM);
    }

  /* careful here: the data-vector contains a COMPLEX8-struct, so
     we better don't rely on memory-alignment and read into a
     REAL4 array first */
  if (fread( rawdata, 2*SFTlen * sizeof(REAL4), 1, fp) != 1) {
    fclose(fp);
    LALCDestroyVector (status->statusPtr, &(outputSFT->data));
    LALFree (outputSFT);
    LALFree (rawdata);
    ABORT (status, SFTFILEIOH_EREAD, SFTFILEIOH_MSGEREAD);
  }
  fclose(fp);

  /* do endian byte-swapping on the data if required */
  if ( swapEndian )
    endian_swap ( (CHAR*)rawdata, sizeof(REAL4), 2 * SFTlen);

  /* is this correct  ? */
  /* FIXME: double-check this !! */
  renorm = 1.0 * SFTlen / header.length;
    
  /* let's re-normalize and fill data into output-vector */
  for (i=0; i < SFTlen; i++)
    {
      outputSFT->data->data[i].re = renorm * rawdata[2 * i];
      outputSFT->data->data[i].im = renorm * rawdata[2 * i + 1];
    }    
  
  LALFree (rawdata);

  /* now fill in the header-info */
  outputSFT->deltaF  = deltaF;
  outputSFT->f0      = f0 + offset0 * deltaF;
  outputSFT->epoch.gpsSeconds     = header.gpsSeconds;
  outputSFT->epoch.gpsNanoSeconds = header.gpsNanoSeconds;
  
  /* that's it: return */
  *sft = outputSFT;

  DETATCHSTATUSPTR (status);
  RETURN(status);

} /* LALReadSFTfile() */

/***********************************************************************
 * internal helper functions
 ***********************************************************************/
/* a little endian-swapper needed for SFT reading/writing */
static void 
endian_swap(CHAR * pdata, size_t dsize, size_t nelements)
{
  UINT4 i, j, indx;
  CHAR tempbyte;

  if (dsize <= 1) return;

  for (i=0; i<nelements; i++)
    {
      indx = dsize;
      for (j=0; j<dsize/2; j++)
	{
	  tempbyte = pdata[j];
	  indx = indx - 1;
	  pdata[j] = pdata[indx];
	  pdata[indx] = tempbyte;
	}
      
      pdata = pdata + dsize;
    }
  
  return;

} /* endian swap */

/*----------------------------------------------------------------------
 * do correct byte-swapping for a whole SFT-header
 * 
 * currently only SFT-v1.0 is supported
 *----------------------------------------------------------------------*/
void
ByteSwapSFTHeader (LALStatus *stat, SFTHeader *header)
{
  INITSTATUS (stat, "LALByteSwapSFTHeader", SFTFILEIOC);

  /* we have to do this one by one, because memory-alignment
   * is not guaranteed !! 
   */
  endian_swap ((CHAR*)&(header->version), sizeof(header->version), 1);

  if (header->version > 1) {
    ABORT (stat, SFTFILEIOH_EVERSION, SFTFILEIOH_MSGEVERSION);
  }

  endian_swap ((CHAR*)&(header->gpsSeconds), sizeof(header->gpsSeconds), 1);
  endian_swap ((CHAR*)&(header->gpsNanoSeconds), sizeof(header->gpsNanoSeconds), 1);
  endian_swap ((CHAR*)&(header->timeBase), sizeof(header->timeBase), 1);
  endian_swap ((CHAR*)&(header->fminBinIndex), sizeof(header->fminBinIndex), 1);
  endian_swap ((CHAR*)&(header->length), sizeof(header->length), 1);

  RETURN (stat);

} /* ByteSwapSFTHeader() */
