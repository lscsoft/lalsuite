/*-----------------------------------------------------------------------
 *
 * File Name: SFTbin.c
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

/************************************ <lalVerbatim file="SFTbinCV">
Author: Sintes, A.M., Krishnan, B.
$Id$
************************************* </lalVerbatim> */

/* <lalLaTeX>  *******************************************************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\section{MOdule \texttt{SFTbin.c}}
\label{ss:SFTbin.c}
Routines for reading SFT binary files

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Prototypes}
\vspace{0.1in}
\input{SFTbinD}
\index{\verb&LALReadSFTbinHeader1()&}
\index{\verb&LALReadSFTtype()&}
\index{\verb&ReadCOMPLEX8SFTbinData1()&}
\index{\verb&ReadCOMPLEX16SFTbinData1()&}
\idx{ReadSFTfile()}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Description}

the output of the periodogram should be properly normalized !!!


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Uses}
\begin{verbatim}
LALHO()
\end{verbatim}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Notes}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\vfill{\footnotesize\input{SFTbinCV}}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

*********************************************** </lalLaTeX> */


#include "./SFTbinIO.h"

NRCSID (SFTBINC, "$Id$");

/*
 * The functions that make up the guts of this module
 */


/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */

/* *******************************  <lalVerbatim file="SFTbinD"> */
void LALReadSFTbinHeader1 (LALStatus  *status,
                     SFTHeader1   *header,
                     CHAR          *fname)
{ /*   *********************************************  </lalVerbatim> */
  
  FILE     *fp = NULL;
  SFTHeader1  header1;
  size_t    errorcode;
  
  /* --------------------------------------------- */
  INITSTATUS (status, "LALReadSFTbinHeader1", SFTBINC);
  ATTATCHSTATUSPTR (status); 
  
  /*   Make sure the arguments are not NULL: */ 
  ASSERT (header, status, SFTBINH_ENULL,  SFTBINH_MSGENULL);
  ASSERT (fname,  status, SFTBINH_ENULL,  SFTBINH_MSGENULL);
  
  /* opening the SFT binary file */
  fp = fopen( fname, "r");
  ASSERT (fp, status, SFTBINH_EFILE,  SFTBINH_MSGEFILE);
  
  /* read the header */
  errorcode = fread( (void *)&header1, sizeof(SFTHeader1), 1, fp);
  ASSERT (errorcode ==1, status, SFTBINH_EHEADER,  SFTBINH_MSGEHEADER);
  ASSERT (header1.endian ==1, status, SFTBINH_EENDIAN, SFTBINH_MSGEENDIAN);

  header->endian     = header1.endian;
  header->gpsSeconds = header1.gpsSeconds;
  header->gpsNanoSeconds = header1.gpsNanoSeconds;
  header->timeBase       = header1.timeBase;
  header->fminBinIndex   = header1.fminBinIndex ;
  header->length     = header1.length ;
  
  fclose(fp);
  
  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);
}

/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/* *******************************  <lalVerbatim file="SFTbinD"> */
void LALReadSFTtype(LALStatus  *status,
                 SFTtype    *sft,    /* asumed  memory is allocated  */
                 CHAR       *fname,
		 INT4       fminBinIndex)
{ /*   *********************************************  </lalVerbatim> */

  FILE        *fp = NULL;
  SFTHeader1  header1;
  size_t     errorcode;
  INT4       offset;
  INT4       length;
  REAL8      deltaF,timeBase;
  
  /* --------------------------------------------- */
  INITSTATUS (status, "LALReadSFTtype", SFTBINC);
  ATTATCHSTATUSPTR (status); 
  
  /*   Make sure the arguments are not NULL: */ 
  ASSERT (sft,   status, SFTBINH_ENULL, SFTBINH_MSGENULL);
  ASSERT (sft->data, status, SFTBINH_ENULL, SFTBINH_MSGENULL);
  ASSERT (fname, status, SFTBINH_ENULL, SFTBINH_MSGENULL);
  
  /* opening the SFT binary file */
  fp = fopen( fname, "r");
  ASSERT (fp, status, SFTBINH_EFILE,  SFTBINH_MSGEFILE);
  
  /* read the header */
  errorcode = fread( (void *)&header1, sizeof(SFTHeader1), 1, fp);
  ASSERT (errorcode == 1,      status, SFTBINH_EHEADER, SFTBINH_MSGEHEADER);
  ASSERT (header1.endian == 1, status, SFTBINH_EENDIAN, SFTBINH_MSGEENDIAN);
  
  /* check that the data we want is in the file and it is correct */
  ASSERT (header1.timeBase > 0.0, status, SFTBINH_EVAL, SFTBINH_MSGEVAL);

  timeBase = header1.timeBase;
  deltaF  = 1.0/timeBase;

  sft->deltaF  = deltaF;
  sft->f0      = fminBinIndex*deltaF;
  sft->epoch.gpsSeconds     = header1.gpsSeconds;
  sft->epoch.gpsNanoSeconds = header1.gpsNanoSeconds;
  
  offset = fminBinIndex - header1.fminBinIndex;
  ASSERT (offset >= 0, status, SFTBINH_EVAL, SFTBINH_MSGEVAL);  

  length=sft->data->length;
    
  if (length > 0){
    COMPLEX8  *dataIn1 = NULL;
    COMPLEX8  *dataIn;
    COMPLEX8  *dataOut;
    INT4      n;
    REAL8     factor;
    
    ASSERT (sft->data->data, status, SFTBINH_ENULL,  SFTBINH_MSGENULL);
    ASSERT (header1.length >= offset+length, status, SFTBINH_EVAL, SFTBINH_MSGEVAL);
    dataIn1 = (COMPLEX8 *)LALMalloc(length*sizeof(COMPLEX8) );
    /* skip offset data points and read the required amount of data */
    ASSERT (0 == fseek(fp, offset * sizeof(COMPLEX8), SEEK_CUR), status, SFTBINH_EVAL, SFTBINH_MSGEVAL); 
    errorcode = fread( (void *)dataIn1, length*sizeof(COMPLEX8), 1, fp);

    dataIn = dataIn1;
    dataOut = sft->data->data;
    n= length;
    
    /* is this correct  ? */
    /* or if the data is normalized properly, there will be no need to consider. */
    factor = ((double) length)/ ((double) header1.length);
    
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

/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/* *******************************  <lalVerbatim file="SFTbinD"> */
void ReadCOMPLEX8SFTbinData1(LALStatus  *status,
                   COMPLEX8SFTData1    *sft,
                   CHAR                *fname)
{ /*   *********************************************  </lalVerbatim> */

  FILE        *fp = NULL;
  SFTHeader1  header1;
  size_t     errorcode;
  INT4       offset;
  INT4       length;
  
  /* --------------------------------------------- */
  INITSTATUS (status, "ReadCOMPLEX8SFTbinData1", SFTBINC);
  ATTATCHSTATUSPTR (status); 
  
  /*   Make sure the arguments are not NULL: */ 
  ASSERT (sft,   status, SFTBINH_ENULL, SFTBINH_MSGENULL);
  ASSERT (fname, status, SFTBINH_ENULL, SFTBINH_MSGENULL);
  
  /* opening the SFT binary file */
  fp = fopen( fname, "r");
  ASSERT (fp, status, SFTBINH_EFILE,  SFTBINH_MSGEFILE);
  
  /* read the header */
  errorcode = fread( (void *)&header1, sizeof(SFTHeader1), 1, fp);
  ASSERT (errorcode == 1,      status, SFTBINH_EHEADER, SFTBINH_MSGEHEADER);
  ASSERT (header1.endian == 1, status, SFTBINH_EENDIAN, SFTBINH_MSGEENDIAN);
  
  /* check that the data we want is in the file and it is correct */
  ASSERT (header1.timeBase > 0.0, status, SFTBINH_EVAL, SFTBINH_MSGEVAL);
  sft->timeBase             = header1.timeBase;
  sft->epoch.gpsSeconds     = header1.gpsSeconds;
  sft->epoch.gpsNanoSeconds = header1.gpsNanoSeconds;
  
  offset = sft->fminBinIndex - header1.fminBinIndex;
  ASSERT (offset >= 0, status, SFTBINH_EVAL, SFTBINH_MSGEVAL);  

  length=sft->length;
  
  if (length > 0){
    COMPLEX8  *dataIn1 = NULL;
    COMPLEX8  *dataIn;
    COMPLEX8  *dataOut;
    INT4      n;
    REAL8     factor;
    
    ASSERT (sft->data, status, SFTBINH_ENULL,  SFTBINH_MSGENULL);
    ASSERT (header1.length >= offset+length, status, SFTBINH_EVAL, SFTBINH_MSGEVAL);
    dataIn1 = (COMPLEX8 *)LALMalloc(length*sizeof(COMPLEX8) );
    /* skip offset data points and read the required amount of data */
    ASSERT (0 == fseek(fp, offset * sizeof(COMPLEX8), SEEK_CUR), status, SFTBINH_EVAL, SFTBINH_MSGEVAL); 
    errorcode = fread( (void *)dataIn1, length*sizeof(COMPLEX8), 1, fp);

    dataIn = dataIn1;
    dataOut = sft->data;
    n= length;
    
    /* is this correct  ? */
    /* or if the data is normalized properly, there will be no need to consider. */
    factor = ((double) length)/ ((double) header1.length);
    
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

/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/* *******************************  <lalVerbatim file="SFTbinD"> */
void ReadCOMPLEX16SFTbinData1(LALStatus  *status,
                   COMPLEX16SFTData1    *sft,
                   CHAR                *fname)
{ /*   *********************************************  </lalVerbatim> */

  FILE        *fp = NULL;
  SFTHeader1  header1;
  size_t     errorcode;
  INT4       offset;
  INT4       length;
  
  /* --------------------------------------------- */
  INITSTATUS (status, "ReadCOMPLEX16SFTbinData1", SFTBINC);
  ATTATCHSTATUSPTR (status); 
  
  /*   Make sure the arguments are not NULL: */ 
  ASSERT (sft,   status, SFTBINH_ENULL, SFTBINH_MSGENULL);
  ASSERT (fname, status, SFTBINH_ENULL, SFTBINH_MSGENULL);
  
  /* opening the SFT binary file */
  fp = fopen( fname, "r");
  ASSERT (fp, status, SFTBINH_EFILE,  SFTBINH_MSGEFILE);
  
  /* read the header */
  errorcode = fread( (void *)&header1, sizeof(SFTHeader1), 1, fp);
  ASSERT (errorcode == 1,      status, SFTBINH_EHEADER, SFTBINH_MSGEHEADER);
  ASSERT (header1.endian == 1, status, SFTBINH_EENDIAN, SFTBINH_MSGEENDIAN);
  
  /* check that the data we want is in the file and it is correct */
  ASSERT (header1.timeBase > 0.0, status, SFTBINH_EVAL, SFTBINH_MSGEVAL);
  sft->timeBase             = header1.timeBase;
  sft->epoch.gpsSeconds     = header1.gpsSeconds;
  sft->epoch.gpsNanoSeconds = header1.gpsNanoSeconds;
  
  offset = sft->fminBinIndex - header1.fminBinIndex;
  ASSERT (offset >= 0, status, SFTBINH_EVAL, SFTBINH_MSGEVAL);
  
  length=sft->length;
  
  if (length > 0){
    COMPLEX16  *dataIn1 = NULL;
    COMPLEX16  *dataIn;
    COMPLEX16  *dataOut;
    INT4      n;
    REAL8     factor;
    
    ASSERT (sft->data, status, SFTBINH_ENULL,  SFTBINH_MSGENULL);
    ASSERT (header1.length >= offset+length, status, SFTBINH_EVAL, SFTBINH_MSGEVAL);
    dataIn1 = (COMPLEX16 *)LALMalloc(length*sizeof(COMPLEX16) );
    /* skip offset data points and read the required amount of data */
    ASSERT (0 == fseek(fp, offset * sizeof(COMPLEX16), SEEK_CUR), status, SFTBINH_EVAL, SFTBINH_MSGEVAL); 
    errorcode = fread( (void *)dataIn1, length*sizeof(COMPLEX16), 1, fp);
    
    dataOut = sft->data;
    dataIn = dataIn1;
    n= length;
    
    /* is this correct  ? */
    /* or if the data is normalized properly, there will be no need to consider. */
    factor = ((double) length)/ ((double) header1.length);
    
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
 

/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/* *******************************  <lalVerbatim file="SFTbinD"> */
void WriteCOMPLEX8SFTbinData1 (LALStatus          *status,
		       COMPLEX8SFTData1   *sft,
		       CHAR               *outfname)
{/*   *********************************************  </lalVerbatim> */
  /* function to write an entire SFT to a file which has been read previously*/
  
  FILE  *fp = NULL;
  COMPLEX8  *inData;
  SFTHeader1  header;
  INT4  i, w;
  REAL4  rePt, imPt;

  /* --------------------------------------------- */
  INITSTATUS (status, "WriteCOMPLEX8SFTbinData1", SFTBINC);
  ATTATCHSTATUSPTR (status);   
 
  /*   Make sure the arguments are not NULL and perform basic checks*/ 
  ASSERT (sft,   status, SFTBINH_ENULL, SFTBINH_MSGENULL);
  ASSERT (sft->length > 0,  status, SFTBINH_EVAL, SFTBINH_MSGEVAL);
  ASSERT (sft->fminBinIndex > 0,  status, SFTBINH_EVAL, SFTBINH_MSGEVAL);
  ASSERT (sft->timeBase > 0,  status, SFTBINH_EVAL, SFTBINH_MSGEVAL);
  ASSERT (sft->data,   status, SFTBINH_ENULL, SFTBINH_MSGENULL);
  ASSERT (outfname, status, SFTBINH_ENULL, SFTBINH_MSGENULL); 

  /* fill in the header information */
  header.endian = 1.0;
  header.gpsSeconds = sft->epoch.gpsSeconds;
  header.gpsNanoSeconds = sft->epoch.gpsNanoSeconds;
  header.timeBase = sft->timeBase;
  header.fminBinIndex = sft->fminBinIndex;
  header.length = sft->length; 


  /* open the file for writing */
  fp = fopen(outfname, "w");
  ASSERT (fp, status, SFTBINH_EFILE,  SFTBINH_MSGEFILE);

  /* write the header*/
  w = fwrite((void* )&header, sizeof(header), 1, fp);
  ASSERT (w == 1, status, SFTBINH_EWRITE, SFTBINH_MSGEWRITE);

  inData = sft->data;
  for ( i = 0; i < header.length; i++){

    rePt = inData[i].re;
    imPt = inData[i].im;

    /* write the real and imaginary parts */
    w = fwrite((void* )&rePt, sizeof(REAL4), 1, fp);
    ASSERT (w == 1, status, SFTBINH_EWRITE, SFTBINH_MSGEWRITE);
    w = fwrite((void* )&imPt, sizeof(REAL4), 1, fp);
    ASSERT (w == 1, status, SFTBINH_EWRITE, SFTBINH_MSGEWRITE);
  }

  fclose(fp);

  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);
}



/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/* *******************************  <lalVerbatim file="SFTbinD"> */
void WriteCOMPLEX16SFTbinData1(LALStatus          *status,
		       COMPLEX16SFTData1   *sft,
		       CHAR               *outfname)
{/*   *********************************************  </lalVerbatim> */
  /* function to write entirely an SFT which has been previously read */
  
  FILE  *fp = NULL;
  COMPLEX16  *inData;
  SFTHeader1 header;
  INT4  i, w;
  REAL4  rePt, imPt;

  /* --------------------------------------------- */
  INITSTATUS (status, "WriteCOMPLEX16SFTbinData1", SFTBINC);
  ATTATCHSTATUSPTR (status);   
 
  /*   Make sure the arguments are not NULL and perform basic checks*/ 
  ASSERT (sft,   status, SFTBINH_ENULL, SFTBINH_MSGENULL);
  ASSERT (sft->length > 0,  status, SFTBINH_EVAL, SFTBINH_MSGEVAL);
  ASSERT (sft->fminBinIndex > 0,  status, SFTBINH_EVAL, SFTBINH_MSGEVAL);
  ASSERT (sft->timeBase > 0,  status, SFTBINH_EVAL, SFTBINH_MSGEVAL);
  ASSERT (sft->data,   status, SFTBINH_ENULL, SFTBINH_MSGENULL);
  ASSERT (outfname, status, SFTBINH_ENULL, SFTBINH_MSGENULL); 

  /* fill in the header information */
  header.endian = 1.0;
  header.gpsSeconds = sft->epoch.gpsSeconds;
  header.gpsNanoSeconds = sft->epoch.gpsNanoSeconds;
  header.timeBase = sft->timeBase;
  header.fminBinIndex = sft->fminBinIndex;
  header.length = sft->length; 


  /* open the file for writing */
  fp = fopen(outfname, "w");
  ASSERT (fp, status, SFTBINH_EFILE,  SFTBINH_MSGEFILE);

  /* write the header*/
  w = fwrite((void* )&header, sizeof(header), 1, fp);
  ASSERT (w == 1, status, SFTBINH_EWRITE, SFTBINH_MSGEWRITE);

  inData = sft->data;
  for ( i = 0; i < header.length; i++){

    rePt = inData[i].re;
    imPt = inData[i].im;

    /* write the real and imaginary parts */
    w = fwrite((void* )&rePt, sizeof(REAL8), 1, fp);
    ASSERT (w == 1, status, SFTBINH_EWRITE, SFTBINH_MSGEWRITE);
    w = fwrite((void* )&imPt, sizeof(REAL8), 1, fp);
    ASSERT (w == 1, status, SFTBINH_EWRITE, SFTBINH_MSGEWRITE);
  }

  fclose(fp);

  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);
}


/*----------------------------------------------------------------------
 * Reinhard draft versions of SFT IO functions, 
 * which combines Alicia's ReadSFTbinHeader1() and ReadSFTtype()
 *----------------------------------------------------------------------*/

/* <lalVerbatim file="SFTbinD"> */
void
LALReadSFTfile (LALStatus *status, 
		SFTtype **sft, 		/* output SFT */
		REAL8 fmin, 		/* lower frequency-limit */
		REAL8 fmax,		/* upper frequency-limit */
		const CHAR *fname)	/* filename */
{ /* </lalVerbatim> */
  FILE     *fp = NULL;
  SFTHeader1  header;		/* SFT file-header version1 */
  REAL8 deltaF, f0, Band;
  UINT4 offset0, offset1, SFTlen;
  SFTtype *outputSFT;
  UINT4 i;
  REAL4 factor;
  

  INITSTATUS (status, "ReadSFTfile", SFTBINC);
  ATTATCHSTATUSPTR (status); 
  
  ASSERT (sft, status, SFTBINH_ENULL,  SFTBINH_MSGENULL);
  ASSERT (*sft == NULL, status, SFTBINH_ENONULL, SFTBINH_MSGENONULL);
  ASSERT (fname,  status, SFTBINH_ENULL,  SFTBINH_MSGENULL);
  
  /* opening the SFT binary file */
  if ( (fp = fopen( fname, "rb")) == NULL) {
    ABORT (status, SFTBINH_EFILE, SFTBINH_MSGEFILE);
  }
  
  /* read the header */
  
  if (fread( &header, sizeof(header), 1, fp) != 1) {
    ABORT (status, SFTBINH_EHEADER,  SFTBINH_MSGEHEADER);
  }
  if (header.endian != 1) {
    ABORT (status, SFTBINH_EENDIAN, SFTBINH_MSGEENDIAN);
  }

  deltaF = 1.0 / header.timeBase;
  f0 = header.fminBinIndex * deltaF;
  Band = header.length * deltaF;

  /* check that the required frequency-interval is part of the SFT */
  if ( (fmin < f0) || (fmax > f0 + Band) ) {
    ABORT (status, SFTBINH_EFREQBAND, SFTBINH_MSGEFREQBAND);
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
    ABORT (status, SFTBINH_ESEEK, SFTBINH_MSGESEEK);
  }
  
  /* now allocate the SFT to be returned */
  outputSFT = LALMalloc ( sizeof(SFTtype) );
  if (outputSFT == NULL) {
    fclose(fp);
    ABORT (status, SFTBINH_EMEM, SFTBINH_MSGEMEM);
  }
  LALCCreateVector (status->statusPtr, &(outputSFT->data), SFTlen);
  BEGINFAIL (status) {
    fclose(fp);
    LALFree (outputSFT);
  } ENDFAIL (status);

  if (fread( outputSFT->data->data, SFTlen * sizeof(COMPLEX8), 1, fp) != 1) {
    fclose(fp);
    LALCDestroyVector (status->statusPtr, &(outputSFT->data));
    LALFree (outputSFT);
    ABORT (status, SFTBINH_EREAD, SFTBINH_MSGEREAD);
  }
  fclose(fp);

  /* is this correct  ? */
  /* or if the data is normalized properly, there will be no need to consider. */
  factor = 1.0 * SFTlen / header.length;
    
  /* let's re-normalize */
  for (i=0; i < SFTlen; i++)
    {
      outputSFT->data->data[i].re *= factor;
      outputSFT->data->data[i].im *= factor;
    }    
  
  /* now fill in the header-info */
  outputSFT->deltaF  = deltaF;
  outputSFT->f0      = header.fminBinIndex * deltaF;
  outputSFT->epoch.gpsSeconds     = header.gpsSeconds;
  outputSFT->epoch.gpsNanoSeconds = header.gpsNanoSeconds;
  
  /* that's it: return */
  *sft = outputSFT;

  DETATCHSTATUSPTR (status);
  RETURN(status);

} /* LALReadSFTfile() */
