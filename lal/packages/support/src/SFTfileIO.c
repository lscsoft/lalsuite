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
\idx{LALReadSFTdata()}
\idx{LALReadSFTfile()}
\idx{LALReadSFTfiles()}
\idx{LALWriteSFTtoFile()}
\vspace{0.1in}
\input{SFTfileIOD}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Description}

Will be written soon! :)

Watch this space.

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
#include <sys/types.h>
#include <dirent.h>

#include <lal/SFTfileIO.h>

NRCSID (SFTFILEIOC, "$Id$");



/* for private use: a vector of 'strings' */
typedef struct {
  UINT4 length;
  CHAR **data;
} StringVector;

/* prototypes for internal helper functions */
static StringVector *find_files (const CHAR *globdir);
static void DestroyStringVector (StringVector *strings);
static void LALCopySFT (LALStatus *stat, SFTtype *dest, const SFTtype *src);
static void endian_swap(CHAR * pdata, size_t dsize, size_t nelements);

void write_timeSeriesR4 (FILE *fp, const REAL4TimeSeries *series);
void write_timeSeriesR8 (FILE *fp, const REAL8TimeSeries *series);
void dump_SFT (const SFTtype *sft, FILE *fp);
void LALwriteSFTtoXMGR (LALStatus *stat, const SFTtype *sft, const CHAR *fname);

/* number of bytes in SFT-header v1.0 */
static const size_t header_len_v1 = 32;	

/*
 * The functions that make up the guts of this module
 */


/*----------------------------------------------------------------------*/
/* <lalVerbatim file="SFTfileIOD"> */
void LALReadSFTheader (LALStatus  *status,
                       SFTHeader   *header,
		       const CHAR  *fname)
{ /* </lalVerbatim> */
  
  FILE *fp = NULL;
  SFTHeader  header1;
  CHAR *rawheader = NULL;
  CHAR *ptr = NULL;
  REAL8 version;
  BOOLEAN swapEndian = 0;

  INITSTATUS (status, "LALReadSFTHeader", SFTFILEIOC);
  ATTATCHSTATUSPTR (status); 
  
  /*   Make sure the arguments are not NULL: */ 
  ASSERT (header, status, SFTFILEIOH_ENULL,  SFTFILEIOH_MSGENULL);
  ASSERT (fname,  status, SFTFILEIOH_ENULL,  SFTFILEIOH_MSGENULL);
  
  /* opening the SFT binary file */
  fp = fopen( fname, "r");
  if (fp == NULL) {
    ABORT (status, SFTFILEIOH_EFILE,  SFTFILEIOH_MSGEFILE);
  }
  
  /* read version-number */
  if  (fread (&version, sizeof(version), 1, fp) != 1) {
    ABORT (status, SFTFILEIOH_EHEADER,  SFTFILEIOH_MSGEHEADER);
  }

  /* check endian-ness */
  if (version > 1000000) {
    endian_swap ((CHAR*)&version, sizeof(version), 1);
    swapEndian = 1;
  }

  /* check compatibility of version with this function */
  if (version != 1) {
    ABORT (status, SFTFILEIOH_EVERSION, SFTFILEIOH_MSGEVERSION);
  }

  /* big-endian not yet supported */
  if (swapEndian) {
    ABORT (status, SFTFILEIOH_EENDIAN, SFTFILEIOH_MSGEENDIAN);
  }

  /* read the header */
  rawheader = LALCalloc (1, header_len_v1);
  if (rawheader == NULL) {
    ABORT (status, SFTFILEIOH_EMEM, SFTFILEIOH_MSGEMEM);    
  }
  
  rewind (fp);	/* go back to start */
  if (fread( rawheader, header_len_v1, 1, fp) != 1) {
    LALFree (rawheader);
    ABORT (status, SFTFILEIOH_EHEADER,  SFTFILEIOH_MSGEHEADER);
  }
  fclose(fp);

  /* now fill-in the header-struct with the appropriate fields */
  /* NOTE: we have to do it this way, because the struct can have 
   * padding in memory, so the fields are not guaranteed to lie 'close'
   */
  ptr = rawheader;
  header1.version	= *(REAL8*) ptr;
  ptr += sizeof(REAL8);
  header1.gpsSeconds 	= *(INT4*) ptr;
  ptr += sizeof(INT4);
  header1.gpsNanoSeconds= *(INT4*) ptr;
  ptr += sizeof(INT4);
  header1.timeBase      = *(REAL8*) ptr;
  ptr += sizeof(REAL8);
  header1.fminBinIndex  = *(INT4*) ptr;
  ptr += sizeof(INT4);
  header1.length     	= *(INT4*) ptr;

  LALFree (rawheader);

  /* ----- do some consistency-checks on the header-fields: ----- */

  /* 1. version has to be an integer! */
  if ( header1.version - (UINT4)(header1.version) != 0) {
    ABORT (status, SFTFILEIOH_EHEADER,  SFTFILEIOH_MSGEHEADER);
  }

  /* 2. gps_sec and gps_nsec >= 0 */
  if ( (header1.gpsSeconds < 0) || (header1.gpsNanoSeconds <0) ) {
    ABORT (status, SFTFILEIOH_EHEADER,  SFTFILEIOH_MSGEHEADER);
  }

  /* 3. tbase > 0 */
  if ( header1.timeBase <= 0 ) {
    ABORT (status, SFTFILEIOH_EHEADER,  SFTFILEIOH_MSGEHEADER);
  }

  /* 4. fminindex >= 0 */
  if (header1.fminBinIndex < 0) {
    ABORT (status, SFTFILEIOH_EHEADER,  SFTFILEIOH_MSGEHEADER);
  }
  
  /* 5. nsamples >= 0 */
  if (header1.length < 0) {
    ABORT (status, SFTFILEIOH_EHEADER,  SFTFILEIOH_MSGEHEADER);
  }

  /* ok, the SFT-header seems consistent, so let's return it */
  *header = header1;
  
  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* LALReadSFTheader() */



/*----------------------------------------------------------------------
 * Read a given frequency-range from an SFT-file into a memory struct SFTtype
 *
 * NOTE: currently only SFT-spec v1.0 is supported! 
 *
 * NOTE2: this is a convenience wrapper for LALReadSFTheader() and LALReadSFTdata()
 *
 *----------------------------------------------------------------------*/
/* <lalVerbatim file="SFTfileIOD"> */
void
LALReadSFTfile (LALStatus *stat, 
		SFTtype **sft, 		/* output SFT */
		REAL8 fmin, 		/* lower frequency-limit */
		REAL8 fmax,		/* upper frequency-limit */
		const CHAR *fname)	/* filename */
{ /* </lalVerbatim> */
  SFTHeader  header;		/* SFT file-header version1 */
  REAL8 deltaF;
  UINT4 readlen;
  INT4 fminBinIndex, fmaxBinIndex;
  SFTtype *outputSFT = NULL;
  UINT4 i;
  REAL4 renorm;

  INITSTATUS (stat, "LALReadSFTfile", SFTFILEIOC);
  ATTATCHSTATUSPTR (stat); 
  
  ASSERT (sft, stat, SFTFILEIOH_ENULL,  SFTFILEIOH_MSGENULL);
  ASSERT (*sft == NULL, stat, SFTFILEIOH_ENONULL, SFTFILEIOH_MSGENONULL);
  ASSERT (fname,  stat, SFTFILEIOH_ENULL,  SFTFILEIOH_MSGENULL);
  ASSERT (fmin < fmax, stat, SFTFILEIOH_EVAL, SFTFILEIOH_MSGEVAL);

  /* read the header */
  TRY ( LALReadSFTheader (stat->statusPtr, &header, fname), stat);

  /* ----- figure out which data we want to read ----- */
  deltaF = 1.0 / header.timeBase;

  /* find the right frequency-bin and number of bins
   * The rounding here is chosen such that the required 
   * frequency-interval is _guaranteed_ to lie within the 
   * returned range  */
  fminBinIndex = (INT4) floor (fmin * header.timeBase);  /* round this down */
  fmaxBinIndex = (INT4) ceil  (fmax * header.timeBase);  /* round up */

  readlen = (UINT4)(fmaxBinIndex - fminBinIndex) + 1;	/* number of bins to read */

  /* allocate the final SFT to be returned */
  TRY ( LALCreateSFTtype (stat->statusPtr, &outputSFT, readlen), stat);


  /* and read it, using the lower-level function: */
  LALReadSFTdata (stat->statusPtr, outputSFT, fname, fminBinIndex);
  BEGINFAIL (stat) {
    LALDestroySFTtype (stat->statusPtr, &outputSFT);
  } ENDFAIL (stat);


  /***************************************************
   * FIXME: questionable renormalization follows... to be clarified
   ***************************************************/
  renorm = 1.0 * readlen / header.length;
  /* **************************************************/

  /* let's re-normalize and fill data into output-vector */
  if (renorm != 1)
    for (i=0; i < readlen; i++)
      {
	outputSFT->data->data[i].re *= renorm;
	outputSFT->data->data[i].im *= renorm;
      }    
  
  /* that's it: return */
  *sft = outputSFT;

  DETATCHSTATUSPTR (stat);
  RETURN(stat);

} /* LALReadSFTfile() */



/* ----------------------------------------------------------------------
 * read a whole bunch of SFTs at once: given a vector of filenames, return
 * a vector of SFTs
 *
 * the input-glob is a bit special: the file-part is interpreted as a 
 * match-string: so test1/thisdir/SFT  found match all "*SFt* files
 * in test1/thisdir/
 *
 *----------------------------------------------------------------------*/
/* <lalVerbatim file="SFTfileIOD"> */
void
LALReadSFTfiles (LALStatus *stat,
		 SFTVector **sftvect,
		 REAL8 fmin,
		 REAL8 fmax,
		 const CHAR *globdir)
{ /* </lalVerbatim> */

  UINT4 i, numSFTs;
  SFTVector *out = NULL;
  SFTtype *sft = NULL;
  StringVector *fnames;

  INITSTATUS (stat, "LALReadSFTfiles", SFTFILEIOC);
  ATTATCHSTATUSPTR (stat); 
  
  ASSERT (sftvect, stat, SFTFILEIOH_ENULL,  SFTFILEIOH_MSGENULL);
  ASSERT (*sftvect == NULL, stat, SFTFILEIOH_ENONULL, SFTFILEIOH_MSGENONULL);
  ASSERT (globdir,  stat, SFTFILEIOH_ENULL,  SFTFILEIOH_MSGENULL);
  ASSERT (fmin < fmax, stat, SFTFILEIOH_EVAL, SFTFILEIOH_MSGEVAL);


  /* make filelist 
   * NOTE: we don't use glob() as it was reported to fail under condor */
  if ( (fnames = find_files (globdir)) == NULL) {
    ABORT (stat, SFTFILEIOH_EGLOB, SFTFILEIOH_MSGEGLOB);
  }

  numSFTs = fnames->length;

  for (i=0; i < numSFTs; i++)
    {
      LALReadSFTfile (stat->statusPtr, &sft, fmin, fmax, fnames->data[i]);
      BEGINFAIL (stat) {
	if (out) LALDestroySFTVector (stat->statusPtr, &out);
      } ENDFAIL (stat);
      /* this is a bit tricky: first we need to read one SFT to know how
       * many frequency-bins we need */
      if (out == NULL) {
	LALCreateSFTVector (stat->statusPtr, &out, numSFTs, sft->data->length);
      }
	
      /* transfer the returned SFT into the SFTVector 
       * this is a bit complicated by the fact that it's a vector
       * of SFTtypes, not pointers: therefore we need to *COPY* the stuff !
       */
      LALCopySFT (stat->statusPtr, &(out->data[i]), sft);
      BEGINFAIL (stat) {
	LALDestroySFTVector (stat->statusPtr, &out);
	LALDestroySFTtype (stat->statusPtr, &sft);
      } ENDFAIL (stat);

      LALDestroySFTtype (stat->statusPtr, &sft);
      sft = NULL;

    } /* for i < numSFTs */

  DestroyStringVector (fnames);

  *sftvect = out;

  DETATCHSTATUSPTR (stat);
  RETURN (stat);

} /* LALReadSFTfiles () */



/* ----------------------------------------------------------------------*/
/* function to write an entire SFT to a file 
 * 
 * NOTE: currently only SFT-spec v1.0 is supported
 *
 *----------------------------------------------------------------------*/
/* <lalVerbatim file="SFTfileIOD"> */
void
LALWriteSFTtoFile (LALStatus  *status,
		   const SFTtype *sft,
		   const CHAR *outfname)
{ /* </lalVerbatim> */

  FILE  *fp = NULL;
  COMPLEX8  *inData;
  INT4  i;
  UINT4 datalen;
  REAL4  *rawdata;
  CHAR *rawheader, *ptr;
  SFTHeader header;

  INITSTATUS (status, "LALWriteSFTtoFile", SFTFILEIOC);
  ATTATCHSTATUSPTR (status);   
 
  /*   Make sure the arguments are not NULL and perform basic checks*/ 
  ASSERT (sft,   status, SFTFILEIOH_ENULL, SFTFILEIOH_MSGENULL);
  ASSERT (sft->data,  status, SFTFILEIOH_EVAL, SFTFILEIOH_MSGEVAL);
  ASSERT (sft->deltaF > 0, status, SFTFILEIOH_EVAL, SFTFILEIOH_MSGEVAL);
  ASSERT (outfname, status, SFTFILEIOH_ENULL, SFTFILEIOH_MSGENULL); 

  /* fill in the header information */
  header.version = 1.0;
  header.gpsSeconds = sft->epoch.gpsSeconds;
  header.gpsNanoSeconds = sft->epoch.gpsNanoSeconds;
  header.timeBase = 1.0 / sft->deltaF;
  header.fminBinIndex = (INT4) floor (sft->f0 / sft->deltaF + 0.5);	/* round to closest int! */
  header.length = sft->data->length; 

  /* build raw header for writing to disk */
  rawheader = LALCalloc (1, header_len_v1);
  if (rawheader == NULL) {
    ABORT (status, SFTFILEIOH_EMEM, SFTFILEIOH_MSGEMEM);    
  }
  ptr = rawheader;
  *(REAL8*) ptr = header.version;
  ptr += sizeof (REAL8);
  *(INT4*) ptr = header.gpsSeconds;
  ptr += sizeof (INT4);
  *(INT4*) ptr = header.gpsNanoSeconds;
  ptr += sizeof (INT4);
  *(REAL8*) ptr = header.timeBase;
  ptr += sizeof (REAL8);
  *(INT4*) ptr = header.fminBinIndex;
  ptr += sizeof (INT4);
  *(INT4*) ptr = header.length; 

  /* write data into a contiguous REAL4-array */
  datalen = 2 * header.length * sizeof(REAL4);	/* amount of bytes for SFT-data */

  rawdata = LALCalloc (1, datalen);
  if (rawdata == NULL) {
    LALFree (rawheader);
    ABORT (status, SFTFILEIOH_EMEM, SFTFILEIOH_MSGEMEM);    
  }

  inData = sft->data->data;
  for ( i = 0; i < header.length; i++)
    {
      rawdata[2 * i]     = inData[i].re;
      rawdata[2 * i + 1] = inData[i].im;
    } /* for i < length */


  /* open the file for writing */
  fp = fopen(outfname, "wb");
  if (fp == NULL) {
    LALFree (rawheader);
    LALFree (rawdata);
    ABORT (status, SFTFILEIOH_EFILE,  SFTFILEIOH_MSGEFILE);
  }

  /* write the header*/
  if( fwrite( rawheader, header_len_v1, 1, fp) != 1) {
    LALFree (rawheader);
    LALFree (rawdata);
    fclose (fp);
    ABORT (status, SFTFILEIOH_EFILE, SFTFILEIOH_MSGEFILE);
  }
  
  /* write the data */
  if (fwrite( rawdata, datalen, 1, fp) != 1) {
    LALFree (rawheader);
    LALFree (rawdata);
    fclose (fp);
    ABORT (status, SFTFILEIOH_EFILE, SFTFILEIOH_MSGEFILE);
  }

  /* done */
  fclose(fp);
  LALFree (rawheader);
  LALFree (rawdata);


  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* WriteSFTtoFile() */

/* ------------------------------------------------------------
 * this is a function for low-level SFT data-reading:
 * the SFT-data is read starting from fminBinIndex and filled
 * into the pre-allocate vector sft of length N
 *
 * NOTE: !! NO re-normalization is done here!! this remains up
 *       to the caller of this function!!
 *
 *------------------------------------------------------------*/
/* <lalVerbatim file="SFTfileIOD"> */
void
LALReadSFTdata(LALStatus *status,
	       SFTtype    *sft,    /* asumed  memory is allocated  */
	       const CHAR *fname,
	       INT4 fminBinIndex)
{ /* </lalVerbatim> */
  FILE        *fp = NULL;
  SFTHeader  header;
  UINT4 offset, readlen;
  REAL4 *rawdata = NULL;
  REAL8 version;
  BOOLEAN swapEndian = 0;
  UINT4 i;

  INITSTATUS (status, "LALReadSFTtype", SFTFILEIOC);
  ATTATCHSTATUSPTR (status); 
  
  /*   Make sure the arguments are not NULL: */ 
  ASSERT (sft,   status, SFTFILEIOH_ENULL, SFTFILEIOH_MSGENULL);
  ASSERT (sft->data, status, SFTFILEIOH_ENULL, SFTFILEIOH_MSGENULL);
  ASSERT (fname, status, SFTFILEIOH_ENULL, SFTFILEIOH_MSGENULL);
  
  /* Read header */
  TRY ( LALReadSFTheader (status->statusPtr, &header, fname), status);


  /* check that the required frequency-interval is part of the SFT */
  readlen = sft->data->length;
  if ( (fminBinIndex < header.fminBinIndex) 
       || (fminBinIndex + (INT4)readlen > header.fminBinIndex + header.length) ) {
    ABORT (status, SFTFILEIOH_EFREQBAND, SFTFILEIOH_MSGEFREQBAND);
  }

  /* how many frequency-bins to skip */
  offset = fminBinIndex - header.fminBinIndex;

  /* open file for reading */
  if ( (fp = fopen( fname, "rb")) == NULL) {
    ABORT (status, SFTFILEIOH_EFILE, SFTFILEIOH_MSGEFILE);
  }
  /* read version-number and check endianness */
  if (fread (&version, sizeof(version), 1, fp) != 1) {
    fclose (fp);
    ABORT (status, SFTFILEIOH_EFILE, SFTFILEIOH_MSGEFILE);
  }
  if (version > 1000000) {
    endian_swap ((CHAR*)&version, sizeof(version), 1);
    swapEndian = 1;
  }
  /* only SFT-spec v1.0 currently supported! */
  if (version != 1) {	
    fclose (fp);
    ABORT (status, SFTFILEIOH_EVERSION, SFTFILEIOH_MSGEVERSION);
  }
  /* big-endian currently not supported */
  if (swapEndian) {
    fclose (fp);
    ABORT (status, SFTFILEIOH_EENDIAN, SFTFILEIOH_MSGEENDIAN);
  }

  /* skip SFT-header in file */
  rewind (fp);
  if (fseek(fp, header_len_v1, SEEK_SET) != 0) {
    ABORT (status, SFTFILEIOH_EFILE, SFTFILEIOH_MSGEFILE);
  }

  /* skip offset data points to the correct frequency-bin */
  if (fseek(fp, offset * 2 * sizeof(REAL4), SEEK_CUR) != 0) {
    ABORT (status, SFTFILEIOH_EFILE, SFTFILEIOH_MSGEFILE);
  }

  /* ----- prepare memory for data-reading ----- */
  rawdata = LALCalloc (1, 2 * readlen *sizeof(REAL4) );
  if (rawdata == NULL) {
    fclose (fp);
    ABORT (status, SFTFILEIOH_EMEM, SFTFILEIOH_MSGEMEM);
  }

  /* we don't rely on memory-packing, so we read into a REAL4 array first */
  if (fread( rawdata, 2 * readlen * sizeof(REAL4), 1, fp) != 1) {
    fclose(fp);
    LALFree (rawdata);
    ABORT (status, SFTFILEIOH_EFILE, SFTFILEIOH_MSGEFILE);
  }
  fclose(fp);

  /* now fill data into output-vector */
  for (i=0; i < readlen; i++)
    {
      sft->data->data[i].re = rawdata[2 * i];
      sft->data->data[i].im = rawdata[2 * i + 1];
    }    
  
  LALFree (rawdata);

  /* now fill in the header-info */
  strncpy (sft->name, fname, LALNameLength);
  sft->name[LALNameLength - 1 ] = '\0';	/* make sure it's 0-terminated */
  sft->deltaF  			= 1.0 / header.timeBase;
  sft->f0      			= fminBinIndex / header.timeBase;
  sft->epoch.gpsSeconds     	= header.gpsSeconds;
  sft->epoch.gpsNanoSeconds 	= header.gpsNanoSeconds;

  
  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* LALReadSFTdata() */


/*----------------------------------------------------------------------
 *  copy an entire SFT-type into another
 *  we require the destination to have at least as many frequency-bins
 *  as the source, but it can have less..
 *----------------------------------------------------------------------*/
void
LALCopySFT (LALStatus *stat, SFTtype *dest, const SFTtype *src)
{

  INITSTATUS( stat, "LALDestroySFTVector", SFTFILEIOC);

  ASSERT (dest,  stat, SFTFILEIOH_ENULL,  SFTFILEIOH_MSGENULL);
  ASSERT (dest->data,  stat, SFTFILEIOH_ENULL,  SFTFILEIOH_MSGENULL);
  ASSERT (src, stat, SFTFILEIOH_ENULL,  SFTFILEIOH_MSGENULL);
  ASSERT (src->data, stat, SFTFILEIOH_ENULL,  SFTFILEIOH_MSGENULL);

  /* some hard requirements */
  if ( dest->data->length < src->data->length ) {
    ABORT (stat, SFTFILEIOH_ECOPYSIZE, SFTFILEIOH_MSGECOPYSIZE);
  }
  
  /* copy head */
  strcpy (dest->name, src->name);
  dest->epoch = src->epoch;
  dest->f0 = src->f0;
  dest->deltaF = src->deltaF;
  dest->sampleUnits = src->sampleUnits;
  /* copy data */
  memcpy (dest->data->data, src->data->data, dest->data->length * sizeof (dest->data->data[0]));
  
  RETURN (stat);

} /* LALCopySFT() */



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
 * glob() has been reported to fail under condor, so we use our own
 * function to get a filelist from a directory, using a sub-pattern.
 *
 * This is not really a glob, but a "trick": the globdir is interpreted as
 *  "directory/pattern", ie. we actually match "directory/ *pattern*"
 *
 *----------------------------------------------------------------------*/
StringVector *
find_files (const CHAR *globdir)
{
  CHAR *dname, *ptr1, *ptr2;
  CHAR fpattern[512];
  size_t dirlen;
  DIR *dir;
  struct dirent *entry;
  CHAR **filelist = NULL; 
  UINT4 numFiles = 0;
  StringVector *ret = NULL;
  UINT4 j;
  UINT4 namelen;

  ptr1 = strrchr (globdir, '/');
  
  if (ptr1) 
    dirlen = (size_t)(ptr1 - globdir) + 1;
  else
    dirlen = 2;   /* for "." */

  if ( (dname = LALCalloc (1, dirlen)) == NULL) {
    return (NULL);
  }
  if (ptr1) {
    strncpy (dname, globdir, dirlen);
    dname[dirlen-1] = '\0';
  }
  else
    strcpy (dname, ".");
  
  /* copy the rest as a substring for matching */
  ptr2 = ptr1 + 1;
  strcpy (fpattern, ptr2);

  /* now go through the filelist in this directory */
  if ( (dir = opendir(dname)) == NULL)
    {
      LALPrintError ("Can't open data-directory `%s`\n", dname);
      LALFree (dname);
      return (NULL);
    }

  while ( (entry = readdir (dir)) != NULL )
    {
      if ( strstr (entry->d_name, fpattern) ) 	/* found a matching file */
	{
	  numFiles ++;
	  if ( (filelist = LALRealloc (filelist, numFiles * sizeof(CHAR*))) == NULL) {
	    return (NULL);
	  }

	  namelen = strlen(entry->d_name) + strlen(dname) + 2 ;
	  if ( (filelist[ numFiles - 1 ] = LALCalloc (1, namelen)) == NULL) {
	    for (j=0; j < numFiles; j++)
	      LALFree (filelist[j]);
	    LALFree (filelist);
	    return (NULL);
	  }

	  sprintf(filelist[numFiles-1], "%s/%s", dname, entry->d_name);
	}

    } /* while more directory entries */

  closedir (dir);
  LALFree (dname);

  /* ok, did we find anything? */
  if (numFiles == 0)
    return (NULL);

  if ( (ret = LALCalloc (1, sizeof (StringVector) )) == NULL) 
    {
      for (j=0; j<numFiles; j++)
	LALFree (filelist[j]);
      LALFree (filelist);
      return (NULL);
    }

  ret->length = numFiles;
  ret->data = filelist;

  return (ret);

} /* find_files() */


void
DestroyStringVector (StringVector *strings)
{
  UINT4 i;

  for (i=0; i < strings->length; i++)
    LALFree (strings->data[i]);
  
  LALFree (strings->data);
  LALFree (strings);

} /* DestroyStringVector () */


/***********************************************************************
 *
 * the following are
 * DEBUG + testing functions: not meant for general use!! 
 *
 ***********************************************************************/

/* write a time-series into a text-file */
void
write_timeSeriesR4 (FILE *fp, const REAL4TimeSeries *series)
{
  REAL8 timestamp; 
  UINT4 i;

  if (series == NULL)
    {
      printf ("\nempty input!\n");
      return; 
    }

  timestamp = 1.0*series->epoch.gpsSeconds + series->epoch.gpsNanoSeconds * 1.0e-9;

  for( i = 0; i < series->data->length; i++)
    {
      fprintf( fp, "%16.9f %e\n", timestamp, series->data->data[i] );
      timestamp += series->deltaT;
    }

  return;

} /* write_timeSeriesR4() */

void
write_timeSeriesR8 (FILE *fp, const REAL8TimeSeries *series)
{
  REAL8 timestamp; 
  UINT4 i;

  timestamp = 1.0*series->epoch.gpsSeconds + series->epoch.gpsNanoSeconds * 1.0e-9;

  for( i = 0; i < series->data->length; i++)
    {
      fprintf( fp, "%f %e\n", timestamp, series->data->data[i] );
      timestamp += series->deltaT;
    }

  return;

} /* write_timeSeriesR4() */


/* write an SFT in xmgrace-readable format */
void 
LALwriteSFTtoXMGR (LALStatus *stat, const SFTtype *sft, const CHAR *fname)
{
  FILE *fp;

  REAL4 val;
  UINT4 i;
  REAL8 Tsft, freqBand;
  REAL8 f0, df, ff;
  UINT4 set, nsamples;

  const CHAR *xmgrHeader = 
    "@version 50103\n"
    "@xaxis label \"f (Hz)\"\n";

  INITSTATUS( stat, "LALwriteSFTtoXMGR", SFTFILEIOC);
  ATTATCHSTATUSPTR( stat );

  fp = fopen (fname, "w");

  if( !fp ) {
    ABORT (stat, SFTFILEIOH_EFILE, SFTFILEIOH_MSGEFILE);
  }


  f0 = sft->f0;
  df = sft->deltaF;
  nsamples = sft->data->length;
  Tsft = 1.0 / sft->deltaF;
  freqBand = nsamples * df;

  fprintf (fp, xmgrHeader);
  fprintf (fp, "@subtitle \"epoch = (%d s, %d ns), Tsft = %f\"\n", 
	   sft->epoch.gpsSeconds, sft->epoch.gpsNanoSeconds, Tsft);

  set = 0;
  /* Print set header. */
  fprintf( fp, "\n@target G0.S%d\n@type xy\n", set);
  fprintf( fp, "@s%d color (0,0,1)\n", set );
  for (i=0; i < nsamples; i++)
    {
      ff = f0 + i*df;
      val = sft->data->data[i].re;

      fprintf(fp, "%f %e\n", ff, val);
        
    } /* for i < nsamples */

  set ++;
  /* Print set header. */
  fprintf( fp, "\n@target G0.S%d\n@type xy\n", set);
  fprintf( fp, "@s%d color (0,1,0)\n", set );
  for (i=0; i < nsamples; i++)
    {
      ff = f0 + i*df;
      val = sft->data->data[i].im;

      fprintf(fp, "%f %e\n", ff, val);
        
    } /* for i < nsamples */

  set ++;
  /* Print set header. */
  fprintf( fp, "\n@target G0.S%d\n@type xy\n", set);
  fprintf( fp, "@s%d color (1,0,0)\n", set );
  for (i=0; i < nsamples; i++)
    {
      ff = f0 + i*df;
      val = (sft->data->data[i].re * sft->data->data[i].re + sft->data->data[i].im*sft->data->data[i].im) / freqBand;
      
      fprintf(fp, "%f %e\n", ff, val);
        
    } /* for i < nsamples */

  
  fclose(fp);

  DETATCHSTATUSPTR( stat );
  RETURN (stat);
  
} /* write_SFT() */


/* dump an SFT into a text-file */
void dump_SFT (const SFTtype *sft, FILE *fp)
{
  REAL4 valre, valim;
  UINT4 i;
  REAL8 Tsft, freqBand;
  REAL8 f0, df, ff;
  UINT4 nsamples;

  f0 = sft->f0;
  df = sft->deltaF;
  nsamples = sft->data->length;
  Tsft = 1.0 / sft->deltaF;
  freqBand = nsamples * df;

  for (i=0; i < nsamples; i++)
    {
      ff = f0 + i*df;
      valre = sft->data->data[i].re;
      valim = sft->data->data[i].im;
      fprintf(fp, "%f %e %e\n", ff, valre, valim);
        
    } /* for i < nsamples */
  
  return;
  
} /* dump_SFT() */
