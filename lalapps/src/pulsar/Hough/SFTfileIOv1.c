/*
 * Copyright (C) 2010 Karl Wette
 * Copyright (C) 2004, 2005 R. Prix, B. Machenschalk, A.M. Sintes
 *
 * crc64() taken from SFTReferenceLibrary.c Copyright (C) 2004 Bruce Allen
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

/*---------- INCLUDES ----------*/
#include <sys/types.h>
#include <errno.h>
#include <string.h>
#include <strings.h>
#include <ctype.h>

#include <lal/LALStdio.h>
#include <lal/LALString.h>
#include <lal/FileIO.h>
#include <lal/SFTfileIO.h>
#include <lal/StringVector.h>
#include <lal/Sequence.h>
#include <lal/ConfigFile.h>
#include <lal/UserInputParse.h>
#include <lal/LogPrintf.h>
#include <lal/SFTutils.h>

#include "SFTfileIOv1.h"

/** \cond DONT_DOXYGEN */

/*---------- internal prototypes ----------*/
static void endian_swap(CHAR * pdata, size_t dsize, size_t nelements);

/*==================== FUNCTION DEFINITIONS ====================*/

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


/**
 * [DEPRECATED]: Low-level function to read only the SFT-header of a given file.
 *
 * NOTE: don't use! This function is obsolete and SFT-v1-specific and only kept for
 * backwards-compatibility with Hough-codes.
 */
void
LALReadSFTheader (LALStatus  *status,			/**< pointer to LALStatus structure */
		  SFTHeader   *header,	/**< [out] returned header */
		  const CHAR  *fname)	/**< path+filename */
{
  FILE *fp = NULL;
  SFTHeader  header1;
  CHAR *rawheader = NULL;
  CHAR *ptr = NULL;
  CHAR inVersion[8];
  REAL8 version;
  BOOLEAN swapEndian = 0;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /*   Make sure the arguments are not NULL: */
  ASSERT (header, status, SFTFILEIO_ENULL,  SFTFILEIO_MSGENULL);
  ASSERT (fname,  status, SFTFILEIO_ENULL,  SFTFILEIO_MSGENULL);

  /* opening the SFT binary file */
  fp = LALFopen( fname, "rb" );
  if (fp == NULL) {
    ABORT (status, SFTFILEIO_EFILE,  SFTFILEIO_MSGEFILE);
  }

  /* read version-number */
  if  (fread (inVersion, sizeof(inVersion), 1, fp) != 1) {
    fclose (fp);
    if (lalDebugLevel) XLALPrintError ("\nInvalid SFT-file: %s\n\n", fname);
    ABORT (status, SFTFILEIO_EHEADER,  SFTFILEIO_MSGEHEADER);
  }

  /* try version 1.0 */
  version = 1.0;
  /* try swapping the version if it is not equal */
  if(memcmp(inVersion,&version,sizeof(version))){
    endian_swap (inVersion, sizeof(inVersion),1);
    swapEndian = 1;
    /* fail if still not true */
    if(memcmp(inVersion,&version,sizeof(version))){
      fclose (fp);
      if (lalDebugLevel) XLALPrintError ("\nOnly v1-SFTs supported at the moment!: %s\n\n", fname);
      ABORT (status, SFTFILEIO_EHEADER,  SFTFILEIO_MSGEHEADER);
    }
  }

  /* read the whole header */
  rawheader = LALCalloc (1, sizeof(SFTHeader) );
  if (rawheader == NULL) {
    fclose (fp);
    ABORT (status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM);
  }

  rewind (fp);	/* go back to start */
  if (fread( rawheader, sizeof(SFTHeader), 1, fp) != 1) {
    fclose (fp);
    LALFree (rawheader);
    if (lalDebugLevel) XLALPrintError ("\nInvalid SFT-file: %s\n\n", fname);
    ABORT (status, SFTFILEIO_EHEADER,  SFTFILEIO_MSGEHEADER);
  }

  fclose(fp);

  /* now fill-in the header-struct with the appropriate fields */
  /* NOTE: we have to do it this way, because the struct can have
   * padding in memory, so the fields are not guaranteed to lie 'close'
   * Endian-swapping ist done here if necessary
   */
  ptr = rawheader;
  if (swapEndian) endian_swap((CHAR*)ptr,sizeof(REAL8),1);
  memcpy( &header1.version, ptr, sizeof(REAL8) );
  ptr += sizeof(REAL8);
  if (swapEndian) endian_swap((CHAR*)ptr,sizeof(INT4),1);
  memcpy( &header1.gpsSeconds, ptr, sizeof(INT4) );
  ptr += sizeof(INT4);
  if (swapEndian) endian_swap((CHAR*)ptr,sizeof(INT4),1);
  memcpy( &header1.gpsNanoSeconds, ptr, sizeof(INT4) );
  ptr += sizeof(INT4);
  if (swapEndian) endian_swap((CHAR*)ptr,sizeof(REAL8),1);
  memcpy( &header1.timeBase, ptr, sizeof(REAL8) );
  ptr += sizeof(REAL8);
  if (swapEndian) endian_swap((CHAR*)ptr,sizeof(INT4),1);
  memcpy( &header1.fminBinIndex, ptr, sizeof(INT4) );
  ptr += sizeof(INT4);
  if (swapEndian) endian_swap((CHAR*)ptr,sizeof(INT4),1);
  memcpy( &header1.length, ptr, sizeof(INT4) );

  LALFree (rawheader);

  /* ----- do some consistency-checks on the header-fields: ----- */

  /* gps_sec and gps_nsec >= 0 */
  if ( (header1.gpsSeconds < 0) || (header1.gpsNanoSeconds <0) ) {
    if (lalDebugLevel) XLALPrintError ("\nInvalid SFT-file: %s\n\n", fname);
    ABORT (status, SFTFILEIO_EHEADER,  SFTFILEIO_MSGEHEADER);
  }

  /* tbase > 0 */
  if ( header1.timeBase <= 0 ) {
    if (lalDebugLevel) XLALPrintError ("\nInvalid SFT-file: %s\n\n", fname);
    ABORT (status, SFTFILEIO_EHEADER,  SFTFILEIO_MSGEHEADER);
  }

  /* fminindex >= 0 */
  if (header1.fminBinIndex < 0) {
    if (lalDebugLevel) XLALPrintError ("\nInvalid SFT-file: %s\n\n", fname);
    ABORT (status, SFTFILEIO_EHEADER,  SFTFILEIO_MSGEHEADER);
  }

  /* nsamples >= 0 */
  if (header1.length < 0) {
    if (lalDebugLevel) XLALPrintError ("\nInvalid SFT-file: %s\n\n", fname);
    ABORT (status, SFTFILEIO_EHEADER,  SFTFILEIO_MSGEHEADER);
  }

  /* ok, the SFT-header seems consistent, so let's return it */
  *header = header1;

  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* LALReadSFTheader() */


/**
 * [DEPRECATED] This is a function for low-level SFT data-reading:
 * the SFT-data is read starting from fminBinIndex and filled
 * into the pre-allocate vector sft of length N
 *
 * NOTE: !! NO re-normalization is done here!! this remains up
 * to the caller of this function!!
 *
 */
void
LALReadSFTdata(LALStatus *status,			/**< pointer to LALStatus structure */
	       SFTtype    *sft,    /**< [out] output-SFT: assuming memory is allocated  */
	       const CHAR *fname,  /**< path+filename */
	       INT4 fminBinIndex)  /**< minimun frequency-index to read */
{
  FILE        *fp = NULL;
  SFTHeader  header;
  UINT4 offset, readlen;
  REAL4 *rawdata = NULL;
  CHAR inVersion[8];
  REAL8 version;
  BOOLEAN swapEndian = 0;
  UINT4 i;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /*   Make sure the arguments are not NULL: */
  ASSERT (sft,   status, SFTFILEIO_ENULL, SFTFILEIO_MSGENULL);
  ASSERT (sft->data, status, SFTFILEIO_ENULL, SFTFILEIO_MSGENULL);
  ASSERT (fname, status, SFTFILEIO_ENULL, SFTFILEIO_MSGENULL);

  /* Read header */
  TRY ( LALReadSFTheader (status->statusPtr, &header, fname), status);

  /* check that the required frequency-interval is part of the SFT */
  readlen = sft->data->length;
  if ( (fminBinIndex < header.fminBinIndex)
       || (fminBinIndex + (INT4)readlen > header.fminBinIndex + header.length) ) {
    ABORT (status, SFTFILEIO_EFREQBAND, SFTFILEIO_MSGEFREQBAND);
  }

  /* how many frequency-bins to skip */
  offset = fminBinIndex - header.fminBinIndex;

  /* open file for reading */
  if ( (fp = LALFopen( fname, "rb" )) == NULL) {
    ABORT (status, SFTFILEIO_EFILE, SFTFILEIO_MSGEFILE);
  }

  /* read version-number */
  if  (fread (inVersion, sizeof(inVersion), 1, fp) != 1) {
    fclose (fp);
    if (lalDebugLevel) XLALPrintError ("\nInvalid SFT-file: %s\n\n", fname);
    ABORT (status, SFTFILEIO_EHEADER,  SFTFILEIO_MSGEHEADER);
  }

  /* set invalid version */
  version = -1.0;

  if (version < 0) {
    /* try version 1.0 */
    version = 1.0;
    /* try swapping the version if it is not equal */
    if(memcmp(inVersion,&version,sizeof(version))){
      endian_swap (inVersion, sizeof(inVersion),1);
      swapEndian = 1;
      /* set invalid version if still not true */
      if(memcmp(inVersion,&version,sizeof(version))){
	version = -1;
	endian_swap (inVersion, sizeof(inVersion),1);
	swapEndian = 0;
      }
    }
  }

  /* fail if the version is invalid */
  if (version < 0) {
    fclose (fp);
    if (lalDebugLevel) XLALPrintError ("\nInvalid SFT-file: %s\n\n", fname);
    ABORT (status, SFTFILEIO_EHEADER,  SFTFILEIO_MSGEHEADER);
  }

  /* check compatibility of version with this function */
  if (version != 1) {
    fclose (fp);
    ABORT (status, SFTFILEIO_EVERSION, SFTFILEIO_MSGEVERSION);
  }

  /* skip SFT-header in file */
  rewind (fp);
  if (fseek(fp, sizeof(SFTHeader), SEEK_SET) != 0) {
    fclose (fp);
    ABORT (status, SFTFILEIO_EFILE, SFTFILEIO_MSGEFILE);
  }

  /* skip offset data points to the correct frequency-bin */
  if (fseek(fp, offset * 2 * sizeof(REAL4), SEEK_CUR) != 0) {
    fclose (fp);
    ABORT (status, SFTFILEIO_EFILE, SFTFILEIO_MSGEFILE);
  }

  /* ----- prepare memory for data-reading ----- */
  rawdata = LALCalloc (1, 2 * readlen *sizeof(REAL4) );
  if (rawdata == NULL) {
    fclose (fp);
    ABORT (status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM);
  }

  /* we don't rely on memory-packing, so we read into a REAL4 array first */
  if (fread( rawdata, 2 * readlen * sizeof(REAL4), 1, fp) != 1) {
    fclose(fp);
    LALFree (rawdata);
    ABORT (status, SFTFILEIO_EFILE, SFTFILEIO_MSGEFILE);
  }
  fclose(fp);

  /* now fill data into output-vector */
  for (i=0; i < readlen; i++)
    {
      if (swapEndian)
	endian_swap((CHAR*)&rawdata[2*i], sizeof(REAL4), 2);
      sft->data->data[i] = crectf( rawdata[2 * i], rawdata[2 * i + 1] );
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


/**
 * [DEPRECATED] Basic SFT reading-function.
 * Given a filename \a fname and frequency-limits [\a fMin, \a fMax],
 * returns an SFTtype \a sft containing the SFT-data.
 *
 * \note 1) the actual returned frequency-band is
 * <tt>[floor(Tsft * fMin), ceil(Tsft * fMax)] / Tsft</tt>,
 * i.e. we need to round to an integer frequency-bin within the SFT, but
 * the requested frequency-band is guaranteed to be contained in the output
 * (if present in the SFT-file), but can be slightly larger.
 *
 * 2) The special input <tt>fMin=fMax=0</tt> means to read and
 * return the <em>whole</em> frequency-band contained in the SFT-file.
 *
 * 3) Currently only SFTv1 are supported!!
 *
 */
void
LALReadSFTfile (LALStatus *status,			/**< pointer to LALStatus structure */
		SFTtype **sft, 		/**< [out] output SFT */
		REAL8 fMin, 		/**< lower frequency-limit */
		REAL8 fMax,		/**< upper frequency-limit */
		const CHAR *fname)	/**< path+filename */
{
  SFTHeader  header;		/* SFT file-header version1 */
  UINT4 readlen;
  INT4 fminBinIndex, fmaxBinIndex;
  SFTtype *outputSFT = NULL;
  UINT4 i;
  REAL4 renorm;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  ASSERT (sft, status, SFTFILEIO_ENULL,  SFTFILEIO_MSGENULL);
  ASSERT (*sft == NULL, status, SFTFILEIO_ENONULL, SFTFILEIO_MSGENONULL);
  ASSERT (fname,  status, SFTFILEIO_ENULL,  SFTFILEIO_MSGENULL);
  ASSERT (fMin <= fMax, status, SFTFILEIO_EVAL, SFTFILEIO_MSGEVAL);

  /* read the header */
  TRY ( LALReadSFTheader (status->statusPtr, &header, fname), status);

  /* ----- figure out which data we want to read ----- */

  /* special case: fMin==fMax==0 means "read all" */
  if ( (fMin == 0) && (fMax == 0) )
    {
      fminBinIndex = header.fminBinIndex;
      fmaxBinIndex = fminBinIndex + header.length - 1;
    }
  else
    {
      /* find the right frequency-bin and number of bins
       * The rounding here is chosen such that the required
       * frequency-interval is _guaranteed_ to lie within the
       * returned range  */
      fminBinIndex = (INT4) floor (fMin * header.timeBase);  /* round this down */
      fmaxBinIndex = (INT4) ceil  (fMax * header.timeBase);  /* round up */
    }

  readlen = (UINT4)(fmaxBinIndex - fminBinIndex) + 1;	/* number of bins to read */

  /* allocate the final SFT to be returned */
  XLAL_CHECK_LAL ( status, ( outputSFT = XLALCreateSFT ( readlen) ) != NULL, XLAL_EFUNC);


  /* and read it, using the lower-level function: */
  LALReadSFTdata (status->statusPtr, outputSFT, fname, fminBinIndex);
  BEGINFAIL (status) {
    XLALDestroySFT ( outputSFT);
  } ENDFAIL (status);


  /*
   * NOTE: the following renormalization is necessary for v1-SFTs
   * as the data are not multiplied by dt, therefore extracting a
   * sub-band B' of the total band B requires mulitplication of
   * the data by B'/B.
   * SFTv2 will store FFT-data multiplied by dt and then this will
   * not be necessary any more.
   *
   */
  renorm = 1.0 * readlen / header.length;

  /* let's re-normalize and fill data into output-vector */
  if (renorm != 1)
    for (i=0; i < readlen; i++)
      {
	outputSFT->data->data[i] *= ((REAL4) renorm);
      }

  /* that's it: return */
  *sft = outputSFT;

  DETATCHSTATUSPTR (status);
  RETURN(status);

} /* LALReadSFTfile() */


/**
 * [DEPRECATED] Higher-level SFT-reading function to read a whole vector of SFT files
 * and return an SFTvector \a sftvect. The handling of
 * [fMin, fMax] is identical to LALReadSFTfile().
 *
 * \note 1) the file-pattern \a fpattern can use a wide range of glob-patterns.
 *
 * 2) currently the SFTs matching the pattern are required to have the same
 * number of frequency bins, otherwise an error will be returned.
 * (This might be relaxed in the future).
 *
 * 3) This function does <em>not</em> use <tt>glob()</tt> and should therefore
 * be safe even under condor.
 *
 */
void
LALReadSFTfiles (LALStatus *status,			/**< pointer to LALStatus structure */
		 SFTVector **sftvect,	/**< [out] output SFT vector */
		 REAL8 fMin,	       	/**< lower frequency-limit */
		 REAL8 fMax,		/**< upper frequency-limit */
		 UINT4 wingBins,	/**< number of frequency-bins to be added left and right. */
		 const CHAR *fpattern)	/**< path/filepattern */
{

  UINT4 i, numSFTs;
  SFTVector *out = NULL;
  SFTtype *oneSFT = NULL;
  SFTHeader header;
  REAL8 dFreq; 		/* frequency spacing in SFT */
  REAL8 fWing;		/* frequency band to be added as "wings" to the 'physical band' */
  LALStringVector *fnames;
  UINT4 firstlen = 0;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  ASSERT (sftvect, status, SFTFILEIO_ENULL,  SFTFILEIO_MSGENULL);
  ASSERT (*sftvect == NULL, status, SFTFILEIO_ENONULL, SFTFILEIO_MSGENONULL);
  ASSERT (fpattern,  status, SFTFILEIO_ENULL,  SFTFILEIO_MSGENULL);
  ASSERT (fMin <= fMax, status, SFTFILEIO_EVAL, SFTFILEIO_MSGEVAL);

  /* make filelist
   * NOTE: we don't use glob() as it was reported to fail under condor */
  if ( (fnames = XLALFindFiles (fpattern)) == NULL) {
    ABORT (status, SFTFILEIO_EGLOB, SFTFILEIO_MSGEGLOB);
  }

  /* allocate head of output SFT-vector */
  if ( (out = LALCalloc ( 1, sizeof(SFTVector) )) == NULL ) {
    ABORT ( status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM );
  }

  /* read header of first sft to determine Tsft, and therefore dfreq */
  LALReadSFTheader(status->statusPtr, &header, fnames->data[0]);
  BEGINFAIL(status) {
    XLALDestroyStringVector (fnames);
  } ENDFAIL(status);

  numSFTs = fnames->length;
  dFreq = 1.0 / header.timeBase;
  fWing = wingBins * dFreq;

  /* main loop: load all SFTs and put them into the SFTvector */
  for (i=0; i < numSFTs; i++)
    {
      LALReadSFTfile (status->statusPtr, &oneSFT, fMin-fWing, fMax+fWing, fnames->data[i]);
      BEGINFAIL (status) {
	XLALDestroyStringVector (fnames);
	XLALDestroySFT ( oneSFT);
	if (out) XLALDestroySFTVector ( out);
      } ENDFAIL (status);

      if ( !firstlen )
	firstlen = oneSFT->data->length;
      /* make sure all SFTs have same length */
      if ( oneSFT->data->length != firstlen )
	{
	  XLALDestroyStringVector (fnames);
	  XLALDestroySFT ( oneSFT);
	  XLALDestroySFTVector ( out);
	  ABORT (status, SFTFILEIO_EDIFFLENGTH, SFTFILEIO_MSGEDIFFLENGTH);
	} /* if length(thisSFT) != common length */

      if ( XLALAppendSFT2Vector ( out, oneSFT ) != XLAL_SUCCESS ) {
	XLALDestroyStringVector (fnames);
	XLALDestroySFT ( oneSFT);
	XLALDestroySFTVector ( out);
	ABORTXLAL( status );
      }

      XLALDestroySFT ( oneSFT);
      oneSFT = NULL;	/* important for next call of LALReadSFTfile()! */

    } /* for i < numSFTs */

  XLALDestroyStringVector (fnames);

  *sftvect = out;

  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* LALReadSFTfiles () */


/** \endcond */
