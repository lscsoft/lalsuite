/* $Id$ */
/* This is a reference library for the SFT data format
 *
 * Copyright (C) 2004 Bruce Allen <ballen@uwm.edu>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2, or (at your option)
 * any later version.
 *
 * You should have received a copy of the GNU General Public License
 * (for example COPYING); if not, write to the Free
 * Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "SFTReferenceLibrary.h"
#include "config.h"

/*
  The quantity below is: D800000000000000 (base-16) =
  1101100000000000000000000000000000000000000000000000000000000000
  (base-2).  The primitive polynomial is x^64 + x^4 + x^3 + x + 1.
*/
#define POLY64 0xd800000000000000ULL
#define TABLELEN 256
#define BLOCKSIZE 65536

/* some local prototypes */
static unsigned long long crc64(const unsigned char* data, unsigned int length, unsigned long long crc);
static void swap2(char *location);
static void swap4(char *location);
static void swap8(char *location);
static int validate_sizes(void);

const char* ReferenceSFTLibraryVersion(void) {
  return PACKAGE_VERSION;
}

/* The crc64 checksum of M bytes of data at address data is returned
   by crc64(data, M, ~(0ULL)). Call the function multiple times to
   compute the checksum of data made in contiguous chunks, setting
   final argument to the previously accumulated checksum value. */
unsigned long long crc64(const unsigned char* data,
			 unsigned int length,
			 unsigned long long crc) {

  unsigned long long CRCTable[TABLELEN];
  unsigned int i;

  /* is there is no data, simply return previous checksum value */
  if (!length)
    return crc;
  
  /* initialize the CRC table for fast computation.  We could keep
     this table in memory to make the computation faster, but that is
     not re-entrant for multi-threaded code.
  */
  for (i = 0; i < TABLELEN; i++) {
    int j;
    unsigned long long part = i;
    for (j = 0; j < 8; j++) {
      if (part & 1)
        part = (part >> 1) ^ POLY64;
      else
        part >>= 1;
    }
    CRCTable[i] = part;
  }
  
  /* compute the CRC-64 code */
  for (i=0; i<length; i++) {
    unsigned long long temp1 = crc >> 8;
    unsigned long long temp2 = CRCTable[(crc ^ (unsigned long long) data[i]) & 0xff];
    crc = temp1 ^ temp2;
  }
  
  return crc;
}

/* swap 2, 4 or 8 bytes.  Point to low address */
static void swap2(char *location){
  char tmp=*location;
  *location=*(location+1);
  *(location+1)=tmp;
  return;
}
 
static void swap4(char *location){
  char tmp=*location;
  *location=*(location+3);
  *(location+3)=tmp;
  swap2(location+1);
  return;
}
 
static void swap8(char *location){
  char tmp=*location;
  *location=*(location+7);
  *(location+7)=tmp;
  tmp=*(location+1);
  *(location+1)=*(location+6);
  *(location+6)=tmp;
  swap4(location+2);
  return;
}

/* this routine checks that the assumed sizes and structure packing
   conventions of this code are valid.  These checks could also be
   done at compile time with #error statements */
static int validate_sizes(void) {
  if (
      sizeof(char) != 1      ||
      sizeof(int) != 4       ||
      sizeof(long long) != 8 ||
      sizeof(struct headertag1) != 32 ||
      sizeof(struct headertag2) != 48
      )
    return SFTESIZEWRONG;
  
  return SFTNOERROR;
}


/* translate return values from SFT routines into human-readable
   character string error messages */
const char *SFTErrorMessage(int errorcode) {

  switch (errorcode) {
  case 0:
    return "Success";
  case SFTENULLFP:
    return "SFT file pointer is NULL";
  case SFTESEEK:
    return "SFT fseek() failed in stream";
  case SFTEREAD:
    return "SFT fread() failed in stream";
  case SFTEUNKNOWN:
    return "SFT version in header is unknown (not 1 or 2)";
  case SFTEGPSNSEC:
    return "SFT header GPS nsec not in range 0 to 10^9-1";
  case SFTEBADCOMMENT:
    return "SFT comment length not a multiple of 8";
  case SFTEGETSTREAMPOS:
    return "SFT fgetpos() failed: unable to save position of stream";
  case SFTEBADCRC64:
    return "SFT corrupted: CRC checksum in header does not match data";
  case SFTERESTORESTREAMPOS:
    return "SFT fsetpos() failed: unable to restore stream position";
  case SFTENOMEM:
    return "SFT calloc() failed: unable to allocate memory";
  case SFTESIZEWRONG:
    return "SFT sizeof() objects does not match assumptions";
  case SFTEWRITE:
    return "SFT fwrite() failed";
  case SFTENULLPOINTER:
    return "SFT library routine passed a null data pointer";
  case SFTENONE:
    return "SFT file empty (zero length)";
  case SFTEHIDDENCOMMENT:
    return "SFT comment contains data AFTER null character";
  case SFTENONULLINCOMMENT:
    return "SFT comment field is not null-terminated";
  case SFTEGPSNOTINCREASING:
    return "SFT GPS times not increasing between SFT blocks";
  case SFTETBASECHANGES:
    return "SFT time base changes between SFT blocks";
  case SFTEFIRSTINDEXCHANGES:
    return "SFT first frequency index changes between SFT blocks";
  case SFTENSAMPLESCHANGES:
    return "SFT number of data samples changes between SFT blocks";
  case SFTEINSTRUMENTCHANGES:
    return "SFT instrument changes between SFT blocks";
  case SFTEVERSIONCHANGES:
    return "SFT version changes between SFT blocks";
  case SFTETBASENOTPOS:
    return "SFT time base is not positive";
  case SFTEFIRSTINDEXNEG:
    return "SFT first frequency index is negative";
  case SFTENSAMPLESNOTPOS:
    return "SFT number of data samples is not positive";
  case SFTEINSTRUMENTUNKNOWN:
    return "SFT detector not one of A1 B1 E1 G1 H1 H2 K1 L1 N1 O1 P1 T1 V1 V2";
  case SFTEBEFOREDATA:
    return "SFT data requested lies before available data";
  case SFTEAFTERDATA:
    return "SFT data requested lies after available data";
  }
  return "SFT Error Code not recognized";
}

/* return values: see header file
   On return, this routine leaves the file pointer at the end of the SFT
*/
int WriteSFT(FILE *fp,            /* stream to write to */
	     int gps_sec,         /* GPS sec of first sample */
	     int gps_nsec,        /* GPS nsec of first sample */
	     double tbase,        /* time baseline of SFTs */
	     int firstfreqindex,  /* index of first frequency bin included in data (0=DC)*/
	     int nsamples,        /* number of frequency bins to include in SFT */
	     const char *detector,/* channel-prefix defining detector */
	     const char *comment, /* null-terminated comment string to include in SFT */
	     float *data          /* points to nsamples x 2 x floats (Real/Imag)  */
	     ) {
  struct headertag2 header;
  int comment_length, inc, i;
  char pad[7];
  
  /* check that all data types have correct length */
  if (validate_sizes())
    return SFTESIZEWRONG;
  
  /* check that file pointer is valid */
  if (!fp)
    return SFTENULLFP;
  
  /* check that nsec times are sensible */
  if (gps_nsec < 0 || gps_nsec > 999999999)
    return SFTEGPSNSEC;

  /* check that tbase is positive */
  if (tbase<=0.0)
    return SFTETBASENOTPOS;

  /* check that first frequency index is non-negative */
  if (firstfreqindex<0)
    return SFTEFIRSTINDEXNEG;

  /* check that number of samples is 1 or greater */
  if (nsamples<=0)
    return SFTENSAMPLESNOTPOS;

  /* check that detector type is defined and data is present */
  if (!detector || !data)
    return SFTENULLPOINTER;

  /* check that detector type is recognized */
  if (strlen(detector)>2 || unknownDetector(detector))
    return SFTEINSTRUMENTUNKNOWN;

  /* comment length including null terminator to string must be an
     integer multiple of eight bytes. comment==NULL means 'no
     comment'  */
  if (comment) {
    comment_length = strlen(comment) + 1;
    inc = (8 - (comment_length % 8)) % 8;
    for (i=0; i<inc; i++)
      pad[i]=0;
  }
  else {
    comment_length=0;
    inc=0;
  }
  
  /* fill out header */
  header.version        = 2;
  header.gps_sec        = gps_sec;
  header.gps_nsec       = gps_nsec;
  header.tbase          = tbase;
  header.firstfreqindex = firstfreqindex;
  header.nsamples       = nsamples;
  header.crc64          = 0;
  header.detector[0]    = detector[0];
  header.detector[1]    = detector[1];
  header.padding[0]     = 0;
  header.padding[1]     = 0;
  header.comment_length = comment_length+inc;
  
  /* compute CRC of header */
  header.crc64 = crc64((const unsigned char *)&header, sizeof(header), ~(0ULL));
  
  /* compute CRC of comment */
  header.crc64 = crc64((const unsigned char *)comment, comment_length, header.crc64);
  
  /* compute CRC of comment padding */
  header.crc64 = crc64((const unsigned char *)pad, inc, header.crc64);

  /* compute CRC of data */
  header.crc64 = crc64((const unsigned char *)data, nsamples*2*sizeof(float), header.crc64);
  
  /* write the header to file */
  if (1 != fwrite((const void *)&header, sizeof(header), 1, fp))
    return SFTEWRITE;
  
  /* write the comment to file */
  if (comment_length != (int)fwrite((const void *)comment, 1, comment_length, fp))
    return SFTEWRITE;

  /* write comment padding to file */
  if (inc != (int)fwrite((const void *)pad, 1, inc, fp))
    return SFTEWRITE;

  /* write the data to the file.  Data must be packed
     REAL,IMAG,REAL,IMAG,... */
  if (nsamples != (int)fwrite((const void *)data, 2*sizeof(float), nsamples, fp))
    return SFTEWRITE;
  
  return SFTNOERROR;
}

/* On return this routine leaves the stream in the same
   place as when the routine was called */
int ReadSFTHeader(FILE *fp,                  /* stream to read */
		  struct headertag2 *info,   /* address to return header */
		  char **comment,            /* if non-NULL, put pointer to comment */
		  int *swapendian,           /* set nonzero if data in reverse endian order */
		  int validate)              /* validate checksum of the file */
{
  struct headertag2 header,header_unswapped;
  int swap=0;
  int what, version, retval=0;
  fpos_t streamposition;
  char *mycomment=NULL;

  /* check that all data types have correct length */
  if (validate_sizes())
    return SFTESIZEWRONG;

  /* check that pointers are valid */
  if (!fp)
    return SFTENULLFP;

  if (!info || !swapendian)
    return SFTENULLPOINTER;

  /* save stream position */
  if (fgetpos(fp, &streamposition))
    return SFTEGETSTREAMPOS;
  
  /* read in header.  Note that for v1 SFTs this reads 16 bytes past
     the header -- not important since old SFTs are much longer.  Note
     also that the second and third arguments to fread() are reversed
     from conventional practice, so that we can see if the file is
     zero length */
  if (sizeof(header) != (what=fread((void *)&header, 1, sizeof(header), fp))) {
    if (!what && feof(fp))
      retval=SFTENONE;
    else
      retval=SFTEREAD;
    goto error;
  }
  
  /* save an unswapped copy */
  header_unswapped=header;

  /* check endian ordering, and swap if needed */
  if (header.version != 1 && header.version != 2) {
    swap8((char *)&header.version);
    swap4((char *)&header.gps_sec);
    swap4((char *)&header.gps_nsec);
    swap8((char *)&header.tbase);
    swap4((char *)&header.firstfreqindex);
    swap4((char *)&header.nsamples);
    swap8((char *)&header.crc64);
    swap4((char *)&header.comment_length);
    swap = 1;
  }
  
  /* check if header version is recognized */
  if (header.version == 1) {
    version = 1;
    header.crc64          = 0;
    header.detector[0]    = 0;
    header.detector[1]    = 0;
    header.padding[0]     = 0;
    header.padding[1]     = 0;
    header.comment_length = 0;
  }
  else if (header.version == 2)
    version = 2;
  else {
    retval=SFTEUNKNOWN;
    goto error;
  }
  
  if (header.comment_length % 8) {
    retval=SFTEBADCOMMENT;
    goto error;
  }

  /* check that number of samples is 1 or greater.  We do this check
     before calculating CRC since the number of values to CRC check
     must be known first */
  if (header.nsamples<=0){
    retval=SFTENSAMPLESNOTPOS;
    goto error;
  }

  /* validate crc64 checksum ??  Do this BEFORE other checks since if
     a problem occurs it is more likely file corruption than a bad
     SFT */
  if (version !=1 && validate) {
    unsigned long long crc;
    unsigned long long crc64save=header.crc64;
    int total_length = header.comment_length + 2*sizeof(float)*header.nsamples;
    char block[BLOCKSIZE];
    fpos_t streamposition2;
    int i, tocheck, foundhidden=0, comment_left = header.comment_length, foundnull=!comment_left;

    /* save stream position */
    if (fgetpos(fp, &streamposition2)) {
      retval=SFTEGETSTREAMPOS;
      goto error;
    }

    /* compute CRC of header */
    header_unswapped.crc64 = 0ULL;;
    crc = crc64((unsigned char *)&header_unswapped, sizeof(header_unswapped), ~(0ULL));

    /* read data in lengths of BLOCKSIZE, computing CRC */
    while (total_length > 0) {
      /* read either BLOCKSIZE or amount remaining */
      int toread = (BLOCKSIZE < total_length) ? BLOCKSIZE : total_length;
      if (toread != (int)fread(block, 1, toread, fp)) {
	retval=SFTEREAD;
	goto error;
      }      
      total_length -= toread;
      crc = crc64((unsigned char *)block, toread, crc);

      /* check to see if comment contains NULL termination character
	 and no hidden characters */
      tocheck = (comment_left < toread) ? comment_left : toread;
      for (i=0; i<tocheck; i++) {
	if (!block[i])
	  foundnull=1;
	if (foundnull && block[i])
	  foundhidden=1;
      }
      comment_left -= tocheck;
    }
     
    /* check that checksum is consistent */
    if (crc != crc64save) {
#ifdef SFTDEVEL
      printf("CRC64 computes as %llu\n", crc);
      printf("CRC64 in SFT is   %llu\n", crc64save);
#endif
      retval=SFTEBADCRC64;
      goto error;
    }

    /* check that comment has correct NULL termination */
    if (!foundnull) {
      retval=SFTENONULLINCOMMENT;
      goto error;
    }

    /* check that comment has no hidden characters */
    if (foundhidden) {
      retval=SFTEHIDDENCOMMENT;
      goto error;
    }

    /* return to position just after header to read comment if desired */
    if (fsetpos(fp, &streamposition2)) {
      retval=SFTERESTORESTREAMPOS;
      goto error;
    }
  }

  /* check that time stamps are in a valid range */
  if (header.gps_nsec < 0 || header.gps_nsec > 999999999) {
    retval=SFTEGPSNSEC;
    goto error;
  }
  
  /* check that tbase is positive */
  if (header.tbase<=0.0) {
    retval=SFTETBASENOTPOS;
    goto error;
  }

  /* check that first frequency index is non-negative */
  if (header.firstfreqindex<0){
    retval=SFTEFIRSTINDEXNEG;
    goto error;
  }

  /* check that detector type is known */
  if (header.version!=1 && unknownDetector(header.detector))
    return SFTEINSTRUMENTUNKNOWN;

  /* if user has asked for comment, store it */
  if (comment && header.comment_length) {
    int i;
    
    /* first allocate memory for storing the comment */
    if (!(mycomment=calloc(header.comment_length, 1))) {
      retval=SFTENOMEM;
      goto error;
    }
    
    /* now read the comment into memory */
    if (header.comment_length != (int)fread(mycomment, 1, header.comment_length, fp)) {
      free(mycomment);
      mycomment = NULL;
      retval=SFTEREAD;
      goto error;
    }
    
    /* check that the comment is null-terminated and contains no
       messages 'hidden' after a NULL character */
    retval=SFTENONULLINCOMMENT;
    for (i=0; i<header.comment_length; i++) {
      if (!mycomment[i])
	retval=0;
      if (!retval && mycomment[i]) {
	retval=SFTEHIDDENCOMMENT;
	break;
      }
    }  
    if (retval) {
      free(mycomment);
      mycomment = NULL;
      goto error;
    }
  }
  
 error:
  /* restore stream pointer to the correct position */
  if (fsetpos(fp, &streamposition))
    retval=SFTERESTORESTREAMPOS;
  
  if (retval)
    return retval;

  /* no errors found: now set return values and return */
  *swapendian=swap;
  
  if (comment)
    /* note: mycomment may be NULL if there is NO comment */
    *comment=mycomment;

  /* return header-info */
  *info=header;

  return SFTNOERROR;
}



int ReadSFTData(FILE *fp,                /* data file.  Position left unchanged on return */
		float *data,             /* location where data should be written */
		int firstbin,            /* first frequency bin to read from data set */
		int nsamples,            /* number of frequency bin samples to retrieve */
		char **comment,          /* if non-NULL, will contain pointer to comment string */
		struct headertag2 *info  /* if non-NULL, will contain header information */
		) {
  fpos_t streamposition;
  int retval=0, swapendian;
  int seekforward;
  struct headertag2 myinfo;
  
  /* check for null pointers */
  if (!fp)
    return SFTENULLFP;

  if (!data && nsamples)
    return SFTENULLPOINTER;

  /* save position of stream */
  if (fgetpos(fp, &streamposition))
    return SFTEGETSTREAMPOS;

  /* read header of SFT */
  if ((retval=ReadSFTHeader(fp, &myinfo, comment, &swapendian, 0)))
    goto error;

  /* sanity checks -- do we ask for data before the first frequency bin */
  if (firstbin<myinfo.firstfreqindex) {
    retval=SFTEBEFOREDATA;
    goto error;
  }

  /* warning, after this subtraction firstbin contains the offset to
     the first requested bin! */
  firstbin -= myinfo.firstfreqindex;

  /* sanity checks -- do we ask for data after the last frequency bin */
  if ((firstbin+nsamples)> myinfo.nsamples) {
    retval=SFTEAFTERDATA;
    goto error;
  }

  /* seek to start of data (skip comment if present) */
  if (myinfo.version == 1)
    seekforward=sizeof(struct headertag1)+firstbin*2*sizeof(float);
  else
    seekforward=sizeof(struct headertag2)+myinfo.comment_length+firstbin*2*sizeof(float);
  if (fseek(fp, seekforward, SEEK_CUR)) {
    retval=SFTESEEK;
    goto error2;
  }

  /* read in the data */
  if (nsamples != (int)fread((void *)data, 2*sizeof(float), nsamples, fp)) {
    retval=SFTEREAD;
    goto error2;
  }

  /* byte reverse the data if necessary */
  if (swapendian) {
    int i;
    for (i=0; i<2*nsamples; i++)
      swap4((char *)(data+i));
  }
  
 error2:
  /* free storage for the comemnt string */
  if (retval && comment && *comment) {
    free(*comment);
    *comment=NULL;
  }
  
 error:
  /* return to starting position in stream */
  if (fsetpos(fp, &streamposition))
    return SFTERESTORESTREAMPOS;

  /* if no errors, return header information */
  if (info && !retval)
    *info=myinfo;
  
  return retval;
}

/* This routine returns zero if the two headers contain consistent
   information, else an error code if they are not consistent */
int CheckSFTHeaderConsistency(struct headertag2 *headerone, /* pointer to earlier header */
			      struct headertag2 *headertwo  /* pointer to later header */
			      ) {
  /* check for null pointer */
  if (!headerone || !headertwo)
    return SFTENULLPOINTER;
  
  /* Same version number */
  if (headerone->version != headertwo->version)
    return SFTEVERSIONCHANGES;

  /* GPS times increasing */
  if (headerone->gps_sec >headertwo->gps_sec || (headerone->gps_sec==headertwo->gps_sec && headerone->gps_nsec>=headertwo->gps_nsec))
    return SFTEGPSNOTINCREASING;
  
  /* Time base the same */
  if (headerone->tbase != headertwo->tbase)
    return SFTETBASECHANGES;
  
  /* First frequency index the same */
  if (headerone->firstfreqindex != headertwo->firstfreqindex)
    return SFTEFIRSTINDEXCHANGES;
  
  /* Number of samples the same */
  if (headerone->nsamples != headertwo->nsamples)
    return SFTENSAMPLESCHANGES;

  /* check for identical detectors */
  if ( (headerone->detector[0] != headertwo->detector[0]) || (headerone->detector[1] != headertwo->detector[1]) )
    return SFTEINSTRUMENTCHANGES;

  return SFTNOERROR;
}

/* check that channel-prefix defines a 'known' detector.  The list of
 * known detectors implemented here for now follows the list in
 * Appendix D of LIGO-T970130-F-E:
 *
 * returns 0 if known, error code otherwise */
int unknownDetector (const char *detector) {
  int i;
  const char *knownDetectors[] = {
    "A1",       /* ALLEGRO */
    "B1",       /* NIOBE */
    "E1",       /* EXPLORER */
    "G1",       /* GEO_600 */
    "H1",       /* LHO_4k */
    "H2",       /* LHO_2k */
    "K1",       /* ACIGA */
    "L1",       /* LLO_4k */
    "N1",       /* Nautilus */
    "O1",       /* AURIGA */
    "P1",       /* CIT_40 */
    "T1",       /* TAMA_300 */
    "V1",       /* Virgo_CITF */
    "V2",       /* Virgo (3km) */
    NULL
  };

  if (!detector)
    return SFTENULLPOINTER;

  for (i=0; knownDetectors[i]; i++) {
    if (knownDetectors[i][0]==detector[0] && knownDetectors[i][1]==detector[1])
      return 0;
  }
  
  return SFTEINSTRUMENTUNKNOWN;

} /* uknownDetector() */
