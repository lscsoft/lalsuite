/*
 * Copyright (C) 2010, 2014, 2016, 2020, 2022 Karl Wette
 * Copyright (C) 2010 Chris Messenger
 * Copyright (C) 2004--2006, 2010 Reinhard Prix
 * Copyright (C) 2004, 2005 Bernd Machenschalk
 * Copyright (C) 2004, 2005 Alicia Sintes
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with with program; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA  02110-1301  USA
 */

/**
 * \file
 * \ingroup SFTfileIO_h
 * \brief Internal SFT types and functions
 */

/*---------- includes ----------*/

#include <lal/SFTfileIO.h>

/*---------- SFT versions ----------*/

/**
 * \name SFT versions
 *
 * - Version 1 (--2004):
 *   - No formal specification.
 *   - No longer supported.
 * - Version 2 (2004--2022):
 *   - Specified in https://dcc.ligo.org/LIGO-T040164-v1/public.
 *   - Added header fields for detector name, checksum, comment string.
 *   - Incompatible with version 1.
 * - Version 3 (2022--):
 *   - Specified in https://dcc.ligo.org/LIGO-T040164-v2/public.
 *   - Modified header to record window type and parameter.
 *   - Compatible with version 2.
 */
/** @{ */

#define MIN_SFT_VERSION 2
#define MAX_SFT_VERSION 3

/** @} */

/*---------- constants ----------*/

#define TRUE    1
#define FALSE   0

static const REAL8 fudge_up   = 1 + 10 * LAL_REAL8_EPS;	// about ~1 + 2e-15
static const REAL8 fudge_down = 1 - 10 * LAL_REAL8_EPS;	// about ~1 - 2e-15

/** size of blocks allocated for SFT data. For Einstein\@home SFTs this should be set to 8000 (externally) */
#ifndef SFTFILEIO_REALLOC_BLOCKSIZE
#define SFTFILEIO_REALLOC_BLOCKSIZE 100
#endif

/*---------- macros ----------*/

#if defined(__GNUC__)
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

#define GPS2REAL8(gps) (1.0 * (gps).gpsSeconds + 1.e-9 * (gps).gpsNanoSeconds )
#define GPSEQUAL(gps1,gps2) (((gps1).gpsSeconds == (gps2).gpsSeconds) && ((gps1).gpsNanoSeconds == (gps2).gpsNanoSeconds))
#define GPSZERO(gps) (((gps).gpsSeconds == 0) && ((gps).gpsNanoSeconds == 0))

/*---------- internal types ----------*/

/* NOTE: the locator is implemented as an OPAQUE type in order to enforce encapsulation
 * of the actual physical storage of SFTs and to ease future extensions of the interface.
 * DO NOT TRY TO USE THIS TYPE OUTSIDE OF THIS MODULE!!
 */
struct tagSFTLocator
{
  CHAR *fname;		/* name of file containing this SFT */
  long offset;		/* SFT-offset with respect to a merged-SFT */
  UINT4 isft;           /* index of SFT this locator belongs to, used only in XLALLoadSFTs() */
};

/*---------- internal prototypes ----------*/

// These functions are defined in SFTtypes.c

REAL8 TSFTfromDFreq ( REAL8 dFreq );

int compareSFTdesc(const void *ptr1, const void *ptr2);
int compareSFTloc(const void *ptr1, const void *ptr2);
int compareDetNameCatalogs ( const void *ptr1, const void *ptr2 );
int compareSFTepoch(const void *ptr1, const void *ptr2);

// These functions are defined in SFTnaming.c

int build_sft_windowspec ( UINT2 *windowspec, CHAR (*windowspec_str)[9], const char *window_type, REAL8 window_param );
int parse_sft_windowspec ( const UINT2 windowspec, const char **window_type, REAL8 *window_param );
int parse_sft_windowspec_str ( const CHAR *windowspec_str, CHAR (*window_type)[32], REAL8 *window_param );

// These functions are defined in SFTfileIO.c

void endian_swap(CHAR * pdata, size_t dsize, size_t nelements);

FILE * fopen_SFTLocator ( const struct tagSFTLocator *locator );

int read_SFTversion_from_fp ( UINT4 *version, BOOLEAN *need_swap, FILE *fp );
int read_sft_header_from_fp (FILE *fp, SFTtype *header, UINT4 *version, UINT8 *crc64, UINT2 *SFTwindowspec, BOOLEAN *swapEndian, CHAR **SFTcomment, UINT4 *numBins );
UINT4 read_sft_bins_from_fp ( SFTtype *ret, UINT4 *firstBinRead, UINT4 firstBin2read, UINT4 lastBin2read , FILE *fp );

BOOLEAN has_valid_crc64 (FILE *fp );

// These functions are defined in SFTReferenceLibrary.c

unsigned long long crc64(const unsigned char* data, unsigned int length, unsigned long long crc);
