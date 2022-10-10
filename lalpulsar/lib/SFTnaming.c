/*
 * Copyright (C) 2010, 2014, 2016, 2022 Karl Wette
 * Copyright (C) 2010 Chris Messenger
 * Copyright (C) 2004, 2005 Reinhard Prix, Bernd Machenschalk, Alicia Sintes
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

/*---------- includes ----------*/

#include <ctype.h>

#include <lal/LISAspecifics.h>

#include "SFTinternal.h"

/*========== function definitions ==========*/

/**
 * Parses valid CW detector names and prefixes. 'name' input can be either a valid detector name or prefix
 *
 * \return XLAL_SUCCESS if 'name' is a valid detector name/prefix, an XLAL error code on failure
 *
 * If passed a non-NULL pointer 'prefix', will be set to the Allocated prefix string (2 characters+0)
 * for valid detectors, NULL otherwise.
 *
 * If passed a non-NULL pointer 'lalCachedIndex', will set to index >= 0 into the
 * lalCachedDetectors[] array if found there, or -1 if it's one of the "CW special" detectors
 *
 * \note This should be the *only* function defining valid CW detector names and prefixes
 *
 * \note 'prefix' is allocated here and needs to be free'ed by caller!
 *
 * \note Unless 'exactMatch' is true, if first two characters of 'name' match a valid detector prefix,
 * that is accepted, which allows passing longer strings beginning with a detector prefix ('H1:blabla')
 * without extra hassle
 */
int
XLALFindCWDetector ( CHAR** prefix, INT4 *lalCachedIndex, const CHAR *name, const BOOLEAN exactMatch )
{
  XLAL_CHECK ( name != NULL, XLAL_EINVAL );
  XLAL_CHECK ( strlen ( name ) >= 2, XLAL_EINVAL );	// need at least a full prefix 'letter+number'

  // ----- first check if 'name' corresponds to one of our 'CW special' detectors (LISA and X-ray satellites)
  const CHAR *specialDetectors[] =
    {
      "Z1",	  /* LISA effective IFO 1 */
      "Z2",	  /* LISA effective IFO 2 */
      "Z3",	  /* LISA effective IFO 3 */
      "Z4",	  /* LISA effective IFO 2 minus 3 */
      "Z5",	  /* LISA effective IFO 3 minus 1 */
      "Z6",	  /* LISA effective IFO 1 minus 2 */
      "Z7",	  /* LISA pseudo TDI A */
      "Z8",	  /* LISA pseudo TDI E */
      "Z9",	  /* LISA pseudo TDI T */
      NULL
    };
  for ( UINT4 i = 0; specialDetectors[i] != NULL; i ++ )
    {
      if ( strncmp ( specialDetectors[i], name, 2 ) == 0 )
        {
          if ( exactMatch && strcmp( name + strlen ( specialDetectors[i] ), "" ) != 0 ) {
            XLAL_ERROR ( XLAL_EINVAL, "Trailing garbage '%s' after detector prefix\n", name + strlen ( specialDetectors[i] ) );
          }
          if ( prefix != NULL ) {
            (*prefix) = XLALStringDuplicate ( specialDetectors[i] );
          }
          if ( lalCachedIndex != NULL ) {
            (*lalCachedIndex) = -1;
          }
          return XLAL_SUCCESS;
        }
    } // for i < len(specialDetectors)

  // ----- if not found, go through list of 'official' cached lalsuite detectors
  UINT4 numLALDetectors = sizeof(lalCachedDetectors) / sizeof(lalCachedDetectors[0]);
  for ( UINT4 i = 0; i < numLALDetectors; i ++)
    {
      const CHAR *prefix_i = lalCachedDetectors[i].frDetector.prefix;
      const CHAR *name_i   = lalCachedDetectors[i].frDetector.name;
      if ( strncmp ( prefix_i, name, 2 ) == 0 )
        {
          if ( exactMatch && strcmp( name + strlen ( prefix_i ), "" ) != 0 ) {
            XLAL_ERROR ( XLAL_EINVAL, "Trailing garbage '%s' after detector prefix\n", name + strlen ( prefix_i ) );
          }
          if ( prefix != NULL ) {
            (*prefix) = XLALStringDuplicate ( prefix_i );
          }
          if ( lalCachedIndex != NULL ) {
            (*lalCachedIndex) = i;
          }
          return XLAL_SUCCESS;
        }
      if ( strncmp ( name, name_i, strlen ( name_i ) ) == 0 )
        {
          if ( exactMatch && strcmp( name + strlen ( name_i ), "" ) != 0 ) {
            XLAL_ERROR ( XLAL_EINVAL, "Trailing garbage '%s' after detector prefix\n", name + strlen ( name_i ) );
          }
          if ( prefix != NULL ) {
            (*prefix) = XLALStringDuplicate ( prefix_i );
          }
          if ( lalCachedIndex != NULL ) {
            (*lalCachedIndex) = i;
          }
          return XLAL_SUCCESS;
        }

    } // for i < numLALDetectors

  if ( prefix != NULL ) {
    (*prefix) = NULL;
  }
  XLAL_ERROR ( XLAL_EINVAL, "Unknown detector name '%s'\n", name );

} // XLALFindCWDetector()


/**
 * Determine if 'name' is a valid detector name or prefix
 *
 * \return TRUE if 'name' is a valid detector name/prefix, FALSE otherwise
 * (or if an error is encountered; this function does not raise XLAL errors)
 */
BOOLEAN
XLALIsValidCWDetector ( const CHAR *name )
{
  int UNUSED errnum;
  int retn;
  XLAL_TRY_SILENT( retn = XLALFindCWDetector ( NULL, NULL, name, 1 ), errnum );
  return ( retn == XLAL_SUCCESS );
}


/**
 * Find the valid CW detector prefix. 'name' input can be either a valid detector name or prefix
 *
 * \note The returned string is allocated here!
 */
CHAR *
XLALGetChannelPrefix ( const CHAR *name )
{
  XLAL_CHECK_NULL ( name != NULL, XLAL_EINVAL );

  CHAR *prefix = NULL;
  XLAL_CHECK_NULL ( XLALFindCWDetector ( &prefix, NULL, name, 0 ) == XLAL_SUCCESS, XLAL_EINVAL, "Unknown detector name '%s'\n", name );

  return prefix;
} // XLALGetChannelPrefix()


/**
 * Find the site geometry-information 'LALDetector' for given a detector name (or prefix).
 *
 * \note The LALDetector struct is allocated here and needs to be free'ed by caller!
 */
LALDetector *
XLALGetSiteInfo ( const CHAR *name )
{
  XLAL_CHECK_NULL ( name != NULL, XLAL_EINVAL );

  const INT4 numLALDetectors = sizeof(lalCachedDetectors) / sizeof(lalCachedDetectors[0]);

  // first turn the free-form 'detector-name' into a well-defined channel-prefix, and find lalCachedDetector[] index
  INT4 lalCachedIndex = -1;
  CHAR *prefix;
  XLAL_CHECK_NULL ( XLALFindCWDetector ( &prefix, &lalCachedIndex, name, 0 ) == XLAL_SUCCESS, XLAL_EFUNC );

  LALDetector *site;
  XLAL_CHECK_NULL ( ( site = XLALCalloc ( 1, sizeof( *site) )) != NULL, XLAL_ENOMEM );

  switch ( prefix[0] )
    {
    case 'Z':       // create dummy-sites for LISA
      XLAL_CHECK_NULL ( XLALcreateLISA ( site, prefix[1] ) == XLAL_SUCCESS, XLAL_EFUNC, "Failed to created LISA detector 'name=%s, prefix=%s'\n", name, prefix );
      break;

    default:
      XLAL_CHECK_NULL ( (lalCachedIndex >= 0) && (lalCachedIndex < numLALDetectors), XLAL_EFAILED, "Internal inconsistency found (for 'name=%s, prefix=%s')\n", name, prefix );
      (*site) = lalCachedDetectors[lalCachedIndex];
      break;
    } /* switch channel[0] */

  XLALFree ( prefix );

  return site;

} // XLALGetSiteInfo()


/**
 * Return the 'official' file name for a given SFT, folllowing the SFT-v2 naming convention
 * LIGO-T040164-01 https://dcc.ligo.org/cgi-bin/DocDB/ShowDocument?docid=27385, namely
 *
 * name = S-D-G-T.sft
 * where
 * S = Source: upper-case single letter site designation 'G', 'H', 'L', 'V', ...
 * D = description: a free-form string of alphanumerics and {_, +, #}
 * G = GPS start time of first SFT in seconds (9- or 10-digit number)
 * T = total time interval covered by the data in this file
 *
 * furthermore, the v2-spec uses the following convention for the description field 'D':
 * D = numSFTs_IFO_SFTtype[_Misc]
 * where
 * numSFTs : number of SFTs in the file
 * IFO     : 2-character detector name, eg 'G1', 'H1', 'H2', 'L1', 'V1', ...
 * SFTtype : SFT-timebase, in the form '[T]SFT', where [T] is the SFT-duration in seconds, eg "1800SFT"
 * Misc    : optional string providing additional information
 */
char *
XLALOfficialSFTFilename ( char site,		//!< site-character 'G', 'H', 'L', ...
                          char channel,	//!< channel character '1', '2', ...
                          UINT4 numSFTs,	//!< number of SFTs in SFT-file
                          UINT4 Tsft,		//!< time-baseline in (integer) seconds
                          UINT4 GPS_start,	//!< GPS seconds of first SFT start time
                          UINT4 Tspan,		//!< total time-spanned by all SFTs in seconds
                          const char *Misc	//!< [in] optional 'Misc' entry in the SFT 'D' field (can be NULL)
                          )
{
  if ( Misc != NULL ) {
    XLAL_CHECK_NULL ( XLALCheckValidDescriptionField ( Misc ) == XLAL_SUCCESS, XLAL_EINVAL );
  }

  // ----- S
  char S[2] = { site, 0 };

  // ----- D
  char D[512];
  char IFO[2] = { site, channel };
  size_t written = snprintf ( D, sizeof(D), "%d_%c%c_%dSFT%s%s", numSFTs, IFO[0], IFO[1], Tsft, Misc ? "_" : "", Misc ? Misc : "" );
  XLAL_CHECK_NULL ( written < sizeof(D), XLAL_EINVAL, "Description field length of %zu exceeds buffer length of %zu characters\n", written, sizeof(D)-1 );

  // ----- G
  char G[11];
  written = snprintf ( G, sizeof(G), "%09d", GPS_start );
  XLAL_CHECK_NULL ( written < sizeof(G), XLAL_EINVAL, "GPS seconds %d exceed buffer length of %zu characters\n", GPS_start, sizeof(G)-1 );

  // ----- T
  char T[10];
  written = snprintf ( T, sizeof(T), "%d", Tspan );
  XLAL_CHECK_NULL ( written < sizeof(T), XLAL_EINVAL, "Tspan=%d s exceed buffer length of %zu characters\n", Tspan, sizeof(T)-1 );

  // S-D-G-T.sft
  size_t len = strlen(S) + 1 + strlen(D) + 1 + strlen(G) + 1 + strlen(T) + 4 + 1;
  char *filename;
  XLAL_CHECK_NULL ( (filename = XLALCalloc ( 1, len )) != NULL, XLAL_ENOMEM );

  written = snprintf ( filename, len, "%s-%s-%s-%s.sft", S, D, G, T );
  XLAL_CHECK_NULL ( written < len, XLAL_EFAILED, "Miscounted string-length, expected %zu characters but got %zu\n", len - 1, written );

  return filename;

} // XLALGetOfficialName4SFT()


/**
 * Return the 'official' file name for a given SFT, folllowing the SFT-v2 naming convention
 * LIGO-T040164-01 https://dcc.ligo.org/cgi-bin/DocDB/ShowDocument?docid=27385,
 * see also XLALOfficialSFTFilename() for details.
 */
char *
XLALGetOfficialName4SFT ( const SFTtype *sft,	//!< [in] input SFT to generate name for
                          const char *Misc	//!< [in] optional 'Misc' entry in the SFT 'D' field (can be NULL)
                          )
{
  XLAL_CHECK_NULL ( sft != NULL, XLAL_EINVAL );


  UINT4 Tsft = (UINT4) round ( 1.0 / sft->deltaF );

  /* calculate sft 'duration' -- may be different from timebase if nanosecond of sft-epoch is non-zero */
  UINT4 Tspan = Tsft;
  if ( sft->epoch.gpsNanoSeconds > 0) {
    Tspan += 1;
  }

  char *filename;
  XLAL_CHECK_NULL ( (filename = XLALOfficialSFTFilename ( sft->name[0], sft->name[1], 1, Tsft, sft->epoch.gpsSeconds, Tspan, Misc )) != NULL, XLAL_EFUNC );

  return filename;

} // XLALGetOfficialName4SFT()


/**
 * Return the 'official' file name for a given SFT-vector written into a single "merged SFT-file",
 * folllowing the SFT-v2 naming convention
 * LIGO-T040164-01 https://dcc.ligo.org/cgi-bin/DocDB/ShowDocument?docid=27385,
 * see also XLALOfficialSFTFilename() for details.
 */
char *
XLALGetOfficialName4MergedSFTs ( const SFTVector *sfts,	//!< [in] input SFT vector to generate name for
                                 const char *Misc	//!< [in] optional 'Misc' entry in the SFT 'D' field (can be NULL)
                                 )
{
  XLAL_CHECK_NULL ( sfts != NULL, XLAL_EINVAL );
  XLAL_CHECK_NULL ( sfts->length > 0, XLAL_EINVAL );

  UINT4 numSFTs = sfts->length;
  SFTtype *sftStart       = &(sfts->data[0]);
  SFTtype *sftEnd         = &(sfts->data[numSFTs-1]);
  LIGOTimeGPS *epochStart = &(sftStart->epoch);
  LIGOTimeGPS *epochEnd   = &(sftEnd->epoch);

  const char *name = sftStart->name;
  UINT4 Tsft = (UINT4) round ( 1.0 / sftStart->deltaF );

  /* calculate time interval covered -- may be different from timebase if nanosecond of sft-epochs are non-zero */
  UINT4 Tspan = epochEnd->gpsSeconds - epochStart->gpsSeconds + Tsft;
  if ( epochStart->gpsNanoSeconds > 0) {
    Tspan += 1;
  }
  if ( epochEnd->gpsNanoSeconds > 0) {
    Tspan += 1;
  }

  char *filename;
  XLAL_CHECK_NULL ( (filename = XLALOfficialSFTFilename ( name[0], name[1], numSFTs, Tsft, epochStart->gpsSeconds, Tspan, Misc )) != NULL, XLAL_EFUNC );

  return filename;

} // XLALGetOfficialName4MergedSFTs()


/**
 * Check whether given string qualifies as a valid 'description' field of a FRAME (or SFT)
 * filename, according to  LIGO-T010150-00-E "Naming Convention for Frame Files which are to be Processed by LDAS",
 * LIGO-T040164-01 at https://dcc.ligo.org/LIGO-T040164-x0/public
 *
 */
int
XLALCheckValidDescriptionField ( const char *desc )
{
  XLAL_CHECK ( desc != NULL, XLAL_EINVAL );

  size_t len = strlen ( desc );

  if ( len == 1 && isupper(desc[0]) ) {
    XLAL_ERROR ( XLAL_EINVAL, "Single uppercase description reserved for class-1 raw frames!\n" );
  }

  for ( UINT4 i=0; i < len; i ++ )
    {
      int c = desc[i];
      if ( !isalnum(c) && (c!='_') && (c!='+') && (c!='#') ) {	// all the valid characters allowed
        XLAL_ERROR ( XLAL_EINVAL, "Invalid chacter '%c' found, only alphanumeric and ['_', '+', '#'] are allowed\n", c );
      }
    } // for i < len

  return XLAL_SUCCESS;

} // XLALCheckValidDescriptionField()
