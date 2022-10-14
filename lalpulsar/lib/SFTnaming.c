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

#define _XOPEN_SOURCE   // for isascii()

#include <ctype.h>

#include <lal/LISAspecifics.h>

#include "SFTinternal.h"

/*---------- macros ----------*/

#define LAST_ELEMENT(buf) ( (buf)[XLAL_NUM_ELEM(buf) - 1] )

#define STRCPYCHK(name, dest, src) do { \
    const char *const STRCPYCHK_p = (src); \
    XLAL_CHECK( strlen(STRCPYCHK_p) < sizeof(dest), XLAL_ESIZE, "%s '%s' does not fit in a buffer of size %zu", name, STRCPYCHK_p, sizeof(dest) ); \
    strncpy((dest), STRCPYCHK_p, sizeof(dest)); \
  } while(0)

/*---------- constants ----------*/

/** SFT window specification; see \cite SFT-spec */
static const struct {
  const char *window_type;
  const char *short_name;
  UINT2 A;
  UINT2 B;
} windowspec_table[] = {
  /* zero-parameter windows */
  { .window_type = "unknown",     .short_name = "UNKN", .B = 0 },
  { .window_type = "rectangular", .short_name = "RECT", .B = 1 },
  { .window_type = "hann",        .short_name = "HANN", .B = 2 },
  /* one-parameter windows */
  { .window_type = "tukey",       .short_name = "TKEY", .A = 1 },
};

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
 * Build an SFT file name from the given specification
 */
char *XLALBuildSFTFilenameFromSpec(
  const SFTFilenameSpec *spec   /**< [in] SFT filename specification */
  )
{

  // check input
  XLAL_CHECK_NULL( spec != NULL, XLAL_EFAULT );

  // check specification
  XLAL_CHECK_NULL( LAST_ELEMENT(spec->path) == 0, XLAL_EINVAL,
                   "'path' is not a null-terminated string" );
  XLAL_CHECK_NULL( LAST_ELEMENT(spec->extn) == 0, XLAL_EINVAL,
                   "'extn' is not a null-terminated string" );
  XLAL_CHECK_NULL( spec->numSFTs > 0, XLAL_EINVAL,
                   "'numSFTs' must be strictly positive" );
  XLAL_CHECK_NULL( LAST_ELEMENT(spec->detector) == 0, XLAL_EINVAL,
                   "'detector' is not a null-terminated string" );
  XLAL_CHECK_NULL( XLALIsValidCWDetector( spec->detector ), XLAL_EINVAL,
                   "'%s' is not a valid 2-character detector prefix", spec->detector );
  XLAL_CHECK_NULL( spec->SFTtimebase > 0, XLAL_EINVAL,
                   "'SFTtimebase' must be strictly positive" );
  XLAL_CHECK_NULL( LAST_ELEMENT(spec->window_type) == 0, XLAL_EINVAL,
                   "'window_type' is not a null-terminated string" );
  XLAL_CHECK_NULL( LAST_ELEMENT(spec->privMisc) == 0, XLAL_EINVAL,
                   "'privMisc' is not a null-terminated string" );
  XLAL_CHECK_NULL( spec->pubObsRun == 0 || strlen(spec->privMisc) == 0, XLAL_EINVAL,
                   "Public SFTs (with pubObsRun=%u) cannot include a private description '%s'",
                   spec->pubObsRun, spec->privMisc );
  XLAL_CHECK_NULL( strlen(spec->privMisc) == 0 || XLALCheckValidDescriptionField(spec->privMisc) == XLAL_SUCCESS, XLAL_EINVAL,
                   "Private description '%s' contains disallowed characters", spec->privMisc );
  XLAL_CHECK_NULL( LAST_ELEMENT(spec->pubObsKind) == 0, XLAL_EINVAL,
                   "'pubObsKind' is not a null-terminated string" );
  XLAL_CHECK_NULL( spec->pubObsRun == 0 ||
                   strcmp(spec->pubObsKind, "RUN") == 0 ||
                   strcmp(spec->pubObsKind, "AUX") == 0 ||
                   strcmp(spec->pubObsKind, "SIM") == 0 ||
                   strcmp(spec->pubObsKind, "DEV") == 0, XLAL_EINVAL,
                   "'pubObsKind=%s' must be one of 'RUN'|'AUX'|'SIM'|'DEV'", spec->pubObsKind );
  XLAL_CHECK_NULL( spec->pubObsRun == 0 || spec->pubVersion > 0, XLAL_EINVAL,
                   "Public SFTs (with pubObsRun=%u) must include a positive version 'pubVersion'",
                   spec->pubObsRun );
  XLAL_CHECK_NULL( LAST_ELEMENT(spec->pubChannel) == 0, XLAL_EINVAL,
                   "'pubChannel' is not a null-terminated string" );
  XLAL_CHECK_NULL( spec->pubObsRun == 0 || strlen(spec->pubChannel) > 0, XLAL_EINVAL,
                   "Public SFTs (with pubObsRun=%u) must include a channel name 'pubChannel'",
                   spec->pubObsRun );
  XLAL_CHECK_NULL( spec->pubObsRun == 0 || strlen(spec->window_type) > 0, XLAL_EINVAL,
                   "Public SFTs (with pubObsRun=%u) must include a window function 'window_type'",
                   spec->pubObsRun );
  XLAL_CHECK_NULL( spec->nbFirstBinFreq == 0 || spec->nbBinWidthFreq > 0, XLAL_EINVAL,
                   "Narrow-band SFTs (with nbFirstBinFreq>0) must have nbBinWidthFreq>0" );
  XLAL_CHECK_NULL( spec->gpsStart > 0, XLAL_EINVAL,
                   "'gpsStart' must be strictly positive" );
  XLAL_CHECK_NULL( spec->SFTspan > 0, XLAL_EINVAL,
                   "'SFTspan' must be strictly positive" );

  char windowspec_str[9];
  if ( spec->pubObsRun > 0 ) {
    XLAL_CHECK_NULL( build_sft_windowspec( NULL, &windowspec_str, spec->window_type, spec->window_beta ) == XLAL_SUCCESS, XLAL_EINVAL,
                     "'%s' is not a valid SFT window function", spec->window_type );
  }

  // SFT filename string
  char *fname = NULL;

  // append SFT path, if any
  if (strlen(spec->path) > 0) {
    fname = XLALStringAppendFmt(
      fname,
      "%s/",
      spec->path
      );
    XLAL_CHECK_NULL( fname != NULL, XLAL_EFUNC );
  }

  // append SFT filename fields <S>, <D1>, <D2>, <D3> (see \cite SFT-spec)
  fname = XLALStringAppendFmt(
    fname,
    "%c-%d_%s_%dSFT",
    spec->detector[0],
    spec->numSFTs,
    spec->detector,
    spec->SFTtimebase
    );
  XLAL_CHECK_NULL( fname != NULL, XLAL_EFUNC );

  // append optional SFT filename field <D4> (see \cite SFT-spec)
  if ( spec->pubObsRun == 0 ) {   // private SFT

    if (strlen(spec->privMisc) > 0) {   // miscellaneous description

      // append private SFT miscellaneous description
      fname = XLALStringAppendFmt(
        fname,
        "_%s",
        spec->privMisc
        );
      XLAL_CHECK_NULL( fname != NULL, XLAL_EFUNC );

    }

  } else {   // public SFT

    // strip prefix and disallowed characters from channel name
    char pubChannel[XLAL_NUM_ELEM(spec->pubChannel)];
    strcpy(pubChannel, spec->pubChannel);
    char *pubChannel_no_prefix = strchr(pubChannel, ':');
    if (!pubChannel_no_prefix) {
      pubChannel_no_prefix = pubChannel;
    }
    XLAL_CHECK_NULL( XLALStringKeepChars(pubChannel_no_prefix, isascii) != NULL, XLAL_EFUNC );
    XLAL_CHECK_NULL( XLALStringKeepChars(pubChannel_no_prefix, isalnum) != NULL, XLAL_EFUNC );

    // append public SFT descriptor
    fname = XLALStringAppendFmt(
      fname,
      "_O%d%s+V%d+C%s+W%s",
      spec->pubObsRun,
      spec->pubObsKind,
      spec->pubVersion,
      pubChannel_no_prefix,
      windowspec_str
      );
    XLAL_CHECK_NULL( fname != NULL, XLAL_EFUNC );

  }

  // append optional SFT filename field <D5> (see \cite SFT-spec)
  if (spec->nbFirstBinFreq > 0) {   // narrow-band SFTs

    // append narrow-band SFT descriptor
    fname = XLALStringAppendFmt(
      fname,
      "_NBF%04dHz%dW%04dHz%d",
      spec->nbFirstBinFreq,
      spec->nbFirstBinRem,
      spec->nbBinWidthFreq,
      spec->nbBinWidthRem
      );
    XLAL_CHECK_NULL( fname != NULL, XLAL_EFUNC );

  }

  // append optional SFT filename fields <G>, <T> (see \cite SFT-spec), and extension
  fname = XLALStringAppendFmt(
    fname,
    "-%09d-%d.%s",
    spec->gpsStart,
    spec->SFTspan,
    strlen(spec->extn) > 0 ? spec->extn : "sft"
    );
  XLAL_CHECK_NULL( fname != NULL, XLAL_EFUNC );

  return fname;

} // XLALBuildSFTFilenameFromSpec(()


/**
 * Parse a SFT file path and return its specification
 */
int XLALParseSFTFilenameIntoSpec(
  SFTFilenameSpec *spec,        /**< [out] SFT filename specification */
  const char *SFTpath           /**< [in] SFT file path */
  )
{

  // check input
  XLAL_CHECK( spec != NULL, XLAL_EFAULT );
  XLAL_CHECK( SFTpath != NULL, XLAL_EFAULT );
  XLALPrintInfo("%s(): SFTpath='%s'\n",
                __func__, SFTpath);

  // clear specification
  memset(spec, 0, sizeof(*spec));

  // local copy of path
  char localSFTpath[8192];
  STRCPYCHK( "SFT filename and path", localSFTpath, SFTpath );

  // split SFT filename
  char *SFTfname = strrchr(localSFTpath, '/');
  if ( !SFTfname ) {
    SFTfname = localSFTpath;
  } else {
    *SFTfname++ = '\0';
    STRCPYCHK( "SFT path", spec->path, localSFTpath );
  }
  XLALPrintInfo("%s(): spec->path='%s' SFTfname='%s'\n",
                __func__, spec->path, SFTfname);

  // split SFT filename <S>, <D>, <G>, <T> fields (see \cite SFT-spec)
  char *S, *D, *G, *T;
  {
    char *p = SFTfname;
    S = XLALStringToken(&p, "-", 1);
    XLAL_CHECK( S != NULL, XLAL_EINVAL,
                "SFT file path '%s' contains no <S> field", SFTpath );
    D = XLALStringToken(&p, "-", 1);
    XLAL_CHECK( D != NULL, XLAL_EINVAL,
                "SFT file path '%s' contains no <D> field", SFTpath );
    G = XLALStringToken(&p, "-", 1);
    XLAL_CHECK( G != NULL, XLAL_EINVAL,
                "SFT file path '%s' contains no <G> field", SFTpath );
    T = XLALStringToken(&p, ".", 1);
    XLAL_CHECK( T != NULL, XLAL_EINVAL,
                "SFT file path '%s' contains no <T> field", SFTpath );
    XLALPrintInfo("%s(): S='%s' D='%s' G='%s' T='%s' extn='%s'\n",
                  __func__, S, D, G, T, p);
    STRCPYCHK( "SFT extension", spec->extn, p );
  }

  // split SFT filename <D> field into <D1>, <D2>, <D3>[, <D4>, <D5>] fields (see \cite SFT-spec)
  char *D1, *D2, *D3, *D4, *D5;
  {
    char *p = D;
    D1 = XLALStringToken(&p, "_", 1);
    XLAL_CHECK( D1 != NULL, XLAL_EINVAL,
                "SFT file path '%s' contains no <D1> field", SFTpath );
    D2 = XLALStringToken(&p, "_", 1);
    XLAL_CHECK( D2 != NULL, XLAL_EINVAL,
                "SFT file path '%s' contains no <D2> field", SFTpath );
    D3 = XLALStringToken(&p, "_", 1);
    XLAL_CHECK( D3 != NULL, XLAL_EINVAL,
                "SFT file path '%s' contains no <D3> field", SFTpath );
    D4 = XLALStringToken(&p, "_", 1);
    D5 = XLALStringToken(&p, "_", 1);
    XLALPrintInfo("%s(): D1='%s' D2='%s' D3='%s' D4='%s' D5='%s' p='%s'\n",
                  __func__, D1, D2, D3, D4, D5, p);
  }

  // parse required SFT filename fields <S>, <D1>, <D2>, <D3>, <G>, <T>
  XLAL_CHECK( sscanf( D1, "%d", &spec->numSFTs ) == 1, XLAL_EFUNC,
              "Could not parse 'numSFTs' from field D1='%s'", D1 );
  XLAL_CHECK( spec->numSFTs > 0, XLAL_EINVAL,
              "'numSFTs' must be strictly positive" );
  STRCPYCHK( "SFT detector prefix", spec->detector, D2 );
  XLAL_CHECK( XLALIsValidCWDetector( spec->detector ), XLAL_EINVAL,
              "'%s' is not a valid 2-character detector prefix", spec->detector );
  XLAL_CHECK( S[0] == spec->detector[0] && S[1] == '\0', XLAL_EINVAL,
              "Inconsistent site/detector fields S='%s' and D2='%s'", S, D2 );
  XLAL_CHECK( sscanf( D3, "%dSFT", &spec->SFTtimebase ) == 1, XLAL_EFUNC,
              "Could not parse 'SFTtimebase' from field D3='%s'", D3 );
  XLAL_CHECK( spec->SFTtimebase > 0, XLAL_EINVAL,
              "'SFTtimebase' must be strictly positive" );
  XLAL_CHECK( sscanf( G, "%d", &spec->gpsStart ) == 1, XLAL_EFUNC,
              "Could not parse 'gpsStart' from field G='%s'", G );
  XLAL_CHECK( spec->gpsStart > 0, XLAL_EINVAL,
              "'gpsStart' must be strictly positive" );
  XLAL_CHECK( sscanf( T, "%d", &spec->SFTspan ) == 1, XLAL_EFUNC,
              "Could not parse 'SFTspan' from field T='%s'", T );
  XLAL_CHECK( spec->SFTspan > 0, XLAL_EINVAL,
              "'SFTspan' must be strictly positive" );

  // parse optional SFT filename field <D4>
  if (D4 != NULL) {
    if (strchr(D4, '+')) {   // public SFT
      char *D41, *D42, *D43, *D44;
      char *p = D4;
      D41 = XLALStringToken(&p, "+", 1);
      XLAL_CHECK( D41 != NULL, XLAL_EINVAL,
                  "SFT file path '%s' contains no <D41> field", SFTpath );
      D42 = XLALStringToken(&p, "+", 1);
      XLAL_CHECK( D42 != NULL, XLAL_EINVAL,
                  "SFT file path '%s' contains no <D42> field", SFTpath );
      D43 = XLALStringToken(&p, "+", 1);
      XLAL_CHECK( D43 != NULL, XLAL_EINVAL,
                  "SFT file path '%s' contains no <D43> field", SFTpath );
      D44 = XLALStringToken(&p, "+", 1);
      XLAL_CHECK( D44 != NULL, XLAL_EINVAL,
                  "SFT file path '%s' contains no <D44> field", SFTpath );
      XLALPrintInfo("%s(): D41='%s' D42='%s' D43='%s' D44='%s' p='%s'\n",
                    __func__, D41, D42, D43, D44, p);
      char buf[8192];
      XLAL_CHECK( sscanf( D41, "O%d%255s", &spec->pubObsRun, buf ) == 2, XLAL_EINVAL,
                  "Could not parse public SFT fields 'pubObsRun' and 'pubObsKind' from field D41='%s'", D41 );
      STRCPYCHK( "Public SFT field 'pubObsKind'", spec->pubObsKind, buf );
      XLAL_CHECK( sscanf( D42, "V%d", &spec->pubVersion ) == 1, XLAL_EFUNC,
                  "Could not parse public SFT field 'pubVersion' from field D42='%s'", D42 );
      XLAL_CHECK( sscanf( D43, "C%255s", buf ) == 1, XLAL_EINVAL,
                  "Could not parse public SFT field 'pubChannel' from field D43='%s'", D43 );
      STRCPYCHK( "Public SFT field 'pubChannel'", spec->pubChannel, buf );
      XLAL_CHECK( sscanf( D44, "W%255s", buf ) == 1, XLAL_EINVAL,
                  "Could not parse public SFT fields 'window_type' and 'window_beta' from field D44='%s'", D44 );
      XLAL_CHECK( parse_sft_windowspec_str( buf, &spec->window_type, &spec->window_beta ) == XLAL_SUCCESS, XLAL_EINVAL,
                  "'%s' does not specify a valid SFT window function", buf );
    } else { // private SFT
      XLAL_CHECK( XLALCheckValidDescriptionField(D4) == XLAL_SUCCESS, XLAL_EINVAL,
                  "Private description '%s' contains disallowed characters", D4 );
      STRCPYCHK( "Private SFT miscellaneous description", spec->privMisc, D4 );
    }
  }

  // parse optional SFT filename field <D5>
  if (D5 != NULL) {   // narrow-band SFTs
    XLAL_CHECK( sscanf( D5, "NBF%dHz%dW%dHz%d",
                        &spec->nbFirstBinFreq, &spec->nbFirstBinRem,
                        &spec->nbBinWidthFreq, &spec->nbBinWidthRem
                  ) == 4, XLAL_EINVAL,
                "Could not parse narrow-band SFT field D5='%s'", D5 );
  }

  return XLAL_SUCCESS;

} // XLALParseSFTFilenameIntoSpec()


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


/**
 * Build an SFT 2-byte 'windowspec' or filename field 'windowspec_str' for the window given by 'window_type' and 'window_beta'
 */
int build_sft_windowspec ( UINT2 *windowspec, CHAR (*windowspec_str)[9], const char *window_type, REAL8 window_beta )
{

  // check input
  XLAL_CHECK ( window_type != NULL, XLAL_EFAULT );

  // default to unknown window
  if (windowspec) {
    *windowspec = windowspec_table[0].B;
  }
  if (windowspec_str) {
    XLAL_CHECK ( snprintf( *windowspec_str, sizeof(*windowspec_str), "%s", windowspec_table[0].short_name ) < (int) sizeof(*windowspec_str), XLAL_ESYS );
  }

  // lookup windows
  for (size_t w = 0; w < XLAL_NUM_ELEM(windowspec_table); ++w) {
    if ( XLALStringCaseCompare(window_type, windowspec_table[w].window_type) == 0 ) {

      // build windowspec and windowspec_str
      if ( windowspec_table[w].A == 0 ) {   /* zero-parameter windows */

        XLAL_CHECK ( window_beta <= 0, XLAL_EINVAL, "Invalid window_beta=%0.10g for SFT window type '%s'", window_beta, window_type );

        if (windowspec) {
          *windowspec = windowspec_table[w].B;
        }

        if (windowspec_str) {
          XLAL_CHECK ( snprintf( *windowspec_str, sizeof(*windowspec_str), "%s", windowspec_table[w].short_name ) < (int) sizeof(*windowspec_str), XLAL_ESYS );
        }

      } else {                              /* one-parameter windows */

        REAL8 Breal = window_beta * 5000;
        long B = lrint(Breal);
        XLAL_CHECK ( B == ((long) Breal), XLAL_ETOL, "SFT window_beta=%0.10g cannot be exactly represented by an integer [0,5000]", window_beta );
        XLAL_CHECK ( 0 <= B && B <= 5000, XLAL_ERANGE, "SFT window_beta=%0.10g is out of range [0.0,1.0] (B=%ld)", window_beta, B );

        if (windowspec) {
          *windowspec = windowspec_table[w].A * 5001 + B;
        }

        if (windowspec_str) {
          XLAL_CHECK ( snprintf( *windowspec_str, sizeof(*windowspec_str), "%s%ld", windowspec_table[w].short_name, B ) < (int) sizeof(*windowspec_str), XLAL_ESYS );
        }

      }

      break;

    }
  }

  return XLAL_SUCCESS;

}


/**
 * Parse an SFT 2-byte 'windowspec' into a window name 'window_type' and possible parameter 'window_beta'
 */
int parse_sft_windowspec ( const UINT2 windowspec, const char **window_type, REAL8 *window_beta )
{

  // parse windowspec
  const UINT2 A = windowspec / 5001;
  const UINT2 B = windowspec - 5001 * A;

  // lookup windows
  for (size_t w = 0; w < XLAL_NUM_ELEM(windowspec_table); ++w) {
    if ( A == windowspec_table[w].A ) {

      // build window_type and window_beta
      if ( A == 0 && B == windowspec_table[w].B ) {   /* zero-parameter windows */

        if (window_type) {
          *window_type = windowspec_table[w].window_type;
        }

        if (window_beta) {
          *window_beta = 0;
        }

        return XLAL_SUCCESS;

      } else if ( A > 0 ) {                              /* one-parameter windows */

        if (window_type) {
          *window_type = windowspec_table[w].window_type;
        }

        if (window_beta) {
          *window_beta = ((REAL8) B) / 5000;
        }

        return XLAL_SUCCESS;

      }

    }
  }

  XLAL_ERROR ( XLAL_EINVAL, "Illegal SFT windowspec=%d : windowspecA=%d, windowspecB=%d",
               windowspec, A, B );

}


/**
 * Parse an SFT filename field 'windowspec_str' into a window name 'window_type' and possible parameter 'window_beta'
 */
int parse_sft_windowspec_str ( const CHAR *windowspec_str, CHAR (*window_type)[32], REAL8 *window_beta )
{

  // check input
  XLAL_CHECK ( windowspec_str != NULL, XLAL_EFAULT );

  // parse windowspec_str
  char short_name[5];
  long B = 0;
  int parsed = sscanf( windowspec_str, "%4s%li", short_name, &B );
  XLAL_CHECK ( 0 < parsed, XLAL_EINVAL, "Could not parse SFT filename field windowspec_str='%s'", windowspec_str );
  XLAL_CHECK ( 0 <= B && B <= 5000, XLAL_ERANGE, "SFT window_beta=%0.10g is out of range [0.0,1.0] (B=%ld)", ((REAL8) B) / 5000, B );

  // lookup windows
  for (size_t w = 0; w < XLAL_NUM_ELEM(windowspec_table); ++w) {
    if ( strcmp(short_name, windowspec_table[w].short_name) == 0 ) {

      // build window_type
      if (window_type) {
        XLAL_CHECK ( snprintf( *window_type, sizeof(*window_type), "%s", windowspec_table[w].window_type ) < (int) sizeof(*window_type), XLAL_ESYS );
      }

      // build window_beta
      if (window_beta) {
        if ( windowspec_table[w].A == 0 ) {   /* zero-parameter windows */
          *window_beta = 0;
        } else {                              /* one-parameter windows */
          *window_beta = ((REAL8) B) / 5000;
        }
      }

      return XLAL_SUCCESS;

    }
  }

  XLAL_ERROR ( XLAL_EINVAL, "Illegal SFT windowspec_str='%s'", windowspec_str );

}
