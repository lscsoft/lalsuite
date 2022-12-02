/*
 * Copyright (C) 2004, 2005, 2013 Reinhard Prix
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
 *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *  MA  02110-1301  USA
 */


/**
 * \author Reinhard Prix
 * \date 2005, 2013
 * \file
 * \ingroup lalpulsar_bin_SFTTools
 * \brief Code to dump various SFT-info in human-readable form to stdout.
 */

/* ---------- includes ---------- */
#include "config.h"

#include <lal/UserInput.h>
#include <lal/SFTfileIO.h>
#include <lal/PulsarDataTypes.h>
#include <lal/LALString.h>
#include <lal/LALPulsarVCSInfo.h>

/*---------- DEFINES ----------*/
#define MIN(x,y) ((x) < (y) ? (x) : (y))

/*----- Macros ----- */

/*---------- internal types ----------*/
/*---------- empty initializers ---------- */
/*---------- Global variables ----------*/

/* User variables */
typedef struct
{
  CHAR *SFTfiles;
  BOOLEAN headerOnly;
  BOOLEAN dataOnly;
  BOOLEAN timestampsOnly;
  UINT4 Nmax;
} UserVariables_t;

/*---------- internal prototypes ----------*/
int XLALprintDescriptor ( const SFTDescriptor *ptr );
int XLALprintHeader ( const SFTtype *header );
int XLALprintData ( const SFTtype *sft );
int XLALprintTimestamps ( const SFTCatalog *catalog, const UINT4 Nmax );

int XLALReadUserInput ( int argc, char *argv[], UserVariables_t *uvar );

/*==================== FUNCTION DEFINITIONS ====================*/

/*----------------------------------------------------------------------
 * main function
 *----------------------------------------------------------------------*/
int
main(int argc, char *argv[])
{
  /* register all our user-variable */
  UserVariables_t XLAL_INIT_DECL(uvar);
  XLAL_CHECK ( XLALReadUserInput ( argc, argv, &uvar ) == XLAL_SUCCESS, XLAL_EFUNC );

  SFTCatalog *catalog;
  XLAL_CHECK ( (catalog = XLALSFTdataFind ( uvar.SFTfiles, NULL )) != NULL, XLAL_EFUNC, "No SFTs matched your --SFTfiles query\n" );
  if ( XLALUserVarWasSet(&uvar.Nmax) && ( uvar.Nmax > catalog->length ) ) {
    XLALPrintWarning("Nmax=%u requested but SFT catalog has only %u entries. Returning them all.", uvar.Nmax, catalog->length);
  }

  if ( uvar.headerOnly  )
    {
      for ( UINT4 i=0; i < MIN(uvar.Nmax,catalog->length); i ++ )
        {
          SFTDescriptor *ptr = &(catalog->data[i]);

          XLALprintHeader ( &(ptr->header) );
          XLALprintDescriptor ( ptr );

        } // for i < catalog->length
    } // if --headerOnly
  else if ( uvar.timestampsOnly )
    {
      XLAL_CHECK ( XLALprintTimestamps ( catalog, uvar.Nmax ) == XLAL_SUCCESS, XLAL_EFUNC );
    } // if --timestampsOnly
  else // if --dataOnly or data+headers [default]
    {
      SFTVector *sfts;
      XLAL_CHECK ( (sfts = XLALLoadSFTs ( catalog, -1, -1 )) != NULL, XLAL_EFUNC );

      BOOLEAN segmentedSFTs = 0;
      if ( sfts->length < catalog->length ) {
        segmentedSFTs = 1;
        printf ("\n");
        printf ("%%%% Frequency-segmented SFTs: len(catalog) = %d, len(sft-vector) = %d\n", catalog->length, sfts->length );
        if ( !uvar.dataOnly ) {
          printf ("%%%% For segmented SFTS we can't output SFT 'descriptors' + data. Use --headerOnly if you need the descriptor fields\n\n");
        }
      } // if we're dealing with 'segmented SFTs': currently can't map them to catalog-descriptors

      for ( UINT4 i=0; i < MIN(uvar.Nmax,sfts->length); i ++ )
        {
          SFTtype *sft_i = &(sfts->data[i]);;

          if ( !uvar.dataOnly ) {
            XLALprintHeader ( sft_i );
            if ( !segmentedSFTs ) {	// skip descriptor fields for segmented SFTs (as we can't map them to SFTs)
              XLALprintDescriptor ( &(catalog->data[i]) );
            }
          } // output header info

          XLALprintData ( sft_i );

        } // for i < num_sfts

      XLALDestroySFTVector ( sfts );
    } /* if !headerOnly */

  /* free memory */
  XLALDestroySFTCatalog (catalog);
  XLALDestroyUserVars();

  LALCheckMemoryLeaks();

  return 0;
} /* main */

int
XLALprintDescriptor ( const SFTDescriptor *desc )
{
  XLAL_CHECK ( desc != NULL, XLAL_EINVAL );

  printf ("%%%% ----- Descriptor:\n");
  printf ( "Locator:     '%s'\n", XLALshowSFTLocator ( desc->locator ) );
  printf ( "SFT version: %d\n", desc->version );
  printf ( "numBins:     %d\n", desc->numBins );
  printf ( "crc64:       %" LAL_UINT8_FORMAT "\n", desc->crc64 );
  printf ( "window:      %s(%0.10g)\n", desc->window_type, desc->window_param );
  if ( desc->comment ) {
    if ( strchr( desc->comment, '\n' ) ) {
      printf ( "comment: ==========\n%s\nend comment: ------\n", desc->comment );
    } else {
      printf ( "comment:     %s\n", desc->comment );
    }
  } else {
    printf ( "comment:     <none>\n" );
  }

  return XLAL_SUCCESS;

} // XLALprintDescriptor()


int
XLALprintHeader ( const SFTtype *header )
{
  XLAL_CHECK ( header != NULL, XLAL_EINVAL );

  printf ("%%%% ----- Header:\n");
  printf ( "Name:        '%s'\n", header->name );
  printf ( "epoch:       [%d, %d]\n", header->epoch.gpsSeconds, header->epoch.gpsNanoSeconds );
  printf ( "f0:          %.9f\n", header->f0 );
  printf ( "deltaF:      %.9g\n", header->deltaF );
  if ( header->data ) {
    printf ("numBins:     %d\n", header->data->length );
  }

  return XLAL_SUCCESS;
} // XLALprintHeader()

int
XLALprintData ( const SFTtype *sft )
{
  XLAL_CHECK ( (sft != NULL) && (sft->data != NULL), XLAL_EINVAL );

  printf ("%%%% ----- Data x(f):\n");
  printf ("%%%% Frequency f[Hz]     Real(x)           Imag(x) \n");
  for ( UINT4 k=0; k < sft->data->length; k ++ ) {
            printf ( "%.9f      % 6e  % 6e  \n", sft->f0 + k * sft->deltaF, crealf(sft->data->data[k]), cimagf(sft->data->data[k]) );
  }

  return XLAL_SUCCESS;

} // XLALprintData()

// output timestamps in a format parseable by XLALReadTimestampsFile(), ie.
// "The timestamps file is of the format: <repeated lines of the form "seconds nano-seconds">
// allowing for '%#' as comments, which are ignored."
int
XLALprintTimestamps ( const SFTCatalog *catalog, const UINT4 Nmax )
{
  XLAL_CHECK ( (catalog != NULL) && (catalog->data != NULL), XLAL_EINVAL );

  // print fully descriptive commented header-log for reproducibility:
  char *version;
  XLAL_CHECK_MAIN ( (version = XLALVCSInfoString(lalPulsarVCSInfoList, 0, "%% ")) != NULL, XLAL_EFUNC );
  char *cmdline;
  XLAL_CHECK ( (cmdline = XLALUserVarGetLog ( UVAR_LOGFMT_CMDLINE )) != NULL, XLAL_EFUNC );

  char *logstr = NULL;
  XLAL_CHECK ( (logstr = XLALStringAppend ( logstr, version )) != NULL, XLAL_EFUNC );
  XLALFree ( version );
  XLAL_CHECK ( (logstr = XLALStringAppend ( logstr, "%% cmdline: " )) != NULL, XLAL_EFUNC );
  XLAL_CHECK ( (logstr = XLALStringAppend ( logstr, cmdline )) != NULL, XLAL_EFUNC );
  XLALFree ( cmdline );
  XLAL_CHECK ( (logstr = XLALStringAppend ( logstr, "\n" )) != NULL, XLAL_EFUNC );
  printf ( "%s", logstr );
  XLALFree ( logstr );

  // print timestamps
  printf ( "%%%% Timestamps: seconds nano-seconds\n" );
  for ( UINT4 i=0; i < MIN(Nmax,catalog->length); i ++ )
    {
      const SFTtype *header = &(catalog->data[i].header);
      printf ( "%10d %09d\n", header->epoch.gpsSeconds, header->epoch.gpsNanoSeconds );
    }

  return XLAL_SUCCESS;

} // XLALprintTimestamps()

int
XLALReadUserInput ( int argc, char *argv[], UserVariables_t *uvar )
{
  /* set a few defaults */
  uvar->Nmax = LAL_UINT4_MAX;

  XLALRegisterUvarMember(	SFTfiles,	STRING,  'i', REQUIRED, "File-pattern for input SFTs. Possibilities are:\n"
                                " - '<SFT file>;<SFT file>;...', where <SFT file> may contain wildcards\n - 'list:<file containing list of SFT files>'");
  XLALRegisterUvarMember(	headerOnly,	BOOLEAN, 'H', OPTIONAL, "Output only SFT headers");
  XLALRegisterUvarMember(	dataOnly,	BOOLEAN, 'd', OPTIONAL, "Output only SFT data, no header info");
  XLALRegisterUvarMember(	timestampsOnly,	BOOLEAN, 't', OPTIONAL, "Output only timestamps, in timestamps-file format");
  XLALRegisterUvarMember(	Nmax,	UINT4, 'N', OPTIONAL, "When run on multiple SFTs, exit after this many");

  /* read cmdline & cfgfile  */
  BOOLEAN should_exit = 0;
  XLAL_CHECK( XLALUserVarReadAllInput( &should_exit, argc, argv, lalPulsarVCSInfoList ) == XLAL_SUCCESS, XLAL_EFUNC );
  if ( should_exit ) {
    exit (1);
  }

  // ---------- sanity input checks ----------
  XLAL_CHECK ( UVAR_SET3( headerOnly, dataOnly, timestampsOnly ) <= 1, XLAL_EINVAL, "Contradictory input: at most *one* of --headerOnly, --dataOnly or --timestampsOnly allowed\n" );
  XLAL_CHECK ( uvar->Nmax > 0, XLAL_EINVAL );

  return XLAL_SUCCESS;

} // XLALReadUserInput()
