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
 *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 *  MA  02111-1307  USA
 */


/**
 * \author Reinhard Prix
 * \date 2005, 2013
 * \file
 * \ingroup pulsarApps
 * \brief Code to dump various SFT-info in human-readable form to stdout.
 */

/* ---------- includes ---------- */
#include <lalapps.h>

#include <lal/UserInput.h>
#include <lal/SFTfileIO.h>
#include <lal/SFTutils.h>

/*---------- DEFINES ----------*/

/*----- Macros ----- */
#define INIT_MEM(x) memset(&(x), 0, sizeof((x)))

/*---------- internal types ----------*/
/*---------- empty initializers ---------- */
/*---------- Global variables ----------*/

/* User variables */
typedef struct
{
  BOOLEAN help;
  CHAR *SFTfiles;
  BOOLEAN headerOnly;
  BOOLEAN noHeader;
} UserVariables_t;

/*---------- internal prototypes ----------*/
int XLALprintDescriptor ( const SFTDescriptor *ptr );
int XLALprintHeader ( const SFTtype *header );
int XLALprintData ( const SFTtype *sft );

int XLALReadUserInput ( int argc, char *argv[], UserVariables_t *uvar );

/*==================== FUNCTION DEFINITIONS ====================*/

/*----------------------------------------------------------------------
 * main function
 *----------------------------------------------------------------------*/
int
main(int argc, char *argv[])
{
  /* register all our user-variable */
  UserVariables_t uvar; INIT_MEM ( uvar );
  XLAL_CHECK ( XLALReadUserInput ( argc, argv, &uvar ) == XLAL_SUCCESS, XLAL_EFUNC );

  SFTConstraints constraints; INIT_MEM ( constraints );
  CHAR detector[2] = "??";	/* allow reading v1-SFTs without detector-info */
  constraints.detector = detector;
  SFTCatalog *catalog;
  XLAL_CHECK ( (catalog = XLALSFTdataFind ( uvar.SFTfiles, &constraints )) != NULL, XLAL_EFUNC, "No SFTs matched your --SFTfiles query\n" );

  if ( uvar.headerOnly )
    {
      for ( UINT4 i=0; i < catalog->length; i ++ )
        {
          SFTDescriptor *ptr = &(catalog->data[i]);

          XLALprintHeader ( &(ptr->header) );
          XLALprintDescriptor ( ptr );

        } // for i < catalog->length
    } // if header
  else
    {
      SFTVector *sfts;
      XLAL_CHECK ( (sfts = XLALLoadSFTs ( catalog, -1, -1 )) != NULL, XLAL_EFUNC );

      BOOLEAN segmentedSFTs = 0;
      if ( sfts->length < catalog->length ) {
        segmentedSFTs = 1;
        printf ("\n");
        printf ("%%%% Frequency-segmented SFTs: len(catalog) = %d, len(sft-vector) = %d\n", catalog->length, sfts->length );
        if ( !uvar.noHeader ) {
          printf ("%%%% For segmented SFTS we can't output SFT 'descriptors' + data. Use --headerOnly if you need the descriptor fields\n\n");
        }
      } // if we're dealing with 'segmented SFTs': currently can't map them to catalog-descriptors

      for ( UINT4 i=0; i < sfts->length; i ++ )
        {
          SFTtype *sft_i = &(sfts->data[i]);;

          if ( !uvar.noHeader ) {
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
  printf ( "comment:     %s\n", (desc->comment)?(desc->comment) : "<none>" );

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


int
XLALReadUserInput ( int argc, char *argv[], UserVariables_t *uvar )
{
  XLALregBOOLUserStruct ( 	help,		'h', UVAR_HELP,     "Print this help/usage message");
  XLALregSTRINGUserStruct (	SFTfiles,	'i', UVAR_REQUIRED, "File-pattern for input SFTs");
  XLALregBOOLUserStruct (	headerOnly,	'H', UVAR_OPTIONAL, "Output only header-info");
  XLALregBOOLUserStruct (	noHeader,	'n', UVAR_OPTIONAL, "Output only data, no header");

  /* read cmdline & cfgfile  */
  XLAL_CHECK ( XLALUserVarReadAllInput ( argc, argv ) == XLAL_SUCCESS, XLAL_EFUNC );

  if (uvar->help) { 	/* help requested: we're done */
    exit (0);
  }

  XLAL_CHECK ( !(uvar->headerOnly && uvar->noHeader), XLAL_EINVAL, "Contradictory input --headerOnly --noHeader\n" );

  return XLAL_SUCCESS;
}
