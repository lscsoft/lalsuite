 /*
 * Copyright (C) 2020 Pep Covas, David Keitel
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
 * \file
 * \ingroup lalapps_WriteSFTsfromSFDBs
 * \author Pep Covas, David Keitel
 * \brief Read data in the form of Virgo Collaboration's SFDBs, convert to SFTs in the usual LALSuite format, and write those to disk.
 */

 #include <sys/stat.h>
 #include <sys/types.h>
 #include <LALAppsVCSInfo.h>
 #include <lal/UserInput.h>
 #include <lal/PulsarDataTypes.h>
 #include <lal/SFTfileIO.h>
 #include <lal/SFTutils.h>
 #include <lal/LogPrintf.h>

typedef struct
{
    REAL8 fmin;
    REAL8 fmax;
    CHAR *file_pattern;
    CHAR *timeStampsStarting;
    CHAR *timeStampsFinishing;
    CHAR *outSFTnames;
} UserInput_t;

int initUserVars ( UserInput_t *uvar );

int main(int argc, char *argv[]) {

    UserInput_t XLAL_INIT_DECL(uvar);

    /* register all user-variables */
    XLAL_CHECK_MAIN ( initUserVars ( &uvar ) == XLAL_SUCCESS, XLAL_EFUNC );

    /* read cmdline & cfgfile  */
    BOOLEAN should_exit = 0;
    XLAL_CHECK_MAIN( XLALUserVarReadAllInput( &should_exit, argc, argv, lalAppsVCSInfoList ) == XLAL_SUCCESS, XLAL_EFUNC );
    if ( should_exit ) {
      exit (1);
    }

    MultiSFTVector *inputSFTs = NULL;
    XLAL_CHECK_MAIN ( (inputSFTs = XLALReadSFDB(uvar.fmin, uvar.fmax, uvar.file_pattern, uvar.timeStampsStarting, uvar.timeStampsFinishing)) != NULL, XLAL_EFUNC );

    LALStringVector *fnames;
    XLAL_CHECK_MAIN ( (fnames = XLALFindFiles (uvar.outSFTnames)) != NULL, XLAL_EFUNC, "Failed to find filelist matching pattern '%s'.\n\n", uvar.outSFTnames );


    for ( UINT4 k = 0; k < inputSFTs->length; k++) {

        FILE *fpSingleSFT = NULL;
        const CHAR *filenameOut = fnames->data[k];
        XLAL_CHECK_MAIN ( ( fpSingleSFT = fopen ( filenameOut, "wb" )) != NULL, XLAL_EIO, "Failed to open singleSFT file '%s' for writing\n", filenameOut );

        UINT4 numSFTs = inputSFTs->data[k]->length;
        for ( UINT4 i=0; i < numSFTs; i++ )
        {   
            const SFTtype *thisSFT = inputSFTs->data[k]->data + i;

            XLAL_CHECK_MAIN ( XLALWriteSFT2fp( thisSFT, fpSingleSFT, NULL ) == XLAL_SUCCESS, XLAL_EFUNC,  "XLALWriteSFT2fp() failed to write SFT to '%s'!\n", filenameOut );
                 
        } /* for i < numSFTs */

        if ( fpSingleSFT ) fclose ( fpSingleSFT );

    }

    XLALDestroyMultiSFTVector (inputSFTs); 
    XLALDestroyStringVector(fnames);

    return 0;

}


int initUserVars ( UserInput_t *uvar ) {
   
   XLAL_CHECK ( uvar != NULL, XLAL_EINVAL );

   uvar->fmin = 0;
   uvar->fmax = 0;
   uvar->file_pattern = NULL;
   uvar->timeStampsStarting = NULL;
   uvar->timeStampsFinishing = NULL;
   uvar->outSFTnames = NULL;

   /* now register all our user-variable */
   XLALRegisterUvarMember( file_pattern,            STRING, 'i', REQUIRED, "String of SFDB files (possibly from more than one detector, separated by a ;)");
   XLALRegisterUvarMember( timeStampsStarting,      STRING, 's', OPTIONAL, "File(s) containing the starting timestamps of science segments (possibly from more than one detector, separated by a ;)");
   XLALRegisterUvarMember( timeStampsFinishing,     STRING, 'f', OPTIONAL, "File(s) containing the starting timestamps of science segments (possibly from more than one detector, separated by a ;)");
   XLALRegisterUvarMember( outSFTnames,         STRING, 'd', REQUIRED, "Output file(s) (possibly from more than one detector, separated by a ;)");
   XLALRegisterUvarMember(   fmin,               REAL8, 0, REQUIRED, "Lowest frequency to extract from SFTs");
   XLALRegisterUvarMember(   fmax,               REAL8, 0, REQUIRED, "Highest frequency to extract from SFTs");

   return XLAL_SUCCESS;
 
} /* initUserVars() */
