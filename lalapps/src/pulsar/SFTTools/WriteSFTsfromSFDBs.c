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
    UINT4 numifo;
    REAL8 fmin;
    REAL8 fmax;
    CHAR *file_pattern;
    CHAR *timeStampsStarting;
    CHAR *timeStampsFinishing;
    CHAR *outputSingleSFT;
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
    inputSFTs =  XLALReadSFDB(uvar.fmin, uvar.fmax, uvar.file_pattern, uvar.timeStampsStarting, uvar.timeStampsFinishing);


    for ( UINT4 k = 0; k < uvar.numifo; k++) {

        FILE *fpSingleSFT = NULL;
        CHAR   fileSigma[ 500 ];
        strcpy ( fileSigma, uvar.outputSingleSFT);
        if (uvar.numifo>1) {
            if (k==0) strcat ( fileSigma, "_H1");  // Assuming that the returned SFTs are alphabetically sorted in the first index
            if (k==1) strcat ( fileSigma, "_L1");
            if (k==2) strcat ( fileSigma, "_V1");
        }
        XLAL_CHECK ( ( fpSingleSFT = fopen ( fileSigma, "wb" )) != NULL, XLAL_EIO, "Failed to open singleSFT file '%s' for writing\n", uvar.outputSingleSFT );

        UINT4 numSFTs = inputSFTs->data[k]->length;
        for ( UINT4 i=0; i < numSFTs; i++ )
        {   
            const SFTtype *thisSFT = inputSFTs->data[k]->data + i;

            /*const CHAR *sft_comment;
            CHAR *new_comment;
            UINT4 comment_len = 0;
            comment_len = strlen ( add_comment ) + 10;
            sft_comment = oneSFTCatalog.data->comment;
            if ( sft_comment ) comment_len += strlen ( sft_comment );
            XLAL_CHECK_MAIN ( ( new_comment  = LALCalloc (1, comment_len )) != NULL, XLAL_ENOMEM );
            if ( sft_comment ) {
              strcpy ( new_comment, sft_comment );
              strcat ( new_comment, ";\n");
            }
            strcat ( new_comment, add_comment );*/
        

        
            // if user asked for single-SFT output, add this SFT to the open file
            //if ( uvar.outputSingleSFT ) XLAL_CHECK ( XLAL_SUCCESS == XLALWriteSFT2fp( &(thisSFT->data[0]), fpSingleSFT, new_comment ), XLAL_EFUNC,  "XLALWriteSFT2fp() failed to write SFT to '%s'!\n", uvar.outputSingleSFT );
            XLAL_CHECK ( XLAL_SUCCESS == XLALWriteSFT2fp( thisSFT, fpSingleSFT, NULL ), XLAL_EFUNC,  "XLALWriteSFT2fp() failed to write SFT to '%s'!\n", uvar.outputSingleSFT );
                
            //XLALFree ( new_comment );
 
        } /* for i < numSFTs */

        if ( fpSingleSFT ) fclose ( fpSingleSFT );

    }

    XLALDestroyMultiSFTVector (inputSFTs); 

    return 0;

}


int initUserVars ( UserInput_t *uvar ) {
   
   XLAL_CHECK ( uvar != NULL, XLAL_EINVAL );

   uvar->numifo = 0;
   uvar->fmin = 0;
   uvar->fmax = 0;
   uvar->file_pattern = NULL;
   uvar->timeStampsStarting = NULL;
   uvar->timeStampsFinishing = NULL;
   uvar->outputSingleSFT = NULL;

   /* now register all our user-variable */
   XLALRegisterUvarMember( file_pattern,            STRING, 'i', REQUIRED, "File-pattern for input SFDBs");
   XLALRegisterUvarMember( timeStampsStarting,      STRING, 's', OPTIONAL, "File-pattern for timestamp file with starting times of SCIENCE segments");
   XLALRegisterUvarMember( timeStampsFinishing,     STRING, 'f', OPTIONAL, "File-pattern for timestamp file with ending times of SCIENCE segments");
   XLALRegisterUvarMember( outputSingleSFT,         STRING, 'd', REQUIRED, "File-pattern for output file");
   XLALRegisterUvarMember(   fmin,               REAL8, 0, REQUIRED, "Lowest frequency to extract from SFTs");
   XLALRegisterUvarMember(   fmax,               REAL8, 0, REQUIRED, "Highest frequency to extract from SFTs");
   XLALRegisterUvarMember(   numifo,             UINT4, 0,  REQUIRED, "Number of detectors");

   return XLAL_SUCCESS;
 
} /* initUserVars() */
