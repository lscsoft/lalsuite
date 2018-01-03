/*
*  Copyright (C) 2014 Evan Goetz
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

#include <lal/UserInput.h>
#include "templates.h"

typedef struct
{
   REAL8 Pmin;
   REAL8 Pmax;
   REAL8 dfmin;
   REAL8 dfmax;
   REAL8 Tsft;
   REAL8 SFToverlap;
   REAL8 Tobs;
   INT4 maxVectorLength;
   INT4 minTemplateLength;
   INT4 maxTemplateLength;
   INT4 vectorMath;
   BOOLEAN exactflag;
   CHAR *filename;
} UserVariables_t;

INT4 InitUserVars2(UserVariables_t *uvar, int argc, char *argv[]);

INT4 main(int argc, char *argv[])
{
   UserVariables_t XLAL_INIT_DECL(uvar);
   XLAL_CHECK ( InitUserVars2(&uvar, argc, argv) == XLAL_SUCCESS, XLAL_EFUNC );

   //Set vectormath
   //if (uvar.vectorMath==1) XLAL_CHECK( XLALVectorDeviceSet(VECTORDEVICE_SSE) == XLAL_SUCCESS, XLAL_EFUNC );
   //else if (uvar.vectorMath==2) XLAL_CHECK( XLALVectorDeviceSet(VECTORDEVICE_AVX) == XLAL_SUCCESS, XLAL_EFUNC );
   //else XLAL_CHECK( XLALVectorDeviceSet(VECTORDEVICE_FPU) == XLAL_SUCCESS, XLAL_EFUNC );

   TwoSpectTemplateVector *vector = NULL;
   XLAL_CHECK( (vector = generateTwoSpectTemplateVector(uvar.Pmin, uvar.Pmax, uvar.dfmin, uvar.dfmax, uvar.Tsft, uvar.SFToverlap, uvar.Tobs, uvar.maxVectorLength, uvar.minTemplateLength, uvar.maxTemplateLength, uvar.vectorMath, uvar.exactflag)) != NULL, XLAL_EFUNC );

   if (XLALUserVarWasSet(&uvar.filename)) XLAL_CHECK( writeTwoSpectTemplateVector(vector, uvar.filename) == XLAL_SUCCESS, XLAL_EFUNC );

   destroyTwoSpectTemplateVector(vector);
   XLALDestroyUserVars();
   
   return 0;
}

INT4 InitUserVars2(UserVariables_t *uvar, int argc, char *argv[])
{
   XLAL_CHECK ( uvar != NULL, XLAL_EINVAL, "Invalid NULL input 'uvar'\n");
   XLAL_CHECK ( argv != NULL, XLAL_EINVAL, "Invalid NULL input 'argv'\n");

   uvar->Tsft = 1800;
   uvar->SFToverlap = 900;

   XLALRegisterUvarMember(  Pmin,              REAL8, 0 , REQUIRED, "Minimum period");
   XLALRegisterUvarMember(  Pmax,              REAL8, 0 , REQUIRED, "Maximum period");
   XLALRegisterUvarMember(  dfmin,             REAL8, 0 , REQUIRED, "Minimum modulation depth");
   XLALRegisterUvarMember(  dfmax,             REAL8, 0 , REQUIRED, "Maximum modulation depth");
   XLALRegisterUvarMember(  Tsft,              REAL8, 0 , OPTIONAL, "SFT coherence length");
   XLALRegisterUvarMember(  SFToverlap,        REAL8, 0 , OPTIONAL, "SFT overlap in second");
   XLALRegisterUvarMember(  Tobs,              REAL8, 0 , REQUIRED, "Total observation time");
   XLALRegisterUvarMember(   minTemplateLength, INT4, 0 , REQUIRED, "Minimum number of pixels in templates");
   XLALRegisterUvarMember(   maxTemplateLength, INT4, 0 , REQUIRED, "Maximum number of pixels in tempaltes");
   XLALRegisterUvarMember(   maxVectorLength,   INT4, 0 , REQUIRED, "Maximum vector length");
   XLALRegisterUvarMember(   vectorMath,        INT4, 0 , OPTIONAL, "Vector math flag: 0 = no SSE/AVX, 1 = SSE, 2 = AVX");
   XLALRegisterUvarMember(  exactflag,         BOOLEAN, 0 , OPTIONAL, "Flag to specify using exact templates");
   XLALRegisterUvarMember(filename,          STRING, 0 , OPTIONAL, "Filename of output file (if not specified, the vector is destroyed upon exit)");

   BOOLEAN should_exit = 0;
   XLAL_CHECK( XLALUserVarReadAllInput( &should_exit, argc, argv ) == XLAL_SUCCESS, XLAL_EFUNC );
   if ( should_exit ) exit (1);

   return XLAL_SUCCESS;
}
