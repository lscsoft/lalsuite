/*
*  Copyright (C) 2014 Reinhard Prix
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

#include <lal/XLALError.h>
#include <lal/LALMalloc.h>
#include <lal/AVFactories.h>
#include <lal/StreamInput.h>
#include <lal/LogPrintf.h>
#include <lal/FileIO.h>
#include <lal/ConfigFile.h>

// some basic consistency checks of the (XLAL) FileIO and Configfile module

int
main(int argc, char *argv[])
{
  XLAL_CHECK ( argc == 1, XLAL_EINVAL, "No input arguments allowed.\n" );
  XLAL_CHECK ( argv != NULL, XLAL_EINVAL );

  char testFile[] = "earth00-19-DE405.dat.gz";

  char *testFilePath;
  XLAL_CHECK ( (testFilePath = XLALFileResolvePathLong ( testFile, TEST_DATA_DIR )) != NULL, XLAL_EINVAL );

 // read gzipped ephemeris file once with XLALCHARReadSequence() and once with XLALFileLoad()
  CHARSequence *sequence = NULL;
  char *data = NULL;
  REAL8 tic;
  // ----- 1. XLALCHARReadSequence() -----
  FILE *fp0;
  tic = XLALGetTimeOfDay();
  XLAL_CHECK ( (fp0 = LALFopen ( testFilePath, "rb" )) != NULL, XLAL_EIO, "Failed to fopen('%s', 'rb')\n", testFilePath );
  XLAL_CHECK ( XLALCHARReadSequence( &sequence, fp0 ) == XLAL_SUCCESS, XLAL_EFUNC );
  LALFclose ( fp0 );
  REAL8 time_CHARReadSequence = XLALGetTimeOfDay() - tic;

  // ----- 2. XLALFileLoad() -----
  tic = XLALGetTimeOfDay();
  XLAL_CHECK ( (data = XLALFileLoad ( testFilePath )) != NULL, XLAL_EFUNC );
  REAL8 time_FileLoad = XLALGetTimeOfDay() - tic;

  // and compare them against each other
  XLAL_CHECK ( strcmp ( data, sequence->data ) == 0, XLAL_EFAILED, "Strings read by XLALCHARReadSequence() and XLALFileLoadData() differ!\n" );

  // time XLALCreateTokenList()
  TokenList *tokens = NULL;
  tic = XLALGetTimeOfDay();
  XLAL_CHECK ( XLALCreateTokenList ( &tokens, data, "\n") == XLAL_SUCCESS, XLAL_EFUNC );
  REAL8 time_CreateTokenList = XLALGetTimeOfDay() - tic;

  // time ParseDatafileContent() [ = 'cleanConfig()' + 'XLALCreateTokenList()' ]
  LALParsedDataFile *content = NULL;
  tic = XLALGetTimeOfDay();
  XLAL_CHECK ( XLALParseDataFileContent ( &content, data ) == XLAL_SUCCESS, XLAL_EFUNC );
  REAL8 time_ParseDataFileContent = XLALGetTimeOfDay() - tic;

  // output results
  XLALPrintInfo ( "'earth00-19-DE405.dat.gz' read and uncompressed by XLALCHARReadSequence() and XLALFileLoadData() is identical!\n" );
  XLALPrintInfo ( "time XLALCHARReadSequence():     %f s\n", time_CHARReadSequence );
  XLALPrintInfo ( "time XLALFileLoad():             %f s\n", time_FileLoad );
  XLALPrintInfo ( "time XLALCreateTokenList():      %f s\n", time_CreateTokenList );
  XLALPrintInfo ( "time XLALParseDataFileContent(): %f s\n", time_ParseDataFileContent );

  // free remaining memory
  XLALFree ( data );
  XLALDestroyCHARVector ( sequence );
  XLALDestroyTokenList ( tokens );
  XLALDestroyParsedDataFile ( content );
  XLALFree ( testFilePath );

  LALCheckMemoryLeaks();

  return XLAL_SUCCESS;

} // main()
