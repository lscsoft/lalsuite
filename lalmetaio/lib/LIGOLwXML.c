/*
*  Copyright (C) 2007 Andres C. Rodriguez, Sukanta Bose, Alexander Dietz, Duncan Brown, Jolien Creighton, Kipp Cannon, Lisa M. Goggin, Patrick Brady, Robert Adam Mercer, Saikat Ray-Majumder, Anand Sengupta, Stephen Fairhurst, Xavier Siemens, Sean Seader, Thomas Cokelaer
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

/*-----------------------------------------------------------------------
 *
 * File Name: LIGOLwXML.c
 *
 * Author: Brown, D. A.
 *
 *-----------------------------------------------------------------------
 */

/**
 * \author Brown, D. A.
 * \file
 * \ingroup lalmetaio_general
 *
 * \brief Routines to write LIGO metadata database structures to LIGO lightweight XML files.
 *
 * ### Description ###
 *
 * The routine \c XLALCreateLIGOLwXMLFileName creates a name for a  LIGO lightweight XML file that is in accordance with the specifications of document T050017.
 *
 * ### Algorithm ###
 *
 * None.
 *
 * ### Uses ###
 *
 * <tt>fopen()</tt>
 * <tt>XLALFilePrintf()</tt>
 * <tt>fclose()</tt>
 *
 * ### Notes ###
 *
 */


#include <lal/FileIO.h>
#include <lal/LALMalloc.h>
#include <lal/LIGOLwXML.h>
#include <lal/XLALError.h>
#include <LIGOLwXMLHeaders.h>


/**
 * Open an XML file for writing.  The return value is a pointer to a new
 * LIGOLwXMLStream file handle or NULL on failure.
 */
LIGOLwXMLStream *
XLALOpenLIGOLwXMLFile (
    const char *path
)
{
  LIGOLwXMLStream *new;

  /* malloc a new XML file handle */

  new = XLALMalloc( sizeof( *new ) );
  if ( ! new )
    XLAL_ERROR_NULL( XLAL_EFUNC );

  /* fopen() the underlying C file */

  new->fp = XLALFileOpen( path, "w" );
  if ( ! new->fp )
  {
    XLALFree(new);
    XLAL_ERROR_NULL( XLAL_EFUNC );
  }

  /* write the XML header */

  if ( XLALFilePuts( LAL_LIGOLW_XML_HEADER, new->fp ) < 0 )
  {
    XLALFileClose( new->fp );
    XLALFree( new );
    XLAL_ERROR_NULL( XLAL_EIO );
  }

  /* done */

  return new;
}


/**
 * Close an XML stream.  On failure the stream is left in an undefined
 * state, and has not been free()'ed.  Sorry.
 */

int
XLALCloseLIGOLwXMLFile (
  LIGOLwXMLStream *xml
)
{
  if ( xml )
  {
    if ( XLALFilePuts( LAL_LIGOLW_XML_FOOTER, xml->fp ) < 0 )
      /* can't write XML footer */
      XLAL_ERROR( XLAL_EIO );
    if ( XLALFileClose( xml->fp ) )
      /* fclose() on the underlying C file failed */
      XLAL_ERROR( XLAL_EFUNC );
  }

  XLALFree( xml );

  return 0;
}

/**
 * Creates a XML filename accordingly to document T050017
 */

int XLALCreateLIGODataFileName(
        char* filename,
        size_t size,
        const char* dataSource,
        const char* dataDescription,
        const LIGOTimeGPS* gpsStartTime,
        const LIGOTimeGPS* gpsEndTime,
        const char* extension
)
{
     INT4 gpsDuration;

     /* check input structures */
     if (!filename || !dataSource || !dataDescription ||
	 !gpsStartTime || !gpsEndTime || !extension)
          XLAL_ERROR(XLAL_EFAULT);

     /* check the correctnes of the input strings */
     if ( strchr(dataSource, '-') || strchr(dataDescription, '-'))
     {
          filename = NULL;
          XLALPrintError("the input character strings contain invalid"
			 " dashes ('-').");
          XLAL_ERROR(XLAL_EINVAL);
      }

      /* calculate the GPS duration */
      gpsDuration = gpsEndTime->gpsSeconds - gpsStartTime->gpsSeconds;
      if (gpsEndTime->gpsNanoSeconds > 0) ++gpsDuration;

      /* and here put it all together */
      snprintf( filename, size, "%s-%s-%d-%d.%s",
		   dataSource, dataDescription, gpsStartTime->gpsSeconds,
		   gpsDuration, extension );

      /* return success */
      return 0;
}
