/*
 * Copyright (C) 2006 S.Fairhurst, B. Krishnan, L.Santamaria
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

/** \file NRWaveIO.c
 *  \ingroup NRWaveIO
 *  \author S.Fairhurst, B.Krishnan, L.Santamaria
 * 
 *  \brief Functions for reading/writing numerical relativity waveforms
 *
 * $Id$ 
 *
 */

#include <lal/LALStdio.h>
#include <lal/FileIO.h>
#include <lal/NRWaveIO.h>



NRCSID( NRWAVEIOC, "$Id$");


/** Reads a numerical relativity waveform given a filename and a value of the 
    total mass for setting the timescale */
REAL4TimeVectorSeries *
XLALReadNRWave( REAL8  mass,       /**< Value of total mass for setting time scale */
		CHAR   *filename   /**< File containing numrel waveform */) 
{

  INT4 r, count;
  FILE *fp=NULL;
  REAL4TimeVectorSeries *out=NULL;
  REAL4VectorSequence *data=NULL;
  REAL4Vector *timeVec=NULL;

  CHAR str1[16], str2[16],str3[16], str4[16];
  REAL4 tmp1, tmp2, tmp3;


  fp = fopen(filename, "r");
  if (!fp) {
    XLAL_ERROR_NULL( "XLALReadNRWave", XLAL_ENAME );
  }

  /* read the first comment line which is expected to be  
     something like  "# time h+ hx"  -- do we need this?*/ 
  r = fscanf(fp, "%s%s%s%s\n", str1, str2, str3, str4);
  
  /* count number of lines */
  count = 0;
  do 
    {
      r = fscanf(fp,"%f%f%f\n", &tmp1, &tmp2, &tmp3);
      /* make sure the line has the right number of entries or is EOF */
      if (r==3) count++;
    } while ( r != EOF);
  rewind(fp);  

  /* read the first line again */  
  r = fscanf(fp, "%s%s%s%s\n", str1, str2, str3, str4);

  /* allocate memory */
  out = LALCalloc(1, sizeof(*out));
  if (!out) {
    XLAL_ERROR_NULL( "XLALReadNRWave", XLAL_ENOMEM );
  }
  strcpy(out->name,filename);
  out->f0 = 0;
  
  data =  XLALCreateREAL4VectorSequence (2, count);
  if (!data) {
    XLAL_ERROR_NULL( "XLALReadNRWave", XLAL_ENOMEM );
  }
  
  timeVec = XLALCreateREAL4Vector ( count );
  if (!timeVec) {
    XLAL_ERROR_NULL( "XLALReadNRWave", XLAL_ENOMEM );
  }

  /* loop over file again */
  count = 0;
  do 
    {
      r = fscanf(fp,"%f%f%f\n", &tmp1, &tmp2, &tmp3);
      /* make sure the line has the right number of entries or is EOF */
      if (r==3) {
	data->data[count] = tmp2;
	data->data[data->length + count] = tmp3;
	count++;
      }

    } while ( r != EOF);


  /* need additional consistency check on timeVec*/
  out->deltaT = LAL_MTSUN_SI * mass * ( timeVec->data[1] - timeVec->data[0]);

  fclose(fp);

  XLALDestroyREAL4Vector (timeVec);

  out->data = data;
  
  /* normal exit */	
  return out;
  
} /* XLALReadNRWave() */


