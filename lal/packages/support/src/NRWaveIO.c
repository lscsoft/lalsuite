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
    total mass for setting the timescale.  The output waveform is scaled corresponding
    to a distance of 1Mpc.
*/
REAL4TimeVectorSeries *
XLALReadNRWave( REAL4  mass,       /**< Value of total mass for setting time scale */
		CHAR   *filename   /**< File containing numrel waveform */) 
{

  INT4 r, count;
  FILE *fp=NULL;
  REAL4TimeVectorSeries *out=NULL;
  REAL4VectorSequence *data=NULL;
  REAL4Vector *timeVec=NULL;
  REAL4 massMpc;

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


  /* mass in Mpc -- we multiply h(t) by this factor */
  massMpc = LAL_MRSUN_SI / ( LAL_PC_SI * 1.0e6);

  /* loop over file again */
  count = 0;
  do 
    {
      r = fscanf(fp,"%f%f%f\n", &tmp1, &tmp2, &tmp3);
      /* make sure the line has the right number of entries or is EOF */
      if (r==3) {
	timeVec->data[count] = tmp1;
	data->data[count] = massMpc * tmp2;
	data->data[data->vectorLength + count] = massMpc * tmp3;
	count++;
      }
      
    } while ( r != EOF);
  

  /* scale time */
  out->deltaT = LAL_MTSUN_SI * mass * ( timeVec->data[1] - timeVec->data[0]);
  /* need additional consistency check on timeVec to 
     make sure it is spaced uniformly */

  fclose(fp);

  XLALDestroyREAL4Vector (timeVec);

  out->data = data;
  
  /* normal exit */	
  return out;
  
} /* XLALReadNRWave() */



/** Function for reading a numerical relativity metadata file.
    It returns a list of numrel wave parameters.  It uses 
    LALParseDataFile() for reading the data file.  This automatically
    takes care of removing comment lines starting with # and other details. 
 */
void
LALNRDataFind( LALStatus *status,
	       NRWaveCatalog *out, /**< list of numrel metadata */
	       CHAR   *filename /**< File with metadata information */)
{
  LALParsedDataFile *cfgdata=NULL;
  UINT4 k, numWaves;

  INITSTATUS (status, "LALNRDataFind", NRWAVEIOC);
  ATTATCHSTATUSPTR (status); 

  ASSERT (filename != NULL, status, NRWAVEIO_ENULL, NRWAVEIO_MSGENULL );
  ASSERT ( out != NULL, status, NRWAVEIO_ENULL, NRWAVEIO_MSGENULL );

  TRY( LALParseDataFile ( status->statusPtr, &cfgdata, filename), status);
  numWaves = cfgdata->lines->nTokens; /*number of waves */

  /* allocate memory for output catalog */
  out->length = numWaves;
  out->data = LALCalloc(1, out->length * sizeof(NRWaveMetaData));

  /* now get wave parameters from each line of data */
  for (k = 0; k < numWaves; k++) {
    TRY(LALGetSingleNRMetaData( status->statusPtr, out->data + k, cfgdata->lines->tokens[k]), status);
  }

  TRY( LALDestroyParsedDataFile (status->statusPtr, &cfgdata), status);

  DETATCHSTATUSPTR(status);
  RETURN(status);
}



/** Parse a single string to fill the NRWaveMetaData structure */
void
LALGetSingleNRMetaData( LALStatus             *status,
			NRWaveMetaData  *out, /**< Meta data struct to be filled */
			const CHAR            *cfgstr /**< config string containing the data for a single NR wave*/) 
{
  REAL4 tmpR[7];
  INT4  tmpI[2];
  INT4 test;

  INITSTATUS (status, "LALGetSingleNRMetaData", NRWAVEIOC);
  ATTATCHSTATUSPTR (status); 

  ASSERT (cfgstr != NULL, status, NRWAVEIO_ENULL, NRWAVEIO_MSGENULL );
  ASSERT (out != NULL, status, NRWAVEIO_ENULL, NRWAVEIO_MSGENULL );

  test = sscanf(cfgstr, "%f%f%f%f%f%f%f%d%d%s", tmpR, tmpR+1, tmpR+2, tmpR+3, tmpR+4, tmpR+5, 
	 tmpR+6, tmpI, tmpI+1, out->filename);

  /* Check the metadata file format */
  if ( test != 10) {
    /* there must be exactly 10 data entries 
       -- massratio, spin1[3], spin2[3], l, m, filename */
    ABORT( status, NRWAVEIO_EFORMAT, NRWAVEIO_MSGEFORMAT );
  }

  /*** need a few more format checks here ***/

  /* copy values to out */
  out->massRatio = tmpR[0];
  out->spin1[0] = tmpR[1];
  out->spin1[1] = tmpR[2];
  out->spin1[2] = tmpR[3];
  out->spin2[0] = tmpR[4];
  out->spin2[1] = tmpR[5];
  out->spin2[2] = tmpR[6];

  out->mode[0] = tmpI[0];
  out->mode[1] = tmpI[1];

  DETATCHSTATUSPTR(status);
  RETURN(status);
}
