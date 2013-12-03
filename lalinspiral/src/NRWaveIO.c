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

/**
 * \file NRWaveIO.c
 * \ingroup NRWaveIO_h
 * \author S.Fairhurst, B.Krishnan, L.Santamaria
 *
 * \brief Functions for reading/writing numerical relativity waveforms
 *
 */

#include <lal/LALStdio.h>
#include <lal/FileIO.h>
#include <lal/NRWaveIO.h>
#include <lal/NRWaveInject.h>

/**
 * Functionfor reading the numrel waveform -- just returns the numrel
 * data as it is without any rescaling of time or amplitude
 */

void LALReadNRWave_raw(LALStatus *status,	/**< pointer to LALStatus structure */
		       REAL4TimeVectorSeries **out, /**< [out] output time series for h+ and hx */
		       const CHAR  *filename        /**< [in] File containing numrel waveform */)
{

  UINT4 length, k, r;
  REAL4TimeVectorSeries *ret=NULL;
  REAL4VectorSequence *data=NULL;
  REAL4Vector *timeVec=NULL;
  LALParsedDataFile *cfgdata=NULL;
  REAL4 tmp1, tmp2, tmp3;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /* some consistency checks */
  ASSERT (filename != NULL, status, NRWAVEIO_ENULL, NRWAVEIO_MSGENULL );
  ASSERT ( out != NULL, status, NRWAVEIO_ENULL, NRWAVEIO_MSGENULL );
  ASSERT ( *out == NULL, status, NRWAVEIO_ENONULL, NRWAVEIO_MSGENONULL );

  TRY( LALParseDataFile ( status->statusPtr, &cfgdata, filename), status);
  length = cfgdata->lines->nTokens; /*number of data points */


  /* allocate memory */
  ret = LALCalloc(1, sizeof(*ret));
  if (!ret) {
    ABORT( status, NRWAVEIO_ENOMEM, NRWAVEIO_MSGENOMEM );
  }
  strcpy(ret->name,filename);
  ret->f0 = 0;

  data =  XLALCreateREAL4VectorSequence (2, length);
  if (!data) {
    ABORT( status, NRWAVEIO_ENOMEM, NRWAVEIO_MSGENOMEM );
  }

  timeVec = XLALCreateREAL4Vector (length);
  if (!timeVec) {
    ABORT( status, NRWAVEIO_ENOMEM, NRWAVEIO_MSGENOMEM );
  }

  /* now get the data */
  for (k = 0; k < length; k++) {
    r = sscanf(cfgdata->lines->tokens[k], "%f%f%f", &tmp1, &tmp2, &tmp3);

    /* Check the data file format */
    if ( r != 3) {
      /* there must be exactly 3 data entries -- time, h+, hx */
      ABORT( status, NRWAVEIO_EFORMAT, NRWAVEIO_MSGEFORMAT );
    }

    timeVec->data[k] = tmp1;
    data->data[k] = tmp2;
    data->data[data->vectorLength + k] = tmp3;

  }

  /*  scale time */
  ret->deltaT = timeVec->data[1] - timeVec->data[0];

  /* might also want to go through timeVec to make sure it is evenly spaced */


  ret->data = data;
  (*out) = ret;

  XLALDestroyREAL4Vector (timeVec);
  TRY( LALDestroyParsedDataFile (status->statusPtr, &cfgdata), status);

  DETATCHSTATUSPTR(status);
  RETURN(status);

} /* LALReadNRWave() */



void LALReadNRWave_raw_real8(LALStatus *status,	/**< pointer to LALStatus structure */
		       REAL8TimeVectorSeries **out, /**< [out] output time series for h+ and hx */
		       const CHAR  *filename        /**< [in] File containing numrel waveform */)
{

  UINT4 length, k, r;
  REAL8TimeVectorSeries *ret=NULL;
  REAL8VectorSequence *data=NULL;
  REAL8Vector *timeVec=NULL;
  LALParsedDataFile *cfgdata=NULL;
  REAL8 tmp1, tmp2, tmp3;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /* some consistency checks */
  ASSERT (filename != NULL, status, NRWAVEIO_ENULL, NRWAVEIO_MSGENULL );
  ASSERT ( out != NULL, status, NRWAVEIO_ENULL, NRWAVEIO_MSGENULL );
  ASSERT ( *out == NULL, status, NRWAVEIO_ENONULL, NRWAVEIO_MSGENONULL );

  TRY( LALParseDataFile ( status->statusPtr, &cfgdata, filename), status);
  length = cfgdata->lines->nTokens; /*number of data points */


  /* allocate memory */
  ret = LALCalloc(1, sizeof(*ret));
  if (!ret) {
    ABORT( status, NRWAVEIO_ENOMEM, NRWAVEIO_MSGENOMEM );
  }
  strcpy(ret->name,filename);
  ret->f0 = 0;

  data =  XLALCreateREAL8VectorSequence (2, length);
  if (!data) {
    ABORT( status, NRWAVEIO_ENOMEM, NRWAVEIO_MSGENOMEM );
  }

  timeVec = XLALCreateREAL8Vector (length);
  if (!timeVec) {
    ABORT( status, NRWAVEIO_ENOMEM, NRWAVEIO_MSGENOMEM );
  }

  /* now get the data */
  for (k = 0; k < length; k++) {
    r = sscanf(cfgdata->lines->tokens[k], "%lf%lf%lf", &tmp1, &tmp2, &tmp3);

    /* Check the data file format */
    if ( r != 3) {
      /* there must be exactly 3 data entries -- time, h+, hx */
      ABORT( status, NRWAVEIO_EFORMAT, NRWAVEIO_MSGEFORMAT );
    }

    timeVec->data[k] = tmp1;
    data->data[k] = tmp2;
    data->data[data->vectorLength + k] = tmp3;

  }

  /*  scale time */
  ret->deltaT = timeVec->data[1] - timeVec->data[0];

  /* might also want to go through timeVec to make sure it is evenly spaced */


  ret->data = data;
  (*out) = ret;

  XLALDestroyREAL8Vector (timeVec);
  TRY( LALDestroyParsedDataFile (status->statusPtr, &cfgdata), status);

  DETATCHSTATUSPTR(status);
  RETURN(status);

} /* LALReadNRWave() */


/**
 * Reads a numerical relativity waveform given a filename and a value of the
 * total mass for setting the timescale.  The output waveform is scaled corresponding
 * to a distance of 1Mpc.
 */
void LALReadNRWave(LALStatus *status,		/**< pointer to LALStatus structure */
		   REAL4TimeVectorSeries **out, /**< [out] output time series for h+ and hx */
		   const REAL4  mass,           /**< [in] Value of total mass for setting time scale */
		   const CHAR  *filename        /**< [in] File containing numrel waveform */)
{

  UINT4 length, k, r;
  REAL4TimeVectorSeries *ret=NULL;
  REAL4VectorSequence *data=NULL;
  REAL4Vector *timeVec=NULL;
  REAL4 massMpc;
  LALParsedDataFile *cfgdata=NULL;
  REAL4 tmp1, tmp2, tmp3;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /* some consistency checks */
  ASSERT (filename != NULL, status, NRWAVEIO_ENULL, NRWAVEIO_MSGENULL );
  ASSERT ( mass > 0, status, NRWAVEIO_EVAL, NRWAVEIO_MSGEVAL );
  ASSERT ( out != NULL, status, NRWAVEIO_ENULL, NRWAVEIO_MSGENULL );
  ASSERT ( *out == NULL, status, NRWAVEIO_ENONULL, NRWAVEIO_MSGENONULL );

  TRY( LALParseDataFile ( status->statusPtr, &cfgdata, filename), status);
  length = cfgdata->lines->nTokens; /*number of data points */


  /* allocate memory */
  ret = LALCalloc(1, sizeof(*ret));
  if (!ret) {
    ABORT( status, NRWAVEIO_ENOMEM, NRWAVEIO_MSGENOMEM );
  }
  strcpy(ret->name,filename);
  ret->f0 = 0;

  data =  XLALCreateREAL4VectorSequence (2, length);
  if (!data) {
    ABORT( status, NRWAVEIO_ENOMEM, NRWAVEIO_MSGENOMEM );
  }

  timeVec = XLALCreateREAL4Vector (length);
  if (!timeVec) {
    ABORT( status, NRWAVEIO_ENOMEM, NRWAVEIO_MSGENOMEM );
  }

  /* mass in Mpc -- we multiply h(t) by this factor */
  massMpc = mass * LAL_MRSUN_SI / ( LAL_PC_SI * 1.0e6);

  /* now get the data */
  for (k = 0; k < length; k++) {
    r = sscanf(cfgdata->lines->tokens[k], "%f%f%f", &tmp1, &tmp2, &tmp3);

    /* Check the data file format */
    if ( r != 3) {
      /* there must be exactly 3 data entries -- time, h+, hx */
      ABORT( status, NRWAVEIO_EFORMAT, NRWAVEIO_MSGEFORMAT );
    }

    timeVec->data[k] = tmp1;
    data->data[k] = massMpc * tmp2;
    data->data[data->vectorLength + k] = massMpc * tmp3;

  }

  /*  scale time */
  ret->deltaT = LAL_MTSUN_SI * mass * ( timeVec->data[1] - timeVec->data[0]);

  /* go through timeVec to make sure it is uniformly spaced */


  ret->data = data;
  (*out) = ret;

  XLALDestroyREAL4Vector (timeVec);
  TRY( LALDestroyParsedDataFile (status->statusPtr, &cfgdata), status);

  DETATCHSTATUSPTR(status);
  RETURN(status);

} /* LALReadNRWave() */



/**
 * Function for reading a numerical relativity metadata file.
 * It returns a list of numrel wave parameters.  It uses
 * LALParseDataFile() for reading the data file.  This automatically
 * takes care of removing comment lines starting with # and other details.
 */
void
LALNRDataFind( LALStatus *status,   /**< pointer to LALStatus structure */
	       NRWaveCatalog *out,  /**< [out] list of numrel metadata */
	       const CHAR *dir,     /**< [in] directory with data files */
	       const CHAR *filename /**< [in] File with metadata information */)
{
  LALParsedDataFile *cfgdata=NULL;
  UINT4 k, numWaves;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  ASSERT (filename != NULL, status, NRWAVEIO_ENULL, NRWAVEIO_MSGENULL );
  ASSERT ( out != NULL, status, NRWAVEIO_ENULL, NRWAVEIO_MSGENULL );
  ASSERT ( dir != NULL, status, NRWAVEIO_ENULL, NRWAVEIO_MSGENULL );

  TRY( LALParseDataFile ( status->statusPtr, &cfgdata, filename), status);
  numWaves = cfgdata->lines->nTokens; /*number of waves */

  /* allocate memory for output catalog */
  out->length = numWaves;
  out->data = LALCalloc(1, out->length * sizeof(NRWaveMetaData));
  if ( out->data == NULL) {
    ABORT( status, NRWAVEIO_ENOMEM, NRWAVEIO_MSGENOMEM );
  }

  /* now get wave parameters from each line of data */
  for (k = 0; k < numWaves; k++) {
    TRY(LALGetSingleNRMetaData( status->statusPtr, out->data + k, dir, cfgdata->lines->tokens[k]), status);
  }

  TRY( LALDestroyParsedDataFile (status->statusPtr, &cfgdata), status);

  DETATCHSTATUSPTR(status);
  RETURN(status);
}



/** Parse a single string to fill the NRWaveMetaData structure */
void
LALGetSingleNRMetaData( LALStatus       *status, /**< pointer to LALStatus structure */
			NRWaveMetaData  *out,   /**< [out] Meta data struct to be filled */
			const CHAR      *dir,   /**< [in] data directory */
			const CHAR      *cfgstr /**< [in] config string containing the data for a single NR wave*/)
{
  REAL4 tmpR[7];
  INT4  tmpI[2];
  INT4 test;
  CHAR tmpStr[512];

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  ASSERT (cfgstr != NULL, status, NRWAVEIO_ENULL, NRWAVEIO_MSGENULL );
  ASSERT (out != NULL, status, NRWAVEIO_ENULL, NRWAVEIO_MSGENULL );
  ASSERT (dir != NULL, status, NRWAVEIO_ENULL, NRWAVEIO_MSGENULL );

  test = sscanf(cfgstr, "%f%f%f%f%f%f%f%d%d%s", tmpR, tmpR+1, tmpR+2, tmpR+3, tmpR+4, tmpR+5,
	 tmpR+6, tmpI, tmpI+1, tmpStr);

  /* Check the metadata file format */
  if ( test != 10) {
    /* there must be exactly 10 data entries
       -- massratio, spin1[3], spin2[3], l, m, filename */
    ABORT( status, NRWAVEIO_EFORMAT, NRWAVEIO_MSGEFORMAT );
  }

  /* the mass ratio must be positive */
  ASSERT (tmpR[0] >= 0, status, NRWAVEIO_EVAL, NRWAVEIO_MSGEVAL );

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

  strcpy(out->filename, dir);
  strcat(out->filename, "/");
  strcat(out->filename, tmpStr);

  DETATCHSTATUSPTR(status);
  RETURN(status);
}

/** Put the main functionalities of nr_wave.c together */
void
LALAddStrainModes(
  LALStatus              *status,	/**< pointer to LALStatus structure */
  REAL4TimeVectorSeries  **outStrain,   /**< h+, hx data       */
  NRWaveCatalog          *nrCatalog,    /**< NR wave metadata struct        */
  INT4                   modeLlo,      /**< contains modeLlo and modeLhi   */
  INT4                   modeLhi,      /**< modeLhi                        */
  const SimInspiralTable       *thisInj     /**< injection                      */)
{
  INT4 modeL, modeM;
  REAL4TimeVectorSeries *sumStrain=NULL;
  REAL4TimeVectorSeries *tempStrain=NULL;
  NRWaveMetaData thisMetaData;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /* loop over l values */
  for ( modeL = modeLlo; modeL <= modeLhi; modeL++ )
    {
      /* loop over m values */
      for ( modeM = -modeL; modeM <= modeL; modeM++ )
        {
          /* find nearest matching numrel waveform */
          XLALFindNRFile( &thisMetaData, nrCatalog, thisInj, modeL, modeM );

          /* read numrel waveform */
          TRY( LALReadNRWave( status->statusPtr, &tempStrain, thisInj->mass1 +
                thisInj->mass2, thisMetaData.filename ), status );

          /* compute the h+ and hx for given inclination and coalescence phase*/
          XLALOrientNRWave( tempStrain, thisMetaData.mode[0],
			    thisMetaData.mode[1], thisInj->inclination, thisInj->coa_phase );

	  if (sumStrain == NULL) {

	    sumStrain = LALCalloc(1, sizeof(*sumStrain));

	    sumStrain->data =  XLALCreateREAL4VectorSequence(2, tempStrain->data->vectorLength);
	    sumStrain->deltaT = tempStrain->deltaT;
	    sumStrain->f0 = tempStrain->f0;
	    sumStrain->sampleUnits = tempStrain->sampleUnits;

	    memset(sumStrain->data->data,0.0,2*tempStrain->data->vectorLength*sizeof(REAL4));

	    sumStrain = XLALSumStrain( sumStrain, tempStrain );
	  }

	  else {
	    sumStrain = XLALSumStrain( sumStrain, tempStrain );
	  }

	  /*fprintf(stdout, "\nInjecting NR waveform from file %s", thisMetaData.filename);*/

          /* clear memory for strain */
          XLALDestroyREAL4VectorSequence ( tempStrain->data );
          LALFree( tempStrain );
          tempStrain = NULL;

        } /* end loop over modeM values */

      } /* end loop over modeL values */

  /*fprintf(stdout, "\nNR injections done\n");*/

  *outStrain = sumStrain;

  DETATCHSTATUSPTR(status);
  RETURN(status);

}




/** Main driver funtion for doing Numerical Relativity Injections */
void LALDriveNRInject( LALStatus *status,	/**< pointer to LALStatus structure */
		       REAL4TimeSeries *injData, /**< The time series to inject into */
		       SimInspiralTable *injections, /**< The list of injections to perform */
		       NumRelInjectParams *params /**< Parameters */
		       )
{

  REAL4TimeVectorSeries *sumStrain = NULL;
  SimInspiralTable *thisInj    = NULL;   /* current injection */

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /* loop over injections */
  for ( thisInj = injections; thisInj; thisInj = thisInj->next )
    {

      TRY( LALAddStrainModes(status->statusPtr, &sumStrain, params->nrCatalog,
				  params->modeLlo, params->modeLhi, thisInj), status);

      TRY( LALInjectStrainGW( status->statusPtr, injData, sumStrain, thisInj, params->ifo, params->dynRange), status);

      XLALDestroyREAL4VectorSequence ( sumStrain->data );
      LALFree( sumStrain );
      sumStrain = NULL;

    } /* end loop over injections */


  DETATCHSTATUSPTR(status);
  RETURN(status);

}


