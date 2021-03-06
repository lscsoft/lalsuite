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
 *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *  MA  02110-1301  USA
 */

#ifndef _NRWAVEIO_H  	/* Double-include protection. */
#define _NRWAVEIO_H

/* includes */
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/ConfigFile.h>
#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/LALStdio.h>
#include <lal/FileIO.h>
#include <lal/StreamInput.h>
#include <lal/LIGOMetadataTables.h>
/*  #include <lal/NRWaveInject.h> */

#ifdef  __cplusplus   /* C++ protection. */
extern "C" {
#endif

/**
 * \defgroup NRWaveIO_h Header NRWaveIO.h
 * \ingroup lalinspiral_UNCLASSIFIED
 * \author S.Fairhurst, B. Krishnan, L.Santamaria
 *
 * \brief Module for reading/writing Numrel waveforms
 *
 */
/** @{ */

/** \name Error codes */
/** @{ */
#define NRWAVEIO_ENULL 	  1	/**< Null pointer */
#define NRWAVEIO_EFILE 	  2	/**< Error in file-IO */
#define NRWAVEIO_ENONULL  3	/**< Not a Null pointer */
#define NRWAVEIO_ENOMEM   4	/**< Memory ellocation error */
#define NRWAVEIO_EVAL 	  5	/**< Invalid value */
#define NRWAVEIO_EFORMAT  6	/**< Meta data file format incorrect */
/** @} */

/** \cond DONT_DOXYGEN */
#define NRWAVEIO_MSGENULL 	"Null pointer"
#define NRWAVEIO_MSGEFILE 	"Error in file-IO"
#define NRWAVEIO_MSGENONULL 	"Not a Null pointer"
#define NRWAVEIO_MSGENOMEM 	"Memory ellocation error"
#define NRWAVEIO_MSGEVAL  	"Invalid value"
#define NRWAVEIO_MSGEFORMAT     "Meta data file format incorrect"
/** \endcond */


/**
 * Struct containing metadata information about a
 * single numerical relativity waveform.  This information
 * will be read from a metadata file. It is expected that
 * more elements will be added to this struct as required.
 */
typedef struct
tagNRWaveMetaData
{
  REAL8 massRatio; /**< Mass ratio m1/m2 where we assume m1 >= m2*/
  REAL8 spin1[3];  /**< Spin of m1 */
  REAL8 spin2[3];  /**< Spin of m2 */
  INT4  mode[2];   /**< l and m values */
  CHAR  filename[LALNameLength]; /**< filename where data is stored */
  /*   NumRelGroup group; */
}
NRWaveMetaData;


/** List of numrel waveform metadata */
typedef struct
tagNRWaveCatalog
{
  UINT4           length; /**< Number of waveforms */
  NRWaveMetaData  *data;  /**< List of waveforms */
}
NRWaveCatalog;

typedef struct
tagNumRelInjectParams
{
  NRWaveCatalog *nrCatalog;
  INT4 modeLlo;
  INT4 modeLhi;
  CHAR *ifo;
  REAL8 dynRange;
} NumRelInjectParams;

void LALReadNRWave_raw(LALStatus *status, REAL4TimeVectorSeries **out, const CHAR  *filename);

void LALReadNRWave_raw_real8(LALStatus *status, REAL8TimeVectorSeries **out, const CHAR  *filename);

void LALReadNRWave(LALStatus *status, REAL4TimeVectorSeries **out, const REAL4  mass, const CHAR  *filename);

void LALNRDataFind( LALStatus *status, NRWaveCatalog *out, const CHAR *dir, const CHAR *filename);

void LALGetSingleNRMetaData( LALStatus *status, NRWaveMetaData *data, const CHAR *dir, const CHAR *cfgstr);

void LALAddStrainModes( LALStatus *status, REAL4TimeVectorSeries  **outStrain,
		     NRWaveCatalog *nrCatalog, INT4 modeLlo, INT4 modeLhi, const SimInspiralTable *thisInj);

void LALDriveNRInject( LALStatus *status, REAL4TimeSeries *injData, SimInspiralTable *injections, NumRelInjectParams *params );

/** @} */ /* end:NRWaveIO_h */

#ifdef  __cplusplus
}                /* Close C++ protection */
#endif

#endif     /* Close double-include protection _NRWAVEIO_H */
