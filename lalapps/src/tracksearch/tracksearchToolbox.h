/**** <lalVerbatim file="TSDatgenHV"> *********
Author: Torres. C
$ID: tracksearch.h,v 1.0 2004/04/14 02:00:00 cristina Exp $
***** </lalVerbatim> **********************************/

#ifndef TRACKSEARCHTOOLBOX_H
#define TRACKSEARCHTOOLBOX_H

#include "tracksearch.h"
#include "tracksearchAverager.h"
#include <unistd.h>
#include <errno.h>

#define TSAPERMS 0666

/*
 * Error codes
 */
#define TRACKSEARCHTOOLBOXC_ENORM                          0
#define TRACKSEARCHTOOLBOXC_EFAIL                          1
#define TRACKSEARCHTOOLBOXC_EHEAD                          4
#define TRACKSEARCHTOOLBOXC_EDATA                          8
#define TRACKSEARCHTOOLBOXC_ERHEAD                         16
#define TRACKSEARCHTOOLBOXC_ERDATA                         32
#define TRACKSEARCHTOOLBOXC_EALLOC                          2
#define TRACKSEARCHTOOLBOXC_EWRITE                         64


#define  TRACKSEARCHTOOLBOXC_EMSGNORM            "TSToolbox Normal Exit"
#define  TRACKSEARCHTOOLBOXC_EMSGFAIL            "TSToolbox Subroutine Fail"
#define  TRACKSEARCHTOOLBOXC_EMSGHEAD            "TSToolbox Error writing candidate file header."
#define  TRACKSEARCHTOOLBOXC_EMSGDATA            "TSToolbox Error writing candidate file data."
#define  TRACKSEARCHTOOLBOXC_EMSGRHEAD           "TSToolbox Error reading candidate file header."
#define  TRACKSEARCHTOOLBOXC_EMSGRDATA           "TSToolbox Error reading candidate file data."
#define  TRACKSEARCHTOOLBOXC_EMSGALLOC           "TSToolbox Memory allocation problem."
#define  TRACKSEARCHTOOLBOXC_EMSGWRITE           "TSToolbox Error creating file for writing."

/*
 * Misc functions for lalapps tracksearch related software
 */

void
LALappsTSAReadMapFile(LALStatus*,
		      TSAMap**,
		      CHARVector*);

void
LALappsTSAWriteMapFile(LALStatus*,
		       TSAMap*,
		       CHARVector*);

void
LALappsTSACreateMap(LALStatus                    *status,
		    TSAMap                      **map,
		    TrackSearchMapMarkingParams  *imageBorders,
		    CreateTimeFreqIn             *createParams);

void
LALappsTSADestroyMap(LALStatus *status,
		     TSAMap    **map);

void
LALappsTSAWritePGM(LALStatus  *status,
		   TSAMap     *map,
		   CHARVector *UserNameOverride);

void
LALappsTSassert(UINT4          err,
		INT4           code,
		const CHAR    *msg);

/*
 * Bubble sort the cache list in RAM
 */
void
LALappsTSASortCache(LALStatus         *status,
		    TSAcache          *inputCache,
		    UINT4              ignoreMissing);
void
LALappsTSALoadCacheFile(LALStatus    *status,
			CHARVector   *filename,
			TSAcache    **mapCache);

void
LALappsTSADestroyCache(LALStatus       *status,
		       TSAcache       **mapCache);

void
LALappsDetermineFilename(LALStatus                    *status,
			 TrackSearchMapMarkingParams  imageBorders,
			 CHARVector                 **thisFilename, 
			 const CHAR*                  myFileExt);
void
LALappsCreateR4FromR8TimeSeries(LALStatus               *status,
				REAL4TimeSeries        **R4TS,
				REAL8TimeSeries         *R8TS);

void LALappsPSD_Check(REAL8TimeSeries*);

void print_real4tseries(const REAL4TimeSeries *fseries, 
			const char *file);
void print_real8tseries(const REAL8TimeSeries *fseries, 
			const char *file);
void print_real4fseries(const REAL4FrequencySeries *fseries,
			const char *file);
void print_real8fseries(const REAL8FrequencySeries *fseries, 
			const char *file);
void print_complex8fseries(const COMPLEX8FrequencySeries *fseries, 
			   const char *file);
void print_complex8_RandC_fseries(const COMPLEX8FrequencySeries *fseries, 
				  const char *file);
void print_lalUnit(LALUnit unit,const char *file);
#endif
 
