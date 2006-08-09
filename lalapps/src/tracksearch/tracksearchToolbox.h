/**** <lalVerbatim file="TSDatgenHV"> *********
Author: Torres. C
$ID: tracksearch.h,v 1.0 2004/04/14 02:00:00 charlie Exp $
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

#define  TRACKSEARCHTOOLBOXC_EMSGNORM            "TSToolbox Normal Exit"
#define  TRACKSEARCHTOOLBOXC_EMSGFAIL            "TSToolbox Subroutine Fail"

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
#endif
