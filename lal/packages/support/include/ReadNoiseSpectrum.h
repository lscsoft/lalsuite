#ifndef _READNOISESPECTRUMH_H
#define _READNOISESPECTRUMH_H

#include <lal/LALDatatypes.h>
#include <lal/Date.h>
#include <lal/TwoDMesh.h>

#ifdef  __cplusplus
extern "C" {
#pragma }
#endif


NRCSID (READNOISESPECTRUMH, "$Id$");

/* <lalLaTeX> 
\newpage\subsection*{Error codes} 
</lalLaTeX> */
/* <lalErrTable> */
#define LALREADNOISESPECTRUMH_ENULL 1
#define LALREADNOISESPECTRUMH_ENNUL 2
#define LALREADNOISESPECTRUMH_EALOC 3
#define LALREADNOISESPECTRUMH_EOPEN 4
#define LALREADNOISESPECTRUMH_EFCLO 5
#define LALREADNOISESPECTRUMH_EPARS 8

#define LALREADNOISESPECTRUMH_MSGENULL "Null pointer"
#define LALREADNOISESPECTRUMH_MSGENNUL "Non-null pointer"
#define LALREADNOISESPECTRUMH_MSGEALOC "Memory allocation error"
#define LALREADNOISESPECTRUMH_MSGEOPEN "Error opening file"
#define LALREADNOISESPECTRUMH_MSGEFCLO "Error closing file"
#define LALREADNOISESPECTRUMH_MSGEPARS "Error parsing spectrum file"
/* </lalErrTable> */

#define LALREADNOISESPECTRUM_MAXLINELENGTH 2048

void LALReadNoiseSpectrum(LALStatus *stat, REAL4FrequencySeries *spectrum, CHAR *fname);

#ifdef  __cplusplus
#pragma {
}
#endif

#endif /* _READNOISESPECTRUMH_H */
