#ifndef INJSGNL_H
#define INJSGNL_H

/*
 *
 * Routine to inject a signal with parameters read from a LIGOLw-format file.
 *
 */

#include <lal/LALDatatypes.h>

typedef enum { burst_inject, inspiral_inject } inject_type;
int inject_signal( REAL4TimeSeries *series, int injectSignalType, 
    const char *injectFile, const char *calCacheFile, REAL4 responseScale );

#endif /* INJSGNL_H */
