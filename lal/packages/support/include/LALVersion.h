#ifndef _LALVERSION_H
#define _LALVERSION_H

#include <lal/LALDatatypes.h>

#ifdef __cplusplus
extern "C" {
#endif

NRCSID( LALVERSIONH, "$Id$" );

#define LALVERSIONH_ENULL 1
#define LALVERSIONH_ESIZE 2
#define LALVERSIONH_ESPRN 4
#define LALVERSIONH_ESHRT 8

#define LALVERSIONH_MSGENULL "Null string pointer."
#define LALVERSIONH_MSGESIZE "Zero string size."
#define LALVERSIONH_MSGESPRN "Error in snprintf."
#define LALVERSIONH_MSGESHRT "String too short."

extern const char *lalVersion;
extern const int   lalVersionMajor;
extern const int   lalVersionMinor;
extern const char *lalConfigureArgs;
extern const char *lalConfigureDate;

#define LALVersionRequired( major, minor )      \
  ( LAL_VERSION_MAJOR > ( major ) ||            \
    ( LAL_VERSION_MAJOR == ( major ) && LAL_VERSION_MINOR >= ( minor ) ) )

void
LALVersion( LALStatus *status, CHAR *message, UINT4 size, INT4 verbose );

#ifdef __cplusplus
}
#endif

#endif _LALVERSION_H
