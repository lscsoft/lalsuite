/* $Id$ */

#ifndef WrapperInterfaceH
#define WrapperInterfaceH

#include <wrapperInterfaceDatatypes.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifdef LDAS_BUILD
#define LDAS_EXTERN extern
#else
#define LDAS_EXTERN
#endif

/* Declarations for search functions( USE WHEN NEW SO IS READY ) */
LDAS_EXTERN INT4 initSearch( CHAR** initStatus, InitParams* initParams );
LDAS_EXTERN INT4 conditionData( CHAR** conditionStatus, inPut* data,
    SearchParams* searchParams );
LDAS_EXTERN INT4 applySearch( CHAR** searchStatus, inPut* input,
    SearchOutput* output, SearchParams* searchParams );
LDAS_EXTERN INT4 freeOutput( CHAR** freeStatus, SearchOutput* output );
LDAS_EXTERN INT4 finalizeSearch( CHAR** finalizeStatus );

#undef LDAS_EXTERN

#ifdef __cplusplus
}
#endif

#endif /* WrapperInterfaceH */
