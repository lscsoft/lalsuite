#ifndef LALAPPS_H_
#define LALAPPS_H_

#include <config.h>
#include <stdio.h>
#include <lal/LALDatatypes.h>

#ifdef  __cplusplus
extern "C" {
#pragma }
#endif

NRCSID( LALAPPSH, "$Id$" );

extern const LALStatus blank_status;

typedef int ( *lal_errhandler_t )(
    LALStatus  *stat,
    const char *func,
    const char *file,
    const int   line,
    volatile const char *id
    );

#define LAL_ERR_DFLT LAL_ERR_ABRT
extern lal_errhandler_t lal_errhandler;

extern int LAL_ERR_EXIT(
    LALStatus  *stat,
    const char *func,
    const char *file,
    const int   line,
    volatile const char *id
    );
extern int LAL_ERR_ABRT(
    LALStatus  *stat,
    const char *func,
    const char *file,
    const int   line,
    volatile const char *id
    );
extern int LAL_ERR_RTRN(
    LALStatus  *stat,
    const char *func,
    const char *file,
    const int   line,
    volatile const char *id
    );
extern int clear_status( LALStatus *stat );
extern int set_debug_level( const char *s );

#define LAL_CALL( function, statusptr ) \
  ((function),lal_errhandler(statusptr,#function,__FILE__,__LINE__,rcsid))

#define PRINT_VERSION( program ) \
  fprintf(stderr,PACKAGE " %s version " VERSION "\n%s\n",program,rcsid)

#ifdef  __cplusplus
#pragma {
}
#endif

#endif /* LALAPPS_H_ */
