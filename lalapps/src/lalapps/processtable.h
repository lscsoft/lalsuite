#ifndef PROCESSTABLE_H_
#define PROCESSTABLE_H_

#include <lal/LALDatatypes.h>
#include <lal/LIGOMetadataTables.h>

#ifdef  __cplusplus
extern "C" {
#endif

void
populate_process_table (
    LALStatus           *status,
    ProcessTable        *ptable,
    CHAR                *program,
    CHAR                *cvs_revision,
    CHAR                *cvs_source,
    CHAR                *cvs_date
    );

#ifdef  __cplusplus
}
#endif

#endif
