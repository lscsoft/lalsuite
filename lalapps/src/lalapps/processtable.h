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
    const CHAR          *program,
    const CHAR          *cvs_revision,
    const CHAR          *cvs_source,
    const CHAR          *cvs_date
    );

#ifdef  __cplusplus
}
#endif

#endif
