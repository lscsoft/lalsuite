#include <lal/LALDatatypes.h>
#include <lal/LIGOMetadataTables.h>

void
populate_process_table (
    LALStatus           *status,
    ProcessTable        *ptable,
    CHAR                *program,
    CHAR                *cvs_revision,
    CHAR                *cvs_source,
    CHAR                *cvs_date
    );
