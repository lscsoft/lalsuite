/* Extras for BOINC compilation of HierarchicalSearch
   Author: Bernd Machenschalk
*/

#include <sys/types.h>
#include <unistd.h>
#include <stdio.h>

/* BOINC includes */
#include "filesys.h"

#include <lal/LALError.h>
#include <lal/LALRCSID.h>
NRCSID(HSBOINCEXTRASHRCSID,"$Id$");
#include "FstatToplist.h"

/* linking proper functions to the hooks in HierarchicalSearch.c */

#define SHOW_PROGRESS show_progress
#define fopen boinc_fopen

/* not implemented yet

#define GET_CHECKPOINT init_and_read_checkpoint
#define HS_CHECKPOINTING
#define INSERT_INTO_FSTAT_TOPLIST add_candidate_and_checkpoint
#define CLEANUP_CHECKPOINTING

so for now we define only dummies: */

#define GET_CHECKPOINT(filename)
#define CLEANUP_CHECKPOINTING
#undef  HS_CHECKPOINTING
#define INSERT_INTO_FSTAT_TOPLIST insert_into_fstat_toplist

#ifdef  __cplusplus
extern "C" {
#endif

extern LALStatus *global_status;

/* function prototypes, they are defined in boinc_extras.c */

/* allows the App to register another output file to be put into the
   zip archive that is sent back to the server */
extern void register_output_file(char*filename);

/* show progress of the App.
   NOTE: This also set the count & total (skypos) for checkpointing */
extern void show_progress(double rac, double dec, UINT4 count, UINT4 total);

/* inits checkpointing for the toplist and reads the last checkpoint if present
   This expects all passed variables (toplist, total, count) to be already initialized.
   The variables are modified only if a previous checkpoint was found.
   If *cptname (name of the checkpoint file) is NULL,
   the name is constructed by appending ".cpt" to the output filename.
   The FILE* should be the one that checpointed_fopen() above has returned. */
extern void init_and_read_checkpoint(toplist_t*toplist,
				     UINT4*total, UINT4*count,
				     char*outputname, char*cptname);

/* This corresponds to insert_into_fstat_toplist().
   It inserts a candidate into the toplist, updates the file
   and "compacts" it if necessary (i.e. bytes > maxsize).
   NOTE that the toplist parameter is just a dummy to make the interface
        compatible to insert_into_fstat_toplist(). The operation is
        actually performed on the toplist passed to the least recent call
        of init_and_read_checkpoint(), which, however, should be the same
        in all reasonable cases. */
extern int add_candidate_and_checkpoint (toplist_t*toplist, FstatOutputEntry cand);

/* does the final (compact) write of the file and cleans up checkpointing stuff
   The checkpoint file remains there in case something goes wrong during the rest */
extern void write_and_close_checkpointed_file (void);



/* the main() function of HierarchicalSerach.c becomes the extern MAIN(),
   the real main() function of the BOINC App is defined in boinc_extras.c
*/
extern int MAIN(int,char**);

#ifdef  __cplusplus
}
#endif

