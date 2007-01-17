/* Extras for BOINC compilation of HierarchicalSearch
   Author: Bernd Machenschalk
*/

/* BOINC includes */
#include <sys/types.h>
#include <unistd.h>
#include "filesys.h"

#include <lal/LALError.h>
#include <lal/LALRCSID.h>
NRCSID(HSBOINCEXTRASHRCSID,"$Id$");
#include "FstatToplist.h"

/* linking proper functions to the hooks in HierarchicalSearch.c */
/* not implemented yet at all

#define SET_CHECKPOINT(filename,rac,dec,tpl_count,tpl_total)\
        set_checkpoint(filename,rac,dec,tpl_count,tpl_total)
#define GET_CHECKPOINT(filename)\
        get_checkpoint(filename)

so for now we define only dummies: */

#define SET_CHECKPOINT(filename,rac,dec,count,total) filename = filename;
#define GET_CHECKPOINT(filename)
#define REMOVE_CHECKPOINT(filename)


#define SHOW_PROGRESS(rac,dec,count,total)\
        show_progress(rac,dec,count,total)

#define fopen boinc_fopen

#ifdef  __cplusplus
extern "C" {
#endif

extern LALStatus *global_status;

/* function prototypes, they are defined in boinc_extras.c */

/* allows the App to register another output file to be put into the
   zip archive that is sent back to the server */
extern void register_output_file(char*filename);

/* show progress of the App.
   This also set the count & total (skypos) for checkpointing */
extern void show_progress(double rac, double dec, long count, long total);

/* This function is meant as a "dummy" replacement for fopen(output_file)
   It doesn't actually return a FILE*, but just a (casted) pointer to
   the filename. */
extern FILE* checpointed_fopen(cher*filename,char*dummy);

/* inits checkpointing for the toplist and reads the last checkpoint if present
   This expects all passed variables (toplist, total, count) to be already initialized.
   The variables are modified only if a previous checkpoint was found.
   If *cptname (name of the checkpoint file) is NULL,
   the name is constructed by appending ".cpt" to the output filename.
   The FILE* should be the one that checpointed_fopen() above has returned. */
extern void init_and_read_checkpoint(FILE* fp, toplist_t*toplist,
				     unsigned long*total, unsigned long*count,
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
   The checkpoint file remains there in case something goes wrong during the rest
   The interface of the function is that of fclose(), so it can be used to replace
   the usual fclose(output_filepointer) */
extern int write_and_close_checkpointed_file (FILE*fp);



/* the main() function of HierarchicalSerach.c becomes the extern MAIN(),
   the real main() function of the BOINC App is defined in boinc_extras.c
*/
extern int MAIN(int,char**);

#ifdef  __cplusplus
}
#endif

