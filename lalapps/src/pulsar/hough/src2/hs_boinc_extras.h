/* Extras for BOINC compilation of HierarchicalSearch
   Author: Bernd Machenschalk
*/

#include <lal/LALRCSID.h>
NRCSID(HSBOINCEXTRASHRCSID,"$Id$");

/* BOINC includes */
#include <sys/types.h>
#include <unistd.h>
#include "filesys.h"

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
extern "C"{
#endif

/* function prototypes, they are defined in boinc_extras.c */
extern void set_checkpoint(char*filename,double rac,double dec, long count, long total);
extern void get_checkpoint(char*filename);
extern void show_progress(double rac, double dec, long count, long total);

/* the main() function of HierarchicalSerach.c becomes the main_hierarchical_search()
   the real main() function of the BOINC App is defined in boinc_extras.c
*/
extern int MAIN(int,char**);

#ifdef  __cplusplus
}
#endif

