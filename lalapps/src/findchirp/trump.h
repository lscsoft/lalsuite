#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <lal/Thresholds.h>
#include "metaio.h"
#include <lal/LALStdlib.h>

#define MAXSTR       256
#define MAXDATA  10
#define MAXVETODATA  5
#define MAXINJECTDATA  6
#define MAXCANDDATA  5
#define DT 10.01
#define INJECTIONS 1
#define TRIGGERS   0

#define RESPONSEC_ENORM  0
#define RESPONSEC_ESUB   1
#define RESPONSEC_EARG   2
#define RESPONSEC_EVAL   3
#define RESPONSEC_EFILE  4
#define RESPONSEC_EINPUT 5
#define RESPONSEC_EMEM   6

#define RESPONSEC_MSGENORM  "Normal exit"
#define RESPONSEC_MSGESUB   "Subroutine failed"
#define RESPONSEC_MSGEARG   "Error parsing arguments"
#define RESPONSEC_MSGEVAL   "Input argument out of valid range"
#define RESPONSEC_MSGEFILE  "Could not open file"
#define RESPONSEC_MSGEINPUT "Error reading file"
#define RESPONSEC_MSGEMEM   "Out of memory"

/* Usage format string. */
#define USAGE "Usage: %s [options]\n"\
        "--help                 Print this help message\n" \
        "[--veto vetofile]      File with veto metadata\n" \
        "[--trigger trigfile]   File with trigger/injection metadata\n" \
        "[--times file safety]  File with t dt and how much data per seg is ignored\n"\
        "[--coincidence file1 file2 window] Two files for coincidence and dt\n"

typedef struct
tagtimeWindow{
    double                start_time;   /* gps seconds */
    double                end_time;     /* gps seconds */
    struct tagtimeWindow *next_time_window;
}
timeWindow;

typedef struct
tagvetoParams{
    char       name[256];                   /* A veto name                    */
    char       filename[256];               /* XML filename                   */
    char       table_column[256];           /* Name of database column        */
    int        i_table_column;              /* An index TBD of column         */
    float      threshold;                   /* col > thresh => veto           */ 
    double     minusdtime;                  /* veto t_veto - dt < t ....      */
    double     plusdtime;                   /*  ...... < t_veto + dt          */
    timeWindow *vwindows; 
    struct tagvetoParams *next_veto;
}
vetoParams;

typedef struct
tagcandParams{
    char       name[256];
    char       injectfile[256];
    char       triggerfile[256];
    char       injepochs[256];
    float      snr_threshold;
    float      chi_threshold;
    double     dtime;
}
candParams;

typedef struct
tagcandEvent{
    double time;
    float  snr;
    float  chisq;
    struct tagcandEvent *next_event;
}
candEvent;


int buildVetoTimes( vetoParams *thisentry);
int buildEventList( candEvent **eventhead, timeWindow *vwindows, candParams candidates,
        int injectflag);
int resolveVetoTimes( timeWindow **vwindows, vetoParams *thisentry);
int build2DHistogram(candEvent *eventhead, const char *outputfile, int **hist, int nbins, float minsnr, float maxchisq);
int computeUL(const char *outputfile, int **triggerHistogram,
        int **injectHistogram, int numbins, float minsnr, float maxchisq,
        int ninject, float time_analyzed);
