#include <lal/LIGOMetadataTables.h>
#include "metaio.h"

#define INJECTIONS 1
#define TRIGGERS   0

#define EVENTUTILSH_ENULLP 1
#define EVENTUTILSH_EGETRO 2
#define EVENTUTILSH_EFILE  3
#define EVENTUTILSH_MSGENULLP "Null pointer"
#define EVENTUTILSH_MSGEGETRO "Error getting row from table"
#define EVENTUTILSH_MSGEFILE  "Error opening file"

typedef struct
tagSnglInspiralIndex
{
    INT4 ifoIndex;
    INT4 searchIndex;
    INT4 channelIndex;
    INT4 end_timeIndex;
    INT4 end_time_nsIndex;
    INT4 impulse_timeIndex;
    INT4 impulse_time_nsIndex;
    INT4 template_durationIndex;
    INT4 event_durationIndex;
    INT4 amplitudeIndex;
    INT4 eff_distanceIndex;
    INT4 coa_phaseIndex;
    INT4 mass1Index;
    INT4 mass2Index;
    INT4 mchirpIndex;
    INT4 etaIndex;
    INT4 tau0Index;
    INT4 tau2Index;
    INT4 tau3Index;
    INT4 tau4Index;
    INT4 tau5Index;
    INT4 ttotalIndex;
    INT4 snrIndex;
    INT4 chisqIndex;
    INT4 chisq_dofIndex;
    INT4 sigmasqIndex;
}
SnglInspiralIndex;

typedef struct
tagSimInspiralIndex
{
    INT4 geocent_end_timeIndex;
    INT4 geocent_end_time_nsIndex;
    INT4 end_time_gmstIndex;
    INT4 sourceIndex;
    INT4 mtotalIndex;
    INT4 etaIndex;
    INT4 distanceIndex;
    INT4 longitudeIndex;
    INT4 latitudeIndex;
    INT4 inclinationIndex;
    INT4 coa_phaseIndex;
    INT4 polarizationIndex;
    INT4 eff_dist_hIndex;
    INT4 eff_dist_lIndex;
    INT4 simulation_idIndex;
}
SimInspiralIndex;

    

typedef struct
tagSearchSummaryIndex
{
    INT4 process_idIndex;
    INT4 shared_objectIndex;
    INT4 lalwrapper_cvs_tagIndex;
    INT4 lal_cvs_tagIndex;
    INT4 commentIndex;
    INT4 in_start_timeIndex;
    INT4 in_start_time_nsIndex;
    INT4 in_end_timeIndex;
    INT4 in_end_time_nsIndex;
    INT4 out_start_timeIndex;
    INT4 out_start_time_nsIndex;
    INT4 out_end_timeIndex;
    INT4 out_end_time_nsIndex;
    INT4 neventsIndex;
    INT4 nnodesIndex;
}
SearchSummaryIndex;

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
tagsnglIFO{
    char *cfgfile;
    char ifoname[2];
    int  dqbit;
    struct tagcandParams candidates;
    struct tagvetoParams *veto_entry;
    struct tagtimeWindow *awindows;
    struct tagtimeWindow *vwindows;
    struct tagcandEvent *eventhead;
    struct tagcandEvent *Ieventhead;
    double start_time;
    double end_time;
    float time_analyzed; 
    float safety;
}
snglIFO;

typedef struct
tagmultiInspiral{
    int dqbit;
    float snr[8];
    double time[8];
    float effDistance[8];
    float mchirp[8];
    struct tagmultiInspiral *next_event;
}
multiInspiral;

typedef struct
tagtimeWindow{
    double                start_time;   /* gps seconds */
    double                end_time;     /* gps seconds */
    double                snr;          /* loudest snr in window */
    double                ratio;        /* ASQ/VETO ratio above which is veto */
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
    double     ratio;                   /* ASQ/VETO ratio above which is veto */
    timeWindow *vwindows; 
    struct tagvetoParams *next_veto;
}
vetoParams;

typedef struct
tagcandEvent{
    double time;
    float  snr;
    float  chisq;
    float  eff_distance;
    float  mchirp;
    float  mass1;
    float  mass2;
    int significance;
    int candidate;
    int coincident;
    struct tagcandEvent *next_event;
}
candEvent;


int buildVetoTimes( vetoParams *thisentry);
int buildEventList( candEvent **eventhead, timeWindow *vwindows, candParams candidates,
        int injectflag, int maxflag, float calfudge);
int resolveVetoTimes( timeWindow **vwindows, vetoParams *thisentry);
int build2DHistogram(candEvent *eventhead, const char *outputfile, int **hist, int nbins, float minsnr, float maxchisq);
int computeUL(const char *outputfile, int **triggerHistogram,
        int **injectHistogram, int numbins, float minsnr, float maxchisq,
        int ninject, float time_analyzed);
int buildDataQaulity(int **coincident_times, snglIFO *ifo, int numIFO,
        double *dummyStart, double *dummyEnd);
int buildMultiInspiralEvents(multiInspiral **multInspEv, int *coincident_times,
        snglIFO *ifo, int numIFO, int injectflag, double dummyStart, float delm,
        float distance, float coincidence_window);


void
buildSnglInspiralIndex(
        LALStatus             *status,
        const MetaioParseEnv   triggerEnv,
        SnglInspiralIndex        *params
        );

void
getSnglInspiralEvent(
        LALStatus             *status,
        const MetaioParseEnv   triggerEnv,
        SnglInspiralTable        *burstEvent,
        SnglInspiralIndex        *params
        );

void
LALClusterSnglInspiralTable (
	      LALStatus         *status,
              SnglInspiralTable *inspiralEvent,
              INT4              dtime,
	      INT4		clusterchoice
	      );

void readInspiralTriggers( 
        LALStatus *status, 
        SnglInspiralTable **eventList, 
        const CHAR *fname
        );
