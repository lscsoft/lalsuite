#define INJECTIONS 1
#define TRIGGERS   0

#define EVENTUTILSH_ENULLP 1
#define EVENTUTILSH_MSGENULLP "null pointer"
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
