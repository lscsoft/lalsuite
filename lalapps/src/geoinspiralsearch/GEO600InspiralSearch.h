/* <lalVerbatim file="GEO600InspiralSearchHV"> 

Author: Sathyaprakash, B.S.

</lalVerbatim> */

/* <lalLaTeX>

\section{Header \texttt{GEO600InspiralSearch.h}}
\label{s:GEO600InspiralSearch.h}

Header file for GEO600 inspiral search functions.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/GEO600InspiralSearch.h>
\end{verbatim}

\noindent This header file covers routines that are used in 
GEO600's inspiral search using an MPI code.

</lalLaTeX> */

#ifndef _GEO600INSPIRALSEARCH_H 
#define _GEO600INSPIRALSEARCH_H 

#ifdef  __cplusplus
extern "C" {
#endif

NRCSID( GEO600INSPIRALSEARCHH, "$Id$" );

/* <lalLaTeX>

\subsection*{Error codes}

</lalLaTeX>  */

/* <lalErrTable> */

#define GEO600INSPIRALSEARCHH_ENULL 1
#define GEO600INSPIRALSEARCHH_EMEM 2
#define GEO600INSPIRALSEARCHH_ECHOICE 4
#define GEO600INSPIRALSEARCHH_EDIV0 8
#define GEO600INSPIRALSEARCHH_ESIZE 16
#define GEO600INSPIRALSEARCHH_MSGENULL "Arguments contained an unexpected null pointer"
#define GEO600INSPIRALSEARCHH_MSGEMEM "Memory allocation error"
#define GEO600INSPIRALSEARCHH_MSGECHOICE "Invalid choice for an input parameter"
#define GEO600INSPIRALSEARCHH_MSGEDIV0 "Division by zero"
#define GEO600INSPIRALSEARCHH_MSGESIZE "Invalid input size"

/* </lalErrTable> */

/* <lalLaTeX>

\section*{Structures}
\input{GEO600InspiralSearchHS}
</lalLaTeX>  */

#include <lal/LALNoiseModels.h>

typedef enum
{
	realPSD, syntheticPSD
}
SimulationType;

typedef enum
{
	LoudestEvent, AllTemplates, AllEvents
}
FindEventsType;

/*  <lalVerbatim file="GEO600InspiralSearchHS"> */
typedef struct
tagInspiralSearchParams
{
   char  *channelName;      /* Name of the channel to search in */

   INT4  baseDirectory;     /* Base directory */
   UINT4 NumTrials;         /* Number of trials */
   UINT4 startGPSTime;      /* GPS time at which to begin search */
   UINT4 currentGPSTime;    /* GPS time at which to end search */
   UINT4 endGPSTime;        /* GPS time at which to end search */
   UINT4 dataSetDuration;   /* Duration of each dataSet in integer power of 2 */
   UINT4 dataSetLength;     /* Lenght of each dataSet in integer power of 2 */
   UINT4 samplingRate;      /* The rate at which to re-sample data */
                            /* New rate must be an integer power of 2, */
                            /* less than or equal to the original rate; */
                            /* if samplingRate is set to zero, default */
                            /* value will be used */
   UINT4 paddingFactor;     /* Factor by which dataSet.length is larger than */
                            /* the longest template; must be an int power of 2 */
   UINT4 numDataSetsForPSD; /* Number of data sets to be used to compute PSD */

   REAL8 fLo;               /* Lower frequency cutoff for band pass data */
   REAL8 fHi;               /* Upper frequency cutoff for band pass data */
   REAL8 df;                /* Frequency resolution */
   SimulationType PSDType;  /* Type of simulation required: With realPSD or syntheticPSD */
   FindEventsType findEventsType; /* Whether to search for loudest event or all events */
} 
InspiralSearchParams;
/*  </lalVerbatim>  */
/*  <lalLaTeX> \idx[Type]{InspiralSearchParams} </lalLaTeX>  */


/*  <lalVerbatim file="GEO600InspiralSearchHS"> */
typedef struct
tagRunningAveragePSDIn
{
   InspiralSearchParams searchParams;
   RealFFTPlan          *fwdp;
   RealFFTPlan          *revp;
   REAL4Vector          *window;
   REAL8                 norm;
} 
RunningAveragePSDIn;
/*  </lalVerbatim>  */
/*  <lalLaTeX> \idx[Type]{InspiralSearchParams} </lalLaTeX>  */

/*  <lalLaTeX>
\vfill{\footnotesize\input{GEO600InspiralSearchHS}}
</lalLaTeX>  */

/* Function prototypes */

/* <lalLaTeX>
\newpage\input{GEO600InspiralSearchMasterC}
</lalLaTeX>  */

void 
GEO600InspiralSearchMasterDB
   (
   LALStatus *status
   );

/* <lalLaTeX>
\newpage\input{GEO600InspiralSearchSlaveC}
</lalLaTeX>  */
void 
GEO600InspiralSearchSlave
   (
   LALStatus *status
   );

/* <lalLaTeX>
\newpage\input{GEO600InspiralSearchFuncsC}
</lalLaTeX>  */

void 
ReadInspiralSearchParams
   (
   LALStatus *status, 
   InspiralSearchParams *searchParams
   );

void 
ReadInspiralTemplateBankParams
   (
   LALStatus            *status, 
   InspiralCoarseBankIn *coarseIn
   );

void 
ReadInspiralFindEventsIn
   (
   LALStatus            *status, 
   InspiralFindEventsIn *eventsin
   );

void 
ReadRandomInspiralSignalParams
   (
   LALStatus              *status, 
   RandomInspiralSignalIn *randIn,
   InspiralCoarseBankIn   *coarseIn
   );


/* <lalLaTeX>
\newpage\input{WriteRandomInspiralSignalParamsC}
</lalLaTeX>  */
void 
WriteRandomInspiralSignalParams
   (
   LALStatus              *status, 
   FILE                   *logfile,
   RandomInspiralSignalIn *randIn
   );



/* <lalLaTeX>
\newpage\input{WriteInspiralEventsC}
</lalLaTeX>  */
void
WriteInspiralEvents
   (
   LALStatus          *status, 
   InspiralEventsList *eventList,
   INT4               nEvents,
   FILE               *EventsFile
   );



/* <lalLaTeX>
\newpage\input{WriteInspiralLoudestEventC}
</lalLaTeX>  */
void
WriteInspiralLoudestEvent
   (
   LALStatus          *status, 
   InspiralEventsList *eventList,
   UINT4              currentGPSTime,
   FILE               *EventsFile
   );



/* <lalLaTeX>
\newpage\input{WriteInspiralSearchParamsC}
</lalLaTeX>  */
void 
WriteInspiralSearchParams
   (
   LALStatus            *status, 
   FILE                 *logfile,
   InspiralSearchParams *searchParams
   );

/* <lalLaTeX>
\newpage\input{WriteInspiralTemplateBankParamsC}
</lalLaTeX>  */
void 
WriteInspiralTemplateBankParams
   (
   LALStatus            *status, 
   FILE                 *logfile,
   InspiralCoarseBankIn *coarseIn
   );


/* <lalLaTeX>
\newpage\input{WriteInspiralFindEventsInC}
</lalLaTeX>  */
void 
WriteInspiralFindEventsIn
   (
   LALStatus            *status, 
   FILE                 *logfile,
   InspiralFindEventsIn *eventsin
   );

/* 
 * This function down-samples an input dataset at the desired
 * sample-rate by low-pass filtering of the given dataset.
 * Low-pass filtering is done in the Fourier domain
 */

/* <lalLaTeX>
\newpage\input{RunningAverageDoublePSDC}
</lalLaTeX>  */
void 
RunningAverageDoublePSD
   (
   LALStatus           *status, 
   REAL8Vector         *psd, 
   INT4                *ret,
   RunningAveragePSDIn *averagePSDIn
   );

/* <lalLaTeX>
\newpage\input{WriteInspiralTemplateListC}
</lalLaTeX>  */
void 
WriteInspiralTemplateList 
   (
   LALStatus *status, 
   FILE *logfile, 
   INT4 nlist, 
   InspiralTemplateList *list
   );

/* <lalLaTeX>
\newpage\input{GEO600InspiralSearch}
</lalLaTeX>  */

/*  <lalVerbatim file="GEO600InspiralSearchHS"> */
  typedef struct LALDBConnect
{
    char *host;
    char *database;
    char *user;
    char *password;
    INT4 insert_id;
}
LALDBConnect;
/*  </lalVerbatim>  */
/*  <lalLaTeX> \idx[Type]{LALDBConnect} </lalLaTeX>  */


/*  <lalVerbatim file="GEO600InspiralSearchHS"> */
typedef struct LALProcessTable
{
       int creator_db;
       char *program;
       char *version;
       char *cvs_repository;
       int cvs_entry_time;
       char *comment;
       int is_online;
       char *node;
       char *username;
       int unix_procid;
       char* start_time;
       char* end_time;
       int jobid;
       char *domain;
       char *process_id;
       int param_set;
       char *ifos;
}
LALProcessTable;
/*  </lalVerbatim>  */
/*  <lalLaTeX> \idx[Type]{LALProcessTable} </lalLaTeX>  */


/*  <lalVerbatim file="GEO600InspiralSearchHS"> */
typedef struct LALProcessParamsTable
{
    int creator_db;
    char *program;
    int process_id;
    char *param;
    char *type;
    char *value;
}
LALProcessParamsTable;
/*  </lalVerbatim>  */
/*  <lalLaTeX> \idx[Type]{LALProcessParamsTable} </lalLaTeX>  */



/*  <lalVerbatim file="GEO600InspiralSearchHS"> */
typedef struct LALSnglInspiralTable
{
       int creator_db;
       int process_id;
       char *filter_id;
       char *ifo;
       char *search;
       int end_time;
       int end_time_ns;
       int impulse_time;
       int impulse_time_ns;
       double amplitude;
       double eff_distance;
       double coa_phase;
       double mass1;
       double mass2;
       double mchirp;
       double eta;
       double tau0;
       double tau2;
       double tau3;
       double tau4;
       double tau5;
       double ttotal;
       double snr;
       double chisq;
       int chisq_dof;
       double sigmasq;
       char *event_id;
       char *correlationEvent;
}
LALSnglInspiralTable;
/*  </lalVerbatim>  */
/*  <lalLaTeX> \idx[Type]{LALSnglInspiralTable} </lalLaTeX>  */

/* <lalLaTeX>
\newpage\input{WriteInspiralEventsC}
</lalLaTeX>  */

void
WriteInspiralEventsDB
(
   LALStatus          *status, 
   InspiralEventsList *eventList,
   INT4               nEvents,
   LALDBConnect       conn
);




/* <lalLaTeX>
\newpage\input{WriteInspiralLoudestEventC}
</lalLaTeX>  */
void
WriteInspiralLoudestEventDB
(
   LALStatus          *status, 
   InspiralEventsList *eventList,
   UINT4              currentGPSTime,
   LALDBConnect       conn
);



void SetUpDBConnectionInfo(LALDBConnect *conn);

void SetUpProcessTableInfo(LALProcessTable *procTable);

void WriteToLALProcessTable(LALProcessTable table, LALDBConnect *conn); 

void WriteToLALProcessParamsTable(LALProcessParamsTable table, LALDBConnect conn);

void AppendLALProcessTable(LALDBConnect conn);


#ifdef  __cplusplus
}
#endif

#endif /* _GEO600INSPIRALSEARCH_H */
