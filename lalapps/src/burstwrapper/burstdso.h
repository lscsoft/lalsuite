#ifndef BURSTDSOH
#define BURSTDSOH

#include "datacondAPI/DatacondCaller.h"
#include <lal/StdBurstSearch.h>
#include <lal/LALDatatypes.h>
#include <lal/LIGOMetadataTables.h>

#include "FrameL.h"
#include <lal/FrameStream.h>

#ifdef  __cplusplus
extern "C" {
#endif

  /*****************************************************************/

  /* Define DEBUGBURST for debug info to be sent to disk */
  /*
#define DEBUGBURST 
  */

  /* big endian when not linux */
#ifdef linux
#undef WORDS_BIGENDIAN
#endif

  /* max number of triggers */
#define MAXNBURSTS 50000

  /* max size of a blob object */
#define MAXBLOB 800000

  /* names of the channels used by the wrapper */
#define GW_STRAIN_DATA "GW_STRAIN_DATA"
#define GW_STRAIN_PSD "GW_STRAIN_PSD"

  /* name of the frame cache file */
#define CACHEFILENAME "FrCacheFile"

  /* abort macro for LAL functions */
#define SABORT(I,M) fprintf(stderr,"ERROR: %i\n%s\n%s:%i\n",I,M,__FILE__,__LINE__); return 1

  /* macro to check LAL status pointer */
#define RCHECKSTATUSPTR( statusptr ) \
  if ( (statusptr)->statusPtr->statusCode ) \
  { \
    SETSTATUS( statusptr, -1, "Recursive error" ); \
    DETATCHSTATUSPTR(status); \
    RETURN (status); \
  } else (void)(0)

  /*
#define LALCALL (x) ((x), \
  if ( stat.statusPtr->statusCode ) \
  { \
    fprintf(stderr,"LAL error: %s, line %i\n",__FILE__,__LINE__); \
    return 1; \
  } else (void)(0)
  */

  /* Assert/Abort macros for LAL functions */
#define SASSERT(p1,p3,p4) if(!(p1)) {SABORT(p3,p4);}

#define RASSERT(p1,p2,p3,p4) if(!(p1)) {ABORT(p2,p3,p4);}

#define TASSERT(p1,p2,p3,p4) if(!(p1)) {SABORT(p3,p4);}


  /* linked list for burst parameters */
  typedef struct tagBurstParameterList {
    
    struct tagBurstParameterList *next;

    BurstParameter param;

  } BurstParameterList;


  /* implemented ETGs */
  typedef enum {TFCLUSTERS,SLOPE,POWER} BurstETG;

  /* burst search structure */
  typedef struct tagBurstSearchParams {

    CHAR channel[LIGOMETA_CHANNEL_MAX]; /* channel name */
    
    UINT4 Ndata; /* number of data points */

    BurstParameter waveforms; /* linked list of waveform names */

    UINT4 Nwaveforms; /* number of waveforms */

    REAL4TimeVectorSeries *wave_pc; /* input data */

    BurstParameter waveform_amplitudes; /* linked list of amplitudes */

    UINT4 Nwaveform_amplitude; /* number of amplitudes */

    BurstParameter alpha; /* linked list of alpha */

    UINT4 Nalpha; /* numebr of alpha */

    BurstParameter delta; /* linked list of delta */

    UINT4 Ndelta; /* number of delta */

    BurstParameter psi; /* linked list of psi */

    UINT4 Npsi; /* number of psi */

    UINT4 MAXinj; 

    UINT4 Ninj; /* number of injections */

    BurstParameter injTimes; /* linked list of injection times */

    UINT4 Ninj_per_segment; /* number of injection times */

    BurstETG ETG; /* which ETG */

    BurstParameterList ETGparams; /* linked list of ETG parameters; each parameter can be a linked list */

    BOOLEAN ETGparamsLOCK;

    UINT4 Nparams; /* number of ETG parameters */

    UINT4 *NNparams; /* number of element of each ETG parameter */

    BurstOutputDataSegment data; /* data for Parameter Estimation function */

    BurstOutputParameters oparams; /* parameters for PEst function */

    BOOLEAN searchMaster; /* Am I searchMaster? */

    void (*ETGfun)(LALStatus *, EventIDColumn *, REAL4TimeVectorSeries *, BurstParameter *); /* pointer to LAL function implementing the ETG */

    INT4 outputMethod; /* flag for output method */

    BOOLEAN noLALBurstOutput; /* flag for disabling the PEst function */

    CHAR dsoOutput[1024]; /* output file */

    BOOLEAN allowNoResponse; /* flag to avoid using response files */

  } BurstSearchParams;


  void ParseParameter( 
		      LALStatus *stat,
		      BurstParameter *params,
		      CHAR *argv,
		      UINT4 *N
		      );

/******** <lalErrTable file="TFCInitSearchHErrTab"> ********/
#define INITSEARCHH_ENULL     1
#define INITSEARCHH_ENNUL     2
#define INITSEARCHH_EALOC     3
#define INITSEARCHH_EARGS     4
#define INITSEARCHH_ESPOS      5
#define INITSEARCHH_ENAN      6
#define INITSEARCHH_EAOR      7
#define INITSEARCHH_EPOS      8
#define INITSEARCHH_EMMF      9
#define INITSEARCHH_ENAM      10
#define INITSEARCHH_ETAB      11
#define INITSEARCHH_EIN      12
#define INITSEARCHH_ERANK     13
#define INITSEARCHH_EMAXBLOB 14
#define INITSEARCHH_EETG 15
#define INITSEARCHH_EETGDATA 16
#define INITSEARCHH_ENB 17
#define INITSEARCHH_EUI 18
#define INITSEARCHH_EFILE 19

#define INITSEARCHH_MSGENULL     "Null pointer"
#define INITSEARCHH_MSGENNUL     "Non-null pointer"
#define INITSEARCHH_MSGEALOC     "Memory allocation error"
#define INITSEARCHH_MSGEARGS     "Wrong number of arguments"
#define INITSEARCHH_MSGESPOS      "Parameter <= 0 expected to be strictly positive"
#define INITSEARCHH_MSGENAN       "Numerical parameter is NaN"
#define INITSEARCHH_MSGEAOR      "Parameter expected to be in [0,1] is out of range"
#define INITSEARCHH_MSGEPOS      "Parameter < 0 when expected to be positive or zero"
#define INITSEARCHH_MSGEMMF      "Parameter minf > maxf"
#define INITSEARCHH_MSGENAM       "Channel/IFO name is zero length"
#define INITSEARCHH_MSGETAB       "Unknown table name"
#define INITSEARCHH_MSGEIN       "invalid input"
#define INITSEARCHH_MSGERANK     "rank of node"
#define INITSEARCHH_MSGEMAXBLOB "Maximum allowed BLOB size (1 Mo) exceeded"
#define INITSEARCHH_MSGEETG "ETG function generated a SIGSEGV or SIGBUS signal!"
#define INITSEARCHH_MSGEETGDATA "ETG function modifies its input!"
#define INITSEARCHH_MSGENB "Maximum number of burst exceeded!"
#define INITSEARCHH_MSGEUI "Unimplemented feature!"
#define INITSEARCHH_MSGEFILE "File error!"

/******** </lalErrTable> ********/


  /*****************************************************************/

typedef enum {Unknown, TimeSeries, FrequencySeries, Other1D, TimeFrequency, Wavelets, MultiDimensional} ObjectType;

void
LALFrGetSeriesType(
    LALStatus   *status,
    LALTYPECODE *output,
    ObjectType *objtype,
    FrChanIn    *chanin,
    FrStream    *stream
    );

int OutputSymbols(char *algorithms, 
		 int *Nsymbols,
		  datacond_symbol_type **symbols);

int getFrameCache(char *fQuery, 
		  char *dataserver);

int getFrameData(char *fQuery, 
		 char *dataserver,
		 int *Nsymbols,
		 datacond_symbol_type **symbols);

int getNonFrameData(char *rFiles, 
		    char **algorithms, 
		    int *Nsymbols,
		    datacond_symbol_type **symbols);

int InitSearch(char *filterParams, 
	       BurstSearchParams *bparams);

int ConditionData(int Nsymbols,
		  datacond_symbol_type *symbols,
		  char *algorithms);

int ReConditionData(int Nsymbols,
		    datacond_symbol_type *symbols,
		    char *algorithms,
		    BurstSearchParams *bparams);

int RunSearch(BurstSearchParams *bparams,
	      double f0, 
	      double f1
	      );

int ExtractResponse( COMPLEX8FrequencySeries *respPtr, 
		     int Nsymbols, 
		     datacond_symbol_type *symbols, 
		     char *ifo );

int symbolsfree(datacond_symbol_type *symbols);

int bparamsfree(BurstSearchParams *bparams); 

#ifdef  __cplusplus
}
#endif

#endif
