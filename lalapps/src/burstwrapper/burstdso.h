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

#define DEBUGBURST 

#ifdef linux
#undef WORDS_BIGENDIAN
#endif

#define MAXNBURSTS 50000
#define MAXBLOB 800000

#define GW_STRAIN_DATA "GW_STRAIN_DATA"
#define GW_STRAIN_PSD "GW_STRAIN_PSD"

#define SABORT(I,M) fprintf(stderr,"ERROR: %i\n%s\n%s:%i\n",I,M,__FILE__,__LINE__); return 1

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

#define SASSERT(p1,p3,p4) if(!(p1)) {SABORT(p3,p4);}

#define RASSERT(p1,p2,p3,p4) if(!(p1)) {ABORT(p2,p3,p4);}

#define TASSERT(p1,p2,p3,p4) if(!(p1)) {SABORT(p3,p4);}


  typedef struct tagBurstParameterList {
    
    struct tagBurstParameterList *next;

    BurstParameter param;

  } BurstParameterList;


  typedef enum {TFCLUSTERS,SLOPE,POWER} BurstETG;


  typedef struct tagBurstSearchParams {

    CHAR channel[LIGOMETA_CHANNEL_MAX];
    
    UINT4 Ndata;

    BurstParameter waveforms;

    UINT4 Nwaveforms;

    REAL4TimeVectorSeries *wave_pc;

    BurstParameter waveform_amplitudes;

    UINT4 Nwaveform_amplitude;

    BurstParameter alpha;

    UINT4 Nalpha;

    BurstParameter delta;

    UINT4 Ndelta;

    BurstParameter psi;

    UINT4 Npsi;

    UINT4 MAXinj;

    UINT4 Ninj;

    BurstParameter injTimes;

    UINT4 Ninj_per_segment;

    BurstETG ETG;

    BurstParameterList ETGparams;

    BOOLEAN ETGparamsLOCK;

    UINT4 Nparams;

    UINT4 *NNparams;

    BurstOutputDataSegment data;

    BurstOutputParameters oparams;

    BOOLEAN searchMaster;

    void (*ETGfun)(LALStatus *, EventIDColumn *, REAL4TimeVectorSeries *, BurstParameter *);

    INT4 outputMethod;

    BOOLEAN noLALBurstOutput;

    CHAR dsoOutput[1024];

    BOOLEAN allowNoResponse;

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
