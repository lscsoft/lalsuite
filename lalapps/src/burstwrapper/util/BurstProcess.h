#ifndef BURSTPROCESSH
#define BURSTPROCESSH

#include <stdlib.h>
#include <sys/types.h>
#include <regex.h>

/*
#include <lal/EPSearch.h>
#include <BurstProcessTypeDef.h>
*/

#include <Coincidences.h>

#define regex(s1_,s2_) regexec(s1_,s2_,0,NULL,0)

/******************************************************************/
typedef struct tagBurstSegInfo {

  char *params; /* comma separated list of job parameters */
  char *info;   /* comma separated list of info from segment */

  struct tagBurstSegInfo *next;

} BurstSegInfo;

typedef struct tagBurstSegParams {

  char *params;

  int t0s;
  int t0ns;

} BurstSegParams;

int BurstProcess(
		 unsigned int *nback, /* number of background segments */
		 BurstSegInfo *backData, /* params,info for each segment */ 
		 int (*backFun)(char **, SnglBurstTableC *, void *, BurstSegParams *), /* function to process each segment */
		 void *backParams0,
		 size_t backSize,
		 unsigned int *nfore,
		 BurstSegInfo *foreData,
		 int (*foreFun)(char **, SnglBurstTableC *, void *, BurstSegParams *),
		 void *foreParams0,
		 size_t foreSize,
		 char **files,
		 int nfiles
		 );

int BurstProcessSub(
		 unsigned int *nback, /* number of background segments */
		 BurstSegInfo *backData, /* params,info for each segment */ 
		 int (*backFun)(char **, SnglBurstTableC *, void *, BurstSegParams *), /* function to process each segment */
		 void *backParams0,
		 size_t backSize,
		 unsigned int *nfore,
		 BurstSegInfo *foreData,
		 int (*foreFun)(char **, SnglBurstTableC *, void *, BurstSegParams *),
		 void *foreParams0,
		 size_t foreSize,
		 char **files,
		 int nfiles,
		 char *match
		 );

int BurstProcess2(
		  unsigned int *nback, /* number of background segments */
		  BurstSegInfo *backData, /* params,info for each segment */ 
		  int (*backFun)(char **, SnglBurstTableC *, SnglBurstTableC *, void *, BurstSegParams *, BurstSegParams *), /* function to process each segment */
		  void *backParams0,
		  size_t backSize,
		  unsigned int *nfore,
		  BurstSegInfo *foreData,
		  int (*foreFun)(char **, SnglBurstTableC *, SnglBurstTableC *, void *, BurstSegParams *, BurstSegParams *),
		  void *foreParams0,
		  size_t foreSize,
		  char **files1,
		  int nfiles1,
		  char **files2,
		  int nfiles2,
		  int mixETGParams,
		  char *IFO1,
		  char *IFO2
		  );

/******************************************************************/
typedef struct tagBackFunNumberOfUnclusteredEventsParams {
  int nBursts;
  int nSegments;

  int Nbands;
  double Fbands[10];
  int nBurstsBands[10];

} BackFunNumberOfUnclusteredEventsParams;

int BackFunNumberOfUnclusteredEvents(
				     char **info,
				     SnglBurstTableC *input,
				     void *parameters,
				     BurstSegParams *bparams
				     );


/******************************************************************/
typedef struct tagForeFunIsDetectedParams {

  int dType;
  
  int nDetected;
  int nInjected;

  double twin;
  double toff;

  int doFCut;
  double f0, f1;

} ForeFunIsDetectedParams;

int ForeFunIsDetected(
		      char **info,
		      SnglBurstTableC *input,
		      void *parameters,
		      BurstSegParams *bparams
		      );

/******************************************************************/
typedef struct tagEstimationErrorsParams {
  
  int nDetected;
  int nInjected;

  double dt, dt2;
  double df, df2;
  double dur, dur2;
  double bw, bw2;
  double h, h2;
  double conf, conf2;
  double snr, snr2;

  double bw0;
  double dur0;
  double f0;

  double twin;
  double toff;

  int outfile;

  int reportAll;

} EstimationErrorsParams;

int EstimationErrors(
		     char **info,
		     SnglBurstTableC *input,
		     void *parameters,
		     BurstSegParams *bparams
		     );

/******************************************************************/

typedef struct tagBackFunNumberOfUnclusteredEventsC2Params {

  Coincidence2Params *cparams;

  double nBursts;
  int nSegments; 

  double dtmin;
  double dtmax;
  double dt;

} BackFunNumberOfUnclusteredEventsC2Params;


int BackFunNumberOfUnclusteredEventsC2(
				       char **, 
				       SnglBurstTableC *, 
				       SnglBurstTableC *, 
				       void *, 
				       BurstSegParams *, 
				       BurstSegParams *
				       );

/******************************************************************/
typedef struct tagForeFunIsDetectedC2Params {
  
  Coincidence2Params *cparams;

  int nDetected;
  int nInjected;

  double twin;
  double toff;

} ForeFunIsDetectedC2Params;

int ForeFunIsDetectedC2(
			char **info,
			SnglBurstTableC *input1,
			SnglBurstTableC *input2,
			void *parameters,
			BurstSegParams *bparams1,
			BurstSegParams *bparams2
		      );
#endif
