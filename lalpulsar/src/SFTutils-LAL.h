#ifdef _SFTUTILS_H  /* Only include via SFTutils.h */
#ifndef SWIG /* exclude from SWIG interface */
/** \cond DONT_DOXYGEN */

/*----- Error-codes -----*/

#define SFTUTILS_ENULL 		1
#define SFTUTILS_ENONULL	2
#define SFTUTILS_EMEM		3
#define SFTUTILS_EINPUT		4
#define SFTUTILS_EFUNC		6

#define SFTUTILS_MSGENULL 	"Arguments contained an unexpected null pointer"
#define SFTUTILS_MSGENONULL	"Output pointer is not NULL"
#define SFTUTILS_MSGEMEM	"Out of memory"
#define SFTUTILS_MSGEINPUT	"Invald input parameter"
#define SFTUTILS_MSGEFUNC	"Sub-routine failed"

void LALCreateSFTtype (LALStatus *status, SFTtype **sft, UINT4 SFTlen);
void LALCreateSFTVector (LALStatus *status, SFTVector **sftvect, UINT4 numSFTs, UINT4 SFTlen);
void LALCreateMultiSFTVector ( LALStatus *status, MultiSFTVector **out, UINT4 length, UINT4Vector *numsft );

void upsampleMultiSFTVector (LALStatus *, MultiSFTVector *inout, UINT4 upsample, UINT4 Dterms);
void upsampleSFTVector (LALStatus *, SFTVector *inout, UINT4 upsample, UINT4 Dterms);

void LALDestroySFTtype (LALStatus *status, SFTtype **sft);
void LALDestroySFTVector (LALStatus *status, SFTVector **sftvect);
void LALDestroyPSDVector (LALStatus *status, PSDVector **vect);
void LALDestroyMultiSFTVector (LALStatus *status, MultiSFTVector **multvect);
void LALDestroyMultiPSDVector (LALStatus *status, MultiPSDVector **multvect);

void LALCopySFT (LALStatus *status, SFTtype *dest, const SFTtype *src);

void LALSubtractSFTVectors (LALStatus *, SFTVector **outVect, const SFTVector *inVect1, const SFTVector *inVect2 );
void LALLinearlyCombineSFTVectors (LALStatus *, SFTVector **outVect, SFTVector **inVects, const COMPLEX16Vector *weights, const CHAR *outName);
void LALAppendSFT2Vector (LALStatus *, SFTVector *vect, const SFTtype *sft );

void LALCreateTimestampVector (LALStatus *status, LIGOTimeGPSVector **vect, UINT4 len);
void LALDestroyTimestampVector (LALStatus *status, LIGOTimeGPSVector **vect);

void LALMakeTimestamps (LALStatus *, LIGOTimeGPSVector **timestamps, const LIGOTimeGPS tStart, REAL8 duration, REAL8 Tsft);

void LALComputeNoiseWeights  (LALStatus *status, REAL8Vector *weightV, const SFTVector *sftVect, INT4 blkSize, UINT4 excludePercentile);
void LALComputeMultiNoiseWeights  (LALStatus *status, MultiNoiseWeights **weightsV, const MultiPSDVector *multipsd, UINT4 blocksRngMed, UINT4 excludePercentile);
void LALDestroyMultiNoiseWeights  (LALStatus *status, MultiNoiseWeights **weights);

void LALGetSFTtimestamps (LALStatus *, LIGOTimeGPSVector **timestamps, const SFTVector *sfts );

/** \endcond */
#endif /* SWIG */
#endif /* _SFTUTILS_H */
