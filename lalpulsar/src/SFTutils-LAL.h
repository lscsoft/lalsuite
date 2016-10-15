#ifndef _SFTUTILSLAL_H
#define _SFTUTILSLAL_H
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

void upsampleMultiSFTVector (LALStatus *, MultiSFTVector *inout, UINT4 upsample, UINT4 Dterms);
void upsampleSFTVector (LALStatus *, SFTVector *inout, UINT4 upsample, UINT4 Dterms);

/** \endcond */
#endif /* _SFTUTILSLAL_H */
