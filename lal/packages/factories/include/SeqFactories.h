/*-----------------------------------------------------------------------
 *
 * File Name: SeqFactories.h
 *
 * Author: Finn, L. S.
 *
 * Revision: $Id$
 *
 *-----------------------------------------------------------------------
 *
 * NAME
 * SeqFactories.h
 *
 * SYNOPSIS
 * #include "SeqFactories.h"
 *
 * DESCRIPTION 
 * Provides prototype and status code information for use
 * of CreateSequence, CreateVectorSequence, CreateArraySequence, 
 * DestroySequence, DestroyVectorSequence and DestroyArraySequence
 *
 * DIAGNOSTICS
 *
 *----------------------------------------------------------------------- 
 */

#ifndef _SEQFACTORIES_H
#define _SEQFACTORIES_H

#ifndef _LALDATATYPES_H
#include "LALDatatypes.h"
#ifndef _LALDATATYPES_H
#define _LALDATATYPES_H
#endif
#endif

#ifndef _AVFACTORIES_H
#include "AVFactories.h"
#ifndef _AVFACTORIES_H
#define _AVFACTORIES_H
#endif
#endif

NRCSID (SEQFACTORIESH, "$Id$");

#define CREATESEQ_ESLENGTH  1
#define CREATESEQ_EVLENGTH  2
#define CREATESEQ_EVPTR    4
#define CREATESEQ_EUPTR    8
#define CREATESEQ_EMALLOC  16

#define CREATESEQ_MSGESLENGTH   "Illegal length" /* ESLENGTH */
#define CREATESEQ_MSGEVLENGTH   "Illegal vector length" /* EVLENGTH */
#define CREATESEQ_MSGEVPTR    "seq == NULL" /* EVPTR */
#define CREATESEQ_MSGEUPTR    "*seq != NULL" /* EUPTR */
#define CREATESEQ_MSGEMALLOC   "Malloc failure" /* EMALLOC */

#define DESTROYSEQ_EVPTR  1
#define DESTROYSEQ_EUPTR  2
#define DESTROYSEQ_EDPTR  8

#define DESTROYSEQ_MSGEVPTR "seq == NULL" /* EVPTR */
#define DESTROYSEQ_MSGEUPTR "*seq == NULL" /* EUPTR */
#define DESTROYSEQ_MSGEDPTR "(*seq)->data == NULL" /* EDPTR */

#define CREATEVECSEQ_ESLENGTH    1
#define CREATEVECSEQ_EVPTR      2
#define CREATEVECSEQ_EUPTR      4
#define CREATEVECSEQ_EMALLOC    8
#define CREATEVECSEQ_EVLENGTH  16
#define CREATEVECSEQ_EINPTR    32

#define CREATEVECSEQ_MSGESLENGTH   "Illegal sequence length" /* ESLENGTH */
#define CREATEVECSEQ_MSGEVPTR     "seq == NULL" /* EVPTR */
#define CREATEVECSEQ_MSGEUPTR     "*seq != NULL" /* EUPTR */
#define CREATEVECSEQ_MSGEMALLOC   "Malloc failure" /* EMALLOC */
#define CREATEVECSEQ_MSGEVLENGTH  "Illegal vector length" /* EVLENGTH */
#define CREATEVECSEQ_MSGEINPTR    "in == NULL" /* EINPTR */

#define DESTROYVECSEQ_EVPTR  1
#define DESTROYVECSEQ_EUPTR  2
#define DESTROYVECSEQ_EDPTR  8

#define DESTROYVECSEQ_MSGEVPTR "seq == NULL" /* EVPTR */
#define DESTROYVECSEQ_MSGEUPTR "*seq == NULL" /* EUPTR */
#define DESTROYVECSEQ_MSGEDPTR "(*seq)->data == NULL" /* EDPTR */

typedef struct tagCreateVectorSequenceIn {
  INT4 length;
  INT4 vectorLength;
} CreateVectorSequenceIn;

/*
 * 9. Functions Declarations (i.e., prototypes).
 */

void CreateSequence(Status *, REAL4Sequence **, INT4);
void CHARCreateSequence(Status *, CHARSequence **, INT4);
void I2CreateSequence(Status *, INT2Sequence **, INT4);
void I4CreateSequence(Status *, INT4Sequence **, INT4);
void I8CreateSequence(Status *, INT8Sequence **, INT4);
void SCreateSequence(Status *, REAL4Sequence **, INT4);
void DCreateSequence(Status *, REAL8Sequence **, INT4);
void CCreateSequence(Status *, COMPLEX8Sequence **, INT4);
void ZCreateSequence(Status *, COMPLEX16Sequence **, INT4);

void DestroySequence(Status *, REAL4Sequence **);
void CHARDestroySequence(Status *, CHARSequence **);
void I2DestroySequence(Status *, INT2Sequence **);
void I4DestroySequence(Status *, INT4Sequence **);
void I8DestroySequence(Status *, INT8Sequence **);
void SDestroySequence(Status *, REAL4Sequence **);
void DDestroySequence(Status *, REAL8Sequence **);
void CDestroySequence(Status *, COMPLEX8Sequence **);
void ZDestroySequence(Status *, COMPLEX16Sequence **);

void CreateVectorSequence(Status *, 
                             REAL4VectorSequence **,
			     CreateVectorSequenceIn *);
void CHARCreateVectorSequence(Status *, 
                             CHARVectorSequence **,
			     CreateVectorSequenceIn *);
void I2CreateVectorSequence(Status *, 
			     INT2VectorSequence **,
			     CreateVectorSequenceIn *);
void I4CreateVectorSequence(Status *, 
			     INT4VectorSequence **,
			     CreateVectorSequenceIn *);
void I8CreateVectorSequence(Status *, 
			     INT8VectorSequence **,
			     CreateVectorSequenceIn *);
void SCreateVectorSequence(Status *, 
			     REAL4VectorSequence **,
			     CreateVectorSequenceIn *);
void DCreateVectorSequence(Status *, 
			     REAL8VectorSequence **,
			     CreateVectorSequenceIn *);
void CCreateVectorSequence(Status *, 
			     COMPLEX8VectorSequence **, 
			     CreateVectorSequenceIn *);
void ZCreateVectorSequence(Status *, 
			     COMPLEX16VectorSequence **, 
			     CreateVectorSequenceIn *);

void DestroyVectorSequence (Status *, 
                             REAL4VectorSequence **);
void CHARDestroyVectorSequence (Status *, 
                             CHARVectorSequence **);
void I2DestroyVectorSequence(Status *, 
			     INT2VectorSequence **);
void I4DestroyVectorSequence(Status *, 
			     INT4VectorSequence **);
void I8DestroyVectorSequence(Status *, 
			     INT8VectorSequence **);
void SDestroyVectorSequence(Status *, 
			     REAL4VectorSequence **);
void DDestroyVectorSequence(Status *, 
			     REAL8VectorSequence **);
void CDestroyVectorSequence(Status *, 
			     COMPLEX8VectorSequence **);
void ZDestroyVectorSequence(Status *, 
			     COMPLEX16VectorSequence **);


#endif

