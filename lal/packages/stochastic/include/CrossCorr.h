/*-----------------------------------------------------------------------
 *
 * File Name:  CrossCorr.h
 *
 * Author: Steve Drasco
 *
 * Revision: $Id$
 *
 *-----------------------------------------------------------------------
 *
 * NAME
 * CrossCorr.h
 *
 * SYNOPSIS
 * #include <lal/CcrossCorr.h>
 *
 * DESCRIPTION
 * Error codes, typedefs, and protypes for the function LALCrossCorr()
 *
 * DIAGNOSTICS
 *
 *-----------------------------------------------------------------------
 */

#ifndef _CROSSCORR_H
#define _CROSSCORR_H

#include <lal/LALStdlib.h>

#ifdef  __cplusplus
extern "C" {
#endif

NRCSID (CROSSCORRH, "$Id$");

typedef struct tagCCIn {		/* input structure                  		*/
	REAL4TimeSeries	*g1;		/* output of detector 1				*/
	REAL4TimeSeries	*g2;		/* output of detector 2				*/
	REAL4Vector	*QmaxTilde;	/* optimized kernel                 		*/
	INT2		plan;		/* 1 -> estimate plan, otherwise measure plan 	*/
} CCIn;

#define CROSSCORR_EIN		1
#define CROSSCORR_EOUT		2
#define CROSSCORR_ENULL		4
#define CROSSCORR_ESIZE1	8
#define CROSSCORR_ESIZE2	16


#define CROSSCORR_MSGEIN	"Pointer to input structure must be non-NULL"
#define CROSSCORR_MSGEOUT	"Pointer to output value must be non-NULL"
#define CROSSCORR_MSGENULL      "The pointers to the detector output vectors must be non-NULL"
#define CROSSCORR_MSGESIZE1	"Detector output vectors must have equal length"
#define CROSSCORR_MSGESIZE2	"Kernel vector and detector output vectors must have equal length"

#ifdef  __cplusplus
}
#endif

#endif
