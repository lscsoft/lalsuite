/*----------------------------------------------------------------------- 
 * 
 * File Name: FindChirpSP.h
 *
 * Author: Brown, D. A.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#ifndef _FINDCHIRPSP_H
#define _FINDCHIRPSP_H

#include <lal/LALDatatypes.h>
#include <lal/RealFFT.h>
#include <lal/DataBuffer.h>
#include <lal/LALInspiral.h>
#include <lal/FindChirp.h>

#ifdef  __cplusplus
extern "C" {
#endif

#define FINDCHIRPSP_ENULL 1
#define FINDCHIRPSP_ENNUL 2
#define FINDCHIRPSP_EALOC 3
#define FINDCHIRPSP_ENUMZ 4
#define FINDCHIRPSP_ESEGZ 5
#define FINDCHIRPSP_EMISM 6
#define FINDCHIRPSP_EDELT 7
#define FINDCHIRPSP_EFLOW 8
#define FINDCHIRPSP_EDYNR 9
  
#define FINDCHIRPSP_MSGENULL "Null pointer"
#define FINDCHIRPSP_MSGENNUL "Non-null pointer"
#define FINDCHIRPSP_MSGEALOC "Memory allocation error"
#define FINDCHIRPSP_MSGENUMZ "Invalid number of segments"
#define FINDCHIRPSP_MSGESEGZ "Invalid number of points in segments"
#define FINDCHIRPSP_MSGEMISM "Mismatch between number of points in segments"
#define FINDCHIRPSP_MSGEDELT "deltaT is zero or negative"
#define FINDCHIRPSP_MSGEFLOW "Low frequency cutoff is negative"
#define FINDCHIRPSP_MSGEDYNR "Dynamic range scaling is zero or negative"



typedef struct
tagFindChirpSPDataParams
{
  REAL4Vector                  *ampVec;
  RealFFTPlan                  *fwdPlan;
  RealFFTPlan                  *invPlan;
  REAL4Vector                  *vVec;
  REAL4Vector                  *wVec;
  COMPLEX8Vector               *wtildeVec;
  REAL4Vector                  *tmpltPowerVec;
  REAL4                         deltaT;
  REAL4                         fLow;
  REAL4                         dynRange;
  UINT4                         invSpecTrunc;
}
FindChirpSPDataParams;

typedef struct
tagFindChirpSPTmpltParams
{
  REAL4                         deltaT;
  REAL4                         fLow;
  REAL4                         dynRange;
  REAL4Vector                  *xfacVec;
}
FindChirpSPTmpltParams;


void
LALFindChirpSPDataInit (
    LALStatus                  *status,
    FindChirpSPDataParams     **output,
    FindChirpInitParams        *params
    );

void
LALFindChirpSPData (
    LALStatus                  *status,
    FindChirpSegmentVector     *fcSegVec,
    DataSegmentVector          *dataSegVec,
    FindChirpSPDataParams      *params
    );

void
LALFindChirpSPDataFinalize (
    LALStatus                  *status,
    FindChirpSPDataParams     **output
    );


void
LALFindChirpSPTemplateInit (
    LALStatus                  *status,
    FindChirpSPTmpltParams    **output,
    FindChirpInitParams        *params
    );

void
LALFindChirpSPTemplate (
    LALStatus                  *status,
    FindChirpTemplate          *fcTmplt,
    InspiralTemplate           *tmplt,
    FindChirpSPTmpltParams     *params
    );

void
LALFindChirpSPTemplateFinalize (
    LALStatus                  *status,
    FindChirpSPTmpltParams    **output
    );


#ifdef  __cplusplus
}
#endif

#endif /* _FINDCHIRPSP_H */
