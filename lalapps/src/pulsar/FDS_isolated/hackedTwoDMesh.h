
#ifndef _HACKEDTWODMESH_H
#define _HACKEDTWODMESH_H

#include <lal/LALStdlib.h>
#include <lal/TwoDMesh.h>

#ifdef __cplusplus
extern "C" {
#pragma }
#endif

#include <lal/LALStdlib.h>

NRCSID(HTWODMESHH,"$Id$");


typedef struct h_tagTwoDMeshNode {
  REAL8 x, y;
  REAL8 dx;
  REAL8 dy[2];
  struct h_tagTwoDMeshNode *next;
  struct h_tagTwoDMeshNode *subMesh;
  UINT4 nSub;
} h_TwoDMeshNode;


typedef struct  {
  REAL8 domain[2];
  void (*getRange)( LALStatus *, REAL8 [2], REAL8, void *);
  void *rangeParams;
  void (*getMetric)( LALStatus *, REAL8 [3], REAL8 [2], void *);
  void *metricParams;
  REAL8 mThresh;
  REAL8 widthMaxFac;
  REAL8 widthRetryFac;
  UINT4 maxColumns;
  UINT4 nIn;
  UINT4 nOut;
} h_TwoDMeshParamStruc;

typedef struct  {
  REAL8 domain[2];
  REAL8 leftRange[2];  
  REAL8 rightRange[2];
  REAL8 leftClip[2];
  REAL8 rightClip[2];
  BOOLEAN tooWide;
} h_TwoDColumnParamStruc;


void
hackedLALCreateTwoDMesh( LALStatus          *stat,
			 h_TwoDMeshNode       **mesh,
			 h_TwoDMeshParamStruc *params );

void
hackedLALTwoDMesh( LALStatus          *stat,
		   h_TwoDMeshNode       **tail,
		   h_TwoDMeshParamStruc *params );

void
hackedLALTwoDColumn( LALStatus            *stat,
		     h_TwoDMeshNode         **tail,
		     h_TwoDColumnParamStruc *column,
		     h_TwoDMeshParamStruc   *params );

void
hackedLALDestroyTwoDMesh( LALStatus    *stat,
			  h_TwoDMeshNode **mesh,
			  UINT4        *nFree );



#ifdef __cplusplus
#pragma {
}
#endif

#endif /* _TWODMESH_H */
