
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


void
hackedLALCreateTwoDMesh( LALStatus          *stat,
			 TwoDMeshNode       **mesh,
			 TwoDMeshParamStruc *params );

void
hackedLALTwoDMesh( LALStatus          *stat,
		   TwoDMeshNode       **tail,
		   TwoDMeshParamStruc *params );

void
hackedLALTwoDColumn( LALStatus            *stat,
		     TwoDMeshNode         **tail,
		     TwoDColumnParamStruc *column,
		     TwoDMeshParamStruc   *params );



#ifdef __cplusplus
#pragma {
}
#endif

#endif /* _TWODMESH_H */
