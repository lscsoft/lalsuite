#include <math.h>
#include "LALStdlib.h"
#include "Inspiral.h"
#include "FindRoot.h"

NRCSID (INSPIRALVELOCITYC, "$Id$"); 

void InspiralVelocity(Status *status,
		      REAL8 *v,
		      TofVIn *ak)
{
  DFindRootIn rootIn;
  void *funcParams;


  INITSTATUS (status, "InspiralVelocity", INSPIRALVELOCITYC);
  ATTATCHSTATUSPTR(status);

  ASSERT (v, status, INSPIRALVELOCITY_ENULL, INSPIRALVELOCITY_MSGENULL);
  ASSERT (ak, status, INSPIRALVELOCITY_ENULL, INSPIRALVELOCITY_MSGENULL);

  rootIn.function = TofV;
  rootIn.xmax = ak->vlso;
  rootIn.xmin = ak->v0/2.;
  rootIn.xacc = 1.0e-8;



  funcParams = (void *) ak;


     if (ak->t==ak->t0) {
     *v = ak->v0;
  DETATCHSTATUSPTR(status);
  RETURN(status);
  }


  DBisectionFindRoot(status->statusPtr, v, &rootIn, funcParams);
  CHECKSTATUSPTR(status);




  DETATCHSTATUSPTR(status);
  RETURN(status);


}
