#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/Inspiral.h>
#include <lal/FindRoot.h>

NRCSID (INSPIRALVELOCITYC, "$Id$"); 

void LALInspiralVelocity(LALStatus *status,
		      REAL8 *v,
		      TofVIn *ak)
{
  DFindRootIn rootIn;
  void *funcParams;


  INITSTATUS (status, "LALInspiralVelocity", INSPIRALVELOCITYC);
  ATTATCHSTATUSPTR(status);

  ASSERT (v, status, INSPIRALVELOCITY_ENULL, INSPIRALVELOCITY_MSGENULL);
  ASSERT (ak, status, INSPIRALVELOCITY_ENULL, INSPIRALVELOCITY_MSGENULL);

  rootIn.function = LALTofV;
  rootIn.xmax = ak->vlso;
  rootIn.xmin = ak->v0/2.;
  rootIn.xacc = 1.0e-8;



  funcParams = (void *) ak;


     if (ak->t==ak->t0) {
     *v = ak->v0;
  DETATCHSTATUSPTR(status);
  RETURN(status);
  }


  LALDBisectionFindRoot(status->statusPtr, v, &rootIn, funcParams);
  CHECKSTATUSPTR(status);




  DETATCHSTATUSPTR(status);
  RETURN(status);


}
