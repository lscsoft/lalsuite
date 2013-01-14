#ifndef _LALPSIRD_H
#define _LALPSIRD_H

#include <math.h>

#include <lal/Units.h>
#include <lal/LALInspiral.h>
#include <lal/SeqFactories.h>

#include <lal/LALAdaptiveRungeKutta4.h>

#ifdef  __cplusplus
extern "C" {
#pragma }
#endif

/* use error codes above 1024 to avoid conflicts with GSL */
#define LALPSIRDPN_TEST_ENERGY		1025
#define LALPSIRDPN_TEST_OMEGADOT	1026
#define LALPSIRDPN_TEST_OMEGANAN	1028
#define LALPSIRDPN_TEST_OMEGAMATCH      1029
#define LALPSIRDPN_TEST_OMEGANONPOS     1031
#define LALPSIRDPN_TEST_OMEGACUT        1032

#ifdef  __cplusplus
#pragma {
}
#endif

#endif /* _LALPSPININSPIRALRD_H */
