/********************************** <lalVerbatim file="PowerHV">
Author: Brady, P
$Id$
**************************************************** </lalVerbatim> */

#ifndef _POWER_H
#define _POWER_H

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <unistd.h>
#include <lal/EPData.h>
#include <lal/ExcessPower.h>

/* turn these on to enable debugging tools */
#undef POWER_SO_DEBUG              /* print debug information               */
#undef POWER_SO_ATTACH_DDD         /* sleep for a few seconds to attach ddd */
#undef POWER_SO_GRAPHING_FUNCS     /* enable some printing tools            */
#define IFO_STRING_LENGTH 2

#ifdef  __cplusplus
extern "C" {
#endif

NRCSID (POWERH, "power $Id$");

#ifdef  __cplusplus
}
#endif

#endif /* _POWER_H  */




