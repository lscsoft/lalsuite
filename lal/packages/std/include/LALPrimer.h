/* <lalVerbatim> */
#ifndef _LALPRIMER_H  /* Double-include protection. */
#define _LALPRIMER_H

#include "LALDatatypes.h"

#ifdef  __cplusplus   /* C++ protection. */
extern "C" {
#endif


NRCSID( LALPRIMERH, "$Id$" );

#define LALPRIMERH_ENULL 1
#define LALPRIMERH_EDIV0 2
#define LALPRIMERH_MSGENULL "Null pointer"
#define LALPRIMERH_MSGEDIV0 "Division by zero"

void
REAL4Invert( Status *stat, REAL4 *output, REAL4 input );

void
REAL4Divide( Status *stat, REAL4 *output, REAL4 numer, REAL4 denom);


#ifdef  __cplusplus
}
#endif  /* C++ protection. */
#endif  /* Double-include protection. */
/* </lalVerbatim> */
