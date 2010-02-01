

#ifndef __IHS_H__
#define __IHS_H__

#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include "TwoSpect.h"

typedef struct
{
   REAL4Vector *maxima;
   INT4Vector *locations;
   INT4 columns;
} ihsMaximaStruct;

typedef struct
{
   REAL4 ihs;
   INT4 loc;
} ihsVals;

typedef struct
{
   REAL4Vector *ihsfar;
   REAL4Vector *ihsdistMean;
   REAL4Vector *ihsdistSigma;
} ihsfarStruct;


ihsMaximaStruct *new_ihsMaxima(ffdataStruct *ffdata, INT4 columns);
void free_ihsMaxima(ihsMaximaStruct *data);
void runIHS(ihsMaximaStruct *out, ffdataStruct *in, INT4 columns);

ihsVals * new_ihsVals(void);
void free_ihsVals(ihsVals *ihsvals);
void incHarmSum(ihsVals *out, REAL4Vector *in);

ihsfarStruct * new_ihsfarStruct(INT4 columns);
void free_ihsfarStruct(ihsfarStruct *ihsfarstruct);
void genIhsFar(ihsfarStruct *out, ffdataStruct *ffdata, INT4 columns, REAL4 threshold);

REAL4Vector * ihsSums(REAL4Vector *ihss, INT4 cols);

REAL4 ihsFOM(REAL4Vector *ihss, INT4Vector *locs, REAL4Vector *expect);
REAL4 ihsLoc(REAL4Vector *ihss, INT4Vector *locs, REAL4Vector *expect);


#endif



