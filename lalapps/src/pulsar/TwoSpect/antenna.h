
#ifndef __ANTENNA_H__
#define __ANTENNA_H__

#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/DetResponse.h>
#include <lal/LALInitBarycenter.h>
#include <lal/Velocity.h>

EphemerisData * new_Ephemeris(CHAR *earth_ephemeris, CHAR *sun_ephemeris);
void free_Ephemeris(EphemerisData *ephemdata);
void initEphemeris(EphemerisData *ephemdata);

INT4Vector * CompBinShifts(REAL4 freq, REAL4Vector *velocities, REAL4 Tcoh, REAL4 dopplerMultiplier);

REAL4Vector * CompAntennaPatternWeights(REAL4 ra, REAL4 dec, REAL8 t0, REAL4 Tcoh, REAL8 Tobs, LALDetector det);
REAL4Vector * CompAntennaVelocity(REAL4 ra, REAL4 dec, REAL8 t0, REAL4 Tcoh, REAL8 Tobs, LALDetector det, EphemerisData *edat);

#endif

