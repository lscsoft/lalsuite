#include <lal/DopplerScan.h>
#include <lal/LALInitBarycenter.h>
#include <lal/UserInput.h>
#include <lal/LALString.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_rng.h>
#include "../TwoSpectSpecFunc.h"

typedef struct
{
   BOOLEAN help;		/**< Print this help/usage message */
   REAL8 Tsft;
   REAL8 SFToverlap;
   REAL8 t0;
   REAL8 Tobs;
   REAL8 cosi;
   REAL8 psi;
   REAL8 alpha;
   REAL8 delta;
   INT4 skylocations;
   LALStringVector *IFO;
   CHAR *outfilename;
   CHAR *ephemEarth;
   CHAR *ephemSun;
} UserVariables_t;

INT4 InitUserVars(UserVariables_t *uvar, int argc, char *argv[]);

int main(int argc, char *argv[])
{
   UserVariables_t XLAL_INIT_DECL(uvar);
   XLAL_CHECK ( InitUserVars(&uvar, argc, argv) == XLAL_SUCCESS, XLAL_EFUNC );
   
   MultiLALDetector *detectors = NULL;
   XLAL_CHECK( (detectors = XLALMalloc(sizeof(MultiLALDetector))) != NULL, XLAL_ENOMEM );
   detectors->length = uvar.IFO->length;
   for (UINT4 ii=0; ii<detectors->length; ii++) {
      if (strcmp("H1", uvar.IFO->data[ii])==0) {
         detectors->sites[ii] = lalCachedDetectors[LAL_LHO_4K_DETECTOR]; //H1
      } else if (strcmp("L1",uvar.IFO->data[ii])==0) {
         detectors->sites[ii] = lalCachedDetectors[LAL_LLO_4K_DETECTOR]; //L1
      } else if (strcmp("V1", uvar.IFO->data[ii])==0) {
         detectors->sites[ii] = lalCachedDetectors[LAL_VIRGO_DETECTOR];  //V1
      } else if (strcmp("H2", uvar.IFO->data[ii])==0) {
         detectors->sites[ii] = lalCachedDetectors[LAL_LHO_2K_DETECTOR]; //H2
      } else if (strcmp("H2r", uvar.IFO->data[ii])==0) {
         LALDetector H2 = lalCachedDetectors[LAL_LHO_2K_DETECTOR]; //H2 rotated
         H2.frDetector.xArmAzimuthRadians -= 0.25*LAL_PI;
         H2.frDetector.yArmAzimuthRadians -= 0.25*LAL_PI;
         memset(&(H2.frDetector.name), 0, sizeof(CHAR)*LALNameLength);
         snprintf(H2.frDetector.name, LALNameLength, "%s", "LHO_2k_rotatedPiOver4");
         XLAL_CHECK( (XLALCreateDetector(&(detectors->sites[ii]), &(H2.frDetector), LALDETECTORTYPE_IFODIFF)) != NULL, XLAL_EFUNC );
      } else {
         XLAL_ERROR(XLAL_EINVAL, "Not using valid interferometer! Expected 'H1', 'H2', 'H2r' (rotated H2), 'L1', or 'V1' not %s.\n", uvar.IFO->data[ii]);
      }
   }

   EphemerisData *edat = NULL;
   XLAL_CHECK( (edat = XLALInitBarycenter(uvar.ephemEarth, uvar.ephemSun)) != NULL, XLAL_EFUNC );

   LIGOTimeGPS tStart;
   XLALGPSSetREAL8 ( &tStart, uvar.t0 );
   XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC, "XLALGPSSetREAL8 failed\n" );

   MultiLIGOTimeGPSVector *multiTimestamps = NULL;
   XLAL_CHECK( (multiTimestamps = XLALMakeMultiTimestamps(tStart, uvar.Tobs, uvar.Tsft, uvar.SFToverlap, detectors->length)) != NULL, XLAL_EFUNC );

   LIGOTimeGPS refTime = multiTimestamps->data[0]->data[0];
   
   MultiDetectorStateSeries *multiStateSeries = NULL;
   XLAL_CHECK( (multiStateSeries = XLALGetMultiDetectorStates(multiTimestamps, detectors, edat, uvar.SFToverlap)) != NULL, XLAL_EFUNC );

   gsl_rng *rng = NULL;
   XLAL_CHECK( (rng = gsl_rng_alloc(gsl_rng_mt19937)) != NULL, XLAL_EFUNC );
   gsl_rng_set(rng, 0);

   FILE *OUTPUT;
   XLAL_CHECK( (OUTPUT = fopen(uvar.outfilename,"w")) != NULL, XLAL_EIO, "Output file %s could not be opened\n", uvar.outfilename );

   for (INT4 n=0; n<uvar.skylocations; n++) {
      SkyPosition skypos;
      if (XLALUserVarWasSet(&(uvar.alpha)) && XLALUserVarWasSet(&(uvar.delta)) && n==0) {
         skypos.longitude = uvar.alpha;
         skypos.latitude = uvar.delta;
         skypos.system = COORDINATESYSTEM_EQUATORIAL;
      } else {
         skypos.longitude = LAL_TWOPI*gsl_rng_uniform(rng);
         skypos.latitude = LAL_PI*gsl_rng_uniform(rng) - LAL_PI_2;
         skypos.system = COORDINATESYSTEM_EQUATORIAL;
      }

      REAL8 cosi0, psi0;
      if (XLALUserVarWasSet(&(uvar.cosi)) && n==0) cosi0 = uvar.cosi;
      else cosi0 = 2.0*gsl_rng_uniform(rng) - 1.0;
      if (XLALUserVarWasSet(&(uvar.psi)) && n==0) psi0 = uvar.psi;
      else psi0 = LAL_PI*gsl_rng_uniform(rng);

      MultiAMCoeffs *multiAMcoefficients = NULL;
      XLAL_CHECK( (multiAMcoefficients = XLALComputeMultiAMCoeffs(multiStateSeries, NULL, skypos)) != NULL, XLAL_EFUNC );

      MultiSSBtimes *multissb = NULL;
      XLAL_CHECK( (multissb = XLALGetMultiSSBtimes(multiStateSeries, skypos, refTime, SSBPREC_RELATIVISTICOPT)) != NULL, XLAL_EFUNC );

      REAL8 frequency0 = 1000.0 + (gsl_rng_uniform(rng)-0.5)/uvar.Tsft;
      REAL8 frequency = 1000.0;

      for (UINT4 ii=0; ii<multiAMcoefficients->data[0]->a->length; ii++) {
         REAL4 Fplus0 = multiAMcoefficients->data[0]->a->data[ii]*cos(2.0*psi0) + multiAMcoefficients->data[0]->b->data[ii]*sin(2.0*psi0);
         REAL4 Fcross0 = multiAMcoefficients->data[0]->b->data[ii]*cos(2.0*psi0) - multiAMcoefficients->data[0]->a->data[ii]*sin(2.0*psi0);
         REAL4 Fplus1 = multiAMcoefficients->data[1]->a->data[ii]*cos(2.0*psi0) + multiAMcoefficients->data[1]->b->data[ii]*sin(2.0*psi0);
         REAL4 Fcross1 = multiAMcoefficients->data[1]->b->data[ii]*cos(2.0*psi0) - multiAMcoefficients->data[1]->a->data[ii]*sin(2.0*psi0);
         COMPLEX16 RatioTerm0 = crect(0.5*Fplus1*(1.0+cosi0*cosi0), Fcross1*cosi0)/crect(0.5*Fplus0*(1.0+cosi0*cosi0), Fcross0*cosi0);

         REAL4 detPhaseArg = 0.0, detPhaseMag = 0.0;
         BOOLEAN loopbroken = 0;
         for (INT4 jj=0; jj<50 && !loopbroken; jj++) {
            REAL4 psi = 0.02*jj*LAL_PI;
            Fplus0 = multiAMcoefficients->data[0]->a->data[ii]*cos(2.0*psi) + multiAMcoefficients->data[0]->b->data[ii]*sin(2.0*psi);
            Fcross0 = multiAMcoefficients->data[0]->b->data[ii]*cos(2.0*psi) - multiAMcoefficients->data[0]->a->data[ii]*sin(2.0*psi);
            Fplus1 = multiAMcoefficients->data[1]->a->data[ii]*cos(2.0*psi) + multiAMcoefficients->data[1]->b->data[ii]*sin(2.0*psi);
            Fcross1 = multiAMcoefficients->data[1]->b->data[ii]*cos(2.0*psi) - multiAMcoefficients->data[1]->a->data[ii]*sin(2.0*psi);
            for (INT4 kk=0; kk<51 && !loopbroken; kk++) {
               //REAL4 cosi = 1.0 - 2.0*0.02*kk;
               REAL4 cosi;
               if (cosi0<0.0) cosi = -0.02*kk;
               else cosi = 0.02*kk;
               COMPLEX16 complexnumerator = crect(0.5*Fplus1*(1.0+cosi*cosi), Fcross1*cosi);
               COMPLEX16 complexdenominator = crect(0.5*Fplus0*(1.0+cosi*cosi) , Fcross0*cosi);
               if (cabs(complexdenominator)>1.0e-7) {
                  COMPLEX16 complexval = complexnumerator/complexdenominator;
                  detPhaseMag += fmin(cabs(complexval), 10.0);
                  detPhaseArg += gsl_sf_angle_restrict_pos(carg(complexval));
               } else {
                  loopbroken = 1;
                  detPhaseMag = 0.0;
                  detPhaseArg = 0.0;
               }
            }
         }
         detPhaseMag /= 2550.0;
         detPhaseArg /= 2550.0;

         REAL8 timediff0 = multissb->data[0]->DeltaT->data[ii] - 0.5*uvar.Tsft*multissb->data[0]->Tdot->data[ii];
         REAL8 timediff1 = multissb->data[1]->DeltaT->data[ii] - 0.5*uvar.Tsft*multissb->data[1]->Tdot->data[ii];
         REAL8 tau = timediff1 - timediff0;
         REAL8 freqshift0 = -LAL_TWOPI*tau*frequency0;
         REAL8 freqshift = -LAL_TWOPI*tau*frequency;

         REAL8 nearestFrequency = round(multissb->data[0]->Tdot->data[ii]*frequency*uvar.Tsft)/uvar.Tsft;

         COMPLEX16 dirichlet0 = conj(DirichletKernelLargeNHann((multissb->data[1]->Tdot->data[ii]*frequency0-nearestFrequency)*uvar.Tsft)/DirichletKernelLargeNHann((multissb->data[0]->Tdot->data[ii]*frequency0-nearestFrequency)*uvar.Tsft));
         COMPLEX16 dirichlet = conj(DirichletKernelLargeNHann((multissb->data[1]->Tdot->data[ii]*frequency-nearestFrequency)*uvar.Tsft)/DirichletKernelLargeNHann((multissb->data[0]->Tdot->data[ii]*frequency-nearestFrequency)*uvar.Tsft));

         COMPLEX16 signal0 = 0.5*crect(0.5*Fplus0*(1.0+cosi0*cosi0), Fcross0*cosi0)*cpolar(1.0, LAL_TWOPI*frequency0*0.5*uvar.Tsft+LAL_TWOPI*frequency0*(multissb->data[0]->DeltaT->data[ii]+multissb->data[0]->Tdot->data[ii]*0.5*uvar.Tsft))*uvar.Tsft*DirichletKernelLargeNHann((multissb->data[0]->Tdot->data[ii]*frequency0-nearestFrequency)*uvar.Tsft);
         COMPLEX16 signal1 = 0.5*crect(0.5*Fplus1*(1.0+cosi0*cosi0), Fcross1*cosi0)*cpolar(1.0, LAL_TWOPI*frequency0*0.5*uvar.Tsft+LAL_TWOPI*frequency0*(multissb->data[1]->DeltaT->data[ii]+multissb->data[1]->Tdot->data[ii]*0.5*uvar.Tsft))*uvar.Tsft*DirichletKernelLargeNHann((multissb->data[1]->Tdot->data[ii]*frequency0-nearestFrequency)*uvar.Tsft);
      
         fprintf(OUTPUT, "%g %g %g %g %g %g %g %g %g %g %g %g %g %g\n", cabs(signal0), gsl_sf_angle_restrict_pos(carg(signal0)), cabs(signal1), gsl_sf_angle_restrict_pos(carg(signal1)), cabs(RatioTerm0), gsl_sf_angle_restrict_pos(carg(RatioTerm0)), detPhaseMag, detPhaseArg, freqshift0, freqshift, cabs(dirichlet0), gsl_sf_angle_restrict_pos(carg(dirichlet0)), cabs(dirichlet), gsl_sf_angle_restrict_pos(carg(dirichlet)));

         //fprintf(OUTPUT, "%g %g %g %g\n", cabs(signal0), gsl_sf_angle_restrict_pos(carg(signal0)), cabs(signal1), gsl_sf_angle_restrict_pos(carg(signal1)));
      }

      XLALDestroyMultiAMCoeffs(multiAMcoefficients);
      XLALDestroyMultiSSBtimes(multissb);
   }

   fclose(OUTPUT);
   gsl_rng_free(rng);
   XLALDestroyMultiDetectorStateSeries(multiStateSeries);
   XLALDestroyMultiTimestamps(multiTimestamps);
   XLALDestroyEphemerisData(edat);
   XLALFree(detectors);
   XLALDestroyUserVars();
}

INT4 InitUserVars(UserVariables_t *uvar, int argc, char *argv[])
{
   XLAL_CHECK ( uvar != NULL, XLAL_EINVAL, "Invalid NULL input 'uvar'\n");
   XLAL_CHECK ( argv != NULL, XLAL_EINVAL, "Invalid NULL input 'argv'\n");

   uvar->ephemEarth = XLALStringDuplicate("earth00-19-DE405.dat.gz");
   uvar->ephemSun = XLALStringDuplicate("sun00-19-DE405.dat.gz");
   uvar->outfilename = XLALStringDuplicate("output.dat");
   uvar->Tsft = 1800;
   uvar->SFToverlap = 900;
   uvar->skylocations = 1;

   XLALregBOOLUserStruct(  help,        'h', UVAR_HELP     , "Print this help/usage message");
   XLALregREALUserStruct(  Tsft,         0 , UVAR_OPTIONAL , "SFT coherence time");
   XLALregREALUserStruct(  SFToverlap,   0 , UVAR_OPTIONAL , "SFT overlap in seconds, usually Tsft/2");
   XLALregREALUserStruct(  t0,           0 , UVAR_OPTIONAL , "GPS start time of the search");
   XLALregREALUserStruct(  Tobs,         0 , UVAR_OPTIONAL , "Duration of the search (in seconds)");
   XLALregREALUserStruct(  cosi,         0 , UVAR_OPTIONAL , "Cosine of NS inclinaiont angle");
   XLALregREALUserStruct(  psi,          0 , UVAR_OPTIONAL , "Polarization angle of GW");
   XLALregREALUserStruct(  alpha,        0 , UVAR_OPTIONAL , "Right ascension of source (in radians)");
   XLALregREALUserStruct(  delta,        0 , UVAR_OPTIONAL , "Declination of source (in radians)");
   XLALregINTUserStruct(   skylocations, 0 , UVAR_OPTIONAL , "Number of sky locations");
   XLALregLISTUserStruct(  IFO,          0 , UVAR_REQUIRED , "CSV list of detectors, eg. \"H1,H2,L1,G1, ...\" ");
   XLALregSTRINGUserStruct(outfilename,  0 , UVAR_OPTIONAL , "Output filename");
   XLALregSTRINGUserStruct(ephemEarth,   0 , UVAR_OPTIONAL , "Earth ephemeris file");
   XLALregSTRINGUserStruct(ephemSun,     0 , UVAR_OPTIONAL , "Sun ephemeris file");

   XLAL_CHECK( XLALUserVarReadAllInput(argc, argv) == XLAL_SUCCESS, XLAL_EFUNC );

   if ( uvar->help ) exit (0);

   return XLAL_SUCCESS;
}
