#include <lal/DopplerScan.h>
#include <lal/LALInitBarycenter.h>
#include <lal/UserInput.h>
#include <lal/LALString.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_rng.h>
#include "../TwoSpectSpecFunc.h"

typedef struct
{
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
   BOOLEAN unrestrictedCosi;
   BOOLEAN rectWindow;
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

      REAL8 frequency = 1000.0;
      REAL8 frequency0 = frequency + (gsl_rng_uniform(rng)-0.5)/uvar.Tsft;

      for (UINT4 ii=0; ii<multiAMcoefficients->data[0]->a->length; ii++) {
         REAL4 Fplus0 = multiAMcoefficients->data[0]->a->data[ii]*cos(2.0*psi0) + multiAMcoefficients->data[0]->b->data[ii]*sin(2.0*psi0);
         REAL4 Fcross0 = multiAMcoefficients->data[0]->b->data[ii]*cos(2.0*psi0) - multiAMcoefficients->data[0]->a->data[ii]*sin(2.0*psi0);
         REAL4 Fplus1 = multiAMcoefficients->data[1]->a->data[ii]*cos(2.0*psi0) + multiAMcoefficients->data[1]->b->data[ii]*sin(2.0*psi0);
         REAL4 Fcross1 = multiAMcoefficients->data[1]->b->data[ii]*cos(2.0*psi0) - multiAMcoefficients->data[1]->a->data[ii]*sin(2.0*psi0);
         COMPLEX16 RatioTerm0 = crect(0.5*Fplus1*(1.0+cosi0*cosi0), Fcross1*cosi0)/crect(0.5*Fplus0*(1.0+cosi0*cosi0), Fcross0*cosi0);  //real det-sig ratio term

         REAL4 detPhaseArg = 0.0, detPhaseMag = 0.0;
         BOOLEAN loopbroken = 0;
         for (INT4 jj=0; jj<16 && !loopbroken; jj++) {
            REAL4 psi = 0.0625*jj*LAL_PI;
            Fplus0 = multiAMcoefficients->data[0]->a->data[ii]*cos(2.0*psi) + multiAMcoefficients->data[0]->b->data[ii]*sin(2.0*psi);
            Fcross0 = multiAMcoefficients->data[0]->b->data[ii]*cos(2.0*psi) - multiAMcoefficients->data[0]->a->data[ii]*sin(2.0*psi);
            Fplus1 = multiAMcoefficients->data[1]->a->data[ii]*cos(2.0*psi) + multiAMcoefficients->data[1]->b->data[ii]*sin(2.0*psi);
            Fcross1 = multiAMcoefficients->data[1]->b->data[ii]*cos(2.0*psi) - multiAMcoefficients->data[1]->a->data[ii]*sin(2.0*psi);
            for (INT4 kk=0; kk<21 && !loopbroken; kk++) {
               REAL4 cosi = 1.0 - 2.0*0.05*kk;
               if (!uvar.unrestrictedCosi) {
                  if (cosi0<0.0) cosi = -0.05*kk;
                  else cosi = 0.05*kk;
               }
               COMPLEX16 complexnumerator = crect(0.5*Fplus1*(1.0+cosi*cosi), Fcross1*cosi);
               COMPLEX16 complexdenominator = crect(0.5*Fplus0*(1.0+cosi*cosi) , Fcross0*cosi);
               if (cabs(complexdenominator)>1.0e-6) {
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
         detPhaseMag /= 336.0;
         detPhaseArg /= 336.0;
         COMPLEX16 RatioTerm = cpolar(detPhaseMag, detPhaseArg);

         //Bin of interest
         REAL8 signalFrequencyBin = round(multissb->data[0]->Tdot->data[ii]*frequency0*uvar.Tsft) - frequency*uvar.Tsft;  //estimated nearest freq in ref IFO

         REAL8 timediff0 = multissb->data[0]->DeltaT->data[ii] - 0.5*uvar.Tsft*multissb->data[0]->Tdot->data[ii];
         REAL8 timediff1 = multissb->data[1]->DeltaT->data[ii] - 0.5*uvar.Tsft*multissb->data[1]->Tdot->data[ii];
         REAL8 tau = timediff1 - timediff0;
         REAL8 freqshift0 = -LAL_TWOPI*tau*frequency0;  //real freq shift
         REAL8 freqshift = -LAL_TWOPI*tau*(round(multissb->data[0]->Tdot->data[ii]*frequency0*uvar.Tsft)/uvar.Tsft);    //estimated freq shift
         COMPLEX16 phaseshift0 = cpolar(1.0, freqshift0);
         COMPLEX16 phaseshift = cpolar(1.0, freqshift);

         REAL8 delta0_0 = (multissb->data[0]->Tdot->data[ii]*frequency0-frequency)*uvar.Tsft - signalFrequencyBin;
         REAL8 delta1_0 = (multissb->data[1]->Tdot->data[ii]*frequency0-frequency)*uvar.Tsft - signalFrequencyBin;
         REAL8 realSigBinDiff = round(delta0_0) - round(delta1_0);
         delta1_0 += realSigBinDiff;
         REAL8 delta0 = round(multissb->data[0]->Tdot->data[ii]*frequency0*uvar.Tsft)*(multissb->data[0]->Tdot->data[ii] - 1.0);
         REAL8 delta1 = round(multissb->data[0]->Tdot->data[ii]*frequency0*uvar.Tsft)*(multissb->data[1]->Tdot->data[ii] - 1.0);
         REAL8 estSigBinDiff = round(delta0) - round(delta1);
         delta1 += estSigBinDiff;
         COMPLEX16 dirichlet0;
         if (!uvar.rectWindow) {
            if (fabsf((REAL4)delta1_0)<(REAL4)1.0e-6) {
               if (fabsf((REAL4)delta0_0)<(REAL4)1.0e-6) dirichlet0 = crect(1.0, 0.0);
               else if (fabsf((REAL4)(delta0_0*delta0_0-1.0))<(REAL4)1.0e-6) dirichlet0 = crect(-2.0, 0.0);
               else if (fabsf((REAL4)(delta0_0-roundf(delta0_0)))<(REAL4)1.0e-6) {
                  dirichlet0 = crect(0.0, 0.0);
                  continue;
               }
	       else dirichlet0 = -LAL_PI*crectf(cos(LAL_PI*delta0_0), -sin(LAL_PI*delta0_0))*delta0_0*(delta0_0*delta0_0 - 1.0)/sin(LAL_PI*delta0_0);
            }
            else if (fabsf((REAL4)(delta1_0*delta1_0-1.0))<(REAL4)1.0e-6) {
               if (fabsf((REAL4)delta0_0)<(REAL4)1.0e-6) dirichlet0 = crect(-0.5, 0.0);
               else if (fabsf((REAL4)(delta0_0*delta0_0-1.0))<(REAL4)1.0e-6) dirichlet0 = crect(1.0, 0.0);
               else if (fabsf((REAL4)(delta0_0-roundf(delta0_0)))<(REAL4)1.0e-6) {
                  dirichlet0 = crect(0.0, 0.0);
                  continue;
               }
	       else dirichlet0 = -LAL_PI_2*crectf(-cos(LAL_PI*delta0_0), sin(LAL_PI*delta0_0))*delta0_0*(delta0_0*delta0_0 - 1.0)/sin(LAL_PI*delta0_0);
            }
            else if (fabsf((REAL4)delta0_0)<(REAL4)1.0e-6) dirichlet0 = -LAL_1_PI*crectf(cos(LAL_PI*delta1_0), sin(LAL_PI*delta1_0))*sin(LAL_PI*delta1_0)/(delta1_0*(delta1_0*delta1_0-1.0));
            else if (fabsf((REAL4)(delta0_0*delta0_0-1.0))<(REAL4)1.0e-6) dirichlet0 = LAL_2_PI*crectf(cos(LAL_PI*delta1_0), sin(LAL_PI*delta1_0))*sin(LAL_PI*delta1_0)/(delta1_0*(delta1_0*delta1_0-1.0));
            else if (fabsf((REAL4)(delta0_0-roundf(delta0_0)))<(REAL4)1.0e-6) {
               dirichlet0 = crect(0.0, 0.0);
               continue;
            }
            else dirichlet0 = sin(LAL_PI*delta1_0)/sin(LAL_PI*delta0_0)*(delta0_0*(delta0_0*delta0_0-1.0))/(delta1_0*(delta1_0*delta1_0-1.0))*cpolarf(1.0,LAL_PI*(delta1_0-delta0_0));
         } else {
            if (fabsf((REAL4)delta1_0)<(REAL4)1.0e-6) {
               if (fabsf((REAL4)delta0_0)<(REAL4)1.0e-6) dirichlet0 = crect(1.0, 0.0);
               else if (fabsf((REAL4)(delta0_0-roundf(delta0_0)))<(REAL4)1.0e-6) {
                  dirichlet0 = crect(0.0, 0.0);
                  continue;
               }
               else dirichlet0 = conj(1.0/((cpolar(1.0,LAL_TWOPI*delta0_0)-1.0)/crect(0.0,LAL_TWOPI*delta0_0)));
            }
            else if (fabsf((REAL4)delta0_0)<(REAL4)1.0e-6) dirichlet0 = conj((cpolar(1.0,LAL_TWOPI*delta1_0)-1.0)/crect(0.0,LAL_TWOPI*delta1_0));
            else if (fabsf((REAL4)(delta0_0-roundf(delta0_0)))<(REAL4)1.0e-6) {
               dirichlet0 = crect(0.0, 0.0);
               continue;
            }
            else dirichlet0 = conj(((cpolar(1.0,LAL_TWOPI*delta1_0)-1.0)/crect(0.0,LAL_TWOPI*delta1_0))/((cpolar(1.0,LAL_TWOPI*delta0_0)-1.0)/crect(0.0,LAL_TWOPI*delta0_0)));
         }

         COMPLEX16 dirichlet;
         if (!uvar.rectWindow) {
            if (fabsf((REAL4)delta1)<(REAL4)1.0e-6) {
               if (fabsf((REAL4)delta0)<(REAL4)1.0e-6) dirichlet = crect(1.0, 0.0);
               else if (fabsf((REAL4)(delta0*delta0-1.0))<(REAL4)1.0e-6) dirichlet = crect(-2.0, 0.0);
               else if (fabsf((REAL4)(delta0-roundf(delta0)))<(REAL4)1.0e-6) dirichlet = crect(0.0, 0.0);
               else dirichlet = -LAL_PI*crectf(cos(LAL_PI*delta0), -sin(LAL_PI*delta0))*delta0*(delta0*delta0 - 1.0)/sin(LAL_PI*delta0);
            }
            else if (fabsf((REAL4)(delta1*delta1-1.0))<(REAL4)1.0e-6) {
               if (fabsf((REAL4)delta0)<(REAL4)1.0e-6) dirichlet = crect(-0.5, 0.0);
               else if (fabsf((REAL4)(delta0*delta0-1.0))<(REAL4)1.0e-6) dirichlet = crect(1.0, 0.0);
               else if (fabsf((REAL4)(delta0-roundf(delta0)))<(REAL4)1.0e-6) dirichlet = crect(0.0, 0.0);
               else dirichlet = -LAL_PI_2*crectf(-cos(LAL_PI*delta0), sin(LAL_PI*delta0))*delta0*(delta0*delta0 - 1.0)/sin(LAL_PI*delta0);
            }
            else if (fabsf((REAL4)delta0)<(REAL4)1.0e-6) dirichlet = -LAL_1_PI*crectf(cos(LAL_PI*delta1), sin(LAL_PI*delta1))*sin(LAL_PI*delta1)/(delta1*(delta1*delta1-1.0));
            else if (fabsf((REAL4)(delta0*delta0-1.0))<(REAL4)1.0e-6) dirichlet = LAL_2_PI*crectf(cos(LAL_PI*delta1), sin(LAL_PI*delta1))*sin(LAL_PI*delta1)/(delta1*(delta1*delta1-1.0));
            else if (fabsf((REAL4)(delta0-roundf(delta0)))<(REAL4)1.0e-6) dirichlet = crect(0.0, 0.0);
            else dirichlet = sin(LAL_PI*delta1)/sin(LAL_PI*delta0)*(delta0*(delta0*delta0-1.0))/(delta1*(delta1*delta1-1.0))*cpolarf(1.0,LAL_PI*(delta1-delta0));
         } else {
            if (fabsf((REAL4)delta1)<(REAL4)1.0e-6) {
               if (fabsf((REAL4)delta0)<(REAL4)1.0e-6) dirichlet = crect(1.0, 0.0);
               else if (fabsf((REAL4)(delta0-roundf(delta0)))<(REAL4)1.0e-6) dirichlet = crect(0.0, 0.0);
               else dirichlet = conj(1.0/((cpolar(1.0,LAL_TWOPI*delta0)-1.0)/crect(0.0,LAL_TWOPI*delta0)));
            }
            else if (fabsf((REAL4)delta0)<(REAL4)1.0e-6) dirichlet = conj((cpolar(1.0,LAL_TWOPI*delta1)-1.0)/crect(0.0,LAL_TWOPI*delta1));
            else if (fabsf((REAL4)(delta0-roundf(delta0)))<(REAL4)1.0e-6) dirichlet = crect(0.0, 0.0);
            else dirichlet = conj(((cpolar(1.0,LAL_TWOPI*delta1)-1.0)/crect(0.0,LAL_TWOPI*delta1))/((cpolar(1.0,LAL_TWOPI*delta0)-1.0)/crect(0.0,LAL_TWOPI*delta0)));
         }
         dirichlet = cpolar(1.0, carg(dirichlet));
         if (fabs(floor(delta0)-floor(delta1))>=1.0) dirichlet *= -1.0;

         COMPLEX16 realRatio = RatioTerm0*phaseshift0*conj(dirichlet0);
         COMPLEX16 estRatio = RatioTerm*phaseshift*conj(dirichlet);

         fprintf(OUTPUT, "%g %g %g %g\n", cabs(realRatio), gsl_sf_angle_restrict_pos(carg(realRatio)), cabs(estRatio), gsl_sf_angle_restrict_pos(carg(estRatio)));
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
   uvar->unrestrictedCosi = 0;
   uvar->rectWindow = 0;

   XLALRegisterUvarMember(  Tsft,             REAL8, 0 , OPTIONAL , "SFT coherence time");
   XLALRegisterUvarMember(  SFToverlap,       REAL8, 0 , OPTIONAL , "SFT overlap in seconds, usually Tsft/2");
   XLALRegisterUvarMember(  t0,               REAL8, 0 , OPTIONAL , "GPS start time of the search");
   XLALRegisterUvarMember(  Tobs,             REAL8, 0 , OPTIONAL , "Duration of the search (in seconds)");
   XLALRegisterUvarMember(  cosi,             REAL8, 0 , OPTIONAL , "Cosine of NS inclinaiont angle");
   XLALRegisterUvarMember(  psi,              REAL8, 0 , OPTIONAL , "Polarization angle of GW");
   XLALRegisterUvarMember(  alpha,            REAL8, 0 , OPTIONAL , "Right ascension of source (in radians)");
   XLALRegisterUvarMember(  delta,            REAL8, 0 , OPTIONAL , "Declination of source (in radians)");
   XLALRegisterUvarMember(   skylocations,     INT4, 0 , OPTIONAL , "Number of sky locations");
   XLALRegisterUvarMember(  IFO,              STRINGVector, 0 , REQUIRED , "CSV list of detectors, eg. \"H1,H2,L1,G1, ...\" ");
   XLALRegisterUvarMember(outfilename,      STRING, 0 , OPTIONAL , "Output filename");
   XLALRegisterUvarMember(ephemEarth,       STRING, 0 , OPTIONAL , "Earth ephemeris file");
   XLALRegisterUvarMember(ephemSun,         STRING, 0 , OPTIONAL , "Sun ephemeris file");
   XLALRegisterUvarMember(  unrestrictedCosi, BOOLEAN, 0 , OPTIONAL , "Marginalize over cos(iota) from -1 to 1");
   XLALRegisterUvarMember(  rectWindow,       BOOLEAN, 0 , OPTIONAL , "Use rectangular window function instead of Hann windowing");

   BOOLEAN should_exit = 0;
   XLAL_CHECK( XLALUserVarReadAllInput( &should_exit, argc, argv ) == XLAL_SUCCESS, XLAL_EFUNC );
   if ( should_exit ) exit (1);

   return XLAL_SUCCESS;
}
