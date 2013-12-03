#include <stdio.h>
#include <stdlib.h>

#include <lal/Units.h>
#include <lal/LALAdaptiveRungeKutta4.h>
#include <lal/FindRoot.h>
#include <lal/SeqFactories.h>
#include <lal/LALSimInspiral.h>
#include <lal/LALSimIMR.h>
#include <lal/Date.h>
#include <lal/TimeSeries.h>
#include <gsl/gsl_sf_gamma.h>

//#include <lal/LALSimulation.h>
//#include <lal/LALSimIMREOBNRv2.h>

#define UNUSED(expr) do { (void)(expr); } while (0)

int main (int argc, char *argv[])
{

  REAL8TimeSeries *hplus  = NULL;
  REAL8TimeSeries *hcross = NULL;

  REAL8 phiC = 0.;
  REAL8 deltaT = 1./16384.;
  REAL8 m1     = 75. * LAL_MSUN_SI;
  REAL8 m2     = 25. * LAL_MSUN_SI;
  REAL8 fMin   = 10;
  REAL8 r      = 1.0e6 * LAL_PC_SI;
  REAL8 i      = 0./*LAL_PI / 6.*/;
  //REAL8 spin1[3], spin2[3];

  /*memset( spin1, 0, sizeof(spin1) );
  memset( spin2, 0, sizeof(spin2) );
  spin1[2] = 0.;
  spin2[2] = 0.;*/


  UNUSED(argc);
  UNUSED(argv);
  
  if ( XLALSimIMREOBNRv2DominantMode( &hplus, &hcross, phiC, deltaT, m1, m2, fMin, r, i) == XLAL_FAILURE )
  {
    fprintf( stderr, "The waveform generation function has failed!!\n" );
    exit(1);
  }

  FILE *out = fopen( "eobHPlusCross.dat", "w" );

  for (unsigned int j = 0; j < hplus->data->length; j++ )
  {
    fprintf( out, "%e %e %e\n", j * deltaT / ((m1+m2)/LAL_MSUN_SI*LAL_MTSUN_SI ), hplus->data->data[j], hcross->data->data[j] );
  }
  fclose( out );

  return 0;
}
