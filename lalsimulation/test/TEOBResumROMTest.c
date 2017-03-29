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
#include <lal/LALgetopt.h>

static void print_usage(char *program)
{
  fprintf(stderr,
          "%s [options]\n"\
          "The following options are recognized.  Options not surrounded in []\n"\
          "are required. Defaults are shown in brackets\n", program );
  fprintf(stderr,
          " [--help ]                 display this message\n"\
          " [--mass1 ] mass1          mass1 (1.4)\n"\
          " [--mass2 ] mass2          mass2 (1.4)\n"\
          " [--lambda1 ] lambda1      lambda1 (5000)\n"\
          " [--lambda2 ] lambda2      lambda2 (5000)\n"\
          " [--f-lower ] fmin         start frequency (30)\n"\
          " [--sample-rate ] srate    deltaT = 1/srate (4096)\n"\
          " [--no-output ]            do not write out data to ./rom_out.txt. For performance check\n"\
          " [--distance] dist         distance in Mpc (100)\n\n");
}

#define UNUSED(expr) do { (void)(expr); } while (0)

int main (int argc, char *argv[])
{

  REAL8 mass1 = 1.4 ;
  REAL8 mass2 = 1.4 ;
  REAL8 lambda1 = 5000.0 ;
  REAL8 lambda2 = 5000.0 ;
  REAL8 fMin = 30.0 ;
  REAL8 srate = 4096.0 ;
  REAL8 distance = 100.0 ;
  int no_out = 0 ;

  /* LALgetopt arguments */
  struct LALoption long_options[] =
  {
    {"help",                    no_argument,       0,                'h'},
    {"mass1",                   required_argument, 0,                'm'},
    {"mass2",                   required_argument, 0,                'M'},
    {"lambda1",                 required_argument, 0,                'l'},
    {"lambda2",                 required_argument, 0,                'L'},
    {"f-lower",                 required_argument, 0,                'f'},
    {"sample-rate",             required_argument, 0,                's'},
    {"distance",                required_argument, 0,                'd'},
    {"no-output",               no_argument,       &no_out,          1},
    {0, 0, 0, 0}
  };
  int c;

  while ( 1 )
  {
    /* LALgetopt_long stores long option here */
    int option_index = 0;

    c = LALgetopt_long_only( argc, argv,"h:m:M:l:L:f:s:d:*", long_options, &option_index );

    /* detect the end of the options */
    if ( c == - 1 )
    {
      break;
    }

    switch ( c )
    {
      case 0:
        /* if this option set a flag, do nothing else now */
        if ( long_options[option_index].flag != 0 )
        {
          break;
        }
        else
        {
          fprintf( stderr, "error parsing option %s with argument %s\n",
                   long_options[option_index].name, LALoptarg );
          exit( 1 );
        }
        break;
      case 'm':
        mass1 = atof( LALoptarg );
        break;
      case 'M':
        mass2 = atof( LALoptarg );
        break;
      case 'l':
        lambda1 = atof( LALoptarg );
        break;
      case 'L':
        lambda2 = atof( LALoptarg );
        break;
      case 'f':
        fMin = atof( LALoptarg );
        break;
      case 'd':
        distance = atof( LALoptarg );
        break;
      case 's':
        srate = atof( LALoptarg );
        break;
      case 'h':
        print_usage(argv[0]);
        exit( 0 );
        break;
      case '?':
        print_usage(argv[0]);
        exit( 1 );
        break;
      default:
        fprintf( stderr, "unknown error while parsing options\n" );
        print_usage(argv[0]);
        exit( 1 );
    }
  }


  REAL8TimeSeries *hplus  = NULL;
  REAL8TimeSeries *hcross = NULL;

  REAL8 phiC = 0.;
  REAL8 deltaT = 1./srate;
  REAL8 r      = distance*1.0e6 * LAL_PC_SI;
  REAL8 i      = 0.;

  UNUSED(argc);
  UNUSED(argv);

  if ( XLALSimInspiralTEOBResumROM(&hplus,&hcross,phiC,deltaT,fMin,0.0,r,i,mass1*LAL_MSUN_SI,mass2*LAL_MSUN_SI,lambda1,lambda2) == XLAL_FAILURE )
  {
    fprintf( stderr, "The waveform generation function has failed!!\n" );
    exit(1);
  }

  if (no_out == 0) {
    double t0 = XLALGPSGetREAL8(&hplus->epoch);

    int n=0;
    int N=(int) hplus->data->length;
    FILE *out = fopen("./rom_out.txt","w");
    for (n=0;n<N;n++){
      fprintf(out,"%.9e %.9e %.9e\n",t0 + n * hplus->deltaT,hplus->data->data[n],hcross->data->data[n]);
    }
    fclose(out);
  }

  return 0;
}
