/* Code to check that the function to read TEMPO-style par files return the expected values
 *
 * Matthew Pitkin 21/01/2014
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/ReadPulsarParFile.h>
#include <lal/LALString.h>

#define PARFILE "test.par" /* define par file name */

typedef struct tagParamData{
  const CHAR *name;  /* parameter name as given by the conventions (see param variable in TEMPOs mxprt.f file */
  const CHAR *val;   /* parameter value for output to par file */
  const CHAR *valcheck; /* if parameter will be converted to SI units, this is the conversion (other wise
                           * set to the same value as val) given to 5 d.p. */
  const CHAR *sigma; /* standard deviation on the parameter as read from the .par file */
  const CHAR *sigmacheck; /* check of read in sigma converted to SI units */
  const CHAR *fitFlag; /* add a TEMPO-style fitting flag to some parameters */
}ParamData;

#define NUMPARS 54

/* setup a collection of allowed .par file parameters */
ParamData p[NUMPARS] =
{
  /* frequency parameters */
  { "F0",       "123.00000",    "123.00000",        "1.00000e-08",  "1.00000e-08",  "1" },
  { "F1",       "-1.23000e-10", "-1.23000e-10",     "1.00000e-12",  "1.00000e-12",  "1" },
  { "F2",       "1.23000e-15",  "1.23000e-15",      "1.00000e-17",  "1.00000e-17",  "1" },
  { "F3",       "-1.23000e-20", "-1.23000e-20",     "1.00000e-22",  "1.00000e-22",  "1" },
  { "F4",       "1.23000e-25",  "1.23000e-25",      "1.00000e-27",  "1.00000e-27",  "1" },
  { "F5",       "-1.23000e-30", "-1.23000e-30",     "1.00000e-32",  "1.00000e-32",  " " },
  { "F6",       "1.23000e-35",  "1.23000e-35",      "1.00000e-37",  "1.00000e-37",  "1" },
  { "F7",       "-1.23000e-40", "-1.23000e-40",     "1.00000e-42",  "1.00000e-42",  "1" },
  { "F8",       "1.23000e-45",  "1.23000e-45",      "1.00000e-47",  "1.00000e-47",  "1" },
  { "F9",       "-1.23000e-50", "-1.23000e-50",     "1.00000e-52",  "1.00000e-52",  " " },
  { "F10",      "1.23000e-67",  "1.23000e-67",      "1.00000e-68",  "1.00000e-68",  " " },
  { "PEPOCH",   "54321.0",      "870652748.81600",  NULL,           NULL,           " " },
  /* name parameters */
  { "PSRJ",     "J0000+0000",   "J0000+0000",       NULL,           NULL,           " " },
  { "PSRB",     "B0000+0000",   "B0000+0000",       NULL,           NULL,           " " },
  { "PSR",      "0000+0000",    "0000+0000",        NULL,           NULL,           " " },
  { "NAME",     "J0000+0000",   "J0000+0000",       NULL,           NULL,           " " },
  /* position parameters */
  { "RAJ",      "12:34:56.7",   "3.29407",          "1.23",         "8.94481e-05",  "1" },
  { "ra",       "12:34:56.7",   "3.29407",          "1.23",         "8.94481e-05",  " " }, /* use lower case input name*/
  { "DECJ",     "-12:34:56.7",  "-0.21960",         "1.23",         "5.96321e-06",  " " },
  { "DEC",      "-12:34:56.7",  "-0.21960",         "1.23",         "5.96321e-06",  " " },
  { "POSEPOCH", "54321.0",      "870652748.81600",  NULL,           NULL         ,  " " },
  { "PMRA",     "12.345",       "1.89654e-15",      "1.23",         "1.88963e-16",  " " },
  { "EPHEM",    "DE405",        "DE405",            NULL,           NULL         ,  " " },
  { "PMDEC",    "-12.345",      "-1.89654e-15",     "1.23",         "1.88963e-16",  "1" },
  { "DIST",     "123.00000",    "3.79538e+21",      "12.30000",     "3.79538e+20",  " " },
  /* glitch parameters */
  { "GLEP_1",   "54321.0",      "870652748.81600",  "0.00123",      "106.27200",    " " },
  { "GLEP_2",   "54322.0",      "870739148.81600",  "0.00123",      "106.27200",    " " },
  { "GLPH_1",   "0.43453",      "0.43453",          "0.00453",      "0.00453",      "1" },
  { "GLPH_2",   "0.93453",      "0.93453",          "0.01453",      "0.01453",      "1" },
  { "GLF0_1",   "2.50000e-8",   "2.50000e-08",      "4.50000e-09",  "4.50000e-09",  " " },
  { "GLF0_2",   "3.50000e-8",   "3.50000e-08",      "2.50000e-09",  "2.50000e-09",  " " },
  { "GLF0D_1",  "5.50000e-6",   "5.50000e-06",      "7.50000e-09",  "7.50000e-09",  " " },
  { "GLF0D_2",  "6.50000e-6",   "6.50000e-06",      "9.50000e-09",  "9.50000e-09",  " " },
  { "GLTD_1",   "234.500000",   "20260800.00000",   "2.00000",      "172800.00000", " " },
  { "GLTD_2",   "23.450000",    "2026080.00000",    "1.00000",      "86400.00000",  " " },
  /* binary parameters */
  { "BINARY",   "BT",           "BT",               NULL,           NULL         ,  " " },
  { "OM",       "123.45",       "2.15461",          "1.23",         "0.02147"    ,  " " },
  { "A1",       "12.34500",     "12.34500",         "1.23000e-04",  "1.23000e-04",  "1" },
  { "ECC",      "1.23400e-05",  "1.23400e-05",      "1.23400e-08",  "1.23400e-08",  " " },
  { "PB",       "1.23",         "106272.00000",     "0.00123",      "106.27200"  ,  " " },
  { "T0",       "54321.0",      "870652748.81600",  "0.00123",      "106.27200"  ,  " " },
  { "TASC",     "54321.0",      "870652748.81600",  "0.00123",      "106.27200"  ,  " " },
  { "EPS1",     "1.23400e-05",  "1.23400e-05",      "1.23400e-08",  "1.23400e-08",  "1" },
  { "EPS2",     "1.23400e-05",  "1.23400e-05",      "1.23400e-08",  "1.23400e-08",  "1" },
  { "FB0",      "1.23400e-05",  "1.23400e-05",      "1.23400e-14",  "1.23400e-14",  "1" },
  { "FB2",      "1.23400e-09",  "1.23400e-09",      "1.23400e-18",  "1.23400e-18",  " " },
  { "FB1",      "1.23400e-09",  "1.23400e-09",      "1.23400e-18",  "1.23400e-18",  "1" },
  { "EDOT",     "1.23400e-05",  "1.23400e-17",      "1.23400e-18",  "1.23400e-18",  " " },
  /* GW parameters */
  { "H0",       "1.23000e-22",  "1.23000e-22",      "1.23000E-23",  "1.23000e-23",  " " }, /* input exponent as E */
  { "COSIOTA",  "-0.12300",     "-0.12300",         "0.00123",      "0.00123"    ,  " " },
  { "PHI0",     "1.23000",      "1.23000",          "0.12300",      "0.12300"    ,  " " },
  { "PSI",      "-0.12300",     "-0.12300",         "0.01230",      "0.01230"    ,  " " },
  { "APLUS",    "1.23000e-22",  "1.23000e-22",      "1.23000D-23",  "1.23000e-23",  " " }, /* input exponent as D */
  { "ACROSS",   "1.23000e-22",  "1.23000e-22",      "1.23000D-23",  "1.23000e-23",  " " }, /* input exponent as D */
};


int main( void ){
  INT4 i = 0;
  PulsarParameters *pars = NULL;

  /* output par file */
  FILE *fp = NULL;
  if ( (fp = fopen(PARFILE, "w")) == NULL ){
    XLAL_ERROR( XLAL_EIO, "Error... problem writing parfile!\n" );
  }

  for( i=0; i<NUMPARS; i++ ){
    if ( p[i].sigma != NULL ){ /* write out value and error */
      fprintf(fp, "%s\t%s\t%s\t%s\n", p[i].name, p[i].val, p[i].fitFlag, p[i].sigma);
    }
    else{ /* just write out value */
      fprintf(fp, "%s\t%s\n", p[i].name, p[i].val);
    }
  }

  fclose(fp);

  /* read in par file */
  pars = XLALReadTEMPOParFileNew( PARFILE );

  /* check read-in parameters against originals */
  for ( i=0; i<NUMPARS; i++ ){
    CHAR outval[256];
    CHAR outsigma[256];
    CHAR outfitflag[256];

    /* see if the parameter can be split into an alphabet and numerical value (as is the case
     * for e.g. FB, which is held as a vector */
    CHAR *namecopy = NULL, *token = NULL;
    namecopy = XLALStringDuplicate( p[i].name );
    /* get part of the name before the numerical/underscore delimiter */
    if ( strchr( p[i].name, '_' ) && strncmp( p[i].name, "WAVE", 4 ) ) { /* delimiter by underscores and not a FITWAVES parameter */
      token = strtok(namecopy, "_");
    }
    else{ token = strtok(namecopy, "0123456789"); }

    if ( !PulsarCheckParam( pars, p[i].name ) && !PulsarCheckParam( pars, token ) ){
      XLAL_ERROR( XLAL_EFAILED, "Error... parameter %s does not exist in the read-in data!\n", p[i].name );
    }

    /* value is a vector */
    if ( PulsarCheckParam( pars, token ) && PulsarGetParamType( pars, token ) == PULSARTYPE_REAL8Vector_t ){
     /* get value and convert to string */
      if ( !strchr( p[i].valcheck, 'e' ) ){ /* number doesn't contain an exponent */
        sprintf(outval, "%.5lf", PulsarGetREAL8VectorParamIndividual(pars, p[i].name));
      }
      else{ /* number does contain an exponent */
        sprintf(outval, "%.5le", PulsarGetREAL8VectorParamIndividual(pars, p[i].name));
      }

      if ( p[i].sigmacheck != NULL ){
        /* get error and convert to string */
        if ( !strchr( p[i].sigmacheck, 'e' ) ){ /* number doesn't contain an exponent */
          sprintf(outsigma, "%.5lf", PulsarGetREAL8VectorParamErrIndividual(pars, p[i].name));
        }
        else{
          sprintf(outsigma, "%.5le", PulsarGetREAL8VectorParamErrIndividual(pars, p[i].name));
        }

        UINT4 *fitFlag = PulsarGetParamFitFlag(pars, token);

        if ( fitFlag[atoi(p[i].name+strlen(token))] == 1 ){ sprintf(outfitflag, "1"); }
        else if ( fitFlag[atoi(p[i].name+strlen(token))] == 0 ){ sprintf(outfitflag, " "); }
        else {
          XLAL_ERROR( XLAL_EFAILED, "Error... fit flag incorrect for %s.\n", p[i].name );
        }
      }
    }
    /* value is a double */
    else if( PulsarGetParamType( pars, p[i].name ) == PULSARTYPE_REAL8_t ){
      /* get value and convert to string */
      if ( !strchr( p[i].valcheck, 'e' ) ){ /* number doesn't contain an exponent */
        sprintf(outval, "%.5lf", PulsarGetREAL8Param(pars, p[i].name));
      }
      else{ /* number does contain an exponent */
        sprintf(outval, "%.5le", PulsarGetREAL8Param(pars, p[i].name));
      }

      if ( p[i].sigmacheck != NULL ){
        /* get error and convert to string */
        if ( !strchr( p[i].sigmacheck, 'e' ) ){ /* number doesn't contain an exponent */
          sprintf(outsigma, "%.5lf", PulsarGetREAL8ParamErr(pars, p[i].name));
        }
        else{
          sprintf(outsigma, "%.5le", PulsarGetREAL8ParamErr(pars, p[i].name));
        }

        UINT4 *fitFlag = PulsarGetParamFitFlag(pars, p[i].name);

        if ( fitFlag[0] == 1 ){ sprintf(outfitflag, "1"); }
        else if ( fitFlag[0] == 0 ){ sprintf(outfitflag, " "); }
        else {
          XLAL_ERROR( XLAL_EFAILED, "Error... fit flag incorrect for %s.\n", p[i].name );
        }
      }
    }
    /* value is a string */
    else if ( PulsarGetParamType( pars, p[i].name ) == PULSARTYPE_string_t ) {
      CHAR *out = XLALStringDuplicate(PulsarGetStringParam(pars, p[i].name));
      sprintf(outval, "%s", out);
      XLALFree(out);
    }

    /* compare returned value with input value */
    if ( strcmp(outval, p[i].valcheck) != 0 ){
      XLAL_ERROR( XLAL_EFAILED, "Error... parameter %s does not match input (%s cf. %s)!\n", p[i].name, outval,
                  p[i].valcheck );
    }

    if( p[i].sigma != NULL ){
      /* compare returned value with input value */
      if ( strcmp(outsigma, p[i].sigmacheck) != 0 ){
        XLAL_ERROR( XLAL_EFAILED, "Error... parameter sigma %s does not match input (%s cf. %s)!\n", p[i].name,
                    outsigma, p[i].sigmacheck );
      }

      if ( strcmp(outfitflag, p[i].fitFlag) != 0 ){
        XLAL_ERROR( XLAL_EFAILED, "Error... parameter %s fit flag does not match input!\n", p[i].fitFlag );
      }
    }

    XLALFree( namecopy );
  }

  /* remove par file */
  if ( remove(PARFILE) != 0 ){
    XLAL_ERROR( XLAL_EIO, "Error... problem removing parfile!\n" );
  }

  fprintf(stderr, "Successfully read in .par file values!\n");

  PulsarFreeParams( pars );

  LALCheckMemoryLeaks();

  return XLAL_SUCCESS;
}
