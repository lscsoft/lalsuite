/*********************************************************************************/
/*     Cosmological functions for cosmic string burst computation (header)       */
/*                                                                               */
/*                  Jolien Creighton, Irit Maor, Xavier Siemens                  */
/*                                                                               */
/*                         UWM/Caltech - September 2006                          */
/*********************************************************************************/

char cvsinfo[]  = "$Id$" "$Name$";

#ifndef CS_COSMO_H
#define CS_COSMO_H
#include <stdlib.h>

typedef struct cs_cosmo_functions {
	double *phit;
	double *phiA;
	double *phiV;
	double *z;
	double  zmin;
	double  dlnz;
	size_t  n;
} cs_cosmo_functions_t;

extern const char * cosmoname;
extern double zeq;
extern double H0;
cs_cosmo_functions_t cs_cosmo_functions( double *z, size_t n );
cs_cosmo_functions_t cs_cosmo_functions_alloc( double zmin, double dlnz, size_t n );
#define cs_cosmo_functions_free(cosmofns) \
	((void)(free((cosmofns).z),free((cosmofns).phit),free((cosmofns).phiA),free((cosmofns).phiV)))
#endif /* CS_COSMO_H */
