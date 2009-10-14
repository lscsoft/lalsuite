/*
*  Copyright (C) 2007 David Chin
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

/*
 * Author: David Chin <dwchin@umich.edu> +1-734-709-9119
 * $Id$
 */

/* 
 * Makes a grid given appropriate params
 */
 
#include "make_gridding.h"
#include "skygrid.h"
#include "util.h"

extern int lalDebugLevel;

/* 
 * private functions 
 */
static double rad_to_deg(double rad);
static double rad_to_hr(double rad);
static UINT4 modulo_sum(UINT4 a, UINT4 b, UINT4 mod);
static UINT4 irr_grid_ra(LALStatus *s, REAL8Vector *ra_grid, 
                         REAL8 earth_phi);
static REAL8 earth_location(LALStatus *s, LIGOTimeGPS *p_gps, 
                            EphemerisData *p_eph);
/* Mollweide gridding functions */ 
static void fprintf_laldvector_hor(FILE *f, REAL8Vector *v, char *delim,
                                   REAL8 scale);
static void fprintf_laldvector_ver(FILE *f, REAL8Vector *v, REAL8 scale);



void 
init_gridding(gridding_t *g)
{
  g->gps.gpsSeconds = 0;
  g->gps.gpsNanoSeconds = 0;
  g->ra = NULL;
  g->ra_irr = NULL;
  g->dec = NULL;
} /* END: init_gridding() */

void 
make_gridding(LALStatus *s, gridding_t *g, 
              UINT4 num_ra, gridding_geom_t ra_geom, 
              UINT4 num_dec, gridding_geom_t dec_geom,
              EphemerisData *e, LIGOTimeGPS   *gps )
{
  UINT4 i, j;
  REAL8 pi_num_ra = (REAL8)LAL_PI/(REAL8)num_ra;
  REAL8 earth_phi; /* azi. position of Earth in solar system barycenter */
  
  if (dec_geom == DETRESP_VARGRID)
  {
    fprintf(stderr, "gridding_t does not support variable Dec gridding\n");
    exit(11);
  }

  g->gps.gpsSeconds = gps->gpsSeconds;
  g->gps.gpsNanoSeconds = gps->gpsNanoSeconds;
  g->ra_geom = ra_geom;
  g->dec_geom = dec_geom;


  LALDCreateVector(s, &(g->dec), num_dec);
  
  if (g->dec_geom == DETRESP_REGGRID)
  {
    for (i = 0; i < g->dec->length; ++i)
      g->dec->data[i] = (REAL8)LAL_PI * 
        (-0.5 + (0.5 + i)/(REAL8)(g->dec->length));
  }
  else if (g->dec_geom == DETRESP_IRRGRID)
  {
    for (i = 0; i < g->dec->length; ++i)
      g->dec->data[i] = asin(-1. + 
                             (REAL8)(1 + 2*i)/(REAL8)(g->dec->length));
  }
 
  
  if (g->ra_geom != DETRESP_VARGRID)
  {
    LALDCreateVector(s, &(g->ra), num_ra);
    
    if (g->ra_geom == DETRESP_REGGRID)
    {
      for (i = 0; i < g->ra->length; ++i)
        g->ra->data[i] = pi_num_ra * (1. + 2.*(REAL8)i);
    }
    else if (g->ra_geom == DETRESP_IRRGRID)
    {
      earth_phi = earth_location(s, gps, e);
      
      if (num_ra != irr_grid_ra(s, g->ra, earth_phi))
      {
        fprintf(stderr, "Error in irr_grid_ra()\n");
        exit(13);
      }
    }
  }
  else
  {
    g->ra_irr = (REAL8Vector**)xcalloc(sizeof(REAL8Vector**), 
                                       g->dec->length);
    
    for (j = 0; j < g->dec->length; ++j)
    {
      /* number of RA grid points goes as cos(Dec) */
      i = (UINT4)rint(fabs(cos(g->dec->data[j]) * num_ra)) + 2;
      LALDCreateVector(s, &(g->ra_irr[j]), i);
      
      if (lalDebugLevel & 7)
        printf("dec = % e, cos(dec) = % e, num_ra = %d\n", 
               g->dec->data[j], cos(g->dec->data[j]), i);
    
      /* the gridding at each Dec is regular */
      for (i = 0; i < (g->ra_irr[j])->length; ++i)
        (g->ra_irr[j])->data[i] = (REAL8)LAL_PI * (1. + 2.*(REAL8)i)
                                  / (REAL8)((g->ra_irr[j])->length);
    }
  }
} /* END: init_gridding() */


/*
 * deallocate gridding
 */
void 
cleanup_gridding(LALStatus *s, gridding_t *g)
{
  UINT4 j;
  
  if (g->ra_geom == DETRESP_VARGRID)
  {
    for (j = 0; j < g->dec->length; ++j)
    {
      LALDDestroyVector(s, &(g->ra_irr[j]));
      g->ra_irr[j] = NULL;
    }
    free(g->ra_irr);
  }
  else
  {
    LALDDestroyVector(s, &(g->ra));
  }
  
  LALDDestroyVector(s, &(g->dec));

  init_gridding(g);
} /* END: cleanup_gridding() */


/*
 * Fill gridding with zeros
 */
void 
zero_gridding(LALStatus *s, gridding_t *g)
{
  UINT4 i;
  
  if (g->ra == NULL || g->dec == NULL)
  {
    exit(11);
  }
  
  for (i = 0; i < g->ra->length; ++i)
    g->ra->data[i] = 0.;
    
  for (i = 0; i < g->dec->length; ++i)
    g->dec->data[i] = 0.;
} /* END: zero_gridding() */


void 
print_gridding(gridding_t *g, char *fn, gridding_printmode_t mode)
{
  UINT4 i, j;
  FILE *outfile = NULL;
  
  if (fn != NULL)
    outfile = xfopen(fn, "wo");
  else
    outfile = stdout;
  
  switch (mode)
  {
    case DETRESP_HUMANREAD:
      fprintf(outfile, "GPS = %d:%d\n", g->gps.gpsSeconds,
              g->gps.gpsNanoSeconds);
      switch (g->ra_geom)
      {
        case DETRESP_REGGRID:
          fprintf(outfile, "ra_geom =\tREGULAR\n");
          break;
        case DETRESP_IRRGRID:
          fprintf(outfile, "ra_geom =\tIRREGULAR\n");
          break;
        case DETRESP_VARGRID:
          fprintf(outfile, "ra_geom =\tVARIABLE\n");
          break;
        default:
          fprintf(outfile, "ra_geom =\t???\n");
      }
      
      switch (g->dec_geom)
      {
        case DETRESP_REGGRID:
          fprintf(outfile, "dec_geom =\tREGULAR\n");
          break;
        case DETRESP_IRRGRID:
          fprintf(outfile, "dec_geom =\tIRREGULAR\n");
          break;
        case DETRESP_VARGRID:
          fprintf(outfile, "dec_geom =\tVARIABLE\n");
          break;
        default:
          fprintf(outfile, "dec_geom =\t???\n");
      }
      
      fprintf(outfile, "Declination (deg):\n");
      fprintf_laldvector_ver(outfile, g->dec, 180./(REAL8)LAL_PI);
    
      fprintf(outfile, "\n");
        
      fprintf(outfile, "Right Ascenscion (h):\n");
      if (g->ra_geom == DETRESP_REGGRID || g->ra_geom == DETRESP_IRRGRID)
      {
        fprintf_laldvector_ver(outfile, g->ra, 12./(REAL8)LAL_PI);
      }
      else
      {
        for (i = 0; i < g->dec->length; ++i)
          fprintf_laldvector_hor(outfile, g->ra_irr[i], ", ", 
                                 12./(REAL8)LAL_PI);
      }
      break; /* DETRESP_HUMANREAD */
      
    case DETRESP_XYPAIRS_ASCII:
      if (g->ra_geom != DETRESP_VARGRID)
      {
        for (i = 0; i < g->dec->length; ++i)
          for (j = 0; j < g->ra->length; ++j)
            fprintf(outfile, "% 20.14e\t% 20.14e\n", 
                    rad_to_hr(g->ra->data[j]), 
                    rad_to_deg(g->dec->data[i]));
      }
      else
      {
        for (i = 0; i < g->dec->length; ++i)
          for (j = 0; j < (g->ra_irr[i])->length; ++j)
            fprintf(outfile, "% 20.14e\t% 20.14e\n",
                    rad_to_hr((g->ra_irr[i])->data[j]), 
                    rad_to_deg(g->dec->data[i]));
      }
      break; /* DETRESP_XYPAIRS_ASCII */
      
    case DETRESP_XYPAIRS_BIN:
      /* FIXME */
      break; /* DETRESP_XYPAIRS_BIN */
  }
  
  if (outfile != stdout)
    fclose(outfile);  
} /* END: print_gridding() */


void
print_ra_grid(gridding_t *g, char *fn)
{
  UINT4 i, j;
  FILE *outfile = NULL;

  if (fn != NULL)
    outfile = xfopen(fn, "wo");
  else
    outfile = stdout;

  for (i = 0; i < g->ra->length; ++i)
  {
    for (j = 0; j < g->dec->length; ++j)
    {
      fprintf(outfile, "% 20.14e\t", g->ra->data[i]);
    }
    fprintf(outfile, "\n");
  }
  fclose(outfile);
} /* END: print_ra_grid() */


void
print_dec_grid(gridding_t *g, char *fn)
{
  UINT4 i, j;
  FILE *outfile = NULL;

  if (fn != NULL)
    outfile = xfopen(fn, "wo");
  else
    outfile = stdout;

  for (i = 0; i < g->ra->length; ++i)
  {
    for (j = 0; j < g->dec->length; ++j)
    {
      fprintf(outfile, "% 20.14e\t", g->dec->data[j]);
    }
    fprintf(outfile, "\n");
  }
  fclose(outfile);
} /* END: print_ra_grid() */



static 
double 
rad_to_deg(double rad)
{
  return (rad * 180. / (double)LAL_PI);
} 


static 
double
rad_to_hr(double rad)
{
  return (fmod(rad * 12. / (double)LAL_PI, 24.));
}

static
UINT4
modulo_sum(UINT4 a, UINT4 b, UINT4 mod)
{
  return (UINT4)fmod((double)a+(double)b, (double)mod);
} /* END: modulo_sum() */


/* fills in an irregular grid 
 * FIXME: we can deal only with an even number of grid points */
static
UINT4 
irr_grid_ra(LALStatus *s, REAL8Vector *ra_grid, REAL8 earth_phi)
{
  UINT4 num_grid = ra_grid->length;
  UINT4 num_grid_4 = num_grid / 4;
  UINT4 num_grid_2 = num_grid_4 * 2;
  UINT4 num_grid_3_4 = num_grid_4 * 3;
  UINT4 i;
  REAL8 a, b;
  
  if (num_grid != 0)
  {
    if (lalDebugLevel > 1)
    {
      fprintf(stderr, "num_grid_4   = %d\n", num_grid_4);
      fprintf(stderr, "num_grid_2   = %d\n", num_grid_2);
      fprintf(stderr, "num_grid_3_4 = %d\n", num_grid_3_4);
      fprintf(stderr, "num_grid     = %d\n", num_grid);
    }
    
    a = 2. / (REAL8)num_grid;
    b = 2. * a;

    for (i = 0; i < num_grid_2; ++i)
      ra_grid->data[i] = asin(-1. + a + b*(REAL8)i) + earth_phi;
      
    for (i = 0; i < num_grid_2; ++i)
      ra_grid->data[num_grid_2 + i] = ra_grid->data[i] + (REAL8)LAL_PI;
  }
  
  if (lalDebugLevel & 7)
  {
    printf("earth_phi = %e\n", earth_phi);
    printf("irr_grid_ra (rad):\n");
    for (i = 0; i < num_grid; ++i)
      printf("\t% e\n", ra_grid->data[i]);
  }
    
  return num_grid;
} /* END: irr_grid_ra() */


/* returns azimuthal angle of Earth's location in ICRS */
static
REAL8 
earth_location(LALStatus *s, LIGOTimeGPS *p_gps, EphemerisData *p_eph)
{
  EarthState earth;

  LALBarycenterEarth(s, &earth, p_gps, p_eph);
  
  if (lalDebugLevel & 7)
  {
    printf("earth.posNow = %e, %e, %e\n", earth.posNow[0],
           earth.posNow[1], earth.posNow[2]);
  }
  
  return atan2(earth.posNow[1], earth.posNow[0]);
} /* END: earth_location() */

static 
void 
fprintf_laldvector_hor(FILE *f, REAL8Vector *v, char *delim, REAL8 scale)
{
  UINT4 i;
  
  for (i = 0; i < v->length - 1; ++i)
    fprintf(f, "% 20.14e%s", v->data[i] * scale, delim);
  fprintf(f, "% 20.14e\n", v->data[i] * scale);
    
  return;
} /* END: fprintf_laldvector_hor() */


static 
void 
fprintf_laldvector_ver(FILE *f, REAL8Vector *v, REAL8 scale)
{
  UINT4 i;

  for (i = 0; i < v->length; ++i)
    fprintf(f, "\t% 20.14e\n", v->data[i] * scale);
    
  return;
} /* END: fprintf_laldvector_ver() */

