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

static
UINT4
modulo_sum(UINT4 a, UINT4 b, UINT4 mod)
{
  return (UINT4)fmod((double)a+(double)b, (double)mod);
}

/* returns azimuthal angle of Earth's location in ICRS */
static
REAL8 
earth_location(LALStatus *s, LIGOTimeGPS *p_gps, EphemerisData *p_eph)
{
  EarthState earth;

  LALBarycenterEarth(s, &earth, p_gps, p_eph);
  
  return atan2(earth.posNow[1], earth.posNow[0]);
}


/* fills in an irregular grid */
static
UINT4 
irr_grid_ra(LALStatus *s, REAL8Vector *ra_grid, REAL8 earth_phi)
{
  REAL8Vector *g = NULL;
  UINT4 num_grid = ra_grid->length;
  UINT4 num_grid_4 = num_grid / 4;
  UINT4 num_grid_2 = num_grid_4 * 2;
  UINT4 num_grid_3_4 = num_grid_4 * 3;
  UINT4 i, offset;
  REAL8 a, b;
  REAL8 twopi_num_grid = (REAL8)LAL_TWOPI/(REAL8)num_grid;

  LALDCreateVector(s, &g, num_grid);
  
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
  
    g->data[num_grid_4] = asin(-1. + a + b*(REAL8)num_grid_4);
    g->data[num_grid_2] = asin(-1. + a + b*(REAL8)num_grid_2);
    for (i = 0; i < num_grid_4; ++i)
    {
      g->data[i] = asin(-1. + a + b*(REAL8)i);
      g->data[num_grid_2 - i - 1] = -g->data[i];
      g->data[num_grid_2 + i] = g->data[i] + (REAL8)LAL_PI;
      g->data[num_grid - i - 1] = -g->data[i] + (REAL8)LAL_PI;
    }
  }
  
  if (lalDebugLevel > 2)
  {
    printf("irr_grid_ra (sin(ra)):\n");
    for (i = 0; i < num_grid; ++i)
      printf("\t% e\n", sin(g->data[i]));
  }
  
  /* shift the computed grid above, and fill in the RA grid */
  offset = (UINT4)rint(earth_phi / twopi_num_grid);
  for (i = 0; i < num_grid; ++i)
    ra_grid->data[modulo_sum(i, offset, num_grid)] = g->data[i];
  
  LALDDestroyVector(s, &g);
  
  return num_grid;
}


void 
init_gridding(gridding_t *g)
{
  g->ra = NULL;
  g->ra_irr = NULL;
  g->dec = NULL;
} /* END: init_gridding() */

void 
make_gridding(LALStatus *s, gridding_t *g, 
              UINT4 num_ra, gridding_geom_t ra_geom, 
              UINT4 num_dec, gridding_geom_t dec_geom,
              LIGOTimeGPS *gps)
{
  UINT4 i, j;
  REAL8 pi_num_ra = (REAL8)LAL_PI/(REAL8)num_ra;
  EphemerisData ephemerides;
  REAL8         earth_phi; /* azi. position of Earth in solar
                              system barycenter */
  
  if (dec_geom == DETRESP_VARGRID)
  {
    fprintf(stderr, "gridding_t does not support variable Dec gridding\n");
    exit(11);
  }
  
  g->ra_geom = ra_geom;
  g->dec_geom = dec_geom;
  

  LALDCreateVector(s, &(g->dec), num_dec);
  
  if (g->dec_geom == DETRESP_REGGRID)
  {
    for (i = 0; i < g->dec->length; ++i)
      g->dec->data[i] = (REAL8)LAL_PI * (-0.5 + (0.5 + i)/(REAL8)g->dec->length);
  }
  else if (g->dec_geom == DETRESP_IRRGRID)
  {
    for (i = 0; i < g->dec->length; ++i)
      g->dec->data[i] = asin(-1. + (REAL8)(1 + 2*i)/(REAL8)g->dec->length);
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
      /* gridding depends on Earth's position in orbit */ 
      init_ephemeris(s, &ephemerides);
      
      earth_phi = earth_location(s, gps, &ephemerides);
      
      if (num_ra != irr_grid_ra(s, g->ra, earth_phi))
      {
        fprintf(stderr, "Error in irr_grid_ra()\n");
        exit(13);
      }
      
      cleanup_ephemeris(s, &ephemerides);
    }
  }
  else
  {
    g->ra_irr = (REAL8Vector**)xcalloc(sizeof(REAL8Vector**), g->dec->length);
    
    for (j = 0; j < g->dec->length; ++j)
    {
      i = (UINT4)rint(sin(g->dec->data[j]));
      LALDCreateVector(s, &(g->ra_irr[j]), i);
    
      /* the gridding at each Dec is regular */
      for (i = 0; i < (g->ra_irr[j])->length; ++i)
        (g->ra_irr[j])->data[i] = pi_num_ra * (1. + 2.*(REAL8)i);
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
print_gridding(gridding_t *g, char *fn)
{
  UINT4 i, j;
  FILE *outfile = NULL;
  
  if (fn != NULL)
    outfile = xfopen(fn, "wo");
  else
    outfile = stdout;
    
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
  for (i = 0; i < g->dec->length; ++i)
    fprintf(outfile, "\t% 20.14e\n", g->dec->data[i]/(double)LAL_PI*180.);
  fprintf(outfile, "\n");
    
  fprintf(outfile, "Right Ascencion (h):\n");
  if (g->ra_geom == DETRESP_REGGRID || g->ra_geom == DETRESP_IRRGRID)
  {
    for (i = 0; i < g->ra->length; ++i)
      fprintf(outfile, "\t% 20.14e\n", g->ra->data[i]/(double)LAL_PI * 12.);
  }
  else
  {
    fprintf(outfile, "RSN\n");
  }

  fclose(outfile);  
} /* END: print_gridding() */

