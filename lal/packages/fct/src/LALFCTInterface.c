#include <lal/LALStdlib.h>
#include <lal/LALFCTInterface.h>
#include <lal/LALMalloc.h>
#include <config.h>

#define LALFCT

#include "fct.c"
#include "fct_fft.c"

NRCSID(LALFCTINTERFACEC, "$Id:");

struct
tagLALFCTPlan
{
    fct_plan* fctPlan;
    fct_status* status;
    UINT4 num_data_cubes;
    BOOLEAN fct_initialised;
}
/* LALFCTPlan */;

void LALCreateFCTPlan(LALStatus* const status,
                      LALFCTPlan** plan_ptr,
                      const LALCreateFCTPlanInput* const in)
{
    INT4 i = 0;
    
    INITSTATUS(status, "LALCreateFCTPlan", LALFCTINTERFACEC);
    
    /* No harm in checking this each time */
    ASSERT((sizeof(fct_real) == sizeof(REAL4)), status,
           LALFCTINTERFACEH_EDATATYPE, LALFCTINTERFACEH_MSGEDATATYPE);
    ASSERT(*plan_ptr == 0, status,
           LALFCTINTERFACEH_ENNUL, LALFCTINTERFACEH_MSGENNUL);
    ASSERT(in != 0, status,
           LALFCTINTERFACEH_ENULL, LALFCTINTERFACEH_MSGENULL);

    ASSERT(in->data_length > 0, status,
           LALFCTINTERFACEH_ENPOS, LALFCTINTERFACEH_MSGENPOS);

    /* No point doing a 1-d FCT */
    ASSERT(in->number_of_dimensions > 1, status,
           LALFCTINTERFACEH_ENPOS, LALFCTINTERFACEH_MSGENPOS);
    ASSERT(in->dimension_0_stride > 0, status,
           LALFCTINTERFACEH_ENPOS, LALFCTINTERFACEH_MSGENPOS);

    for (i = 1; i < in->number_of_dimensions; ++i)
    {
        ASSERT(in->phase_fn[i-1] != 0, status,
               LALFCTINTERFACEH_ENULL, LALFCTINTERFACEH_MSGENULL);
    }

    /* Set up the fct malloc hooks */
    fct_malloc_hook = LALMalloc;
    fct_calloc_hook = LALCalloc;
    fct_free_hook = LALFree;
    
    /* Also set the FFTW hooks *unless* they've already been set */
    if (fftw_malloc_hook == 0)
    {
	fftw_malloc_hook = LALMalloc;
    }
    
    if (fftw_free_hook == 0)
    {
	fftw_free_hook = LALFree;
    }

    *plan_ptr = LALCalloc(1, sizeof(**plan_ptr));
    
    (*plan_ptr)->status = LALCalloc(1, sizeof(*(*plan_ptr)->status));

    (*plan_ptr)->fctPlan = fct_init_plan(in->data_length,
                                         in->number_of_dimensions,
                                         in->dimension_0_stride,
                                         (*plan_ptr)->status);

    ASSERT((*plan_ptr)->status->fct_errno == 0, status,
    LALFCTINTERFACEH_EINTERNAL, fct_strerror((*plan_ptr)->status->fct_errno));

    for (i = 1; i < in->number_of_dimensions; ++i)
    {
        fct_add_phase_function((*plan_ptr)->fctPlan, i, in->phase_fn[i-1],
                               (*plan_ptr)->status);
    }

    (*plan_ptr)->num_data_cubes = 0;
    (*plan_ptr)->fct_initialised = 0;
}

void LALDestroyFCTPlan(LALStatus* const status,
                       LALFCTPlan** plan_ptr)
{
    INITSTATUS(status, "LALDestroyFCTPlan", LALFCTINTERFACEC);

    ASSERT(*plan_ptr != 0, status,
           LALFCTINTERFACEH_ENULL, LALFCTINTERFACEH_MSGENULL);

    fct_destroy_plan((*plan_ptr)->fctPlan, (*plan_ptr)->status);

    ASSERT((*plan_ptr)->status->fct_errno == 0, status,
     LALFCTINTERFACEH_EINTERNAL, fct_strerror((*plan_ptr)->status->fct_errno));

    LALFree((*plan_ptr)->status);
    (*plan_ptr)->status = 0;

    LALFree(*plan_ptr);
    *plan_ptr = 0;
}

void LALFCTSetUnits(LALStatus* const status,
                    const LALFCTSetUnitsInput* const in,
                    LALFCTPlan* const plan)
{
    INITSTATUS(status, "LALFCTSetUnits", LALFCTINTERFACEC);

    ASSERT(plan != 0, status,
           LALFCTINTERFACEH_ENULL, LALFCTINTERFACEH_MSGENULL);
    ASSERT(in != 0, status,
           LALFCTINTERFACEH_ENULL, LALFCTINTERFACEH_MSGENULL);

    fct_set_units(plan->fctPlan, in->offset, in->delta, plan->status);

    ASSERT(plan->status->fct_errno == 0, status,
           LALFCTINTERFACEH_EINTERNAL, fct_strerror(plan->status->fct_errno));

    plan->fct_initialised = 0;
}

void LALFCTSetMaxSegments(LALStatus* const status,
                          const INT4 max_segments,
                          LALFCTPlan* const plan)
{
    INITSTATUS(status, "LALFCTSetMaxSegments", LALFCTINTERFACEC);

    ASSERT(plan != 0, status,
           LALFCTINTERFACEH_ENULL, LALFCTINTERFACEH_MSGENULL);
    ASSERT(max_segments > 0, status,
           LALFCTINTERFACEH_ENPOS, LALFCTINTERFACEH_MSGENPOS);

    fct_set_max_segments(plan->fctPlan, max_segments, plan->status);

    ASSERT(plan->status->fct_errno == 0, status,
           LALFCTINTERFACEH_EINTERNAL, fct_strerror(plan->status->fct_errno));

    plan->fct_initialised = 0;
}

void LALFCTSetDataCubes(LALStatus* const status,
                       const LALFCTSetDataCubesInput* const in,
                       LALFCTPlan* const plan)
{
    UINT4 i = 0;

    INITSTATUS(status, "LALFCTSetDataCubes", LALFCTINTERFACEC);

    ASSERT(plan != 0, status,
           LALFCTINTERFACEH_ENULL, LALFCTINTERFACEH_MSGENULL);
    ASSERT(in != 0, status,
           LALFCTINTERFACEH_ENULL, LALFCTINTERFACEH_MSGENULL);
    ASSERT(in->data_cube != 0, status,
           LALFCTINTERFACEH_ENULL, LALFCTINTERFACEH_MSGENULL);
    ASSERT(in->num_data_cubes > 0, status,
           LALFCTINTERFACEH_ENPOS, LALFCTINTERFACEH_MSGENPOS);

    for (i = 0; i < plan->num_data_cubes; ++i)
    {
        INT4 j = 0;
        for (j = 0; j < plan->fctPlan->number_of_dimensions; ++j)
        {
            ASSERT(in->data_cube[i].start_locations[j] >= 0, status,
                   LALFCTINTERFACEH_ECUBE, LALFCTINTERFACEH_MSGECUBE);
            ASSERT(in->data_cube[i].end_locations[j]
                       > in->data_cube[i].start_locations[j], status,
                   LALFCTINTERFACEH_ECUBE, LALFCTINTERFACEH_MSGECUBE);
            ASSERT(in->data_cube[i].stride[j] > 0, status,
                   LALFCTINTERFACEH_ENPOS, LALFCTINTERFACEH_MSGENPOS);
        }
    }

    for (i = 0; i < plan->num_data_cubes; ++i)
    {
        fct_remove_data_cube(plan->fctPlan, plan->status);
        ASSERT(plan->status->fct_errno == 0, status,
            LALFCTINTERFACEH_EINTERNAL, fct_strerror(plan->status->fct_errno));
    }

    plan->num_data_cubes = in->num_data_cubes;

    for (i = 0; i < plan->num_data_cubes; ++i)
    {
        fct_add_data_cube(plan->fctPlan,
                          in->data_cube[i].start_locations,
                          in->data_cube[i].end_locations,
                          in->data_cube[i].stride,
                          0 /* mask */,
                          FCT_SPECIFY_RANGES /* mode */,
                          plan->status);
        ASSERT(plan->status->fct_errno == 0, status,
            LALFCTINTERFACEH_EINTERNAL, fct_strerror(plan->status->fct_errno));
    }

    plan->fct_initialised = 0;
}

void LALFCTOutputDataSize(LALStatus* const status,
                          UINT8* const output_data_size,
                          LALFCTPlan* const plan)
{
    INITSTATUS(status, "LALFCTOutputDataSize", LALFCTINTERFACEC);

    ASSERT(plan != 0, status,
           LALFCTINTERFACEH_ENULL, LALFCTINTERFACEH_MSGENULL);
    ASSERT(output_data_size != 0, status,
           LALFCTINTERFACEH_ENULL, LALFCTINTERFACEH_MSGENULL);
    ASSERT(plan->num_data_cubes > 0, status,
           LALFCTINTERFACEH_ENPOS, LALFCTINTERFACEH_MSGENPOS);

#ifndef FCT_INIT_DEPRECATED
    if (plan->fct_initialised == 0)
    {
        fct_init(plan->fctPlan);
        plan->fct_initialised = 1;
    }
#endif
    
    /*
      Divide by two since the output vector in LAL is an array
      of complex numbers but the output vector in FCT is an
      array of real numbers (real/imag components)
    */
    *output_data_size = fct_output_data_size(plan->fctPlan, plan->status)/2;
    ASSERT(plan->status->fct_errno == 0, status,
           LALFCTINTERFACEH_EINTERNAL, fct_strerror(plan->status->fct_errno));
}

void LALFCTCalculate(LALStatus* const status,
                     COMPLEX8Vector* const out,
                     const COMPLEX8Vector* const in,
                     LALFCTPlan* const plan)
{
    UINT8 output_data_size = 0;

    INITSTATUS(status, "LALFCTCalculate", LALFCTINTERFACEC);

    ASSERT(plan != 0, status,
           LALFCTINTERFACEH_ENULL, LALFCTINTERFACEH_MSGENULL);
    ASSERT(in != 0, status,
           LALFCTINTERFACEH_ENULL, LALFCTINTERFACEH_MSGENULL);
    ASSERT(out != 0, status,
           LALFCTINTERFACEH_ENULL, LALFCTINTERFACEH_MSGENULL);

    ASSERT((INT8)(in->length) == (INT8)(plan->fctPlan->data_length), status,
           LALFCTINTERFACEH_ESIZE, LALFCTINTERFACEH_MSGESIZE);

    ASSERT(plan->num_data_cubes > 0, status,
           LALFCTINTERFACEH_ENPOS, LALFCTINTERFACEH_MSGENPOS);

#ifndef FCT_INIT_DEPRECATED
    if (plan->fct_initialised == 0)
    {
        fct_init(plan->fctPlan);
        plan->fct_initialised = 1;
    }
#endif
    
    output_data_size = fct_output_data_size(plan->fctPlan, plan->status)/2;

    ASSERT((out->length == output_data_size), status,
           LALFCTINTERFACEH_ESZMM, LALFCTINTERFACEH_MSGESZMM);

    fct_calculate(plan->fctPlan,
                  (fct_real*) in->data,
                  (fct_real*) out->data,
                  plan->status);

    ASSERT(plan->status->fct_errno == 0, status,
           LALFCTINTERFACEH_EINTERNAL, fct_strerror(plan->status->fct_errno));
}
