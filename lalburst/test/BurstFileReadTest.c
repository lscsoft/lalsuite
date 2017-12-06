#include <lal/LALSimBurst.h> 
#include <lal/GenerateBurst.h> 
#include <lal/LIGOLwXMLBurstRead.h> 
#include <lal/SnglBurstUtils.h> 
#include <lal/TimeSeries.h>

int main(void){

    SimBurst* sim_burst;
    REAL8TimeSeries *hp, *hx;
    int retcode;
    size_t i;

    /* Get the injection */
    printf("Reading injection from file\n");
    sim_burst = XLALSimBurstTableFromLIGOLw("ad_hoc_test.xml.gz", NULL, NULL);

    /* Do the injection */
    printf("Generating injection from file\n");
    retcode = XLALGenerateSimBurst(&hp, &hx, sim_burst, 1.0/16384);

    for (i=0; i<hp->data->length; i++) {
        printf("%lu: %f %f\n", i, hp->data->data[i], hx->data->data[i]);
    }

    /* Clean up */
    XLALDestroySimBurstTable(sim_burst);
    XLALDestroyREAL8TimeSeries(hp);
    XLALDestroyREAL8TimeSeries(hx);

    return retcode;
}
