#include <lal/LALSimInspiralWaveformFlags.h>
#include <lal/LALStatusMacros.h>
#include <stdio.h>

int main(int argc , char *argv[])
{
    /* Set lalDebugLevel to print all info, warnings, errors */
    lalDebugLevel = 7;

    (void) argc;
    (void) argv;

    /* Create a new struct */
    printf("Creating new WaveformFlags with default values...\n");
    LALSimInspiralWaveformFlags *test=XLALSimInspiralCreateWaveformFlags();
    /* Print values of flags and check if default */
    printf("Get spinO: %i\n",XLALSimInspiralGetSpinOrder(test));
    printf("spinO default? %i\n",XLALSimInspiralSpinOrderIsDefault(XLALSimInspiralGetSpinOrder(test)));
    printf("Get tideO: %i\n",XLALSimInspiralGetTidalOrder(test));
    printf("tideO default? %i\n",XLALSimInspiralTidalOrderIsDefault(XLALSimInspiralGetTidalOrder(test)));
    printf("Get FrameAxis: %i\n",XLALSimInspiralGetFrameAxis(test));
    printf("FrameAxis default? %i\n",XLALSimInspiralFrameAxisIsDefault(XLALSimInspiralGetFrameAxis(test)));
    printf("Get ModesChoice: %i\n",XLALSimInspiralGetModesChoice(test));
    printf("ModesChoice default? %i\n",XLALSimInspiralModesChoiceIsDefault(XLALSimInspiralGetModesChoice(test)));

    /* Change individual fields */
    printf("Set spinO=2.5PN using enum member...\n");
    XLALSimInspiralSetSpinOrder(test,LAL_SIM_INSPIRAL_SPIN_ORDER_25PN);
    printf("Get spinO: %i\n",XLALSimInspiralGetSpinOrder(test));
    printf("Set spinO=2PN using integer...\n");
    XLALSimInspiralSetSpinOrder(test,4);
    printf("Get spinO: %i\n",XLALSimInspiralGetSpinOrder(test));
    printf("Set tideO=5PN using integer...\n");
    XLALSimInspiralSetTidalOrder(test,10);
    printf("Get tideO: %i\n",XLALSimInspiralGetTidalOrder(test));
    printf("Set tideO=5PN using enum member...\n");
    XLALSimInspiralSetTidalOrder(test,LAL_SIM_INSPIRAL_TIDAL_ORDER_6PN);
    printf("Get tideO: %i\n",XLALSimInspiralGetTidalOrder(test));

    /* Recheck for default values */
    printf("spinO default? %i\n",XLALSimInspiralSpinOrderIsDefault(XLALSimInspiralGetSpinOrder(test)));
    printf("tideO default? %i\n",XLALSimInspiralTidalOrderIsDefault(XLALSimInspiralGetTidalOrder(test)));
    printf("FrameAxis default? %i\n",XLALSimInspiralFrameAxisIsDefault(XLALSimInspiralGetFrameAxis(test)));
    printf("ModesChoice default? %i\n",XLALSimInspiralModesChoiceIsDefault(XLALSimInspiralGetModesChoice(test)));

    printf("Destroying LALSimInspiralWaveformFlags struct\n");
    XLALSimInspiralDestroyWaveformFlags(test);

    return 0;
}
