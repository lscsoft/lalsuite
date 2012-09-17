#include <lal/LALSimInspiralTestGRParams.h>
#include <lal/LALStatusMacros.h>

int main(int argc , char *argv[])
{
    /* Set lalDebugLevel to print all info, warnings, errors */
    lalDebugLevel = 7;

    (void) argc;
    (void) argv;
    int errnum;

    /* Create a new struct */
    printf("Creating a new struct with one parameter...\n");
    LALSimInspiralTestGRParam *test=NULL;
    if( XLALSimInspiralAddTestGRParam(&test,"first_param",10.0) != XLAL_SUCCESS)
        XLAL_ERROR(XLAL_EFUNC);
    printf("value:%lf\n",XLALSimInspiralGetTestGRParam(test, "first_param"));

    /* Add a second parameter to the struct */
    printf("Adding a second parameter to the struct...\n");
    if( XLALSimInspiralAddTestGRParam(&test,"second_param",20.0) != XLAL_SUCCESS)
        XLAL_ERROR(XLAL_EFUNC);
    printf("Printing the struct after adding second parameter...\n");
    if( XLALSimInspiralPrintTestGRParam(stderr,test) != XLAL_SUCCESS)
        XLAL_ERROR(XLAL_EFUNC);
    printf("Get second parameter and print its value...\n");
    printf("value:%lf\n",XLALSimInspiralGetTestGRParam(test, "second_param"));

    /* Set first parameter to new value */
    printf("Setting first parameter to new value...\n");
    if( XLALSimInspiralSetTestGRParam(test, "first_param",1000.0)
            != XLAL_SUCCESS )
        XLAL_ERROR(XLAL_EFUNC);
    printf("new value:%lf\n",
            XLALSimInspiralGetTestGRParam(test, "first_param"));

    /* Check for existing and non-existent parameters */
    printf("Checking for existing first parameter and non-existent third parameter...\n");
    printf("first_param:%d second_param:%d third_param:%d\n",
            XLALSimInspiralTestGRParamExists(test, "first_param"),
            XLALSimInspiralTestGRParamExists(test, "second_param"),
            XLALSimInspiralTestGRParamExists(test, "third_param"));
    printf("Now add a third parameter...\n");
    if( XLALSimInspiralAddTestGRParam(&test,"third_param",12.0) != XLAL_SUCCESS )
        XLAL_ERROR(XLAL_EFUNC);
    printf("first_param:%d second_param:%d third_param:%d\n",
            XLALSimInspiralTestGRParamExists(test, "first_param"),
            XLALSimInspiralTestGRParamExists(test, "second_param"),
            XLALSimInspiralTestGRParamExists(test, "third_param"));

    /* Print the params struct */
    printf("We print the struct as it appears now...\n");
    if( XLALSimInspiralPrintTestGRParam(stderr,test) != XLAL_SUCCESS )
        XLAL_ERROR(XLAL_EFUNC);

    /* Try to add a parameter that already exists */
    printf("Trying to add a parameter that already exists...\n");
    XLAL_TRY( XLALSimInspiralAddTestGRParam(&test,"third_param",12.0), errnum );
    printf("This throws the above message and error code %d\n", errnum);

    /* Destroy the params struct */
    XLALSimInspiralDestroyTestGRParam(test);

    return 0;
}
