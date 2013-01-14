/*
 * Copyright (C) 2008 Jordi Burguet-Castell, Xavier Siemens.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with with program; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 * MA 02111-1307 USA
 */

/*
 * Read the filters file. To be used for calibration.
 */

#include <complex.h>
#include <lal/ConfigFile.h>   /* to use LALParseDataFile() */
#include <lal/AVFactories.h>  /* to use XLALCreateREAL8Vector() */

#include <lal/ReadFiltersFile.h>


static int check_checksum(const char *filterfile, const char *claimed_chksum);

/*
 * A common check in the code when reading the filters file.
 */
#define CHECK(VAR, NAME)                                                \
    do {                                                                \
        if (strcmp(VAR, NAME) != 0) {                                   \
            fprintf(stderr,                                             \
                    "ERROR: Line (%s) of file %s is not properly terminated " \
                    "by '%s' marker!\n\n", thisline, filterfile, NAME); \
            XLALDestroyParsedDataFile(Filters);                         \
            return -1;                                                  \
        }                                                               \
    } while (0)



/*
 * Read the filters file and fill the appropiate structures.
 *
 * InputData must be allocated before calling this function. After
 * calling this function, many of its values will be initialized.
 */
int XLALReadFiltersFile(const char *filterfile, StrainIn *InputData)
{
    LALParsedDataFile *Filters = NULL;  /* preparsed contents of Filters-file */
    int numlines, i;  /* total number of lines and line counter */
    int n, l;         /* counters */
    CHAR *thisline;
    char paramname[64];         /* name of the parameter */
    char sensingstr[8], usfstr[18], delaystr[6];  /* filters labels */
    char aastr[10], servostr[6], awstr[14];
    int NCinv, NA, ND, NAW;     /* number of points in filter */
    int err = 0;  /* error code */
    REAL8 re, im;

    err = XLALParseDataFile(&Filters, filterfile);
    if (err) {
        fprintf(stderr, "Error parsing data file %s\n", filterfile);
        return err;
    }

    numlines = Filters->lines->nTokens; /* how many lines of data */

    /* Check that file is not empty */
    if (numlines == 0) {
        fprintf(stderr, "File %s has no contents!\n", filterfile);
        XLALDestroyParsedDataFile(Filters);
        return -2;
    }

    /**------------------------------------------------------------------**/
    /* Read checksum */
    i = 0;                                  /* start with first line */
    thisline = Filters->lines->tokens[i];   /* get line i */

    if (strncmp(thisline, "SHA-1 checksum", 14) != 0) {
        fprintf(stderr, "ERROR: file %s does not contain a checksum.\n",
                filterfile);
        XLALDestroyParsedDataFile(Filters);
        return -1;
    }

    {
        char *claimed_chksum = (char *) &thisline[16];
        claimed_chksum[40] = '\0';

        strncpy(InputData->filter_chksum, claimed_chksum,
                sizeof InputData->filter_chksum);
        if (check_checksum(filterfile, claimed_chksum) != 0) {
            fprintf(stderr, "ERROR: bad checksum for file %s\n", filterfile);
            XLALDestroyParsedDataFile(Filters);
            return -1;
        }
    }

    /**------------------------------------------------------------------**/
    /* Read VC info (CVS information for the moment) */
    thisline = Filters->lines->tokens[++i];   /* get next line */
    strncpy(InputData->filter_vc_info, thisline, sizeof(InputData->filter_vc_info));

    /**------------------------------------------------------------------**/
    /* Read parameters that used to be passed manually */
    /* (Do, Go, Wo, f) */

    /* Frequency of the calibration line */
    thisline = Filters->lines->tokens[++i];   /* get next line */
    sscanf(thisline, "%" LAL_REAL8_FORMAT " %s", &InputData->f, paramname);
    CHECK(paramname, "CAL_LINE_FREQ");

    /* Open loop gain at cal line freq */
    thisline = Filters->lines->tokens[++i];   /* get next line */
    sscanf(thisline, "%" LAL_REAL8_FORMAT " %s", &re, paramname);
    CHECK(paramname, "OLG_RE");

    thisline = Filters->lines->tokens[++i];   /* get next line */
    sscanf(thisline, "%" LAL_REAL8_FORMAT " %s", &im, paramname);
    CHECK(paramname, "OLG_IM");

    InputData->Go = re + I * im;

    /* Whitening filter at cal line freq */
    thisline = Filters->lines->tokens[++i];   /* get next line */
    sscanf(thisline, "%" LAL_REAL8_FORMAT " %s", &re, paramname);
    CHECK(paramname, "WHITENER_RE");

    thisline = Filters->lines->tokens[++i];   /* get next line */
    sscanf(thisline, "%" LAL_REAL8_FORMAT " %s", &im, paramname);
    CHECK(paramname, "WHITENER_IM");

    InputData->Wo = re + I * im;

    /* Digital filter (servo) at cal line freq */
    thisline = Filters->lines->tokens[++i];   /* get next line */
    sscanf(thisline, "%" LAL_REAL8_FORMAT " %s", &re, paramname);
    CHECK(paramname, "SERVO_RE");

    thisline = Filters->lines->tokens[++i];   /* get next line */
    sscanf(thisline, "%" LAL_REAL8_FORMAT " %s", &im, paramname);
    CHECK(paramname, "SERVO_IM");

    InputData->Do = re + I * im;

    /**------------------------------------------------------------------**/
    /* Read sensing function info */
    thisline = Filters->lines->tokens[++i];   /* get next line */
    sscanf(thisline, "%s", sensingstr);
    CHECK(sensingstr, "SENSING");

    thisline = Filters->lines->tokens[++i];   /* get next line */
    sscanf(thisline, "%" LAL_INT4_FORMAT " %s", &InputData->CinvUSF, usfstr);
    CHECK(usfstr, "UPSAMPLING_FACTOR");

    /* Check upsamplig factor USF = 1, or positive and mod 2=0 */
    if (! (InputData->CinvUSF == 1 ||
           (InputData->CinvUSF > 1 && InputData->CinvUSF % 2 == 0)) ) {
        fprintf(stderr, "ERROR: bad upsampling factor %d", InputData->CinvUSF);
        XLALDestroyParsedDataFile(Filters);
        return -1;
    }

    /* Read Delay */
    thisline = Filters->lines->tokens[++i];   /* get next line */
    sscanf(thisline, "%" LAL_INT4_FORMAT " %s", &InputData->CinvDelay, delaystr);
    CHECK(delaystr, "DELAY");

    /* Read number of sensing filters and their orders */
    thisline = Filters->lines->tokens[++i];   /* get next line */
    NCinv = strtol(thisline, &thisline, 10);
    CHECK(thisline, " FILTER_ORDER");

    /* Allocate inverse sensing funtion filters */
    InputData->Cinv = (REAL8IIRFilter *)LALMalloc(sizeof(REAL8IIRFilter));

    /* Allocate inverse sensing function filter */
    InputData->Cinv->directCoef = XLALCreateREAL8Vector(NCinv);
    InputData->Cinv->recursCoef = XLALCreateREAL8Vector(NCinv);
    InputData->Cinv->history    = XLALCreateREAL8Vector(NCinv-1);

    for (l=0; l < NCinv; l++)    InputData->Cinv->directCoef->data[l] = 0.0;
    for (l=0; l < NCinv; l++)    InputData->Cinv->recursCoef->data[l] = 0.0;
    for (l=0; l < NCinv-1; l++)  InputData->Cinv->history->data[l]    = 0.0;

    for (n = 0; n < NCinv; n++) {   /* read direct coeffs */
        thisline = Filters->lines->tokens[++i];   /* get next line */
        InputData->Cinv->directCoef->data[n] = strtod(thisline, NULL);
    }

    /**--------------------------------------------------------------------**/
    /* Read in servo */
    thisline = Filters->lines->tokens[++i];   /* get next line */
    sscanf(thisline, "%s", servostr);
    CHECK(servostr, "SERVO" );

    /* Read number of sensing filters and their orders */
    thisline = Filters->lines->tokens[++i];   /* get next line */
    ND = strtol(thisline, &thisline, 10);
    CHECK(thisline, " FILTER_ORDER");

    /* Allocate inverse sensing funtion filters */
    InputData->D = (REAL8IIRFilter *)LALMalloc(sizeof(REAL8IIRFilter));

    /* Allocate inverse sensing function filter */
    InputData->D->directCoef = XLALCreateREAL8Vector(ND);
    InputData->D->recursCoef = XLALCreateREAL8Vector(ND);
    InputData->D->history    = XLALCreateREAL8Vector(ND-1);

    for (l=0; l < ND; l++)    InputData->D->directCoef->data[l] = 0.0;
    for (l=0; l < ND; l++)    InputData->D->recursCoef->data[l] = 0.0;
    for (l=0; l < ND-1; l++)  InputData->D->history->data[l]    = 0.0;

    for (n = 0; n < ND; n++) {   /* read direct coeffs */
        thisline = Filters->lines->tokens[++i];   /* get next line */
        InputData->D->directCoef->data[n] = strtod(thisline, NULL);
    }

    /**--------------------------------------------------------------------**/
    /* Read in actuation */
    thisline = Filters->lines->tokens[++i];
    sscanf(thisline, "%s", aastr);
    CHECK(aastr, "ACTUATION");

    /* Read number of sensing filters and their orders */
    thisline = Filters->lines->tokens[++i];
    NA = strtol(thisline, &thisline, 10);
    CHECK(thisline, " FILTER_ORDER");

    /* Allocate inverse sensing funtion filters */
    InputData->A = (REAL8IIRFilter *)LALMalloc(sizeof(REAL8IIRFilter));

    /* Allocate inverse sensing function filter */
    InputData->A->directCoef = XLALCreateREAL8Vector(NA);
    InputData->A->recursCoef = XLALCreateREAL8Vector(NA);
    InputData->A->history    = XLALCreateREAL8Vector(NA-1);

    for (l=0; l < NA; l++)    InputData->A->directCoef->data[l] = 0.0;
    for (l=0; l < NA; l++)    InputData->A->recursCoef->data[l] = 0.0;
    for (l=0; l < NA-1; l++)  InputData->A->history->data[l]    = 0.0;

    for (n = 0; n < NA; n++) {   /* read direct coeffs */
        thisline = Filters->lines->tokens[++i];
        InputData->A->directCoef->data[n] = strtod(thisline, NULL);
    }

    /**-------------------------------------------------------------------**/
    /* Read in antiwhitening filter */
    thisline = Filters->lines->tokens[++i];
    sscanf(thisline,"%s", awstr);
    CHECK(awstr, "ANTIWHITENING");

    /* Read number of sensing filters and their orders */
    thisline = Filters->lines->tokens[++i];
    NAW = strtol(thisline, &thisline, 10);
    CHECK(thisline, " FILTER_ORDER");

    /* Allocate inverse sensing funtion filters */
    InputData->AW = (REAL8IIRFilter *)LALMalloc(sizeof(REAL8IIRFilter));

    /* Allocate inverse sensing function filter */
    InputData->AW->directCoef = XLALCreateREAL8Vector(NAW);
    InputData->AW->recursCoef = XLALCreateREAL8Vector(NAW);
    InputData->AW->history    = XLALCreateREAL8Vector(NAW-1);

    for (l=0; l < NAW; l++)    InputData->AW->directCoef->data[l] = 0.0;
    for (l=0; l < NAW; l++)    InputData->AW->recursCoef->data[l] = 0.0;
    for (l=0; l < NAW-1; l++)  InputData->AW->history->data[l]    = 0.0;

    for (n = 0; n < NAW; n++) {   /* read direct coeffs */
        thisline = Filters->lines->tokens[++i];   /* get next line */
        InputData->AW->directCoef->data[n] = strtod(thisline, NULL);
    }

    /**--------------------------------------------------------------------**/
    XLALDestroyParsedDataFile(Filters);

    return 0;
}


/*
 * Free memory reserved for the filters.
 */
void XLALDestroyFiltersFile(StrainIn* InputData)
{
    int i;

    REAL8IIRFilter *filters[4];

    filters[0] = InputData->Cinv;
    filters[1] = InputData->A;
    filters[2] = InputData->AW;
    filters[3] = InputData->D;

    for (i = 0; i < 4; i++) {
        XLALDestroyREAL8Vector(filters[i]->directCoef);
        XLALDestroyREAL8Vector(filters[i]->recursCoef);
        XLALDestroyREAL8Vector(filters[i]->history);
        LALFree(filters[i]);
    }

}


/*
 * Check the filters file checksum.
 */
static int check_checksum(const char *filterfile, const char *claimed_chksum)
{
    FILE *fpipe;
    char command[1024];
    char actual_chksum[1024];

    /* Actual checksum */
    snprintf(command, sizeof command, "tail -n +3 %s | sha1sum", filterfile);
    if ( !(fpipe = (FILE*) popen(command, "r")) ) {
        perror("Problem in popen when checking checksum.");
        return -1;
    }
    if ( !(fgets(actual_chksum, sizeof actual_chksum, fpipe)) ) {
        perror("Problem in fgets when checking checksum.");
        return -1;
    }
    pclose(fpipe);

    /* Compare! */
    return strncmp(claimed_chksum, actual_chksum, 40);
}
