/*----------------------------------------------------------------------- 
 * 
 * File Name: DateString.c 
 * 
 * Author: David Chin <dwchin@umich.edu>
 * 
 * Revision: $Id$
 * 
 *----------------------------------------------------------------------- 
 * 
 * NAME 
 * DateString
 * 
 * SYNOPSIS
 *
 * LALDateString(): Returns a formatted string for the date and time in
 *               ISO 8601 format, given an LALDate.
 * 
 * DESCRIPTION
 *
 * LALDateString():
 *      Inputs:   (LALDate *) time
 *
 *      Outputs:  (CHAR *) formatted string date and time
 * 
 * 
 * DIAGNOSTICS 
 * (Abnormal termination conditions, error and warning codes summarized 
 * here. More complete descriptions are found in documentation.)
 *
 * CALLS
 * (list of LLAL, LDAS, other non-system functions/procedures called. 
 * 
 * NOTES
 * See: http://www.iso.ch/markete/8601.pdf -- Official ISO Document
 *      http://www.cl.cam.ac.uk/~mgk25/iso-time.html -- overview
 * 
 *-----------------------------------------------------------------------
 */
#include <lal/LALRCSID.h>

NRCSID (DATESTRINGC, "$Id$");

#include <lal/Date.h>


void
LALDateString (LALStatus     *status,
               CHARVector    *timestamp,
               const LALDate *date)
{
    INITSTATUS (status, "LALDateString", DATESTRINGC);

    /*
     * Check pointer to input variable
     */
    ASSERT (date != (LALDate *)NULL, status,
            DATESTRING_ENULLINPUT, DATESTRING_MSGENULLINPUT);

    /*
     * Check pointer to output variable.
     */
    ASSERT (timestamp != (CHARVector *)NULL, status,
            DATESTRING_ENULLOUTPUT, DATESTRING_MSGENULLOUTPUT);

    /*
     * Check that timestamp buffer is large enough
     */
    ASSERT (timestamp->length >= 26, status,
            DATESTRING_EBUFFTOOSMALL, DATESTRING_MSGEBUFFTOOSMALL);

    /*
     * Use strftime (3) to form ISO8601-format time stamp, plus day name
     */
    strftime(timestamp->data, timestamp->length,
             "%Y-%m-%d %H:%M:%S UTC %a", &(date->unixDate));

    RETURN (status);
} /* END LALDateString() */

    
