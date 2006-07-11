/*
 * stochastic.h - SGWB Standalone Analysis Pipeline
 *                      - header file
 *
 * Tania Regimbau <regimbau@obs-nice>
 *
 * $Id$
 */

#ifndef _STOCHASTIC2_H
#define _STOCHASTIC2_H

#ifdef  __cplusplus
extern "C" {
#endif

NRCSID (STOCHASTIC2H, "$Id$" );


void parseOptions(INT4 argc, CHAR *argv[]);
void displayUsage(INT4 exitcode);

#ifdef  __cplusplus
}
#endif

#endif /* _STOCHASTIC_H */
