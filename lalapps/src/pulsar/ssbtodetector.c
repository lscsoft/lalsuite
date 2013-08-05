/*
*  Copyright (C) 2013 Matthew Pitkin
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

/**
 * \file
 * \ingroup pulsarApps
 * \author Matt Pitkin
 *
 * \brief Program to convert an input MJD, or GPS, time at the solar system barycentre to a GPS time at a detector
 *
 * This code will take in an MJD time in TDB, or a GPS-style time in TDB at the solar system barycentre and convert
 * it to a GPS time at a given detector for a given sky position. The input detectors can be any of the standard
 * acronyms for gravitational wave detectors, but can also include several radio telescopes (these use position
 * information from TEMPO2):
 *   - the Robert C. Byrd Greeen Bank Telescope (GBT),
 *   - the Parkes Telescope (PKS),
 *   - the Lovell Telescope at Jodrell Bank (JBO),
 *   - the Arecibo Telescope (AO),
 *   - the Effelsberg 100m Radio Telescope (EFF),
 *   - the Nancay Decimetre Radio Telescope (NRT),
 *   - the Mount Pleasant Radio Observatory, Hobart (HOB),
 *   - the Hartebeesthoek Radio Astronomy Observatory (HART),
 *   - the Very Large Array (VLA),
 *   - and the Westerbork Synthesis Radio Telescope (WSRT).
 */

#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>

#include <lal/LALInitBarycenter.h>
#include <lal/GeneratePulsarSignal.h>
#include <lal/BinaryPulsarTiming.h>
#include <lal/SFTutils.h>
#include <lal/LALString.h>

/* define the radio telescope positions */
typedef enum{
  GBT = 0, /* Robert C. Byrd GBT */
  PKS,     /* Parkes */
  JBO,     /* JBO Lovell telescope */
  AO,      /* Arecibo Telescope */
  EFF,     /* Effelsberg */
  NRT,     /* Nancay */
  HOB,     /* Hobart */
  HART,    /* Hartebeesthoek */
  VLA,     /* Very Large Array */
  WSRT,    /* Westerbork */
  NUMSCOPES
} Scopes;

/* create some detector locations for radio telescopes - positions taken from TEMPO2 observatories.dat file */
REAL8 scopelocations[NUMSCOPES][3] = { { 882589.65, -4924872.32, 3943729.348 },     /* GBT */
                                       { -4554231.5, 2816759.1, -3454036.3 },       /* Parkes */
                                       { 3822626.04, -154105.65, 5086486.04 },      /* JBO */
                                       { 2390490.0, -5564764.0, 1994727.0 },        /* Arecibo Telescope */
                                       { 4033949.5, 486989.4, 4900430.8 },          /* Effelsberg */
                                       { 4324165.81, 165927.11, 4670132.83},        /* Nancay */
                                       { -3950077.96, 2522377.31, -4311667.52 },    /* Hobart */
                                       { 5085442.780, 2668263.483, -2768697.034 },  /* Hartebeesthoek */
                                       { -1601192.0, -5041981.4, 3554871.4 },       /* Very Large Array */
                                       { 3828445.659, 445223.600, 5064921.5677 } }; /* Westerbork */
#define USAGE \
"Usage: %s [options]\n\n"\
" --help (-h)              display this message\n"\
" --mjd (-m)               a MJD time in TDB at the solar system barycenter\n"\
" --gps (-g)               a GPS equivalent time in TDB at the solar system barycenter\n"\
" --ra (-r)                the source right ascension (e.g. 15:21:34.76)\n"\
" --dec (-d)               the source declination (e.g. -09:53:12.36)\n"\
" --telescope (-t)         a detector acronym:\n\
                            GW detectors:\n\
                              H1/H2 - LIGO Hanford,\n\
                              L1 - LIGO Livingston,\n\
                              V1 - Virgo,\n\
                              G1 - GEO600,\n\
                              T1 - TAMA300,\n\
                            Radio Telescopes:\n\
                              GBT - the Robert C. Byrd Green Bank Telescope,\n\
                              PKS - the Parkes Telescope,\n\
                              JBO - the Lovell Telescope at Jodrell Bank,\n\
                              AO - the Arecibo Telescope,\n\
                              EFF - the Effelsberg 100m Radio Telescope,\n\
                              NRT - the Nancay Decimetre Radio Telescope,\n\
                              HOB - the Mount Pleasant Radio Observatory, Hobart,\n\
                              HART - the Hartebeesthoek Radio Astronomy Observatory,\n\
                              VLA - the Very Large Array,\n\
                              WSRT - the Westerbork Synthesis Radio Telescope.\n"\
"\n"

int main( int argc, char *argv[] ){
  double mjdtime = 0.;
  char *det = NULL, *ra = NULL, *dec = NULL;
  LIGOTimeGPS gpsin;
  LIGOTimeGPS gpsout;
  double gpstime = 0.;
  char earth[1024], sun[1024];
  PulsarSignalParams params;
  EphemerisData *edat;
  Scopes sn = NUMSCOPES;

  struct option long_options[] =
  {
    { "help",      no_argument,       0, 'h' },
    { "telescope", required_argument, 0, 't' },
    { "ra",        required_argument, 0, 'r' },
    { "dec",       required_argument, 0, 'd' },
    { "mjd",       required_argument, 0, 'm' },
    { "gps",       required_argument, 0, 'g' },
    { 0, 0, 0, 0 }
  };

  CHAR args[] = "ht:r:d:m:";
  CHAR *program = argv[0];

  /* get input arguments */
  while(1){
    int option_index = 0;
    int c;

    c = getopt_long( argc, argv, args, long_options, &option_index );
    if ( c == -1 ) /* end of options */
      break;

    switch(c){
      case 0: /* if option set a flag, nothing else to do */
        if ( long_options[option_index].flag )
          break;
        else
          fprintf(stderr, "Error parsing option %s with argument %s\n", long_options[option_index].name, optarg );
      case 'h': /* help message */
        fprintf(stderr, USAGE, program);
        exit(0);
      case 't': /* the detector */
        det = XLALStringDuplicate( optarg );
        break;
      case 'r': /* the right ascension */
        ra = XLALStringDuplicate( optarg );
        break;
      case 'd': /* the declination */
        dec = XLALStringDuplicate( optarg );
        break;
      case 'm': /* the mjd time */
        mjdtime = atof(optarg);
        break;
      case 'g': /* the gps time */
        gpstime = atof(optarg);
      case '?':
        fprintf(stderr, "Unknown error while parsing options\n" );
        exit(0);
      default:
        fprintf(stderr, "Unknown error while parsing options\n" );
        exit(0);
    }
  }

  if ( mjdtime <= 0. && gpstime <= 0. ){
    fprintf(stderr, "Error... input MJD or GPS time is not sensible!\n");
    exit(1);
  }

  if (  mjdtime > 0. && gpstime > 0. ){
    fprintf(stderr, "Error... required either an input MJD to or an input GPS time!\n");
    exit(1);
  }

  if ( ra == NULL ){
    fprintf(stderr, "Error... no right ascension given!\n");
    exit(1);
  }

  if ( dec == NULL ){
    fprintf(stderr, "Error... no declination given!\n");
    exit(1);
  }

  if ( det == NULL ){
    fprintf(stderr, "Error... no telescope has been given!\n");
    exit(1);
  }

  /* set detector/telescope */
  if( !strcmp( det, "GBT" ) )      { sn = GBT; }
  else if( !strcmp( det, "PKS" ) ) { sn = PKS; }
  else if( !strcmp( det, "JBO" ) ) { sn = JBO; }
  else if( !strcmp( det, "AO" ) )  { sn = AO; }
  else if( !strcmp( det, "EFF" ) ) { sn = EFF; }
  else if( !strcmp( det, "NRT" ) ) { sn = NRT; }
  else if( !strcmp( det, "HOB" ) ) { sn = HOB; }
  else if( !strcmp( det, "HART" ) ){ sn = HART; }
  else if( !strcmp( det, "VLA" ) ) { sn = VLA; }
  else if( !strcmp( det, "WSRT" ) ){ sn = WSRT; }

  if( sn != NUMSCOPES ){
    LALDetector *site = NULL;
    site = (LALDetector *)XLALMalloc( sizeof(LALDetector) );
    //memcpy(site->location, scopelocations[sn], 3*sizeof(REAL8));
    site->location[0] = scopelocations[sn][0];
    site->location[1] = scopelocations[sn][1];
    site->location[2] = scopelocations[sn][2];
    params.site = site;
  }
  else{ params.site = XLALGetSiteInfo( det ); } /* try a GW detector */

  params.pulsar.position.latitude = XLALdmsToRads( dec );
  params.pulsar.position.longitude = XLALhmsToRads( ra );
  params.pulsar.position.system = COORDINATESYSTEM_EQUATORIAL;

  /* get the Earth and Sun ephemeris files - note yo may have to change this for different systems */
  sprintf(earth, "%s/share/lalpulsar/earth00-19-DE405.dat.gz", getenv("LALPULSAR_PREFIX"));
  sprintf(sun, "%s/share/lalpulsar/sun00-19-DE405.dat.gz", getenv("LALPULSAR_PREFIX"));

  /* double check that the files exist */
  if( fopen(sun, "r") == NULL || fopen(earth, "r") == NULL ){
    fprintf(stderr, "Error... ephemeris files not, or incorrectly, defined!\n");
    exit(1);
  }

  /* set up ephemeris files */
  edat = XLALInitBarycenter( earth, sun );
  params.ephemerides = edat;

  /* convert MJD time to GPS format (remaining at the SSB) */
  if ( mjdtime > 0. ) { gpstime = XLALTTMJDtoGPS( mjdtime ); }

  int ephemstart = 630720013; /* GPS time of Jan 1, 2000, 00:00:00 UTC */
  int ephemend = 1261872015;  /* GPS time of Jan 1, 2020, 00:00:00 UTC */

  if( gpstime < ephemstart || gpstime > ephemend ){
    fprintf(stderr, "Time (GPS %.9lf) is outside the ephemeris file ranges!\n", gpstime);
    exit(1);
  }

  /* put into LIGOTimeGPS format */
  XLALGPSSetREAL8( &gpsin, gpstime );

  /* convert time at SSB to GPS time at detector */
  if ( XLALConvertSSB2GPS( &gpsout, gpsin, &params ) != XLAL_SUCCESS ){
    fprintf(stderr, "Problem converting time!\n");
    exit(1);
  }

  //fprintf(stdout, "%.9lf\n", XLALGPSGetREAL8( &gpsout ) );
  fprintf(stdout, "%d.%09d\n", gpsout.gpsSeconds, gpsout.gpsNanoSeconds );

  return 0;
}
