/*
*  Copyright (C) 2007 Jolien Creighton, Teviet Creighton
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
\author Creighton, T. D.
\file
\ingroup pulsarTODO

\heading{Module \ref PulsarCatInput.c}

Parses a catalogue of pulsar data.

\heading{Prototypes}





\heading{Description}

The routine <tt>LALReadPulsarCatHead()</tt> takes a set of tokens
<tt>*list</tt> (as generated from an input line by
<tt>LALCreateTokenList()</tt>, and determines which tokens correspond to
the pulsar catalogue fields enumerated by <tt>indx[]</tt>: each element
of <tt>indx[]</tt> stores the number of the corresponding token.  Each
token should be a field name corresponding to one of the enumeration
constants in \c PulsarCatIndex (e.g.\ <tt>"RAJ"</tt>,
<tt>"POSEPOCH"</tt>, <tt>"F"</tt>), or <tt>"e"</tt> to represent an
uncertainty in the preceding field.  Unrecognized tokens are
ignored. If any pusar catalogue field does not have a corresponding
token in <tt>indx[]</tt> is set to \f$-1\f$.

The routine <tt>LALReadPulsarCatLine()</tt> takes the input
<tt>*line</tt>, splits it into whitespace-separated tokens, and parses
each token into the corresponding field of <tt>*node</tt>, using
<tt>indx[]</tt> to determine which tokens correspond to which fields.
In general, each field has a distinct parsing algorithm, as specified
below:
<dl>
<dt>NAME</dt><dd> A standard B1950 or J2000 pulsar name (e.g.\
<tt>"B0021-72C"</tt>, <tt>"J0024-7203U"</tt>), copied directly into
<tt>node->bname</tt> or <tt>node->jname</tt>.</dd>
<dt>RAJ</dt><dd> J2000 right ascencion in the form
<tt>"</tt>hours\c :minutes\c :seconds<tt>"</tt>, where hours is a
signed integer, minutes an unsigned integer, and seconds is an
unsigned floating-point number in normal place-index notation (i.e.\
integral part, optional decimal place, and optional fractional part;
no exponential notation).</dd>
<dt>DECJ</dt><dd> J2000 declination in the form
<tt>"</tt>degrees\c :minutes\c :seconds<tt>"</tt>, where degrees
is a signed integer, minutes an unsigned integer, and seconds is an
unsigned floating-point number in normal place-index notation (i.e.\
integral part, optional decimal place, and optional fractional part;
no exponential notation).</dd>
<dt>PMRA</dt><dd> Right ascension component of proper motion in
milliarcseconds per year, as a floating-point number (any notation).</dd>
<dt>PMDEC</dt><dd> Declination component of proper motion in
milliarcseconds per year, as a floating-point number (any notation).</dd>
<dt>POSEPOCH</dt><dd> Epoch of position/proper motion measurements
in Julian days, as a floating-point number (any notation).  If the
number is less than 2~million, then it is assumed that the actual
Julian day is 2~million plus the number given.</dd>
<dt>F</dt><dd> The pulsar spin frequency in Hz, as a floating-point
number (any notation).</dd>
<dt>F1</dt><dd> The first derivative of the pulsar spin frequency
in Hz\f${}^2\f$, as a floating-point number (any notation).</dd>
<dt>F2</dt><dd> The pulsar spin frequency in Hz\f${}^3\f$, as a
floating-point number (any notation).</dd>
<dt>PEPOCH</dt><dd> Epoch of frequency and frequency-derivative
measurements in Julian days, as a floating-point number (any
notation).  If the number is less than 2~million, then it is assumed
that the actual Julian day is 2~million plus the number given.</dd>
<dt>e</dt><dd> Uncertainty in any of the preceding quantities.
This is given as an unsigned integer corresponding to the uncertainty
in the last 1 or 2 significant digits of the corresponding quantity.
Thus, the parsing routine for that quantity is also responsible for
reading its uncertainty, accounting for the number of significant
digits.</dd>
</dl>
An asterisk <tt>*</tt> in any field means that the quantity is not
measured.  In most cases this means it will be treated as zero.

\heading{Uses}
\code
lalDebugLevel
LALWarning()                  LALStringToU2()
XLALGPSLeapSeconds()
LALDCreateVector()            LALDDestroyVector()
\endcode

\heading{Notes}

At present, <tt>LALLeapSecs()</tt> fails for times prior to the GPS
epoch (1980-01-06 00:00:00, or JD2444244.5000), which prevents earlier
Julian dates from being converted into \c LIGOTimeGPS structures.
Pulsar catalogue entries with epochs earlier than this date will cause
<tt>LALReadPulsarCatLine()</tt> to fail.



*/

#include <math.h>
#include <ctype.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/StringInput.h>
#include <lal/AVFactories.h>
#include <lal/Date.h>
#include <lal/LALBarycenter.h>
#include <lal/SkyCoordinates.h>
#include <lal/PulsarCat.h>

/* Local routine to parse a token into a REAL8 number, with additional
   precision information stored in codes: codes[0] is -1 for no number
   parsed, 0 for a normal representable output, 1 for a subnormal
   number, 2 for off-scale low, and 3 for off-scale high.  codes[1] is
   the total number of characters read from *string.  codes[2] is the
   exponent of the least-significant digit read; it is not guaranteed
   to be meaningful for off-scale output. */
static void
LALParseREAL8( LALStatus  *stat,
	       REAL8      *output,
	       const CHAR *string,
	       INT4       codes[3] )
{
  INT2 sign = 1;          /* overall sign of number */
  REAL8 mantissa = 0.0;   /* value normalised to range [0.1,1) */
  INT8 offset = 0;        /* log(normalisation factor) */
  INT4 digits = 0;        /* number of significant digits read */
  REAL8 currentExp = 1.0; /* running count of 10^(digits) */
  INT4 exponent = 0;      /* power of 10 */
  INT2 expSign = 1;       /* sign of exponent */
  UINT4 i = 0;            /* index */

  INITSTATUS(stat);

  /* Check input parameters. */
  ASSERT( output, stat, PULSARCATH_ENUL, PULSARCATH_MSGENUL );
  ASSERT( string, stat, PULSARCATH_ENUL, PULSARCATH_MSGENUL );
  if ( codes ) {
    codes[1] = codes[2] = 0;
    codes[0] = -1;
  }

  /* Skip leading spaces and read sign. */
  while ( isspace( string[i] ) )
    i++;
  if ( string[i] == '-' ) {
    sign = -1.0;
    i++;
  } else if ( string[i] == '+' )
    i++;

  /* Preliminary sanity check: Is there any parseable number? */
  if ( !isdigit( string[i] ) && ( string[i] != '.' ||
				  isdigit( string[i+1] ) ) ) {
    RETURN( stat );
  }

  /* Parse whole part of mantissa. */
  while ( string[i] == '0' )
    i++;
  while ( isdigit( string[i] ) ) {
    offset++;
    digits++;
    currentExp *= 10.0;
    mantissa += (REAL8)( string[i++] - '0' )/currentExp;
  }

  /* Parse fractional part of mantissa. */
  if ( string[i] == '.' ) {
    i++;
    if ( digits == 0 ) {
      while ( string[i] == '0' ) {
	offset--;
	i++;
      }
    }
    while ( isdigit( string[i] ) ) {
      digits++;
      currentExp *= 10.0;
      mantissa += (REAL8)( string[i++] - '0' )/currentExp;
    }
  }

  /* Parse exponent. */
  if ( string[i] == 'e' || string[i] == 'E' ||
       string[i] == 'd' || string[i] == 'D' ) {
    i++;
    if ( string[i] == '-' ) {
      expSign = -1.0;
      i++;
    } else if ( string[i] == '+' )
      i++;
    while ( isdigit( string[i] ) &&
	    exponent + offset*expSign < 214748363 ) {
      exponent *= 10;
      exponent += (INT4)( string[i++] - '0' );
    }
    while ( isdigit( string[i] ) )
      i++;
  }

  /* Number has been read, so set certain output codes. */
  if ( codes ) {
    codes[1] = i;
    codes[2] = expSign*exponent + offset - digits;
  }

  /* Compute output, truncating out-of-range numbers. */
  if ( expSign > 0 ) {
    if ( exponent > 309 ||
	 ( exponent == 309 &&
	   mantissa > 0.17976931348623157 ) ) {
      *output = sign*LAL_REAL8_MAX;
      if ( codes )
	codes[0] = 3;
    } else {
      *output = 10.0*sign*mantissa
	*pow( 10.0, (REAL8)( exponent + offset - 1 ) );
      if ( codes )
	codes[0] = 0;
    }
  } else {
    if ( exponent > 309 ||
	 ( exponent == 309 &&
	   mantissa < 0.22204460492503131 ) ) {
      *output = sign*LAL_REAL8_MIN;
      if ( codes )
	codes[0] = 2;
    } else if ( exponent > 309 ||
		( exponent == 309 &&
		  mantissa < 0.22204460492503131 ) ) {
      *output = 0.1*sign*mantissa
	*pow( 10.0, (REAL8)( -exponent + offset + 1 ) );
      if ( codes )
	codes[0] = 1;
    } else {
      *output = 0.1*sign*mantissa
	*pow( 10.0, (REAL8)( -exponent + offset + 1 ) );
      if ( codes )
	codes[0] = 0;
    }
  }

  /* Everything has been set in the output. */
  RETURN( stat );
}


/* Local routine to parse a token of the form
   "[-hhh...]h:mm:ss[.sss...]"  into a REAL8 number of seconds, with
   additional precision information stored in codes: codes[0] is -1
   for no number parsed, 0 for a normal representable output, 1 for a
   subnormal number, 2 for off-scale low, and 3 for off-scale high.
   codes[1] is the total number of characters read from *string.
   codes[2] is the exponent of the least-significant digit read; it is
   not guaranteed to be meaningful for off-scale output. */
static void
LALParseREAL8HMS( LALStatus  *stat,
		  REAL8      *output,
		  const CHAR *string,
		  INT4       codes[3] )
{
  INT2 sign = 1.0;   /* overall sign of number */
  UINT4 digits = 0;  /* number of digits read after decimal point */
  REAL8 value = 0.0; /* absolute value of number */
  REAL8 norm = 1.0;  /* running count of 10^(digits) */
  UINT4 i = 0;       /* index */

  INITSTATUS(stat);

  /* Check input parameters. */
  ASSERT( output, stat, PULSARCATH_ENUL, PULSARCATH_MSGENUL );
  ASSERT( string, stat, PULSARCATH_ENUL, PULSARCATH_MSGENUL );
  if ( codes ) {
    codes[1] = codes[2] = 0;
    codes[0] = -1;
  }

  /* Skip leading spaces. */
  while ( isspace( string[i] ) )
    i++;
  if ( string[i] == '-' ) {
    sign = -1.0;
    i++;
  } else if ( string[i] == '+' )
    i++;
  if ( !isdigit( string[i] ) ) {
    RETURN( stat );
  }

  /* Read hours. */
  while ( isdigit( string[i] ) ) {
    value *= 10.0;
    value += (REAL8)( string[i++] - '0' );
  }

  /* Check format of remaining part. */
  if ( string[i] != ':' || !isdigit( string[i+1] ) ||
       !isdigit( string[i+2] ) || string[i+3] != ':' ||
       !isdigit( string[i+4] ) || !isdigit( string[i+5] ) ) {
    RETURN( stat );
  }

  /* Read minutes and integral part of seconds. */
  value = 3600.0*value + 600.0*(REAL8)( string[i+1] - '0' ) +
    60.0*(REAL8)( string[i+2] - '0' ) +
    10.0*(REAL8)( string[i+4] - '0' ) + (REAL8)( string[i+5] - '0' );
  i += 6;

  /* Parse fractional part of seconds. */
  if ( string[i] == '.' ) {
    i++;
    while ( isdigit( string[i] ) ) {
      digits++;
      norm *= 10.0;
      value += (REAL8)( string[i++] - '0' )/norm;
    }
  }

  /* Number has been read, so set output codes and return. */
  if ( value >= LAL_REAL8_MAX ) {
    value = LAL_REAL8_MAX;
    if ( codes )
      codes[0] = 3;
  } else if ( value <= LAL_REAL8_MIN ) {
    value = LAL_REAL8_MIN;
    if ( codes )
      codes[0] = 2;
  } else if ( value <= LAL_REAL8_MIN ) {
    value = LAL_REAL8_MIN;
    if ( codes )
      codes[0] = 1;
  } else if ( codes )
    codes[0] = 0;
  *output = sign*value;
  if ( codes ) {
    codes[1] = i;
    codes[2] = -digits;
  }
  RETURN( stat );
}



void
LALReadPulsarCatHead( LALStatus *stat,
		      INT4      indx[PULSARCATINDEX_NUM],
		      TokenList *list )
{
  UINT4 i;                /* an index */
  INITSTATUS(stat);

  /* Check that required input exists. */
  ASSERT( list, stat, PULSARCATH_ENUL, PULSARCATH_MSGENUL );

  /* Set default values. */
  for ( i = 0; i < PULSARCATINDEX_NUM; i++ )
    indx[i] = -1;

  /* Parse token list and assign indx[]. */
  for ( i = 0; i < list->nTokens; i++ ) {
    if ( !strcmp( list->tokens[i], "NAME" ) )
      indx[PULSARCATINDEX_NAME] = i;
    else if ( !strcmp( list->tokens[i], "RAJ" ) ) {
      indx[PULSARCATINDEX_RAJ] = i;
      if ( !strcmp( list->tokens[i+1], "e" ) )
	indx[PULSARCATINDEX_RAJERR] = ++i;
    }
    else if ( !strcmp( list->tokens[i], "DECJ" ) ) {
      indx[PULSARCATINDEX_DECJ] = i;
      if ( !strcmp( list->tokens[i+1], "e" ) )
	indx[PULSARCATINDEX_DECJERR] = ++i;
    }
    else if ( !strcmp( list->tokens[i], "PMRA" ) ) {
      indx[PULSARCATINDEX_PMRA] = i;
      if ( !strcmp( list->tokens[i+1], "e" ) )
	indx[PULSARCATINDEX_PMRAERR] = ++i;
    }
    else if ( !strcmp( list->tokens[i], "PMDEC" ) ) {
      indx[PULSARCATINDEX_PMDEC] = i;
      if ( !strcmp( list->tokens[i+1], "e" ) )
	indx[PULSARCATINDEX_PMDECERR] = ++i;
    }
    else if ( !strcmp( list->tokens[i], "POSEPOCH" ) )
      indx[PULSARCATINDEX_POSEPOCH] = i;
    else if ( !strcmp( list->tokens[i], "F" ) ) {
      indx[PULSARCATINDEX_F] = i;
      if ( !strcmp( list->tokens[i+1], "e" ) )
	indx[PULSARCATINDEX_FERR] = ++i;
    }
    else if ( !strcmp( list->tokens[i], "F1" ) ) {
      indx[PULSARCATINDEX_F1] = i;
      if ( !strcmp( list->tokens[i+1], "e" ) )
	indx[PULSARCATINDEX_F1ERR] = ++i;
    }
    else if ( !strcmp( list->tokens[i], "F2" ) ) {
      indx[PULSARCATINDEX_F2] = i;
      if ( !strcmp( list->tokens[i+1], "e" ) )
	indx[PULSARCATINDEX_F2ERR] = ++i;
    }
    else if ( !strcmp( list->tokens[i], "PEPOCH" ) )
      indx[PULSARCATINDEX_PEPOCH] = i;
    else if ( !strcmp( list->tokens[i], "Dist" ) )
      indx[PULSARCATINDEX_Dist] = i;
  }
  RETURN( stat );
}



void
LALReadPulsarCatLine( LALStatus     *stat,
		      PulsarCatNode *node,
		      TokenList     *list,
		      INT4          indx[PULSARCATINDEX_NUM] )
{
  INT4 i;       /* an index */
  CHAR *endptr; /* pointer to end of a given token */

  INITSTATUS(stat);
  ATTATCHSTATUSPTR( stat );

  /* Check that pointers exist. */
  ASSERT( node, stat, PULSARCATH_ENUL, PULSARCATH_MSGENUL );
  ASSERT( list, stat, PULSARCATH_ENUL, PULSARCATH_MSGENUL );


  /* Parse pulsar name. */
  i = indx[PULSARCATINDEX_NAME];
  if ( i >= 0 && i < (INT4)list->nTokens ) {
    if ( list->tokens[i][0] == 'B' ) {
      memcpy( node->bname, list->tokens[i], 9*sizeof(CHAR) );
      node->bname[9] = '\0';
    } else if ( list->tokens[i][0] == 'J' ) {
      memcpy( node->jname, list->tokens[i], 11*sizeof(CHAR) );
      node->jname[11] = '\0';
    } else {
      LALWarning( stat, "Unrecognized pulsar name" );
    }
  }


  /* Parse right ascension. */
  i = indx[PULSARCATINDEX_RAJ];
  if ( i >= 0 && i < (INT4)list->nTokens && list->tokens[i][0] != '*' ) {
    REAL8 seconds; /* seconds of right ascension */
    INT4 codes[3]; /* return codes from LALParseREAL8HMS() */

    /* Read number of seconds. */
    TRY( LALParseREAL8HMS( stat->statusPtr, &seconds, list->tokens[i],
			   codes ), stat );
    if ( codes[0] == -1 || list->tokens[i][codes[1]] != '\0' ) {
      ABORT( stat, PULSARCATH_EPARSE, PULSARCATH_MSGEPARSE );
    }

    /* Convert to radians. */
    node->pos.system = COORDINATESYSTEM_EQUATORIAL;
    node->pos.longitude = (LAL_PI/43200.0)*seconds;

    /* Compute uncertainty. */
    i = indx[PULSARCATINDEX_RAJERR];
    if ( i >= 0 && i < (INT4)list->nTokens && list->tokens[i][0] != '*' ) {
      UINT2 err;      /* uncertainty */
      TRY( LALStringToU2( stat->statusPtr, &err, list->tokens[i],
			  &endptr ), stat );
      if ( endptr == list->tokens[i] || *endptr != '\0' ) {
	ABORT( stat, PULSARCATH_EPARSE, PULSARCATH_MSGEPARSE );
      }
      node->dpos.system = COORDINATESYSTEM_EQUATORIAL;
      node->dpos.longitude = (LAL_PI/43200.0)*(REAL8)( err )
	*pow( 10.0, (REAL8)( codes[2] ) );
    }
  }


  /* Parse declination. */
  i = indx[PULSARCATINDEX_DECJ];
  if ( i >= 0 && i < (INT4)list->nTokens && list->tokens[i][0] != '*' ) {
    REAL8 seconds; /* arcseconds of declination */
    INT4 codes[3]; /* return codes from LALParseREAL8HMS() */

    /* Read number of seconds. */
    TRY( LALParseREAL8HMS( stat->statusPtr, &seconds, list->tokens[i],
			   codes ), stat );
    if ( codes[0] == -1 || list->tokens[i][codes[1]] != '\0' ) {
      ABORT( stat, PULSARCATH_EPARSE, PULSARCATH_MSGEPARSE );
    }

    /* Convert to radians. */
    node->pos.system = COORDINATESYSTEM_EQUATORIAL;
    node->pos.latitude = (LAL_PI/648000.0)*seconds;

    /* Compute uncertainty. */
    i = indx[PULSARCATINDEX_DECJERR];
    if ( i >= 0 && i < (INT4)list->nTokens && list->tokens[i][0] != '*' ) {
      UINT2 err;      /* uncertainty */
      TRY( LALStringToU2( stat->statusPtr, &err, list->tokens[i],
			  &endptr ), stat );
      if ( endptr == list->tokens[i] || *endptr != '\0' ) {
	ABORT( stat, PULSARCATH_EPARSE, PULSARCATH_MSGEPARSE );
      }
      node->dpos.system = COORDINATESYSTEM_EQUATORIAL;
      node->dpos.latitude = (LAL_PI/648000.0)*(REAL8)( err )
	*pow( 10.0, (REAL8)( codes[2] ) );
    }
  }


  /* Parse right ascension proper motion. */
  i = indx[PULSARCATINDEX_PMRA];
  if ( i >= 0 && i < (INT4)list->nTokens && list->tokens[i][0] != '*' ) {
    REAL8 pmra;    /* proper motion in milliarcseconds per year */
    INT4 codes[3]; /* return codes from LALParseREAL8() */

    /* Parse proper motion. */
    TRY( LALParseREAL8( stat->statusPtr, &pmra, list->tokens[i],
			codes ), stat );
    if ( codes[0] == -1 || list->tokens[i][codes[1]] != '\0' ) {
      ABORT( stat, PULSARCATH_EPARSE, PULSARCATH_MSGEPARSE );
    }

    /* Convert to radians per second. */
    node->pm.system = COORDINATESYSTEM_EQUATORIAL;
    node->pm.longitude = (LAL_PI/6.48e8)/(LAL_YRTROP_SI)*pmra;

    /* Compute uncertainty. */
    i = indx[PULSARCATINDEX_PMRAERR];
    if ( i >= 0 && i < (INT4)list->nTokens && list->tokens[i][0] != '*' ) {
      UINT2 err;      /* uncertainty */
      TRY( LALStringToU2( stat->statusPtr, &err, list->tokens[i], &endptr ),
	   stat );
      if ( endptr == list->tokens[i] || *endptr != '\0' ) {
	ABORT( stat, PULSARCATH_EPARSE, PULSARCATH_MSGEPARSE );
      }
      node->dpm.system = COORDINATESYSTEM_EQUATORIAL;
      node->dpm.longitude = (LAL_PI/6.48e8)/(LAL_YRTROP_SI)
	*(REAL8)( err )*pow( 10.0, (REAL8)( codes[2] ) );
    }
  }


  /* Parse declination proper motion. */
  i = indx[PULSARCATINDEX_PMDEC];
  if ( i >= 0 && i < (INT4)list->nTokens && list->tokens[i][0] != '*' ) {
    REAL8 pmdec;   /* proper motion in milliarcseconds per year */
    INT4 codes[3]; /* return codes from LALParseREAL8() */

    /* Parse proper motion. */
    TRY( LALParseREAL8( stat->statusPtr, &pmdec, list->tokens[i],
			codes ), stat );
    if ( codes[0] == -1 || list->tokens[i][codes[1]] != '\0' ) {
      ABORT( stat, PULSARCATH_EPARSE, PULSARCATH_MSGEPARSE );
    }

    /* Convert to radians per second. */
    node->pm.system = COORDINATESYSTEM_EQUATORIAL;
    node->pm.latitude = (LAL_PI/6.48e8)/(LAL_YRTROP_SI)*pmdec;

    /* Compute uncertainty. */
    i = indx[PULSARCATINDEX_PMDECERR];
    if ( i >= 0 && i < (INT4)list->nTokens && list->tokens[i][0] != '*' ) {
      UINT2 err;      /* uncertainty */
      TRY( LALStringToU2( stat->statusPtr, &err, list->tokens[i], &endptr ),
	   stat );
      if ( endptr == list->tokens[i] || *endptr != '\0' ) {
	ABORT( stat, PULSARCATH_EPARSE, PULSARCATH_MSGEPARSE );
      }
      node->dpm.system = COORDINATESYSTEM_EQUATORIAL;
      node->dpm.latitude = (LAL_PI/6.48e8)/(LAL_YRTROP_SI)
	*(REAL8)( err )*pow( 10.0, (REAL8)( codes[2] ) );
    }
  }


  /* Parse position epoch. */
  i = indx[PULSARCATINDEX_POSEPOCH];
  if ( i >= 0 && i < (INT4)list->nTokens && list->tokens[i][0] != '*' ) {
    REAL8 jday;        /* Julian day */
    INT4 codes[3];     /* return codes from LALParseREAL8() */
    INT8 gpsNan;       /* GPS nanoseconds */
    INT4 leap1, leap2; /* number of leap seconds to date */
    LIGOTimeGPS epoch; /* GPS epoch */

    /* Parse Julian date. */
    TRY( LALParseREAL8( stat->statusPtr, &jday, list->tokens[i],
			codes ), stat );
    if ( codes[0] == -1 || list->tokens[i][codes[1]] != '\0' ) {
      ABORT( stat, PULSARCATH_EPARSE, PULSARCATH_MSGEPARSE );
    }
    if ( jday < 1.0e6 )
      jday += 2.4e6;

    /* Convert Julian days to GPS nanoseconds. */
    gpsNan = (INT8)( ( jday - 2444244.5 )*(8.64e13L) );
    XLALINT8NSToGPS( &epoch, gpsNan );
    leap2 = 0;
    do {
      leap1 = leap2;
      leap2 = XLALGPSLeapSeconds ( epoch.gpsSeconds );
      epoch.gpsSeconds += leap2 - leap1;
    } while ( leap2 != leap1 );
    node->posepoch = epoch;
  }


  /* Parse frequency and its derivatives. */
  i = indx[PULSARCATINDEX_F];
  if ( i >= 0 && i < (INT4)list->nTokens && list->tokens[i][0] != '*' ) {
    UINT4 length = 1;  /* number of frequency terms */
    REAL8 f[3], df[3]; /* frequency terms and uncertainties */
    INT4 codes[3];     /* return codes from LALParseREAL8() */
    memset( f, 0, 3*sizeof(REAL8) );
    memset( df, 0, 3*sizeof(REAL8) );

    /* First read higher derivatives. */
    i = indx[PULSARCATINDEX_F1];
    if ( i >= 0 && i < (INT4)list->nTokens && list->tokens[i][0] != '*' ) {
      length = 2;
      i = indx[PULSARCATINDEX_F2];
      if ( i >= 0 && i < (INT4)list->nTokens && list->tokens[i][0] != '*' ) {
	length = 3;

	/* Parse F2. */
	TRY( LALParseREAL8( stat->statusPtr, f+2, list->tokens[i],
			    codes ), stat );
	if ( codes[0] == -1 || list->tokens[i][codes[1]] != '\0' ) {
	  ABORT( stat, PULSARCATH_EPARSE, PULSARCATH_MSGEPARSE );
	}

	/* Compute uncertainty. */
	i = indx[PULSARCATINDEX_F2ERR];
	if ( i >= 0 && i < (INT4)list->nTokens && list->tokens[i][0] != '*' ) {
	  UINT2 err;      /* uncertainty */
	  TRY( LALStringToU2( stat->statusPtr, &err, list->tokens[i],
			      &endptr ), stat );
	  if ( endptr == list->tokens[i] || *endptr != '\0' ) {
	    ABORT( stat, PULSARCATH_EPARSE, PULSARCATH_MSGEPARSE );
	  }
	  df[2] = (REAL8)( err )*pow( 10.0, (REAL8)( codes[2] ) );
	}
      }

      /* Parse F1. */
      i = indx[PULSARCATINDEX_F1];
      TRY( LALParseREAL8( stat->statusPtr, f+1, list->tokens[i],
			  codes ), stat );
      if ( codes[0] == -1 || list->tokens[i][codes[1]] != '\0' ) {
	ABORT( stat, PULSARCATH_EPARSE, PULSARCATH_MSGEPARSE );
      }

      /* Compute uncertainty. */
      i = indx[PULSARCATINDEX_F1ERR];
      if ( i >= 0 && i < (INT4)list->nTokens && list->tokens[i][0] != '*' ) {
	UINT2 err;      /* uncertainty */
	TRY( LALStringToU2( stat->statusPtr, &err, list->tokens[i],
			    &endptr ), stat );
	if ( endptr == list->tokens[i] || *endptr != '\0' ) {
	  ABORT( stat, PULSARCATH_EPARSE, PULSARCATH_MSGEPARSE );
	}
	df[1] = (REAL8)( err )*pow( 10.0, (REAL8)( codes[2] ) );
      }
    }

    /* Parse F. */
    i = indx[PULSARCATINDEX_F];
    TRY( LALParseREAL8( stat->statusPtr, f, list->tokens[i],
			codes ), stat );
    if ( codes[0] == -1 || list->tokens[i][codes[1]] != '\0' ) {
      ABORT( stat, PULSARCATH_EPARSE, PULSARCATH_MSGEPARSE );
    }

    /* Compute uncertainty. */
    i = indx[PULSARCATINDEX_FERR];
    if ( i >= 0 && i < (INT4)list->nTokens && list->tokens[i][0] != '*' ) {
      UINT2 err;      /* uncertainty */
      TRY( LALStringToU2( stat->statusPtr, &err, list->tokens[i],
			  &endptr ), stat );
      if ( endptr == list->tokens[i] || *endptr != '\0' ) {
	ABORT( stat, PULSARCATH_EPARSE, PULSARCATH_MSGEPARSE );
      }
      df[0] = (REAL8)( err )*pow( 10.0, (REAL8)( codes[2] ) );
    }

    /* Allocate space and fill arrays. */
    if ( node->f ) {
      TRY( LALDDestroyVector( stat->statusPtr, &(node->f) ), stat );
    }
    if ( node->df ) {
      TRY( LALDDestroyVector( stat->statusPtr, &(node->df) ), stat );
    }
    TRY( LALDCreateVector( stat->statusPtr, &(node->f), length ),
	 stat );
    TRY( LALDCreateVector( stat->statusPtr, &(node->df), length ),
	 stat );
    memcpy( node->f->data, f, length*sizeof(REAL8) );
    memcpy( node->df->data, df, length*sizeof(REAL8) );
  }


  /* Parse frequency epoch. */
  i = indx[PULSARCATINDEX_PEPOCH];
  if ( i >= 0 && i < (INT4)list->nTokens && list->tokens[i][0] != '*' ) {
    REAL8 jday;        /* Julian day */
    INT4 codes[3];     /* return codes from LALParseREAL8() */
    INT8 gpsNan;       /* GPS nanoseconds */
    INT4 leap1, leap2; /* number of leap seconds to date */
    LIGOTimeGPS epoch; /* GPS epoch */

    /* Parse Julian date. */
    TRY( LALParseREAL8( stat->statusPtr, &jday, list->tokens[i],
			codes ), stat );
    if ( codes[0] == -1 || list->tokens[i][codes[1]] != '\0' ) {
      ABORT( stat, PULSARCATH_EPARSE, PULSARCATH_MSGEPARSE );
    }
    if ( jday < 1.0e6 )
      jday += 2.4e6;

    /* Convert Julian days to GPS nanoseconds. */
    gpsNan = (INT8)( ( jday - 2444244.5 )*(8.64e13L) );
    XLALINT8NSToGPS( &epoch, gpsNan );
    leap2 = 0;
    do {
      leap1 = leap2;
      leap2 = XLALGPSLeapSeconds ( epoch.gpsSeconds );
      epoch.gpsSeconds += leap2 - leap1;
    } while ( leap2 != leap1 );
    node->fepoch = epoch;
  }


  /* Parse distance. */
  i = indx[PULSARCATINDEX_Dist];
  if ( i >= 0 && i < (INT4)list->nTokens && list->tokens[i][0] != '*' ) {
    REAL8 dist;    /* distance in kpc */
    INT4 codes[3]; /* return codes from LALParseREAL8() */

    TRY( LALParseREAL8( stat->statusPtr, &dist, list->tokens[i],
			codes ), stat );
    if ( codes[0] == -1 || list->tokens[i][codes[1]] != '\0' ) {
      ABORT( stat, PULSARCATH_EPARSE, PULSARCATH_MSGEPARSE );
    }
    node->dist = 1000.0*LAL_PC_SI*dist;
  }


  /* Done. */
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}
