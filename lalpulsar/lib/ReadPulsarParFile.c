/*
*  Copyright (C) 2013 Jolien Creighton, Matt Pitkin
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
*  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
*  MA  02110-1301  USA
*/

/**
 * \author Matt Pitkin
 * \date 2013
 * \file
 * \ingroup lalpulsar_general
 * \brief Functions to read TEMPO pulsar parameter files
 *
   Functions for reading pulsar parameters from TEMPO .par files.

   # Prototypes



   # Description

   Radio astronomers fit pulsar parameters using TEMPO(2) which will output
   the parameters in a <tt>.par</tt> file. The values allowed in this file can be
   found in the TEMPO documentation. Two function are available to extract these
   parameters from the <tt>.par</tt> files:
   <ul>
   <li>\c XLALReadTEMPOParFileOrig - this will read parameters into a \c BinaryPulsarParams structure and set
   any unused parameters to zero or \c NULL. To use this you must know the correct parameter name within
   the structure. This is deprecated in favour of \c XLALReadTEMPOParFile and is no longer maintained</li>
   <li>\c XLALReadTEMPOParFile - this reads the parameters into a linked list \c PulsarParameters structure, from which the
   parameters can be accessed using the appropriate access function. These use a hash table to quick look-up.
   The parameters are assigned names, which are used as the hash table keys, which are fully uppercase
   versions of the TEMPO parameter names.
   </ul>
   All parameters read in are converted into SI units.

   Functions are is also included which converts a string containing the right ascension or
   declination in the format <tt>ddd/hh:mm:ss.s</tt> or <tt>ddd/hhmmss.s</tt>
   (as is given in the <tt>.par</tt> file) into a \c REAL8 value in
   radians.

   # Notes

*/

#include <config.h>
#include <string.h>
#include <math.h>
#include <locale.h>
#include <ctype.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include <lal/LALConstants.h>
#include <lal/LALStdlib.h>
#include <lal/LALString.h>
#include <lal/ComputeFstat.h>
#include <lal/TranslateAngles.h>
#include <lal/TranslateMJD.h>
#include <lal/LALHashFunc.h>
#include <lal/ReadPulsarParFile.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif


#define DAYSTOSECS 86400.0 /* number of seconds in an SI day */

size_t PulsarTypeSize[5] = {
  sizeof( UINT4 ),
  sizeof( REAL8 ),
  sizeof( REAL8Vector * ),
  sizeof( CHAR ) *PULSAR_PARNAME_MAX,
  sizeof( void * )
};


/** Array for conversion from lowercase to uppercase */
static const CHAR a2A[256] = {
  ['a'] = 'A', ['b'] = 'B', ['c'] = 'C', ['d'] = 'D', ['e'] = 'E', ['f'] = 'F', ['g'] = 'G', ['h'] = 'H',
  ['i'] = 'I', ['j'] = 'J', ['k'] = 'K', ['l'] = 'L', ['m'] = 'M', ['n'] = 'N', ['o'] = 'O', ['p'] = 'P',
  ['q'] = 'Q', ['r'] = 'R', ['s'] = 'S', ['t'] = 'T', ['u'] = 'U', ['v'] = 'V', ['w'] = 'W', ['x'] = 'X',
  ['y'] = 'Y', ['z'] = 'Z'
};


/** \brief Convert string to uppercase */
static void strtoupper( CHAR *s )
{
  /* convert everything to uppercase */
  CHAR c;
  for ( ; *s; ++s ) {
    if ( ( c = a2A[( int )( *s )] ) ) {
      *s = c;
    }
  }
}

/* hash functions copied from LALInference.c */
typedef struct taghash_elem {
  const char *name;
  PulsarParam *itemPtr;
} hash_elem;

static void *new_elem( const char *name, PulsarParam *itemPtr )
{
  hash_elem e = { .name = name, .itemPtr = itemPtr };
  return memcpy( XLALMalloc( sizeof( e ) ), &e, sizeof( e ) );
};

static void del_elem( void *elem )
{
  XLALFree( elem );
}


/* Compute a hash value based on element */
static UINT8 PulsarHash( const void *elem );
static UINT8 PulsarHash( const void *elem )
{
  if ( !elem ) {
    XLAL_ERROR( XLAL_EINVAL );
  }
  size_t len = strnlen( ( ( const hash_elem * )elem )->name, PULSAR_PARNAME_MAX );
  return ( XLALCityHash64( ( ( const hash_elem * )elem )->name, len ) );
}


static int PulsarHashElemCmp( const void *elem1, const void *elem2 );
static int PulsarHashElemCmp( const void *elem1, const void *elem2 )
{
  if ( !elem1 || !elem2 ) {
    XLAL_ERROR( XLAL_EINVAL );
  }
  return ( strncmp( ( ( const hash_elem * )elem1 )->name, ( ( const hash_elem * )elem2 )->name, PULSAR_PARNAME_MAX ) );
}


/** \brief Get a pointer to a parameter of a given name from a \c PulsarParameters structure
 *
 * Note this function can only be used internally.
 */
static PulsarParam *PulsarGetParamItemSlow( const PulsarParameters *pars, const CHAR *name )
{
  /* (this function is only to be used internally) */
  /* Returns pointer to item for given item name.  */
  if ( pars == NULL ) {
    return NULL;
  }
  if ( pars->nparams == 0 ) {
    return NULL;
  }

  PulsarParam *this = pars->head;

  /* loop through all values checking the name */
  while ( this != NULL ) {
    if ( !strcmp( this->name, name ) ) {
      break;
    } else {
      this = this->next;
    }
  }

  return ( this );
}


/** \brief Get a pointer to a parameter of a given name from a \c PulsarParameters structure
 *
 * This function will return a pointer to the parameter. It will initially try and use the parameter's hash name,
 * otherwise it will use \c PulsarGetParamItemSlow to loop through all parameters.
 *
 * Note this function can only be used internally.
 */
static PulsarParam *PulsarGetParamItem( const PulsarParameters *pars, const CHAR *name )
{
  /* (this function is only to be used internally) */
  /* Returns pointer to item for given item name.  */
  hash_elem tmp;
  const hash_elem *match = NULL;

  CHAR upperName[PULSAR_PARNAME_MAX];
  XLALStringCopy( upperName, name, PULSAR_PARNAME_MAX );
  strtoupper( upperName );

  if ( pars == NULL ) {
    return NULL;
  }
  if ( pars->nparams == 0 ) {
    return NULL;
  }
  if ( !pars->hash_table ) {
    return PulsarGetParamItemSlow( pars, upperName );
  }

  tmp.name = upperName;
  XLALHashTblFind( pars->hash_table, &tmp, ( const void ** )&match );
  if ( !match ) {
    return NULL;
  }

  return match->itemPtr;
}


void *PulsarGetParam( const PulsarParameters *pars, const CHAR *name )
{
  /* Return the value of variable name from the pars structure */
  PulsarParam *item = PulsarGetParamItem( pars, name );
  if ( !item ) {
    XLAL_ERROR_NULL( XLAL_EFAILED, "Entry \"%s\" not found.", name );
  }

  return ( item->value );
}


PulsarParamType PulsarGetParamType( const PulsarParameters *pars, const char *name )
{
  return PulsarGetParamItem( pars, name )->type;
}


REAL8 PulsarGetREAL8Param( const PulsarParameters *pars, const CHAR *name )
{
  /* check type is a REAL8 */
  if ( PulsarGetParamType( pars, name ) == PULSARTYPE_REAL8_t ) {
    return *( REAL8 * )PulsarGetParam( pars, name );
  } else {
    XLAL_ERROR_REAL8( XLAL_EINVAL, "Used wrong type for required parameter" );
  }
}


REAL8 PulsarGetREAL8ParamOrZero( const PulsarParameters *pars, const CHAR *name )
{
  /* if parameter is there return it otherwise return zero */
  return PulsarCheckParam( pars, name ) ? PulsarGetREAL8Param( pars, name ) : 0.;
}


UINT4 PulsarGetUINT4Param( const PulsarParameters *pars, const CHAR *name )
{
  /* check type is a UINT4 */
  if ( PulsarGetParamType( pars, name ) == PULSARTYPE_UINT4_t ) {
    return *( UINT4 * )PulsarGetParam( pars, name );
  } else {
    XLAL_ERROR( XLAL_EINVAL, "Used wrong type for required parameter" );
  }
}


UINT4 PulsarGetUINT4ParamOrZero( const PulsarParameters *pars, const CHAR *name )
{
  /* if parameter is there return it otherwise return zero */
  return PulsarCheckParam( pars, name ) ? PulsarGetUINT4Param( pars, name ) : 0;
}


const CHAR *PulsarGetStringParam( const PulsarParameters *pars, const CHAR *name )
{
  /* check type is a string */
  if ( PulsarGetParamType( pars, name ) == PULSARTYPE_string_t ) {
    CHAR *rvalue = *( CHAR ** )PulsarGetParam( pars, name );
    return rvalue;
  } else {
    XLAL_ERROR_NULL( XLAL_EINVAL, "Used wrong type for required parameter" );
  }
}


const REAL8Vector *PulsarGetREAL8VectorParam( const PulsarParameters *pars, const CHAR *name )
{
  /* check type is a REAL8Vector */
  if ( PulsarGetParamType( pars, name ) == PULSARTYPE_REAL8Vector_t ) {
    return *( REAL8Vector ** )PulsarGetParam( pars, name );
  } else {
    XLAL_ERROR_NULL( XLAL_EINVAL, "Used wrong type for required parameter" );
  }
}


REAL8 PulsarGetREAL8VectorParamIndividual( const PulsarParameters *pars, const CHAR *name )
{
  /* split the input name into the parameter name and an index e.g. FB0 into FB and 0, or GLEP_1 */
  CHAR *namecopy = NULL, *token;
  namecopy = XLALStringDuplicate( name );
  REAL8 val = 0.;

  if ( strchr( name, '_' ) ) { /* check for underscore delimiting numbers */
    token = strtok( namecopy, "_" );
  } else { /* name delimited by just numbers */
    token = strtok( namecopy, "0123456789" );
  }

  const REAL8Vector *vpars = PulsarGetREAL8VectorParam( pars, token );

  /* get the index value of the parameter name */
  INT4 idx = -1;
  INT4 loc = strlen( token );
  if ( name[loc] == '_' ) {
    loc++;
  }
  if ( sscanf( name + loc, "%d", &idx ) != 1 ) {
    XLAL_ERROR_REAL8( XLAL_EINVAL, "Input parameter (%s) problem", name );
  }

  /* for glitch parameters GL, or WAVE parameters the index starts at 1, so shift it */
  if ( !strncmp( name, "GL", 2 ) || !strncmp( name, "WAVE", 4 ) ) {
    idx--;
  }
  if ( idx < 0 || ( UINT4 )idx > vpars->length - 1 ) {
    XLAL_ERROR_REAL8( XLAL_EINVAL, "Input parameter index %d is wrong", idx );
  }

  val = vpars->data[idx];

  XLALFree( namecopy );

  return val;
}


void *PulsarGetParamErr( const PulsarParameters *pars, const CHAR *name )
{
  /* Return the error value of variable name from the pars structure */
  PulsarParam *item = PulsarGetParamItem( pars, name );
  if ( !item ) {
    XLAL_ERROR_NULL( XLAL_EFAILED, "Entry \"%s\" not found.", name );
  }

  return ( item->err );
}


const UINT4 *PulsarGetParamFitFlag( const PulsarParameters *pars, const CHAR *name )
{
  /* return the fit flag value */
  PulsarParam *item = PulsarGetParamItem( pars, name );
  if ( !item ) {
    XLAL_ERROR_NULL( XLAL_EFAILED, "Entry \"%s\" not found.", name );
  }

  if ( item->fitFlag == NULL ) {
    return NULL;
  } else {
    if ( item->fitFlag->data == NULL ) {
      return NULL;
    } else {
      return item->fitFlag->data;
    }
  }
}


const UINT4Vector *PulsarGetParamFitFlagAsVector( const PulsarParameters *pars, const CHAR *name )
{
  /* return the fit flag value(s) */
  PulsarParam *item = PulsarGetParamItem( pars, name );
  if ( !item ) {
    XLAL_ERROR_NULL( XLAL_EFAILED, "Entry \"%s\" not found.", name );
  }

  if ( item->fitFlag == NULL ) {
    return NULL;
  } else {
    return item->fitFlag;
  }
}


REAL8 PulsarGetREAL8ParamErr( const PulsarParameters *pars, const CHAR *name )
{
  /* Return the error value of variable name from the pars structure */
  /* check type is a REAL8 */
  if ( PulsarGetParamType( pars, name ) == PULSARTYPE_REAL8_t ) {
    void *val = PulsarGetParamErr( pars, name );
    if ( val == NULL ) {
      return 0.;  /* if no error is present return 0 */
    } else {
      return *( REAL8 * )val;
    }
  } else {
    XLAL_ERROR_REAL8( XLAL_EINVAL, "Used wrong type for required parameter" );
  }
}


const REAL8Vector *PulsarGetREAL8VectorParamErr( const PulsarParameters *pars, const CHAR *name )
{
  /* Return the error value of variable name from the pars structure */
  /* check type is a REAL8 */
  if ( PulsarGetParamType( pars, name ) == PULSARTYPE_REAL8Vector_t ) {
    void *val = PulsarGetParamErr( pars, name );
    if ( val == NULL ) {
      return NULL;  /* if no error is present return NULL */
    } else {
      return *( REAL8Vector ** )val;
    }
  } else {
    XLAL_ERROR_NULL( XLAL_EINVAL, "Used wrong type for required parameter" );
  }
}


REAL8 PulsarGetREAL8VectorParamErrIndividual( const PulsarParameters *pars, const CHAR *name )
{
  /* split the input name into the parameter name and an index e.g. FB0 into FB and 0*/
  CHAR *namecopy = NULL, *token;
  namecopy = XLALStringDuplicate( name );
  REAL8 val = 0.;

  if ( strchr( name, '_' ) ) { /* check for underscore delimiting numbers */
    token = strtok( namecopy, "_" );
  } else { /* name delimited by just numbers */
    token = strtok( namecopy, "0123456789" );
  }

  const REAL8Vector *vpars = PulsarGetREAL8VectorParamErr( pars, token );

  /* get the index value of the parameter name */
  INT4 idx = -1;
  INT4 loc = strlen( token );
  if ( name[loc] == '_' ) {
    loc++;
  }
  if ( sscanf( name + loc, "%d", &idx ) != 1 ) {
    XLAL_ERROR_REAL8( XLAL_EINVAL, "Input parameter (%s) problem", name );
  }

  /* for glitch parameters GL, or WAVE parameters the index starts at 1, so shift it */
  if ( !strncmp( name, "GL", 2 ) || !strncmp( name, "WAVE", 4 ) ) {
    idx--;
  }
  if ( idx < 0 || ( UINT4 )idx > vpars->length - 1 ) {
    XLAL_ERROR_REAL8( XLAL_EINVAL, "Input parameter index %d is wrong", idx );
  }

  val = vpars->data[idx];

  XLALFree( namecopy );

  return val;
}


void PulsarAddParam( PulsarParameters *pars, const CHAR *name, void *value, PulsarParamType type )
{
  /* Add the variable name with type type and value value to pars */
  /* If variable already exists, it will over-write the current value if type compatible*/
  PulsarParam *old = NULL;

  /* create the hash table if it does not exist */
  if ( !pars->hash_table ) {
    pars->hash_table = XLALHashTblCreate( del_elem, PulsarHash, PulsarHashElemCmp );
  }

  /* Check input value is accessible */
  if ( !value ) {
    XLAL_ERROR_VOID( XLAL_EFAULT, "Unable to access value through null pointer; trying to add \"%s\".", name );
  }

  /* Check the name doesn't already exist */
  if ( PulsarCheckParam( pars, name ) ) {
    old = PulsarGetParamItem( pars, name );

    /* If the type is incompatible, it should be removed */
    if ( old->type != type || old->type == PULSARTYPE_REAL8Vector_t ) {
      PulsarRemoveParam( pars, name );
    } else {
      PulsarSetParam( pars, name, value );
      return;
    }
  }

  /* If we get this far it is safe to create a new node for this variable */
  PulsarParam *new = XLALMalloc( sizeof( PulsarParam ) );
  memset( new, 0, sizeof( PulsarParam ) );
  if ( new ) {
    new->value = ( void * )XLALMalloc( PulsarTypeSize[type] );
    new->fitFlag = NULL;

    /* for REAL8 and REAL8Vector parameters initialise errors and fit flags to zero  */
    if ( type == PULSARTYPE_REAL8Vector_t ) {
      REAL8Vector *tmpvector = *( REAL8Vector ** )value;
      UINT4 dlength = tmpvector->length;

      new->fitFlag = XLALCreateUINT4Vector( dlength );
      memset( new->fitFlag->data, 0, sizeof( UINT4 )*dlength ); // initialise to zero

      REAL8Vector *err = XLALCreateREAL8Vector( dlength );
      memset( err->data, 0, sizeof( REAL8 )*dlength ); // initialise to zero

      new->err = ( void * )XLALMalloc( PulsarTypeSize[type] );
      memcpy( new->err, ( void * )&err, PulsarTypeSize[type] );
    } else if ( type == PULSARTYPE_REAL8_t ) {
      new->fitFlag = XLALCreateUINT4Vector( 1 );
      memset( new->fitFlag->data, 0, sizeof( UINT4 ) ); // initialise to zero

      new->err = ( void * )XLALMalloc( sizeof( UINT4Vector * ) );
      REAL8 err = 0.;
      memcpy( new->err, ( void * )&err, PulsarTypeSize[type] );
    }
  }

  if ( new == NULL || new->value == NULL ) {
    XLAL_ERROR_VOID( XLAL_ENOMEM, "Unable to allocate memory for list item." );
  }

  CHAR upperName[PULSAR_PARNAME_MAX];
  XLALStringCopy( upperName, name, PULSAR_PARNAME_MAX );
  strtoupper( upperName );

  XLALStringCopy( new->name, upperName, PULSAR_PARNAME_MAX );
  new->type = type;
  memcpy( new->value, value, PulsarTypeSize[type] );
  new->next = pars->head;
  pars->head = new;
  hash_elem *elem = new_elem( new->name, new );
  XLALHashTblAdd( pars->hash_table, ( void * )elem );
  pars->nparams++;
}


void PulsarAddREAL8Param( PulsarParameters *pars, const CHAR *name, REAL8 value )
/* Typed version of PulsarAddParam for REAL8 values.*/
{
  PulsarAddParam( pars, name, ( void * )&value, PULSARTYPE_REAL8_t );
}


void PulsarAddUINT4Param( PulsarParameters *pars, const CHAR *name, UINT4 value )
/* Typed version of PulsarAddParam for UINT4 values.*/
{
  PulsarAddParam( pars, name, ( void * )&value, PULSARTYPE_UINT4_t );
}


void PulsarAddREAL8VectorParam( PulsarParameters *pars, const CHAR *name, const REAL8Vector *value )
/* Typed version of PulsarAddParam for REAL8Vector values.*/
{
  REAL8Vector *fval = NULL;
  fval = XLALCreateREAL8Vector( value->length );
  memcpy( fval->data, value->data, sizeof( REAL8 )*value->length );
  PulsarAddParam( pars, name, ( void * )&fval, PULSARTYPE_REAL8Vector_t );
}


void PulsarAddStringParam( PulsarParameters *pars, const CHAR *name, const CHAR *value )
/* Typed version of PulsarAddParam for string values.*/
{
  CHAR *sval = NULL;
  sval = XLALMalloc( PulsarTypeSize[PULSARTYPE_string_t] );
  XLALStringCopy( sval, value, PulsarTypeSize[PULSARTYPE_string_t] );
  PulsarAddParam( pars, name, ( void * )&sval, PULSARTYPE_string_t );
}


/* Check for existence of name */
int PulsarCheckParam( const PulsarParameters *pars, const CHAR *name )
{
  /* convert name to uppercase */
  CHAR upperName[PULSAR_PARNAME_MAX];
  XLALStringCopy( upperName, name, PULSAR_PARNAME_MAX );
  strtoupper( upperName );

  if ( PulsarGetParamItem( pars, upperName ) ) {
    return 1;
  } else {
    return 0;
  }
}


void PulsarClearParams( PulsarParameters *pars )
{
  /* Free all variables inside the linked list, leaving only the head struct */
  PulsarParam *this, *next = NULL;

  if ( !pars ) {
    return;
  }

  this = pars->head;

  if ( this ) {
    next = this->next;
  }

  while ( this ) {
    if ( this->type == PULSARTYPE_REAL8Vector_t ) {
      if ( this->value ) {
        XLALDestroyREAL8Vector( *( REAL8Vector ** )this->value );
      }
      if ( this->err ) {
        XLALDestroyREAL8Vector( *( REAL8Vector ** )this->err );
      }
    } else if ( this->type == PULSARTYPE_string_t ) {
      if ( this->value ) {
        XLALFree( *( CHAR ** )this->value );
      }
    }
    if ( this->fitFlag ) {
      XLALDestroyUINT4Vector( this->fitFlag );
    }
    XLALFree( this->value );
    XLALFree( this->err );
    XLALFree( this );
    this = next;
    if ( this ) {
      next = this->next;
    }
  }
  pars->head = NULL;
  pars->nparams = 0;

  if ( pars->hash_table ) {
    XLALHashTblDestroy( pars->hash_table );
  }
  pars->hash_table = NULL;
}


/* remove a given parameter */
void PulsarRemoveParam( PulsarParameters *pars, const CHAR *name )
{
  PulsarParam *this;

  if ( !pars ) {
    XLAL_ERROR_VOID( XLAL_EFAULT );
  }

  CHAR upperName[PULSAR_PARNAME_MAX];
  XLALStringCopy( upperName, name, PULSAR_PARNAME_MAX );
  strtoupper( upperName );

  this = pars->head;
  PulsarParam *parent = NULL;

  while ( this ) {
    if ( !strcmp( this->name, upperName ) ) {
      break;
    } else {
      parent = this;
      this = this->next;
    }
  }
  if ( !this ) {
    XLAL_PRINT_WARNING( "Entry \"%s\" not found.", upperName );
    return;
  }

  if ( !parent ) {
    pars->head = this->next;
  } else {
    parent->next = this->next;
  }
  /* Remove from hash table */
  hash_elem elem;
  elem.name = this->name;
  XLALHashTblRemove( pars->hash_table, ( void * )&elem );
  if ( this->type == PULSARTYPE_REAL8Vector_t ) {
    if ( this->value ) {
      XLALDestroyREAL8Vector( *( REAL8Vector ** )this->value );
    }
    if ( this->err ) {
      XLALDestroyREAL8Vector( *( REAL8Vector ** )this->err );
    }
  } else if ( this->type == PULSARTYPE_string_t ) {
    if ( this->value ) {
      XLALFree( *( CHAR ** )this->value );
    }
  }
  if ( this->fitFlag ) {
    XLALDestroyUINT4Vector( this->fitFlag );
  }
  XLALFree( this->value );
  XLALFree( this->err );
  this->value = NULL;
  this->err = NULL;
  this->fitFlag = NULL;
  XLALFree( this );

  this = NULL;
  pars->nparams--;
}


/* Set the value of parameter name in the pars structure to value */
void PulsarSetParam( PulsarParameters *pars, const CHAR *name, const void *value )
{
  PulsarParam *item;

  /* convert name to uppercase */
  CHAR upperName[PULSAR_PARNAME_MAX];
  XLALStringCopy( upperName, name, PULSAR_PARNAME_MAX );
  strtoupper( upperName );

  item = PulsarGetParamItem( pars, upperName );
  if ( !item ) {
    XLAL_ERROR_VOID( XLAL_EINVAL, "Entry \"%s\" not found.", upperName );
  }
  memcpy( item->value, value, PulsarTypeSize[item->type] );
}


/* Set the value of parameter name error in the pars structure to value */
void PulsarSetParamErr( PulsarParameters *pars, const CHAR *name, void *value, const UINT4 *fitFlag, UINT4 len )
{
  PulsarParam *item;

  /* convert name to uppercase */
  CHAR upperName[PULSAR_PARNAME_MAX];
  XLALStringCopy( upperName, name, PULSAR_PARNAME_MAX );
  strtoupper( upperName );

  item = PulsarGetParamItem( pars, upperName );
  if ( !item ) {
    XLAL_ERROR_VOID( XLAL_EINVAL, "Entry \"%s\" not found.", name );
  }

  if ( item->err != NULL ) {
    // Remove previous error values if present
    if ( item->type == PULSARTYPE_REAL8Vector_t ) {
      XLALDestroyREAL8Vector( *( REAL8Vector ** )item->err );
    }
    XLALFree( item->err );
  }

  if ( ( item->err = ( void * )XLALMalloc( PulsarTypeSize[PulsarGetParamType( pars, name )] ) ) == NULL ) {
    XLAL_ERROR_VOID( XLAL_ENOMEM, "Unable to allocate memory for list item." );
  }

  // Remove previous fit flag if present
  if ( item->fitFlag != NULL ) {
    XLALDestroyUINT4Vector( item->fitFlag );

    // set the fit flag values
    item->fitFlag = NULL;
    item->fitFlag = XLALCreateUINT4Vector( len );
    memcpy( item->fitFlag->data, fitFlag, len * sizeof( UINT4 ) );
  }

  // set error values
  memcpy( item->err, value, PulsarTypeSize[item->type] );
}


/* Set the value of parameter name error in the pars structure to value */
void PulsarSetREAL8ParamErr( PulsarParameters *pars, const CHAR *name, REAL8 value, UINT4 fitFlag )
{
  PulsarSetParamErr( pars, name, ( void * )&value, &fitFlag, 1 );
}


/* Set the value of parameter name error in the pars structure to value */
void PulsarSetREAL8VectorParamErr( PulsarParameters *pars, const CHAR *name, const REAL8Vector *value, const UINT4 *fitFlag )
{
  REAL8Vector *evals = NULL;
  evals = XLALCreateREAL8Vector( value->length );
  memcpy( evals->data, value->data, sizeof( REAL8 )*value->length );

  PulsarSetParamErr( pars, name, ( void * )&evals, fitFlag, value->length );
}


/* free memory */
void PulsarFreeParams( PulsarParameters *pars )
{
  if ( !pars ) {
    return;
  }

  PulsarClearParams( pars );
  if ( pars->hash_table ) {
    XLALHashTblDestroy( pars->hash_table );
  }
  XLALFree( pars );
}


/* create a copy of the parameters */
void PulsarCopyParams( PulsarParameters *origin, PulsarParameters *target )
{
  /* Check that the source and origin differ */
  if ( origin == target ) {
    return;
  }

  if ( !origin ) {
    XLAL_ERROR_VOID( XLAL_EFAULT, "Unable to access origin pointer." );
  }

  /* Make sure the structure is initialised */
  if ( !target ) {
    XLAL_ERROR_VOID( XLAL_EFAULT, "Unable to copy to uninitialised PulsarParameters structure." );
  }

  PulsarParam *this, *next = NULL;

  this = origin->head;

  if ( this ) {
    next = this->next;
  }

  /* first clear target variables */
  PulsarClearParams( target );

  while ( this ) {
    if ( this->type == PULSARTYPE_REAL8Vector_t ) {
      const REAL8Vector *old = *( REAL8Vector ** )this->value;
      REAL8Vector *new = XLALCreateREAL8Vector( old->length );
      if ( new ) {
        memcpy( new->data, old->data, new->length * sizeof( REAL8 ) );
      } else {
        XLAL_ERROR_VOID( XLAL_ENOMEM, "Unable to copy vector!\n" );
      }

      PulsarAddREAL8VectorParam( target, this->name, ( const REAL8Vector * )new );

      if ( this->err ) {
        const REAL8Vector *olderr = *( REAL8Vector ** )this->err;
        REAL8Vector *newerr = XLALCreateREAL8Vector( olderr->length );
        if ( newerr ) {
          memcpy( newerr->data, olderr->data, newerr->length * sizeof( REAL8 ) );
        } else {
          XLAL_ERROR_VOID( XLAL_ENOMEM, "Unable to copy vector!\n" );
        }

        const UINT4 *oldfit = PulsarGetParamFitFlag( origin, this->name );
        UINT4 *newfit = ( UINT4 * )XLALCalloc( newerr->length, sizeof( UINT4 ) );

        if ( newfit ) {
          memcpy( newfit, oldfit, newerr->length * sizeof( UINT4 ) );
        } else {
          XLAL_ERROR_VOID( XLAL_ENOMEM, "Unable to copy vector!\n" );
        }

        PulsarSetREAL8VectorParamErr( target, this->name, ( const REAL8Vector * )newerr, ( const UINT4 * )newfit );
      }
    } else if ( this->type == PULSARTYPE_string_t ) {
      CHAR *parstr = XLALStringDuplicate( *( CHAR ** )this->value );
      PulsarAddStringParam( target, this->name, ( const CHAR * )parstr );
    } else {
      PulsarAddParam( target, this->name, this->value, this->type );

      if ( ( this->type == PULSARTYPE_REAL8_t ) && this->err ) {
        REAL8 err = *( REAL8 * )this->err;

        const UINT4 *oldfit = PulsarGetParamFitFlag( origin, this->name );
        UINT4 newfit = oldfit[0];

        PulsarSetREAL8ParamErr( target, this->name, err, newfit );
      }
    }

    this = next;
    if ( this ) {
      next = this->next;
    }
  }

  return;
}


/* create functions for converting par file input e.g. from strings to floats, converting times, converting units */
enum {
  CONVFLOAT = 0,
  CONVINT,
  CONVSTRING,
  CONVHMS,
  CONVDMS,
  CONVMJD,
  CONVBINUNITS
};


#define DEFINE_CONV_FACTOR_FUNCTION( name, convfactor, type ) \
  void ParConv ## name ( const CHAR *inval, void *out ){ \
    CHAR *in = XLALStringDuplicate( inval ); \
    if ( type != CONVSTRING && type != CONVINT ) { \
      REAL8 x; \
      if ( type == CONVFLOAT || type == CONVBINUNITS ){ \
        /* check for exponents expressed as 'D/d' rather than 'E/e' and switch */ \
        for ( INT4 i = 0; i < (INT4)strlen(in); i++ ) { \
          if ( in[i] == 'D' || in[i] == 'd' ) { in[i] = 'e'; } \
        }  \
        x = atof(in); \
        if ( type == CONVBINUNITS ){ \
          /* TEMPO2 checks if this is > 1e-7 then it's in units of 1e-12, so needs converting */ \
          if ( fabs( x ) > 1e-7 ) { x *= 1.e-12; } \
          /* some of the parameter files in the ATNF catalogue have values of EDOT that are stupidly large e.g. O(1e33). \
           * These can cause the time delay routines to fail, so if values of EDOT are greater than 10000 ignore them and \
           * set it to zero */ \
          if ( x > 10000. ) { x = 0.; } \
        } \
        else{ x *= convfactor; } \
      } \
      else if ( type == CONVHMS ) { /* convert to radians from hh:mm:ss.s format */ \
        XLALTranslateHMStoRAD( &x, in );                              \
      } \
      else if ( type == CONVDMS ) { /* convert to radians from hh:mm:ss.s format */ \
        XLALTranslateDMStoRAD( &x, in );                                  \
      } \
      else if ( type == CONVMJD ) { /* convert an MJD to a GPS time */ \
         x = XLALTTMJDtoGPS( atof(in) ); \
      } \
      memcpy(out, &x, sizeof(REAL8) ); \
    } \
    else if ( type == CONVSTRING ){ memcpy( out, in, strlen(in)+1 ); } \
    else if ( type == CONVINT ) { /* convert to an integer */ \
      UINT4 x; \
      x = (UINT4)atoi(in); \
      memcpy(out, &x, sizeof(UINT4) ); \
    } \
    XLALFree( in ); \
  }


DEFINE_CONV_FACTOR_FUNCTION( ToString, 0, CONVSTRING )
DEFINE_CONV_FACTOR_FUNCTION( ToFloat, 1., CONVFLOAT ) /* just convert string to float */
DEFINE_CONV_FACTOR_FUNCTION( ToInt, 0, CONVINT ) /* just convert string to float */
DEFINE_CONV_FACTOR_FUNCTION( DegsToRads, LAL_PI_180, CONVFLOAT ) /* convert degrees to radians */
DEFINE_CONV_FACTOR_FUNCTION( MasPerYrToRadPerSec, LAL_PI_180 / ( 3600.e3 * 365.25 * 86400. ), CONVFLOAT ) /* convert milliarcseconds/yr to radians/s */
DEFINE_CONV_FACTOR_FUNCTION( SecsToRads, LAL_TWOPI / ( 24.*3600. ), CONVFLOAT ) /* convert seconds to radians */
DEFINE_CONV_FACTOR_FUNCTION( ArcsecsToRads, LAL_PI_180 / 3600., CONVFLOAT ) /* convert arcseconds to radians */
DEFINE_CONV_FACTOR_FUNCTION( MasToRads, LAL_PI_180 / 3600.0e3, CONVFLOAT ) /* convert milliarcseconds to radians */
DEFINE_CONV_FACTOR_FUNCTION( InvArcsecsToInvRads, 3600. / LAL_PI_180, CONVFLOAT ) /* convert 1/arcsec to 1/rads */
DEFINE_CONV_FACTOR_FUNCTION( DaysToSecs, DAYSTOSECS, CONVFLOAT ) /* convert days to seconds */
DEFINE_CONV_FACTOR_FUNCTION( KpcToMetres, LAL_PC_SI * 1.e3, CONVFLOAT ) /* convert kiloparsecs to metres */
DEFINE_CONV_FACTOR_FUNCTION( BinaryUnits, 1.e-12, CONVBINUNITS ) /* convert certain binary units as defined in TEMPO2 with factor */
DEFINE_CONV_FACTOR_FUNCTION( MJDToGPS, 0, CONVMJD ) /* convert from MJD to GPS time */
DEFINE_CONV_FACTOR_FUNCTION( DegPerYrToRadPerSec, LAL_PI_180 / ( 365.25 * DAYSTOSECS ), CONVFLOAT ) /* convert degs/year to rads/s */
DEFINE_CONV_FACTOR_FUNCTION( SolarMassToKg, LALPULSAR_TEMPO2_MSUN_SI, CONVFLOAT ) /* convert solar masses to kg */
DEFINE_CONV_FACTOR_FUNCTION( RAToRads, 0, CONVHMS ) /* convert right ascension to radians */
DEFINE_CONV_FACTOR_FUNCTION( DecToRads, 0, CONVDMS ) /* convert declination to radians */
DEFINE_CONV_FACTOR_FUNCTION( MicrosecToSec, 1.e-6, CONVFLOAT ) /* convert microseconds to seconds */

/** Function type definition for a conversion function */
typedef void ( *ParamConversionFunc )( const CHAR *in, void *out );

/** A strcuture to contain all possible pulsar parameters that can be read in from a par file, and define
 * the conversion function and type used for each */
typedef struct tagParConversion {
  CHAR name[PULSAR_PARNAME_MAX];    /** Parameter name */
  ParamConversionFunc convfunc;     /** Conversion function from string to required value */
  ParamConversionFunc converrfunc;  /** Conversion function for error from string to value */
  PulsarParamType ptype;            /** Parameter type */
} ParConversion;


#define NUM_PARS 132 /* number of allowed parameters */

/** Initialise conversion structure with most allowed TEMPO2 parameter names and conversion functions
 * (convert all read in parameters to SI units where necessary). See http://arxiv.org/abs/astro-ph/0603381 and
 * http://arxiv.org/abs/astro-ph/0603381 for parameter definitions.
 *
 * If requiring a new parameter to be read in in should be added to this structure and \c NUM_PARS should be
 * incremented.
 */
ParConversion pc[NUM_PARS] = {
  { .name = "F", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8Vector_t }, /* vector containing frequency and derivatives (Hz, Hz/s, Hz/s^2 ...) */
  { .name = "P", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8Vector_t }, /* vector containing period and derivatives (s, s/s, s/s^2 ...) */
  { .name = "DIST", .convfunc = ParConvKpcToMetres, .converrfunc = ParConvKpcToMetres, .ptype = PULSARTYPE_REAL8_t }, /* distance to pulsar in metres */
  { .name = "PX", .convfunc = ParConvMasToRads, .converrfunc = ParConvMasToRads, .ptype = PULSARTYPE_REAL8_t }, /* parallax (converted to radians) */
  { .name = "DM", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8_t }, /* dispersion measure */
  { .name = "DM1", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8_t }, /*first derivative of the dispersion measure */

  /* position parameters */
  { .name = "RA", .convfunc = ParConvRAToRads, .converrfunc = ParConvSecsToRads, .ptype = PULSARTYPE_REAL8_t }, /* right ascension (converted to radians) */
  { .name = "RAJ", .convfunc = ParConvRAToRads, .converrfunc = ParConvSecsToRads, .ptype = PULSARTYPE_REAL8_t }, /* right ascension (converted to radians) */
  { .name = "DEC", .convfunc = ParConvDecToRads, .converrfunc = ParConvArcsecsToRads, .ptype = PULSARTYPE_REAL8_t }, /* declination (converted to radians) */
  { .name = "DECJ", .convfunc = ParConvDecToRads, .converrfunc = ParConvArcsecsToRads, .ptype = PULSARTYPE_REAL8_t }, /* declination (converted to radians) */
  { .name = "PMRA", .convfunc = ParConvMasPerYrToRadPerSec, .converrfunc = ParConvMasPerYrToRadPerSec, .ptype = PULSARTYPE_REAL8_t }, /* proper motion in right ascension (converted to radians/s) */
  { .name = "PMDEC", .convfunc = ParConvMasPerYrToRadPerSec, .converrfunc = ParConvMasPerYrToRadPerSec, .ptype = PULSARTYPE_REAL8_t }, /* proper motion in declination (converted to radians/s) */
  { .name = "ELONG", .convfunc = ParConvDegsToRads, .converrfunc = ParConvDegsToRads, .ptype = PULSARTYPE_REAL8_t }, /* ecliptic longitude (converted from degs to rads) */
  { .name = "ELAT", .convfunc = ParConvDegsToRads, .converrfunc = ParConvDegsToRads, .ptype = PULSARTYPE_REAL8_t }, /* ecliptic latitude (converted from degs to rads) */
  { .name = "PMELONG", .convfunc = ParConvMasPerYrToRadPerSec, .converrfunc = ParConvMasPerYrToRadPerSec, .ptype = PULSARTYPE_REAL8_t }, /* proper motion in ecliptic longitude (converted to radians/s) */
  { .name = "PMELAT", .convfunc = ParConvMasPerYrToRadPerSec, .converrfunc = ParConvMasPerYrToRadPerSec, .ptype = PULSARTYPE_REAL8_t }, /* proper motion in ecliptic latitude (converted to radians/s) */
  { .name = "BETA", .convfunc = ParConvDegsToRads, .converrfunc = ParConvDegsToRads, .ptype = PULSARTYPE_REAL8_t }, /* an alias for ELONG also giving ecliptic latitude (converted from degs to rads) */
  { .name = "LAMBDA", .convfunc = ParConvDegsToRads, .converrfunc = ParConvDegsToRads, .ptype = PULSARTYPE_REAL8_t }, /* an alias for ELAT also giving ecliptic longitude (converted from degs to rads) */
  { .name = "PMBETA", .convfunc = ParConvMasPerYrToRadPerSec, .converrfunc = ParConvMasPerYrToRadPerSec, .ptype = PULSARTYPE_REAL8_t }, /* proper motion in ecliptic latitude (converted to radians/s) */
  { .name = "PMLAMBDA", .convfunc = ParConvMasPerYrToRadPerSec, .converrfunc = ParConvMasPerYrToRadPerSec, .ptype = PULSARTYPE_REAL8_t }, /* proper motion in ecliptic longitude (converted to radians/s) */

  /* epoch parameters */
  { .name = "PEPOCH", .convfunc = ParConvMJDToGPS, .converrfunc = ParConvDaysToSecs, .ptype = PULSARTYPE_REAL8_t }, /* period epoch (saved as GPS time) */
  { .name = "POSEPOCH", .convfunc = ParConvMJDToGPS, .converrfunc = ParConvDaysToSecs, .ptype = PULSARTYPE_REAL8_t }, /* position epoch (saved as GPS time) */
  { .name = "DMEPOCH", .convfunc = ParConvMJDToGPS, .converrfunc = ParConvDaysToSecs, .ptype = PULSARTYPE_REAL8_t }, /* dispersion measure epoch (saved as GPS time) */

  /* glitch parameters */
  { .name = "GLEP", .convfunc = ParConvMJDToGPS, .converrfunc = ParConvDaysToSecs, .ptype = PULSARTYPE_REAL8Vector_t }, /* glitch epochs (saved as GPS times) */
  { .name = "GLPH", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8Vector_t }, /* glitch phase offsets (rotational phase) */
  { .name = "GLF0", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8Vector_t }, /* glitch frequency offsets (Hz) */
  { .name = "GLF1", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8Vector_t }, /* glitch frequency derivative offsets (Hz/s) */
  { .name = "GLF2", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8Vector_t }, /* glitch second frequency derivative offsets (Hz/s^2) */
  { .name = "GLF0D", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8Vector_t }, /* glitch decaying frequency component (Hz) */
  { .name = "GLTD", .convfunc = ParConvDaysToSecs, .converrfunc = ParConvDaysToSecs, .ptype = PULSARTYPE_REAL8Vector_t }, /* glitch decay time (seconds) */

  /* string parameters */
  { .name = "NAME", .convfunc = ParConvToString, .converrfunc = NULL, .ptype = PULSARTYPE_string_t }, /* pulsar name */
  { .name = "PSR", .convfunc = ParConvToString, .converrfunc = NULL, .ptype = PULSARTYPE_string_t }, /* pulsar name */
  { .name = "PSRJ", .convfunc = ParConvToString, .converrfunc = NULL, .ptype = PULSARTYPE_string_t }, /* pulsar J name */
  { .name = "PSRB", .convfunc = ParConvToString, .converrfunc = NULL, .ptype = PULSARTYPE_string_t }, /* pulsar B name */

  /* binary parameters */
  { .name = "BINARY", .convfunc = ParConvToString, .converrfunc = NULL, .ptype = PULSARTYPE_string_t }, /* binary model type */
  { .name = "A1", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8_t }, /* projected semi-major axis (light seconds) */
  { .name = "OM", .convfunc = ParConvDegsToRads, .converrfunc = ParConvDegsToRads, .ptype = PULSARTYPE_REAL8_t }, /* angle of periastron (convert degrees to radians) */
  { .name = "ECC", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8_t }, /* eccentricity */
  { .name = "PB", .convfunc = ParConvDaysToSecs, .converrfunc = ParConvDaysToSecs, .ptype = PULSARTYPE_REAL8_t }, /* orbital period (convert days to seconds) */
  { .name = "T0", .convfunc = ParConvMJDToGPS, .converrfunc = ParConvDaysToSecs, .ptype = PULSARTYPE_REAL8_t }, /* time of perisatron (GPS) */
  { .name = "TASC", .convfunc = ParConvMJDToGPS, .converrfunc = ParConvDaysToSecs, .ptype = PULSARTYPE_REAL8_t }, /* time ascending noise for ELL1 model (GPS) */
  { .name = "EPS1", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8_t }, /* e*sin(w0) for ELL1 model */
  { .name = "EPS2", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8_t }, /* e*cos(w0) for ELL1 model */
  { .name = "GAMMA", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8_t }, /* relativistic parameter */
  { .name = "OMDOT", .convfunc = ParConvDegPerYrToRadPerSec, .converrfunc = ParConvDegPerYrToRadPerSec, .ptype = PULSARTYPE_REAL8_t }, /* angle of periastron time derivative (degs/year converted to rad/s) */
  { .name = "XDOT", .convfunc = ParConvBinaryUnits, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8_t }, /* project semi-major axis time derivative (light sec/sec) */
  { .name = "PBDOT", .convfunc = ParConvBinaryUnits, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8_t }, /* period time derivative */
  { .name = "EDOT", .convfunc = ParConvBinaryUnits, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8_t }, /* eccentricity time derivative (1/s) */
  { .name = "EPS1DOT", .convfunc = ParConvBinaryUnits, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8_t },
  { .name = "EPS2DOT", .convfunc = ParConvBinaryUnits, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8_t },
  { .name = "XPBDOT", .convfunc = ParConvBinaryUnits, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8_t },
  { .name = "SINI", .convfunc = ParConvToString, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_string_t }, /* sin of inclination angle */
  { .name = "MTOT", .convfunc = ParConvSolarMassToKg, .converrfunc = ParConvSolarMassToKg, .ptype = PULSARTYPE_REAL8_t }, /* total system mass (convert solar masses to kg) */
  { .name = "M2", .convfunc = ParConvSolarMassToKg, .converrfunc = ParConvSolarMassToKg, .ptype = PULSARTYPE_REAL8_t }, /* binary companion mass (convert solar masses to kg) */
  { .name = "DR", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8_t },
  { .name = "DTHETA", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8_t },
  { .name = "SHAPMAX", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8_t }, /* used for DDS model */

  /* multiple system terms for BT1P and BT2P models */
  { .name = "A1_2", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8_t },
  { .name = "A1_3", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8_t },
  { .name = "OM_2", .convfunc = ParConvDegsToRads, .converrfunc = ParConvDegsToRads, .ptype = PULSARTYPE_REAL8_t },
  { .name = "OM_3", .convfunc = ParConvDegsToRads, .converrfunc = ParConvDegsToRads, .ptype = PULSARTYPE_REAL8_t },
  { .name = "PB_2", .convfunc = ParConvDaysToSecs, .converrfunc = ParConvDaysToSecs, .ptype = PULSARTYPE_REAL8_t },
  { .name = "PB_3", .convfunc = ParConvDaysToSecs, .converrfunc = ParConvDaysToSecs, .ptype = PULSARTYPE_REAL8_t },
  { .name = "T0_2", .convfunc = ParConvMJDToGPS, .converrfunc = ParConvDaysToSecs, .ptype = PULSARTYPE_REAL8_t },
  { .name = "T0_3", .convfunc = ParConvMJDToGPS, .converrfunc = ParConvDaysToSecs, .ptype = PULSARTYPE_REAL8_t },
  { .name = "ECC_2", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8_t },
  { .name = "ECC_3", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8_t },
  { .name = "FB", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8Vector_t }, /* orbital frequency components for the BTX model */

  /* Aberration delay parameters */
  { .name = "A0", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8_t },
  { .name = "B0", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8_t },

  /* Kopeikin terms */
  { .name = "D_AOP", .convfunc = ParConvInvArcsecsToInvRads, .converrfunc = NULL, .ptype = PULSARTYPE_REAL8_t }, /* inverse arsecs converted to inverse radians */
  { .name = "KIN", .convfunc = ParConvDegsToRads, .converrfunc = NULL, .ptype = PULSARTYPE_REAL8_t },
  { .name = "KOM", .convfunc = ParConvDegsToRads, .converrfunc = NULL, .ptype = PULSARTYPE_REAL8_t },

  /* FITWAVES parameters */
  { .name = "WAVE_OM", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8_t }, /* fundamental frequency (Hz) of sinusoids used for whitening */
  { .name = "WAVEEPOCH", .convfunc = ParConvMJDToGPS, .converrfunc = NULL, .ptype = PULSARTYPE_REAL8_t }, /* frequency epoch (converted from MJD to GPS) */
  { .name = "WAVESIN", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8Vector_t }, /* amplitudes of the sine terms for the kth sinusoids */
  { .name = "WAVECOS", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8Vector_t }, /* amplitudes of the cosine terms for the kth sinusoids */

  /* TEMPO parameters */
  { .name = "EPHEM", .convfunc = ParConvToString, .converrfunc = NULL, .ptype = PULSARTYPE_string_t }, /* ephemeris type e.g. DE405 */
  { .name = "UNITS", .convfunc = ParConvToString, .converrfunc = NULL, .ptype = PULSARTYPE_string_t }, /* TEMPO2 units e.g. TDB */
  { .name = "START", .convfunc = ParConvMJDToGPS, .converrfunc = NULL, .ptype = PULSARTYPE_REAL8_t }, /* start of observations */
  { .name = "FINISH", .convfunc = ParConvMJDToGPS, .converrfunc = NULL, .ptype = PULSARTYPE_REAL8_t }, /* end of observations */
  { .name = "NTOA", .convfunc = ParConvToInt, .converrfunc = NULL, .ptype = PULSARTYPE_UINT4_t }, /* number of TOAs in observation */
  { .name = "TRES", .convfunc = ParConvMicrosecToSec, .converrfunc = NULL, .ptype = PULSARTYPE_REAL8_t }, /* timing residual (convert microseconds to seconds) */
  { .name = "CLK", .convfunc = ParConvToString, .converrfunc = NULL, .ptype = PULSARTYPE_string_t }, /* The observatory clock */
  { .name = "T2CMETHOD", .convfunc = ParConvToString, .converrfunc = NULL, .ptype = PULSARTYPE_string_t }, /* Method for transforming from terrestrial to celestial frame */
  { .name = "TIMEEPH", .convfunc = ParConvToString, .converrfunc = NULL, .ptype = PULSARTYPE_string_t }, /* Which time ephemeris to use */

  /* GW parameters */
  { .name = "H0", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8_t }, /* gravitational wave amplitude */
  { .name = "APLUS", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8_t }, /* plus polarisation component of GW amplitude */
  { .name = "ACROSS", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8_t }, /* cross polarisation component of GW amplitude */
  { .name = "PHI0", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8_t }, /* gravitational wave initial phase (radians) */
  { .name = "PSI", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8_t }, /* gravitational wave polarisation angle (radians) */
  { .name = "COSIOTA", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8_t }, /* cosine of source inclination angle */
  { .name = "IOTA", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8_t }, /* the source inclination angle in radians */
  { .name = "C22", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8_t }, /* GW amplitude of C22 component */
  { .name = "C21", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8_t }, /* GW amplitude of C21 component */
  { .name = "PHI22", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8_t }, /* initial phase of C22 component (radians) */
  { .name = "PHI21", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8_t }, /* initial phase of C21 component (radians) */
  { .name = "CGW", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8_t }, /* speed of gravitational waves as a fraction of the speed of light */
  { .name = "LAMBDAPIN", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8_t }, /* parameters from http://uk.arxiv.org/abs/0909.4035 */
  { .name = "COSTHETA", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8_t },
  { .name = "THETA", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8_t },
  { .name = "I21", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8_t },
  { .name = "I31", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8_t },
  { .name = "Q22", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8_t }, /* the l=m=2 mass quadrupole in kg m^2 */

  /* GW non-GR polarisation mode amplitude parameters (assuming emission at twice the rotation frequency) */
  { .name = "HPLUS", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8_t }, /* GW tensor plus polarisation amplitude */
  { .name = "HCROSS", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8_t }, /* GW tensor cross polarisation amplitude */
  { .name = "PSITENSOR", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8_t }, /* GW tensor angle polarisation */
  { .name = "PHI0TENSOR", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8_t }, /* initial phase for tensor modes */
  { .name = "HSCALARB", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8_t }, /* GW scalar breathing mode polarisation amplitude */
  { .name = "HSCALARL", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8_t }, /* GW scalar longitudinal polarisation amplitude */
  { .name = "PSISCALAR", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8_t }, /* GW scalar angle polarisation */
  { .name = "PHI0SCALAR", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8_t }, /* initial phase for scalar modes */
  { .name = "HVECTORX", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8_t }, /* GW vector x-mode amplitude */
  { .name = "HVECTORY", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8_t }, /* GW vector y-mode amplitude */
  { .name = "PSIVECTOR", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8_t }, /* GW vector angle polarisation */
  { .name = "PHI0VECTOR", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8_t }, /* GW vector polarisation initial phase */

  /* GW non-GR polarisation mode amplitude parameters (assuming emission at the rotation frequency) */
  { .name = "H0_F", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8_t }, /* GR 21 polarisation amplitude */
  { .name = "HPLUS_F", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8_t }, /* GW tensor plus polarisation amplitude */
  { .name = "HCROSS_F", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8_t }, /* GW tensor cross polarisation amplitude */
  { .name = "PSITENSOR_F", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8_t }, /* GW tensor angle polarisation */
  { .name = "PHI0TENSOR_F", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8_t }, /* initial phase for tensor modes */
  { .name = "HSCALARB_F", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8_t }, /* GW scalar breathing mode polarisation amplitude */
  { .name = "HSCALARL_F", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8_t }, /* GW scalar longitudinal polarisation amplitude */
  { .name = "PSISCALAR_F", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8_t }, /* GW scalar angle polarisation */
  { .name = "PHI0SCALAR_F", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8_t }, /* initial phase for scalar modes */
  { .name = "HVECTORX_F", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8_t }, /* GW vector x-mode amplitude */
  { .name = "HVECTORY_F", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8_t }, /* GW vector y-mode amplitude */
  { .name = "PSIVECTOR_F", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8_t }, /* GW vector angle polarisation */
  { .name = "PHI0VECTOR_F", .convfunc = ParConvToFloat, .converrfunc = ParConvToFloat, .ptype = PULSARTYPE_REAL8_t }, /* GW vector polarisation initial phase */

  /* transient signal parameters */
  { .name = "TRANSIENTWINDOWTYPE", .convfunc = ParConvToString, .converrfunc = NULL, .ptype = PULSARTYPE_string_t }, /* window type, currently "RECT" (rectangular window) or "EXP" (exponentially decaying window) */
  { .name = "TRANSIENTSTARTTIME", .convfunc = ParConvMJDToGPS, .converrfunc = NULL, .ptype = PULSARTYPE_REAL8_t }, /* transient start time in MJD */
  { .name = "TRANSIENTTAU", .convfunc = ParConvDaysToSecs, .converrfunc = NULL, .ptype = PULSARTYPE_REAL8_t } /* transient decay time scale (duration for "RECT" and decay constant for "EXP") */
};

/** \brief Parse a single line from a pulsar parameter file
 *
 * This will parse a line from the TEMPO-style pulsar parameter file containing the
 * parameter given by \c name. The parameter will be added to the \c par structure.
 *
 * This function can only be used internally.
 */
static INT4 ParseParLine( PulsarParameters *par, const CHAR *name, FILE *fp )
{
  INT4 nread = 0; /* number of values read from line */
  CHAR str[PULSAR_PARNAME_MAX];
  /* three potential values on the line */
  CHAR str1[PULSAR_PARNAME_MAX], str2[PULSAR_PARNAME_MAX], str3[PULSAR_PARNAME_MAX];
  INT4 i = 0, tmpnum = 0;
  CHAR *nname = NULL; /* duplicate of name */

  if ( par == NULL ) {
    XLAL_ERROR( XLAL_EINVAL, "Error... PulsarParameter structure is not initialised!\n" );
  }
  if ( name == NULL ) {
    XLAL_ERROR( XLAL_EINVAL, "Error... parameter name is not set!\n" );
  }
  if ( fp == NULL ) {
    XLAL_ERROR( XLAL_EINVAL, "Error... parameter file pointer is NULL!\n" );
  }

  /* parse the line from the fp pointer */
  if ( fgets( str, PULSAR_PARNAME_MAX, fp ) == NULL ) {
    XLAL_PRINT_WARNING( "No value for parameter %s", name );
    return XLAL_SUCCESS;
  }

  /* scan the line for values */
  nread = sscanf( str, "%s %s %s", str1, str2, str3 );

  if ( !nread ) {
    XLAL_PRINT_WARNING( "No value for parameter %s", name );
    return XLAL_SUCCESS;
  }

  /* check for parameters with more than one possible name */
  if ( !strcmp( name, "E" ) ) {
    nname = XLALStringDuplicate( "ECC" );
  } else if ( !strcmp( name, "E_1" ) ) {
    nname = XLALStringDuplicate( "ECC_1" );
  } else if ( !strcmp( name, "E_2" ) ) {
    nname = XLALStringDuplicate( "ECC_2" );
  } else if ( !strcmp( name, "E_3" ) ) {
    nname = XLALStringDuplicate( "ECC_3" );
  } else {
    nname = XLALStringDuplicate( name );
  }

  /* perform parameter dependent inputs */
  for ( i = 0; i < NUM_PARS; i++ ) {
    /* this has to be hard-coded for the F, WAVE, FB and glitch (GL) vector parameters */

    /* check for parameters requiring placement in vectors */
    INT4 isfreq = ( !strncmp( nname, "F", 1 ) && !strcmp( "F", pc[i].name ) && sscanf( nname + strlen( "F" ), "%d",  &tmpnum ) == 1 );
    INT4 isperiod = ( !strncmp( nname, "P", 1 ) && !strcmp( "P", pc[i].name ) && sscanf( nname + strlen( "P" ), "%d",  &tmpnum ) == 1 );
    INT4 iswave = ( ( !strncmp( nname, "WAVE", 4 ) && ( !strcmp( "WAVESIN", pc[i].name ) || !strcmp( "WAVECOS", pc[i].name ) ) ) && strcmp( "WAVE_OM", nname ) && strcmp( "WAVEEPOCH", nname ) );
    INT4 isfb = ( !strncmp( nname, "FB", 2 ) && !strcmp( "FB", pc[i].name ) );
    INT4 isgl = ( !strncmp( nname, "GL", 2 ) && strstr( nname, pc[i].name ) && !strncmp( nname, pc[i].name, strlen( nname ) - 2 ) );

    if ( !strcmp( nname, pc[i].name ) || isfreq || isperiod || iswave || isfb || isgl ) {
      UINT4 num = 0, nsize = 0;

      if ( pc[i].convfunc == NULL ) {
        XLAL_PRINT_WARNING( "No conversion function for parameter %s. Skipping parameter.\n", pc[i].name );
        return XLAL_SUCCESS;
      }

      /* add parameter */
      if ( isfreq || isperiod || isfb ) { /* add frequency/period and frequency/period derivative values, or FB values */
        REAL8Vector *ptr = NULL, *eptr = NULL;
        UINT4 *fitFlag = NULL;

        if ( isfb ) {
          if ( strlen( nname ) > strlen( "FB" ) ) {
            if ( sscanf( nname + strlen( "FB" ), "%d",  &num ) != 1 ) {
              XLAL_ERROR( XLAL_EINVAL, "Error...problem reading %s number from par file.\n", nname );
            }
          }
        } else if ( isfreq ) {
          if ( strlen( nname ) > strlen( "F" ) ) {
            if ( sscanf( nname + strlen( "F" ), "%d", &num ) != 1 ) {
              XLAL_ERROR( XLAL_EINVAL, "Error...problem reading %s number from par file.\n", nname );
            }
          }
        } else {
          if ( strlen( nname ) > strlen( "P" ) ) {
            if ( sscanf( nname + strlen( "P" ), "%d", &num ) != 1 ) {
              XLAL_ERROR( XLAL_EINVAL, "Error...problem reading %s number from par file.\n", nname );
            }
          }
        }

        void *val = ( void * )XLALMalloc( PulsarTypeSize[PULSARTYPE_REAL8_t] );

        pc[i].convfunc( str1, val );

        nsize = num + 1;
        if ( PulsarCheckParam( par, pc[i].name ) ) {
          // count any previous values
          const REAL8Vector *cptr = PulsarGetREAL8VectorParam( par, pc[i].name );
          nsize = ( cptr->length > num + 1 ) ? cptr->length : ( num + 1 );
        }

        // pointer for parameter values
        ptr = XLALCreateREAL8Vector( nsize );
        memset( ptr->data, 0, sizeof( REAL8 )*nsize );

        // pointers for error and fit flags
        eptr = XLALCreateREAL8Vector( nsize );
        memset( eptr->data, 0, sizeof( REAL8 )*nsize );       // initialise with zeros
        fitFlag = ( UINT4 * )XLALCalloc( nsize, sizeof( UINT4 ) ); // set fit values to zero (no-fit) by default

        if ( PulsarCheckParam( par, pc[i].name ) ) {
          // copy any previous values
          const REAL8Vector *cptr = PulsarGetREAL8VectorParam( par, pc[i].name );
          memcpy( ptr->data, cptr->data, cptr->length * sizeof( REAL8 ) );

          /* copy any error values, as updating REAL8Vector parameters destroys the original PulsarParam structure */
          const REAL8Vector *ceptr = PulsarGetREAL8VectorParamErr( par, pc[i].name );
          const UINT4 *ff = PulsarGetParamFitFlag( par, pc[i].name );

          // copy any previous error values
          memcpy( eptr->data, ceptr->data, sizeof( REAL8 )*cptr->length );
          memcpy( fitFlag, ff, sizeof( UINT4 )*cptr->length );
        }

        ptr->data[num] = *( REAL8 * )val;

        // NOTE: this will destroy the previous held PulsarParam structure for this parameter
        PulsarAddREAL8VectorParam( par, pc[i].name, ( const REAL8Vector * )ptr );

        // (re-)add errors
        PulsarSetREAL8VectorParamErr( par, pc[i].name, ( const REAL8Vector * )eptr, ( const UINT4 * )fitFlag );

        XLALFree( val );
        XLALDestroyREAL8Vector( ptr );
        XLALDestroyREAL8Vector( eptr );
        XLALFree( fitFlag );
      } else if ( iswave && nread == 2 ) { /* add WAVE values */
        REAL8Vector *ptr1 = NULL, *ptr2 = NULL;

        if ( strlen( nname ) > strlen( "WAVE" ) ) {
          if ( sscanf( nname + strlen( "WAVE" ), "%d",  &num ) != 1 ) {
            XLAL_ERROR( XLAL_EINVAL, "Error...problem reading %s number from par file.\n", nname );
          }
        } else {
          break;
        }

        num--; /* WAVE values start from WAVE1, so subtract 1 from the num for the vector index */

        void *val1 = ( void * )XLALMalloc( PulsarTypeSize[PULSARTYPE_REAL8_t] );
        void *val2 = ( void * )XLALMalloc( PulsarTypeSize[PULSARTYPE_REAL8_t] );

        pc[i].convfunc( str1, val1 );
        pc[i].convfunc( str2, val2 );

        nsize = num + 1;
        if ( PulsarCheckParam( par, "WAVESIN" ) && PulsarCheckParam( par, "WAVECOS" ) ) {
          // count number of values
          const REAL8Vector *cptr = PulsarGetREAL8VectorParam( par, pc[i].name );

          nsize = ( cptr->length > num + 1 ) ? cptr->length : ( num + 1 );
        }

        // pointers for parameter values
        ptr1 = XLALCreateREAL8Vector( nsize );
        ptr2 = XLALCreateREAL8Vector( nsize );
        memset( ptr1->data, 0, sizeof( REAL8 )*nsize );
        memset( ptr2->data, 0, sizeof( REAL8 )*nsize );

        if ( PulsarCheckParam( par, "WAVESIN" ) && PulsarCheckParam( par, "WAVECOS" ) ) {
          const REAL8Vector *cptr1 = PulsarGetREAL8VectorParam( par, "WAVESIN" );
          const REAL8Vector *cptr2 = PulsarGetREAL8VectorParam( par, "WAVECOS" );

          // copy previous values
          memcpy( ptr1->data, cptr1->data, ptr1->length * sizeof( REAL8 ) );
          memcpy( ptr2->data, cptr2->data, ptr2->length * sizeof( REAL8 ) );
        }

        ptr1->data[num] = *( REAL8 * )val1;
        ptr2->data[num] = *( REAL8 * )val2;

        PulsarAddREAL8VectorParam( par, "WAVESIN", ( const REAL8Vector * )ptr1 );
        PulsarAddREAL8VectorParam( par, "WAVECOS", ( const REAL8Vector * )ptr2 );

        XLALFree( val1 );
        XLALFree( val2 );
        XLALDestroyREAL8Vector( ptr1 );
        XLALDestroyREAL8Vector( ptr2 );

        /* there are no errors on the wave parameters, so set defaults of zero and then break */
        UINT4 *fitFlag = ( UINT4 * )XLALCalloc( nsize, sizeof( UINT4 ) );
        REAL8Vector *eptr = XLALCreateREAL8Vector( nsize );
        memset( eptr->data, 0, sizeof( REAL8 )*nsize );
        PulsarSetREAL8VectorParamErr( par, "WAVESIN", ( const REAL8Vector * )eptr, ( const UINT4 * )fitFlag );
        PulsarSetREAL8VectorParamErr( par, "WAVECOS", ( const REAL8Vector * )eptr, ( const UINT4 * )fitFlag );
        XLALDestroyREAL8Vector( eptr );
        XLALFree( fitFlag );

        break;
      } else if ( isgl ) { /* get glitch parameters */
        REAL8Vector *ptr = NULL, *eptr = NULL;
        UINT4 *fitFlag = NULL;

        if ( sscanf( strstr( nname, "_" ) + 1, "%d", &num ) != 1 ) {
          XLAL_ERROR( XLAL_EINVAL, "Error...problem reading %s number from par file.\n", nname );
        }

        num--; /* GL values start from e.g. GLF0_1, so subtract 1 from the num for the vector index */

        void *val = ( void * )XLALMalloc( PulsarTypeSize[PULSARTYPE_REAL8_t] );

        pc[i].convfunc( str1, val );

        nsize = num + 1;
        if ( PulsarCheckParam( par, pc[i].name ) ) {
          // count previous values
          const REAL8Vector *cptr = PulsarGetREAL8VectorParam( par, pc[i].name );
          nsize = ( cptr->length > num + 1 ) ? cptr->length : ( num + 1 );
        }

        // pointer for parameter values
        ptr = XLALCreateREAL8Vector( nsize );
        memset( ptr->data, 0, sizeof( REAL8 )*nsize );

        // pointers for error and fit flags
        eptr = XLALCreateREAL8Vector( nsize );
        memset( eptr->data, 0, sizeof( REAL8 )*nsize );       // initialise with zeros
        fitFlag = ( UINT4 * )XLALCalloc( nsize, sizeof( UINT4 ) ); // set fit values to zero (no-fit) by default

        if ( PulsarCheckParam( par, pc[i].name ) ) {
          // copy any previous values
          const REAL8Vector *cptr = PulsarGetREAL8VectorParam( par, pc[i].name );
          memcpy( ptr->data, cptr->data, cptr->length * sizeof( REAL8 ) );

          /* copy any error values, as updating REAL8Vector parameters destroys the original PulsarParam structure */
          const REAL8Vector *ceptr = PulsarGetREAL8VectorParamErr( par, pc[i].name );
          const UINT4 *ff = PulsarGetParamFitFlag( par, pc[i].name );

          // copy any previous error values
          memcpy( eptr->data, ceptr->data, sizeof( REAL8 )*cptr->length );
          memcpy( fitFlag, ff, sizeof( UINT4 )*cptr->length );
        }

        ptr->data[num] = *( REAL8 * )val;

        PulsarAddREAL8VectorParam( par, pc[i].name, ( const REAL8Vector * )ptr );

        // (re-)add errors
        PulsarSetREAL8VectorParamErr( par, pc[i].name, ( const REAL8Vector * )eptr, ( const UINT4 * )fitFlag );

        XLALFree( val );
        XLALDestroyREAL8Vector( ptr );
        XLALDestroyREAL8Vector( eptr );
        XLALFree( fitFlag );
      } else {
        if ( pc[i].ptype != PULSARTYPE_string_t ) {
          void *val = ( void * )XLALMalloc( PulsarTypeSize[pc[i].ptype] );
          pc[i].convfunc( str1, val );
          PulsarAddParam( par, pc[i].name, val, pc[i].ptype );
          XLALFree( val );
        } else {
          CHAR *val = XLALStringDuplicate( str1 );
          PulsarAddStringParam( par, pc[i].name, ( const CHAR * )val );
          XLALFree( val );
        }
      }

      /* check for error values */
      if ( ( nread == 2 && strcmp( str2, "1" ) ) || ( nread == 3 && !strcmp( str2, "1" ) ) ) {
        if ( pc[i].converrfunc == NULL ) {
          XLAL_PRINT_WARNING( "No conversion function for parameter %s error. No error being set.\n", pc[i].name );
          return XLAL_SUCCESS;
        }

        void *val = ( void * )XLALMalloc( PulsarTypeSize[PULSARTYPE_REAL8_t] );
        UINT4 isFit = 0;

        /* get the fit flag */
        if ( isfreq || isperiod || isfb || isgl ) { /* do this for REAL8Vector parameters */
          const REAL8Vector *ptr = PulsarGetREAL8VectorParamErr( par, pc[i].name );
          const UINT4 *fitFlag = PulsarGetParamFitFlag( par, pc[i].name );

          XLAL_CHECK( ptr != NULL, XLAL_EFUNC, "No default errors set for parameter %s", pc[i].name );

          // copy of parameter errors and fit flags
          REAL8Vector *eptr = XLALCreateREAL8Vector( ptr->length );
          UINT4 *ff = ( UINT4 * )XLALMalloc( sizeof( UINT4 ) * ptr->length );

          memcpy( eptr->data, ptr->data, sizeof( REAL8 )*ptr->length );
          memcpy( ff, fitFlag, sizeof( UINT4 )*ptr->length );

          if ( nread == 2 ) {
            pc[i].converrfunc( str2, val );
          } else {
            if ( !strcmp( str2, "1" ) ) {
              isFit = 1;  /* a fit flag is set to one */
            }
            pc[i].converrfunc( str3, val );
          }

          eptr->data[num] = *( REAL8 * )val;
          ff[num] = isFit;

          PulsarSetREAL8VectorParamErr( par, pc[i].name, ( const REAL8Vector * )eptr, ( const UINT4 * )ff );

          XLALDestroyREAL8Vector( eptr );
          XLALFree( ff );
        } else {
          if ( nread == 2 ) {
            pc[i].converrfunc( str2, val );
          } else {
            if ( !strcmp( str2, "1" ) ) {
              isFit = 1;
            }
            pc[i].converrfunc( str3, val );
          }

          PulsarSetREAL8ParamErr( par, pc[i].name, *( REAL8 * )val, isFit );
        }

        XLALFree( val );
      }

      break;
    }
  }

  XLALFree( nname );

  return XLAL_SUCCESS;
}


/* read in the pulsar parameter file */
PulsarParameters *XLALReadTEMPOParFile( const CHAR *pulsarAndPath )
{
  FILE *fp = NULL;
  CHAR str[PULSAR_PARNAME_MAX]; /* string to contain first value on line */

  PulsarParameters *par = XLALCalloc( sizeof( *par ), 1 );

  /* open file */
  if ( ( fp = fopen( pulsarAndPath, "r" ) ) == NULL ) {
    XLAL_PRINT_ERROR( "Error... Cannot open .par file %s\n", pulsarAndPath );
    XLALFree( par );
    XLAL_ERROR_NULL( XLAL_EIO );
  }

  int nread = 0; /* number of parameters read from a line in the par file */
  int UNUSED c;

  /* read in the par file */
  while ( !feof( fp ) ) {
    /* Read in a line from the parameter file */
    nread = fscanf( fp, "%s", str );
    if ( nread == 1 ) {
      /* check for comment line and skip to end of line */
      if ( str[0] == '#' ) {
        c = fscanf( fp, "%*[^\n]" );
        continue;
      }

      CHAR upperName[PULSAR_PARNAME_MAX];
      XLALStringCopy( upperName, str, PULSAR_PARNAME_MAX );
      strtoupper( upperName );

      if ( XLAL_SUCCESS != ParseParLine( par, upperName, fp ) ) {
        XLAL_PRINT_WARNING( "Parameter \"%s\" could not be successfully parsed from par file", str );
      }
    }
  }

  /* check for linked parameters SINI and KIN */
  if ( PulsarCheckParam( par, "SINI" ) ) {
    CHAR *sini = XLALStringDuplicate( PulsarGetStringParam( par, "SINI" ) );
    strtoupper( sini );

    REAL8 sinid;

    PulsarRemoveParam( par, "SINI" );

    if ( !strcmp( sini, "KIN" ) ) {
      if ( PulsarCheckParam( par, "KIN" ) ) {
        sinid = sin( PulsarGetREAL8Param( par, "KIN" ) );
        PulsarAddParam( par, "SINI", &sinid, PULSARTYPE_REAL8_t );
      } else {
        XLAL_PRINT_ERROR( "Error... KIN not set in .par file %s\n", pulsarAndPath );
        PulsarClearParams( par );
        XLALFree( par );
        XLAL_ERROR_NULL( XLAL_EIO );
      }
    } else {
      sinid = atof( sini );
      PulsarAddParam( par, "SINI", &sinid, PULSARTYPE_REAL8_t );
    }

    XLALFree( sini );
  }

  fclose( fp );

  return par;
}

/* NOTE: Convert this function to be more like readParfile.C in TEMPO2 - read
 * in a line at a time using fgets and make each parameter a structure */
void
XLALReadTEMPOParFileOrig( BinaryPulsarParams *output,
                          CHAR      *pulsarAndPath )
{
  XLAL_PRINT_DEPRECATION_WARNING( "XLALReadTEMPOParFile" );

  FILE *fp = NULL;
  CHAR val[500][40]; /* string array to hold all the read in values
                        500 strings of max 40 characters is enough */
  INT4 i = 0, j = 1, k;
  int UNUSED c;

  if ( output == ( BinaryPulsarParams * )NULL ) {
    XLAL_ERROR_VOID( XLAL_EFAULT );
  }

  output->name = NULL;
  output->jname = NULL;
  output->bname = NULL;

  output->model = NULL; /* set binary model to null - in case not a binary */

  /* set all output params to zero*/
  output->e = 0.0;    /* orbital eccentricity */
  output->Pb = 0.0;   /* orbital period (days) */
  output->w0 = 0.0;   /* longitude of periastron (deg) */
  output->x = 0.0;    /* projected semi-major axis/speed of light (light secs) */
  output->T0 = 0.0;   /* time of orbital periastron as measured in TDB (MJD) */

  output->e2 = 0.0;
  output->Pb2 = 0.0;
  output->w02 = 0.0;
  output->x2 = 0.0;
  output->T02 = 0.0;

  output->e3 = 0.0;
  output->Pb3 = 0.0;
  output->w03 = 0.0;
  output->x3 = 0.0;
  output->T03 = 0.0;

  output->xpbdot = 0.0; /* (10^-12) */

  output->eps1 = 0.0;     /* e*sin(w) */
  output->eps2 = 0.0;     /* e*cos(w) */
  output->eps1dot = 0.0;
  output->eps2dot = 0.0;
  output->Tasc = 0.0; /* time of the ascending node (used rather than T0) */

  output->fb = NULL;
  output->fbErr = NULL;
  output->nfb = 0;

  output->wdot = 0.0; /* precesion of longitude of periastron w = w0 + wdot(tb-T0) (degs/year) */
  output->gamma = 0.0; /* gravitational redshift and time dilation parameter (s)*/
  output->Pbdot = 0.0; /* rate of change of Pb (dimensionless 10^-12) */
  output->xdot = 0.0; /* rate of change of x(=asini/c) - optional (10^-12)*/
  output->edot = 0.0; /* rate of change of e (10^-12)*/

  output->s = 0.0;    /* Shapiro 'shape' parameter sin i */
  output->sstr = NULL;

  output->shapmax = 0.;

  /*output.r=0.0; Shapiro 'range' parameter */
  output->dr = 0.0;
  output->dth = 0.0;  /* (10^-6) */
  output->a0 = 0.0;
  output->b0 = 0.0; /* abberation delay parameters */

  output->M = 0.0;   /* M = m1 + m2 (m1 = pulsar mass, m2 = companion mass) */
  output->m2 = 0.0;  /* companion mass */

  output->f0 = 0.0;
  output->f1 = 0.0;
  output->f2 = 0.0;
  output->f3 = 0.0;
  output->f4 = 0.0;
  output->f5 = 0.0;
  output->f6 = 0.0;
  output->f7 = 0.0;
  output->f8 = 0.0;
  output->f9 = 0.0;
  output->f10 = 0.0;

  output->waveSin = NULL;
  output->waveCos = NULL;
  output->wave_om = 0.0;
  output->waveepoch = 0.0;
  output->nwaves = 0;

  output->ra = 0.0;
  output->dec = 0.0;
  output->pmra = 0.0;
  output->pmdec = 0.0;

  output->px = 0.;  /* parallax (mas) */
  output->dist = 0.; /* distance (kpc) */

  output->DM = 0.;  /* dispersion measure */
  output->DM1 = 0.; /* first derivative of dispersion measure */

  output->daop = 0.;
  output->daopset = 0;
  output->kin = 0.;
  output->kinset = 0;
  output->kom = 0.;
  output->komset = 0;

  /* set all errors on params to zero */
  output->raErr = 0.0;
  output->decErr = 0.0;
  output->pmraErr = 0.0;
  output->pmdecErr = 0.0;

  output->posepoch = 0.0;
  output->pepoch = 0.0;

  output->posepochErr = 0.0;
  output->pepochErr = 0.0;

  output->startTime = 0.0;
  output->finishTime = INFINITY;

  output->xpbdotErr = 0.0; /* (10^-12) */

  output->eps1Err = 0.0;      /* e*sin(w) */
  output->eps2Err = 0.0;      /* e*cos(w) */
  output->eps1dotErr = 0.0;
  output->eps2dotErr = 0.0;
  output->TascErr = 0.0;  /* time of the ascending node (used rather than T0) */

  output->wdotErr = 0.0; /* precesion of longitude of periastron w = w0 + wdot(tb-T0) (degs/year) */
  output->gammaErr = 0.0; /* gravitational redshift and time dilation parameter (s)*/
  output->PbdotErr = 0.0; /* rate of change of Pb (dimensionless 10^-12) */
  output->xdotErr = 0.0; /* rate of change of x(=asini/c) - optional (10^-12)*/
  output->edotErr = 0.0; /* rate of change of e (10^-12)*/

  output->sErr = 0.0;   /* Shapiro 'shape' parameter sin i */
  output->shapmaxErr = 0.;

  /*output->rErr=0.0;  Shapiro 'range' parameter */
  output->drErr = 0.0;
  output->dthErr = 0.0; /* (10^-6) */
  output->a0Err = 0.0;
  output->b0Err = 0.0;  /* abberation delay parameters */

  output->MErr = 0.0;   /* M = m1 + m2 (m1 = pulsar mass, m2 = companion mass) */
  output->m2Err = 0.0;  /* companion mass */

  output->f0Err = 0.0;
  output->f1Err = 0.0;
  output->f2Err = 0.0;
  output->f3Err = 0.0;
  output->f4Err = 0.0;
  output->f5Err = 0.0;
  output->f6Err = 0.0;
  output->f7Err = 0.0;
  output->f8Err = 0.0;
  output->f9Err = 0.0;
  output->f10Err = 0.0;

  output->eErr = 0.0;
  output->w0Err = 0.0;
  output->PbErr = 0.0;
  output->xErr = 0.0;
  output->T0Err = 0.0;

  output->e2Err = 0.0;
  output->w02Err = 0.0;
  output->Pb2Err = 0.0;
  output->x2Err = 0.0;
  output->T02Err = 0.0;

  output->e3Err = 0.0;
  output->w03Err = 0.0;
  output->Pb3Err = 0.0;
  output->x3Err = 0.0;
  output->T03Err = 0.0;

  output->pxErr = 0.;
  output->distErr = 0.;

  output->DMErr = 0.;
  output->DM1Err = 0.;

  output->h0 = 0.;
  output->Q22 = 0.;
  output->cosiota = 0.;
  output->iota = 0.;
  output->psi = 0.;
  output->phi0 = 0.;
  output->Aplus = 0.;
  output->Across = 0.;
  output->I21 = 0.;
  output->I31 = 0.;
  output->lambdapin = 0.;
  output->costheta = 0.;
  output->theta = 0.;
  output->C22 = 0.;
  output->C21 = 0.;
  output->phi22 = 0.;
  output->phi21 = 0.;

  output->hPlus = 0.;
  output->hCross = 0.;
  output->psiTensor = 0.;
  output->phi0Tensor = 0.;
  output->hScalarB = 0.;
  output->hScalarL = 0.;
  output->psiScalar = 0.;
  output->phi0Scalar = 0.;
  output->hVectorX = 0.;
  output->hVectorY = 0.;
  output->psiVector = 0.;
  output->phi0Vector = 0.;

  output->h0Err = 0.;
  output->Q22Err = 0.;
  output->cosiotaErr = 0.;
  output->iotaErr = 0.;
  output->psiErr = 0.;
  output->phi0Err = 0.;
  output->AplusErr = 0.;
  output->AcrossErr = 0.;
  output->I21Err = 0.;
  output->I31Err = 0.;
  output->lambdapinErr = 0.;
  output->costhetaErr = 0.;
  output->thetaErr = 0.;
  output->C22Err = 0.;
  output->C21Err = 0.;
  output->phi22Err = 0.;
  output->phi21Err = 0.;
  output->hPlusErr = 0.;
  output->hCrossErr = 0.;
  output->psiTensorErr = 0.;
  output->phi0TensorErr = 0.;
  output->hScalarBErr = 0.;
  output->hScalarLErr = 0.;
  output->psiScalarErr = 0.;
  output->phi0ScalarErr = 0.;
  output->hVectorXErr = 0.;
  output->hVectorYErr = 0.;
  output->psiVectorErr = 0.;
  output->phi0VectorErr = 0.;

  output->wave_omErr = 0.0;

  output->cgw = 1.; /* initialise the GW speed to be the speed of light */
  output->cgwErr = 0.;

  output->units = NULL;
  output->ephem = NULL;

  if ( ( fp = fopen( pulsarAndPath, "r" ) ) == NULL ) {
    XLAL_PRINT_ERROR( "Error... Cannot open .par file %s\n", pulsarAndPath );
    XLAL_ERROR_VOID( XLAL_EIO );
  }

  /* read all the pulsar data into the string array */
  while ( !feof( fp ) ) {
    /* make sure val[i] is clear first */
    sprintf( val[i], "%s", "" );

    c = fscanf( fp, "%s", val[i] );

    /* if line starts with a '#' then skip to end of line */
    if ( val[i][0] == '#' ) {
      /* skip to the end of the line */
      c = fscanf( fp, "%*[^\n]" );
      if ( feof( fp ) ) {
        break;
      }
      continue;
    }

    i++;
  }

  k = i; /* k is the end number */
  i = 0; /* reset i */

  /* set pulsar values for output */
  /* in .par files first column will param name, second will be param value,
     if third is defined it will be an integer to tell TEMPO whether to fit
     the param or not (don't need this), fourth will be the error on the
     param (in same units as the param) */

  /* convert all epochs given in MJD in .par files to secs in TDB  */
  while ( 1 ) {
    j = i;
    if ( !strcmp( val[i], "NAME" ) || !strcmp( val[i], "name" ) ) {
      output->name = XLALStringDuplicate( val[i + 1] );
      j++;
    } else if ( !strcmp( val[i], "PSRJ" ) || !strcmp( val[i], "psrj" ) ) {
      output->jname = XLALStringDuplicate( val[i + 1] );
      j++;
    } else if ( !strcmp( val[i], "PSRB" ) || !strcmp( val[i], "psrb" ) ) {
      output->bname = XLALStringDuplicate( val[i + 1] );
      j++;
    } else if ( !strcmp( val[i], "ra" ) || !strcmp( val[i], "RA" ) || !strcmp( val[i], "RAJ" ) ) {
      /* this can be in form hh:mm:ss.ss or hhmmss.ss */
      XLALTranslateHMStoRAD( &output->ra, val[i + 1] );
      j++;

      /* only try to get error if one exists */
      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        /* assuming at the moment that error is in arcsec */
        output->raErr = LAL_TWOPI * atof( val[i + 3] ) / ( 24.0 * 60.0 * 60.0 );
        j += 2;
      }
    } else if ( !strcmp( val[i], "dec" ) || !strcmp( val[i], "DEC" ) || !strcmp( val[i], "DECJ" ) ) {
      XLALTranslateDMStoRAD( &output->dec, val[i + 1] );
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        /* assuming at the moment that error is in arcsec */
        output->decErr = LAL_TWOPI * atof( val[i + 3] ) / ( 360.0 * 60.0 * 60.0 );
        j += 2;
      }
    } else if ( !strcmp( val[i], "pmra" ) || !strcmp( val[i], "PMRA" ) ) {
      /* convert pmra from mas/year to rads/sec */
      output->pmra = LAL_PI_180 * atof( val[i + 1] ) / ( 60.0 * 60.0 * 1000.*365.25 * 86400. );
      j++;
      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->pmraErr =
          LAL_PI_180 * atof( val[i + 3] ) / ( 60.0 * 60.0 * 1000.*365.25 * 86400. );
        j += 2;
      }
    } else if ( !strcmp( val[i], "pmdec" ) || !strcmp( val[i], "PMDEC" ) ) {
      /* convert pmdec from mas/year to rads/sec */
      output->pmdec = LAL_PI_180 * atof( val[j + 1] ) / ( 60.0 * 60.0 * 1000.*365.25 * 86400. );
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->pmdecErr =
          LAL_PI_180 * atof( val[i + 3] ) / ( 60.0 * 60.0 * 1000.*365.25 * 86400. );
        j += 2;
      }
    } else if ( !strcmp( val[i], "pepoch" ) || !strcmp( val[i], "PEPOCH" ) ) {
      output->pepoch = XLALTTMJDtoGPS( atof( val[i + 1] ) ); /* convert all epochs to
        from MJD to GPS seconds in TDB */
      j++;

    } else if ( !strcmp( val[i], "posepoch" ) || !strcmp( val[i], "POSEPOCH" ) ) {
      output->posepoch = XLALTTMJDtoGPS( atof( val[i + 1] ) );
      j++;
      /* position epoch in GPS seconds TDB */
    } else if ( !strcmp( val[i], "start" ) || !strcmp( val[i], "START" ) ) {
      output->startTime = XLALTTMJDtoGPS( atof( val[i + 1] ) ); /* convert all epochs to
        from MJD to GPS seconds in TDB */
      j++;

    } else if ( !strcmp( val[i], "finish" ) || !strcmp( val[i], "FINISH" ) ) {
      output->finishTime = XLALTTMJDtoGPS( atof( val[i + 1] ) );
      j++;
      /* position epoch in GPS seconds TDB */
    } else if ( !strcmp( val[i], "f0" ) || !strcmp( val[i], "F0" ) ) {
      /* in .par files exponents sometimes shown as D/d rather than e/E
         need way to check this as atof will not convert D (but will
         work for e/E (if a d/D is present atof will convert the number
         before the d/D but not the exponent */
      CHAR *loc;

      output->f0 = atof( val[i + 1] );
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        /* check if exponent contains e/E or d/D or neither */
        if ( ( loc = strstr( val[i + 3], "D" ) ) != NULL || ( loc = strstr( val[i + 3], "d" ) ) != NULL ) {
          output->f0Err = atof( val[i + 3] ) * pow( 10, atof( loc + 1 ) );
        } else {
          output->f0Err = atof( val[i + 3] );
        }
        j += 2;
      }
    } else if ( !strcmp( val[i], "f1" ) || !strcmp( val[i], "F1" ) ) {
      CHAR *loc;

      /* check if exponent contains e/E or d/D or neither */
      if ( ( loc = strstr( val[i + 1], "D" ) ) != NULL || ( loc = strstr( val[i + 1], "d" ) ) != NULL ) {
        output->f1 = atof( val[i + 1] ) * pow( 10, atof( loc + 1 ) );
      } else {
        output->f1 = atof( val[i + 1] );
      }
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        /* check if exponent contains e/E or d/D or neither */
        if ( ( loc = strstr( val[i + 3], "D" ) ) != NULL || ( loc = strstr( val[i + 3], "d" ) ) != NULL ) {
          output->f1Err = atof( val[i + 3] ) * pow( 10, atof( loc + 1 ) );
        } else {
          output->f1Err = atof( val[i + 3] );
        }
        j += 2;
      }
    } else if ( !strcmp( val[i], "f2" ) || !strcmp( val[i], "F2" ) ) {
      CHAR *loc;

      /* check if exponent contains e/E or d/D or neither */
      if ( ( loc = strstr( val[i + 1], "D" ) ) != NULL || ( loc = strstr( val[i + 1], "d" ) ) != NULL ) {
        output->f2 = atof( val[i + 1] ) * pow( 10, atof( loc + 1 ) );
      } else {
        output->f2 = atof( val[i + 1] );
      }
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        /* check if exponent contains e/E or d/D or neither */
        if ( ( loc = strstr( val[i + 3], "D" ) ) != NULL || ( loc = strstr( val[i + 3], "d" ) ) != NULL ) {
          output->f2Err = atof( val[i + 3] ) * pow( 10, atof( loc + 1 ) );
        } else {
          output->f2Err = atof( val[i + 3] );
        }
        j += 2;
      }
    } else if ( !strcmp( val[i], "f3" ) || !strcmp( val[i], "F3" ) ) {
      CHAR *loc;

      /* check if exponent contains e/E or d/D or neither */
      if ( ( loc = strstr( val[i + 1], "D" ) ) != NULL || ( loc = strstr( val[i + 1], "d" ) ) != NULL ) {
        output->f3 = atof( val[i + 1] ) * pow( 10, atof( loc + 1 ) );
      } else {
        output->f3 = atof( val[i + 1] );
      }
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        /* check if exponent contains e/E or d/D or neither */
        if ( ( loc = strstr( val[i + 3], "D" ) ) != NULL || ( loc = strstr( val[i + 3], "d" ) ) != NULL ) {
          output->f3Err = atof( val[i + 3] ) * pow( 10, atof( loc + 1 ) );
        } else {
          output->f3Err = atof( val[i + 3] );
        }
        j += 2;
      }
    } else if ( !strcmp( val[i], "f4" ) || !strcmp( val[i], "F4" ) ) {
      CHAR *loc;

      /* check if exponent contains e/E or d/D or neither */
      if ( ( loc = strstr( val[i + 1], "D" ) ) != NULL || ( loc = strstr( val[i + 1], "d" ) ) != NULL ) {
        output->f4 = atof( val[i + 1] ) * pow( 10, atof( loc + 1 ) );
      } else {
        output->f4 = atof( val[i + 1] );
      }
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        /* check if exponent contains e/E or d/D or neither */
        if ( ( loc = strstr( val[i + 3], "D" ) ) != NULL || ( loc = strstr( val[i + 3], "d" ) ) != NULL ) {
          output->f4Err = atof( val[i + 3] ) * pow( 10, atof( loc + 1 ) );
        } else {
          output->f4Err = atof( val[i + 3] );
        }
        j += 2;
      }
    } else if ( !strcmp( val[i], "f5" ) || !strcmp( val[i], "F5" ) ) {
      CHAR *loc;

      /* check if exponent contains e/E or d/D or neither */
      if ( ( loc = strstr( val[i + 1], "D" ) ) != NULL || ( loc = strstr( val[i + 1], "d" ) ) != NULL ) {
        output->f5 = atof( val[i + 1] ) * pow( 10, atof( loc + 1 ) );
      } else {
        output->f5 = atof( val[i + 1] );
      }
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        /* check if exponent contains e/E or d/D or neither */
        if ( ( loc = strstr( val[i + 3], "D" ) ) != NULL || ( loc = strstr( val[i + 3], "d" ) ) != NULL ) {
          output->f5Err = atof( val[i + 3] ) * pow( 10, atof( loc + 1 ) );
        } else {
          output->f5Err = atof( val[i + 3] );
        }
        j += 2;
      }
    } else if ( !strcmp( val[i], "f6" ) || !strcmp( val[i], "F6" ) ) {
      CHAR *loc;

      /* check if exponent contains e/E or d/D or neither */
      if ( ( loc = strstr( val[i + 1], "D" ) ) != NULL || ( loc = strstr( val[i + 1], "d" ) ) != NULL ) {
        output->f6 = atof( val[i + 1] ) * pow( 10, atof( loc + 1 ) );
      } else {
        output->f6 = atof( val[i + 1] );
      }
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        /* check if exponent contains e/E or d/D or neither */
        if ( ( loc = strstr( val[i + 3], "D" ) ) != NULL || ( loc = strstr( val[i + 3], "d" ) ) != NULL ) {
          output->f6Err = atof( val[i + 3] ) * pow( 10, atof( loc + 1 ) );
        } else {
          output->f6Err = atof( val[i + 3] );
        }
        j += 2;
      }
    } else if ( !strcmp( val[i], "f7" ) || !strcmp( val[i], "F7" ) ) {
      CHAR *loc;

      /* check if exponent contains e/E or d/D or neither */
      if ( ( loc = strstr( val[i + 1], "D" ) ) != NULL || ( loc = strstr( val[i + 1], "d" ) ) != NULL ) {
        output->f7 = atof( val[i + 1] ) * pow( 10, atof( loc + 1 ) );
      } else {
        output->f7 = atof( val[i + 1] );
      }
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        /* check if exponent contains e/E or d/D or neither */
        if ( ( loc = strstr( val[i + 3], "D" ) ) != NULL || ( loc = strstr( val[i + 3], "d" ) ) != NULL ) {
          output->f7Err = atof( val[i + 3] ) * pow( 10, atof( loc + 1 ) );
        } else {
          output->f7Err = atof( val[i + 3] );
        }
        j += 2;
      }
    } else if ( !strcmp( val[i], "f8" ) || !strcmp( val[i], "F8" ) ) {
      CHAR *loc;

      /* check if exponent contains e/E or d/D or neither */
      if ( ( loc = strstr( val[i + 1], "D" ) ) != NULL || ( loc = strstr( val[i + 1], "d" ) ) != NULL ) {
        output->f8 = atof( val[i + 1] ) * pow( 10, atof( loc + 1 ) );
      } else {
        output->f8 = atof( val[i + 1] );
      }
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        /* check if exponent contains e/E or d/D or neither */
        if ( ( loc = strstr( val[i + 3], "D" ) ) != NULL || ( loc = strstr( val[i + 3], "d" ) ) != NULL ) {
          output->f8Err = atof( val[i + 3] ) * pow( 10, atof( loc + 1 ) );
        } else {
          output->f8Err = atof( val[i + 3] );
        }
        j += 2;
      }
    } else if ( !strcmp( val[i], "f9" ) || !strcmp( val[i], "F9" ) ) {
      CHAR *loc;

      /* check if exponent contains e/E or d/D or neither */
      if ( ( loc = strstr( val[i + 1], "D" ) ) != NULL || ( loc = strstr( val[i + 1], "d" ) ) != NULL ) {
        output->f9 = atof( val[i + 1] ) * pow( 10, atof( loc + 1 ) );
      } else {
        output->f9 = atof( val[i + 1] );
      }
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        /* check if exponent contains e/E or d/D or neither */
        if ( ( loc = strstr( val[i + 3], "D" ) ) != NULL || ( loc = strstr( val[i + 3], "d" ) ) != NULL ) {
          output->f9Err = atof( val[i + 3] ) * pow( 10, atof( loc + 1 ) );
        } else {
          output->f9Err = atof( val[i + 3] );
        }
        j += 2;
      }
    } else if ( !strcmp( val[i], "f10" ) || !strcmp( val[i], "F10" ) ) {
      CHAR *loc;

      /* check if exponent contains e/E or d/D or neither */
      if ( ( loc = strstr( val[i + 1], "D" ) ) != NULL || ( loc = strstr( val[i + 1], "d" ) ) != NULL ) {
        output->f10 = atof( val[i + 1] ) * pow( 10, atof( loc + 1 ) );
      } else {
        output->f10 = atof( val[i + 1] );
      }
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        /* check if exponent contains e/E or d/D or neither */
        if ( ( loc = strstr( val[i + 3], "D" ) ) != NULL || ( loc = strstr( val[i + 3], "d" ) ) != NULL ) {
          output->f10Err = atof( val[i + 3] ) * pow( 10, atof( loc + 1 ) );
        } else {
          output->f10Err = atof( val[i + 3] );
        }
        j += 2;
      }
    } else if ( !strcmp( val[i], "WAVE_OM" ) || !strcmp( val[i], "wave_om" ) ) {
      output->wave_om = atof( val[i + 1] );
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->wave_omErr = atof( val[i + 3] );
        j += 2;
      }
    } else if ( !strcmp( val[i], "WAVEEPOCH" ) || !strcmp( val[i], "waveepoch" ) ) {
      output->waveepoch = XLALTTMJDtoGPS( atof( val[i + 1] ) );
      j++;
    } else if ( strstr( val[i], "WAVE" ) != NULL || strstr( val[i], "wave" ) != NULL ) {
      INT4 wnum = 0;

      if ( sscanf( val[i] + 4, "%d", &wnum ) != 1 ) {
        fprintf( stderr, "Error reading WAVE number from par file\n" );
        abort();
      }

      if ( wnum > output->nwaves ) {
        output->nwaves = wnum;
        output->waveSin = XLALRealloc( output->waveSin, wnum * sizeof( REAL8 ) );
        output->waveCos = XLALRealloc( output->waveCos, wnum * sizeof( REAL8 ) );
      }

      output->waveSin[wnum - 1] = atof( val[i + 1] );
      output->waveCos[wnum - 1] = atof( val[i + 2] );

      j++;
    } else if ( !strcmp( val[i], "binary" ) || !strcmp( val[i], "BINARY" ) ) {
      output->model = XLALStringDuplicate( val[i + 1] );
      j++;
    } else if ( !strcmp( val[i], "ephem" ) || !strcmp( val[i], "EPHEM" ) ) {
      output->ephem = XLALStringDuplicate( val[i + 1] );
      j++;
    } else if ( !strcmp( val[i], "units" ) || !strcmp( val[i], "UNITS" ) ) {
      output->units = XLALStringDuplicate( val[i + 1] );
      j++;
    } else if ( !strcmp( val[i], "a1" ) || !strcmp( val[i], "A1" ) ) {
      output->x = atof( val[i + 1] );
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->xErr = atof( val[i + 3] );
        j += 2;
      }
    } else if ( !strcmp( val[i], "e" ) || !strcmp( val[i], "E" ) || !strcmp( val[i], "ECC" ) ||
                !strcmp( val[i], "ecc" ) ) {
      output->e = atof( val[i + 1] );
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->eErr = atof( val[i + 3] );
        j += 2;
      }
    } else if ( !strcmp( val[i], "pb" ) || !strcmp( val[i], "PB" ) ) {
      output->Pb = atof( val[i + 1] ) * DAYSTOSECS;
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->PbErr = atof( val[i + 3] ) * DAYSTOSECS;
        j += 2;
      }
    } else if ( !strcmp( val[i], "om" ) || !strcmp( val[i], "OM" ) ) {
      output->w0 = atof( val[i + 1] ) * LAL_PI_180; /* convert radians to seconds */
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->w0Err = atof( val[i + 3] ) * LAL_PI_180;
        j += 2;
      }
    } else if ( !strcmp( val[i], "T0" ) ) {
      output->T0 = XLALTTMJDtoGPS( atof( val[i + 1] ) );
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->T0Err = atof( val[i + 3] ) * DAYSTOSECS; /* convert to seconds */
        j += 2;
      }
    } else if ( !strcmp( val[i], "Tasc" ) || !strcmp( val[i], "TASC" ) ) {
      output->Tasc = XLALTTMJDtoGPS( atof( val[i + 1] ) );
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->TascErr = atof( val[i + 3] ) * DAYSTOSECS; /* convert to seconds; */
        j += 2;
      }
    } else if ( !strcmp( val[i], "eps1" ) || !strcmp( val[i], "EPS1" ) ) {
      output->eps1 = atof( val[i + 1] );
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->eps1Err = atof( val[i + 3] );
        j += 2;
      }
    } else if ( !strcmp( val[i], "eps2" ) || !strcmp( val[i], "EPS2" ) ) {
      output->eps2 = atof( val[i + 1] );
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->eps2Err = atof( val[i + 3] );
        j += 2;
      }
    } else if ( !strcmp( val[i], "eps1dot" ) || !strcmp( val[i], "EPS1DOT" ) ) {
      output->eps1dot = atof( val[i + 1] );
      /* TEMPO2 checks if this is > 1e-7 then it's in units of 1e-12, so needs converting */
      if ( fabs( output->eps1dot ) > 1e-7 ) {
        output->eps1dot *= 1.e-12;
      }
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->eps1dotErr = atof( val[i + 3] );
        j += 2;
      }
    } else if ( !strcmp( val[i], "eps2dot" ) || !strcmp( val[i], "EPS2DOT" ) ) {
      output->eps2dot = atof( val[i + 1] );
      /* TEMPO2 checks if this is > 1e-7 then it's in units of 1e-12, so needs converting */
      if ( fabs( output->eps2dot ) > 1e-7 ) {
        output->eps2dot *= 1.e-12;
      }
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->eps2dotErr = atof( val[i + 3] );
        j += 2;
      }
    } else if ( !strcmp( val[i], "xpbdot" ) || !strcmp( val[i], "XPBDOT" ) ) {
      output->xpbdot = atof( val[i + 1] );
      /* TEMPO2 checks if this is > 1e-7 then it's in units of 1e-12, so needs converting */
      if ( fabs( output->xpbdot ) > 1e-7 ) {
        output->xpbdot *= 1.e-12;
      }
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->xpbdotErr = atof( val[i + 3] );
        j += 2;
      }
    } else if ( !strcmp( val[i], "omdot" ) || !strcmp( val[i], "OMDOT" ) ) {
      output->wdot = atof( val[i + 1] ) * LAL_PI_180 / ( 365.25 * DAYSTOSECS ); /* convert degs/years to rads/sec */
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->wdotErr = atof( val[i + 3] ) * LAL_PI_180 / ( 365.25 * DAYSTOSECS );
        j += 2;
      }
    } else if ( !strcmp( val[i], "pbdot" ) || !strcmp( val[i], "PBDOT" ) ) {
      output->Pbdot = atof( val[i + 1] );
      /* TEMPO2 checks if this is > 1e-7 then it's in units of 1e-12, so needs converting */
      if ( fabs( output->Pbdot ) > 1e-7 ) {
        output->Pbdot *= 1.e-12;
      }
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->PbdotErr = atof( val[i + 3] );
        j += 2;
      }
    } else if ( !strcmp( val[i], "xdot" ) || !strcmp( val[i], "XDOT" ) ) {
      output->xdot = atof( val[i + 1] );
      /* TEMPO2 checks if this is > 1e-7 then it's in units of 1e-12, so needs converting */
      if ( fabs( output->xdot ) > 1e-7 ) {
        output->xdot *= 1.e-12;
      }
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->xdotErr = atof( val[i + 3] );
        j += 2;
      }
    } else if ( !strcmp( val[i], "edot" ) || !strcmp( val[i], "EDOT" ) ) {
      output->edot = atof( val[i + 1] );
      /* TEMPO2 checks if this is > 1e-7 then it's in units of 1e-12, so needs converting */
      if ( fabs( output->edot ) > 1e-7 ) {
        output->edot *= 1.e-12;
      }
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->edotErr = atof( val[i + 3] );
        j += 2;
      }

      /* some of the parameter files in the ATNF catalogue have values
         of EDOT that are stupidly large e.g. O(1e33). These can cause
         the time delay routines to fail, so if values of EDOT are
         greater than 10000 ignore them and set it to zero */
      if ( output->edot > 10000 ) {
        output->edot = 0.;
        output->edotErr = 0.;
      }
    } else if ( !strcmp( val[i], "gamma" ) || !strcmp( val[i], "GAMMA" ) ) {
      output->gamma = atof( val[i + 1] );
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->gammaErr = atof( val[i + 3] );
        j += 2;
      }
    } else if ( !strcmp( val[i], "sini" ) || !strcmp( val[i], "SINI" ) ) {
      output->sstr = XLALStringDuplicate( val[i + 1] );
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->sErr = atof( val[i + 3] );
        j += 2;
      }
    } else if ( !strcmp( val[i], "mtot" ) || !strcmp( val[i], "MTOT" ) ) {
      output->M = atof( val[i + 1] ) * LALPULSAR_TEMPO2_MSUN_SI;
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->MErr = atof( val[i + 3] ) * LALPULSAR_TEMPO2_MSUN_SI;
        j += 2;
      }
    } else if ( !strcmp( val[i], "m2" ) || !strcmp( val[i], "M2" ) ) {
      output->m2 = atof( val[i + 1] ) * LALPULSAR_TEMPO2_MSUN_SI;
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->m2Err = atof( val[i + 3] ) * LALPULSAR_TEMPO2_MSUN_SI;
        j += 2;
      }
    } else if ( !strcmp( val[i], "a0" ) || !strcmp( val[i], "A0" ) ) {
      output->a0 = atof( val[i + 1] );
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->a0Err = atof( val[i + 3] );
        j += 2;
      }
    } else if ( !strcmp( val[i], "b0" ) || !strcmp( val[i], "B0" ) ) {
      output->b0 = atof( val[i + 1] );
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->b0Err = atof( val[i + 3] );
        j += 2;
      }
    } else if ( !strcmp( val[i], "dr" ) || !strcmp( val[i], "DR" ) ) {
      output->dr = atof( val[i + 1] );
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->drErr = atof( val[i + 3] );
        j += 2;
      }
    } else if ( !strcmp( val[i], "dtheta" ) || !strcmp( val[i], "DTHETA" ) ) {
      output->dth = atof( val[i + 1] );
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->dthErr = atof( val[i + 3] );
        j += 2;
      }
    } else if ( !strcmp( val[i], "shapmax" ) || !strcmp( val[i], "SHAPMAX" ) ) {
      output->shapmax = atof( val[i + 1] );
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->shapmaxErr = atof( val[i + 3] );
        j += 2;
      }
    }

    /* parameters for Kopeikin terms */
    else if ( !strcmp( val[i], "D_AOP" ) || !strcmp( val[i], "d_aop" ) ) {
      /* convert into 1/rads (factor from T2model.C in TEMPO2 */
      output->daop = atof( val[i + 1] ) * 3600.0 / LAL_PI_180;
      output->daopset = 1;
      j++;
    } else if ( !strcmp( val[i], "KIN" ) || !strcmp( val[i], "kin" ) ) {
      output->kin = atof( val[i + 1] ) * LAL_PI_180; /* convert degs to rads */
      output->kinset = 1;
      j++;
    } else if ( !strcmp( val[i], "KOM" ) || !strcmp( val[i], "kom" ) ) {
      output->kom = atof( val[i + 1] ) * LAL_PI_180; /* convert degs to rads */
      output->komset = 1;
      j++;
    }

    /* parameters for distance */
    else if ( !strcmp( val[i], "px" ) || !strcmp( val[i], "PX" ) ) { /* parallax */
      /* convert from mas to rads (factor from T2model.C in TEMPO2) */
      output->px = atof( val[i + 1] ) * LAL_PI_180 / 3600.0e3;
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->pxErr = atof( val[i + 3] );
        j += 2;
      }
    } else if ( !strcmp( val[i], "dist" ) || !strcmp( val[i], "DIST" ) ) { /* distance */
      output->dist = atof( val[i + 1] ); /* in kpc */
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->distErr = atof( val[i + 3] );
        j += 2;
      }
    }

    /* dispersion measure parameters */
    else if ( !strcmp( val[i], "dm" ) || !strcmp( val[i], "DM" ) ) {
      output->DM = atof( val[i + 1] );
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->DMErr = atof( val[i + 3] );
        j += 2;
      }
    } else if ( !strcmp( val[i], "dm1" ) || !strcmp( val[i], "DM1" ) ) {
      output->DM1 = atof( val[i + 1] );
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->DM1Err = atof( val[i + 3] );
        j += 2;
      }
    }

    /* add parameters extra orbital parameters for the BT1P and BT2P models */
    else if ( !strcmp( val[i], "a1_2" ) || !strcmp( val[i], "A1_2" ) ) {
      output->x2 = atof( val[i + 1] );
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->x2Err = atof( val[i + 3] );
        j += 2;
      }
    } else if ( !strcmp( val[i], "e_2" ) || !strcmp( val[i], "E_2" ) ||
                !strcmp( val[i], "ECC_2" )      || !strcmp( val[i], "ecc_2" ) ) {
      output->e2 = atof( val[i + 1] );
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->e2Err = atof( val[i + 3] );
        j += 2;
      }
    } else if ( !strcmp( val[i], "pb_2" ) || !strcmp( val[i], "PB_2" ) ) {
      output->Pb2 = atof( val[i + 1] ) * DAYSTOSECS;
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->Pb2Err = atof( val[i + 3] ) * DAYSTOSECS;
        j += 2;
      }
    } else if ( !strcmp( val[i], "om_2" ) || !strcmp( val[i], "OM_2" ) ) {
      output->w02 = atof( val[i + 1] ) * LAL_PI_180;
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->w02Err = atof( val[i + 3] ) * LAL_PI_180;
        j += 2;
      }
    } else if ( !strcmp( val[i], "T0_2" ) ) {
      output->T02 = XLALTTMJDtoGPS( atof( val[i + 1] ) );
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->T02Err = atof( val[i + 3] ) * DAYSTOSECS; /* convert to seconds */
        j += 2;
      }
    } else if ( !strcmp( val[i], "a1_3" ) || !strcmp( val[i], "A1_3" ) ) {
      output->x3 = atof( val[i + 1] );
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->x3Err = atof( val[i + 3] );
        j += 2;
      }
    } else if ( !strcmp( val[i], "e_3" ) || !strcmp( val[i], "E_3" ) ||
                !strcmp( val[i], "ECC_3" )      || !strcmp( val[i], "ecc_3" ) ) {
      output->e3 = atof( val[i + 1] );
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->e3Err = atof( val[i + 3] );
        j += 2;
      }
    } else if ( !strcmp( val[i], "pb_3" ) || !strcmp( val[i], "PB_3" ) ) {
      output->Pb3 = atof( val[i + 1] ) * DAYSTOSECS;
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->Pb3Err = atof( val[i + 3] ) * DAYSTOSECS;
        j += 2;
      }
    } else if ( !strcmp( val[i], "om_3" ) || !strcmp( val[i], "OM_3" ) ) {
      output->w03 = atof( val[i + 1] ) * LAL_PI_180;
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->w03Err = atof( val[i + 3] ) * LAL_PI_180;
        j += 2;
      }
    } else if ( !strcmp( val[i], "T0_3" ) ) {
      output->T03 = XLALTTMJDtoGPS( atof( val[i + 1] ) );
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->T03Err = atof( val[i + 3] ) * DAYSTOSECS; /* convert to seconds */
        j += 2;
      }
    }
    /* orbital frequency coefficients for BTX model (up to 12 FB coefficients), but
       only one orbit at the moment i.e. only a two body system */
    else if ( val[i][0] == 'F' && val[i][1] == 'B' ) {
      INT4 fbnum = 0;
      CHAR *loc;

      if ( strlen( val[i] ) == 2 ) {
        fbnum = 0;  /* only one coefficient */
      } else {
        if ( sscanf( val[i] + 2, "%d", &fbnum ) != 1 ) {
          fprintf( stderr, "Error reading FB value from par file\n" );
          abort();
        }
      }

      /* add to number of coefficients */
      if ( output->nfb < fbnum + 1 ) {
        output->fb = XLALRealloc( output->fb, ( fbnum + 1 ) * sizeof( REAL8 ) );
        output->fbErr = XLALRealloc( output->fbErr, ( fbnum + 1 ) * sizeof( REAL8 ) );
        output->nfb = fbnum + 1;
      }

      /* check if exponent contains e/E or d/D or neither */
      if ( ( loc = strstr( val[i + 1], "D" ) ) != NULL || ( loc = strstr( val[i + 1], "d" ) ) != NULL ) {
        output->fb[fbnum] = atof( val[i + 1] ) * pow( 10, atof( loc + 1 ) );
      } else {
        output->fb[fbnum] = atof( val[i + 1] );
      }
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        /* check if exponent contains e/E or d/D or neither */
        if ( ( loc = strstr( val[i + 3], "D" ) ) != NULL || ( loc = strstr( val[i + 3], "d" ) ) != NULL ) {
          output->fbErr[fbnum] = atof( val[i + 3] ) * pow( 10, atof( loc + 1 ) );
        } else {
          output->fbErr[fbnum] = atof( val[i + 3] );
        }
        j += 2;
      }
    }
    /* read in pulsar gravitational wave parameters */
    else if ( !strcmp( val[i], "h0" ) || !strcmp( val[i], "H0" ) ) {
      output->h0 = atof( val[i + 1] );
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->h0Err = atof( val[i + 3] );
        j += 2;
      }
    } else if ( !strcmp( val[i], "q22" ) || !strcmp( val[i], "Q22" ) ) {
      output->Q22 = atof( val[i + 1] );
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->Q22Err = atof( val[i + 3] );
        j += 2;
      }
    } else if ( !strcmp( val[i], "cosiota" ) || !strcmp( val[i], "COSIOTA" ) ||
                !strcmp( val[i], "ciota" ) || !strcmp( val[i], "CIOTA" ) ) {
      output->cosiota = atof( val[i + 1] );
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->cosiotaErr = atof( val[i + 3] );
        j += 2;
      }
    } else if ( !strcmp( val[i], "iota" ) || !strcmp( val[i], "IOTA" ) ) {
      output->iota = atof( val[i + 1] );
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->iotaErr = atof( val[i + 3] );
        j += 2;
      }
    } else if ( !strcmp( val[i], "psi" ) || !strcmp( val[i], "PSI" ) ) {
      output->psi = atof( val[i + 1] );
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->psiErr = atof( val[i + 3] );
        j += 2;
      }
    } else if ( !strcmp( val[i], "phi0" ) || !strcmp( val[i], "PHI0" ) ) {
      output->phi0 = atof( val[i + 1] );
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->phi0Err = atof( val[i + 3] );
        j += 2;
      }
    } else if ( !strcmp( val[i], "aplus" ) || !strcmp( val[i], "APLUS" ) ) {
      output->Aplus = atof( val[i + 1] );
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->AplusErr = atof( val[i + 3] );
        j += 2;
      }
    } else if ( !strcmp( val[i], "across" ) || !strcmp( val[i], "ACROSS" ) ) {
      output->Across = atof( val[i + 1] );
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->AcrossErr = atof( val[i + 3] );
        j += 2;
      }
    } else if ( !strcmp( val[i], "i21" ) || !strcmp( val[i], "I21" ) ) {
      output->I21 = atof( val[i + 1] );
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->I21Err = atof( val[i + 3] );
        j += 2;
      }
    } else if ( !strcmp( val[i], "i31" ) || !strcmp( val[i], "I31" ) ) {
      output->I31 = atof( val[i + 1] );
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->I31Err = atof( val[i + 3] );
        j += 2;
      }
    } else if ( !strcmp( val[i], "lambdapin" ) || !strcmp( val[i], "LAMBDAPIN" ) ) {
      output->lambdapin = atof( val[i + 1] );
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->lambdapinErr = atof( val[i + 3] );
        j += 2;
      }
    } else if ( !strcmp( val[i], "costheta" ) || !strcmp( val[i], "COSTHETA" ) ) {
      output->costheta = atof( val[i + 1] );
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->costhetaErr = atof( val[i + 3] );
        j += 2;
      }
    } else if ( !strcmp( val[i], "theta" ) || !strcmp( val[i], "THETA" ) ) {
      output->theta = atof( val[i + 1] );
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->thetaErr = atof( val[i + 3] );
        j += 2;
      }
    } else if ( !strcmp( val[i], "c22" ) || !strcmp( val[i], "C22" ) ) {
      output->C22 = atof( val[i + 1] );
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->C22Err = atof( val[i + 3] );
        j += 2;
      }
    } else if ( !strcmp( val[i], "c21" ) || !strcmp( val[i], "C21" ) ) {
      output->C21 = atof( val[i + 1] );
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->C21Err = atof( val[i + 3] );
        j += 2;
      }
    } else if ( !strcmp( val[i], "phi22" ) || !strcmp( val[i], "PHI22" ) ) {
      output->phi22 = atof( val[i + 1] );
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->phi22Err = atof( val[i + 3] );
        j += 2;
      }
    } else if ( !strcmp( val[i], "phi21" ) || !strcmp( val[i], "phi21" ) ) {
      output->phi21 = atof( val[i + 1] );
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->phi21Err = atof( val[i + 3] );
        j += 2;
      }
    } else if ( !strcmp( val[i], "hplus" ) || !strcmp( val[i], "HPLUS" ) ) {
      output->hPlus = atof( val[i + 1] );
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->hPlus = atof( val[i + 3] );
        j += 2;
      }
    } else if ( !strcmp( val[i], "hcross" ) || !strcmp( val[i], "HCROSS" ) ) {
      output->hCross = atof( val[i + 1] );
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->hCrossErr = atof( val[i + 3] );
        j += 2;
      }
    } else if ( !strcmp( val[i], "psitensor" ) || !strcmp( val[i], "PSITENSOR" ) ) {
      output->psiTensor = atof( val[i + 1] );
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->psiTensorErr = atof( val[i + 3] );
        j += 2;
      }
    } else if ( !strcmp( val[i], "phi0tensor" ) || !strcmp( val[i], "PHI0TENSOR" ) ) {
      output->phi0Tensor = atof( val[i + 1] );
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->phi0TensorErr = atof( val[i + 3] );
        j += 2;
      }
    } else if ( !strcmp( val[i], "hscalarb" ) || !strcmp( val[i], "HSCALARB" ) ) {
      output->hScalarB = atof( val[i + 1] );
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->hScalarBErr = atof( val[i + 3] );
        j += 2;
      }
    } else if ( !strcmp( val[i], "hscalarl" ) || !strcmp( val[i], "HSCALARL" ) ) {
      output->hScalarL = atof( val[i + 1] );
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->hScalarLErr = atof( val[i + 3] );
        j += 2;
      }
    } else if ( !strcmp( val[i], "psiscalar" ) || !strcmp( val[i], "PSISCALAR" ) ) {
      output->psiScalar = atof( val[i + 1] );
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->psiScalarErr = atof( val[i + 3] );
        j += 2;
      }
    } else if ( !strcmp( val[i], "phi0scalar" ) || !strcmp( val[i], "PHI0SCALAR" ) ) {
      output->phi0Scalar = atof( val[i + 1] );
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->phi0ScalarErr = atof( val[i + 3] );
        j += 2;
      }
    } else if ( !strcmp( val[i], "hvectorx" ) || !strcmp( val[i], "HVECTORX" ) ) {
      output->hVectorX = atof( val[i + 1] );
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->hVectorXErr = atof( val[i + 3] );
        j += 2;
      }
    } else if ( !strcmp( val[i], "hvectory" ) || !strcmp( val[i], "HVECTORY" ) ) {
      output->hVectorY = atof( val[i + 1] );
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->hVectorYErr = atof( val[i + 3] );
        j += 2;
      }
    } else if ( !strcmp( val[i], "psivector" ) || !strcmp( val[i], "PSIVECTOR" ) ) {
      output->psiVector = atof( val[i + 1] );
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->psiVectorErr = atof( val[i + 3] );
        j += 2;
      }
    } else if ( !strcmp( val[i], "phi0vector" ) || !strcmp( val[i], "PHI0VECTOR" ) ) {
      output->phi0Vector = atof( val[i + 1] );
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->phi0VectorErr = atof( val[i + 3] );
        j += 2;
      }
    } else if ( !strcmp( val[i], "cgw" ) || !strcmp( val[i], "CGW" ) ) {
      output->cgw = atof( val[i + 1] );
      j++;

      if ( atoi( val[i + 2] ) == 1 && i + 2 < k ) {
        output->cgwErr = atof( val[i + 3] );
        j += 2;
      }
    }

    if ( j == i ) {
      i++;
    } else {
      i += ( j - i );
    }

    if ( i >= k ) {
      break;
    }
  }

  /*fprintf(stderr, "Have I got to the end of LALReadPARFile.\n");*/
  fclose( fp );

  /* check linked parameters */
  if ( output->sstr != NULL ) {
    if ( !strcmp( output->sstr, "KIN" ) || !strcmp( output->sstr, "kin" ) ) {
      if ( output->kinset ) {
        output->s = sin( output->kin );
      } else {
        XLAL_PRINT_ERROR( "Error... KIN not set in .par file %s\n", pulsarAndPath );
        XLAL_ERROR_VOID( XLAL_EIO );
      }
    } else {
      output->s = atof( output->sstr );
    }
  }
}


/* function to print out to screen all the pulsar parameters and there associated errors */
void PrintPulsarParameters( BinaryPulsarParams params )
{
  fprintf( stderr, "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n" );
  fprintf( stderr, "PULSAR %s :\n", params.name );
  fprintf( stderr, "sky position:\tra %.7lf +/- %.3le rads, dec %.7lf +/- %.3le rads\n", params.ra,
           params.raErr, params.dec, params.decErr );
  if ( params.pmra != 0. || params.pmdec != 0. )
    fprintf( stderr, "proper motion:\tra %.4le +/- %.1le rads/s, dec %.4le +/- %.1le rads/s\n",
             params.pmra, params.pmraErr, params.pmdec, params.pmdecErr );
  if ( params.pepoch != 0. || params.posepoch != 0. )
    fprintf( stderr, "epochs:\tperiod %lf (GPS), position %lf (GPS)\n", params.pepoch,
             params.posepoch );
  fprintf( stderr, "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n\n" );

  fprintf( stderr, "Frequency parameters\n" );
  if ( params.f0 != 0. ) {
    fprintf( stderr, "\tf0 = %.10lf +/- %.3le (Hz)\n", params.f0, params.f0Err );
  }
  if ( params.f1 != 0. ) {
    fprintf( stderr, "\tf1 = %.5le +/- %.3le (Hz/s)\n", params.f1, params.f1Err );
  }
  if ( params.f2 != 0. ) {
    fprintf( stderr, "\tf1 = %.5le +/- %.3le (Hz/s^2)\n", params.f2, params.f2Err );
  }
  /* print binary parameters */
  if ( params.model != NULL ) {
    fprintf( stderr, "\nBinary parameters:\tmodel %s\n", params.model );

    fprintf( stderr, "Keplarian parameters:-\n" );
    if ( params.Pb != 0. ) {
      fprintf( stderr, "\tperiod = %lf +/- %.3le (days)\n", params.Pb, params.PbErr );
    }
    if ( params.x != 0. )
      fprintf( stderr, "\tprojected semi-major axis = %lf +/- %.3le (light sec)\n", params.x,
               params.xErr );
    if ( params.e != 0. ) {
      fprintf( stderr, "\teccentricity = %lf +/- %.3le\n", params.e, params.eErr );
    }
    if ( params.w0 != 0. )
      fprintf( stderr, "\tlongitude of periastron = %lf +/- %.3lf (degs)\n", params.w0,
               params.w0Err );
    if ( params.T0 != 0. ) {
      fprintf( stderr, "\ttime of periastron = %lf +/- %.3lf (GPS)\n", params.T0, params.T0Err );
    }
    if ( params.Tasc != 0. )
      fprintf( stderr, "\ttime of ascending node = %lf +/- %.3lf (GPS)\n", params.Tasc,
               params.TascErr );
    if ( params.eps1 != 0. )
      fprintf( stderr, "\tfirst Laplace-Lagrange parameter (eps1) = %le +/- %.3le\n", params.eps1,
               params.eps1Err );
    if ( params.eps2 != 0. )
      fprintf( stderr, "\tsecond Laplace-Lagrange parameter (eps1) = %le +/- %.3le\n", params.eps2,
               params.eps2Err );
    if ( params.eps2 != 0. )
      fprintf( stderr, "\tsecond Laplace-Lagrange parameter (eps1) = %le +/- %.3le\n", params.eps2,
               params.eps2Err );

    /*fprintf(stderr, "Post-Newtonian parameters:-\n");
    if(params.gamma != 0.)
      fprintf(stderr, "\tGravitational redshift parameter = %le +/- %.3le\n", params.gamma,
    params.gammaErr);*/

  }
}


LALStringVector *XLALReadTEMPOCorFile( REAL8Array *cormat, CHAR *corfile )
{
  FILE *fp = NULL;
  CHAR *firstline = XLALStringDuplicate( "" );
  CHAR onechar[2];
  INT4 i = 0, numPars = 0, c = 1, sl = 0;
  LALStringVector *tmpparams = NULL; /* temporary parameter names */
  LALStringVector *params = NULL;
  UINT4Vector *dims = NULL;

  /* check the file exists */
  if ( access( corfile, F_OK ) != 0 ) {
    XLAL_PRINT_ERROR( "Error... correlation matrix file does not exist!\n" );
    XLAL_ERROR_NULL( XLAL_EIO );
  }

  /* open file */
  if ( ( fp = fopen( corfile, "r" ) ) == NULL ) {
    XLAL_PRINT_ERROR( "Error... cannot open correlation matrix file!\n" );
    XLAL_ERROR_NULL( XLAL_EIO );
  }

  /* read in first line of the file */
  while ( !strchr( fgets( onechar, 2, fp ), '\n' ) ) {
    firstline = XLALStringAppend( firstline, onechar );
  }

  sl = strlen( firstline );

  /* count the number of parameters */
  for ( i = 0; i < sl; i++ ) {
    /* use isspace as delimiters could be unknown generic whitespace */
    if ( !isspace( firstline[i] ) ) {
      if ( c ) {
        numPars++;
        c = 0;
      }
    } else {
      c = 1;
    }
  }

  /* parse the line and put into the params vector */
  rewind( fp ); /* rewind to start of the file */
  for ( i = 0; i < numPars; i++ ) {
    CHAR tmpStr[128];

    if ( fscanf( fp, "%s", tmpStr ) == EOF ) {
      XLAL_PRINT_ERROR( "Error... Problem reading first line of correlation\
 matrix!\n" );
      XLAL_ERROR_NULL( XLAL_EIO );
    }

    tmpparams = XLALAppendString2Vector( tmpparams, tmpStr );

    /* convert some parameter names to a more common convention */
    if ( !XLALStringCaseCompare( tmpStr, "RAJ" ) ) { /* convert RAJ to ra */
      params = XLALAppendString2Vector( params, "RA" );
    } else if ( !XLALStringCaseCompare( tmpStr, "DECJ" ) ) { /* convert DECJ to dec */
      params = XLALAppendString2Vector( params, "DEC" );
    } else {
      params = XLALAppendString2Vector( params, tmpStr );
    }
  }

  dims = XLALCreateUINT4Vector( 2 );
  dims->data[0] = numPars;
  dims->data[1] = numPars;

  /* set the correlation matrix to the correct size */
  cormat = XLALResizeREAL8Array( cormat, dims );

  /* read through covariance values */
  for ( i = 0; i < numPars; i++ ) {
    CHAR tmpStr[128];
    INT4 j = 0;

    if ( fscanf( fp, "%s", tmpStr ) == EOF ) {
      XLAL_PRINT_ERROR( "Error... problem reading in correlation matrix!\n" );
      XLAL_ERROR_NULL( XLAL_EIO );
    }

    if ( strcmp( tmpStr, tmpparams->data[i] ) ) {
      XLAL_PRINT_ERROR( "Error... problem reading in correlation matrix. \
Parameters not in consistent order!\n" );
      XLAL_ERROR_NULL( XLAL_EIO );
    }

    for ( j = 0; j < i + 1; j++ ) {
      REAL8 tmpval = 0.;

      if ( fscanf( fp, "%lf", &tmpval ) == EOF ) {
        XLAL_PRINT_ERROR( "Error... problem reading in correlation matrix!\n" );
        XLAL_ERROR_NULL( XLAL_EIO );
      }

      /* if off diagonal values are +/-1 set to +/- 0.99999 */
      if ( j != i && fabs( tmpval ) == 1. ) {
        tmpval *= 0.99999;
      }

      cormat->data[i * numPars + j] = tmpval;

      /* set opposite elements */
      if ( j != i ) {
        cormat->data[j * numPars + i] = tmpval;
      }
    }
  }

  XLALDestroyStringVector( tmpparams );

  return params;
}


/* functions for converting times given in Terrestrial time TT or TDB in MJD to
times in GPS - this is important for epochs given in .par files which are in
TDB. TT and GPS are different by a factor of 51.184 secs, this is just the
historical factor of 32.184 secs between TT and TAI (International Atomic Time)
and the other 19 seconds come from the leap seonds added between the TAI and
UTC up to the point of definition of GPS time at UTC 01/01/1980 (see
http://www.stjarnhimlen.se/comp/time.html for details) */

/* a very good paper describing the tranforms between different time systems
and why they are necessary can be found in Seidelmann and Fukushima, A&A 265
(1992) http://ukads.nottingham.ac.uk/abs/1992A%26A...265..833S */

/** This function converts a MJD format time corrected to Terrestrial Time (TT)
 * into an equivalent GPS time */
REAL8 XLALTTMJDtoGPS( REAL8 MJD )
{

  LIGOTimeGPS GPS;
  REAL8 mjdInt, mjdFrac;
  mjdFrac = modf( MJD, &mjdInt );
  XLAL_CHECK_REAL8( XLALTranslateMJDTTtoGPS( &GPS, ( INT4 )mjdInt, mjdFrac ) != NULL, XLAL_EFUNC );

  return XLALGPSGetREAL8( &GPS );
}

/** If you have an MJD arrival time on the Earth then this will convert it to
 * the equivalent GPS time in TDB (see Table 1 of Seidelmann and Fukushima,
 * Astronomy & Astrophysics, 265, 833-838 (1992).
 *
 * Note that XLALBarycenter performs these TDBtoTT corrections (i.e. the
 * Einstein delay) when correcting a GPS time on the Earth to TDB. Also, for
 * TEMPO produced pulsar epochs given in MJD these are already in the TDB
 * system and an equivalent GPS time in the TDB can be calculated just using
 * \c XLALTTMJDtoGPS.
 */
REAL8 XLALTDBMJDtoGPS( REAL8 MJD )
{
  REAL8 GPS;
  REAL8 T, TDBtoTT;

  /* Check not before the start of GPS time */
  XLAL_CHECK_REAL8( MJD >= GPS0MJD, XLAL_EDOM, "Input MJD time %.1f is not in range, must be > %.1f.\n", MJD, GPS0MJD );

  /* use factors from Table 1 of Seidelmann and Fukushima, Astronomy &
   * Astrophysics, 265, 833-838 (1992) where TDB = TDT + P
   * and:
   * P = 0.0016568 sin(35999.37 degs x T + 357.5 degs) +
         0.0000224 sin(32964.5 degs x T + 246.0 degs) +
         0.0000138 sin(71998.7 degs x T + 355.0 degs) +
         0.0000048 sin(3034.9 degs x T + 25.0 degs) +
         0.0000047 sin(34777.3 degs x T + 230.0 degs)
   * and T is the elapsed time from J2000 (which has a Julian day date of
   * JD 2451545.0) in Julian centuries.*/
  T = MJD + ( XLAL_MJD_REF - XLAL_EPOCH_J2000_0_JD );
  T /= 36525.; /* covert days to Julian centuries */

  /* time diff in seconds (the Einstein delay) */
  TDBtoTT = 0.0016568 * sin( ( 35999.37 * T + 357.5 ) * LAL_PI_180 ) +
            0.0000224 * sin( ( 32964.5 * T +  246.0 ) * LAL_PI_180 ) +
            0.0000138 * sin( ( 71998.7 * T +  355.0 ) * LAL_PI_180 ) +
            0.0000048 * sin( ( 3034.9 * T + 25.0 ) * LAL_PI_180 ) +
            0.0000047 * sin( ( 34777.3 * T + 230.0 ) * LAL_PI_180 );

  /* convert TDB to TT (TDB-TDBtoTT) and then convert TT to GPS */
  /* there is the magical number factor of 32.184 + 19 leap seconds to the
   * start of GPS time */
  GPS = ( MJD - GPS0MJD ) * 86400. - GPS_TDT - TDBtoTT;

  return GPS;
}


/** If you have an MJD arrival time on the Earth then this will convert it to
 * the equivalent GPS time in TCB (see Table 1 of Seidelmann and Fukushima,
 * Astronomy & Astrophysics, 265, 833-838, 1992).
 *
 * Note that for default TEMPO2 produced pulsar epochs given in MJD these are
 * already in the TCB system and an equivalent GPS time in the TCB can be
 * calculated just using \c XLALTTMJDtoGPS. */
REAL8 XLALTCBMJDtoGPS( REAL8 MJD )
{
  REAL8 GPS;
  REAL8 Tdiff;
  REAL8 TCBtoTDB;

  /* Check not before the start of GPS time (MJD 44244) */
  XLAL_CHECK_REAL8( MJD >= GPS0MJD, XLAL_EDOM, "Input MJD time %.1f is not in\
 range, must be > %.1f.\n", MJD, GPS0MJD );

  /* from Seidelmann and Fukushima we have a linear drift term:
   * TCB - TDB = 1.550506e-8 x (JD - 2443144.5) x 86400
   */
  Tdiff = ( MJD + XLAL_MJD_REF - 2443144.5 ) * 86400.;
  TCBtoTDB = 1.550506e-8 * Tdiff;

  /* convert from TDB to GPS */
  GPS = XLALTDBMJDtoGPS( MJD );

  /* add extra factor as the MJD was really in TCB not TDB) */
  GPS -= TCBtoTDB;

  return GPS;
}
