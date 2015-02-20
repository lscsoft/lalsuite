// Copyright (C) 2015 Reinhard Prix
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with with program; see the file COPYING. If not, write to the
// Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
// MA  02111-1307  USA
//

// ---------- INCLUDES ----------
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>

#include <config.h>

#include <lal/LALConstants.h>

#define IN_VECTORMATH
#include <lal/VectorMath.h>

// ---------- local DEFINES ----------
#if defined(HAVE_SSE)
#define VM_HAVE_SSE 1
#else
#define VM_HAVE_SSE 0
#endif
#if defined(HAVE_AVX)
#define VM_HAVE_AVX 1
#else
#define VM_HAVE_AVX 0
#endif

// ----- Macros -----
// ---------- internal types ----------

// ---------- Global variables ----------
// default VectoDevice is the best available one:
#if defined(HAVE_AVX)
static VectorDevice_type currentVectorDevice = VECTORDEVICE_AVX;
#elif defined(HAVE_SSE)
static VectorDevice_type currentVectorDevice = VECTORDEVICE_SSE;
#else
static VectorDevice_type currentVectorDevice = VECTORDEVICE_FPU;
#endif

static const struct {
  const char *const name;
  BOOLEAN available;
} allVectorDevices[VECTORDEVICE_END] = {
  [VECTORDEVICE_FPU]	= {"FPU",	1 },
  [VECTORDEVICE_SSE] 	= {"SSE", 	VM_HAVE_SSE },
  [VECTORDEVICE_AVX] 	= {"AVX", 	VM_HAVE_AVX }
} ;

//==================== FUNCTION DEFINITIONS ====================*/
///
/// Set Vector-math device at runtime, return error if invalid or unavailable
///
int
XLALVectorDeviceSet ( VectorDevice_type device )
{
  XLAL_CHECK ( (device > VECTORDEVICE_START) && (device < VECTORDEVICE_END), XLAL_EDOM );
  XLAL_CHECK ( allVectorDevices[device].available, XLAL_EINVAL, "Sorry, vector device '%s' not available in this build\n", allVectorDevices[device].name );

  currentVectorDevice = device;
  return XLAL_SUCCESS;
} // XLALVectorDeviceSet()

///
/// Return internal number of currently-selected vector device
///
VectorDevice_type
XLALVectorDeviceGet ( void )
{
  return currentVectorDevice;
} // XLALVectorDeviceGet()

///
/// Return 1 or 0 depending on whether the given vector device is available in this code
///
int
XLALVectorDeviceIsAvailable ( VectorDevice_type device )
{
  XLAL_CHECK ( (device > VECTORDEVICE_START) && (device < VECTORDEVICE_END), XLAL_EDOM );
  return allVectorDevices[device].available;
} // XLALVectorDeviceIsAvailable()

///
/// Provide human-readable names for the different vector device ids in #VectorDevice_type
///
const CHAR *
XLALVectorDeviceName ( VectorDevice_type device )
{
  return allVectorDevices[device].name;
}// XLALVectorDeviceName()

///
/// Return pointer to a static help string enumerating all (available) #VectorDevice_type options.
/// Also indicates which is the currently-selected one.
///
const CHAR *
XLALVectorDeviceHelpString ( void )
{
  static int firstCall = 1;
  static CHAR helpstr[1024];
  if ( firstCall )
    {
      CHAR buf[1024];
      strncpy (helpstr, "Available vector devices: (", sizeof(helpstr));
      UINT4 len = strlen(helpstr);
      const CHAR *separator = "";
      for (UINT4 i = VECTORDEVICE_START + 1; i < VECTORDEVICE_END; i++ )
        {
          if ( ! allVectorDevices[i].available ) {
            continue;
          }
          snprintf ( buf, sizeof(buf), "%s%s", separator, allVectorDevices[i].name );
          separator="|";
          if ( i == currentVectorDevice ) {
            strncat ( buf, "(=default)", sizeof(buf) - strlen(buf) - 1 );
          }
          len += strlen(buf);
          XLAL_CHECK_NULL ( len < sizeof(helpstr), XLAL_EBADLEN, "VectorDevice help-string exceeds buffer length (%lu)\n", sizeof(helpstr) );
          strcat ( helpstr, buf );
        } // for i < VECTORDEVICE_END

      strcat(helpstr, ") ");
      firstCall = 0;

    } // if firstCall

  return helpstr;
} // XLALVectorDeviceHelpString()

///
/// Parse a given string into an #VectorDevice_type if valid and available,
/// return error otherwise.
///
int
XLALVectorDeviceParseString ( VectorDevice_type *device,//!< [out] Parsed #VectorDevice_type enum
                              const char *s		//!< [in] String to parse
                              )
{
  XLAL_CHECK ( s != NULL, XLAL_EINVAL );
  XLAL_CHECK ( device != NULL, XLAL_EINVAL );

  // find matching VectorDevice name string
  for (int i = VECTORDEVICE_START + 1; i < VECTORDEVICE_END; i++ )
    {
      if ( (allVectorDevices[i].name != NULL) && (strcmp ( s, allVectorDevices[i].name ) == 0) )
        {
          if ( allVectorDevices[i].available )
            {
              (*device) = i;
              return XLAL_SUCCESS;
            }
          else
            {
              XLAL_ERROR ( XLAL_EINVAL, "Chosen VectorDevice '%s' valid but unavailable in this binary\n", s );
            }
        } // if found matching VectorDevice
    } // for i < VECTORDEVICE_END

  XLAL_ERROR ( XLAL_EINVAL, "Unknown VectorDevice '%s'\n", s );

} // XLALVectorDeviceParseString()


// ---------- High-level vector maths functions ----------
/// Vector version of 1-output math functions to take advantage of parallel-math devices (SSE,AVX,GPU...)

#define DEFINE_FN_1OUT(funcf)                                           \
  int XLALVector##funcf ( REAL4Vector *out, const REAL4Vector *in ) {   \
    switch ( currentVectorDevice ) {                                    \
    case VECTORDEVICE_FPU:                                              \
      return XLALVector##funcf##_FPU ( out, in );                       \
      break;                                                            \
    case VECTORDEVICE_SSE:                                              \
      return XLALVector##funcf##_SSE ( out, in );                       \
      break;                                                            \
    case VECTORDEVICE_AVX:                                              \
      return XLALVector##funcf##_AVX ( out, in );                       \
      break;                                                            \
    default:                                                            \
      XLAL_ERROR ( XLAL_EINVAL, "Invalid device selected '%d' unknown! [coding error]\n", currentVectorDevice ); \
      break;                                                            \
    } /* end: switch currentVectorDevice */                             \
    return XLAL_SUCCESS;                                                \
  } /* XLALVector##funcf() */

#define DEFINE_FN_2OUT(funcf)                                           \
  int XLALVector##funcf ( REAL4Vector *out, REAL4Vector *out2, const REAL4Vector *in ) { \
    switch ( currentVectorDevice ) {                                    \
    case VECTORDEVICE_FPU:                                              \
      return XLALVector##funcf##_FPU ( out, out2, in );                 \
      break;                                                            \
    case VECTORDEVICE_SSE:                                              \
      return XLALVector##funcf##_SSE ( out, out2, in );                 \
      break;                                                            \
    case VECTORDEVICE_AVX:                                              \
      return XLALVector##funcf##_AVX ( out, out2, in );                 \
      break;                                                            \
    default:                                                            \
      XLAL_ERROR ( XLAL_EINVAL, "Invalid device selected '%d' unknown! [coding error]\n", currentVectorDevice ); \
      break;                                                            \
    } /* end: switch currentVectorDevice */                             \
    return XLAL_SUCCESS;                                                \
  } /* XLALVector##funcf() */

DEFINE_FN_1OUT(Sinf);
DEFINE_FN_1OUT(Cosf);
DEFINE_FN_1OUT(Expf);
DEFINE_FN_1OUT(Logf);

DEFINE_FN_2OUT(SinCosf);
DEFINE_FN_2OUT(SinCosf2PI);

// -------------------- our own failsafe aligned memory handling --------------------
///
/// Create a special REAL4 Vector with n-byte aligned memory 'data' array.
///
/// This does not rely on posix_memalign() being available, and should compile+run everywhere.
/// Use XLALDestroyREAL4VectorAligned() to free this.
///
REAL4VectorAligned *
XLALCreateREAL4VectorAligned ( UINT4 length, UINT4 align )
{
  REAL4VectorAligned *ret;
  XLAL_CHECK_NULL ( (ret = XLALMalloc ( sizeof(*ret))) != NULL, XLAL_ENOMEM );

  ret->length = length;
  UINT4 paddedLength = length + align - 1;
  XLAL_CHECK_NULL ( (ret->data0 = XLALMalloc ( paddedLength * sizeof(REAL4) )) != NULL, XLAL_ENOMEM );

  size_t remBytes = ((size_t)ret->data0) % align;
  size_t offsetBytes = (align - remBytes) % align;
  ret->data = (void*)(((char*)ret->data0) + offsetBytes);

  XLAL_CHECK_NULL ( isMemAligned(ret->data,align), XLAL_EFAULT, "Failed to allocate %zd-byte aligned memory. Must be a coding error.\n", (size_t)align );

  return ret;
} // XLALCreateREAL4VectorAligned()

///
/// Destroy special n-byte aligned  REAL4VectorAligned
///
void
XLALDestroyREAL4VectorAligned ( REAL4VectorAligned *in )
{
  if ( !in ) {
    return;
  }
  if ( in->data0 ) {
    XLALFree ( in->data0 );
  }

  XLALFree ( in );

  return;

} // XLALDestroyREAL4VectorAligned()
