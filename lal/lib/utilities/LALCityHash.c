// Copyright (c) 2011 Google, Inc.
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
//
// CityHash, by Geoff Pike and Jyrki Alakuijala
//
// This file provides CityHash64() and related functions.
//
// It's probably possible to create even faster hash functions by
// writing a program that systematically explores some of the space of
// possible hash functions, by using SIMD instructions, or by
// compromising on hash quality.

// Converted to C by John Veitch

#include <config.h>
#include <lal/LALHashFunc.h>

#include <stdlib.h>  // for size_t.
#include <stdint.h>
#include <string.h>  // for memcpy and memset

typedef struct tagUINT16 {UINT8 first,second;} UINT16;

static inline UINT8 Uint128Low64(const UINT16 *x) { return x->first; }
static inline UINT8 Uint128High64(const UINT16 *x) { return x->second; }

/** Hash 128 input bits down to 64 bits of output.
 * This is intended to be a reasonably good hash function.
 */
static inline UINT8 Hash128to64(const UINT16 *x) {
  // Murmur-inspired hashing.
  const UINT8 kMul = 0x9ddfea08eb382d69ULL;
  UINT8 a = (Uint128Low64(x) ^ Uint128High64(x)) * kMul;
  a ^= (a >> 47);
  UINT8 b = (Uint128High64(x) ^ a) * kMul;
  b ^= (b >> 47);
  b *= kMul;
  return b;
}

static UINT8 UNALIGNED_LOAD64(const char *p) {
  UINT8 result;
  memcpy(&result, p, sizeof(result));
  return result;
}

static UINT4 UNALIGNED_LOAD32(const char *p) {
  UINT4 result;
  memcpy(&result, p, sizeof(result));
  return result;
}

#if defined(_WIN32) && !defined(__CYGWIN__)

#include <stdlib.h>
#define bswap_32(x) _byteswap_ulong(x)
#define bswap_64(x) _byteswap_UINT8(x)

#elif defined(__APPLE__)

// Mac OS X / Darwin features
#include <libkern/OSByteOrder.h>
#define bswap_32(x) OSSwapInt32(x)
#define bswap_64(x) OSSwapInt64(x)

#elif defined(__sun) || defined(sun)

#include <sys/byteorder.h>
#define bswap_32(x) BSWAP_32(x)
#define bswap_64(x) BSWAP_64(x)

#elif defined(__FreeBSD__)

#include <sys/endian.h>
#define bswap_32(x) bswap32(x)
#define bswap_64(x) bswap64(x)

#elif defined(__OpenBSD__)

#include <sys/types.h>
#define bswap_32(x) swap32(x)
#define bswap_64(x) swap64(x)

#elif defined(__NetBSD__)

#include <sys/types.h>
#include <machine/bswap.h>
#if defined(__BSWAP_RENAME) && !defined(__bswap_32)
#define bswap_32(x) bswap32(x)
#define bswap_64(x) bswap64(x)
#endif

#else

#include <byteswap.h>

#endif

#ifdef WORDS_BIGENDIAN
#define UINT4_in_expected_order(x) (bswap_32(x))
#define UINT8_in_expected_order(x) (bswap_64(x))
#else
#define UINT4_in_expected_order(x) (x)
#define UINT8_in_expected_order(x) (x)
#endif

#if !defined(LIKELY)
#if HAVE_BUILTIN_EXPECT
#define LIKELY(x) (__builtin_expect(!!(x), 1))
#else
#define LIKELY(x) (x)
#endif
#endif

static UINT8 Fetch64(const char *p) {
  return UINT8_in_expected_order(UNALIGNED_LOAD64(p));
}

static UINT4 Fetch32(const char *p) {
  return UINT4_in_expected_order(UNALIGNED_LOAD32(p));
}

// Some primes between 2^63 and 2^64 for various uses.
static const UINT8 k0 = 0xc3a5c85c97cb3127ULL;
static const UINT8 k1 = 0xb492b66fbe98f273ULL;
static const UINT8 k2 = 0x9ae16a3b2f90404fULL;

// Magic numbers for 32-bit hashing.  Copied from Murmur3.
static const UINT4 c1 = 0xcc9e2d51;
static const UINT4 c2 = 0x1b873593;

// A 32-bit to 32-bit integer hash copied from Murmur3.
static UINT4 fmix(UINT4 h)
{
  h ^= h >> 16;
  h *= 0x85ebca6b;
  h ^= h >> 13;
  h *= 0xc2b2ae35;
  h ^= h >> 16;
  return h;
}

static UINT4 Rotate32(UINT4 val, int shift) {
  // Avoid shifting by 32: doing so yields an undefined result.
  return shift == 0 ? val : ((val >> shift) | (val << (32 - shift)));
}

#define SWAP(a,b) {a^=b; b^=a; a^=b;}

#undef PERMUTE3
#define PERMUTE3(a, b, c) do { SWAP(a, b); SWAP(a, c); } while (0)

static UINT4 Mur(UINT4 a, UINT4 h) {
  // Helper from Murmur3 for combining two 32-bit values.
  a *= c1;
  a = Rotate32(a, 17);
  a *= c2;
  h ^= a;
  h = Rotate32(h, 19);
  return h * 5 + 0xe6546b64;
}

static UINT4 Hash32Len13to24(const char *s, size_t len) {
  UINT4 a = Fetch32(s - 4 + (len >> 1));
  UINT4 b = Fetch32(s + 4);
  UINT4 c = Fetch32(s + len - 8);
  UINT4 d = Fetch32(s + (len >> 1));
  UINT4 e = Fetch32(s);
  UINT4 f = Fetch32(s + len - 4);
  UINT4 h = len;
  
  return fmix(Mur(f, Mur(e, Mur(d, Mur(c, Mur(b, Mur(a, h)))))));
}

static UINT4 Hash32Len0to4(const char *s, size_t len) {
  UINT4 b = 0;
  UINT4 c = 9;
  for (size_t i = 0; i < len; i++) {
    signed char v = s[i];
    b = b * c1 + v;
    c ^= b;
  }
  return fmix(Mur(b, Mur(len, c)));
}

static UINT4 Hash32Len5to12(const char *s, size_t len) {
  UINT4 a = len, b = len * 5, c = 9, d = b;
  a += Fetch32(s);
  b += Fetch32(s + len - 4);
  c += Fetch32(s + ((len >> 1) & 4));
  return fmix(Mur(c, Mur(b, Mur(a, d))));
}

UINT4 XLALCityHash32(const char *s, size_t len) {
  if (len <= 24) {
    return len <= 12 ?
    (len <= 4 ? Hash32Len0to4(s, len) : Hash32Len5to12(s, len)) :
    Hash32Len13to24(s, len);
  }
  
  // len > 24
  UINT4 h = len, g = c1 * len, f = g;
  UINT4 a0 = Rotate32(Fetch32(s + len - 4) * c1, 17) * c2;
  UINT4 a1 = Rotate32(Fetch32(s + len - 8) * c1, 17) * c2;
  UINT4 a2 = Rotate32(Fetch32(s + len - 16) * c1, 17) * c2;
  UINT4 a3 = Rotate32(Fetch32(s + len - 12) * c1, 17) * c2;
  UINT4 a4 = Rotate32(Fetch32(s + len - 20) * c1, 17) * c2;
  h ^= a0;
  h = Rotate32(h, 19);
  h = h * 5 + 0xe6546b64;
  h ^= a2;
  h = Rotate32(h, 19);
  h = h * 5 + 0xe6546b64;
  g ^= a1;
  g = Rotate32(g, 19);
  g = g * 5 + 0xe6546b64;
  g ^= a3;
  g = Rotate32(g, 19);
  g = g * 5 + 0xe6546b64;
  f += a4;
  f = Rotate32(f, 19);
  f = f * 5 + 0xe6546b64;
  size_t iters = (len - 1) / 20;
  do {
    a0 = Rotate32(Fetch32(s) * c1, 17) * c2;
    a1 = Fetch32(s + 4);
    a2 = Rotate32(Fetch32(s + 8) * c1, 17) * c2;
    a3 = Rotate32(Fetch32(s + 12) * c1, 17) * c2;
    a4 = Fetch32(s + 16);
    h ^= a0;
    h = Rotate32(h, 18);
    h = h * 5 + 0xe6546b64;
    f += a1;
    f = Rotate32(f, 19);
    f = f * c1;
    g += a2;
    g = Rotate32(g, 18);
    g = g * 5 + 0xe6546b64;
    h ^= a3 + a1;
    h = Rotate32(h, 19);
    h = h * 5 + 0xe6546b64;
    g ^= a4;
    g = bswap_32(g) * 5;
    h += a4 * 5;
    h = bswap_32(h);
    f += a0;
    PERMUTE3(f, h, g);
    s += 20;
  } while (--iters != 0);
  g = Rotate32(g, 11) * c1;
  g = Rotate32(g, 17) * c1;
  f = Rotate32(f, 11) * c1;
  f = Rotate32(f, 17) * c1;
  h = Rotate32(h + g, 19);
  h = h * 5 + 0xe6546b64;
  h = Rotate32(h, 17) * c1;
  h = Rotate32(h + f, 19);
  h = h * 5 + 0xe6546b64;
  h = Rotate32(h, 17) * c1;
  return h;
}

// Bitwise right rotate.  Normally this will compile to a single
// instruction, especially if the shift is a manifest constant.
static UINT8 Rotate(UINT8 val, int shift) {
  // Avoid shifting by 64: doing so yields an undefined result.
  return shift == 0 ? val : ((val >> shift) | (val << (64 - shift)));
}

static UINT8 ShiftMix(UINT8 val) {
  return val ^ (val >> 47);
}

static UINT8 HashLen16(UINT8 u, UINT8 v) {
  UINT16 uv={.first=u,.second=v};
  return Hash128to64(&uv);
}

static UINT8 HashLen16mul(UINT8 u, UINT8 v, UINT8 mul) {
  // Murmur-inspired hashing.
  UINT8 a = (u ^ v) * mul;
  a ^= (a >> 47);
  UINT8 b = (v ^ a) * mul;
  b ^= (b >> 47);
  b *= mul;
  return b;
}

static UINT8 HashLen0to16(const char *s, size_t len) {
  if (len >= 8) {
    UINT8 mul = k2 + len * 2;
    UINT8 a = Fetch64(s) + k2;
    UINT8 b = Fetch64(s + len - 8);
    UINT8 c = Rotate(b, 37) * mul + a;
    UINT8 d = (Rotate(a, 25) + b) * mul;
    return HashLen16mul(c, d, mul);
  }
  if (len >= 4) {
    UINT8 mul = k2 + len * 2;
    UINT8 a = Fetch32(s);
    return HashLen16mul(len + (a << 3), Fetch32(s + len - 4), mul);
  }
  if (len > 0) {
    UCHAR a = s[0];
    UCHAR b = s[len >> 1];
    UCHAR c = s[len - 1];
    UINT4 y = (UINT4)(a) + ((UINT4)(b) << 8);
    UINT4 z = len + ((UINT4)(c) << 2);
    return ShiftMix(y * k2 ^ z * k0) * k2;
  }
  return k2;
}

// This probably works well for 16-byte strings as well, but it may be overkill
// in that case.
static UINT8 HashLen17to32(const char *s, size_t len) {
  UINT8 mul = k2 + len * 2;
  UINT8 a = Fetch64(s) * k1;
  UINT8 b = Fetch64(s + 8);
  UINT8 c = Fetch64(s + len - 8) * mul;
  UINT8 d = Fetch64(s + len - 16) * k2;
  return HashLen16mul(Rotate(a + b, 43) + Rotate(c, 30) + d,
                   a + Rotate(b + k2, 18) + c, mul);
}

// Return a 16-byte hash for 48 bytes.  Quick and dirty.
// Callers do best to use "random-looking" values for a and b.
static UINT16 WeakHashLen32WithSeeds_5(UINT8 w, UINT8 x, UINT8 y, UINT8 z, UINT8 a, UINT8 b)
{
  a += w;
  b = Rotate(b + a + z, 21);
  UINT8 c = a;
  a += x;
  a += y;
  b += Rotate(a, 44);
  UINT16 pair={.first=a+z,.second=b+c};
  return pair;
}

// Return a 16-byte hash for s[0] ... s[31], a, and b.  Quick and dirty.
static UINT16 WeakHashLen32WithSeeds(const char* s, UINT8 a, UINT8 b) {
  return WeakHashLen32WithSeeds_5(Fetch64(s),
                                Fetch64(s + 8),
                                Fetch64(s + 16),
                                Fetch64(s + 24),
                                a,
                                b);
}

// Return an 8-byte hash for 33 to 64 bytes.
static UINT8 HashLen33to64(const char *s, size_t len) {
  UINT8 mul = k2 + len * 2;
  UINT8 a = Fetch64(s) * k2;
  UINT8 b = Fetch64(s + 8);
  UINT8 c = Fetch64(s + len - 24);
  UINT8 d = Fetch64(s + len - 32);
  UINT8 e = Fetch64(s + 16) * k2;
  UINT8 f = Fetch64(s + 24) * 9;
  UINT8 g = Fetch64(s + len - 8);
  UINT8 h = Fetch64(s + len - 16) * mul;
  UINT8 u = Rotate(a + g, 43) + (Rotate(b, 30) + c) * 9;
  UINT8 v = ((a + g) ^ d) + f + 1;
  UINT8 w = bswap_64((u + v) * mul) + h;
  UINT8 x = Rotate(e + f, 42) + c;
  UINT8 y = (bswap_64((v + w) * mul) + g) * mul;
  UINT8 z = e + f + c;
  a = bswap_64((x + z) * mul + y) + b;
  b = ShiftMix((z + a) * mul + d + h) * mul;
  return b + x;
}

UINT8 XLALCityHash64(const char *s, size_t len) {
  if (len <= 32) {
    if (len <= 16) {
      return HashLen0to16(s, len);
    } else {
      return HashLen17to32(s, len);
    }
  } else if (len <= 64) {
    return HashLen33to64(s, len);
  }
  
  // For strings over 64 bytes we hash the end first, and then as we
  // loop we keep 56 bytes of state: v, w, x, y, and z.
  UINT8 x = Fetch64(s + len - 40);
  UINT8 y = Fetch64(s + len - 16) + Fetch64(s + len - 56);
  UINT8 z = HashLen16(Fetch64(s + len - 48) + len, Fetch64(s + len - 24));
  UINT16 v = WeakHashLen32WithSeeds(s + len - 64, len, z);
  UINT16 w = WeakHashLen32WithSeeds(s + len - 32, y + k1, x);
  x = x * k1 + Fetch64(s);
  
  // Decrease len to the nearest multiple of 64, and operate on 64-byte chunks.
  len = (len - 1) & ~((size_t)63);
  do {
    x = Rotate(x + y + v.first + Fetch64(s + 8), 37) * k1;
    y = Rotate(y + v.second + Fetch64(s + 48), 42) * k1;
    x ^= w.second;
    y += v.first + Fetch64(s + 40);
    z = Rotate(z + w.first, 33) * k1;
    v = WeakHashLen32WithSeeds(s, v.second * k1, x + w.first);
    w = WeakHashLen32WithSeeds(s + 32, z + w.second, y + Fetch64(s + 16));
    SWAP(z, x);
    s += 64;
    len -= 64;
  } while (len != 0);
  return HashLen16(HashLen16(v.first, w.first) + ShiftMix(y) * k1 + z,
                   HashLen16(v.second, w.second) + x);
}

UINT8 XLALCityHash64WithSeed(const char *s, size_t len, UINT8 seed) {
  return XLALCityHash64WithSeeds(s, len, k2, seed);
}

UINT8 XLALCityHash64WithSeeds(const char *s, size_t len,
                           UINT8 seed0, UINT8 seed1) {
  return HashLen16(XLALCityHash64(s, len) - seed0, seed1);
}
