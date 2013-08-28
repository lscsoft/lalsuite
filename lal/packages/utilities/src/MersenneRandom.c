/*
*  Copyright (C) 2007 Jolien Creighton
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


#include <lal/LALStdlib.h>
#include <lal/Random.h>

/*  Constants  */
#define N 624
#define M 397
#define MATRIX_A 0x9908b0df     /* constant vector a */
#define UPPER_MASK 0x80000000   /* most significant w-r bits */
#define LOWER_MASK 0x7fffffff   /* least significant r bits */

/*  Predefined Macros  */
#define TEMPERING_MASK_B 0x9d2c5680
#define TEMPERING_MASK_C 0xefc60000
#define TEMPERING_SHIFT_U(y)  (y >> 11)
#define TEMPERING_SHIFT_S(y)  (y << 7)
#define TEMPERING_SHIFT_T(y)  (y << 15)
#define TEMPERING_SHIFT_L(y)  (y >> 18)

/**
 * \defgroup MersenneRandom_c Module MersenneRandom.c
 * \ingroup Random_h
 * \author Tibbits, M M
 *
 * \brief Routine to get a random number based on the Mersenne Twister Algorithm.
 *
 * ### Description ###
 *
 * This was implemented as in the paper listed below.  We have provided two functions,
 * each which may be called multiple times.  One returns a single number and the other
 * returns a vector of length prescribed by the Vector-\>length.
 *
 * Below I have listed the abstract from the paper referenced below:
 *
 * <p>
 * A new algorithm called Mersenne Twister (MT) is proposed for generating uniform
 * pseudoran-dom numbers. For a particular choice of parameters, the algorithm provides
 * a super astronom-ical period of 2 19937 2 1 and 623-dimensional equidistribution up
 * to 32-bit accuracy, while using a working area of only 624 words. This is a new
 * variant of the previously proposed generators, TGFSR, modified so as to admit a
 * Mersenne-prime period. The characteristic polynomial has many terms. The distribution
 * up to v bits accuracy for 1 \# v \# 32 is also shown to be good. An algorithm is also
 * given that checks the primitivity of the characteristic polynomial of MT with
 * computational complexity O( p 2 ) where p is the degree of the polynomial.
 *
 * <p>
 * They implemented this generator in portable C-code. It passed several stringent
 * statistical tests, including diehard. Its speed is comparable to other modern
 * generators. Its merits are due to the efficient algorithms that are unique to
 * polynomial calculations over the two-element field.
 *
 * ### Algorithm ###
 *
 * Please see paper listed below:
 *
 * M. Matsumoto and T. Nishimura,
 * "Mersenne Twister: A 623-Dimensionally Equidistributed Uniform
 * Pseudo-Random Number Generator",
 * ACM Transactions on Modeling and Computer Simulation,
 * Vol. 8, No. 1, January 1998, pp 3--30.
 *
 * Copyright (C) 1997 Makoto Matsumoto and Takuji Nishimura.
 * When you use this, send an email to: matumoto@math.keio.ac.jp
 * with an appropriate reference to your work.
 *
 * ### Notes ###
 *
 * Pulled from distributed source code:
 *
 * A C-program for MT19937: Real number version  (1998/4/6)
 * genrand() generates one pseudorandom real number (double)
 * which is uniformly distributed on [0,1]-interval, for each
 * call. sgenrand(seed) set initial values to the working area
 * of 624 words. Before genrand(), sgenrand(seed) must be called
 * once. (seed is any 32-bit integer except for 0).  Integer
 * generator is obtained by modifying two lines.
 *
 * Coded by Takuji Nishimura, considering the suggestions by
 * Topher Cooper and Marc Rieffel in July-Aug. 1997.
 *
 * Seed value MAY NOT EQUAL ZERO.
 *
 */
/*@{*/

typedef struct
tagGenParams
{
        UINT4           *mt;    /* the array for the state vector  */
        INT4            mti;    /* mti==N+1 means mt[N] is not initialized */
        UINT4           mag01[2];
        BOOLEAN         set;
}
GenParams;

/**
 * \ingroup Random_h
 * \brief This structure contains the parameters necessary for generating the current
 * sequence of Mersenne twiser random numbers (based on the initial seed).
 * \note The contents should not be manually adjusted.
 */
struct
tagMTRandomParams
{
        UINT4           seed;
        INT2            initialized;
        GenParams       *priv;
};


static void sgenrand( unsigned long seed, GenParams *params );
static double genrand( GenParams *params );


void LALCreateMTRandomParams
(
	LALStatus	*status,
	REAL8		seed,
	MTRandomParams	**params
)
{

	INITSTATUS(status);

	ASSERT ( seed != 0, status,  RANDOMH_ESEED, RANDOMH_MSGESEED);

	/*  params must come in allocated  */
	ASSERT ( params != NULL, status, RANDOMH_ENULL, RANDOMH_MSGENULL);

	ASSERT ( *params == NULL, status, RANDOMH_ENULL, RANDOMH_MSGENULL);

	*params = (MTRandomParams*) LALMalloc(sizeof(MTRandomParams));

	ASSERT ( *params, status, RANDOMH_EALOC, RANDOMH_MSGEALOC);

	(*params)->priv = (GenParams*) LALMalloc(sizeof(GenParams));

	ASSERT ( (*params)->priv, status, RANDOMH_EALOC, RANDOMH_MSGEALOC);

	(*params)->priv->mt = (UINT4*) LALMalloc(sizeof(UINT4)*N);

	ASSERT ( (*params)->priv->mt, status, RANDOMH_EALOC, RANDOMH_MSGEALOC);

	(*params)->priv->mti      = N + 1;
	(*params)->priv->mag01[0] = 0x0;
	(*params)->priv->mag01[1] = MATRIX_A;

	(*params)->seed = seed;
	(*params)->initialized = 1;

	RETURN (status);
}


void LALDestroyMTRandomParams
(
	LALStatus	*status,
	MTRandomParams	**params
)
{

	INITSTATUS(status);

	/*  Check input for existence.  */
	/*  params must come in allocated  */
	ASSERT ( (*params), status, RANDOMH_ENULL, RANDOMH_MSGENULL);

	/*  Params must come in initialized  */
	ASSERT ( (*params)->initialized, status, RANDOMH_EINIT, RANDOMH_MSGEINIT);

	LALFree((*params)->priv->mt);
	LALFree((*params)->priv);
	LALFree(*params);

	*params = NULL;

	RETURN (status);
}


void LALMersenneRandom
(
	LALStatus	*status,
	REAL8		*output,
	MTRandomParams	*params
)
{

	INITSTATUS(status);

	/*  Check input for existence.  */
	/*  Output should come in Allocated  */
	ASSERT ( output, status, RANDOMH_ENULL, RANDOMH_MSGENULL);

	/*  params must come in allocated  */
	ASSERT ( params, status, RANDOMH_ENULL, RANDOMH_MSGENULL);

	ASSERT ( params->initialized == 1, status, RANDOMH_EINIT, RANDOMH_MSGEINIT);

	if (params->priv->set == 0)
	{
		sgenrand( params->seed, params->priv );
		params->priv->set = 1;
	}

	*output = genrand( params->priv );

	RETURN (status);
}


void LALMersenneRandomVector
(
	LALStatus	*status,
	REAL8Vector	*output,
	MTRandomParams	*params
)
{
	/*  Variable Declaration  */
	UINT8	iterator;

	INITSTATUS(status);

	/*  Variable Initialization  */
	iterator = 0;


	/*  Check input for existence.  */
	/*  Output should come in Allocated  */
	ASSERT ( output, status, RANDOMH_ENULL, RANDOMH_MSGENULL);

	/*  params must come in allocated  */
	ASSERT ( params, status, RANDOMH_ENULL, RANDOMH_MSGENULL);

	/*  Params must come in initialized  */
	ASSERT ( params->initialized, status, RANDOMH_EINIT, RANDOMH_MSGEINIT);

	/*  Vector length must be greater than zero  */
	ASSERT ( output->length > 0, status, RANDOMH_EZERO, RANDOMH_MSGEZERO);

	if (params->priv->set == 0)
	{
		sgenrand(  params->seed, params->priv );
		params->priv->set = 1;
	}

	for (iterator = 0; iterator < output->length; iterator++ )
	{
		output->data[iterator] = genrand( params->priv );
	}

	RETURN (status);
}



/* A C-program for MT19937: Real number version  (1998/4/6)    */
/*   genrand() generates one pseudorandom real number (double) */
/* which is uniformly distributed on [0,1]-interval, for each  */
/* call. sgenrand(seed) set initial values to the working area */
/* of 624 words. Before genrand(), sgenrand(seed) must be      */
/* called once. (seed is any 32-bit integer except for 0).     */
/* Integer generator is obtained by modifying two lines.       */
/*   Coded by Takuji Nishimura, considering the suggestions by */
/* Topher Cooper and Marc Rieffel in July-Aug. 1997.           */

/* This library is free software; you can redistribute it and/or   */
/* modify it under the terms of the GNU Library General Public     */
/* License as published by the Free Software Foundation; either    */
/* version 2 of the License, or (at your option) any later         */
/* version.                                                        */
/* This library is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of  */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.            */
/* See the GNU Library General Public License for more details.    */
/* You should have received a copy of the GNU Library General      */
/* Public License along with this library; if not, write to the    */
/* Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA   */
/* 02111-1307  USA                                                 */

/* Copyright (C) 1997 Makoto Matsumoto and Takuji Nishimura.       */
/* When you use this, send an email to: matumoto@math.keio.ac.jp   */
/* with an appropriate reference to your work.                     */

/* REFERENCE                                                       */
/* M. Matsumoto and T. Nishimura,                                  */
/* "Mersenne Twister: A 623-Dimensionally Equidistributed Uniform  */
/* Pseudo-Random Number Generator",                                */
/* ACM Transactions on Modeling and Computer Simulation,           */
/* Vol. 8, No. 1, January 1998, pp 3--30.                          */

/* static int mti=N+1; */ /* mti==N+1 means mt[N] is not initialized */

static void sgenrand
(
	unsigned long	seed,
	GenParams	*params
)
{

	/* setting initial seeds to mt[N] using         */
	/* the generator Line 25 of Table 1 in          */
	/* [KNUTH 1981, The Art of Computer Programming */
	/*    Vol. 2 (2nd Ed.), pp102]                  */

	params->mt[0]= seed & 0xffffffff;

	for (params->mti = 1; params->mti < N; params->mti++)
	{
		params->mt[params->mti] = (69069 * params->mt[params->mti-1]) & 0xffffffff;
	}
}

double genrand
(
	GenParams	*params
)
{
	unsigned long	y;

	if (params->mti >= N)
	{
		/* generate N words at one time */
		int kk;

		/* if sgenrand() has not been called, */
		if (params->mti == N+1)
		{
			/* a default initial seed is used   */
			sgenrand(4357, params);
		}

		for ( kk = 0; kk < (N - M); kk++)
		{
			y = (params->mt[kk]&UPPER_MASK)|(params->mt[kk+1]&LOWER_MASK);
			params->mt[kk] = params->mt[kk+M] ^ (y >> 1) ^ params->mag01[y & 0x1];
		}

		for (       ; kk < (N - 1); kk++)
		{
			y = (params->mt[kk]&UPPER_MASK)|(params->mt[kk+1]&LOWER_MASK);
			params->mt[kk] = params->mt[kk+(M-N)] ^ (y >> 1) ^ params->mag01[y & 0x1];
		}

		y = (params->mt[N-1]&UPPER_MASK)|(params->mt[0]&LOWER_MASK);
		params->mt[N-1] = params->mt[M-1] ^ (y >> 1) ^ params->mag01[y & 0x1];

		params->mti = 0;
	}

	y = params->mt[params->mti++];
	y ^= TEMPERING_SHIFT_U(y);
	y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
	y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C;
	y ^= TEMPERING_SHIFT_L(y);

    return ( (double)y * 2.3283064370807974e-10 );
}
/*@}*/
