/************************************ <lalVerbatim file="NRRandomCV">
Author: Tibbits, M M
$Id$
************************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Module \texttt{NRRandom.c}}
\label{s:NRRandom.c}

\subsubsection*{Prototypes}
\input{NRRandomCP}

\subsubsection*{Description}

Above are listed the random number generating functions
extracted from Numerical Recipes in C Second Edition which have been
wrapped in LAL code to ASSERT errors in passing NULL pointers and taking care
of the status variable.

\subsubsection*{Algorithm}
The algorythms are on separate pages of Numerical Recipes in C, and are listed on the
following \\
\begin{tabbing}
pages:\=
\\
\>ran3   - p. 283\\
\>ran4   - p. 303 - 304\\
\>gasdev - p. 289\\
\end{tabbing}

\subsubsection*{Uses}

Call \texttt{NRRan3()} with the appropriate arguments and it will execute the function ran3
from the Numerical Recipes Book.  Call \texttt{NRRan4()} with the appropriate arguments and
it will execute the function ran4 from the Numerical Recipes Book.  Call \texttt{NRGasDev()}
with the appropriate arguments and it will execute the function gasdev from the Numerical 
Recipes Book.

\subsubsection*{Notes}

Code taken from Numerical Recipes in C.  Version 2.08

\vfill{\footnotesize\input{NRRandomCV}}

******************************************************* </lalLaTeX> */

#include <math.h>
#include <lal/NRRandom.h>

NRCSID( NRRANDOMC, "$Id$");

/* <lalVerbatim file="NRRandomCP"> */
void NRRan3
(
	LALStatus	*status,
	REAL4		*output,
	INT4		*seed
)
{	/* </lalVerbatim > */
	/*  Variable Declaration  */

	INITSTATUS( status, "NRRan3" , NRRANDOMC);

	/*  Check input for existence.  */
	/*  Output should come in Allocated  */
	ASSERT ( output, status, NRRANDOMH_ENULL, NRRANDOMH_MSGENULL);

	/*  seed should come in Allocated  */
	ASSERT ( seed, status, NRRANDOMH_ENULL, NRRANDOMH_MSGENULL);

	*output = ran3( seed );

	RETURN( status );
}

/* <lalVerbatim file="NRRandomCP"> */
void NRRan4
(
	LALStatus	*status,
	REAL4		*output,
	INT4		*seed
)
{	/* </lalVerbatim > */
	/*  Variable Declaration  */

	INITSTATUS( status, "NRRan4" , NRRANDOMC);

	/*  Check input for existence.  */
	/*  Output should come in Allocated  */
	ASSERT ( output, status, NRRANDOMH_ENULL, NRRANDOMH_MSGENULL);

	/*  seed should come in Allocated  */
	ASSERT ( seed, status, NRRANDOMH_ENULL, NRRANDOMH_MSGENULL);

	*output = ran4( seed );

	RETURN( status );
}


/* <lalVerbatim file="NRRandomCP"> */
void NRGasDev
(
	LALStatus	*status,
	REAL4		*output,
	INT4		*seed,
	INT2		*param
)
{	/* </lalVerbatim > */
	/*  Variable Declaration  */
	const	RanFunction	whichTest[2]   = { ran3, ran4 };

	INITSTATUS( status, "RAN4" , NRRANDOMC);

	/*  Check input for existence.  */
	/*  Output should come in Allocated  */
	ASSERT ( output, status, NRRANDOMH_ENULL, NRRANDOMH_MSGENULL);

	/*  seed should come in Allocated  */
	ASSERT ( seed, status, NRRANDOMH_ENULL, NRRANDOMH_MSGENULL);

	/*  Output should come in Allocated  */
	ASSERT ( param, status, NRRANDOMH_ENULL, NRRANDOMH_MSGENULL);

	/*  Check params bounds  */
	ASSERT ( (((*param) == 0)||((*param) == 1)) , status, NRRANDOMH_ETEST, NRRANDOMH_MSGETEST);

	*output = gasdev( seed, whichTest[*param] );

	RETURN( status );
}


float ran3(long *idum)
{
	static int inext,inextp;
	static long ma[56];
	static int iff=0;
	long mj,mk;
	int i,ii,k;
	
	if(*idum < 0 || iff == 0)
	{
		iff = 1;
		mj  = MSEED - (*idum < 0? -*idum : *idum);
		mj %=MBIG;
		ma[55]=mj;
		mk  = 1;
		for(i = 1; i <= 54; i++)
		{
			ii     = (21 * i) % 55;
			ma[ii] = mk;
			mk     = mj - mk;
			if(mk < MZ)
				mk += MBIG;
			mj     = ma[ii];
		}
		for(k = 1; k <= 4; k++)
			for(i = 1; i <= 55; i++)
			{
				ma[i] -= ma[1+(i+30)%55];
				if(ma[i]<MZ) 
					ma[i] += MBIG;
			}
		inext=0;
		inextp=31;
		*idum=1;
	}
	if(++inext == 56)
		inext=1;
	if(++inextp == 56)
		inextp=1;
	mj=ma[inext]-ma[inextp];
	if(mj < MZ)
		mj+=MBIG;
	ma[inext]=mj;
	return mj*FAC;
}


void psdes(unsigned long* lword, unsigned long *irword)
{
	unsigned long i,ia,ib,iswap,itmph=0,itmpl=0;
	static unsigned long c1[NITER]={
		0xbaa96887L, 0x1e17b32cL, 0x03bcdc3cL, 0x0f33d1b2L};
	static unsigned long c2[NITER]={
		0x4b0f3b58L, 0xe874f0c3L, 0x6955c5a6L, 0x55a7ca46L};
	
	for(i=0;i<NITER;i++)
	{
		ia=(iswap=(*irword)) ^ c1[i];
		itmpl = ia & 0xffff;
		itmph = ia >> 16;
		ib=itmpl*itmpl+ ~(itmph*itmph);
		*irword=*lword ^ (((ia=(ib>>16) | ((ib & 0xffff) << 16)) ^ c2[i])+itmpl*itmph);
		*lword=iswap;
	}
}

float ran4(long *idum)
{
	unsigned long irword,itemp,lword;
	static long idums =0;

#if defined(vax) || defined(_vax_) || defined(__vax__) || defined(VAX)
	static unsigned long jflone = 0x00004080;
	static unsigned long jflmsk = 0xffff007f;
#else
	static unsigned long jflone = 0x3f800000;
	static unsigned long jflmsk = 0x007fffff;
#endif

	if(*idum<0)
	{
		idums = -(*idum);
		*idum=1;
	}

	irword=(*idum);
	lword=idums;
	psdes(&lword,&irword);
	itemp=jflone | (jflmsk&irword);
	++(*idum);
	return (*(float *)&itemp)-1.0;
}


float gasdev(long *idum, RanFunction gen)
{
	static int iset=0;
	static float gset;
	float fac,rsq,v1,v2;

	if(iset == 0)
	{
		do{
			v1=2.0*gen(idum)-1.0;
			v2=2.0*gen(idum)-1.0;
			rsq=v1*v1+v2*v2;
		}while (rsq >= 1.0 || rsq == 0.0);
		fac  = sqrt(-2.0*log(rsq)/rsq);
		gset = v1*fac;
		iset = 1;
		return v2*fac;
	}
	else
	{
		iset=0;
		return gset;
	}
}
