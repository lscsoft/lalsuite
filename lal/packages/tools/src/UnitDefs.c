/************************************ <lalVerbatim file="UnitDefsCV">
Author: J. T. Whelan <jtwhelan@loyno.edu>
$Id$  
************************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Module \texttt{UnitDefs.c}}
\label{tools:ss:UnitDefs.c}

Defines basic and derived SI units and a function to produce a text
string corresponding to a unit structure.

\subsubsection*{Prototypes}
\input{UnitDefsCP}
\idx{LALUnitAsString()}
\idx{LALParseUnitString()}

\subsubsection*{Description}

\texttt{LALUnitAsString()} converts the unit structure
\texttt{*input} into a text string which is stored in the character
vector \texttt{*output}.  Note that the resulting text string is
expressed solely in terms of the basic units (m, kg, s, A,
K, strain and counts), and is thus not necessarily the most
convenient way to check the units of a quantity.  A better method is
to construct a unit structure containing the expected units, then
compare that to the actual units using \texttt{LALUnitCompare()}.

\texttt{LALParseUnitString()} reconstructs the original
\texttt{LALUnit} structure from the string output by
\texttt{LALUnitAsString()}.  It is very sensitive to the exact format
of the string and is not intended for use in parsing user-entered
strings.

\subsubsection*{Algorithm}

\texttt{LALUnitAsString()} moves through the unit structure, appending
the appropriate text to the string as it goes along.

\texttt{LALParseUnitString()} moves through the input string, one
character at a time, building an \texttt{LALUnit} structure as a it
goes along, so long as it encounters precisely the syntax expected.

\subsubsection*{Uses}
None.

\subsubsection*{Notes}

\subsubsection*{Predefined Units}
\idx[Constant]{lalMeterUnit}
\idx[Constant]{lalKiloGramUnit}
\idx[Constant]{lalSecondUnit}
\idx[Constant]{lalAmpereUnit} 
\idx[Constant]{lalKelvinUnit}
\idx[Constant]{lalStrainUnit} 
\idx[Constant]{lalADCCountUnit} 
\idx[Constant]{lalHertzUnit} 
\idx[Constant]{lalNewtonUnit} 
\idx[Constant]{lalPascalUnit} 
\idx[Constant]{lalJouleUnit} 
\idx[Constant]{lalWattUnit} 
\idx[Constant]{lalCoulombUnit} 
\idx[Constant]{lalVoltUnit} 
\idx[Constant]{lalOhmUnit} 
\idx[Constant]{lalFaradUnit} 
\idx[Constant]{lalWeberUnit} 
\idx[Constant]{lalHenryUnit} 
\idx[Constant]{lalTeslaUnit} 
\idx[Constant]{lalYottaUnit} 
\idx[Constant]{lalZettaUnit} 
\idx[Constant]{lalExaUnit} 
\idx[Constant]{lalPetaUnit} 
\idx[Constant]{lalTeraUnit} 
\idx[Constant]{lalGigaUnit} 
\idx[Constant]{lalMegaUnit} 
\idx[Constant]{lalKiloUnit} 
\idx[Constant]{lalHectoUnit} 
\idx[Constant]{lalDekaUnit} 
\idx[Constant]{lalDeciUnit} 
\idx[Constant]{lalCentiUnit} 
\idx[Constant]{lalMilliUnit} 
\idx[Constant]{lalMicroUnit} 
\idx[Constant]{lalNanoUnit} 
\idx[Constant]{lalPicoUnit} 
\idx[Constant]{lalFemtoUnit} 
\idx[Constant]{lalAttoUnit} 
\idx[Constant]{lalZeptoUnit} 
\idx[Constant]{lalYoctoUnit} 
\idx[Constant]{lalGramUnit} 
\idx[Constant]{lalAttoStrainUnit} 
\idx[Constant]{lalPicoFaradUnit} 

This file also defines a number of \texttt{constant} unit structures
(declared \texttt{extern} in \texttt{Units.h}).  Zeroth is
\texttt{lalDimensionlessUnit}, which is simply a \texttt{LALUnit} 
structure to be associated with a unitless quantity.
First, the relevant fundamental SI units and two custom units of use in 
gravitational wave detection:
\begin{center}
\begin{tabular}{|llll|}
\hline
Constant & Name & Abbr.\ & Physical Quantity \\
\hline
\texttt{lalMeterUnit} & meter & m & length \\
\texttt{lalKiloGramUnit} & kilogram & kg & mass \\
\texttt{lalSecondUnit} & second & s & time \\
\texttt{lalAmpereUnit} & ampere & A & electric current \\
\texttt{lalKelvinUnit} & kelvin & K & thermodynamic temperature \\
\texttt{lalStrainUnit} & strain & $\epsilon$ & gravitational strain \\
\texttt{lalADCCountUnit} & ADC count & count & A-to-D converter counts \\
\hline
\end{tabular}
\end{center}
Next, the named derived units in the SI\cite{tools:Halliday:2001}:
\begin{center}
\begin{tabular}{|llllll|}
\hline
Constant & Name & Abbr.\ & Physical Quantity & Def.\ & Fundamental\\
\hline
\texttt{lalHertzUnit} & hertz & Hz & frequency &
  s$^{-1}$ &  s$^{-1}$ \\
\texttt{lalNewtonUnit} & newton & N & force &
   kg$\cdot$ m/s$^2$ & m kg s$^{-2}$\\
\texttt{lalPascalUnit} & pascal & Pa & pressure &
   N/m$^2$ & m$^{-1}$ kg s$^{-2}$ \\
\texttt{lalJouleUnit} & joule & J & energy &
   N$\cdot$m & m$^2$ kg s$^{-2}$ \\
\texttt{lalWattUnit} & watt & W & power &
   J/s &  m$^2$ kg s$^{-3}$ \\
\texttt{lalCoulombUnit} & coulomb & C & electric charge & A$\cdot$s & s A \\
\texttt{lalVoltUnit} & volt & V & potential &
    W/A & m$^2$ kg s$^{-3}$ A$^{-1}$\\
\texttt{lalOhmUnit} & ohm & $\Omega$ & resistance &
    V/A & m$^2$ kg s$^{-3}$ A$^{-2}$\\
\texttt{lalFaradUnit} & farad & F & capacitance &
    C/V & m$^{-2}$ kg$^{-1}$ s$^4$ A$^2$\\
\texttt{lalWeberUnit} & weber & Wb & magnetic flux &
    V$\cdot$s & m$^2$ kg s$^{-2}$ A$^{-1}$\\
\texttt{lalHenryUnit} & henry & H & inductance &
    V$\cdot$s/A & m$^2$ kg s$^{-2}$ A$^{-2}$\\
\texttt{lalTeslaUnit} & tesla & T & magnetic flux density &
    Wb/m$^2$ & kg s$^{-2}$ A$^{-1}$\\
\hline
\end{tabular}
\end{center}
The powers of ten (SI prefixes)
\begin{center}
\begin{tabular}{|llll|}
\hline
Constant & Prefix & Abbr.\ & Value\\
\hline
\texttt{lalYottaUnit} & yotta & Y & $10^{ 24}$ \\
\texttt{lalZettaUnit} & zetta & Z & $10^{ 21}$ \\
\texttt{lalExaUnit} & exa & E & $10^{ 18}$ \\
\texttt{lalPetaUnit} & peta & P & $10^{ 15}$ \\
\texttt{lalTeraUnit} & tera & T & $10^{ 12}$ \\
\texttt{lalGigaUnit} & giga & G & $10^{  9}$ \\
\texttt{lalMegaUnit} & mega & M & $10^{  6}$ \\
\texttt{lalKiloUnit} & kilo & k & $10^{  3}$ \\
\texttt{lalHectoUnit} & hecto & h & $10^{  2}$ \\
\texttt{lalDekaUnit} & deka & da & $10^{  1}$ \\
\texttt{lalDeciUnit} & deci & d & $10^{ -1}$ \\
\texttt{lalCentiUnit} & centi & c & $10^{ -2}$ \\
\texttt{lalMilliUnit} & milli & m & $10^{ -3}$ \\
\texttt{lalMicroUnit} & micro & $\mu$ & $10^{ -6}$ \\
\texttt{lalNanoUnit} & nano & n & $10^{ -9}$ \\
\texttt{lalPicoUnit} & pico & p & $10^{-12}$ \\
\texttt{lalFemtoUnit} & femto & f & $10^{-15}$ \\
\texttt{lalAttoUnit} & atto & a & $10^{-18}$ \\
\texttt{lalZeptoUnit} & zepto & z & $10^{-21}$ \\
\texttt{lalYoctoUnit} & yocto & y & $10^{-24}$ \\
\hline
\end{tabular}
\end{center}
And finally a couple of convenient scaled units:
\begin{center}
\begin{tabular}{|lllll|}
\hline
Constant & Name & Abbr.\ & Def.\ & Fundamental\\
\hline
\texttt{lalGramUnit} & gram & g &
  $10^{-3}$ kg & $10^{-3}$ kg \\
\texttt{lalAttoStrainUnit} & attostrain & a$\epsilon$ &
  $10^{-18} \epsilon$ & $10^{-18} \epsilon$ \\
\texttt{lalPicoFaradUnit} & picofarad & pF & 
  $10^{-12}$ F & $10^{-12}$ m$^{-2}$ kg$^{-1}$ s$^4$ A$^2$\\
\hline
\end{tabular}
\end{center}

\vfill{\footnotesize\input{UnitDefsCV}}

******************************************************* </lalLaTeX> */ 
/**************************************** <lalLaTeX file="UnitDefsCB">
\bibitem{tools:Halliday:2001}
D.~Halliday, R.~Resnick, and J.~Walker, \textit{Fundamentals of
  Physics}.  (Wiley \& Sons, New York, 2001)
******************************************************* </lalLaTeX> */ 

#include <lal/LALStdlib.h>
#include <string.h>
#include <ctype.h>
#include <lal/Units.h>

NRCSID( UNITDEFSC, "$Id$" );

#define UNITDEFSC_TEMPSIZE 20

/* To convert a units structure to a string repesentation, we need to
 * define the names of the basic units.
 */

const CHAR lalUnitName[LALNumUnits][LALUnitNameSize] = 
{
  "m", "kg", "s", "A", "K", "strain", "count"
};

/*********************************************************
 *                                                       *
 *                 Predefined units                      *
 *                                                       *
 *********************************************************/

/* Predefined constant units make it easier for programmers to specify
 * and compare (using LALUnitCompare) units more easily.  Those given
 * here are an example; more can be added.
 */

/* LALUnitsTest.c will verify the definitions of the derived units,
 * for example using LALUnitRaise, LALUnitMultiply and LALUnitCompare
 * to show that 1 Farad = 1 Coulomb Volt^-1
 */

const LALUnit lalDimensionlessUnit = {  0, { 0, 0, 0, 0, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };

/* Basic Units */
const LALUnit lalMeterUnit         = {  0, { 1, 0, 0, 0, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };
const LALUnit lalKiloGramUnit      = {  0, { 0, 1, 0, 0, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };
const LALUnit lalSecondUnit        = {  0, { 0, 0, 1, 0, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };
const LALUnit lalAmpereUnit        = {  0, { 0, 0, 0, 1, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };
const LALUnit lalKelvinUnit        = {  0, { 0, 0, 0, 0, 1, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };
const LALUnit lalStrainUnit        = {  0, { 0, 0, 0, 0, 0, 1, 0}, { 0, 0, 0, 0, 0, 0, 0} };
const LALUnit lalADCCountUnit      = {  0, { 0, 0, 0, 0, 0, 0, 1}, { 0, 0, 0, 0, 0, 0, 0} };

/* Derived Mechanical Units */
const LALUnit lalHertzUnit         = {  0, { 0, 0,-1, 0, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };
const LALUnit lalNewtonUnit        = {  0, { 1, 1,-2, 0, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };
const LALUnit lalPascalUnit        = {  0, {-1, 1,-2, 0, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };
const LALUnit lalJouleUnit         = {  0, { 2, 1,-2, 0, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };
const LALUnit lalWattUnit          = {  0, { 2, 1,-3, 0, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };

/* Derived Electromagnetic Units */
const LALUnit lalCoulombUnit       = {  0, { 0, 0, 1, 1, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };
const LALUnit lalVoltUnit          = {  0, { 2, 1,-3,-1, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };
const LALUnit lalOhmUnit           = {  0, { 2, 1,-3,-2, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };
const LALUnit lalFaradUnit         = {  0, {-2,-1, 4, 2, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };
const LALUnit lalWeberUnit         = {  0, { 2, 1,-2,-1, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };
const LALUnit lalHenryUnit         = {  0, { 2, 1,-2,-2, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };
const LALUnit lalTeslaUnit         = {  0, { 0, 1,-2,-1, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };

/* Powers of Ten */
const LALUnit lalYottaUnit         = { 24, { 0, 0, 0, 0, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };
const LALUnit lalZettaUnit         = { 21, { 0, 0, 0, 0, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };
const LALUnit lalExaUnit           = { 18, { 0, 0, 0, 0, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };
const LALUnit lalPetaUnit          = { 15, { 0, 0, 0, 0, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };
const LALUnit lalTeraUnit          = { 12, { 0, 0, 0, 0, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };
const LALUnit lalGigaUnit          = {  9, { 0, 0, 0, 0, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };
const LALUnit lalMegaUnit          = {  6, { 0, 0, 0, 0, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };
const LALUnit lalKiloUnit          = {  3, { 0, 0, 0, 0, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };
const LALUnit lalHectoUnit         = {  2, { 0, 0, 0, 0, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };
const LALUnit lalDekaUnit          = {  1, { 0, 0, 0, 0, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };
const LALUnit lalDeciUnit          = { -1, { 0, 0, 0, 0, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };
const LALUnit lalCentiUnit         = { -2, { 0, 0, 0, 0, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };
const LALUnit lalMilliUnit         = { -3, { 0, 0, 0, 0, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };
const LALUnit lalMicroUnit         = { -6, { 0, 0, 0, 0, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };
const LALUnit lalNanoUnit          = { -9, { 0, 0, 0, 0, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };
const LALUnit lalPicoUnit          = {-12, { 0, 0, 0, 0, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };
const LALUnit lalFemtoUnit         = {-15, { 0, 0, 0, 0, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };
const LALUnit lalAttoUnit          = {-18, { 0, 0, 0, 0, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };
const LALUnit lalZeptoUnit         = {-21, { 0, 0, 0, 0, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };
const LALUnit lalYoctoUnit         = {-24, { 0, 0, 0, 0, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };

/* Convenient Scaled Units */
const LALUnit lalGramUnit          = { -3, { 0, 1, 0, 0, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };
const LALUnit lalAttoStrainUnit    = {-18, { 0, 0, 0, 0, 0, 1, 0}, { 0, 0, 0, 0, 0, 0, 0} };
const LALUnit lalPicoFaradUnit     = {-12, {-2,-1, 2, 2, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };

/* Static function to read a number into a character array */
/* returns 0 on success, 1 on failure */
/* leaves *charPtrPtr pointing to first non-digit */
static int readNumber( char temp[], const char **charPtrPtr )
{
  CHAR *tempPtr, *tempStopPtr;

  tempPtr = temp;
  /* make sure we don't fall off end of temporary array */
  tempStopPtr = temp + UNITDEFSC_TEMPSIZE;

  if ( ! isdigit(**charPtrPtr) ) return 1;

  do 
  {
    *tempPtr = **charPtrPtr;
    ++tempPtr, ++*charPtrPtr;
    if (tempPtr >= tempStopPtr) return 1;
  }
  while ( isdigit(**charPtrPtr) );
  *tempPtr = '\0';
  return 0;
}

/* Static function to read a string into a character array */
/* returns 0 on success, 1 on failure */
/* leaves *charPtrPtr pointing to first non-letter */
static int readString( char temp[UNITDEFSC_TEMPSIZE], const char **charPtrPtr )
{
  CHAR *tempPtr, *tempStopPtr;

  tempPtr = temp;
  /* make sure we don't fall off end of temporary array */
  tempStopPtr = temp + UNITDEFSC_TEMPSIZE; 

  if ( ! isalpha(**charPtrPtr) ) return 1;

  do 
  {
    *tempPtr = **charPtrPtr;
    ++tempPtr, ++*charPtrPtr;
    if (tempPtr >= tempStopPtr) return 1;
  }
  while ( isalpha(**charPtrPtr) );
  *tempPtr = '\0';
  return 0;
}

char * XLALUnitAsString( char *string, UINT4 length, const LALUnit *input )
{
  static const char *func = "XLALUnitAsString";
  UINT2        i;
  CHAR         temp[UNITDEFSC_TEMPSIZE];
  INT2         numer;
  CHAR         *charPtr, *charStopPtr;

  if ( ! string || ! input )
    XLAL_ERROR_NULL( func, XLAL_EFAULT );
  if ( ! length )
    XLAL_ERROR_NULL( func, XLAL_EBADLEN );

  charPtr = string;
  charStopPtr = string + length;
  /* Points one past the end of the array! */

  *charPtr = '\0';

  if (input->powerOfTen != 0)
  {
    sprintf(temp, "10^%d", input->powerOfTen);
    if ( charPtr + strlen(temp) >= charStopPtr) 
      XLAL_ERROR_NULL( func, XLAL_EBADLEN );
    strncpy(charPtr, temp, charStopPtr - charPtr);
    charPtr += strlen(temp);
  } /* if (input->powerOfTen != 0) */

  for (i=0; i<LALNumUnits; ++i)
  {
    numer = input->unitNumerator[i];
    if (numer != 0)
    {
      if (charPtr != string) 
      {
	*charPtr = ' ';
	++charPtr;
      } /* if (charPtr != output->data) */
      if (input->unitDenominatorMinusOne[i] == 0) 
      {
	if (numer == 1) 
	{
	  sprintf(temp, "%s", lalUnitName[i]);
	} /* if (numer == 1) */
	else
	{
	  sprintf(temp, "%s^%d", lalUnitName[i], numer);
	}
      } /* if (input->unitDenominatorMinusOne[i] == 0) */
      else {
	sprintf(temp, "%s^%d/%d", lalUnitName[i], numer,
		 input->unitDenominatorMinusOne[i] + 1);
      }
      if ( charPtr + strlen(temp) >= charStopPtr) 
        XLAL_ERROR_NULL( func, XLAL_EBADLEN );
      strncpy(charPtr, temp, charStopPtr - charPtr);
      charPtr += strlen(temp);
    } /* if (numer != 0) */
  }  /* for (i=0; i<LALNumUnits; ++i) */

  return string;
}

/* <lalVerbatim file="UnitDefsCP"> */
void 
LALUnitAsString( LALStatus *status,
		 CHARVector *output,
		 const LALUnit *input )
/* </lalVerbatim> */
{
  INITSTATUS( status, "LALUnitAsString", UNITDEFSC );
  /* ATTATCHSTATUSPTR (status); */

  ASSERT( input != NULL, status, UNITSH_ENULLPIN, UNITSH_MSGENULLPIN );

  ASSERT( output != NULL, status, UNITSH_ENULLPOUT, UNITSH_MSGENULLPOUT );

  ASSERT( output->data != NULL, status, UNITSH_ENULLPD,
	  UNITSH_MSGENULLPD );

  ASSERT( output->length > 0, status,
	  UNITSH_ESTRINGSIZE, UNITSH_MSGESTRINGSIZE );

  if ( ! XLALUnitAsString( output->data, output->length, input ) )
  {
    int code = xlalErrno;
    XLALClearErrno();
    switch ( code )
    {
      case XLAL_EFAULT: /* a NULL pointer was passed to the XLAL function */
        if ( ! input )
        {
          ABORT( status, UNITSH_ENULLPIN, UNITSH_MSGENULLPIN );
        }
        else /* must have been a NULL output data pointer */
        {
          ABORT( status, UNITSH_ENULLPIN, UNITSH_MSGENULLPIN );
        }
      default: /* otherwise must have almost overwritten the string */
        ABORT( status, UNITSH_ESTRINGSIZE, UNITSH_MSGESTRINGSIZE );
    }
  }

  RETURN(status);
}


LALUnit * XLALParseUnitString( LALUnit *output, const char *string )
{
  static const char *func = "XLALParseUnitString";
  UINT2        i;
  INT2         sign;
  CHAR         temp[20];
  const CHAR   *charPtr, *charStopPtr;
  int outputAllocated = 0;

  if ( ! output )
  {
    output = LALMalloc( sizeof( *output ) );
    outputAllocated = 1;
    if ( ! output )
      XLAL_ERROR_NULL( func, XLAL_ENOMEM );
  }

  /* Start with dimensionless (all zeros) and fill in from there */
  *output = lalDimensionlessUnit;

  /* If the string is NULL, it represents dimensionless */
  if ( ! string )
    return output;

  charPtr = string;
  charStopPtr = string + strlen(string); 

  /* Start with dimensionless (all zeros) and fill in from there */
  *output = lalDimensionlessUnit;
  
  /* If the string is empty, it represents dimensionless */
  if (charPtr == charStopPtr)
    return output;

  /* Look for power of ten; note LALUnitsAsString is set up to say
   * "10^1" rather than "10", so need not allow for missing '^'
   */
  if (*charPtr == '1' && *(charPtr+1) == '0' && *(charPtr+2) == '^') 
  {
    charPtr += 3; 
    /* now pointing at first digit of power of ten (or minus sign) */

    if ( *charPtr == '-'  ) 
    {
      sign = -1;
      ++charPtr;
    }
    else 
    {
      sign = 1;
    }
    
    /* read power of ten into temp[]; return value of 1 means failure */
    if ( readNumber( temp, &charPtr ) ) 
    {
      if ( outputAllocated )
        LALFree( output );
      XLAL_ERROR_NULL( func, XLAL_EFAILED );
    }
    /* charPtr now points to one after end of power of ten */
    
    output->powerOfTen = sign*atoi(temp);    

    /* If the power of ten was all there was, return */
    if (*charPtr == '\0')
      return output;

    if ( *charPtr != ' ') 
    {
      if ( outputAllocated )
        LALFree( output );
      XLAL_ERROR_NULL( func, XLAL_EFAILED );
    }
    
    ++charPtr;
  } /* if (*charPtr == '1' && *(charPtr+1) == '0' && *(charPtr+2) == '^') */

  /* charPtr now points to start of first unit */

  /* Read units and exponents, one unit per pass of the following do loop */
  do
  {
    /* read unit name into temp[]; return value of 1 means failure */
    if ( readString( temp, &charPtr ) )
    {
      if ( outputAllocated )
        LALFree( output );
      XLAL_ERROR_NULL( func, XLAL_EFAILED );
    }

    /* charPtr now points to one after end of unit name */

    /* find which unit name this matches */
    for (i=0; strcmp(temp,lalUnitName[i]); ++i)
    {
      if (i>=LALNumUnits) /* didn't find it */
      {
        if ( outputAllocated )
          LALFree( output );
        XLAL_ERROR_NULL( func, XLAL_EFAILED );
      }
    }

    /* Make sure we haven't already read in this unit */
    if ( output->unitNumerator[i] || output->unitDenominatorMinusOne[i] )
    {
      if ( outputAllocated )
        LALFree( output );
      XLAL_ERROR_NULL( func, XLAL_EFAILED );
    }

    if ( *charPtr == ' ' || *charPtr == '\0' )
    { /* We have only one power of the unit */
      output->unitNumerator[i] = 1;
    }
    else if ( *charPtr == '^' )
    {
      ++charPtr;
      /* now points to the first digit of the exponent, or minus sign */

      if ( *charPtr == '-'  ) 
      {
	sign = -1;
	++charPtr;
      }
      else 
      {
	sign = 1;
      }

      /* read exponent numerator into temp[];
	 return value of 1 means failure */
      if ( readNumber( temp, &charPtr ) )
      {
        if ( outputAllocated )
          LALFree( output );
        XLAL_ERROR_NULL( func, XLAL_EFAILED );
      }
      output->unitNumerator[i] = sign * atoi(temp);

      if ( *charPtr == '/' )
      {
	++charPtr;
	/* now points to first digit of denominator */

	/* read exponent denominator into temp[];
	   return value of 1 means failure */
	if ( readNumber( temp, &charPtr ) || temp[0] == '0')
        {
          if ( outputAllocated )
            LALFree( output );
          XLAL_ERROR_NULL( func, XLAL_EFAILED );
        }
	output->unitDenominatorMinusOne[i] = atoi(temp) - 1;
      } /* if ( *charPtr == '/' ) */
    } /* else if ( *charPtr == '^' ) */
    else 
    {
      if ( outputAllocated )
        LALFree( output );
      XLAL_ERROR_NULL( func, XLAL_EFAILED );
    }

    if ( *charPtr == ' ') ++charPtr;

  }
  while ( *charPtr != '\0' );
  
  return output;
}


/* <lalVerbatim file="UnitDefsCP"> */
void 
LALParseUnitString ( LALStatus *status,
		     LALUnit *output,
		     const CHARVector *input )
/* </lalVerbatim> */
{
  CHAR         *charPtr, *charStopPtr;

  INITSTATUS( status, "LALParseUnitString", UNITDEFSC );

  ASSERT( input != NULL, status, UNITSH_ENULLPIN, UNITSH_MSGENULLPIN );

  ASSERT( output != NULL, status, UNITSH_ENULLPOUT, UNITSH_MSGENULLPOUT );

  ASSERT( input->data != NULL, status, UNITSH_ENULLPD, UNITSH_MSGENULLPD );

  /* ensure that there's a '\0' within the input CHARVector */
  charPtr = input->data;
  charStopPtr = charPtr + strlen(input->data); 
  /* Should point to first '\0' in string */
  if (charStopPtr >= charPtr + input->length) 
  {
    ABORT( status, UNITSH_EPARSE, UNITSH_MSGEPARSE );
  }

  /* call the XLAL function */
  output = XLALParseUnitString( output, charPtr );
  if ( ! output ) /* there was a parse error */
  {
    XLALClearErrno();
    ABORT( status, UNITSH_EPARSE, UNITSH_MSGEPARSE );
  }

  RETURN(status);
}
