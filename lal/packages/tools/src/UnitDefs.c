/************************************ <lalVerbatim file="UnitDefsCV">
Author: Whelan, J. T.
$Id$  
************************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Module \texttt{UnitDefs.c}}
\label{ss:UnitDefs.c}

Defines basic and derived SI units and a function to produce a text
string corresponding to a unit structure.

\subsubsection*{Prototypes}
\input{UnitDefsCP}
\index{\texttt{LALUnitAsString()}}

\subsubsection*{Description}

\texttt{LALUnitAsString()} converts the unit structure
\texttt{*input} into a text string which is stored in the character
vector \texttt{*output}.  Note that the resulting text string is
expressed solely in terms of the basic units (m, kg, s, A,
K, strain and counts), and is thus not necessarily the most
convenient way to check the units of a quantity.  A better method is
to construct a unit structure containing the expected units, then
compare that to the actual units using \texttt{LALUnitCompare()}.

\subsubsection*{Algorithm}
Moves through the unit structure, appending
the appropriate text to the string as it goes along.

\subsubsection*{Uses}
None.

\subsubsection*{Notes}

\subsubsection*{Predefined Units}

This file also defines a number of \texttt{constant} unit structures
(declared \texttt{extern} in \texttt{Units.h}).
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
Next, the named derived units in the SI\cite{Halliday:2001}:
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
\bibitem{Halliday:2001}
D.~Halliday, R.~Resnick, and J.~Walker, \textit{Fundamentals of
  Physics}.  (Wiley \& Sons, New York, 2001)
******************************************************* </lalLaTeX> */ 

#ifndef _LALSTDLIB_H
#include "LALStdlib.h"
#ifndef _LALSTDLIB_H
#define _LALSTDLIB_H
#endif
#endif

#ifndef _STRING_H
#include <string.h>
#ifndef _STRING_H
#define _STRING_H
#endif
#endif

#ifndef _UNITS_H
#include "Units.h"
#ifndef _UNITS_H
#define _UNITS_H
#endif
#endif

NRCSID( UNITDEFSC, "$Id$" );

/* To convert a units structure to a string repesentation, we need to
 * define the names of the basic units.
 */

enum { LALUnitNameSize = sizeof("strain") };
/* enum { LALUnitTextSize = sizeof("10^-32768 m^-32768/32767 kg^-32768/32767 "
				"s^-32768/32767 A^-32768/32767 " 
				"K^-32768/32767 strain^-32768/32767 "
				"count^-32768/32767") }; */

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

/* Basic Units */
const LALUnit lalMeterUnit       = {  0, { 1} };
const LALUnit lalKiloGramUnit    = {  0, { 0, 1} };
const LALUnit lalSecondUnit      = {  0, { 0, 0, 1} };
const LALUnit lalAmpereUnit      = {  0, { 0, 0, 0, 1} };
const LALUnit lalKelvinUnit      = {  0, { 0, 0, 0, 0, 1} };
const LALUnit lalStrainUnit      = {  0, { 0, 0, 0, 0, 0, 1} };
const LALUnit lalADCCountUnit    = {  0, { 0, 0, 0, 0, 0, 0, 1} };

/* Derived Mechanical Units */
const LALUnit lalHertzUnit       = {  0, { 0, 0,-1} };
const LALUnit lalNewtonUnit      = {  0, { 1, 1,-2} };
const LALUnit lalPascalUnit      = {  0, {-1, 1,-2} };
const LALUnit lalJouleUnit       = {  0, { 2, 1,-2} };
const LALUnit lalWattUnit        = {  0, { 2, 1,-3} };

/* Derived Electromagnetic Units */
const LALUnit lalCoulombUnit     = {  0, { 0, 0, 1, 1} };
const LALUnit lalVoltUnit        = {  0, { 2, 1,-3,-1} };
const LALUnit lalOhmUnit         = {  0, { 2, 1,-3,-2} };
const LALUnit lalFaradUnit       = {  0, {-2,-1, 4, 2} };
const LALUnit lalWeberUnit       = {  0, { 2, 1,-2,-1} };
const LALUnit lalHenryUnit       = {  0, { 2, 1,-2,-2} };
const LALUnit lalTeslaUnit       = {  0, { 0, 1,-2,-1} };

/* Powers of Ten */
const LALUnit lalYottaUnit       = { 24};
const LALUnit lalZettaUnit       = { 21};
const LALUnit lalExaUnit         = { 18};
const LALUnit lalPetaUnit        = { 15};
const LALUnit lalTeraUnit        = { 12};
const LALUnit lalGigaUnit        = {  9};
const LALUnit lalMegaUnit        = {  6};
const LALUnit lalKiloUnit        = {  3};
const LALUnit lalHectoUnit       = {  2};
const LALUnit lalDekaUnit        = {  1};
const LALUnit lalDeciUnit        = { -1};
const LALUnit lalCentiUnit       = { -2};
const LALUnit lalMilliUnit       = { -3};
const LALUnit lalMicroUnit       = { -6};
const LALUnit lalNanoUnit        = { -9};
const LALUnit lalPicoUnit        = {-12};
const LALUnit lalFemtoUnit       = {-15};
const LALUnit lalAttoUnit        = {-18};
const LALUnit lalZeptoUnit       = {-21};
const LALUnit lalYoctoUnit       = {-24};

/* Convenient Scaled Units */
const LALUnit lalGramUnit        = { -3, { 0, 1} };
const LALUnit lalAttoStrainUnit  = {-18, { 0, 0, 0, 0, 0, 1} };
const LALUnit lalPicoFaradUnit   = {-12, {-2,-1, 2, 2} };

/* <lalVerbatim file="UnitDefsCP"> */
void 
LALUnitAsString (LALStatus *status, CHARVector *output, const LALUnit *input)
/* </lalVerbatim> */
{
  UINT2        i;
  CHAR         temp[20];
  INT2         numer;
  CHAR         *charPtr, *charStopPtr;

  INITSTATUS( status, "LALUnitAsString", UNITDEFSC );
  /* ATTATCHSTATUSPTR (status); */

  ASSERT( input != NULL, status, UNITSH_ENULLPIN, UNITSH_MSGENULLPIN );

  ASSERT( output != NULL, status, UNITSH_ENULLPOUT, UNITSH_MSGENULLPOUT );

  ASSERT( output->data != NULL, status, UNITSH_ENULLPDOUT,
	  UNITSH_MSGENULLPDOUT );

  ASSERT( output->length > 0, status,
	  UNITSH_ESTRINGSIZE, UNITSH_MSGESTRINGSIZE );

  charPtr = output->data;
  charStopPtr = charPtr + output->length; 
  /* Points one past the end of the array! */

  *charPtr = '\0';

  if (input->powerOfTen != 0)
  {
    sprintf(temp, "10^%d", input->powerOfTen);
    if ( charPtr + strlen(temp) >= charStopPtr) 
    {
      ABORT( status, UNITSH_ESTRINGSIZE, UNITSH_MSGESTRINGSIZE );
    }
    strncpy(charPtr, temp, charStopPtr - charPtr);
    charPtr += strlen(temp);
  } /* if (input->powerOfTen != 0) */

  for (i=0; i<LALNumUnits; ++i)
  {
    numer = input->unitNumerator[i];
    if (numer != 0)
    {
      if (charPtr != output->data) 
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
      {
	ABORT( status, UNITSH_ESTRINGSIZE, UNITSH_MSGESTRINGSIZE );
      }
      strncpy(charPtr, temp, charStopPtr - charPtr);
      charPtr += strlen(temp);
    } /* if (numer != 0) */
  }  /* for (i=0; i<LALNumUnits; ++i) */

  /* printf("Units are:\"%s\"\n",output->data);*/

  /* DETATCHSTATUSPTR(status); */
  RETURN(status);
}
