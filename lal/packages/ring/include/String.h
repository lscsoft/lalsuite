#ifndef _STSTRING_H
#define _STSTRING_H

#include <lal/LALDatatypes.h>

NRCSID( STRINGH, "$Id$" );

#ifdef  __cplusplus
extern "C" {
#pragma } /** to match the previous brace **/
#endif

/**** <lalLaTeX>
 * \subsection*{Error conditions}
 **** </lalLaTeX> */
/**** <lalErrTable> */
#define STRINGH_ENULL 01
#define STRINGH_ENNUL 02
#define STRINGH_EALOC 04
#define STRINGH_ETPSZ 05
#define STRINGH_MSGENULL "Null pointer"
#define STRINGH_MSGENNUL "Non-null pointer"
#define STRINGH_MSGEALOC "Memory allocation error"
#define STRINGH_MSGETPSZ "Duration longer than template length"
/**** </lalErrTable> */

/**** <lalLaTeX>
 * 
 * \subsection*{Structures}
 * 
 * \subsubsection*{Type \texttt{StringTemplateInput}}
 * \idx[Type]{StringTemplateInput}
 *
 **** </lalLaTeX> */
typedef struct
tagFreqTemplate
{
  REAL4 frequency;
  REAL4 normalization;
}
FreqTemplate;

typedef struct
tagTimeTemplateInput
{
  REAL4 duration;
}
TimeTemplateInput;

typedef struct
tagStringTemplateBank
{
  UINT4              numTmplt;
  FreqTemplate       *tmplt;
}
StringTemplateBank;
typedef struct
tagStringTemplateBankInput
{
  REAL4 freqPower;
  REAL4 minFrequency;
  REAL4 maxFrequency;
  REAL4 maxMismatch;
  REAL4FrequencySeries *invSpectrum;
}
StringTemplateBankInput;

void
LALComputeTimeStringTemplate(
    LALStatus         *status,
    REAL4TimeSeries   *output,
    TimeTemplateInput *input,
    REAL4 power
    );

void
LALComputeFreqStringTemplate(
    LALStatus	  *status,
    COMPLEX8FrequencySeries *output,
    REAL4 		power
    );

void
LALCreateStringTemplateBank(
    LALStatus              *status,
    StringTemplateBank      **output,
    StringTemplateBankInput  *input
    );

void
LALDestroyStringTemplateBank(
    LALStatus         *status,
    StringTemplateBank **bank
    );


#ifdef  __cplusplus
#pragma { /** to match the next brace **/
}
#endif

#endif /* _STRING_H */
