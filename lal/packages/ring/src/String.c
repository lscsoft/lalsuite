#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/String.h>

NRCSID( STRINGC, "$Id$" );


void
LALComputeTimeStringTemplate(
    LALStatus         *status,
    REAL4TimeSeries   *output,
    TimeTemplateInput *input,
    REAL4 power
    )

{  
  UINT4 n, i, d;
  REAL4 amplitude;
  REAL4 shift;

  INITSTATUS( status, "LALComputeTimeStringTemplate", STRINGC );

  ASSERT( input, status, STRINGH_ENULL, STRINGH_MSGENULL );
  ASSERT( output, status, STRINGH_ENULL, STRINGH_MSGENULL );
  ASSERT( output->data, status, STRINGH_ENULL, STRINGH_MSGENULL );

  if( output->data->length < input->duration)
  {
    ABORT( status, STRINGH_ETPSZ, STRINGH_MSGETPSZ);
  }
    
  n = output->data->length;
  d = input->duration;

  amplitude = 2.0 * pow( d/2.0 * output->deltaT , power+1);
  shift = - amplitude * pow(d/2 * output->deltaT,power);
  
  for (i = 0 ; i <= d/2 ; i++)
  {
    output->data->data[d/2+i] = amplitude 
      * pow( (output->deltaT * (i)), power) 
      + shift;
    output->data->data[d/2-i] = output->data->data[d/2+i];
  }
  for( i = d ; i < n ; i++)
  {
    output->data->data[i] = 0;
  }

  RETURN(status);
}

void LALComputeFreqStringTemplate(
    LALStatus         *status,
    COMPLEX8FrequencySeries *output,
    REAL4		power
    )
{ /* </lalVerbatim> */
  REAL4 freq = output->f0;
  UINT4 i;

  INITSTATUS( status, "LALComputeFreqStringTemplate", STRINGC );

  ASSERT( output, status, STRINGH_ENULL, STRINGH_MSGENULL );
  ASSERT( output->data, status, STRINGH_ENULL, STRINGH_MSGENULL );

  for ( i = 0; i < output->data->length; ++i )
  {
    if( 0==freq)
      output->data->data[i].re = 0;
    else
      output->data->data[i].re = (1.0/13.0)*pow(freq, power);
    output->data->data[i].im = 0;
    freq += output->deltaF;
  }
  RETURN(status);
}


static void LALMakeStringBank( 
    LALStatus *status,
    FreqTemplate *tmplt, 
    StringTemplateBankInput *input,
    UINT4 *finalcount );
static void LALMakeStringBank( 
    LALStatus *status,
    FreqTemplate *tmplt, 
    StringTemplateBankInput *input,
    UINT4 *finalcount )
{
  COMPLEX8FrequencySeries basetemplate;
  REAL4 intRatio = pow((1-input->maxMismatch),-2.0)-1.0;
  REAL4 base = 0;
  REAL4 top = 0;
  REAL4 norm = 0;
  UINT4 k=0;
  UINT4 count=0;
  UINT4 mink = floor((input->minFrequency - input->invSpectrum->f0)
		  / input->invSpectrum->deltaF);
  UINT4 maxk = floor((input->maxFrequency - input->invSpectrum->f0)
		  /input->invSpectrum->deltaF);

  
  
  INITSTATUS( status, "LALMakeStringBank", STRINGC );
  ATTATCHSTATUSPTR( status);
  CHECKSTATUSPTR(status);
  
  /*spectrum with no cutoff */
  strncpy( basetemplate.name, "freq domain string signal", 
      sizeof( basetemplate.name ) );
  basetemplate.f0   = input->invSpectrum->f0;
  basetemplate.deltaF = input->invSpectrum->deltaF;
  basetemplate.data = NULL; 
  LALCCreateVector( status->statusPtr , &basetemplate.data,
      input->invSpectrum->data->length);
  CHECKSTATUSPTR( status );
  
  LALComputeFreqStringTemplate( status->statusPtr, &basetemplate, input->freqPower );
  CHECKSTATUSPTR( status );
    
  while ( k < maxk ) /* while frequency < maxFrequency*/
  {
    top += input->invSpectrum->deltaF * input->invSpectrum->data->data[k]
     * basetemplate.data->data[k].re * basetemplate.data->data[k].re;
    if (k > mink && (!base || top/base > intRatio))
    {
      base += top;
      top = 0;
      if( tmplt)
      {
        tmplt[count].frequency = k * input->invSpectrum->deltaF;
        tmplt[count].normalization = norm;
      }
      ++count;
    }
    norm += input->invSpectrum->deltaF * basetemplate.data->data[k].re;
  k++;
  }
  if (finalcount)
   *finalcount = count;
  LALCDestroyVector(status->statusPtr , &basetemplate.data);
  DETATCHSTATUSPTR(status);
  RETURN( status);
}


/* <lalVerbatim file="RingCP"> */
void
LALCreateStringTemplateBank(
    LALStatus              *status,
    StringTemplateBank      **output,
    StringTemplateBankInput  *input
    )
{ /* </lalVerbatim> */
  StringTemplateBank *bank;

  INITSTATUS( status, "LALCreateStringTemplateBank", STRINGC );
  ATTATCHSTATUSPTR( status);

  ASSERT( input, status, STRINGH_ENULL, STRINGH_MSGENULL );
  ASSERT( output, status, STRINGH_ENULL, STRINGH_MSGENULL );
  ASSERT( ! *output, status, STRINGH_ENNUL, STRINGH_MSGENNUL );

  bank = *output = LALCalloc( 1, sizeof( **output ) );
  if ( ! bank )
  {
    ABORT( status, STRINGH_EALOC, STRINGH_MSGEALOC );
  }
  CHECKSTATUSPTR(status);
  LALMakeStringBank( status->statusPtr, NULL, input, &bank->numTmplt ); /*count number of templates needed*/
  
  bank->tmplt = LALMalloc( bank->numTmplt * sizeof( *bank->tmplt ) );
  if ( ! bank->tmplt )
  {
    LALFree( *output );
    *output = NULL;
    ABORT( status, STRINGH_EALOC, STRINGH_MSGEALOC );
  }

  LALMakeStringBank( status->statusPtr, bank->tmplt, input, NULL );
  DETATCHSTATUSPTR( status);
  RETURN( status );
}


/* <lalVerbatim file="RingCP"> */
void
LALDestroyStringTemplateBank(
    LALStatus         *status,
    StringTemplateBank **bank
    )
{ /* </lalVerbatim> */
  INITSTATUS( status, "LALDestroyStringTemplateBank", STRINGC );
  ASSERT( bank, status, STRINGH_ENULL, STRINGH_MSGENULL );
  ASSERT( *bank, status, STRINGH_ENULL, STRINGH_MSGENULL );
  LALFree( (*bank)->tmplt );
  LALFree( *bank );
  *bank = NULL;
  RETURN( status );
}
