import lalsimulation as lalsim
import lal

def read_liv_params(LALparams):
  logLeff = lalsim.SimInspiralWaveformParamsLookupNonGRLIVLogLambdaEff(LALparams) 
  signOfA = lalsim.SimInspiralWaveformParamsLookupNonGRLIVASign(LALparams)
  alpha = lalsim.SimInspiralWaveformParamsLookupNonGRLIVAlpha(LALparams)
  return logLeff, signOfA, alpha

def set_liv_pars(LALparams, logLeff, Asign, alpha):
  print("Setting LIV parameters")
  lalsim.SimInspiralWaveformParamsInsertNonGRLIVLogLambdaEff(LALparams, logLeff)
  lalsim.SimInspiralWaveformParamsInsertNonGRLIVASign(LALparams, Asign)
  lalsim.SimInspiralWaveformParamsInsertNonGRLIVAlpha(LALparams, alpha)
  return None 

def is_liv_enabled_by_default(LALparams):
  return lalsim.SimInspiralWaveformParamsLookupEnableLIV(LALparams)

def enable_liv(LALparams):
  lalsim.SimInspiralWaveformParamsInsertEnableLIV(LALparams, 1)
  return lalsim.SimInspiralWaveformParamsLookupEnableLIV(LALparams)

if __name__ == '__main__':
  ## Read default LIV parameters
  LALpars = lal.CreateDict()
  lleff, As, alph = read_liv_params(LALpars)
  print("By default: log10lambdaEff = %f, A_sign = %f, alpha = %f"%(lleff, As, alph))

  ## Setting LIV parameters to non-default values
  set_liv_pars(LALpars, 25.,-1.,1.5)

  ## Check if they have been correctly inserted
  lleff, As, alph = read_liv_params(LALpars)
  print("After inserting appropriate values: log10lambdaEff = %f, A_sign = %f, alpha = %f"%(lleff, As, alph))

  ## Check if LIV is enabled by default
  #ret = is_liv_enabled_by_default(LALpars)
  #print(ret)

  if is_liv_enabled_by_default(LALpars):
    print("LIV shouldn't be enabled by default - please check!")
  else:
    enable_liv(LALpars)
    print("LIV is enabled now!")
