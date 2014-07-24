import sys,os,math
from optparse import *
import random

sys.path.append('@PYTHONLIBDIR@')
import glue.pipeline

##############################################################################
usage = """
usage: %prog [options] 

"""

def determine_reduced_combos(combo):
  cp2 = glue.pipeline.DeepCopyableConfigParser()
  cp2.addsection('input')
  ifos = []
  temp = 0
  while temp < len(combo):
    ifos.append(combo[temp:temp+2])
    temp = temp+2
  if ('H1') in ifos:
    cp2.set('input','h1-data','')
  if ('H2') in ifos:
    cp2.set('input','h2-data','')
  if ('L1') in ifos:
    cp2.set('input','l1-data','')
  if ('V1') in ifos:
    cp2.set('input','v1-data','')
  if len(combo) > 5:
    cp2.set('input','two-ifos')
  if len(combo) > 7:
    cp2.set('input','three-ifos')
  redCombos = determine_ifo_combos(cp2)
  return redCombos

def determine_ifo_combos(cp):
  ifoCombos = []
  ifos = []
  if cp.has_option('input','h1-data'):
    ifos.append('H1')
  if cp.has_option('input','h2-data'):
    ifos.append('H2')
  if cp.has_option('input','l1-data'):
    ifos.append('L1')
  if cp.has_option('input','v1-data'):
    ifos.append('V1')
  if len(ifos) > 3 and cp.has_option('input','four-ifos'):
    ifoCombos.append(ifos[0] + ifos[1] + ifos[2] + ifos[3])
  if len(ifos) > 2 and cp.has_option('input','three-ifos'):
    if len(ifos) == 3:
      ifoCombos.append(ifos[0] + ifos[1] + ifos[2])
    elif len(ifos) == 4:
      ifoCombos.append(ifos[0] + ifos[1] + ifos[2])
      ifoCombos.append(ifos[0] + ifos[1] + ifos[3])
      ifoCombos.append(ifos[0] + ifos[2] + ifos[3])
      ifoCombos.append(ifos[1] + ifos[2] + ifos[3])
  if len(ifos) > 1 and cp.has_option('input','two-ifos'):
    if len(ifos) == 2:
      ifoCombos.append(ifos[0] + ifos[1])
    if len(ifos) == 3:
      ifoCombos.append(ifos[0] + ifos[1])
      ifoCombos.append(ifos[0] + ifos[2])
      ifoCombos.append(ifos[1] + ifos[2])
    if len(ifos) == 4:
      ifoCombos.append(ifos[0] + ifos[1])
      ifoCombos.append(ifos[0] + ifos[2])
      ifoCombos.append(ifos[0] + ifos[3])
      ifoCombos.append(ifos[1] + ifos[2])
      ifoCombos.append(ifos[1] + ifos[3])
      ifoCombos.append(ifos[2] + ifos[3])
  if cp.has_option('input','no-h1h2'):
    if 'H1H2' in ifoCombos:
      ifoCombos.remove('H1H2')
  return ifoCombos

def define_mass_characteristics(cp,object):
  massChars = {}
  if object[0:5] == 'mcomp':
    massVals = object.split('_')
    massChars['min_mass1']=str(float(massVals[1]))
    massChars['max_mass1']=str(float(massVals[2]))
    for name,value in cp.items('mcomp'):
      massChars[name] = value
  elif object[0:6] == 'mtotal':
    massVals = object.split('_')
    massChars['min_mtotal']=str(float(massVals[1]))
    massChars['max_mtotal']=str(float(massVals[2]))
    for name,value in cp.items('mtotal'):
      massChars[name] = value
  else:
    for name,value in cp.items(object):
      massChars[name] = value
  return massChars

def add_parent_child(parentJobs,childJob,dagmanParentChild):
  for parent in parentJobs:
    dagmanParentChild += 'PARENT ' + parent + ' CHILD ' + childJob + '\n'
  return dagmanParentChild

def add_ini_options(cp,subFile,prog):
  for opt,value in cp.items(prog):
    cmnd = '--' + opt + ' ' + value
    if cmnd[-1] != ' ':
      cmnd += ' '
    subFile += cmnd
  return subFile

def create_injcut_job(dagman,injInputFile,output,massChars):
  jobName = ''.join([random.choice('ABCDEFGHIJKLMNOPQR') for x in xrange(20)])
  dagman += 'JOB ' + jobName + ' upper_limit.injcut_' + massChars['type']\
            + '.sub\n'
  dagman += 'RETRY ' + jobName + ' 1 \n'
  dagman += 'VARS ' + jobName + ' macroinjfile="' + injInputFile + '"'
  if massChars['type'] == 'gaussian' or massChars['type'] == 'component':
    dagman += ' macrominmass1="' + massChars['min_mass1'] + '"' +\
              ' macrominmass2="' + massChars['min_mass2'] + '"' +\
              ' macromaxmass1="' + massChars['max_mass1'] + '"' +\
              ' macromaxmass2="' + massChars['max_mass2'] + '"'
  elif massChars['type'] == 'mtotal':
    dagman += ' macrominmasstot="' + massChars['min_mtotal'] + '"' +\
              ' macromaxmasstot="' + massChars['max_mtotal'] + '"'
  dagman += ' macrooutfile="' + output + '" \n'
  dagman += 'CATEGORY ' + jobName + ' injcut \n'
  return dagman,jobName

def create_injcut_subfile(cp,executable,logPath,priority,type):
  logFile = ''.join([random.choice('ABCDEFGHIJKLMNOPQR') for x in xrange(8)])
  subFile = 'universe = standard \n'
  subFile += 'executable = ' + executable + ' \n'
  subFile += 'arguments = --injection-file $(macroinjfile) '
  if type == 'gaussian' or type == 'component':
    subFile += '--mass-range-low $(macrominmass1) '
    subFile += '--mass2-range-low $(macrominmass2) '
    subFile += '--mass-range-high $(macromaxmass1) '
    subFile += '--mass2-range-high $(macromaxmass2) '
  elif type == 'mtotal':
    subFile += '--mass-range-low $(macrominmasstot) '
    subFile += '--mass-range-high $(macromaxmasstot) '
  subFile += '--output $(macrooutfile) '
  subFile = add_ini_options(cp,subFile,'injcut')
  subFile = add_ini_options(cp,subFile,'injcut-' + type)
  subFile += '\n'
  subFile += 'log = ' + logPath + '/' + logFile + '.tmp \n'
  subFile += 'error = logs/injcut-$(cluster)-$(process).err \n'
  subFile += 'output = logs/injcut-$(cluster)-$(process).out \n' 
  subFile += 'notification = never \n'
  subFile += 'priority = ' + priority + ' \n'
  subFile += 'queue 1 \n'
  submitFile = open('upper_limit.injcut_' + type + '.sub','w')
  submitFile.write(subFile)
  submitFile.close()

def create_inspinj_job(dagman,massChars,mass,sourceFile,type):
  jobName = ''.join([random.choice('ABCDEFGHIJKLMNOPQR') for x in xrange(20)])
  dagman += 'JOB ' + jobName + ' upper_limit.inspinj_'+type+'.sub\n'
  dagman += 'RETRY ' + jobName + ' 1 \n'
  dagman += 'VARS ' + jobName + ' macrosourcefile="' + sourceFile + '"'+ \
            ' macrominmass1="' + massChars['min_mass1'] + '"' +\
            ' macrominmass2="' + massChars['min_mass2'] + '"' +\
            ' macromaxmass1="' + massChars['max_mass1'] + '"' +\
            ' macromaxmass2="' + massChars['max_mass2'] + '"'
  if type == 'gaussian':
    dagman += ' macromeanmass1="' + massChars['mean_mass1'] + '"' +\
              ' macromeanmass2="' + massChars['mean_mass2'] + '"' +\
              ' macrostdmass1="' + massChars['std_mass1'] + '"' +\
              ' macrostdmass2="' + massChars['std_mass2'] + '"' 
  dagman += ' macromaxmtotal="' + massChars['max_mtotal'] + '"' +\
            ' macrominmtotal="' + massChars['min_mtotal'] + '"' +\
            ' macromassstr="' + mass.upper() + '"' +\
            ' macromaxdist="' + massChars['max_dist'] + '" \n'
  dagman += 'CATEGORY ' + jobName + ' inspinj \n'
  return dagman,jobName

def create_inspinj_subfile(cp,executable,logPath,priority,type):
  logFile = ''.join([random.choice('ABCDEFGHIJKLMNOPQR') for x in xrange(8)])
  subFile = 'universe = standard \n'
  subFile += 'executable = ' + executable + ' \n'
  subFile += 'arguments = --source-file ../$(macrosourcefile) '
  if type =='gaussian':
    subFile += '--user-tag GAUSSIANMASS$(macromassstr) '
  elif type == 'component':
    subFile += '--user-tag COMPONENT$(macromassstr) '
  elif type == 'mtotal':
    subFile += '--user-tag MTOTAL$(macromassstr) '
  subFile += '--min-mass1 $(macrominmass1) '
  subFile += '--min-mass2 $(macrominmass2) '
  subFile += '--max-mass1 $(macromaxmass1) '
  subFile += '--max-mass2 $(macromaxmass2) '
  if type == 'gaussian':
    subFile += '--stdev-mass1 $(macrostdmass1) '
    subFile += '--stdev-mass2 $(macrostdmass2) '
    subFile += '--mean-mass1 $(macromeanmass1) '
    subFile += '--mean-mass2 $(macromeanmass2) '
  subFile += '--max-mtotal $(macromaxmtotal) '
  subFile += '--min-mtotal $(macrominmtotal) '
  subFile += '--max-distance $(macromaxdist) '
  subFile = add_ini_options(cp,subFile,'inspinj')
  subFile = add_ini_options(cp,subFile,'inspinj-' + type)
  subFile += '\n'
  subFile += 'log = ' + logPath + '/' + logFile + '.tmp \n'
  subFile += 'error = ../logs/inspinj-$(cluster)-$(process).err \n'
  subFile += 'output = ../logs/inspinj-$(cluster)-$(process).out \n'
  subFile += 'initialdir = inspinj_files/ \n'
  subFile += 'notification = never \n'
  subFile += 'priority = ' + priority + ' \n'
  subFile += 'queue 1 \n'
  submitFile = open('upper_limit.inspinj_'+type+'.sub','w')
  submitFile.write(subFile)
  submitFile.close()

def create_coirefm_job(dagman,injFile,inputFiles,foundOutput,missedOutput,\
                       summaryOutput):
  jobName = ''.join([random.choice('ABCDEFGHIJKLMNOPQR') for x in xrange(20)])
  dagman += 'JOB ' + jobName + ' upper_limit.coirefm.sub\n'
  dagman += 'RETRY ' + jobName + ' 1 \n'
  dagman += 'VARS ' + jobName + ' macroglob="' + inputFiles + '"'+ \
            ' macrooutput="' + foundOutput + '"' +\
            ' macromissed="' + missedOutput + '"' +\
            ' macrosummary="' + summaryOutput + '"' +\
            ' macroinjfile="' + injFile + '" \n' 
  dagman += 'CATEGORY ' + jobName + ' coirefm \n'
  return dagman,jobName

def create_coirefm_subfile(cp,executable,logPath,priority):
  logFile = ''.join([random.choice('ABCDEFGHIJKLMNOPQR') for x in xrange(8)])
  subFile = 'universe = standard \n'
  subFile += 'executable = ' + executable + ' \n'
  subFile += 'arguments = --glob $(macroglob) '
  subFile += '--output $(macrooutput) '
  subFile += '--missed $(macromissed) '
  subFile += '--summary $(macrosummary) '
  subFile = add_ini_options(cp,subFile,'coire')
  subFile += '--injection-file $(macroinjfile) '
  subFile += '\n'
  subFile += 'log = ' + logPath + '/' + logFile + '.tmp \n'
  subFile += 'error = logs/coirefm-$(cluster)-$(process).err \n'
  subFile += 'output = logs/coirefm-$(cluster)-$(process).out \n'
  subFile += 'notification = never \n'
  subFile += 'priority = ' + priority + ' \n'
  subFile += 'queue 1 \n'
  submitFile = open('upper_limit.coirefm.sub','w')
  submitFile.write(subFile)
  submitFile.close()

def create_numgalaxies_job(dagman, inputZeroFiles,inputSlideFiles,\
       foundInjections,missedInjections,sourceFile,populationFile,\
       outputDir,h2CombinedDist,type,Cals):
  if type[0:5] == 'mcomp':
    massVals = type.split('_')
    minMass1=float(massVals[1])
    maxMass1=float(massVals[2])
    type = 'component'
  elif type[0:6] == 'mtotal':
    massVals = type.split('_')
    minMtotal=float(massVals[1])
    maxMtotal=float(massVals[2])
    type = 'mtotal'
  else:
    type = 'gaussian'
  jobName = ''.join([random.choice('ABCDEFGHIJKLMNOPQR') for x in xrange(20)])
  dagman += 'JOB ' + jobName + ' upper_limit.numgalaxies_'+type+'.sub\n'
  dagman += 'RETRY ' + jobName + ' 1 \n'
  dagman += 'VARS ' + jobName + ' macrozerofiles="' + inputZeroFiles + '"'+ \
            ' macroslidefiles="' + inputSlideFiles + '"' +\
            ' macrofoundinj="' + foundInjections + '"' +\
            ' macromissedinj="' + missedInjections + '"' +\
            ' macrosourcefile="' + sourceFile + '"' +\
            ' macrooutputdir="' + outputDir + '"' +\
            ' macropopfile="' + populationFile + '" '
  if type == 'component':
    dagman += 'macrominmass =" ' + str(minMass1) + '" '
    dagman += 'macromaxmass =" ' + str(maxMass1) + '" '
    dagman += 'macrodm =" ' + str(maxMass1 - minMass1) + '" '
  if type == 'mtotal':
    dagman += 'macrominmass ="' + str(minMtotal) + '" '
    dagman += 'macromaxmass ="' + str(maxMtotal) + '" '
    dagman += 'macrodm ="' + str(maxMtotal - minMtotal) + '" '
  if h2CombinedDist:
    dagman += 'macroh2distoption ="--h2-combined-dist" '
  else:
    dagman += 'macroh2distoption =" " '
  dagman += 'macrohcal = "' + Cals['H'] + '" '
  dagman += 'macrolcal = "' + Cals['L'] + '" '
  dagman += 'macrohdccal = "' + Cals['H_DC'] + '" '
  dagman += 'macroldccal = "' + Cals['L_DC'] + '" '
  dagman += '\n'
  dagman += 'CATEGORY ' + jobName + ' numgalaxies \n'
  return dagman,jobName

def create_numgalaxies_subfile(cp,executable,logPath,priority,type,userTag):
  logFile = ''.join([random.choice('ABCDEFGHIJKLMNOPQR') for x in xrange(8)])
  subFile = 'universe = vanilla \n'
  subFile += 'executable = ' + executable + ' \n'
  subFile += 'arguments = --zero-glob $(macrozerofiles) '
  subFile += '--found-glob $(macrofoundinj) '
  subFile += '--slide-glob $(macroslidefiles) '
  subFile += '--missed-glob $(macromissedinj) '
  subFile += '--source-file $(macrosourcefile) '
  subFile += '--population-glob $(macropopfile) '
  subFile += '--figure-name ' + userTag + ' '
  subFile += '$(macroh2distoption) '
  if not type == 'gaussian':
    subFile += '--cut-inj-by-mass 1 '
    subFile += '--m-low $(macrominmass) '
    subFile += '--m-high $(macromaxmass) '
    subFile += '--m-dm $(macrodm) '
  subFile += '--h-calibration $(macrohcal) '
  subFile += '--l-calibration $(macrolcal) '
  subFile += '--h-dc-calibration $(macrohdccal) '
  subFile += '--l-dc-calibration $(macroldccal) '
  subFile = add_ini_options(cp,subFile,'plotnumgalaxies')
  subFile = add_ini_options(cp,subFile,'plotnumgalaxies-' + type)
  subFile += '\n'
  subFile += 'log = ' + logPath + '/' + logFile + '.tmp \n'
  subFile += 'error = ../../../../logs/numgalaxies-$(cluster)-$(process).err \n'
  subFile += 'output = plotnumgalaxies.out \n'
  subFile += 'initialdir = $(macrooutputdir) \n'
  subFile += 'notification = never \n'
  subFile += 'getenv = true \n'
  subFile += 'priority = ' + priority + ' \n'
  subFile += 'queue 1 \n'
  submitFile = open('upper_limit.numgalaxies_' + type +'.sub','w')
  submitFile.write(subFile)
  submitFile.close()

def create_computeposterior_job(dagman, sourceChars,outputDir,\
       galaxiesFile,timeAnalyzedFile):
  jobName = ''.join([random.choice('ABCDEFGHIJKLMNOPQR') for x in xrange(20)])
  dagman += 'JOB ' + jobName + ' upper_limit.computeposterior.sub\n'
  dagman += 'RETRY ' + jobName + ' 1 \n'
  dagman += 'VARS ' + jobName + ' macrooutputdir="' + outputDir + '"'+ \
            ' macronumgalinput="' + galaxiesFile + '"' +\
            ' macrotimeanalyzed="' + timeAnalyzedFile + '"'
  if sourceChars['type'] == 'gaussian':
    dagman += ' macromasstype="gaussian" \n'
  elif sourceChars['type'] == 'mtotal':
    dagman += ' macromasstype="totalmass" \n'
  elif sourceChars['type'] == 'component':
    dagman += ' macromasstype="componentmass" \n'
  dagman += 'CATEGORY ' + jobName + ' numgalaxies \n'
  return dagman,jobName

def create_computeposterior_subfile(cp,executable,logPath,priority,userTag):
  logFile = ''.join([random.choice('ABCDEFGHIJKLMNOPQR') for x in xrange(8)])
  subFile = 'universe = vanilla \n'
  subFile += 'executable = ' + executable + ' \n'
  subFile += 'arguments = --galaxies-file $(macronumgalinput) '
  subFile += '--time-analyzed-file $(macrotimeanalyzed) '
  subFile += '--mass-region $(macromasstype) '
  subFile += '--figure-name ' + userTag + ' '
  subFile = add_ini_options(cp,subFile,'compute_posterior')
  subFile += '\n'
  subFile += 'log = ' + logPath + '/' + logFile + '.tmp \n'
  subFile += 'error = ../../../../logs/computeposterior-$(cluster)-$(process).err \n'
  subFile += 'output = computeposterior.out \n'
  subFile += 'initialdir = $(macrooutputdir) \n'
  subFile += 'notification = never \n'
  subFile += 'getenv = true \n'
  subFile += 'priority = ' + priority + ' \n'
  subFile += 'queue 1 \n'
  submitFile = open('upper_limit.computeposterior.sub','w')
  submitFile.write(subFile)
  submitFile.close()

def create_plotulvsmass_job(dagman,computepostglob,massRegion,figureName,\
                            initialDir):
  jobName = ''.join([random.choice('ABCDEFGHIJKLMNOPQR') for x in xrange(20)])
  dagman += 'JOB ' + jobName + ' upper_limit.plotulvsmass.sub\n'
  dagman += 'RETRY ' + jobName + ' 1 \n'
  dagman += 'VARS ' + jobName + ' macromassregion="' + massRegion + '"' +\
            ' macrocomputeglob="../' + computepostGlob + '"'+ \
            ' macromassregion="' + massRegion + '"' +\
            ' macrofigurename="' + figureName + '"' +\
            ' macroinitdir="' + initialDir + '" \n'
  dagman += 'CATEGORY ' + jobName + ' corse \n'
  return dagman,jobName

def create_plotulvsmass_subfile(cp,executable,logPath,priority):
  logFile = ''.join([random.choice('ABCDEFGHIJKLMNOPQR') for x in xrange(8)])
  subFile = 'universe = vanilla \n'
  subFile += 'executable = ' + executable + ' \n'
  subFile += 'arguments = '
  subFile += '--computepost-glob $(macrocomputeglob) '
  subFile += '--mass-region $(macromassregion) '
  subFile += '--figure-name $(macrofigurename) '
  subFile = add_ini_options(cp,subFile,'plotulvsmass')
  subFile += '\n'
  subFile += 'log = ' + logPath + '/' + logFile + '.tmp \n'
  subFile += 'error = ../logs/plotulvsmass-$(cluster)-$(process).err \n'
  subFile += 'output = plotulvsmass-$(macrofigurename).out \n'
  subFile += 'notification = never \n'
  subFile += 'initialdir = $(macroinitdir) \n'
  subFile += 'priority = ' + priority + ' \n'
  subFile += 'getenv = true \n'
  subFile += 'queue 1 \n'
  submitFile = open('upper_limit.plotulvsmass.sub','w')
  submitFile.write(subFile)
  submitFile.close()

def create_combineposterior_job(dagman,posteriorFiles,
                                initialDir,figureName,relDir,minMass,\
                                maxMass,massType):
  jobName = ''.join([random.choice('ABCDEFGHIJKLMNOPQR') for x in xrange(20)])
  dagman += 'JOB ' + jobName + ' upper_limit.combineposterior.sub\n'
  dagman += 'RETRY ' + jobName + ' 1 \n'
  dagman += 'VARS ' + jobName + ' macrofigurename="' + figureName + '"' +\
            ' macroinitdir="' + initialDir + '"'
  macroaddposts = '"'
  for file,name in posteriorFiles:
    macroaddposts += ' --add-posterior ' + relDir + file + ',' + name 
  macroaddposts += '"'
  dagman += ' macroaddposts =' + macroaddposts
  dagman += ' macrominmass ="' + minMass + '"'
  dagman += ' macromaxmass ="' + maxMass + '"'
  dagman += ' macroreldir ="' + relDir + '"'
  dagman += ' macromasstype ="' + massType.replace('_',' ') + '"'
  dagman += ' \n'
  dagman += 'CATEGORY ' + jobName + ' corse \n'
  return dagman,jobName

def create_combineposterior_subfile(cp,executable,logPath,priority):
  logFile = ''.join([random.choice('ABCDEFGHIJKLMNOPQR') for x in xrange(8)])
  subFile = 'universe = vanilla \n'
  subFile += 'executable = ' + executable + ' \n'
  subFile += 'arguments = --figure-name $(macrofigurename) '
  subFile += '--min-mass $(macrominmass) '
  subFile += '--max-mass $(macromaxmass) '
  subFile += '--plot-title $(macromasstype) '
  subFile += '$(macroaddposts) '
  subFile = add_ini_options(cp,subFile,'combine_posterior')
  subFile += '\n'
  subFile += 'log = ' + logPath + '/' + logFile + '.tmp \n'
  subFile += 'error =$(macroreldir)/logs/combineposterior-$(cluster)-$(process).err \n'
  subFile += 'output = combineposterior-$(macrofigurename).out \n'
  subFile += 'getenv = true \n'
  subFile += 'notification = never \n'
  subFile += 'initialdir = $(macroinitdir) \n'
  subFile += 'priority = ' + priority + ' \n'
  subFile += 'queue 1 \n'
  submitFile = open('upper_limit.combineposterior.sub','w')
  submitFile.write(subFile)
  submitFile.close()

def parse_command_line():
  """
  Parser function dedicated
  """
  parser = OptionParser( usage=usage, version="%prog CVS $Id$ " )
  parser.add_option("-f", "--config-file",action="store",type="string",\
      metavar=" FILE", help="use configuration file FILE")
  parser.add_option("-I", "--skip-injcut",action="store_true",default=False,\
      help="Do not run injcut jobs (assume they are already run).")
  parser.add_option("-i", "--skip-inspinj",action="store_true",default=False,\
      help="Do not run inspinj jobs (assume they are already run).")
  parser.add_option("-C", "--skip-coire",action="store_true",default=False,\
      help="Do not run coire jobs (assume they are already run).")
  parser.add_option("-n", "--skip-numgalaxies",action="store_true",\
      default=False,\
      help="Do not run plotnumgalaxies jobs (assume they are already run).")
  parser.add_option("-c", "--skip-compute-posterior",action="store_true",\
      default=False,\
      help="Do not run lalapps_compute_posterior jobs (assume already run).")
  parser.add_option("-u", "--skip-ulvsmass",action="store_true",default=False,\
      help="Do not run plotulvsmass jobs on individual ifo combos.") 
  parser.add_option("-p", "--skip-combine-posterior",action="store_true",\
      default=False,\
      help="Do not run pylal_combine_posteriors or generate combined plotulvsmass plots.")
  parser.add_option("-s", "--skip-spin",action="store_true",default=False,\
      help="Do not generate upper limit for spinning sources.")
  parser.add_option("-N", "--skip-nonspin",action="store_true",default=False,\
      help="Do not generate upper limits for non-spinning sources.")
  parser.add_option("-Z", "--combine-only-past-results",action="store_true",\
      default=False,\
      help="Use this flag to generate results from past results only.")
  parser.add_option("-g", "--skip-gaussian",action="store_true",default=False,\
      help="Do not generate upper limit for gaussian sources.")
  parser.add_option("-t", "--skip-total-mass",action="store_true",\
      default=False,\
      help="Do not generate upper limits as a function of total mass.")
  parser.add_option("-m", "--skip-component-mass",action="store_true",\
      default=False,\
      help="Do not generate upper limits as a function of component mass.")


  # options related to input and output
#  parser.add_option("","--output-file",action="store",type="string",\
#      default=None,metavar=" FILE",\
#      help="file to write to")

#  parser.add_option("","--num-slides",action="store",type="int",default=0,\
#      metavar=" NUM_SLIDES",help="number of time slides performed" )

  (options,args) = parser.parse_args()


  return options, sys.argv[1:]

# ============================================================================
# -- get command line arguments
opts, args = parse_command_line()
if not opts.config_file:
  print >> sys.stderr , 'You must specify a config file'
  sys.exit(1)

###################################

cp = glue.pipeline.DeepCopyableConfigParser()
cp.read(opts.config_file)

# Read ini file and determine how to run the code

runInjcut = not opts.skip_injcut
runInspinj = not opts.skip_inspinj
runCoireFM = not opts.skip_coire
runNumgalaxies = not opts.skip_numgalaxies
runComputeposterior = not opts.skip_compute_posterior
runPlotulvsmass = not opts.skip_ulvsmass
runCombinePosterior = not opts.skip_combine_posterior
runSpin = not opts.skip_spin
runNonSpin = not opts.skip_nonspin
runGaussian = not opts.skip_gaussian
runMtotal = not opts.skip_total_mass
runComponent = not opts.skip_component_mass
combineOnlyPast = opts.combine_only_past_results
massCategories = []
sourceTypes = []
gaussianSources = []
componentSources = []
mtotalSources = []
mrangeSources = []
numCategories = {}
H_cal = {}
H_calDC = {}
L_cal = {}
L_calDC = {}
ifoCombos = determine_ifo_combos(cp)
for item in cp.options('ifar-mass-categories'):
  massCategories.append(item)
if runGaussian:
  for item in cp.options('gaussian-types'):
    sourceTypes.append(item)
    gaussianSources.append(item)
if runComponent:
  for item in cp.options('total-mass-ranges'):
    sourceTypes.append(item)
    mtotalSources.append(item)
    mrangeSources.append(item)
if runMtotal:
  for item in cp.options('component-mass-ranges'):
    sourceTypes.append(item)
    componentSources.append(item)
    mrangeSources.append(item)
for item,value in cp.items('H-cal'):
  H_cal[item.upper()] = value
for item,value in cp.items('H-cal-dc'):
  H_calDC[item.upper()] = value
for item,value in cp.items('L-cal'):
  L_cal[item.upper()] = value
for item,value in cp.items('L-cal-dc'):
  L_calDC[item.upper()] = value

nonSpinInjRuns = []
injSeeds = {}
spinInjRuns = []
if runNonSpin:
  temp = cp.items('non-spin-injection-runs')
  for dummy,string in temp:
    name,seed=string.split(',')
    nonSpinInjRuns.append(name)
    injSeeds[name] = seed
  nonspin = nonSpinInjRuns[0]
if runSpin:
  temp = cp.items('spin-injection-runs')
  for dummy,string in temp:
    name,seed=string.split(',')
    spinInjRuns.append(name)
    injSeeds[name] = seed
  spin = spinInjRuns[0]
ihopeDir = cp.get('main','ihope-directory')
ifarDir = cp.get('main','ifar-directory')
gpsStartTime = cp.get('main','gps-start-time')
gpsEndTime = cp.get('main','gps-end-time')
gpsDuration = str(int(gpsEndTime) - int(gpsStartTime))
logPath = cp.get('main','log-path')
priority = cp.get('main','dagman-job-priority')
userTag = cp.get('main','user-tag')
dagman = ''
dagmanParentChild = ''
dagmanMaxJobs = ''

spinTypes = []
if runNonSpin: 
  spinTypes.append('nonspin')
if runSpin:
  spinTypes.append('spin')

# Add old posterior files to ini file from sets
nameNumber = 1000
for dummyTag,info in cp.items('past-posterior-result-sets'):
  pastTag,setDir = info.split(',')
  for type in spinTypes:
    for item in sourceTypes:
      tempTag = 'post' + str(nameNumber)
      nameNumber += 1
      section = 'past-posteriors-' + item + '-' + type
      file = setDir + '/' + type + '/' + item + '/' + pastTag + \
          '_posterior-pdf.txt'
      cp.set(section,tempTag,','.join([pastTag,file])) 

if not os.path.isdir('logs'):
  os.mkdir('logs')

# Run lalapps_injcut jobs
injcutJobs = {}
executable =cp.get('executables','lalapps_injcut')
if runSpin or runNonSpin:
  injcutOutputDirs = {}
  if not os.path.isdir('injcut_files'):
    os.mkdir('injcut_files')
  if runSpin:
    for inj in spinInjRuns:
      injcutOutputDirs[inj]='injcut_files/spin/'
      if not os.path.isdir('injcut_files/spin'):
        os.mkdir('injcut_files/spin')
  if runNonSpin:
    for inj in nonSpinInjRuns:
      injcutOutputDirs[inj]='injcut_files/nonspin/'
      if not os.path.isdir('injcut_files/nonspin'):
        os.mkdir('injcut_files/nonspin')
else:
  runInjcut = False

if runInjcut:
  for source in sourceTypes:
    sourceChars = define_mass_characteristics(cp,source)
    for inj in injcutOutputDirs.keys():
      injInputFile = ihopeDir + '/' + inj + '/HL-INJECTIONS_' + \
          injSeeds[inj] + '_' + inj.upper() + '-' + gpsStartTime + '-' + \
          gpsDuration + '.xml'
      output = injcutOutputDirs[inj] + 'HL-INJECTIONS_' + source.upper()\
               + '_' + inj.upper()\
               + '-' + gpsStartTime + '-' + gpsEndTime + '.xml'
      dagman,injcutJobs[source + inj] = create_injcut_job(dagman, \
               injInputFile,output,sourceChars)
  create_injcut_subfile(cp,executable,logPath,priority,'gaussian')
  create_injcut_subfile(cp,executable,logPath,priority,'component')
  create_injcut_subfile(cp,executable,logPath,priority,'mtotal')
 
# Run lalapps_inspinj jobs
if runInspinj:
  inspinjJobs = {}
  executable =cp.get('executables','lalapps_inspinj')
  if not os.path.isdir('inspinj_files'):
    os.mkdir('inspinj_files')
  for source in sourceTypes:
    sourceChars = define_mass_characteristics(cp,source)
    if sourceChars['type'] == 'gaussian':
      sourceFile = cp.get('inspinj-source-files',source)
      dagman,inspinjJobs[source] = create_inspinj_job(dagman, sourceChars,\
          source,sourceFile,'gaussian')
    if sourceChars['type'] == 'component':
      sourceFile = cp.get('inspinj-source-files','component')
      dagman,inspinjJobs[source] = create_inspinj_job(dagman, sourceChars,\
          source,sourceFile,'component')
    if sourceChars['type'] == 'mtotal':
      sourceFile = cp.get('inspinj-source-files','mtotal')
      dagman,inspinjJobs[source] = create_inspinj_job(dagman, sourceChars,\
          source,sourceFile,'mtotal')
  create_inspinj_subfile(cp,executable,logPath,priority,'gaussian')
  create_inspinj_subfile(cp,executable,logPath,priority,'component')   
  create_inspinj_subfile(cp,executable,logPath,priority,'mtotal')
  
# Run lalapps_coire to determine missed/found
coireFMJobs = {}
executable =cp.get('executables','lalapps_coire')
if runSpin or runNonSpin:
  coireFMOutputDirs = {}
  if not os.path.isdir('coire_found_missed_files'):
    os.mkdir('coire_found_missed_files')
  if runSpin:
    for inj in spinInjRuns:
      coireFMOutputDirs[inj]='coire_found_missed_files/spin/'
  if runNonSpin:
    for inj in nonSpinInjRuns:
      coireFMOutputDirs[inj]='coire_found_missed_files/nonspin/'
else:
  runCoireFM = False

if runCoireFM:
  for source in sourceTypes:
    for combo in ifoCombos:
      for inj in coireFMOutputDirs.keys():
        outputDir = coireFMOutputDirs[inj] + source
        if not os.path.isdir(outputDir):
          os.makedirs(outputDir)
        injFile = injcutOutputDirs[inj] + '/HL-INJECTIONS_' + source.upper()\
               + '_' + inj.upper() + '-' + gpsStartTime + '-' + \
               gpsEndTime + '.xml'
        inputFiles = ifarDir + '/combined_ifar_files/' + inj.upper() + \
                  '/' + combo + '-CORSE_*-' + inj.upper() + \
                  '_COMBINED_IFAR_CAT_3-' + gpsStartTime + '-' + \
                  gpsEndTime + '.xml.gz'
#        inputFiles = ifarDir + '/corse_all_data_files/' + inj.upper() + \
#                 '/' + combo + '_*-CORSE_' + inj.upper() + \
#                 '_*_CAT_3-' + gpsStartTime + '-' + gpsDuration + '.xml*'
        foundOutput = outputDir + '/' + combo +\
               '-CORSE_' + inj.upper() + 'FOUND_CAT_3-' + gpsStartTime +\
               '-' + gpsDuration + '.xml.gz'
        missedOutput = outputDir + '/' + combo + \
               '-CORSE_' + inj.upper() + 'MISSED_CAT_3-' + gpsStartTime \
               + '-' + gpsDuration + '.xml.gz'
        summaryOutput = outputDir + '/' + combo +\
               '-CORSE_' + inj.upper() + 'FOUND_CAT_3-' + gpsStartTime +\
               '-' + gpsDuration + '.txt'
        dagman,coireFMJobs[source+combo+inj] = create_coirefm_job(dagman, \
               injFile, inputFiles, foundOutput,missedOutput,summaryOutput)
        if runInjcut:
          parentJob = [injcutJobs[source+inj]]
          dagmanParentChild = add_parent_child(parentJob,\
               coireFMJobs[source+combo+inj],dagmanParentChild)
  create_coirefm_subfile(cp,executable,logPath,priority)
      
# Run plotnumgalaxies
numgalaxiesJobs = {}
executable = cp.get('executables','plotnumgalaxies')
if runSpin or runNonSpin:
  numgalaxiesOutputDirs = {}
  if not os.path.isdir('plotnumgalaxies_files'):
    os.mkdir('plotnumgalaxies_files')
  for type in spinTypes:
    for source in sourceTypes:
      for combo in ifoCombos:
        if not os.path.isdir('plotnumgalaxies_files/' + type + '/' + \
               source + '/' + combo):
          os.makedirs('plotnumgalaxies_files/' + type + '/' + source +\
                 '/' + combo)
        numgalaxiesOutputDirs[source+combo+type]=\
            'plotnumgalaxies_files/'+type+'/' + source + '/' + combo + '/'
else:
  runNumgalaxies = False

if runNumgalaxies:
  for type in spinTypes:
    for source in sourceTypes:
      sourceChars = define_mass_characteristics(cp,source)
      for combo in ifoCombos:
         outputDir = numgalaxiesOutputDirs[source+combo+type]
         inputZeroFiles = ifarDir + '/combined_ifar_files/exclude_play/' + \
             combo + '-CORSE_ALL_MASSES-exclude_play_COMBINED_IFAR_CAT_3-' +\
             gpsStartTime + '-' + gpsEndTime + '.xml.gz'
#         inputZeroFiles = ifarDir+'/corse_all_data_files/all_data/' + combo \
#             + '*-CORSE_ALL_DATA_*_CAT_3-*.xml.gz'
         inputSlideFiles = ifarDir+'/combined_ifar_files/all_data_slide/' + \
             combo+'-CORSE_ALL_MASSES-all_data_slide_COMBINED_IFAR_CAT_3-' +\
             gpsStartTime + '-' + gpsEndTime + '.xml.gz'
#         inputSlideFiles = ifarDir+'/corse_all_data_files/all_data_slide/' + \
#             combo + '*-CORSE_SLIDE_ALL_DATA_*_CAT_3-*.xml.gz'
         foundInjections = '../../../../' + coireFMOutputDirs[eval(type)] + \
             source + '/' + combo + '-CORSE_*FOUND_CAT_3-' \
             + gpsStartTime + '-' + gpsDuration + '.xml.gz'
         missedInjections = '../../../../' + coireFMOutputDirs[eval(type)] + \
             source + '/' + combo + '-CORSE_*MISSED_CAT_3-'+ gpsStartTime + \
             '-' + gpsDuration + '.xml.gz'
         if sourceChars['type']=='gaussian':
           sourceFile = cp.get('inspinj-source-files',source)
           populationFile = '../../../../inspinj_files/HL-INJECTIONS_1' + \
               '_GAUSSIANMASS' + source.upper() + '-793130413-2548800.xml.gz'
         elif sourceChars['type']=='component':
           sourceFile = cp.get('inspinj-source-files','component')
           populationFile = '../../../../inspinj_files/HL-INJECTIONS_1' + \
               '_COMPONENT' + source.upper() + '-793130413-2548800.xml.gz'
         elif sourceChars['type']=='mtotal':
           sourceFile = cp.get('inspinj-source-files','mtotal')
           populationFile = '../../../../inspinj_files/HL-INJECTIONS_1' + \
               '_MTOTAL' + source.upper() + '-793130413-2548800.xml.gz'
         if not sourceFile[0] == '/':
           sourceFile = '../../../../' + sourceFile
         if cp.has_option('H2-dist-option',combo):
           h2CombinedDist = True
         else:
           h2CombinedDist = False
         Cals = {}
         Cals['H'] = H_cal[combo]
         Cals['L'] = L_cal[combo]
         Cals['H_DC'] = H_calDC[combo]
         Cals['L_DC'] = L_calDC[combo]
         dagman,numgalaxiesJobs[source + combo + type] = \
             create_numgalaxies_job(dagman, \
             inputZeroFiles,inputSlideFiles,foundInjections,\
             missedInjections,sourceFile,populationFile,\
             outputDir,h2CombinedDist,source,Cals)
         parentJobs = []
         if runInjcut:
           for key in injcutJobs.keys():
             if key[0:len(source)] == source:
               parentJobs.append(injcutJobs[key])
         if runInspinj:
           parentJobs.append(inspinjJobs[source])
         if runCoireFM:
           for key in coireFMJobs.keys():
             if key[0:len(source+combo)] == source+combo:
               parentJobs.append(coireFMJobs[key])
         if parentJobs:
           dagmanParentChild = add_parent_child(parentJobs,\
               numgalaxiesJobs[source+combo+type],dagmanParentChild)
  create_numgalaxies_subfile(cp,executable,logPath,priority,'mtotal',userTag)
  create_numgalaxies_subfile(cp,executable,logPath,priority,'gaussian',userTag)
  create_numgalaxies_subfile(cp,executable,logPath,priority,'component',userTag)                 

computeposteriorJobs = {}
executable = cp.get('executables','lalapps_compute_posterior')

if not (runSpin or runNonSpin):
  runComputeposterior = False

if runComputeposterior:
  for type in spinTypes:
    for source in sourceTypes:
      sourceChars = define_mass_characteristics(cp,source)
      for combo in ifoCombos:
        outputDir = numgalaxiesOutputDirs[source+combo+type]
        galaxiesFile = '../../../../' + \
            numgalaxiesOutputDirs[source+combo+type] + '/png-output.ini'
        timeAnalyzedFile = ifarDir + '/corse_all_data_files/exclude_play/' + \
            combo + '_' + combo + '-CORSE_EXCLUDE_PLAY_' + \
            cp.options('ifar-mass-categories')[0] + '_CAT_3-' + gpsStartTime + \
             '-' + gpsDuration + '.txt'
        dagman,computeposteriorJobs[source+combo+type]=\
            create_computeposterior_job(\
            dagman, sourceChars,outputDir,galaxiesFile,timeAnalyzedFile)
        if runNumgalaxies:
          dagmanParentChild = add_parent_child(\
              [numgalaxiesJobs[source+combo+type]],\
              computeposteriorJobs[source+combo+type],dagmanParentChild)
  create_computeposterior_subfile(cp,executable,logPath,priority,userTag)

plotulvsmassJobs = {}
executable = cp.get('executables','plotulvsmass')

if runSpin or runNonSpin:
  if not os.path.isdir('plotulvsmass_files'):
    os.mkdir('plotulvsmass_files')
else:
  runPlotulvsmass = False

if runPlotulvsmass:
  for type in spinTypes:
    for combo in ifoCombos:
      if runMtotal:
        computepostGlob = 'plotnumgalaxies_files/' + type + '/mtotal_*/'\
            + combo + '/' + userTag + '-upper-limit'
        massRegion = 'totalmass'
        dagman,plotulvsmassJobs[combo + type] = \
            create_plotulvsmass_job(dagman,computepostGlob,massRegion,\
            'mtotal_' + type + '_' + combo,'plotulvsmass_files')    
        for item in mtotalSources:
          if runComputeposterior:
            dagmanParentChild = add_parent_child(\
                [computeposteriorJobs[item+combo+type]],\
                plotulvsmassJobs[combo+type],dagmanParentChild)
      if runComponent:
        computepostGlob = 'plotnumgalaxies_files/' + type + '/mcomp_*/'\
            + combo + '/' + userTag + '-upper-limit'
        massRegion = 'componentmass'
        dagman,plotulvsmassJobs[combo + type] = \
            create_plotulvsmass_job(dagman,computepostGlob,massRegion,\
            'mcomp_' + type + '_' + combo,'plotulvsmass_files')
        for item in componentSources:
          if runComputeposterior:
            dagmanParentChild = add_parent_child(\
                [computeposteriorJobs[item+combo+type]],\
                plotulvsmassJobs[combo+type],dagmanParentChild)
  create_plotulvsmass_subfile(cp,executable,logPath,priority)

combineposteriorJobs = {}

executable = cp.get('executables','pylal_combine_posteriors')
executable2 = cp.get('executables','plotulvsmass')

if (runSpin or runNonSpin):
  if not os.path.isdir('combine_posteriors'):
    os.mkdir('combine_posteriors')
  for type in spinTypes:
    for source in sourceTypes:
      if not os.path.isdir('combine_posteriors/' + type + '/' + source):
        os.makedirs('combine_posteriors/' + type + '/' + source)
else:
  runCombinePosterior = False

if runCombinePosterior:
  for type in spinTypes:
    childJobs = []
    for item in sourceTypes:
      posteriorFiles = []
      parentJobs = []
      if len(item.split('_')) == 3:
        massLow = (item.split('_'))[1]
        massHigh = (item.split('_'))[2]
      else:
        sourceChars = define_mass_characteristics(cp,item)
        massLow = sourceChars['min_mtotal']
        massHigh = sourceChars['max_mtotal']
      if not combineOnlyPast:
        initialDir = 'plotnumgalaxies_files/' + type + '/' + item
        figureName = userTag + '_' + item + '_' + type
        figureName += '_combos'
        for combo in ifoCombos:
          posteriorFiles.append((numgalaxiesOutputDirs[item+combo+type]+'/'+\
              userTag + '-posterior-pdf.txt', userTag + '_' + combo))
          if runComputeposterior:
            parentJobs.append(computeposteriorJobs[item+combo+type])
        posteriorFiles.sort()
        dagman,combineposteriorJobs[item + type + 'nopast'] = \
            create_combineposterior_job(dagman,posteriorFiles,\
            initialDir,figureName,'../../../',massLow,massHigh,item.upper())
        if runComputeposterior:
          dagmanParentChild = add_parent_child(parentJobs,\
              combineposteriorJobs[item+type+'nopast'],dagmanParentChild)
      initialDir = 'combine_posteriors/' + type + '/' + item
      figureName = userTag + '_' + item + '_' + type
      for dummyTag,info in cp.items('past-posteriors-' + item + '-' + type):
        option,value = info.split(',')
        posteriorFiles.append((value,option))
      posteriorFiles.sort()
      dagman,combineposteriorJobs[item + type] = \
          create_combineposterior_job(dagman,posteriorFiles,\
          initialDir,figureName,'../../../',massLow,massHigh,item.upper())
      childJobs.append(combineposteriorJobs[item+type])
      if runComputeposterior:
        dagmanParentChild = add_parent_child(parentJobs,\
            combineposteriorJobs[item+type],dagmanParentChild)

    if runMtotal:
      computepostGlob = 'combine_posteriors/' + type + '/mtotal_*/'\
          + userTag + '_*_' + type + '-combined-upper-limit'
      massRegion = 'totalmass'
      dagman,plotulvsmassJobs[type] = create_plotulvsmass_job(dagman,\
                 computepostGlob,massRegion,'mtotal_' + type + '-combined',\
                 'combine_posteriors')
      dagmanParentChild = add_parent_child(childJobs,\
          plotulvsmassJobs[type],dagmanParentChild)
    if runComponent:
      computepostGlob = 'combine_posteriors/' + type + '/mcomp_*/'\
          +  userTag + '_*_' + type + '-combined-upper-limit'
      massRegion = 'componentmass'
      dagman,plotulvsmassJobs[type] = create_plotulvsmass_job(dagman,\
                 computepostGlob,massRegion,'mcomp_' + type + '-combined',\
                 'combine_posteriors')  
      dagmanParentChild = add_parent_child(childJobs,\
          plotulvsmassJobs[type],dagmanParentChild)

  create_combineposterior_subfile(cp,executable,logPath,priority)
  create_plotulvsmass_subfile(cp,executable2,logPath,priority)


dagmanFile = open('upper_limit.dag','w')
dagFile = dagman + '\n' + dagmanParentChild + '\n' + dagmanMaxJobs
dagmanFile.write(dagFile)
dagmanFile.close()

