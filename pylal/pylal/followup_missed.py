# Copyright (C) 2006  Alexander Dietz
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

__author__ = "Darren Woods and Stephen Fairhurst <sfairhurs@gravity.phys.uwm.edu>"
__prog__ = "followup_missed.py"
__title__ = "Followup missed injections"

import os, sys, exceptions, copy
from math import sqrt, pi

from pylab import rcParams, fill, figtext, figure, plot, axes, axis, xlabel, ylabel, title, close, grid, legend

from pylal import SnglInspiralUtils
from pylal import InspiralUtils
from pylal import SimInspiralUtils
from pylal import CoincInspiralUtils
from pylal import SearchSummaryUtils
from pylal import git_version
from glue import lal
from glue import markup
from glue import pipeline
from glue import segments
from glue import segmentsUtils
from glue.markup import oneliner as extra
from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw import utils
import numpy


##########################################################
class FollowupMissed:
  """
  This defines a class for followup missed injections or follow up any
  injections in general. This class essentially creates time-series of triggers
  around the time of a missed injection for the various steps of the pipeline. 

  Usage:

  # first need to initialize the class
  followup = followup_missed.FollowupMissed( cache, opts)

  # later, one can do a followup of a missed injection 'inj', which returns the filename
  # of the html file containing the tables and pictures created
  followuphtml = followup.followup( inj, ifo )
  """


  # -----------------------------------------------------
  def __init__(self, cache, opts):
    """
    Initialize this class and sets up all the cache files.
    @param cache: The cache of all files
    @param opts: The 'opts' structure from the main code
    """
    rcParams.update({'text.usetex': False})

    # Check all the required options
    option_list = ['verbose','followup_exttrig','output_path',\
                   'followup_time_window','followup_tag','prefix',\
                   'suffix','figure_resolution']
    for option in option_list:
      if not hasattr(opts, option):
        raise "Error: The following parameter is required in the "\
              "opts structure: ", option
    

    # defining some colors and the stages
    self.colors = {'H1':'r','H2':'b','L1':'g','V1':'m','G1':'c'}
    self.stageLabels = ['INSPIRAL_FIRST', 'THINCA_FIRST',\
                        'INSPIRAL_SECOND', 'THINCA_SECOND']
    self.orderLabels = copy.deepcopy(self.stageLabels)
    self.orderLabels.extend( [ 'THINCA_SECOND_CAT_1','THINCA_SECOND_CAT_2', \
                               'THINCA_SECOND_CAT_3','THINCA_SECOND_CAT_4'] )

    # set arguments from the options
    self.opts = opts
    self.verbose = opts.verbose
    self.exttrig = opts.followup_exttrig
    self.output_path = opts.output_path
    self.time_window = opts.followup_time_window
    
    # default value for the injection time-window. This value
    # will be taken later from a processParams table (if possible)
    self.injection_window = 0.050

    # initialize a list of images created
    self.fnameList = []

    # counter for the followups
    self.number = 0

    # for the estimated distances
    self.flow = opts.followup_flow
    self.spectrum = createSpectrum( self.flow, \
                                    sampleRate = 4096, \
                                    nPoints = 1048576)

    if self.verbose:
      print "\nStarting initializing the Followup class..."

    # setting the caches
    if opts.followup_tag:
      self.cache = cache.sieve(description = opts.followup_tag)
    else:
      self.cache = cache
      
    # retrieving the caches for the different stages
    self.triggerCache = {}
    for stage in self.stageLabels:
      pattern = stage
      self.triggerCache[stage] = self.cache.sieve(description=pattern)
      if self.opts.verbose:
        print "%d files found for stage %s" % (len(self.triggerCache[stage]), stage)

    # sieve the injections
    self.injectionCache = self.cache.sieve(description = "INJECTION").\
                          sieve(ifos='HL')
    
    # generate a dictionary based on the event-ID
    self.injections=dict()
    for file, entry in zip(self.injectionCache.pfnlist(), self.injectionCache):
      injection_id = self.get_injection_id(cache_entry = entry)
      self.injections[injection_id] = SimInspiralUtils.\
                                        ReadSimInspiralFromFiles( [file], verbose=False )
    if self.verbose:
      print "parsing of cache files for the Followup class done..."

    

    # get the process params table from one of the COIRE files
    coire_file = self.cache.sieve(description = "FOUND").checkfilesexist()[0].pfnlist()[0]
    try:
      doc = SearchSummaryUtils.ReadTablesFromFiles([coire_file],[ lsctables.ProcessParamsTable])
      process_params = table.get_table(doc, lsctables.ProcessParamsTable.\
                                       tableName)
    except IOError: 	    
      sys.stderr.write("ERROR (IOError) while reading process_params table from file %s. "\
                       "Does this file exist and does it contain a search_summary table?\n" %(coire_file))
      raise 	 
    except:
      raise "Error while reading process_params table from file: ", coire_file

    # and retrieve the time window from this file
    found_flag = False
    for tab in process_params:
      if tab.param=='--injection-window':
        found_flag = True
        self.injection_window = float(tab.value)/1000.0
    if not found_flag: 	 
      sys.stderr.write("WARNING: No entry '--injection-window' found in file %s"\
                       "Value used is %.1f ms. If incorrect, please change file at %s\n" %\
                       (coire_file, 1000.0*self.injection_window, __file__))
    if self.verbose:
      print "Injection-window set to %.0f ms" % (1000*self.injection_window)

    # read the veto files
    self.readVetoFiles()
    
    if self.verbose:
      print "Initializing the Followup class done..."


  # -----------------------------------------------------
  def get_injection_id(self, filename=None, url=None, cache_entry=None):
    """
    Extracting the injection ID from the filename, using
    the mechanism as used in lalapps_path2cache. You can specify the filename
    itself, the url or the cache entry. You must not specify more than one input!
    The injection-ID is calculated in the following way (exttrig only):

    The code expects the INSPIRAL and THINCA files in the following scheme (example):
      PREFIX-TAGPART_injections32_77-GPS-DURATION.xml
    The number of the injection run is extracted (INJRUN) as well as the following
    number (INJNUMBER). The injection ID is then calculated as:
       INJID = 100000*INJRUN + INJNUMBER
    so for this example the injectionj ID is 3200077. 
    
    @param file: filename from which the injection ID is extracted
    @param url:  url from which the injection ID is extracted
    @param cache_entry: cache entry from which the injection ID is extracted
    """
    
    # Check that only one input is given
    if cache_entry:
      if filename and url:
        raise "Error in function get_injection_id: Only one input should be "\
              "specified. Now 'cache_entry' and another variable is specified. Check the code."
    
    if cache_entry is None:
      if filename and url:
        raise "Error in function get_injection_id: Only one input should be "\
              "specified. Now 'filename' and 'url' is specified. Check the code."
      
      if filename:
        path, filename = os.path.split(filename.strip())
        url = "file://localhost%s" % os.path.abspath(os.path.join(path, filename))

      try:
        cache_entry = lal.CacheEntry.from_T050017(url)
      except ValueError, e:
        raise "Error while extracting injection ID from file ", filename

    # split the expression into different sub-pieces
    pieces = cache_entry.description.split('_')
    if self.exttrig:

      # its easy for the exttrig case
      injection_id = pieces[-2]+'_'+pieces[-1]
    else:

      # but need to check for the appearance of the CAT suffix else
      index = 0
      for ind, piece in enumerate(pieces):
        if 'CAT' in piece:
          index = ind          
      injection_id = pieces[index-1] 
        
    return injection_id

  # -----------------------------------------------------
  def readVetoFiles( self ):
    """
    Reads the veto segments given by veto-files (if any)
    """
    self.vetodict = dict()

    # loop over the IFO names
    for ifoName in self.colors:

      self.vetodict[ifoName]=None
      # create the attribute name and check if it exists
      attributeName = 'followup_vetofile_'+ifoName.lower()
      if hasattr( self.opts, attributeName):

         # get the filename
         filename = getattr( self.opts, attributeName )

         if filename:
           #read the segment lists
           self.vetodict[ifoName] = segmentsUtils.fromsegwizard(open(filename))

  # -----------------------------------------------------
  def reset( self ):
    """
    Resets the counting number for the time-series plots generated.
    """
    self.number=0

  # -----------------------------------------------------
  def setTag( self, tag ):
    """
    Just sets the tag, called from plotinspmissed (why needed?)
    @param tag: tag for this followup
    """

    self.tag = tag

  # -----------------------------------------------------
  def print_inj( self, inj, injID ):
    """
    Print some useful informations to the screen.
    @param inj:   the current missed injection
    @param injID: the injection ID (used for exttrig only)
    """

    if self.exttrig:
      print "\nInjection details for injection %d with injID %s: " %\
            (self.number, injID)
    else:
      print "\nInjection details for injection %d:" % (self.number)    
      
    print "m1: %.2f  m2:%.2f  | end_time: %d.%d | "\
          "distance: %.2f  eff_dist_h: %.2f eff_dist_l: %.2f" % \
          ( inj.mass1, inj.mass2, inj.geocent_end_time, inj.geocent_end_time_ns,\
            inj.distance, inj.eff_dist_h, inj.eff_dist_l )

  # ----------------------------------------------------
  def savePlot( self, stage ):
    """
    Saves the plots and store them in a seperate fnameList.
    @param stage: the stage this plot belongs to (e.g. INSPIRAL, THINCA,...)
    """
    fname = 'Images/'+self.opts.prefix + "_"+self.tag+"_map-"+\
            stage+"-"+str(self.number) +self.opts.suffix+'.png'
    fname_thumb = InspiralUtils.\
                  savefig_pylal( filename = self.output_path+fname,\
                                 doThumb = True, 
                                 dpi_thumb = self.opts.figure_resolution)
    
    self.fnameList.append( fname ) 
    return fname
 

  # -----------------------------------------------------  
  def findInjection( self, missedInj ):
    """
    Find the injection-ID corresponding to this particular
    missed injection.
    @param missedInj: the missed injection
    """
    
    # injID: the ID (a number)
    # groupInj: a list of SimInspiral tables
    for injID, groupInj in self.injections.iteritems():
      for inj in groupInj:
        if missedInj.geocent_end_time==inj.geocent_end_time and \
               missedInj.geocent_end_time_ns==inj.geocent_end_time_ns:
          return injID
          
    self.print_inj( missedInj, None)
    raise "No injection ID found for the above particular missed Injection " 

    return None

  # -----------------------------------------------------  
  def getTimeTrigger( self, trig ):
    """
    This is a helper function to return a GPS time as one float number
    @param trig: a sngl_inspiral table entry
    """

    return float(trig.end_time) + float(trig.end_time_ns) * 1.0e-9
  

  # -----------------------------------------------------
  def getTimeSim(self, sim, ifo=None ):
    """
    This is a helper function to return a GPS time as one float number
    for a certain IFO. If no IFO is specified the injected geocentric
    time is returned.
    @param sim: a sim_inspiral table entry
    @param ifo: the IFO for which we want the sim time
    """
    
    time=0
    nano=0

    if not ifo:
      time = sim.geocent_end_time
      nano = sim.geocent_end_time_ns
    if ifo:
      time = sim.get(ifo[0].lower()+'_end_time' )
      nano = sim.get(ifo[0].lower()+'_end_time_ns' )    

    return  float(time) + float(nano) * 1.0e-9

  # -----------------------------------------------------
  def highlightVeto( self, timeInjection, segLarge, ifoName, ylims  ):
    """
    Finds the intersection of the drawn triggers with the veto segments
    for this IFO
    """

    if not self.vetodict[ifoName]:
      return
    
    # find intersecting veto segments
    segVeto =  self.vetodict[ifoName] & segments.segmentlist( [segLarge] )

    # draw all intersecting segments
    for seg1 in segVeto:

      # after a tim-shift
      seg=seg1.shift( -timeInjection)
      vetox = [seg[0], seg[1], seg[1], seg[0], seg[0]]
      vetoy = [ylims[0], ylims[0], ylims[1], ylims[1], ylims[0]]
      fill ( vetox, vetoy, 'y', alpha=0.2)  


  # -----------------------------------------------------
  def isThereVeto( self, timeTrigger, ifoName ):
    """
    This function checks if at the time 'timeTrigger' the
    IFO 'ifoName' is vetoed.
    @param timeTrigger: The time to be investigated
    @param ifoName: The name of the IFO to be investigated
    """
    
    if  self.vetodict[ifoName] is None:
      flag = False
    else:
      flag = timeTrigger in self.vetodict[ifoName]

    return flag

  # -----------------------------------------------------
  def getExpectedSNR(self, triggerFiles, inj, number ):
    """
    Investigate template bank and returns exepcted horizon distance
    @param triggerFiles: List of files containing the inspiral triggers
    @param inj: The current injection
    @param ifo: The current IFO
    @param number: The consecutive number for this inspiral followup
    """
    
    # read the inspiral file(s)
    if self.verbose: print "Processing TMPLTBANK triggers from files ",\
       triggerFiles      
    inspiralSumm, massInspiralSumm = InspiralUtils.\
                                     readHorizonDistanceFromSummValueTable(\
                                     triggerFiles, self.verbose)

    # get different informations on the missed injection
    injMass = [inj.mass1, inj.mass2]
    timeInjection = self.getTimeSim( inj )    
    totalMass  = inj.mass1 + inj.mass2
    eta = inj.mass1 * inj.mass2 / totalMass / totalMass

    # Given the masses (and therefore eta), the expected horizon distance(SNR) has to be scaled by this factor
    factor = sqrt(4 * eta)

    output = {}
    # for each ifo that has been found
    for ifo in massInspiralSumm.keys():     
      # loop over the summ value table as long as horizon is not set
      output[ifo] = [] 
      horizon = 0 
      for massNum in range(len(massInspiralSumm[ifo])):
        # looking at the masses and horizon distance (and threshold)
        if horizon > 0:
          break
        for this, horizon  in zip(massInspiralSumm[ifo][massNum].\
                                  getColumnByName('comment'),\
            massInspiralSumm[ifo][massNum].getColumnByName('value').asarray()):
          masses = this.split('_')
          threshold = float(masses[2])
          readTotalMass = float(masses[0])+float(masses[1])
          # check if the current total Mass is greater tan the requested total Mass
          # if so, then keep this horizon and break 
          if readTotalMass > totalMass:
            startTimeSec = \
              massInspiralSumm[ifo][massNum].getColumnByName('start_time').asarray()
            #sanity check 
            if len(startTimeSec)>1: 
              print >> sys.stderr, 'Warning in fct. expectedHorizonDistance: '\
                    'More than 1 file found at particular GPS time. Using the first file.' 
            if startTimeSec[0] > inj.geocent_end_time:
              text= """the start time of the template bank must be less than 
		the end_time of the injection. We found startTime of the bank to be %s and
		 the geocent end time of the injection to be %s""",startTimeSec[0], inj.geocent_end_time
              raise ValueError,  text
            output[ifo] = horizon * factor * threshold / float(getattr(inj, 'eff_dist_'+ifo[0].lower() ))
            break
          #'otherwise, reset horizon distance
          else: horizon = 0           

    return output

  # -----------------------------------------------------
  def putText( self, text ):
    """
    Puts some text into an otherwise empty plot.
    @param text: text to put in the empty plot
    """
    newText = ''
    for i in range( int(len(text)/60.0)+1):
      newText+=text[60*i:60*i+60]+'\n'
    figtext(0.15,0.15, newText)
 
  # -----------------------------------------------------
  def investigateTimeseries(self, triggerFiles, inj,  ifoName, stage, number ):
    """
    Investigate inspiral triggers and create a time-series
    of the SNRs around the injected time
    @param triggerFiles: List of files containing the inspiral triggers
    @param inj:          the current missed injection
    @param ifoName:      the IFO for which the plot is made 
    @param stage:        the name of the stage (FIRST, SECOND)
    @param number:        the consecutive number for this inspiral followup
    """
    
    # read the inspiral file(s)
    if self.verbose:
      print "Processing INSPIRAL triggers from files ", triggerFiles
      
    snglTriggers = SnglInspiralUtils.ReadSnglInspiralFromFiles( \
      triggerFiles , verbose=False)

    # create a figure and initialize some lists
    fig=figure()
    foundSet = set()
    loudest_details = {}
    noTriggersFound = True
    
    if snglTriggers is None:
      # put message on the plot instead
      self.putText( 'No sngl_inspiral triggers in %s' % str(triggerFiles))

    else:
      # selection segment
      timeInjection = self.getTimeSim( inj )
      segSmall =  segments.segment( timeInjection-self.injection_window, \
                                    timeInjection+self.injection_window )
      segLarge =  segments.segment( timeInjection-self.time_window, \
                                    timeInjection+self.time_window )

      # create coincidences for THINCA stage
      coincTriggers = None
      if 'THINCA' in stage:
        coincTriggers = CoincInspiralUtils.coincInspiralTable( snglTriggers, \
                      CoincInspiralUtils.coincStatistic("snr") )
        selectedCoincs = coincTriggers.vetoed( segSmall )
      
      # loop over the IFOs (although this is a plot for IFO 'ifoName')
      for ifo in self.colors.keys():

        # get the singles for this ifo
        snglInspiral = snglTriggers.ifocut(ifo)

        # select a range of triggers
        selectedLarge = snglInspiral.vetoed( segLarge )
        timeLarge = [ self.getTimeTrigger( sel )-timeInjection \
                      for sel in selectedLarge ]

        selectedSmall = snglInspiral.vetoed( segSmall )
        timeSmall = [ self.getTimeTrigger( sel )-timeInjection \
                      for sel in selectedSmall ]

        # use the set of selected coincident triggers in the THINCA stages
        if coincTriggers:
          selectedSmall = selectedCoincs.cluster(2* self.injection_window).getsngls(ifo)
          timeSmall = [ self.getTimeTrigger( sel )-timeInjection \
                        for sel in selectedSmall ]
          
        # skip if no triggers in the large time window
        if len(timeLarge)==0:
          continue
        noTriggersFound = False

        # add IFO to this set; the injection is found for this IFO and stage
        if len(timeSmall)>0:
          foundSet.add(ifo)                  

          # record details of the loudest trigger
          loudest_details[ifo] = {}
          loudest = selectedSmall[selectedSmall.get_column('snr').argmax()]
          loudest_details[ifo]["snr"] = loudest.snr
          loudest_details[ifo]["mchirp"] = loudest.mchirp
          loudest_details[ifo]["eff_dist"] = loudest.eff_distance
          loudest_details[ifo]["chisq"] = loudest.chisq
          loudest_details[ifo]["timeTrigger"] = self.getTimeTrigger( loudest )

          timeTrigger = self.getTimeTrigger( loudest )
          vetoSegs = self.vetodict[ifoName]
          
        # plot the triggers
        plot( timeLarge, selectedLarge.get_column('snr'),\
              self.colors[ifo]+'o', label="_nolegend_")
        plot( timeSmall, selectedSmall.get_column('snr'), \
              self.colors[ifo]+'s', label=ifo)

      # draw the injection times and other stuff
      if noTriggersFound:
        self.putText( 'No triggers/coincidences found within time window')
        
      ylims=axes().get_ylim()
      plot( [0,0], ylims, 'g--', label="_nolegend_")
      plot( [-self.injection_window, -self.injection_window], ylims, 'c:',\
            label="_nolegend_")
      plot( [+self.injection_window, +self.injection_window], ylims, 'c:',\
            label="_nolegend_")

      self.highlightVeto( timeInjection, segLarge, ifoName, ylims  )

      # save the plot
      grid(True)
      legend()

    ylims=axes().get_ylim()
    axis([-self.time_window, +self.time_window, ylims[0], ylims[1]])
    xlabel('time [s]')
    ylabel('SNR')
    title(stage+'_'+str(self.number))    
    fname = self.savePlot( stage )
    close(fig)

    result = {'filename':fname, 'foundset':foundSet, 'loudest_details':loudest_details}
    return result


  # -----------------------------------------------------    
  def select_category(self, trigger_files, category):
    """
    Return a trigger list that contains only files for the choosen category.
    @param triggerList : a list of file names
    @param category: a category number
    @return: a sub list of filename corresponding to the category requested
    """

    # there are two different labels used to denote
    # the categories. THIS NEEDS TO BE UNIFIED
    cat1 = 'CAT_'+str(category)
    cat2 = 'CATEGORY_'+str(category)
    
    if category==1:
      # Category 1 files might not be labelled with a 'CAT_1' suffix.
      # So, for now, all files NOT containing the expression
      # 'CAT' in the filename are supposed to be CAT_1 files
      new_list = [file for file in trigger_files \
                  if 'CAT' not in file or cat1 in file or cat2 in file]
                     
    else:
      cat = 'CAT_'+str(category)
      new_list = [file for file in trigger_files if cat1 in file\
                  or cat2 in file]      
          
    return new_list


  # -----------------------------------------------------
  def followup(self, inj, selectIFO, description = None):
    """
    Do the followup procedure for the missed injection 'inj'
    and create the several time-series for INSPIRAL and THINCA.
    The return value is the name of the created html file.
    @param inj: sim_inspiral table of the injection that needs to be
                followed up
    @param selectIFO: the IFO that is investigated
    @param description: Can be used to sieve further this pattern
                        from the description field.
    """
    
    def fill_table(page, contents ):
      """
      Making life easier...
      """
      page.add('<tr>')
      for content in contents:
        page.add('<td>')
        page.add( str(content) )
        page.add('</td>')
      page.add('</tr>')

   
    # get the ID corresponding to this injection
    injection_id = self.findInjection( inj )

    # increase internal number:
    self.number+=1

    ## create the web-page and add a table
    page = markup.page()
    page.h1("Followup missed injection #"+str(self.number)+" in "+selectIFO )
    page.hr()
    page.add('<table border="3" ><tr><td>')
    page.add('<table border="2" >')          
    fill_table( page, ['<b>parameter','<b>value'] )
    fill_table( page, ['Number', self.number] )
    fill_table( page, ['inj ID', injection_id] )
    fill_table( page, ['mass1', '%.2f'% inj.mass1] )
    fill_table( page, ['mass2', '%.2f'%inj.mass2] )
    fill_table( page, ['mtotal', '%.2f' % (inj.mass1+inj.mass2)] )
    fill_table( page, ['mchirp', '%.2f' % (inj.mchirp)] )
    fill_table( page, ['end_time', inj.geocent_end_time] )
    fill_table( page, ['end_time_ns', inj.geocent_end_time_ns] )    
    fill_table( page, ['distance', '%.1f' % inj.distance] )
    fill_table( page, ['eff_dist_h','%.1f' %  inj.eff_dist_h] )
    fill_table( page, ['eff_dist_l','%.1f' %  inj.eff_dist_l] )
    fill_table( page, ['eff_dist_v','%.1f' %  inj.eff_dist_v] )
    fill_table( page, ['eff_dist_g','%.1f' %  inj.eff_dist_g] )  
    fill_table( page, ['playground','%s' %  pipeline.s2play(inj.geocent_end_time)] )    
    page.add('</table></td>')
    
    # print infos to screen if required
    if self.opts.verbose:
      self.print_inj( inj,  injection_id)

    # sieve the cache for the required INSPIRAL and THINCA files
    invest_dict = {}
    for stage, cache in self.triggerCache.iteritems():

      trig_cache = lal.Cache()
      for c in cache:

        # check the time and the injection ID
        if inj.geocent_end_time in c.segment:
          if self.get_injection_id(url = c.url) == injection_id:
            trig_cache.append( c )

      # create a filelist
      file_list = trig_cache.sieve(description = description).pfnlist()
        
      # check if the pfnlist is empty. `
      if len(file_list)==0:
        print >>sys.stderr, "Error: No files found for stage %s in the "\
              "cache for ID %s and time %d; probably mismatch of a "\
              "pattern in the options. " % \
              ( stage, injection_id, inj.geocent_end_time)        
        continue

      # if the stage if THINCA_SECOND...
      if 'THINCA_SECOND' in stage:

        # ... need to loop over the four categories
        for cat in [1,2,3,4]:
          
          select_list=self.select_category( file_list, cat)
          if len(select_list)==0:
            print "WARNING: No THINCA_SECOND files found for category ", cat
            continue
          
          modstage = stage+'_CAT_' + str(cat)
          invest_dict[modstage] = self.investigateTimeseries( select_list, inj, selectIFO, modstage, self.number )

        #sys.exit(0)
      else:
        invest_dict[stage]=self.investigateTimeseries( file_list, inj, selectIFO, stage, self.number)

      
      
    ## print out the result for this particular injection
    page.add('<td><table border="2" >')
    fill_table( page, ['<b>step','<b>F/M', '<b>Rec. SNR', '<b>Rec. mchirp', \
                      '<b>Rec. eff_dist', '<b>Rec. chisq', '<b>Veto ON/OFF'] )

    # loop over the stages and create the table with
    # the various data in it (when available)
    for stage in self.orderLabels:
      if stage in invest_dict:
        result = invest_dict[stage]

        # Fill in the details of the loudest found coinc.
        #found_ifo=''
        #if "INSPIRAL" in stage or "THINCA" in stage:
        found_ifo=''
        loudest_snr=''
        loudest_mchirp=''
        loudest_eff_dist=''
        loudest_chisq=''
        veto_onoff=''

        # add all the IFO's for this coincident
        for ifo in result['foundset']:
          found_ifo += ifo+' '
          
          # Parameters of the loudest trigger, taken from the
          # 'loudest-details' dictionary, created in 'investigateTimeseries'
          loudest_snr += ifo + ': ' + str(result['loudest_details'][ifo]['snr'])+'<br>'
          loudest_mchirp += ifo + ': ' + str(result['loudest_details'][ifo]['mchirp'])+'<br>'
          loudest_eff_dist += ifo + ': ' + str(result['loudest_details'][ifo]['eff_dist'])+'<br>'
          loudest_chisq += ifo + ': ' + str(result['loudest_details'][ifo]['chisq'])+'<br>'
          
          # Check whether some of the ifo times is vetoed
          timeTrigger = float(result['loudest_details'][ifo]['timeTrigger'])
          if (self.vetodict[ifo]):
            veto = self.isThereVeto (timeTrigger, ifo)
            veto_txt = 'OFF'
            if veto:
              veto_txt = 'ON'              
            veto_onoff+=ifo+': '+veto_txt+'<br>'
          else: 
            veto_onoff+=ifo+': No info<br>'

        # Fill the table whether something is found or not
        if len(result['foundset'])>0:
          fill_table( page, [ stage,  'FOUND in '+found_ifo, 'loudest<br>'+loudest_snr, \
                              'loudest<br>'+loudest_mchirp, 'loudest<br>'+loudest_eff_dist,\
                              'loudest<br>'+loudest_chisq, veto_onoff])
        else:
          fill_table( page, [ stage,  '<font color="red">MISSED'])
          
    page.add('</table>')
    page.add('</td></tr></table><br><br>')


    ## add the pictures to the webpage
    for stage in self.orderLabels:
      if stage in invest_dict:
        result = invest_dict[stage]
      
        ##if stage!="TMPLTBANK":
        if True:
          fname = result['filename']
          page.a(extra.img(src=[fname], width=400, \
                           alt=fname, border="2"), title=fname, href=[ fname ])
          
    # add version information
    page.add('<hr>Page created with %s Version %s' % \
        (__prog__, git_version.verbose_msg))
    
    # and write the html file
    htmlfilename = self.opts.prefix + "_"+selectIFO+"_followup_"+str(self.number) +\
                         self.opts.suffix+'.html'
    file = open(self.opts.output_path+htmlfilename,'w')      
    file.write(page(False))
    file.close()

    # store html file in fnameList
    self.fnameList.append(htmlfilename)

    # supply the output
    return htmlfilename

  # -----------------------------------------------------
  def estimatedDistance( self, mass1, mass2, distTarget):
    """
    Calculates the approx. effective distance (in Mpc) for a
    binary of two objects with masses 'mass1' and 'mass2',
    when the efective distance for a 10/10 binary is 'distTarget'.
    This distance is rescaled given the two masses and the spectrum.
    Example

      estimatedDistance(10.0, 10.0, 100.0)
      should return exactly 100.0 again.
      
    @param mass1: The mass of the first component (Solar masses)
    @param mass2: The mass of the second component (Solar masses)
    @param distTarget: Effective distance for a 10/10 binary
    """

    snrActual = 8.0
    distActual = computeCandleDistance( 10.0, 10.0, \
                                        self.flow, self.spectrum, \
                                        snrActual)

    # rescale the SNR threshold
    snrTarget = distActual / distTarget*snrActual 
    distance = computeCandleDistance( mass1, mass2, \
                                      self.flow, self.spectrum,\
                                      snrTarget)
    return distance

  
# -----------------------------------------------------
def createSpectrum( flow, sampleRate = 4096, nPoints = 1048576):
  """
  Computes the approximate spectrum of LIGO-I.
  @param flow: The lower cutoff frequency
  @param sampleRate: The sample rate used (default: 4096)
  @param nPoints The number of points (default: 1048576)
  """

  def calcPSD( f ):
    f0=150.0
    a = pow(4.49*f/f0, -56.0)
    b = 0.16*pow(f/f0, -4.52)
    c = 0.52
    d = 0.32*pow(f/f0,2)
    return (a+b+c+d)*9.0e-46
  
  chanDeltaT = 1.0/float(sampleRate)
  deltaF =  1.0/( float(nPoints)*chanDeltaT )
  
  spectrum = []
  for k in range( nPoints/2+1):
    f = k*deltaF
    if f<flow:
      spectrum.append(0)
    else:
      spectrum.append( calcPSD(f) )
    
  return spectrum

# -----------------------------------------------------
def computeCandleDistance( candleM1, candleM2, fLow, spectrum = None, \
                           snr=8.0, sampleRate = 4096, nPoints = 1048576):
  """
  Computes the candle distance as computed in inspiralutils.c.
  This code has been tested on 21 Feb 2008 and the derivation between
  these calculated values with the LAL code is less than 2 Percent.
  @param candleM1: The mass of the first component (Solar masses)
  @param candleM2: The mass of the second component (Solar masses)
  @param flow: The lower cutoff frequency
  @param spectrum: The PSD that has been calculated earlier
  @param snr: The SNR at which the distance is calculated (default: 8)
  @param sampleRate: The sample rate used (default: 4096)
  @param nPoints The number of points (default: 1048576)
  """
  # TestCode: testEstimatedDistance.png
  

  LAL_MTSUN_SI = 4.92549095e-6
  chanDeltaT = 1.0/float(sampleRate)
  negativeSevenOverThree = -7.0/3.0
  totalMass = candleM1 + candleM2
  mu = candleM1 * candleM2 / totalMass
  distNorm = 9.5708317e-20 # = 2.0 * LAL_MRSUN_SI / (1.0e6 * LAL_PC_SI )
  a = sqrt( (5.0 * mu) / 96.0 ) * \
      pow( totalMass / ( pi*pi ), 1.0/3.0 ) *\
      pow( LAL_MTSUN_SI / chanDeltaT, -1.0/6.0 )
  sigmaSq = 4.0 * ( chanDeltaT / float(nPoints) ) * \
            distNorm * distNorm * a * a
  
  fmax = 1.0 / (6.0 * sqrt(6.0) * pi * totalMass * LAL_MTSUN_SI)
  deltaF =  1.0/( float(nPoints)*chanDeltaT )

  cut = int( fLow/deltaF )
  kmax = min( int( fmax/deltaF ), len(spectrum) )

  sigmaSqSum = 0
  for k in range( cut, kmax):
    sigmaSqSum += pow( float(k)/float(nPoints),  negativeSevenOverThree ) /\
                  spectrum[k]

  sigmaSq *= sigmaSqSum
  distance = sqrt( sigmaSq ) / snr

  return distance
