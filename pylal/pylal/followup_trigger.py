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
from __future__ import division

__prog__ = "followup_trigger.py"

import os
import sys
import copy
from math import sqrt, pi
import subprocess
import tempfile

import pylab
import numpy

from pylal import SnglInspiralUtils
from pylal import InspiralUtils
from pylal import SimInspiralUtils
from pylal import CoincInspiralUtils
from pylal import SearchSummaryUtils
from pylal import git_version
from pylal import grbsummary
from pylal import viz
from pylal import tools
from glue import lal
from glue import markup
from glue import pipeline
from glue import segments
from glue import segmentsUtils
from glue import iterutils
from glue.markup import oneliner as extra
from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw import utils
from glue.ligolw import ligolw
from glue.ligolw.utils import ligolw_add


# DB content handler for reading xml input files
class ContentHandler(ligolw.LIGOLWContentHandler):
        pass
lsctables.use_in(ContentHandler)


##########################################################
class FollowupTrigger:
  """
  This defines a class for following up a trigger and to create
  the timeseries of found triggers in each stage of the pipeline.

  Usage:

  # first need to initialize the class
  followup = followup_trigger.FollowupTrigger(cache, opts)

  # later, one can do a followup of different objects, like an injection,
  #   a coinc. trigger or a single trigger:
  htmlname1 = followup.from_coinc(coinc)
  htmlname2 = followup.from_sngl(sngl_inspiral)
  htmlname3 = followup.from_missed(missed_inj)
  htmlname4 = followup.from_found(found_inj)
  htmlname5 = followup.from_new_coinc(new_coinc,[sngls])
  htmlname6 = followup.from_new_slide_coinc(new_coinc,[sngls],slideDict,segs)

  # In each case the path to the created html file is returned.
  # In the first call a CoincInspirals table is expected, a SnglInspiral
  # table in the second case and a SimInspiral table in the two last.
  """


  # -----------------------------------------------------
  def __init__(self, cache, opts, use_injections = True,do_slides=False):
    """
    Initialize this class and sets up all the cache files.
    @param cache: The cache of all files
    @param opts: The 'opts' structure from the main code
    @param use_injections: Specifying if injections are being used
                           (if no injections are used and this is set to
                           False, it speeds up the whole initalization...)
    """
    pylab.rcParams.update({'text.usetex': False})

    # Check all the required options
    if not hasattr(opts, 'followup_exttrig'):
      # default: this is a regular search, not a exttrig search
      opts.followup_exttrig = False
    if not hasattr(opts, 'followup_time_window'):
      # default: set the time window for the timeseries to 10 seconds
      opts.followup_time_window = 10.0
    if not hasattr(opts, 'followup_tag'):
      # default: don't specify a followup-tag
      opts.followup_tag = None
    if not hasattr(opts, 'followup_sned'):
      # default: do not incorporate the updating of effective distances
      #          with lalapps_sned
      opts.followup_sned = None

    option_list = ['verbose','followup_exttrig','output_path',\
                   'followup_time_window','prefix',\
                   'suffix','figure_resolution','user_tag',\
                   'followup_tag','followup_sned']
    for option in option_list:
      if not hasattr(opts, option):
        raise "Error: The following parameter is required in the "\
              "opts structure for followup_trigger: ", option
    

    # setting the color definition and the stages of the pipeline
    self.colors = {'H1':'r','H2':'b','L1':'g','V1':'m','G1':'c'}
    self.stageLabels = [('MATCHEDFILTER',['INSPIRAL'])]
    if do_slides:
      self.stageLabels.append(('COINCIDENCE-SLID',['THINCA_1','THINCA_SLIDES']))
    else:
      self.stageLabels.append(('COINCIDENCE',['THINCA_0','THINCA_ZEROLAG']))

    self.orderLabels = ['MATCHEDFILTER']
    if do_slides:
      self.orderLabels.extend( [ 'COINCIDENCE-SLID_CAT_1',\
          'COINCIDENCE-SLID_CAT_2', 'COINCIDENCE-SLID_CAT_3',\
          'COINCIDENCE-SLID_CAT_4', 'COINCIDENCE-SLID_CAT_5'] )
    else:
      self.orderLabels.extend( [ 'COINCIDENCE_CAT_1','COINCIDENCE_CAT_2', \
                                 'COINCIDENCE_CAT_3','COINCIDENCE_CAT_4', \
                                 'COINCIDENCE_CAT_5'] )

    # set arguments from the options
    self.opts = opts
    self.tag = opts.user_tag
    self.verbose = opts.verbose
    self.exttrig = opts.followup_exttrig
    self.output_path = opts.output_path
    self.time_window = opts.followup_time_window
    self.sned = opts.followup_sned

    # Set argument "old_document" to true if the option is specified.
    # Do not crash if "opts.old_document" is not defined in the script
    # calling this method.
    if not hasattr(opts, "old_document"):
      self.old_document = False
    else:
      self.old_document = opts.old_document

    # Choose the cache-file
    if opts.followup_tag is None:
      if opts.verbose:
        print "WARNING: All injection files are considered. Is that desired?"
        print "         Or do you want to select an injection run with --followup-tag INJRUN ?"
      self.cache = cache
    else:
      if opts.verbose:
        print "WARNING: Only a selected injection run is taken into account: ", opts.followup_tag
        print "         Is that the desired behaviour? You might reconsider removing --followup-tag INJRUN "
      self.cache = cache.sieve(description = opts.followup_tag)
    
    # Value for the injection window. This value might be taken from
    # a processParams table (see function get_injection_window)
    self.get_injection_window()

    # initialize a list of images created
    self.fname_list = []

    # counter for the followups
    self.number = 0

    # setting for each time to followup:
    self.followup_time = None
    self.injection_id = None  # only needed if time is an injection
    self.flag_followup = False

    # a dictionary used to create a list of injections to followup
    # taking into account the 'combined' effective distance
    # and vetoes as well
    self.followup_dict = {'inj':lsctables.New(lsctables.SimInspiralTable), 'dist':[], 'type':[]}
    self.vetoed_injections = lsctables.New(lsctables.SimInspiralTable)

    if self.verbose:
      print "\nStarting initializing the Followup class..."
      
    # splitting up the cache for the different stages
    self.trigger_cache = {}
    for stageName, stagePatterns in self.stageLabels:
      sievedFiles = []
      for pattern in stagePatterns:
          sievedFiles.extend(self.cache.sieve(description=pattern))
      self.trigger_cache[stageName] = sievedFiles
      if self.opts.verbose:
        print "%d files found for stage %s" %\
                                (len(self.trigger_cache[stageName]), stageName)


    # generating a dictionary for injection followups
    self.injections = dict()
    if use_injections:
      if self.verbose:
        print "Creating the injection cache..."

      self.injection_cache = self.cache.sieve(description = "INJECTION").\
                             sieve(ifos='HL')
      
      for file, entry in zip(self.injection_cache.pfnlist(), \
                             self.injection_cache):
        injection_id = self.get_injection_id(cache_entry = entry)
        self.injections[injection_id] = SimInspiralUtils.\
                                          ReadSimInspiralFromFiles( \
                                           [file], verbose=False )

    # read the veto files
    self.read_veto_files()
    
    if self.verbose:
      print "Initializing the Followup class done..."

  # -----------------------------------------------------
  def set_sned(self, executable):
    """
    Sets the sned flag to 'flag'
    @param flag: the executable used to calculate the sned
    """
    self.sned = executable

  # -----------------------------------------------------
  def execute_sned(self, injections):
    """
    Makes an external call to lalapps_sned to recalculate
    the effective distances
    @param injections: the list of injections to be converted
    @return: the recalculated injection
    """

    def fix_numrel_columns(sims):
       """	
       Due to a mismatch between our LAL tag and the LAL HEAD, against which
       we compile pylal, there are missing columns that cannot be missing.
       Just put in nonsense values so that glue.ligolw code won't barf.
       """
       for sim in sims:
         sim.numrel_data="nan"
         sim.numrel_mode_max = 0
         sim.numrel_mode_min = 0

    file1 = '.followup_trigger.output.xml'
    file2 = '.followup_trigger.input.xml'
    
    # write out a dummy file, after fixing the columns
    fix_numrel_columns(injections)
    grbsummary.write_rows(injections, lsctables.SimInspiralTable, file1)
  
    # call command for lalapps_sned
    # FIXME: Hardcoded lower cutoff frequency
    command = self.sned+' --f-lower 40.0  --inj-file '+file1+\
              '  --output '+file2
    if self.verbose:
      print "Executing command '%s'" % command
    subprocess.call(command.split())

    # read in the 'converted' injection
    doc = ligolw_add.ligolw_add(ligolw.Document(), [file2])
    inj_sned = lsctables.table.getTablesByName(doc, lsctables.SimInspiralTable.tableName)

    # return a single SimInspiral table
    return inj_sned[0]

  # -----------------------------------------------------
  def setTag(self, tag):
    """
    Setting a tag (non-conformal naming because of backwards compatibality)
    @param tag: the tag to be set
    """
    self.set_tag(tag)

  def set_tag(self, tag):
    """
    Setting a tag
    @param tag: well, its the tag!
    """
    self.tag = tag
    
      
  # -----------------------------------------------------
  def get_injection_window(self):
    """
    Extracting the length of the used injection_window from
    any 'FOUND' file in the cache.
    """

    # FIXME: This value should now be provided to minifollowups!!
    #        Ssipe output does come through here, but will not get the correct
    #        value (as SIRE files are not doing actual injection finding).

    # get the process params table from one of the COIRE files
    found_cache = self.cache.sieve(description = "FOUND")
    if len(found_cache)==0:
      self.injection_window = 0.1
      print "INFO: No FOUND files found, so setting the injection window"\
            " to default value."
      return
      
    coire_file = found_cache.checkfilesexist()[0].pfnlist()[0]
    try:
      doc = SearchSummaryUtils.ReadTablesFromFiles([coire_file],\
                                                   [lsctables.ProcessParamsTable])
      process_params = table.get_table(doc, lsctables.ProcessParamsTable.\
                                       tableName)
    except IOError:
      sys.stderr.write("ERROR (IOError) while reading process_params table from"\
                       " file %s. Does this file exist and does it contain"\
                       " a search_summary table?\n" %(coire_file))
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

    # set the parameter
    if self.verbose:
      print "Injection-window set to %.0f ms" % (1000*self.injection_window)


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
    
    @param filename: filename from which the injection ID is extracted
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

  def check_injection_id(self, cache_entry, tag):
      """
      The above relies on a very specific naming convention, here we check if
      the injection tag is present in the files' description.
      """
      if tag in cache_entry.description:
          return True
      else:
          return False
  
  # -----------------------------------------------------
  def find_injection_id(self, injection):
    """
    Find the injection-ID corresponding to this particular injection.
    @param injection: the injection object
    @return: the injection ID
    """

    if self.injections:
      # injection_id: the injection ID (a number or a string)
      # group_inj: a list of SimInspiral tables
      for injection_id, group_inj in self.injections.iteritems():
        for inj in group_inj:
          if injection.geocent_end_time==inj.geocent_end_time and \
                 injection.geocent_end_time_ns==inj.geocent_end_time_ns:
            return injection_id
          
      raise "No injection ID found for the above particular missed Injection "

    else:

      # cache not presearched, so searching for this particular injection
      if self.verbose:
        print "INFO: Searching for the injection at time %d.%d" % \
              (injection.geocent_end_time, injection.geocent_end_time_ns)
      injection_cache = self.cache.sieve(description = "INJECTION").\
                        sieve(ifos='HL')
      
      for file, entry in zip(injection_cache.pfnlist(), injection_cache):
        injection_id = self.get_injection_id(cache_entry = entry)
        sims = SimInspiralUtils.ReadSimInspiralFromFiles( [file], verbose=False )
        
        for sim in sims:
          # searching just by comparing the times...
          if injection.geocent_end_time==sim.geocent_end_time and \
                 injection.geocent_end_time_ns==sim.geocent_end_time_ns:
            if self.verbose:
              print "... found it: injection_id = ", injection_id
            return injection_id
      

    return None

  # -----------------------------------------------------
  def read_veto_files( self ):
    """
    Reads the veto segments given by veto-files (if any)
    """
    self.vetodict = dict()

    # loop over the IFO names
    for ifoName in self.colors.keys():

      self.vetodict[ifoName]=None
      # create the attribute name and check if it exists
      attributeName = 'followup_vetofile_'+ifoName.lower()
      if hasattr( self.opts, attributeName):

         # get the filename
         filename = getattr( self.opts, attributeName )
         if filename:
           self.vetodict[ifoName] = segmentsUtils.fromsegwizard(open(filename))

  # -----------------------------------------------------
  def reset( self ):
    """
    Resets the counting number for the time-series plots generated.
    """
    self.number=0

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
  def save_plot( self, stage ):
    """
    Saves the plots and store them in a seperate fname_list.
    @param stage: the stage this plot belongs to (e.g. INSPIRAL, THINCA,...)
    """
    fname = 'Images/'+self.opts.prefix + "_"+self.tag+"_map-"+\
            stage+"-"+str(self.number) +self.opts.suffix+'.png'
    fname_thumb = InspiralUtils.\
                  savefig_pylal( filename = self.output_path+fname,\
                                 doThumb = True,
                                 dpi_thumb = self.opts.figure_resolution)
    
    self.fname_list.append( fname )
    return fname
 

  # -----------------------------------------------------
  def get_time_trigger(self, trig):
    """
    This is a helper function to return a GPS time as one float number
    @param trig: a sngl_inspiral table entry
    """
    return float(trig.end_time) + float(trig.end_time_ns) * 1.0e-9


  # -----------------------------------------------------
  def get_effective_snr(self, trig, fac=50):
    if trig.chisq>0:
      return trig.get_effective_snr(fac=fac)
    else:
      return 0.0

  # -----------------------------------------------------
  def get_new_snr(self, trig, index=6.):
    return trig.get_new_snr(index=index)

  # -----------------------------------------------------
  def get_sim_time(self, sim, ifo = None):
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
      time = getattr(sim, ifo[0].lower()+'_end_time' )
      nano = getattr(sim, ifo[0].lower()+'_end_time_ns' )

    return  float(time) + float(nano) * 1.0e-9


  # -----------------------------------------------------
  def is_veto(self, time_trigger, ifo):
    """
    This function checks if there is a veto at the time 'timeTrigger'
    for the IFO 'ifo'.
    @param time_trigger: The time to be investigated
    @param ifo: The name of the IFO to be investigated
    """
    if self.vetodict[ifo] is None:
      return False
    
    return iterutils.any(time_trigger in seg for seg in self.vetodict[ifo])

  # -----------------------------------------------------
  def put_text(self, text):
    """
    Puts some text into an otherwise empty plot.
    @param text: text to put in the empty plot
    """
    
    newText = ''
    for i in range( int(len(text)/60.0)+1):
      newText+=text[60*i:60*i+60]+'\n'
    pylab.figtext(0.15,0.15, newText)
 
  # -----------------------------------------------------
  def create_timeseries(self, trigger_files, stage, number,\
                        slideDict=None):
    """
    Investigate inspiral triggers and create a time-series
    of the SNRs around the injected time
    @param trigger_files: List of files containing the inspiral triggers
    @param stage:        the name of the stage (INSPIRAL_FIRST, THINCA_0_CAT_2)
    @param number:       the consecutive number for this inspiral followup
    @param slideDict: A dictionary of ifo keyed slide times if using slides
    """
    # create the small and large segments for storing triggers
    seg_small =  segments.segment(self.followup_time - self.injection_window, \
                                  self.followup_time + self.injection_window)
    seg_large =  segments.segment(self.followup_time - self.time_window, \
                                  self.followup_time + self.time_window)
    if abs(seg_small) > abs(seg_large):
      err_msg = "Injection window must be smaller than time_window."
      err_msg = "Got injection window = %f and time window = %g." \
                 %(self.injection_window,self.time_window)
      raise ValueError(err_msg)
   
    # read the file(s) and get the desired sngl_inspiral rows
    if self.verbose:
      print "Processing INSPIRAL triggers from files ", trigger_files
      
    # Memory usage here can balloon, so read in files one-by-one and only keep
    # triggers within time_window
    sngls = lsctables.New(lsctables.SnglInspiralTable, \
      columns=lsctables.SnglInspiralTable.loadcolumns)

    for file in trigger_files:
      xmldoc = utils.load_filename(file, verbose=self.verbose, 
                                   contenthandler=ContentHandler)
      try:
        sngl_table = table.get_table(xmldoc,
                                     lsctables.SnglInspiralTable.tableName)
      except ValueError: # Some files have no sngl table. That's okay
        xmldoc.unlink() # Free memory
        continue
      if slideDict: # If time slide, slide the triggers
        for event in sngl_table:
          event.set_end( event.get_end() + slideDict[event.ifo] )
      # Remove triggers not within time window
      sngl_table = sngl_table.vetoed(seg_large)

      # Add to full list
      if sngl_table:
        sngls.extend(sngl_table)
      
      xmldoc.unlink() # Free memory

    # create a figure and initialize some lists
    fig=pylab.figure()
    foundSet = set()
    loudest_details = {}
    no_triggers_found = True

    if len(sngls) == 0:
      self.put_text( 'No triggers/coincidences found within time window')
    else:
      # loop over the IFOs
      for ifo in self.colors.keys():
        # get the singles for this ifo
        sngls_ifo = sngls.ifocut(ifo)

        # select the triggers within a given time window
        selected_large = sngls_ifo
        time_large = [ float(sel.get_end()) - self.followup_time \
                      for sel in selected_large ]
        selected_small = sngls_ifo.vetoed(seg_small)
        time_small = [ float(sel.get_end()) - self.followup_time \
                      for sel in selected_small ]
 
        # skip if no triggers from ifo in the large time window
        if len(time_large) == 0:
          continue
        no_triggers_found = False

        # add IFO to this set; the injection is found for this IFO-stage
        if len(time_small)>0:
          foundSet.add(ifo)

          # record details of the loudest trigger
          loudest = selected_small[selected_small.get_column('snr').argmax()]
          loudest_details[ifo] = {}
          loudest_details[ifo]["snr"] = loudest.snr
          loudest_details[ifo]["mchirp"] = loudest.mchirp
          loudest_details[ifo]["eta"] = loudest.eta
          loudest_details[ifo]["eff_dist"] = loudest.eff_distance
          loudest_details[ifo]["rchisq"] = 0
          loudest_details[ifo]["bank_rchisq"] = 0
          loudest_details[ifo]["auto_rchisq"] = 0
          if loudest.chisq_dof:
            loudest_details[ifo]["rchisq"] = loudest.chisq/(2*loudest.chisq_dof-2) 
          if loudest.bank_chisq_dof:
            loudest_details[ifo]["bank_rchisq"] = loudest.bank_chisq/loudest.bank_chisq_dof
          if loudest.cont_chisq_dof:
            loudest_details[ifo]["auto_rchisq"] = loudest.cont_chisq/loudest.cont_chisq_dof
          loudest_details[ifo]["timeTrigger"] = float(loudest.get_end())
          loudest_details[ifo]["eff_snr"] = self.get_effective_snr(loudest)
          loudest_details[ifo]["new_snr"] = self.get_new_snr(loudest)
          loudest_details[ifo]["end_time"] = loudest.end_time+loudest.end_time_ns*1E-9
          loudest_details[ifo]["trig"] = loudest

        # plot the triggers
        pylab.plot( time_large, selected_large.get_column('snr'),\
              self.colors[ifo]+'o', label="_nolegend_")
        pylab.plot( time_small, selected_small.get_column('snr'), \
              self.colors[ifo]+'s', label=ifo)

        # highlight any veto
        # FIXME: FOR NOW: COMMENTED OUT...
        ## self.highlight_veto(self.followup_time, seg_large, ifo, ylims)

      # draw the injection times and other stuff
      if no_triggers_found:
        self.put_text( 'No triggers/coincidences found within time window')
        
      ylims=pylab.axes().get_ylim()
      pylab.plot([0,0], ylims, 'g--', label="_nolegend_")
      pylab.plot([-self.injection_window, -self.injection_window], ylims, 'c:',\
            label="_nolegend_")
      pylab.plot([+self.injection_window, +self.injection_window], ylims, 'c:',\
            label="_nolegend_")

      # save the plot
      pylab.grid(True)
      pylab.legend()

    ylims=pylab.axes().get_ylim()
    pylab.axis([-self.time_window, +self.time_window, ylims[0], ylims[1]])
    pylab.xlabel('time [s]')
    pylab.ylabel('SNR')
    pylab.title(stage+'_'+str(self.number))
    fname = self.save_plot( stage )
    pylab.close(fig)

    result = {'filename':fname, 'foundset':foundSet, 'loudest_details':loudest_details}
    return result


  # -----------------------------------------------------
  def select_category(self, trigger_files, category):
    """
    Return a trigger list that contains only files for the choosen category.
    @param trigger_files : a list of trigger file names
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
  def fill_table(self, page, contents,header=False,no_wrapping=False):
    """
    Fills contents in a html table
    @param page: the page object describing a html page
    @contents: the contents of the next table item
    """

    page.add('<tr>')
    for content in contents:
      if header:
        tmpString = '<th'
      else:
        tmpString = '<td'
      if no_wrapping:
        tmpString += ' style="white-space: nowrap;"'
      tmpString += '>'
      page.add(tmpString)
      page.add( str(content) )
      if header:
        page.add('</th>')
      else:
        page.add('</td>')
    page.add('</tr>')

 
  # --------------------------------------------
  def create_table_inj(self, inj):
    """
    Creates the first table containing basic properties
    of the injection which is followed up.
    @param inj: an injection table
    """

    if self.sned:
      inj_sned = self.execute_sned([inj])[0]

    # FIXME: Add expected SNR
    # create the web-page and add a table
    page = markup.page()
    page.h1("Followup injection #"+str(self.number))
    page.add('<table border="2" >')
    page.add('<caption><b> Injection parameters </b> </caption>')
    self.fill_table( page, ['<b>parameter','<b>value'] )
    self.fill_table( page, ['Number', self.number] )
    self.fill_table( page, ['inj ID', self.injection_id] )
    self.fill_table( page, ['mass1', '%.2f' % inj.mass1] )
    self.fill_table( page, ['mass2', '%.2f' % inj.mass2] )
    self.fill_table( page, ['mtotal', '%.2f' % (inj.mass1+inj.mass2)] )
    self.fill_table( page, ['mchirp', '%.2f' % (inj.mchirp)] )
    self.fill_table( page, ['eta', '%.2f' % (inj.eta)] )
    self.fill_table( page, ['spin1z', '%.2f' % (inj.spin1z) ] )
    self.fill_table( page, ['spin2z', '%.2f' % (inj.spin2z) ] )
    self.fill_table( page, ['end_time', '%010d' % inj.geocent_end_time] )
    self.fill_table( page, ['end_time_ns', '%09d' % inj.geocent_end_time_ns] )
    self.fill_table( page, ['distance', '%.1f' % inj.distance] )
    for ifo_id in ['h','l','v','g']:
      if self.sned:
        self.fill_table( page, ['eff_dist_%s' % ifo_id, '%5.1f / %5.1f' %  \
	           (eval("inj.eff_dist_%s" % ifo_id), eval("inj_sned.eff_dist_%s" % ifo_id))] )
      else:
        self.fill_table( page, ['eff_dist_%s' % ifo_id, '%5.1f' % eval("inj.eff_dist_%s" % ifo_id)] )
    self.fill_table( page, ['playground','%s' %  pipeline.s2play(inj.geocent_end_time)] )
    page.add('</table></td>')
    
    return page
  
  # --------------------------------------------
  def create_table_sngl(self, trig):
    """
    Creates the first table containing basic properties
    of the trigger which is followed up.
    @param trig: an sngl_inspiral table
    """
    
    # create the web-page and add a table
    page = markup.page()
    page.h1("Followup trigger #"+str(self.number))
    page.add('<table border="2" >')
    self.fill_table(page, ['<b>parameter','<b>value'] )
    self.fill_table(page, ['Number', self.number] )
    self.fill_table(page, ['inj ID', self.injection_id] )
    self.fill_table(page, ['mass1', '%.2f'% trig.mass1] )
    self.fill_table(page, ['mass2', '%.2f'% trig.mass2] )
    self.fill_table(page, ['mtotal', '%.2f' % (trig.mass1+trig.mass2)] )
    self.fill_table(page, ['mchirp', '%.2f' % (trig.mchirp)] )
    self.fill_table(page, ['end_time', '%010d' % trig.end_time] )
    self.fill_table(page, ['end_time_ns', '%09d' % trig.end_time_ns] )
    self.fill_table(page, ['snr', trig.snr])
    if trig.chisq_dof:
      self.fill_table(page, ['rchisq', trig.chisq/(2*trig.chisq_dof-2)])
    if trig.bank_chisq_dof:
      self.fill_table(page, ['bank_rchisq', trig.bank_chisq/trig.bank_chisq_dof])
    if trig.cont_chisq_dof:
      self.fill_table(page, ['auto_rchisq', trig.cont_chisq/trig.cont_chisq_dof])
    self.fill_table(page, ['eff_snr (fac=50)', self.get_effective_snr(trig)])
    self.fill_table(page, ['new_snr', self.get_new_snr(trig)])
    self.fill_table(page, ['eff_distance', '%.1f' % trig.eff_distance] )
    page.add('</table></td><br>')
    page.hr()
    
    return page
  
  # --------------------------------------------
  def create_table_coinc(self, coinc,snglInspirals=None,page=None,
          slideDict=None):
    """
    Creates the first table containing basic properties
    of the coincidence which is followed up.
    @param coinc: an CoincInspiral table
    """

    if slideDict:
      timeSlide = True
    else:
      timeSlide = False
    
    ## create the web-page and add a table
    if not page:
      page = markup.page()
      page.h1("Followup trigger #"+str(self.number))
    page.add('<table border="2">')

    page.add('<caption><b>Coincidence Information</b></caption>')
    if not snglInspirals:
      self.fill_table( page, ['Statistic: ', coinc.stat] )
    else:
      self.fill_table( page, ['Combined FAR: ', coinc.combined_far] )
      self.fill_table( page, ['Uncombined FAR: ', coinc.false_alarm_rate] )
      for i in range(len(snglInspirals)):
        for j in range(i+1,len(snglInspirals)):
          sngl1 = snglInspirals[i]
          sngl2 = snglInspirals[j]
          ifo1 = sngl1.ifo
          ifo2 = sngl2.ifo
          try:
            ethinca = tools.XLALCalculateEThincaParameter(sngl1,sngl2)
          except:
            ethinca = 'Not coincident'
          Name = 'Ethinca distance between ' + ifo1 + ' and ' + ifo2
          self.fill_table( page, [Name + ': ', ethinca])

    page.add('</table><br>')

    page.add('<table border="2" >')
    page.add('<caption><b>Individual IFO Information</b></caption>')
    for ifo in ['H1','H2','L1','V1','G1']:
      trig = None
      if snglInspirals:
        for sngl in snglInspirals:
          if sngl.ifo == ifo:
            trig = sngl
            if timeSlide:
              trig2 = copy.deepcopy(trig)
              trig2.set_end(trig2.get_end() - slideDict[trig2.ifo])
      elif hasattr(coinc,ifo):
        trig = getattr(coinc,ifo)

      if trig:
        page.add('<td><table border="2" >')
    
        self.fill_table( page, ['parameter', ifo], header=True )
        self.fill_table( page, ['Number', self.number] )
        self.fill_table( page, ['inj ID', self.injection_id] )
        self.fill_table( page, ['SNR', '%.3f' % trig.snr] )
        self.fill_table( page, ['Effective SNR (fac=50)', '%.3f' % self.get_effective_snr(trig,fac=50.)] )
        self.fill_table( page, ['New SNR', '%.3f' % self.get_new_snr(trig,index=6.)] )
        if trig.chisq_dof:
          self.fill_table( page, ['Chisq/dof', '%.3f' % (trig.chisq/(2*trig.chisq_dof-2))] )
        if trig.bank_chisq_dof:
          self.fill_table( page, ['Bank chisq/dof', '%.3f' % (trig.bank_chisq/trig.bank_chisq_dof)] )
        if trig.cont_chisq_dof:
          self.fill_table( page, ['Auto chisq/dof', '%.3f' % (trig.cont_chisq/trig.cont_chisq_dof)] )
        self.fill_table( page, ['Rsq duration (s)', '%.4f' % trig.rsqveto_duration] )
        self.fill_table( page, ['''Mass1 (M<sub>&#x2A00;</sub>)''', '%.2f' % trig.mass1] )
        self.fill_table( page, ['''Mass2 (M<sub>&#x2A00;</sub>)''', '%.2f' % trig.mass2] )
        self.fill_table( page, ['''Mtotal (M<sub>&#x2A00;</sub>)''', '%.2f' % (trig.mass1+trig.mass2)] )
        self.fill_table( page, ['''Mchirp (M<sub>&#x2A00;</sub>)''', '%.3f' % trig.mchirp] )
        self.fill_table( page, ['Template duration (s)', '%.3f' % trig.template_duration ] )
        if timeSlide:
          endTime = trig.end_time + 1E-9*trig.end_time_ns
          self.fill_table( page, ['Slid GPS end time', '%.4f' % endTime] )
          slidEndTime = trig2.end_time + 1E-9*trig2.end_time_ns
          self.fill_table( page, ['Unslid end time', '%.4f' % slidEndTime] )
        else:
          endTime = trig.end_time + 1E-9*trig.end_time_ns
          self.fill_table( page, ['GPS end time', '%.3f' % endTime] )
        self.fill_table( page, ['Effective distance (Mpc)', '%.1f' % trig.eff_distance] )
        page.add('</table></td>')

    page.add('</table><br>')
    page.hr()
    
    return page
  
  # --------------------------------------------
  def create_table_time(self, trigger_time):
    """
    Creates the first table containing the time
    of the followup
    @param trigger_time: well, just the trigger time
    """
    
    ## create the web-page and add a table
    page = markup.page()
    page.h1("Followup time around GPS "+str(trigger_time) )
    page.add('<table border="2" >')
    self.fill_table( page, ['Number', self.number] )
    self.fill_table( page, ['Time', trigger_time] )
    page.add('</table></td><br>')
    page.hr()

    return page


  # --------------------------------------------
  def add_table_followup(self, page, invest_dict):
    """
    Adds the table containing specific information on the loudest
    trigger found in the followup region for each IFO.
    @param page: the html page to which to add the table
    @param invest_dict: dictionary containing the stage results
    """
    
    ## print out the result for this particular injection
    page.add('<table border="2" >')
    page.add('<caption><b> Parameters of the loudest (by SNR) recovered single ifo triggers at each stage of the pipeline </b> </caption>')
    self.fill_table( page, ['step', 'F/M', 'SNR', \
                            'Mchirp', 'eta', 'eff_dist', \
                            'rchisq', 'bank_rchisq', 'auto_rchisq', 'eff_snr',\
                            'new_snr', 'end_time', 'ethinca', 'Veto ON/OFF'], header=True )

    # loop over the stages and create the table with
    # the various data in it (when available)
    for stage in self.orderLabels:
      if stage in invest_dict:
        result = invest_dict[stage]

        # Fill in the details of the loudest found coinc.
        found_ifo = ''
        loudest_snr = ''
        loudest_mchirp = ''
        loudest_eta = ''
        loudest_eff_dist = ''
        loudest_rchisq = ''
        loudest_bank_rchisq = ''
        loudest_auto_rchisq = ''
	loudest_effsnr = ''
        loudest_newsnr = ''
        loudest_ethinca = ' '
        loudest_time = ' '
        veto_onoff = ''

        # add all the IFO's for this coincident
        result['foundlist'] = list(result['foundset'])
        for i in range(len(result['foundlist'])):
          ifo = (result['foundlist'])[i]
          found_ifo += ifo+' '
          
          # Parameters of the loudest trigger, taken from the
          # 'loudest-details' dictionary, created in 'create_timeseries'
	  loudest_snr += "%s : %.3f <br>" % \
                         (ifo, result['loudest_details'][ifo]['snr'])
	  loudest_mchirp += "%s : %.3f <br>" % \
                         (ifo, result['loudest_details'][ifo]['mchirp'])
          loudest_eta += "%s : %.3f <br>" % \
                         (ifo, result['loudest_details'][ifo]['eta'])
	  loudest_eff_dist += "%s : %.3f <br>" % \
                         (ifo, result['loudest_details'][ifo]['eff_dist'])
	  loudest_rchisq += "%s : %.3f <br>" % \
                         (ifo, result['loudest_details'][ifo]['rchisq'])
          loudest_bank_rchisq += "%s : %.3f <br>" % \
                         (ifo, result['loudest_details'][ifo]['bank_rchisq'])
          loudest_auto_rchisq += "%s : %.3f <br>" % \
                         (ifo, result['loudest_details'][ifo]['auto_rchisq'])
	  loudest_effsnr += "%s : %.3f <br>" % \
                         (ifo, result['loudest_details'][ifo]['eff_snr'])
	  loudest_newsnr += "%s : %.3f <br>" % \
                         (ifo, result['loudest_details'][ifo]['new_snr'])
          loudest_time += "%s : %.3f <br>" % \
                         (ifo, result['loudest_details'][ifo]['end_time'])
          for j in range(i+1,len(result['foundlist'])):
            ifo2 = (result['foundlist'])[j]
            try:
              ethinca = tools.XLALCalculateEThincaParameter(
                         result['loudest_details'][ifo]['trig'],
                         result['loudest_details'][ifo2]['trig'])
              loudest_ethinca += "%s and %s: %.3f <br>" % \
                         (ifo,ifo2,ethinca)
            except:
              loudest_ethinca += "%s and %s: %s <br>" % \
                         (ifo,ifo2,'Not coincident')
          
          # Check whether some of the ifo times is vetoed
          time_trigger = float(result['loudest_details'][ifo]['timeTrigger'])
          if self.vetodict[ifo]:
            veto = self.is_veto(time_trigger, ifo)
            veto_txt = 'OFF'
            if veto:
              veto_txt = 'ON'
            veto_onoff+=ifo+': '+veto_txt+'<br>'
          else:
            veto_onoff+=ifo+': No info<br>'

        # Fill the table whether something is found or not
        if len(result['foundset'])>0:
          self.fill_table( page, [ stage,  'FOUND in <br>'+found_ifo, \
                                   loudest_snr, \
                                   loudest_mchirp, \
                                   loudest_eta, \
                                   loudest_eff_dist,\
                                   loudest_rchisq, \
                                   loudest_bank_rchisq, \
                                   loudest_auto_rchisq, \
                                   loudest_effsnr, \
                                   loudest_newsnr, \
                                   loudest_time, \
                                   loudest_ethinca, \
				   veto_onoff],no_wrapping=True)
        else:
          self.fill_table( page, [ stage,  '<font color="red">MISSED'])
          
    page.add('</table>')
    page.add('</td></tr></table><br><br>')

    return page


  # -----------------------------------------------------
  def from_coinc(self, coinc, ifo = None, more_infos = False, \
                 injection_id = None,slideDict=None):
    """
    Creates a followup page from a coincident trigger.
    @param coinc: the coincidence to be followed up
    @param ifo: specifies the ifo to be used from the coinc.
    @param more_infos: to have some additional informations
    @param injection_id: Must be specified for exttrig search
                         to specify what injection to use
    """

    if not ifo:
      print "WARNING (not bad): No IFO specified, using data from the "\
            "first IFO in the coincidence "
      for ifo_name in self.colors.keys():
        if hasattr(coinc, ifo_name):
          ifo = ifo_name
          break
    sngl = getattr(coinc, ifo)
    
    # set the time
    self.followup_time = float(sngl.get_end())
 
    # prepare the page
    self.injection_id = injection_id
    page =  self.create_table_coinc(coinc,slideDict=slideDict)
    self.flag_followup = more_infos

    # do the followup
    return self.followup(page,slideDict=slideDict)

  # -----------------------------------------------------
  def from_new_coinc(self, coinc, sngls,\
                 more_infos = False, injection_id = None):
    """
    Creates a followup page from a coincident trigger.
    @param coinc: the coincidence to be followed up
    @param ifo: specifies the ifo to be used from the coinc.
    @param more_infos: to have some additional informations
    @param injection_id: Must be specified for exttrig search
                         to specify what injection to use
    """

    sngl = sngls[0]

    # set the time
    self.followup_time = float(sngl.get_end())

    # prepare the page
    self.injection_id = injection_id
    page =  self.create_table_coinc(coinc,snglInspirals= sngls)
    self.flag_followup = more_infos

    # do the followup
    return self.followup(page)

  # -----------------------------------------------------
  def from_new_slide_coinc(self, coinc, sngls,slideDict,\
                 more_infos = False, injection_id = None):
    """
    Creates a followup page from a slid coincident trigger. This function
    does not yet produce the plots (as I'm not sure how to!) but the
    relevant information to do this (the slide dictionary and the segment list)
    are provided to this function.
    @param coinc: the coincidence to be followed up
    @param ifo: specifies the ifo to be used from the coinc.
    @param more_infos: to have some additional informations
    @param injection_id: Must be specified for exttrig search
                         to specify what injection to use
    """

    sngl = sngls[0]

    # set the time
    self.followup_time = float(sngl.get_end())

    # prepare the page
    self.injection_id = injection_id
    page =  self.create_table_coinc(coinc,snglInspirals= sngls,\
        slideDict=slideDict)
    self.flag_followup = more_infos

    return self.followup(page,slideDict=slideDict)

  # -----------------------------------------------------
  def from_sngl(self, sngl, ifo = None, more_infos = False, \
                injection_id = None):
    """
    Creates a followup page from a single trigger.
    @param sngl: the sngl trigger to be followed up
    @param ifo: NOT USED
    @param more_infos: to have some additional informations
    @param injection_id: Must be specified for exttrig search
                         to specify what injection to use
    """

    # prepare the page
    self.injection_id = injection_id
    page =  self.create_table_sngl(sngl)
    self.flag_followup = more_infos

    # set the time
    self.followup_time = float(sngl.get_end())

    # do the followup
    return self.followup(page)
    
  # -----------------------------------------------------
  def from_missed(self, missed, ifo = None, more_infos = True, \
                  injection_id = None):
    """
    Creates a followup page from a missed injection.
    @param sngl: the missed injection to be followed up
    @param ifo: The ifo whose time is used (geocent if None)
    @param more_infos: to have some additional informations
    @param injection_id: Must be specified for exttrig search
                         to specify what injection to use
    """

    return self.from_injection(missed, ifo = ifo, more_infos = more_infos,\
                               injection_id = injection_id )
    
  # -----------------------------------------------------
  def from_found(self, found, ifo = None, more_infos = False, \
                 injection_id = None,coinc = None, sngls = None):
    """
    Creates a followup page from a found injection.
    @param sngl: the found injection to be followed up
    @param ifo: The ifo whose time is used (geocent if None)
    @param more_infos: to have some additional informations
    @param injection_id: Must be specified for exttrig search
                         to specify what injection to use
    """

    return self.from_injection(found, ifo = ifo, more_infos = more_infos, \
                               injection_id = injection_id,coinc=coinc, \
                               sngls = sngls )
    
  # -----------------------------------------------------
  def from_injection(self, injection, ifo = None, more_infos = True, \
                     injection_id = None,coinc = None, sngls = None):
    """
    Creates a followup page from an injection.
    @param injection: the injection to be followed up
    @param ifo: The ifo whose time is used (geocent if None)
    @param more_infos: to have some additional informations
    @param injection_id: Must be specified for exttrig search
                         to specify what injection to use
    """

    # Set the injection ID if required
    if injection_id:
      self.injection_id = injection_id
    else:
      self.injection_id = self.find_injection_id(injection)

    # prepare the page
    page =  self.create_table_inj(injection)

    if coinc and sngls:
      page =  self.create_table_coinc(coinc,snglInspirals= sngls,page=page)

    self.flag_followup = more_infos
    
    # set the time and do the followup
    self.followup_time = self.get_sim_time(injection, ifo)
    
    # do the followup
    return self.followup(page)


  # -----------------------------------------------------
  def from_time(self, trigger_time, ifo = None, more_infos = False, \
                injection_id = None):
    """
    Creates a followup page from a given time.
    @param trigger_time: the time to be followed up
    @param ifo: NOT USED
    @param injection_id: Must be specified for exttrig search
                         to specify what injection to use
    """

    self.flag_followup = more_infos

    # prepare the page
    page =  self.create_table_time(trigger_time)

    # set the time
    self.followup_time = trigger_time
    self.injection_id = injection_id

    # do the followup
    return self.followup(page)
    

  # -----------------------------------------------------
  def followup(self, page,slideDict = None):
    """
    Central followup procedure, finding corresponding files,
    generating the time-series and creating the output html files
    @param page: The head of the html page created with different informations
    @param slideDict: A dictionary of ifo keyed slide times if using slides
    @return: filename of the created html
    """
  
    # increase internal number:
    self.number+=1
    page.add('<br>')

    # loop over each stage
    invest_dict = {}
    for stage, cache in self.trigger_cache.iteritems():

      # loop over each file in a stage
      trig_cache = lal.Cache()
      for c in cache:

        # check the time and the injection ID
        # Also pick up files +/- 2048s of this trigger to avoid boundary issues
        if (self.followup_time in c.segment) or ((self.followup_time-2048) in c.segment) or ((self.followup_time+2048) in c.segment):
          if not self.injection_id or \
                 (self.injection_id and \
                  self.check_injection_id(c, self.injection_id)):
            trig_cache.append( c )
        
      # check if the pfnlist is empty. `
      file_list = trig_cache.pfnlist()
      if len(file_list)==0:
        print >>sys.stderr, "ERROR: No files found for stage %s in the "\
              "cache for ID %s and time %d; probably mismatch of a "\
              "pattern in the options. " % \
              ( stage, self.injection_id, self.followup_time)
        continue

      # call the function to create the timeseries
      if stage in ('COINCIDENCE','COINCIDENCE-SLID'):
        # ... need to loop over the four categories
        for cat in [1,2,3,4,5]:
          select_list=self.select_category(file_list, cat)
          if len(select_list)==0:
            print "WARNING (not that bad): "\
                  "No THINCA files found for category ", cat
            continue          
          modstage = stage+'_CAT_' + str(cat)
          invest_dict[modstage] = self.create_timeseries(select_list,modstage,\
                                       self.number,slideDict)
      else:
        invest_dict[stage]=self.create_timeseries(file_list, stage, \
                                self.number,slideDict)


    ## add some more followup if required
    if self.flag_followup:
      self.add_table_followup(page, invest_dict)

    ## add the pictures to the webpage
    for stage in self.orderLabels:
      if stage in invest_dict:
        result = invest_dict[stage]
      
        fname = result['filename']
        page.a(extra.img(src=[fname], width=400, \
                         alt=fname, border="2"), title=fname, href=[fname])


    ## add the version of this code
    page.add("<hr>")
    page.add("Figure(s) and data produced with " + __prog__ + ", version " \
              + git_version.verbose_msg)
        
    # and write the html file
    htmlfilename = self.opts.prefix + "_followup_"+str(self.number) +\
                         self.opts.suffix+'.html'
    file = open(self.opts.output_path+htmlfilename,'w')
    file.write(page(False))
    file.close()

    # store html file in fname_list and return filename
    self.fname_list.append(htmlfilename)
    return htmlfilename

  # -----------------------------------------------------
  def find_unvetoed_ifos(self, inj, ifo_list):
    """
    Find the list of unvetoed injections for this particular
    injection (missed or found)
    @param inj: the injection (SimInspiral table)
    @param ifo_list: the ifo list
    """

    unvetoed_ifos = []
    for ifo in ifo_list:
      veto_time = float(getattr(inj, ifo[0].lower()+'_end_time'))
      if not self.is_veto(veto_time, ifo):
        unvetoed_ifos.append(ifo)

    return unvetoed_ifos

  # -----------------------------------------------------
  def get_relevant_distance(self, inj, ifo_list):
    """
    Find the relevant distance for this injection,
    given  the list of IFO's
    """

    # take the second largest effective distance
    list_dist = []
    for ifo in ifo_list:
      list_dist.append(getattr(inj, 'eff_dist_'+ifo[0].lower()))
    list_dist.sort()

    return list_dist[1]
  
  # -----------------------------------------------------
  def populate_followup_list(self, inj, ifo_list):
    """
    Find the relevant effective distance for this injection
    including any possible vetoes.
    """

    # retreive list of unvetoed IFOs
    unvetoed_ifos = self.find_unvetoed_ifos(inj, ifo_list)

    if len(unvetoed_ifos)>=2:
      # if there are 2 or more unvetoed IFO's left
      # a coincidence is possible. Need to followup this case

      # find the type (for plotting issues)
      type = 'ks'
      if len(unvetoed_ifos) == len(ifo_list):
        # no veto for any ifo
        type = 'ro'

      # get the relevant distance
      rel_dist = self.get_relevant_distance(inj, unvetoed_ifos)

      # append to the list of potential injections to be followed up
      self.followup_dict['inj'].append(inj)
      self.followup_dict['dist'].append(rel_dist)
      self.followup_dict['type'].append(type)

    else:
      self.vetoed_injections.append(inj)

  # -----------------------------------------------------
  def followup_missed_injections(self, opts, missed_injections, \
                                 found_injections, \
                                 xvalue = 'mtotal', xlog = 'linear', \
                                 ylog = 'log'):
    """
    This function will do the complete followup on each relevant
    missed injection.
    """


    # FIXME: put this dict into a basic class
    foundSymbol = {'H1H2':'x', 'H1L1':'^', 'H2L1':'+', 'H1V1':'v', \
                   'H2V1':'<','L1V1':'>', 'H1H2L1':'s', 'H1H2V1':'D',\
                   'H2L1V1':'d', 'H1L1V1':'1',\
                   'H1H2L1V1':'p', 'H1':'x', 'H2':'x', 'L1':'x', 'V1':'x'}
  
    # populate the injections to be followed up
    ifo_list =  get_ifo_list(opts.ifo_times)
    for inj in missed_injections:
      self.populate_followup_list(inj, ifo_list)

    # select the first N injections to followup
    if opts.followup_number:
      index = numpy.argsort(self.followup_dict['dist'])
      number_to_followup = min(len(index), opts.followup_number)
    else:
      raise NotImplementedError,"NotImplemented. You only can specify the "\
            "number of injections to be followed up. "

    # Handle the found injections. Distinguish them by the
    # IFO's they were found with
    found_dist = {}
    found_masses = {}
    for type in found_injections.keys():
      temp_dist = None

      # create a list of IFO's from the one-word-list
      sub_ifo_list = get_ifo_list(type)

      # get the distance and the masses
      found_dist[type] = [self.get_relevant_distance(inj, sub_ifo_list)\
                          for inj in found_injections[type]]
      
      found_masses[type], legend = get_data_from_table( \
        found_injections[type], xvalue)
      
    ## create the main image
    pylab.clf()

    # plot the found injections
    for type in found_injections.keys():

      # create a list of IFO's from the one-word-list
      sub_ifo_list = get_ifo_list(type)
      
      # Color the markers blue iff this inj is found in all
      # possible detectors
      col = 'm'
      if len(sub_ifo_list) == len(ifo_list):
        col = 'b'

      # FIXME: foundSymbol not defined
      pylab.plot( found_masses[type], found_dist[type], col+foundSymbol[type], \
            markerfacecolor='None',markeredgecolor=col, \
            label = type, markersize=10, markeredgewidth=1)

    # don't forget to plot the real vetoed ones
    dist = [self.get_relevant_distance(inj, ifo_list)\
            for inj in self.vetoed_injections]
    masses =  viz.readcol(self.vetoed_injections,"mass1")+\
             viz.readcol(self.vetoed_injections,"mass2")
    pylab.plot( masses, dist, 'ks', \
                markerfacecolor='None',label='vetoed',\
                markeredgecolor='k',markersize=10, markeredgewidth=1)

    # plot the missed ones
    missed_masses, legend = get_data_from_table( self.followup_dict['inj'], xvalue)

    pylab.plot( missed_masses, self.followup_dict['dist'], 'ro', \
          markerfacecolor='None',label='(prop) missed',\
          markeredgecolor='r',markersize=10, markeredgewidth=1)


    legy = 'Second smallest eff. distance [Mpc]'
    title_text = legend+' vs '+legy
    pylab.title( opts.title + ' '+title_text+' in '+opts.ifo_times+\
                 ' times', size='x-large')
    pylab.xlabel(legend, size='x-large')
    pylab.ylabel(legy, size='x-large')
    pylab.grid(True)
    pylab.legend()

    # create a mapDist
    mapDict ={'object':'', \
              'text':'Click on a filled circle to go to the followup page',\
              'xCoords':[], 'yCoords':[], 'links': []}

    # coordinates of the actual picture in the image
    # FIXME: These numbers might change for different solutions,
    # and image types.
    boundFigX = [100.0,720.0]
    boundFigY = [540.0, 60.0]

    # choose lin/log for the axes
    pylab.gca().set_xscale(xlog)
    pylab.gca().set_yscale(ylog)

    # limits on the axes for the image
    axes=pylab.gcf().axes[0]
    rangeFigX = axes.get_xlim()
    rangeFigY = numpy.log10( axes.get_ylim() )
  
    # calculating scaling factors
    slopeX = (boundFigX[1]-boundFigX[0])/(rangeFigX[1]-rangeFigX[0])
    interX = boundFigX[1] - slopeX*rangeFigX[1]
    slopeY = (boundFigY[1]-boundFigY[0])/(rangeFigY[1]-rangeFigY[0])
    interY = boundFigY[1] - slopeY*rangeFigY[1]
    

    # loop over all missed injections
    for ind in index[:number_to_followup]:
      
      # get the relevant data for this injection to be followed up
      inj = self.followup_dict['inj'][ind]
      type = self.followup_dict['type'][ind]
      eff_dist = self.followup_dict['dist'][ind]
      eff_dist_log = numpy.log10(eff_dist)
      total_mass = inj.mass1+inj.mass2


      # get the coordinates of this missed injection
      px = int(interX+slopeX*total_mass)
      py = int(interY+slopeY*eff_dist_log)
      
      # highlight the clickable missed ones; depending if it is
      # a real missed one or a possible missed one (with vetoes)
      pylab.plot( [total_mass], [eff_dist], type, markerfacecolor = type[0], \
            markersize=10)

      # do the followup and create the followup page
      followuphtml = self.from_missed(inj)
      
      # add the point and the link to the mapDict
      mapDict['xCoords'].append(px)
      mapDict['yCoords'].append(py)
      mapDict['links'].append(followuphtml)

    return mapDict
  
#######################################################################
def get_data_from_table(table, col_name, ifo = None):
  """
  Retrieves data from a table, including non-table entries
  such as 'mtotal' and 'time'
  @param table: the table containing the data
  @param col_name: the name of the column
  @param ifo: the ifo for which the column is returned
  """

  if col_name=='time':
     data = [t for t in viz.timeindays(viz.readcol(table,"end_time", ifo))]
     legend = "End time (in days)"
  elif col_name == 'mtotal':
     data = viz.readcol(table,"mass1")+viz.readcol( table,"mass2")
     legend = "Total mass"
  else:
     data = viz.readcol(table, col_name)
     legend = xname

  return data, legend

#######################################################################
def getData( table, xname, yname, ifo  ):
  """
  Retrieves data from a table, including non-table entries
  such as 'mtotal' and 'time'
  @param table: the table with the data
  @param xname: the x-value
  @param yname: the y-value
  @param ifo  : the ifo for which the time or distance is returned
  """

  if xname=='time':
    xp= [ t-opts.time_offset \
          for t in viz.timeindays(viz.readcol( table,"end_time", ifo)) ]
    legx = "End time (in days)"
  elif xname == 'mtotal':
    xp = viz.readcol( table,"mass1")+viz.readcol( table,"mass2")
    legx = "Total mass"
  else:
    xp = viz.readcol( table,xname)
    legx = xname
  

  if yname == 'eff_dist':
    yp = viz.readcol( table,"eff_dist", ifo )
    legy = "Effective distance"
  else:
    yp = viz.readcol( table,yname)
    legy = yname

  return xp, yp, legx, legy

####################################################
def get_ifo_list(ifo_times):
  ifo_list = []
  for i in range(len(ifo_times)/2):
    ifo_list.append(ifo_times[2*i:2*i+2])
  return ifo_list
