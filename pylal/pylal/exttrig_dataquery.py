#!/usr/bin/python

import glob
import os
import sys
import subprocess
import shutil
import optparse
import urllib
import ConfigParser
import warnings
warnings.simplefilter('ignore',DeprecationWarning)
warnings.simplefilter('ignore',UserWarning)

itertools = __import__("itertools")  # absolute import of system-wide itertools

import numpy
import matplotlib; matplotlib.use('Agg')
import pylab as plt

from glue.ligolw import ligolw
from glue.ligolw import utils
from glue.ligolw import table
from glue.ligolw import lsctables
from glue import lal
from glue import segments
from glue import segmentsUtils
from glue import pipeline
from glue import iterutils

##############################################################################
# set-up
##############################################################################

def parse_args():
  """
  Parses the command line arguments.
  """
  parser = optparse.OptionParser()
  # files
  parser.add_option("--config-file", help="Specifies the configuration file to use.")
  parser.add_option("--log-file", help="Specifies the log file for any outputs.",default = "default.log")
  parser.add_option("--grb-file", help="Full path to the GRB xml file.")
  # GRB information
  parser.add_option("--time", help="GPS time of the GRB trigger.")
  parser.add_option("--ra",  help="Right ascension of the GRB (degree).")
  parser.add_option("--dec", help="Declination of the GRB (degree).")
  parser.add_option("--name", help="Name of the GRB (e.g. 100122A).")
  # segment search specs
  parser.add_option("--offset", help="Length of time to check in both directions [s].",type=int,default = 2500)
  parser.add_option("--useold", action="store_true", default=False, help="If this option is specified, the old old segment files will be used. Otherwise they will be queries freshly.")
  parser.add_option("--extend", action="store_true", default=False, help="Adds quanta of adjacent segments of data with duration of segment-duration/2.")
  # outputs
  parser.add_option("--make-plots", action = "store_true", default = False, help = "If this option is specified, segment plots will be created.")
  parser.add_option("--make-xml", action = "store_true", default = False, help = "If this option is specified, a xml file will be created.")

  options, arguments = parser.parse_args()

  # checks
  if not options.config_file:
    raise ValueError, "Must enter a configuration file to use."
  if not options.grb_file and (not options.time or not options.name):
    raise ValueError, "Either a valid GRB xml file must be specified or the GPS time and name of the GRB!"
  if options.make_xml:
    if not options.ra or not options.dec:
      raise ValueError, "The right ascension (--ra) and the declination (--dec) is required to be specified, because --make-xml was set."
  if options.grb_file and options.make_xml:
    print "Warning! It does not make sense to specfy the xml file as input and as output. Therefore, the make-xml flag will be set to False"
    options.make_xml = False

  return options, arguments

def plot_segments(segdict, onsource, offsource, centertime, plot_window, output_filename, tag):
  """
  Creates time series plots of segments for each IFO.
  """
  # define hardcoded variables
  color_code = {'H1':'r', 'H2':'b', 'L1':'g', 'V1':'m', 'G1':'k'}

  # create the plot
  fig = plt.figure()
  ax  = fig.add_subplot(111)
  ax.set_xlabel("time (s)")
  ax.set_ylabel("IFO")
  window = None

  # converts from GPS time to relative +/- centertime
  _time_transform = lambda t: t - centertime

  # add content
  ifolist = segdict.keys()
  ifolist.sort()
  for row, ifo in enumerate(ifolist):
    color = color_code[ifo]
    for seg in segdict[ifo]:
      a = _time_transform(seg[0])
      b = _time_transform(seg[1])
      ax.fill([a, b, b, a, a], [row, row, row+1, row+1, row], color)

  # set x limit of plot
  if offsource:
    c = _time_transform(plot_window[0])
    d = _time_transform(plot_window[1])
    window = segments.segment((c, d))

    # plot vertical lines for offsource
    if abs(plot_window) > 0:
      a = _time_transform(offsource[0])
      b = _time_transform(offsource[1])
      ax.axvline(a, color='k', linewidth=2)
      ax.axvline(b, color='k', linewidth=2)

  # plot vertical lines for onsource
  ax.axvline(_time_transform(onsource[0]), color='k', linestyle='--')
  ax.axvline(_time_transform(onsource[1]), color='k', linestyle='--')

  # plot ticks
  ticks = plt.arange(len(ifolist)) + 0.5
  ax.set_yticks(ticks)
  ax.set_yticklabels(ifolist)

  # set y axis as IFOs
  if window is not None:
    ax.set_xlim(window)
  ax.set_ylim((0, len(ifolist)))

  # save file
  ax.set_title('Segments for GRB '+tag)
  fig.savefig(output_filename)
  plt.close(fig)

def check_segment_availability(grb_name, grb_time, query_start, query_end, offset, ifo, segmentName):
  '''
  Searches +/- offset from GRB time to download the latest segment lists then extracts times and puts them into a txt file.
  '''
  args = {'grb_name'    : grb_name,
          'query_start' : query_start,
          'query_end'   : query_end,
          'ifo'         : ifo,
          'segmentName' : segmentName}
  cmd  = "ligolw_segment_query --database --query-segments --include-segments '{segmentName}' --gps-start-time {query_start} --gps-end-time {query_end} > ./segments{ifo}_grb{grb_name}.xml".format(**args)
  print '>>',cmd
  print
  process    = subprocess.Popen([cmd], shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
  output,err = process.communicate()

  # try to open the file
  try:
    doc = utils.load_filename("segments{ifo}_grb{grb_name}.xml".format(**args), contenthandler = lsctables.use_in(ligolw.LIGOLWContentHandler))
  except:
    raise IOError, "Error reading file: segments{ifo}_grb{grb_name}.xml".format(**args)

  # extract the segment list from segment:table and store in a txt file
  segs = table.get_table(doc, "segment")
  seglist = segments.segmentlist(segments.segment(s.start_time, s.end_time) for s in segs)
  segmentsUtils.tosegwizard(file("{ifo}-science_grb{grb_name}.txt".format(**args),'w'),seglist,header = True)

  print ">> %s segments +/-%ds from %ds found:"%(ifo,offset,grb_time)
  for s in segs:
    print "Start:",s.start_time,"End:",s.end_time,"Duration:",s.end_time-s.start_time
  print

  return

def exttrig_dataquery(grb_name, grb_time, grb_ra, grb_dec, offset, config_file, extend=False, useold=False, make_plots=False, make_xml=False):
  '''
  Finds science time of all available IFOs.
  '''
  ##############################################################################
  # get segment duration and minimum amount of science time
  ##############################################################################

  # read the configuration file
  cp = ConfigParser.ConfigParser()
  cp.read(config_file)

  # define hardcoded variables
  basic_ifolist    = ifolist = ['H1','H2','L1','V1']
  catlist          = [1,2,3]
  sensitivity_dict = {"H1": 1, "L1": 2, "H2": 3, "V1": 4, "G1": 5}

  # get segment length from configuration file
  pad_data = int(cp.get('data','pad-data'))
  if cp.has_option('data','segment-duration'):
    blockDuration = segmentDuration = psdDuration = int(cp.get('data','segment-duration'))
  elif cp.has_option('data','segment-length'):
    blockDuration = segmentDuration = psdDuration = int(cp.get('data','segment-length')) /int(cp.get('data','sample-rate'))
  else:
    raise ValueError, "EXIT: Cannot find segment-duration in [data] section of configuration file!"

  # get sample rate
  if cp.has_option('data','sample-rate'):
    sampleRate = int(cp.get('data', 'sample-rate'))
    print ">> Sample rate has been set to: %d"%sampleRate
    print
  else:
    print ">> ERROR: Need to specify sample-rate in [data] section of configuration file in order to calculate inputs for downstream processes."
    sys.exit()

  # if not extend option then need to get block duration
  if not extend:
    if cp.has_option('data','block-duration'):
      blockDuration = int(cp.get('data','block-duration'))
    elif cp.has_option('data','segment-length'):
      s_length  = int(cp.get('data', 'segment-length'))
      s_num     = int(cp.get('data', 'number-of-segments'))
      s_rate    = int(cp.get('data', 'sample-rate'))
      s_overlap = int(cp.get('inspiral', 'segment-overlap'))
      # calculate blockDuration
      blockDuration = ( s_length * s_num - ( s_num - 1 ) * s_overlap ) / s_rate
    else:
      raise ValueError, "EXIT: Cannot find block-duration in [data] section of configuration file! Either set block-duration or use --extend option."

  # calculate the minimum amount of science time need and how the length of quanta to be added on both ends of the analysis time
  minscilength = blockDuration + 2 * pad_data
  quanta       = segmentDuration / 2

  # if extend beyond minscilength; add segments of quanta length to each end of segment
  print ">> Minimum science segment length is: %ss"%minscilength
  print
  if extend:
    print ">> Will extend minimum science segment by quanta of: %ss"%quanta
    print

  ##############################################################################
  # get list of segments for each IFO and put in science txt file
  ##############################################################################

  if not useold:
    # external call to ligolw_segment_query
    query_start = int(grb_time - offset)
    query_end   = int(grb_time + offset)
    for ifo in ifolist:
      if cp.has_option('segments','%s-segments'%ifo.lower()):
	segmentName = cp.get('segments','%s-segments'%ifo.lower())
	check_segment_availability(grb_name, grb_time, query_start, query_end, offset, ifo, segmentName)

  ##############################################################################
  # get veto files
  ##############################################################################

  if not useold:
    # get and read veto definer file
    veto_file_url          = cp.get('exttrig','cvs_veto_definer')
    veto_file_path,headers = urllib.urlretrieve(veto_file_url,os.path.basename(veto_file_url))

    # do ligolw_segments_from_cats
    deltat     = 500
    args       = {'start_time'     : int(grb_time - offset - deltat),
		  'end_time'       : int(grb_time + offset + deltat),
		  'veto_file_path' : veto_file_path}
    cmd        = "ligolw_segments_from_cats --database --veto-file={veto_file_path} --separate-categories --gps-start-time {start_time}  --gps-end-time {end_time} --output-dir=. --individual-results".format(**args)
    print '>>',cmd
    print
    process    = subprocess.Popen([cmd], shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output,err = process.communicate()

    # Rename the veto files for easier handling
    veto_files = glob.glob('./*VETOTIME_CAT*{start_time}*xml'.format(**args))
    for filename in veto_files:
      p = filename.split('-')
      newname = "%s-%s_grb%s.xml"%(p[0], p[1], grb_name)
      shutil.move(filename, newname)

  ##############################################################################
  # look in txt files and find segment with onsource and minscilength
  ##############################################################################

  # create segment that is +/- offset of GRB time
  onsource        = [grb_time - int(cp.get('exttrig','onsource_left')),\
		     grb_time + int(cp.get('exttrig','onsource_right'))]
  onSourceSegment = segments.segment(onsource[0], onsource[1])

  # get segments in science txt files; see if segments length at least minscilength
  # if no then discard them; if yes then put in segdict[ifo] and ifo in ifolist
  basic_segdict = segdict = segments.segmentlistdict()
  for ifo in ifolist:
    # check configuration file
    if not cp.has_option('segments','%s-segments' % ifo.lower()):
      continue
    # find segment with onsource and check it is at least minscilength
    ifo_segfile = '%s-science_grb%s.txt' %(ifo, grb_name)
    if os.path.exists(ifo_segfile):
      tmplist = segmentsUtils.fromsegwizard(open(ifo_segfile))
      try:
	s = tmplist.find(onSourceSegment)
      except ValueError:
	# if onsource not in segments then move on to next IFO
	continue
      if abs(tmplist[s]) >= minscilength:
	segdict[ifo] = segments.segmentlist([tmplist[s]])
      basic_segdict[ifo] = segments.segmentlist([s for s in tmplist])
  ifolist = segdict.keys()

  if len(ifolist) < 2:
    print "EXIT: Less than 2 interferometers have available data!"
    sys.exit()

  ##############################################################################
  # apply vetoes
  ##############################################################################

  # apply
  print ">> Vetoes that overlap with science segments:"
  for ifo in ifolist:
    # flag; True if IFO not vetoed
    cat_flag = True
    for cat in catlist:
      # create list and check for overlaps
      xmlsegfile = "./%s-VETOTIME_CAT%s_grb%s.xml" %(ifo, cat, grb_name)
      if os.path.exists(xmlsegfile) and cat_flag:
	testseg = segments.segment([segdict[ifo][0][0],segdict[ifo][0][1]])
	list_overlaps = []

	# load the content of the veto-file
        xmldoc = utils.load_filename(xmlsegfile, gz = False, contenthandler = lsctables.use_in(ligolw.LIGOLWContentHandler))
	segs = table.get_table(xmldoc, lsctables.SegmentTable.tableName)
	segdefs = table.get_table(xmldoc, lsctables.SegmentDefTable.tableName)

	# create a mapping between the segments and their definitions
	defdict = {}
	for segdef in segdefs:
	  defdict[segdef.segment_def_id] = segdef.name

	# find veto segments that intersect science segment of IFO with onsource
	for seg in segs:
	  s = segments.segment(seg.start_time, seg.end_time)
	  if testseg.intersects(s):
	   id = seg.segment_def_id
	   list_overlaps.append([defdict[id], seg.start_time, seg.end_time])

	# cut veto CAT1 segments out of science segment; CAT1,2,3 veto IFO if in onsource will veto IFO
	for name, segstart, segend in list_overlaps:
	  print "CAT%s IFO %s, Start: %d End: %d because %s"%(cat, ifo, segstart, segend, name)
	  s = segments.segment(segstart, segend)
	  if onSourceSegment.intersects(s):
	    segdict.pop(ifo, None)
	    cat_flag = False
	    break
	if cat == 1:
	  vetoes = segments.segmentlist(segments.segment(s[1], s[2]) for s in list_overlaps)
	  segdict[ifo] -= vetoes

  # get list of IFOs
  ifolist = segdict.keys()

  print

  if len(ifolist) < 2:
    print "EXIT: After vetoes, less than 2 interferometers have available data!"
    sys.exit()

  ##############################################################################
  # determine segment to be analyzed
  ##############################################################################

  # sort from most sensitive to least sensitive
  def sensitivity_cmp(ifo1, ifo2):
    return cmp(sensitivity_dict[ifo1], sensitivity_dict[ifo2])
  ifolist.sort(sensitivity_cmp)

  # compares IFOs and finds the segment to analyze
  # now try getting off-source segments
  # start trying with all IFOs
  # work our way through subsets; beginning with most sensitive combinations
  test_combos = itertools.chain(*itertools.imap(lambda n: iterutils.choices(ifolist, n), xrange(len(ifolist), 1, -1)))
  off_source_segment = None
  the_ifo_combo = []
  for ifo_combo in test_combos:
    # find conincident science time of IFOs
    trial_seglist = segdict.intersection(ifo_combo)
    if abs(trial_seglist) < minscilength:
      print "EXIT: IFOs do not overlap enough for minscilength",abs(trial_seglist)
      sys.exit()
    else:
      pass

      # find segment with grb_time inside
      try:
	super_seg = trial_seglist[trial_seglist.find(onSourceSegment)].contract(pad_data)
      except ValueError:
	print "EXIT: ValueError with super_seg"
	sys.exit()
      if onSourceSegment not in super_seg:
	print "EXIT: onsource not in super_seg"
	sys.exit()

      # find int division of onsource time intervals before and after grb 
      tplus  = (super_seg[1] - onSourceSegment[1])
      tminus = (onSourceSegment[0] - super_seg[0])

      # get minimum number of onsource time intervals in offsource
      tmin = ( minscilength - 2*pad_data - abs(onSourceSegment) )

      # cut to get minscilength
      if tplus + tminus > tmin:
	half_max = tmin // 2
	if tplus < half_max:
	  print ">> Left sticks out so cut it."
	  remainder = tmin - tplus
	  tminus = min(remainder, tminus)
	elif tminus < half_max:
	  print ">> Right sticks out so cut it."
	  remainder = tmin - tminus
	  tplus = min(remainder, tplus)
	else:
	  print ">> Both sides stick out so cut as symmetrically as possible."
	  tminus = half_max
	  tplus = tmin - half_max  # odd trial sticks out on right
      if tplus + tminus < tmin:
	 offsource = None
      temp_segment = segments.segment((onSourceSegment[0] - tminus - pad_data, onSourceSegment[1] + tplus + pad_data))

      if temp_segment is not None:
	offsource = temp_segment
	ifolist   = list(ifo_combo)

	if extend:
	  # extend as many adjacent 128 second blocks as possible
	  begin_time = offsource[0] - quanta * (abs(super_seg[0]-offsource[0])//quanta)
	  end_time   = offsource[1] + quanta * (abs(super_seg[1]-offsource[1])//quanta)
	  offsource  = segments.segment((begin_time,end_time))

	break
  print

  # check length at least minscilength
  if abs(offsource) < minscilength:
    print abs(offsource),minscilength
    print "EXIT: Calculated offsource segment but less than minscilength!"
    sys.exit()

  # check if no detectors can be used then exit
  if len(ifolist) < 2:
    print "EXIT: Calculated offsource segment but less than two IFOs!"
    sys.exit()

  # check edge case
  if abs(offsource[0]-onsource[0]) < pad_data or abs(offsource[1]-onsource[1]) < pad_data:
    print "WARNING: GRB time close to edge of offsource. Its within the padding time."

  # concatenate "H1L1V1", etc.
  ifolist.sort()
  ifotag = "".join(ifolist)
  print ">> Offsource segment for %s GRB is:"%ifotag
  print "Start:",offsource[0],"End:",offsource[1],"Duration:",offsource[1]-offsource[0],"Left:",grb_time-offsource[0],"Right:",offsource[1]-grb_time
  print


  ##############################################################################
  # output
  ##############################################################################

  # write analyse txt files
  for ifo in basic_ifolist:
    if ifo in ifolist:
      analysisFP = open('%s-analyse_grb%s.txt' %(ifo,grb_name),'w')
      analysisFP.write('# seg\t start    \t stop    \t duration\n')
      analysisFP.write('0\t %d\t %d\t %d\n' %(offsource[0],offsource[1],offsource[1]-offsource[0]))
    else:
      analysisFP = open('%s-analyse_grb%s.txt' %(ifo,grb_name),'w')
      analysisFP.write('# seg\t start    \t stop    \t duration\n')

  # calculate blockDuration
  blockDuration = int(abs(offsource[0]-offsource[1])) - 2 * pad_data

  # calculate psdDuration
  # gets largest power of two such that blockDuration/psdDuration = psdRatio
  # could have done a binary & operator that is faster but this is more user-friendly I believe
  min_psdDuration = int(cp.get('exttrig', 'min-psd-length'))
  psdRatio        = int(cp.get('exttrig', 'psd-ratio'))
  psdDuration = 2**int(numpy.log2(blockDuration/psdRatio))
  if psdDuration < min_psdDuration:
    print "EXIT: PSD segment duration is too short. It is %ds but needs to be at least %ds in length."%(psdDuration,min_psdDuration)
    sys.exit()

  # some downstream processes (e.g. lalapps_tmpltbank) cannot handle these inputs
  if cp.has_option('data', 'segment-duration'):
    cp.remove_option('data', 'segment-duration')
    cp.remove_option('data', 'block-duration')

  # some downstream processes (e.g. lalapps_tmpltbank) requires these options to run
  print ">> Using sample rate of %d to calculate inputs for downstream processes."%sampleRate
  print
  segmentLength  = segmentDuration*sampleRate
  segmentCount   = blockDuration/(segmentDuration/2) - 1 # subtract 1 because one segment length is overlapped
  segmentOverlap = segmentLength/2
  cp.set('data', 'segment-length', segmentLength)
  cp.set('data', 'number-of-segments', segmentCount)
  cp.set('inspiral', 'segment-overlap', segmentOverlap)

  # set values for [coh_PTF_inspral] section in configuration file
  cp.set('coh_PTF_inspiral', 'block-duration', blockDuration)
  cp.set('coh_PTF_inspiral', 'segment-duration', segmentDuration)
  cp.set('coh_PTF_inspiral', 'psd-segment-duration', psdDuration)
  cp.set('coh_PTF_inspiral', 'pad-data', pad_data)
  f = open('grb%s.ini'%grb_name,'w')
  cp.write(f)
  f.close()
  print ">> The [data] section of the configuration file has been edited with the following values:"
  print "sample-rate=",sampleRate
  print "segment-length=",segmentLength
  print "number-of-segments=",segmentCount
  print "segment-overlap=",segmentOverlap
  print
  print ">> The [coh_PTF_inspiral] section of the configuration file has been edited with the following values:"
  print "block-duration =",blockDuration
  print "segment-duration =",segmentDuration
  print "psd-segment-duration =",psdDuration
  print "pad-data =",pad_data
  print

  # plot segments
  offSourceSegment = segments.segment(offsource[0], offsource[1])
  plot_window      = segments.segment(grb_time-offset, grb_time+offset)
  plot_segments(basic_segdict, onSourceSegment, offSourceSegment, grb_time, plot_window, "segment_plot_%s.png"%grb_name, grb_name)

  # make xml file
  if make_xml:
    # create a new xml document with an ExtTriggers Table
    xmldoc = ligolw.Document()
    xmldoc.appendChild(ligolw.LIGO_LW())
    tbl = lsctables.New(lsctables.ExtTriggersTable)
    xmldoc.childNodes[-1].appendChild(tbl)

    # set the values we need
    row = lsctables.ExtTriggersTable()
    row.process_id = None
    row.det_alts = None
    row.det_band = None
    row.det_fluence = None
    row.det_fluence_int = None
    row.det_name = None
    row.det_peak = None
    row.det_peak_int = None
    row.det_snr = ''
    row.email_time = 0 
    row.event_dec = float(grb_dec)
    row.event_dec_err = 0.0
    row.event_epoch = ''
    row.event_err_type = ''
    row.event_ra = float(grb_ra)
    row.event_ra_err = 0.0
    row.start_time = grb_time
    row.start_time_ns = 0
    row.event_type = ''
    row.event_z = 0.0
    row.event_z_err = 0.0 
    row.notice_comments = '' 
    row.notice_id = ''
    row.notice_sequence = ''
    row.notice_time = 0
    row.notice_type = ''
    row.notice_url = ''
    row.obs_fov_dec = 0.0
    row.obs_fov_dec_width = 0.0
    row.obs_fov_ra = 0.0
    row.obs_fov_ra_width = 0.0
    row.obs_loc_ele = 0.0
    row.obs_loc_lat = 0.0
    row.obs_loc_long = 0.0
    row.ligo_fave_lho = 0.0
    row.ligo_fave_llo = 0.0
    row.ligo_delay = 0.0
    row.event_number_gcn = 9999
    row.event_number_grb = grb_name
    row.event_status = 0

    # insert into the table and write file
    tbl.extend([row])
    filename = 'grb%s.xml' % grb_name
    utils.write_filename(xmldoc, filename)

  # plot all vetoes
  if make_plots:
    vetodict = segments.segmentlistdict()
    for cat in catlist:
      for ifo in ifolist:
	vetofile = "%s-VETOTIME_CAT%s_grb%s.xml" % (ifo, cat, grb_name)
        xmldoc = utils.load_filename(vetofile, gz = False, contenthandler = lsctables.use_in(ligolw.LIGOLWContentHandler))
	segs = table.get_table(xmldoc, lsctables.SegmentTable.tableName)
	segdefs = table.get_table(xmldoc, lsctables.SegmentDefTable.tableName)
	vetodict[ifo] = segments.segmentlist(segments.segment(s.start_time, s.end_time) for s in segs)

      if vetodict:
	plot_segments(vetodict, onSourceSegment, offSourceSegment, grb_time, plot_window, "veto_plot_CAT%s_%s.png"%(cat,grb_name), "%s CAT%s"%(grb_name, cat))

  # return
  return 'grb%s.ini'%grb_name, ifolist, onSourceSegment, offSourceSegment

if __name__ == "__main__":
  # get the command line arguments
  opts,args = parse_args()

  # get time, RA, DEC and name of GRB; get offset to search from GRB time
  if opts.grb_file:
    xmldoc    = utils.load_filename(opts.grb_file, gz=opts.grb_file.endswith('.gz'), contenthandler = lsctables.use_in(ligolw.LIGOLWContentHandler))
    ext_table = table.get_table(xmldoc,lsctables.ExtTriggersTable.tableName)
    grb_time  = ext_table[0].start_time
    grb_name  = os.path.basename(opts.grb_file)[3:-4]
    grb_ra    = ext_table[0].event_ra
    grb_dec   = ext_table[0].event_dec
  else:
    grb_name = opts.name
    grb_time = int(opts.time)
    grb_ra   = float(opts.ra)
    grb_dec  = float(opts.dec)

  # run
  exttrig_dataquery(grb_name, grb_time, grb_ra, grb_dec, opts.offset, opts.config_file, opts.extend, opts.useold, opts.make_plots, opts.make_xml)
