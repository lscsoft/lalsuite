"""
Utilities for the inspiral plotting functions
"""

from glue import lal
from glue import segments
import socket, os
import sys
import copy
import math

from glue.ligolw import ligolw
from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw import utils
from pylal import SnglInspiralUtils
from pylal import CoincInspiralUtils
from pylal import git_version
from glue import markup
from glue.markup import oneliner as extra_oneliner

# set default color code for inspiral plotting functions
colors = {'G1':'k','H1':'r','H2':'b','L1':'g','V1':'m'}
symbols = {'G1':'Y','H1':'x','H2':'o','L1':'+','V1':'1'}

# set color codes for coincident ifo types
def get_coinc_ifo_colors( ifo_set ):
  """
  Given an ifo set, returns an html color code for plotting.
  """
  # check that ifo_set is set or frozenset
  if not ( isinstance(ifo_set, set) or isinstance(ifo_set, frozenset) ):
    raise ValueError("ifo_set must be of type set or frozenset. "
      "Use lsctables.instrument_set_from_ifos to do this.")

  if ifo_set == set(['H1', 'H2', 'L1', 'V1']):
    return '#F88017' # dark orange
  elif ifo_set == set(['H1', 'H2', 'L1']):
    return '#00FFFF' # cyan
  elif ifo_set == set(['H1', 'L1', 'V1']):
    return '#7D1B7E' # dark orchid
  elif ifo_set == set(['H2', 'L1', 'V1']):
    return '#153E7E' # dodger blue4
  elif ifo_set == set(['H1', 'L1']):
    return '#00FF00' # green
  elif ifo_set == set(['H1', 'V1']):
    return '#6698FF' # sky blue
  elif ifo_set == set(['H2', 'L1']):
    return '#FF0000' # red
  elif ifo_set == set(['H2', 'V1']):
    return '#FF00FF' # magenta
  elif ifo_set == set(['L1', 'V1']):
    return '#254117' # dark green
  else: # other coincs just set to black
    return 'k'


############################################################
class InspiralPage(object):
  """
  This is a class to contain all the bits of a inspiral page
  showing the results of a piece of code.
  """

  # ------------------------------------------------------
  # ------------------------------------------------------
  def __init__(self, options, ifo_times=None, ifo_tag=None, user_tag=None,
               gps_start_time=None, gps_end_time=None):
    """
    Initializes this class with the options.
    @params options: option object from the calling function
    """
    self.opts = options

    # set the version of the code and the caling executable
    self.version = git_version.verbose_msg.replace('\n','<br>')
    self.name = os.path.basename(sys.argv[0])
    self.arguments = sys.argv[1:]

    # use the values from this option object to initialize some variables
    self.initialize()

    # overwrite the default values if needed
    if ifo_times: self.ifo_times = ifo_times
    if ifo_tag: self.ifo_tag = ifo_tag
    if user_tag: self.user_tag = user_tag
    if gps_start_time is not None: self.gps_start_time = gps_start_time
    if gps_end_time is not None: self.gps_end_time = gps_end_time

    # and NOW create the suffix and prefixes
    self.create_affixes()

    # create some empty lists
    self.fname_list = []
    self.tag_list = []
    self.html_footer = ""


  # ------------------------------------------------------
  def initialize(self):
    """
    Extract information from the option structure.
    Does NOT alter the information in the option structure
    """

    if hasattr(self.opts, 'ifo_times'):
      self.ifo_times = self.opts.ifo_times
    else:
      self.ifo_times = None

    if hasattr(self.opts, 'ifo_tag'):
      self.ifo_tag = self.opts.ifo_tag
    else:
      self.ifo_tag = None

    if hasattr(self.opts, 'user_tag'):
      self.user_tag = self.opts.user_tag
    else:
      self.user_tag = None

    if hasattr(self.opts, 'gps_start_time'):
      self.time_string = str(int(self.opts.gps_start_time))+"-" + \
                         str(int(math.ceil(self.opts.gps_end_time))-int(self.opts.gps_start_time))
    else:
      self.time_string = "unspecified-gpstime"

    if hasattr(self.opts,'output_path'):
      self.output_path = self.opts.output_path
      if not self.output_path.endswith('/'):
        self.output_path += '/'
    else:
      self.output_path = './'

    # create the output directory later

  # ------------------------------------------------------
  def create_affixes(self):
    """
    Create the affixes (prefix/suffix) for the image naming
    """
    self.prefix = self.name

    if self.ifo_times:
      self.prefix = self.ifo_times + "-" + self.prefix
    if self.ifo_tag:
      self.prefix = self.prefix + "_" + self.ifo_tag
    if self.user_tag:
      self.prefix = self.prefix + "_" + self.user_tag

    self.suffix = "-"+self.time_string


  # ------------------------------------------------------
  def add_plot(self, plot_fig, fig_description, output_dir = None):
    """
    Add a plot to the page.
    @param plot_fig: handle of the figure
    @param fig_description: descriptive figure text for the filename
    @param output_dir: alternate output directory [optional]
    """

    fname = "Images/" + self.prefix + "_"+ fig_description + self.suffix

    if output_dir:
      fname = output_dir + '/' + fname
    else:
      fname = self.output_path + fname

    filename, filename_thumb = self.savefig(fname, plot_fig)

    self.fname_list.append(filename)
    self.tag_list.append(filename)

  # ------------------------------------------------------
  def savefig(self, filename_base, fig, doThumb=True, dpi = None, dpi_thumb=50):
    """
    Function to create the image file. 
    @param filename_base: basename of the filename (without the .png ending) 
    @param fig: handle to the figure to save
    @param doThumb: save the thumbnail or not (doThumb=True by default)
    @param dpi: resolution of the figure
    @param dpi_thumb: resolution of the thumbnail (dpi=50 by default)
    """

    savefig_kwargs = {}
    if dpi is not None:
      savefig_kwargs["dpi"] = dpi

    # create the image file
    filename = filename_base + '.png'
    fig.savefig(filename, **savefig_kwargs)

    if doThumb:
      filename_thumb = filename_base + '_thumb.png'
      fig.savefig(filename_thumb, dpi=dpi_thumb)
    else:
      filename_thumb = None

    return filename, filename_thumb

  # ------------------------------------------------------
  def write_page(self, infix = None, doThumb = True, 
                 map_list = [], coinc_summ_table = None ):
    """
    Create the pages if output is enabled
    """
    if self.opts.enable_output:
      html_filename = self.create_htmlname(infix)
      self.write_html_output(html_filename, doThumb = doThumb,
                             map_list = map_list, coinc_summ_table = coinc_summ_table,
                             comment=self.html_footer or None)
      self.write_cache_output(html_filename)
      return html_filename

  # ------------------------------------------------------
  def write(self, text):
    """
    Write some text to the standard output AND
    to the page.
    """
    print text
    self.html_footer+=text+'<br>'

  # ------------------------------------------------------
  def create_htmlname(self, infix):
    """
    Create the html filename
    """

    # Any infix to use?
    if infix:
      html_filename = self.prefix + '_' + infix + self.suffix
    else:
      html_filename = self.prefix + self.suffix

    html_filename += ".html"

    if self.output_path:
      html_filename = self.output_path + html_filename

    return html_filename

  # ------------------------------------------------------
  def write_html_output(self, html_filename, doThumb=True, map_list=[],
                        comment=None, coinc_summ_table=None ):
    """
    @param doThumb: Uses the thumbnail file as the sourcs for the images
    @param map_list: A list of dictionaries to create the image maps
    @param comment: A comment that can be added to the page
    @param coinc_summ_table: A CoincSummTable that can be added to the page
    """

    # Initialise the html output file
    page = markup.page()
    try:
      page.init(title=__title__)
    except:
      page.init()

    page.h1(self.name + " results")

    page.p(self.prefix + self.suffix)
    page.hr()

    # open the output file
    html_file = file(html_filename, "w")

    # loop over the contents
    for tag,filename in zip(self.tag_list, self.fname_list):

      # set the correct name for linking (two '//' does not bother)
      fname = "Images/" + os.path.basename(filename)

      # set the thumbnail pictures if required
      if doThumb:
        fname_thumb = fname[:-4] + "_thumb.png"
      else:
        fname_thumb = fname

      # add the image to the page
      page.a(extra_oneliner.img(src=[fname_thumb], width=400,
          alt=tag, border="2"), title=tag, href=[ fname])

    page.add("<hr/>")

    # add maps to this page
    m=0
    for map_dict in map_list:
      m+=1
      page.add( map_dict['text']+'<br>' )
      page.add( '<IMG src="%s" width=800px '
                'usemap="#map%d">' % ( map_dict['object'], m) )
      page.add( '<MAP name="map%d"> <P>' % m )
      n=0
      for px, py, link in zip( map_dict['xCoords'],
                               map_dict['yCoords'],
                               map_dict['links'] ):
        n+=1
        page.add( '<area href="%s" shape="circle" '
                  'coords="%d, %d, 5"> Point%d</a>' %
                  ( link, px, py, n) )
      page.add('</P></MAP></OBJECT><br>')
      page.add("<hr/>")

    # add some extra stuff if needed
    if comment:
      page.add("<div> "+comment+"</div>")
      page.hr()

    if coinc_summ_table:
      page.add(coinc_summ_table)
      page.hr()

    text = self.write_process_params()
    page.add(text)
    html_file.write(page(False))
    html_file.close()


  # ------------------------------------------------------
  def write_cache_output(self, html_filename):
    """
    Write the output cache file of the plotting functions.
    @param: html_filename: the name of the html file
    """

    output_cache_name = self.prefix + self.suffix +'.cache'
    if self.output_path:
      output_cache_name = self.output_path + output_cache_name

    # open the cachefile
    cachefile = open(output_cache_name, 'w')
    cachefile.write(os.path.basename(html_filename) + '\n')

    # loop over each entry in fname_list
    for filename in self.fname_list:
      if filename.endswith('.png'):
        fname = "Images/"+os.path.basename(filename) # set the correct name for linking
      elif filename.endswith('.html'):
        fname = os.path.basename(str(filename)) # set the correct name for linking

      # add the name to the cache file
      cachefile.write(fname + '\n')

    # and close the file
    cachefile.close()

  # ------------------------------------------------------
  def write_process_params(self):
    """
    Returns the version and the full command run 
    """
    text = "Figure(s) produced with '" + self.name + "' with version: <br>" \
           + self.version \
           + '<br>\n<p style="width:80%; color:blue">' + self.name
    for arg in self.arguments:
      text += " " + arg
    text += '</p>'

    return text

def savefig_pylal(filename=None, filename_thumb=None, doThumb=True, dpi=None,
  dpi_thumb=50, fig=None):
  """
  @param filename: filename in which to save the figure
  @param filename_thumb: filename into which save a thumbnail of the figure
  @param doThumb: save the thumbnail or not (doThumb=True by default)
  @param dpi: resolution of the figure
  @param dpi_thumb: resolution of the thumbnail (dpi=50 by default)
  @param fig: the particular figure you wish to save (current figure by
         default)
  @return filename_thumb if a thumbnail was created (computed from filename
          by default)

  """
  import pylab

  # fill in non-trivial defaults
  if fig is None:
    fig = pylab.gcf()
  if dpi is None:
    dpi = pylab.rcParams["savefig.dpi"]
  if doThumb and (filename_thumb is None):
    if filename is None:
      raise ValueError("must provide filename_thumb or filename if doThumb "
        "is True")
    index = filename.rindex('.')
    filename_thumb = filename[0:index] + '_thumb' + filename[index:]

  # save picture into a file
  if filename is not None:
    fig.savefig(filename, dpi=dpi)

  # save thumbnail into a file if requested
  if doThumb:
    fig.savefig(filename_thumb, dpi=dpi_thumb)

  return filename_thumb


def ErrorMessagePlotting(opts, thisplot):
  """

  """
  text = "---Error in "+opts.name+"in plotting functions "+thisplot
  if "chi" in thisplot:
    text += "\n---possible reasons related to chi-square (are you reading first stage triggers ?)"
  print >> sys.stderr, text


def message(opts, text):
  """

  """
  if opts.verbose:
    print text
  return text+'<br>\n'

def set_figure_tag( plot_description, datatype_plotted = '', open_box = True ):
  """
  Returns a string containing a standard naming convention for the tag part of a figure name.
  The convention used is:
  ${plot_description}_${DATATYPE}_PLOTTED(_OPEN_BOX)

  @plot_description: a description of the plot
  @datatype_plotted: datatype that appears in the plot. If specified, '_PLOTTED' will be
   added after the datatype
  @open_box: whether or not looking at the plot constitutes opening the box
   If set to True, _OPEN_BOX will be added to the file name. If set False, no flag
   will be added.
  """
  if datatype_plotted != '':
    datatype_plotted = '_'.join([ '', datatype_plotted, 'PLOTTED' ])
  if open_box is True:
    box_flag = '_OPEN_BOX'
  else:
    box_flag = ''

  return ''.join([ plot_description, datatype_plotted, box_flag ])

def set_figure_name(opts, figure_tag):
  """
  return a string containing a standard output name for pylal 
  plotting functions.
  """
  fname = ''.join([ "Images/", opts.prefix, "_", figure_tag, opts.suffix, ".png" ])

  if opts.output_path is not None:
    fname = opts.output_path + fname

  return fname

def write_coinc_summ_table(tableList = [], commentList = [], stat=None, statTag=None, number=None, format=None,followup = None, followupOpts = None):
  """
  picks out loudest coincident triggers from given CoincInspiralUtils Tables
  and returns info about the coincidences in a html or wiki table 

  @param tableList: a list of CoincInspiralUtils.coincInspiralTables
  @param commentList: comments about each table (e.g., file name)
  @param stat: any CoincInspiralUtils.coincStatistic
  @param statTag: string specifying what stat used
  @param number: number of triggers to list
  @param format: desired output format; can be either 'html' or 'wiki'
  """
  # set format
  if format == 'html':
    tx = '<table border = "1">'
    xt = '</table>'
    thx = '<tr><td colspan=14>'
    rx = '<tr><td>'
    xr = '</td></tr>'
    xccx = '</td><td>'
  elif format == 'wiki':
    tx = ''
    xt = ''
    thx = '||<-14>'
    rx = '||'
    xr = '||\n'
    xccx = '||'
  else:
    raise ValueError, 'unrecognized format; must be either html or wiki'

  # set statTag if not specified
  if statTag is None: statTag = stat.name

  CoincSummTable = ''

  # populate table
  for coincTable, coincComment in zip(tableList,commentList):
    if stat.name == 'far':
      coincTable.sort(descending=False)
    else:
      coincTable.sort()
    rank = 1
    # set table header
    CoincSummTable = CoincSummTable + tx + thx + coincComment + xr
    CoincSummTable = CoincSummTable + \
        rx + ' Rank ' + xccx + ' followup ' + xccx + 'Coinc IFOs' + xccx + \
        statTag + xccx + 'False Alarm Probability' + xccx + ' end_time ' + \
        xccx + ' end_time_ns ' + xccx + ' mass1 ' + xccx + ' mass2 ' + xccx + ' mchirp ' + \
        xccx + ' eta ' + xccx + ' snr ' + xccx + ' chisq ' + xccx + ' effective_snr ' + xr
    for coinc in coincTable:
      if format == 'html':
        CoincSummTable = CoincSummTable + '<tr><td rowspan=' + str(coinc.numifos) + '>' + str(rank) + '</td>'
      elif format == 'wiki':
        CoincSummTable = CoincSummTable + rx + '<|' + str(coinc.numifos) + '>' + str(rank) + xccx
      followupLink = 'None'
      if followup:
        followup.from_coinc( coinc, coinc.get_ifos()[1][0] )
        coinc.get_ifos()[1][0]
        followupFile = followupOpts.prefix \
            + '_followup_' + str(followup.number) + followupOpts.suffix \
            + '.html'
        followupLink = '<a href="./' + followupFile +'"> here </a>'
      if format == 'html':
        CoincSummTable = CoincSummTable + '<td rowspan=' + str(coinc.numifos) + '>' + followupLink + xccx
      elif format == 'wiki':
        CoincSummTable = CoincSummTable + rx + '<|' + str(coinc.numifos) + '>' + followupLink + xccx
      # cycle through info
      for trig in coinc:
        CoincSummTable = CoincSummTable + trig.ifo + xccx + str(coinc.stat) + xccx + str(coinc.fap) + xccx + str(trig.end_time) + \
                xccx + str(trig.end_time_ns) + xccx + str(trig.mass1) + xccx + str(trig.mass2) + xccx + str(trig.mchirp) + \
                xccx + str(trig.eta) + xccx + str(trig.snr) + xccx + str(trig.chisq) + xccx + str(trig.get_effective_snr()) + xr + \
                rx
      CoincSummTable = CoincSummTable + xr
      rank = rank + 1
      if rank > number: break
    CoincSummTable = CoincSummTable + xt

  return CoincSummTable

def write_html_output(opts, args, fnameList, tagLists,
      doThumb=True, mapList = [],
      comment=None, CoincSummTable=None,
      html_tag = '', add_box_flag=False):
  """
  @param opts: The options from the calling code
  @param args: The args from the calling code
  @param fnameList: A list of the filenames
  @param tagLists: A list for the tags, getting added to the links
  @param doThumb: Uses the _thumb file as the sourcs for the images
  @param mapList: A list of dictionaries to create the image maps
  @html_tag: tag to add to html filename
  @add_box_flag: Adds _OPEN_BOX to the html file name if any
   of the files in filelist have "_OPEN_BOX" in their name. Otherwise,
   will add "_CLOSED_BOX" to the file name. These flags go between
   opts.prefix and opts.suffix
  """

  prefix = opts.prefix
  # add the html_tag if desired
  if html_tag != '':
    prefix += '_' + html_tag
  # add the box-flag to the prefix if desired
  if add_box_flag:
    box_flag = ''
    if any(fname for fname in fnameList if 'OPEN_BOX' in fname):
      box_flag ='_OPEN_BOX'
    else:
      box_flag = '_CLOSED_BOX'
    # add the box flag to the prefix
    prefix += box_flag

  # -- the HTML document and output cache file
  # -- initialise the web page calling init_page
  page, extra = init_markup_page(opts)
  page.h1(opts.name + " results")

  page.p(prefix + opts.suffix)
  page.hr()

  # -- filename
  html_filename = prefix + opts.suffix +".html"
  if opts.output_path:
    html_filename = opts.output_path + html_filename
  html_file = file(html_filename, "w")

  # loop over the contents
  for tag,filename in zip(tagLists,fnameList):

    # set the correct name for linking (two '//' does not bother)
    fname = "Images/" + os.path.basename(filename)

      # set the thumbnail pictures if required
    if doThumb:
      fname_thumb = fname[:-4] + "_thumb.png"
    else:
      fname_thumb =fname

    # add the image to tge page
    page.a(extra.img(src=[fname_thumb], width=400,
        alt=tag, border="2"), title=tag, href=[ fname])

  page.add("<hr/>")

  # add maps to this page
  if len(mapList)>0:
    m=0
    for mapDict in mapList:
      m+=1
      page.add( mapDict['text']+'<br>' )
      page.add( '<IMG src="%s" width=800px '
                'usemap="#map%d">' % ( mapDict['object'], m) )
      page.add( '<MAP name="map%d"> <P>' % m )
      n=0
      for px, py, link in zip( mapDict['xCoords'],
                               mapDict['yCoords'],
                               mapDict['links'] ):
        n+=1
        page.add( '<area href="%s" shape="circle" '
                  'coords="%d, %d, 5"> Point%d</a>' %
                  ( link, px, py, n) )
      page.add('</P></MAP></OBJECT><br>')
      page.add("<hr/>")

  if opts.enable_output:
    if comment is not None:
      page.add("<div> "+comment+"</div>")
      page.hr()
    if CoincSummTable is not None:
      page.add(CoincSummTable)
      page.hr()
    text = writeProcessParams(opts.name, opts.version, args)
    page.add(text)
    html_file.write(page(False))
    html_file.close()

  return html_filename


def write_cache_output(opts, html_filename,fnameList):
  """
  write the output cache file of theplotting functions
  """

  output_cache_name = '.'.join([html_filename.rstrip('.html'), 'cache'])
  this = open(output_cache_name, 'w')
  if opts.enable_output:
    this.write(os.path.basename(html_filename) + '\n')
  for filename in fnameList:
    if str(filename).endswith('.png'):
      fname = "Images/"+os.path.basename(filename) # set the correct name for linking
    elif str(filename).endswith('.html'):
      fname = os.path.basename(str(filename)) # set the correct name for linking
    this.write(fname + '\n')
  this.close()


def writeProcessParams(name, version, command):
  """
  Convert input parameters from the process params that the code was called 
  with into a formatted string that can be saved within an other document 
  (e.g., HTML)

  @param name: name of the executable/script
  @param version:version of the executable/script
  @param command: command line arguments from a pylal script
  @return text
  """
  text = "Figure(s) produced with '" + name + "' with version: <br>" \
      + version \
      + '<br>\n<p style="width:80%; color:blue">'+ name
  for arg in command:
    text += " " + arg
  text+='</p>'

  return text

def AddFileToCache(fname, cache):
  """
  Add the given file to the lal.Cache

  @param fname:
  @param cache:
  """
  file_name = fname.split('.')[0].split('-')
  cache.append(lal.CacheEntry( file_name[0], file_name[1],
    segments.segment(int(file_name[2]),
    int(file_name[2]) + int(file_name[3])),
    'file://' + socket.gethostbyaddr(socket.gethostname())[0] + \
    os.getcwd() + '/' + fname))

def GenerateCache(fileList):
  """
  Generate a lal.Cache for the list of files

  @param fileList : a list of files
  @return cache
  """
  cache = lal.Cache()
  for file in fileList:
    AddFileToCache(file, cache)
  return(cache)


class SummValueContentHandler(ligolw.PartialLIGOLWContentHandler):
  """
  Content handler that only reads in the SummValue table
  """
  def __init__(self, xmldoc):
    ligolw.PartialLIGOLWContentHandler.__init__(self, xmldoc, lambda name, attrs: lsctables.IsTableProperties(lsctables.SummValueTable, name, attrs))

try:
  lsctables.use_in(SummValueContentHandler)
except AttributeError:
  # old glue did not allow .use_in().
  # FIXME: remove when we can require the latest version of glue
  pass


def initialise(opts, name, version = None):
  """
  Create suffix and prefix that will be used to name the output files.
  'version' is outdated and not used anymore.

  @param opts : the user arguments (user_tag, gps_end_time and 
  gps_start_time are used).
  @param name: name of the calling function/executable
  @return prefix 
  @return suffix
  """


  # compose prefix
  prefix = name
  try:
    if opts.ifo_times:
      prefix = opts.ifo_times +"-"+ prefix
  except:
    print >> sys.stderr, "--ifo-times option not implemented in the "+name +" executable. skipping..."
    pass
  try:
    if opts.ifo_tag:
      prefix = prefix + "_" + opts.ifo_tag
  except:
    print >> sys.stderr, "--ifo-tag option not implemented in the "+name +" executable. skipping..."
    pass
  try:
    if opts.user_tag:
      prefix = prefix + "_" + opts.user_tag
  except:
    print >> sys.stderr, "--user-tag option not implemented in the "+name +" executable. skipping..."
    pass

  # compose suffix
  try:
    if opts.gps_start_time is not None and opts.gps_end_time is not None:
      suffix = "-"+str(int(opts.gps_start_time))+"-"+str(int(math.ceil(opts.gps_end_time))-int(opts.gps_start_time))
    else:
      suffix = "-unspecified-gpstime"
  except:
    suffix = "-unspecified-gpstime"
    print >> sys.stderr, "--gps-start-time and/or --gps-end-time option not implemented in the " + \
                         name + " executable. skipping..."
    pass

  opts.prefix = prefix
  opts.suffix = suffix
  opts.name = name
  opts.version = git_version.verbose_msg.replace("\n","<br>")

  # make sure output_path is set correctly
  if opts.output_path is not None:
    if not opts.output_path.endswith("/"):
      opts.output_path += "/"

    # create output file if required
    if not os.path.exists(opts.output_path):
      os.mkdir(opts.output_path)

    if not os.path.exists(opts.output_path+"Images"):
      os.mkdir(opts.output_path+"Images")

  else:
    if not os.path.exists("Images"):
      os.mkdir("Images")

  return opts


def init_markup_page( opts):
  """
  Load the markup module, and initialise the HTML document if the opts 
  argument contains enable_ouput option.

  @param  opts : the user arguments 
  @return page 
  @return extra 
  """
  # Initialise the html output file
  if opts.enable_output:
    try:
      from glue import markup
      from glue.markup import oneliner as extra_oneliner
    except:
      raise ImportError("Require markup.py to generate the html page")

    page = markup.page()
    try:
      page.init(title=__title__)
    except:
      page.init()

  return page, extra_oneliner


def readHorizonDistanceFromSummValueTable(fList, verbose=False, contenthandler=SummValueContentHandler):
  """
  read in the SummValueTables from a list of files and return the
  horizon distance versus total mass

  @param fList:   list of input files
  @param verbose: boolean (default False)
  """

  output = {}
  massOutput = {}
  count = 0
  if len(fList) == 0:
    return output

  # for each file in the list
  for thisFile in fList:
    if verbose:
      print str(count+1)+"/"+str(len(fList))+" " + thisFile
    count = count+1
    massNum = 0

    doc = utils.load_filename(thisFile, contenthandler = contenthandler)
    try:
      summ_value_table = table.get_table(doc, lsctables.SummValueTable.tableName)
    except ValueError:
      print "ValueError in readHorizonDistanceFromSummValueTable while reading summvalue table from file ", thisFile
      return output,massOutput

    # if not summ_value table was filled , then simply returns
    if summ_value_table is None:
      return output,massOutput

    # else
    for row in summ_value_table:
      # we should find a name "inspiral_effective_distance"
      if row.name == 'inspiral_effective_distance':
        # it may be that the file read is an inspiral file containing only the BNS infomration
        if (row.comment == '1.40_1.40_8.00') or (row.comment == '1.4_1.4_8'):
          if not output.has_key(row.ifo):
            output[row.ifo] = lsctables.New(lsctables.SummValueTable)
          output[row.ifo].append(row)
        # or a template bank containing a whole list of inspiral_effective_distance
        else:
          if not massOutput.has_key(row.ifo):
            massOutput[row.ifo] = [lsctables.New(lsctables.SummValueTable)]
          if len(massOutput[row.ifo]) < massNum + 1:
            massOutput[row.ifo].append(lsctables.New(lsctables.SummValueTable))
          massOutput[row.ifo][massNum].append(row)
          massNum += 1
  return output,massOutput

