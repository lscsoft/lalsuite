"""
Tools to create a HTML document that summarizes results produced by ihope. 
"""
__author__ = "Thomas Cokelaer <thomas.cokelaer@astro.cf.ac.uk>"
__version__ = "$Revision$"
__date__ = "$Date$"
__Id__ = "$Id$"
__name__ = "write_ihope_page"

import sys
import copy
import os
import re
import glob
import shutil
import StringIO
from optparse import *
from glue import markup
from glue import segments
from glue.ligolw import ligolw
from glue.ligolw import lsctables
from glue.ligolw import table
from glue.ligolw import utils
from glue.ligolw.utils import print_tables
from glue.markup import oneliner as e
import glue.pipeline

class ContentHandler(ligolw.LIGOLWContentHandler):
    pass
lsctables.use_in(ContentHandler)


# ***************************************************************************
def initialize_page(title,style,script):
  """
  A function that returns a markup.py page object with the required html
  header.
  """
  page = markup.page(mode="strict_html")
  page._escape = False
  doctype="""<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">"""
  doctype+="""\n<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">"""
  page.init(title=title, css=style, script=script , doctype=doctype)
  return page

# ***************************************************************************
def add_page_footer(page):
  """
  A function to add a valid xhtml footer to the bottom of the html.
  
  @param page: The markup.py object
  """
  # A placeholder in case something should be added to the footer of all the
  # sub pages.
  return page

# ***************************************************************************
def write_sub_page(opts,section,html_sections,style,script):
  """ 
  A function to write each of the individual sections into a markup.py object.
  
  @param section: The section being generated
  @param html_sections : A list of section titles
  @param style: The css style file sent as location of this file
  @param script: A dict of the java scripts used by the html.  
  """
  global h2_num
  global h3_num
  print "...Processing " + section
  subPage = initialize_page(section,style,script)
  subPage.p(e.br())
  subPage.add("<!-- beginning of a sub section -->")
  subPage.div(class_="contenu")
  id = section.replace(" ","_")
  h2_num += 1
  subPage.h2(str(h2_num) + ' ' + html_sections[section])
  logText(logfile, html_sections[section], "section")
  h3_num = 1
  subPage = write_results(subPage, opts, section)
  subPage.div.close()
  subPage.add("<!-- close div contenu-->")
  subPage.div.close()
  add_page_footer(subPage)
  subHtmlFile = file(opts.physdir + '/' + section + '.html',"w") 
  subHtmlFile.write(subPage(False))
  subHtmlFile.close
  mainpage.div(class_="menuitem")
  mainpage.add('\t<a class="menulink" href="javascript:loadFrame(\'' + \
      section + '.html\');"> ' + str(h2_num) + ': ' + \
      html_sections[section] + '</a>\n' )
  mainpage.div.close()
  return subPage

# ***************************************************************************
def functionId(nFramesUp):
  """ 
  Create a string naming the function n frames up on the stack.
  
  @param nFramesUp: input
  @type nFramesUp: number
  @return: message
  """
  try:
    co = sys._getframe(nFramesUp+1).f_code
    msg = "%s (%s @ %d)" % (co.co_name, co.co_filename, co.co_firstlineno)
    if msg.startswith("?") is False:
      print "-->ERROR in function: " + msg
  except:
    msg="?"

  return msg

# ***************************************************************************
def logText(logfile, text, tag="done"):
  """
  Utility to populate a logfile in HTML format. 
  The third argument is a string. Depending on its value, the text will have 
  different color. 
  
  @param text: a text to be printed
  @type text: string
  @type tag: string 
  @param tag: is in {"done","warning","error"} 
  """ 

  if tag=="warning":
    msg= "<"+tag+">"+text+"</"+tag+">\n"
    logfile.write(msg)
    if opts.verbose is True:
      print >>sys.stdout,text
  elif tag=="error":
    msg = "<"+tag+">"+text
    logfile.write(msg)
    i =1 
    while len(msg)>0:
      msg = functionId(i)
      if msg.startswith("?") is False:
        logfile.write("\n"+msg)
      else :
        logfile.write("</"+tag+">\n")
        msg=[]
      i=i+1
    print >>sys.stderr, text
  else :
    msg = "<"+tag+">"+text+"</"+tag+">\n"
    logfile.write( msg )
    if opts.verbose is True:
      print >>sys.stdout,text

# ***************************************************************************
def patternFoundInFilename(filename, pattern):
  """
   a small function to check that a pattern is contained within a filename.

     >>> filename = "H1H2L1-plotinspmissed_totalMass_versus_eta.png"
     >>> pattern = "plotinspmissed*eta"
     >>> patternfoundInFilename(filename, pattern)

   @param filename: the filename to look at
   @type filename: string
   @param pattern: a pattern which may contain several '*'
   @type pattern: string
   @return: boolean
  """
  tokeep = False
  if pattern is None : return tokeep
  
  # you may have * inside a pattern
  for word in pattern.split('*'):
    if word in filename:
      tokeep = True
    else:
      tokeep = False
      # if one word is not found within filename, we must break and se
      break

  return tokeep

# ***************************************************************************
def make_external_call(\
  command, show_stdout=False, \
  show_command=False, show_error=True):
  """
  Run a command line argument and print informative messages on failure.
  It returns two outputs: the stdout of the command, and its status.  

    >>> make_external_call('cp * /tmp', False, False, True)

  @param command: the command to try
  @type command: string
  @param show_stdout: show the stdout 
  @type show_stdout: boolean
  @param show_command: show the command
  @type show_command: boolean
  @param show_error: show the error if any
  @type show_error: boolean
  @return: the stdout and a status  


  """
  if show_command and opts.verbose is True: 
    print "--- Trying this command :" + command

  stdin, out, err = os.popen3(command)
  pid, status = os.wait()
  this_output = out.read()
  if show_error & status != 0:
    print >>sys.stderr, "External call failed."
    print >>sys.stderr, "  status: %d" % status
    print >>sys.stderr, "  stdout: %s" % this_output
    print >>sys.stderr, "  stderr: %s" % err.read()
    print >>sys.stderr, "  command: %s" % command
    sys.exit(status)
  if show_stdout:
    if this_output[0:1]=='\n': 
      print  this_output[1:]  #first character is \n
    else:
      print this_output

  stdin.close()
  out.close()
  err.close()  
  return this_output, status

# ***************************************************************************
def mkdir( newdir ):
  """
  Create a directory

  @param newdir : name of directory to be created
  @type newdir: string
  """
  if os.path.isdir(newdir): 
    print >>sys.stdout,"WARNING: this directory already exists (" + newdir +")."
    pass
  elif os.path.isfile(newdir):
    raise OSError("a file with the same name as the desired " \
                  "dir, '%s', already exists." % newdir)
  else: 
    os.mkdir(newdir)


# ***************************************************************************
def create_toggle():
  """
  This function is just an alias to create a javascript for the toggle on/off. 

  @return: nothing
  """
  fname = open("toggle.js", "w")
  fname.write("""function afterLoadFrame() {
    $('#iframecontent h3').click(function() { $(this).next().toggle('fast'); });
    $('#iframecontent a[rel="external"]').attr('target','_blank');
    $('#iframecontent input').hide();
    $('#iframecontent p:first').hide();
    }
  function loadFrame(sourceURL) {
    $("#iframecontent").load(sourceURL,{},afterLoadFrame);
    /* Remove the last two arguments to disable toggling from the title. */
    }
""")
  fname.close()

# ***************************************************************************
def write_results(page, opts, section,injection=None):
  """
  This function is just a switch to the function that create each section of
  the HTML document. 
   
    >>> write_results(page, opts, "injection")

  @param page: the HTML document page.
  @type page: markup.py document
  @param opts: the main options variable
  @type opts: write_ihope_page options
  @param section: the current section to switch to 
  @type section: string
  @return: an update of HTML document
  @param injection: The name of the injection when section=injection
  """

  if opts.verbose is True:
    print >>sys.stdout,"--------------------- Creating section "+section
  if section=='general':
    page = write_general(page, opts)
  elif section=='datainfo':
    page = write_datainfo(page, opts)
  elif section=='playsummary':
    page = write_summary(page,opts,thisSearch='playground',injection='allinj',
                         pipedown='pipedown')
  elif section=='fullsummary':
    page = write_summary(page,opts,thisSearch='full_data',injection='allinj',
                         pipedown='pipedown')
  elif section=='play_full_slide_summary':
    page = write_summary(page,opts,thisSearch='full_data_slide',
                         injection='allinj',pipedown='pipedown')
  elif section=='playground': 	
    page = write_analysis(page, opts,thisSearch='playground',
                          pipedown='pipedown')
  elif section=='full_data_slide': 	
    page = write_analysis(page, opts,thisSearch='full_data_slide',
                          pipedown='pipedown')
  elif section in hipecp.options("injections"):
    page = write_injection(page, opts,section,pipedown='pipedown')
  elif section=='allinj':
     page = write_injection(page, opts,section,pipedown='pipedown')
  elif section=='hardware_injections':
    page = write_hw_injection(page, opts)
  elif section=='full_data':
    page=write_analysis(page,opts,thisSearch='full_data',
                        pipedown='pipedown')
  elif section=='upperlimit_play': 	
    page = write_upperlimit(page, opts,thisSearch='playground',pipedown='pipedown')
  elif section=='upperlimit_full':
    page = write_upperlimit(page, opts,thisSearch='full_data',pipedown='pipedown')
  elif section=='summary_files': 
    page = write_summaryFiles(page, opts)
  elif section=='logfile':
    page = write_logfile(page, opts)
  elif section=='about': 
    page = write_about(page, opts)
  else:
    msg = "this section ("+ section +") does not exist. "
    logText(logfile, msg, "error")

  return page


# ***************************************************************************
def write_logfile(page , opts):
  """
  This function creates an entry with the logfile information.
  
  @param page: the HTML document page.
  @type page: markup.py document
  @param opts: the main options variable
  @type opts: write_ihope_page options
  @return: an update of HTML document
  """
  # get the directory of the url
  dir = opts.webdir +'/'
  
  source = dir + __name__ + '.xml'
  
  page.add("<p> Here is a logfile.")
  page = add_href(page, href=source, name=__name__)
  page.add(", which entry should be all in green except if skip options were used.</p>")
  
  return page

def add_href(page, href="", name="click here"):
  try:
    page.add("<a href=\"" + href + "\" rel=\"external\"> " + name + "</a>")
  except: 
    logText(logfile, "Could not create link "+ href, "error")
  return page

# ***************************************************************************

def create_sections(page,opts,sections,tags,captions,imagetags,images_dirs,\
                    titles,comments,configs,htmltags,htmlOnly,htmlAndPlot,\
                    allIfosFlag):
  for section in sections:
      if section in allIfosFlag:
        allIfos = True
      else:
        allIfos = False
      page = heading(page, titles[section])
      if comments[section] is not None:
        page.p(comments[section])
      page = add_config_section(page, configs[section])
      if section in htmlOnly:
        page = html_insert(page,opts, html_tag=htmltags[section],\
          caption=captions[section],directory=images_dirs[section],\
          all_ifos = allIfos)
      elif section in htmlAndPlot:
        page = html_and_plot(page, opts, cachefile_tag=tags[section],\
          caption=captions[section],image_tag=imagetags[section],\
          directory=images_dirs[section],html_tag=htmltags[section],\
          all_ifos = allIfos)
      else:
        page = fom(page, opts, cachefile_tag=tags[section],\
          caption=captions[section],image_tag=imagetags[section],\
          directory=images_dirs[section],all_ifos=allIfos)
      page.div.close() # for each heading

  return page



# ***************************************************************************
def write_general(page,opts):
  """
  Creates the general section. 

  @param page: the HTML document page.
  @type page: markup.py document
  @param opts: the main options variable
  @type opts: write_ihope_page options
  @return: an update of HTML document
  """
  webdir = opts.webdir
  ini = opts.config_file
    
  text = 'This page summarizes the analysis of the data segment from GPS '
  text += 'time %s up to %s' % (opts.gps_start_time,  opts.gps_end_time)
  page.add(text)

  # The search parameters
  if opts.verbose is True: 
    print "Extracting the version and tag of the executables..." 
  
  # todo : get the list of executables  from the ini file
  executables = []
  for item,value in hipecp.items('condor'):
    if item != 'universe':
      executables.append([item,os.path.split(value)[1]])     

  heading(page, "The search used the following resources:")

  page.ol() #starts enumerate
  # first section with version number
  # get the version of each executables
  try:
    page.ul()
    for exe in executables:
      if exe[0] == 'mvsc_train_forest' or exe[0] == 'mvsc_use_forest':
        text = "<b>"+exe[0]+ ': ' + exe[1] + "</b> "
      else:
        text = "<b>"+exe[0]+ ': ' + exe[1] + "</b> "+get_version(exe[1])
      page.li(text)   
    page.ul.close()
    page.li.close()
  except:
    logText(logfile, """Problem with the executable: cannot find them ? """, "warning")
    pass
  
  # second section 
  try:
    page.add("<li>Segment information (see the Data Information section for more details):")
    page.ul()
    for this in get_ifos():
        seg = this +'-SELECTED_SEGS'+opts.txt_suffix
        this_ref = webdir + '/segments/'+seg
        page.li(e.a(seg,href=this_ref,rel="external") )
    page.ul.close()
    page.li.close()
  except:
    logText(logfile, """Problem with the executable: cannot find them ? """, "warning")
    pass

  try:
    # let us give the ihope.ini and log file.
    # we also need to copy some files to the webdir
    cmd = 'cp '+opts.datadir + '/*.pipeline.log  '+opts.physdir
    make_external_call(cmd, opts.debug, opts.debug, True)
    cmd = 'cp '+opts.datadir + '/*.ini  '+opts.physdir
    make_external_call(cmd, opts.debug, opts.debug, True)
    pipelinelogs = glob.glob( opts.physdir + '/*.pipeline.log' )
    iniFiles = glob.glob( opts.physdir + '/*.ini' )
  except:pass

  try:
    page.add("<li>The configuration parameters are contained within the file(s) ")
    page.ul()
    for file in iniFiles:
      cmd = 'mv ' + file + ' ' + file + '.txt'
      make_external_call(cmd, opts.debug, opts.debug, True)
      logname = file.replace( opts.physdir, '' ) 
      page.li( e.a( logname, href=webdir + "/" + logname + '.txt', rel="external" ) )
    page.ul.close()
    page.add("which was/were used by lalapps_ihope as reported in the file(s):")
    page.ul()
    for file in pipelinelogs:
      cmd = 'mv ' + file + ' ' + file + '.txt'
      make_external_call(cmd, opts.debug, opts.debug, True)
      logname = file.replace( opts.physdir, '' )
      page.li( e.a( logname, href=webdir + "/" + logname + '.txt', rel="external" ) )
    page.ul.close()
    page.li.close()
  except:pass
  
  page.ol.close() # end enumerate

  logText(logfile,  "...Get the executables version...")
  page.div.close()

  # The ifo requested
  page = heading(page, "This search concerned the following combination of ifos")
  page = add_config_section(page, "ifo-details")
  page.div.close() #close the main heading

  return page


# ***************************************************************************
def write_table(page, segs, keys,printKeys = True):
  """
  print a table in html format, where keys is the first column and 
  keys to the dictionary "segs"

  @param page: the html document
  @param segs: the segment durations
  @param keys: the segments names
  @type segs: a dictionary which keys are stored in the parameter "keys"
  @type keys: a list of sorted keys corresponding to the segs dictionary keys (not sorted)
  """
  page.table()
  for key in keys:
    if segs.has_key(key):
      page.tr()
      if printKeys:
        page.td(key)
      page.td(segs[key])
      page.tr.close()
  page.table.close()

  return page
# ***************************************************************************
def write_summary(page,opts,thisSearch='playground',injection = 'allinj',
                  pipedown = 'pipedown'):
  """
  Creates the summary section.

  @param page: the HTML document page.
  @type page: markup.py document
  @param opts: the main options variable
  @type opts: write_ihope_page options
  @param thisSearch: either "playground" or "full_data" or "full_data_slide"
  @type thisSearch: string
  @return: an update of HTML document
  """
  webdir = opts.webdir
  ini = opts.config_file
  ifos = get_ifos()
  msg = "This section summarizes the analysis of the "+thisSearch+" data.<br/>"
  page.add(msg)
  if thisSearch=='playground':
    ZLdir='playground_summary_plots/'
    ifartag=['PLAYGROUND','PLAYGROUND','CLOSED_BOX']
  elif thisSearch=='full_data_slide':
    ZLdir='full_data_slide_summary_plots/'
    ifartag=['FULL_DATA','PLAYGROUND','CLOSED_BOX']
  elif thisSearch=='full_data':
    ZLdir='full_data_summary_plots/'
    ifartag=['FULL_DATA','ALL_DATA','OPEN_BOX']
  else:
    print 'Error in write_analysis function: thisSearch is not valid'
    sys.exit(1)
  IJdir = injection + '_summary_plots/'
  PDdir = pipedown + '/'
  if opts.symlink_plots:
    symlink_plot_directory(opts.datadir + ZLdir,\
        opts.physdir + ZLdir)
  else:
    copy_plot_directory(opts.datadir + ZLdir,\
        opts.physdir + ZLdir)
  if opts.symlink_plots:
    symlink_plot_directory(opts.datadir + IJdir,\
        opts.physdir + IJdir)
  else:
    copy_plot_directory(opts.datadir + IJdir,\
        opts.physdir + IJdir)
  if opts.symlink_plots:
    symlink_plot_directory(opts.datadir + PDdir,\
        opts.physdir + PDdir)
  else:
    copy_plot_directory(opts.datadir + PDdir,\
        opts.physdir + PDdir,open_box = opts.full_data)

  # A keyword to identify the section
  # title will be the name of the section.
  # tags is a tag to search for the relevant cache file (will use the first one found)
  # imagetags. if not none, will only pu a subset of images corresponding to the list provided. 
  # captions is the text to be written in the caption
  # Comments: If not none a comment will be appended to the web page
  # configs : is not none, the part of the ini file corresponding to this config name will be written in the web page
  # images_dirs : What directory to look in.

  sections = []
  htmlOnly = []
  htmlAndPlot = []
  allIfosFlag = []
  titles = {}
  tags = {}
  htmltags = {}
  imagetags = {}
  captions = {}
  comments = {}
  configs = {}
  images_dirs = {}

  section = 'inspiralrange'
  sections.append(section)
  allIfosFlag.append(section)
  titles[section] = "Inspiral range plots"
  tags[section] = '*inspiralrange*'
  imagetags[section] = ['range_plot']
  captions[section] = """Inspiral Horizon distance for a \
       (1.4,1.4) solar mass system with SNR=8 (first sub-figure)."""
  comments[section] = None
  configs[section] = None
  images_dirs[section] = ZLdir

  section = 'ifar2'
  sections.append(section)
  htmlAndPlot.append(section)
  titles[section] = "IFAR and loudest event tables with CAT2 vetoes applied (these include the CBC hardware injections)."
  tags[section] = '*plotifar*'+ifartag[0]+'*CAT_2*'+ifartag[1] + '*'
  htmltags[section] = '*'+ifartag[0]+'*CAT_2_VETO_LOUDEST*'+ifartag[1]+'*_EVENTS_BY_COMBINED_FAR*'
  imagetags[section] = ['cumhist_combined_ifar']
  captions[section] = "Combined and uncombined IFAR vs cumulative number of triggers (above). Loudest event tables as a function of combined FAR (below)."
  comments[section] = None
  configs[section] = None
  images_dirs[section] = PDdir

  section = 'ifar3'
  sections.append(section)
  htmlAndPlot.append(section)
  titles[section] = "IFAR and loudest event tables with CAT2 vetoes applied and CBC hardware injections removed."
  tags[section] = '*plotifar*'+ifartag[0]+'*CAT_3*'+ifartag[1] + '*'
  htmltags[section] = '*'+ifartag[0]+'*CAT_3_VETO_LOUDEST*'+ifartag[1]+'*_EVENTS_BY_COMBINED_FAR*'
  imagetags[section] = ['cumhist_combined_ifar']
  captions[section] = "Combined and uncombined IFAR vs cumulative number of triggers (above). Loudest event tables as a function of combined FAR (below)."
  comments[section] = None
  configs[section] = None
  images_dirs[section] = PDdir

  section = 'ifar4'
  sections.append(section)
  htmlAndPlot.append(section)
  titles[section] = "IFAR and loudest event tables with CAT2 and CAT3 vetoes applied and CBC hardware injections removed."
  tags[section] = '*plotifar*'+ifartag[0]+'*CAT_4*'+ifartag[1] + '*'
  htmltags[section] = '*'+ifartag[0]+'*CAT_4_VETO_LOUDEST*'+ifartag[1]+'*_EVENTS_BY_COMBINED_FAR*'
  imagetags[section] = ['cumhist_combined_ifar']
  captions[section] = "Combined and uncombined IFAR vs cumulative number of triggers (above). Loudest event tables as a function of combined FAR (below)."
  comments[section] = None
  configs[section] = None
  images_dirs[section] = PDdir

  section = 'inspmissed'
  sections.append(section)
  htmlAndPlot.append(section)
  titles[section] = "Found and missed plots (second stage) and the 10 closest missed injections (CAT2 and CAT3 vetoes applied and CBC hardware injections removed)."
  tags[section] = '*ligolw_cbc_plotfm_fm_dist_v_param_ALLINJ*CAT_4*'
  htmltags[section] = '*'+injection.upper()+'_CAT_4_VETO_CLOSEST_MISSED_INJECTIONS_SUMMARY*'
  imagetags[section] = ['cbc_plotfm_fm_dist_v_param']
  captions[section] = 'Found and missed injections. Stars are injections found with zero combined far; circles are injections found with non-zero combined far (their color corresponds to their far). Red crosses are missed injections; vetoed injections are excluded. The decisive distance is the second smallest effective distance. The table below shows the 10 closest missed injections (i.e., the 10 closest red crosses). For a table of the non-zero FAR injections (the colored circles), please go to the "All Injections Combined" section and click on "Loudest injections that are not louder than the loudest slide."'
  comments[section] = None
  configs[section] = None
  images_dirs[section] = PDdir
  
  page = create_sections(page,opts,sections,tags,captions,imagetags,\
                    images_dirs,titles,comments,configs,htmltags,htmlOnly,\
                    htmlAndPlot,allIfosFlag)

  return page




# ***************************************************************************
def write_datainfo(page,opts):
  """
  Creates the datainfo section.

  @param page: the HTML document page.
  @type page: markup.py document
  @param opts: the main options variable
  @type opts: write_ihope_page options
  @return: an update of HTML document
  """
  # first, get some information 
  webdir = opts.webdir
  datadir = opts.datadir
  ini = opts.config_file
  ifos = get_ifos()

  page.add("<p>This section summarizes the data which was available for analysis.")
  page.add("All of the segment files are linked in")
  page = add_href(page, href= webdir + "segments", name="this directory")
  page.add(".<p>")

  page.add("Segments were generated using the following information:")
  page = add_config_section(page, "segments")

  liveTimes = {}
  liveTimes[1] = {}
  segs = get_segments_tag("SELECTED_SEGS")
  for ifo in get_ifos():
    liveTimes[1][ifo] = float(segs[ifo][2])/86400

  # for each category, we loop over each ifo and create a table
  # with the time analysed.
  i=0
  catStrings = ["VETOTIME_CAT1"]
  categories = (hipecp.get('segments','veto-categories')).split(',')[:-1]
  for num in categories:
    catStrings.append("VETOTIME_CAT" + str(num))
  for cat in (catStrings):
    vetSegs = get_segments_tag(cat)
    i=i+1
    if i > 2:
      cumVetSegs = get_segments_tag("COMBINED_CAT_"+str(i)+"_VETO_SEGS")
    else:
      cumVetSegs = vetSegs
    msg = "Category " +str(i) + " veto segments (flags and time)"
    page = heading(page, msg) 

    readTest = False
    for ifo in get_ifos():
      if (not vetSegs.has_key(ifo)) or (not cumVetSegs.has_key(ifo)):
        readTest = True

    if readTest:
      page.add("Problem reading files at category " + cat)
      page.div.close()
      continue

    msg = "Single detector live times:"
    page.add(msg)
    key = 0
    keys = []
    data = {}
    keys.append(key)
    data[key] = ['Ifo','Live time(days)','Duty cycle (percent)']
    if i != 1:
      liveTimes[i] = {}
    for ifo in get_ifos():
      key += 1
      if i == 2:
        liveTimes[i][ifo] = liveTimes[i-1][ifo] - float(vetSegs[ifo][2])/86400
      elif i > 2:
        liveTimes[i][ifo] = liveTimes[1][ifo] - float(cumVetSegs[ifo][2])/86400
      keys.append(key)
      data[key]=[ifo,str(round(liveTimes[i][ifo],3)),str(round(liveTimes[i][ifo]*100/opts.daydur,3))]
    page = write_table(page, data,keys,printKeys = False)

    msg = "Coincident live times (only available if pipedown has finished)"
    page.add(msg)
    key = 0
    keys = []
    data = {}
    keys.append(key)
    data[key] = ['Coincidence','Coincidence time(days)','Percent of duration']
    totTime = 0
    for coinc in get_ifo_coinc():
      key += 1
      summGlob = datadir + '/pipedown/' + coinc + '-FULL_DATA'
      if i != 1:
        summGlob += '_CAT_' + str(i) + '_VETO'
      summGlob += '_LOUDEST_ALL_DATA_EVENTS_BY_COMBINED_FAR_SUMMARY-*.xml'
      summFiles = glob.glob(summGlob)
      if len(summFiles) == 1:
        xmldoc = utils.load_filename(summFiles[0], contenthandler=ContentHandler)
        loudest_event_table = table.get_table(xmldoc, 'loudest_events')
        if len (loudest_event_table):
          time = loudest_event_table[0].all_data_duration__Px_days_xP_
        else:
          time = 0
        keys.append(key)
        data[key] = [coinc,str(round(time,3)),str(round(time*100/opts.daydur,3))]  
        totTime += time
    key += 1
    if totTime != 0:
      keys.append(key)
      data[key] = ['Total',str(round(totTime,3)),str(round(totTime*100/opts.daydur,3))]
    page = write_table(page, data,keys,printKeys = False)
    for iterafsaf in [1]:
      msg = "This category includes the following flags : "
      page.add(msg)
      key = 0
      keys = []
      formFlags = {}
      keys.append(str(key))
      formFlags[str(key)] = ['IFO','Flag','Version','Start padding (s)', 'End padding (s)','Start time','End time']
      key += 1
      seg_filename = datadir +'/segments/' + hipecp.get("segments", "veto-def-file")
      xmldoc = utils.load_filename( seg_filename, gz = seg_filename.endswith(".gz"), contenthandler=ContentHandler )
      veto_def_tab = table.get_table(xmldoc,'veto_definer')
      output = StringIO.StringIO()
      tableList = ['veto_definer']
      columnList = ['ifo', 'name', 'version', 'start_pad', 'end_pad', 'category', 'start_time', 'end_time']
      print_tables.print_tables(xmldoc, output, 'wiki', tableList = tableList, columnList = columnList,
        round_floats = False, format_links = False, print_table_names = False )
      flags = output.getvalue().split('\n') 
      output.close()
      for ifo in ifos:
        for line in flags:
          if line == '' or line == '||ifo||name||version||start pad||end pad||category||start time||end time||':
            continue
          flag = line.split('||')
          flag.pop(0)
          flag.pop(-1)
          if int(flag[5]) != i:
            continue
          if flag[0] != ifo:
            continue
          flag.pop(6)
          keys.append(str(key))
          formFlags[str(key)] = flag
          key += 1

      page = write_table(page, formFlags,keys,printKeys = False)

      if i != 1:
        msg="Cumulative time vetoed at category "+str(i)+" according to txt files in segments directory. Does not include CAT 1 vetoed time."
      else:
        msg = "Time vetoed at category 1 according to *" + cat + "*txt"
      page.add(msg)

      vetKeys = ("segments","H1","H2","L1","G1","V1","T1")
      page = write_table(page, cumVetSegs, vetKeys) 
      page.div.close()

  for tag in ["SCIENCE_SEGMENTS"] :
    page = heading(page, tag +" summary ")
    if tag == "SCIENCE_SEGMENTS":
      msg = "The amount of science mode time returned by the segment database (without CAT1 flags applied) is:"

    page.add(msg +" <br/>")
    try:
      segs = get_segments_tag(tag)
      keys = ("segments","H1","H2","L1","G1","V1","T1") 
      page = write_table(page, segs,keys)
      page.div.close() 
    except:
      msg = "Problems parsing category veto segment list "+tag
      logText(logfile, msg,"warning")
      page.div.close() 
  
  return page


# ***************************************************************************
def write_upperlimit(page, opts, thisSearch='playground', pipedown='pipedown'):
  """
  Creates an upper limit section
  
  @param page: the HTML document page.
  @type page: markup.py document
  @param opts: the main options variable
  @type opts: write_ihope_page options
  @return: an update of HTML document
  """
  webdir = opts.webdir
  ini = opts.config_file
  ifos = get_ifos()

  if thisSearch=='playground':
    ifartag=['PLAYGROUND','PLAYGROUND','CLOSED_BOX']
    FARtext="""a user-specified FAR value (1/livetime by default)."""
  elif thisSearch=='full_data_slide':
    ifartag=['FULL_DATA','PLAYGROUND','CLOSED_BOX']
    FARtext="""an unknown FAR: upper limit for full data slide is not well-defined!""" #HERE BE DRAGONS
  elif thisSearch=='full_data':
    ifartag=['FULL_DATA','ALL_DATA','OPEN_BOX']
    FARtext="""the combined FAR of the loudest full data event."""
  else:
    print 'Error in write_analysis function: thisSearch is not valid'
    sys.exit(1)
  PDdir = pipedown + '/'

  page.add("This section gives the search sensitivity and upper limits using "+FARtext)
  page.add("<br/>")

  if opts.symlink_plots:
    symlink_plot_directory(opts.datadir + PDdir,\
        opts.physdir + PDdir)
  else:
    copy_plot_directory(opts.datadir + PDdir,\
        opts.physdir + PDdir,open_box = opts.full_data)

  sections = []
  htmlOnly = []
  htmlAndPlot = []
  allIfosFlag = []
  titles = {}
  tags = {}
  htmltags = {}
  imagetags = {}
  captions = {}
  comments = {}
  configs = {}
  images_dirs = {}

  section= 'upperlimitTotalMassRangeplots'
  sections.append(section)
  htmlAndPlot.append(section)
  titles[section] = "Sensitive Distance as a Function of Total Mass"
  tags[section] = '*sink_by_total_mass_'+ifartag[0]+'_CAT_4_VETO*'
  htmltags[section] = '*sink_by_total_mass_combined_upper_limit_'+ifartag[0]+'_CAT_4_VETO*'
  imagetags[section] = ['distance']
  captions[section] = """Mean sensitive distance for the CBC search pipeline as a function of total mass. The sensitive distance is determined by computing the cube root of the mean distance^3 for injections found with a FAR below """+FARtext+""" (CAT2 and CAT3 vetoes applied and CBC hardware injections removed)."""
  comments[section] = None
  configs[section] = 'search-volume'
  images_dirs[section] = PDdir

  section= 'upperlimitChirpMassRangeplots'
  sections.append(section)
  htmlAndPlot.append(section)
  titles[section] = "Sensitive Distance as a Function of Chirp Mass"
  tags[section] = '*sink_by_chirp_mass_'+ifartag[0]+'_CAT_4_VETO*'
  htmltags[section] = '*sink_by_chirp_mass_combined_upper_limit_'+ifartag[0]+'_CAT_4_VETO*'
  imagetags[section] = ['distance']
  captions[section] = """Mean sensitive distance for the CBC search pipeline as a function of chirp mass. The sensitive distance is determined by computing the cube root of the mean distance^3 for injections found with a FAR below """+FARtext+""" (CAT2 and CAT3 vetoes applied and CBC hardware injections removed)."""
  comments[section] = None
  configs[section] = 'search-volume'
  images_dirs[section] = PDdir

  section= 'upperlimitM1M2Rangeplots'
  sections.append(section)
  htmlAndPlot.append(section)
  titles[section] = "Sensitive Distance as a Function of Component Masses."
  tags[section] = '*sink_by_component_mass_'+ifartag[0]+'_CAT_4_VETO*'
  htmltags[section] = '*sink_by_component_mass_combined_upper_limit_'+ifartag[0]+'_CAT_4_VETO*'
  imagetags[section] = ['distance']
  captions[section] = """Mean sensitive distance for the CBC search pipeline as a function of the component masses. The sensitive distance is determined by computing the cube root of the mean distance^3 for injections found with a FAR below """+FARtext+""" (CAT2 and CAT3 vetoes applied and CBC hardware injections removed)."""
  comments[section] = None
  configs[section] = 'search-volume'
  images_dirs[section] = PDdir


  section= 'combupperlimitVmtotal'
  sections.append(section)
  htmlAndPlot.append(section)
  titles[section] = "Combined Upper Limits as a Function of Total Mass."
  tags[section] = '*sink_by_total_mass_'+ifartag[0]+'_CAT_4_VETO*'
  htmltags[section] = '*sink_by_total_mass_combined_upper_limit_'+ifartag[0]+'_CAT_4_VETO*'
  imagetags[section] = ['upper_limit_plot']
  captions[section] = """Combined upper limits on the CBC rate (in mergers/Mpc^3/yr) as a function of total mass. The plot shows the 90&#37; upper limits arising both from the prior and the posterior (CAT2 and CAT3 vetoes applied and CBC hardware injections removed)."""
  comments[section] = None
  configs[section] = 'upper-limit'
  images_dirs[section] = PDdir

  section= 'combupperlimitVmchirp'
  sections.append(section)
  htmlAndPlot.append(section)
  titles[section] = "Combined Upper Limits as a Function of Chirp Mass."
  tags[section] = '*sink_by_chirp_mass_'+ifartag[0]+'_CAT_4_VETO*'
  htmltags[section] = '*sink_by_chirp_mass_combined_upper_limit_'+ifartag[0]+'_CAT_4_VETO*'
  imagetags[section] = ['upper_limit_plot']
  captions[section] = """Combined upper limits on the CBC rate (in mergers/Mpc^3/yr) as a function of chirp mass. The plot shows the 90&#37; upper limits arising both from the prior and the posterior (CAT2 and CAT3 vetoes applied and CBC hardware injections removed)."""
  comments[section] = None
  configs[section] = 'upper-limit'
  images_dirs[section] = PDdir

  section= 'combupperlimitVm1m2'
  sections.append(section)
  htmlAndPlot.append(section)
  titles[section] = "Combined Upper Limits as a Function of Component Masses."
  tags[section] = '*sink_by_component_mass_'+ifartag[0]+'_CAT_4_VETO*'
  htmltags[section] = '*sink_by_component_mass_combined_upper_limit_'+ifartag[0]+'_CAT_4_VETO*'
  imagetags[section] = ['upper_limit_plot']
  captions[section] = """Combined upper limits on the CBC rate (in mergers/Mpc^3/yr) as a function of the component masses. The plot shows the 90&#37; upper limits arising both from the prior and the posterior (CAT2 and CAT3 vetoes applied and CBC hardware injections removed)."""
  comments[section] = None
  configs[section] = 'upper-limit'
  images_dirs[section] = PDdir

  page = create_sections(page,opts,sections,tags,captions,imagetags,\
                    images_dirs,titles,comments,configs,htmltags,htmlOnly,\
                    htmlAndPlot,allIfosFlag)

  return page

#***************************************************************************


def write_analysis(page, opts, thisSearch='playground',pipedown='pipedown'):
  """
  Creates the playground or full_data section. It uses the same function 
  because except the name and time analysed, the figures of merits are the same. 
  @param page: the HTML document page.
  @type page: markup.py document
  @param opts: the main options variable
  @type opts: write_ihope_page options
  @param thisSearch: either "playground" or "full_data" or "full_data_slide"
  @type thisSearch: string
  @return: an update of HTML document
  """
  webdir = opts.webdir
  ini = opts.config_file
  ifos = get_ifos()
  msg = "This section summarizes the analysis of the "+thisSearch+" data.<br/>"
  page.add(msg)
  if thisSearch=='playground':
    images_dir='playground_summary_plots/'
    ifartag=['PLAYGROUND','PLAYGROUND','CLOSED_BOX']
  elif thisSearch=='full_data_slide':
    images_dir='full_data_slide_summary_plots/'
    ifartag=['FULL_DATA','PLAYGROUND','CLOSED_BOX']
  elif thisSearch=='full_data':
    images_dir='full_data_summary_plots/'
    ifartag=['FULL_DATA','ALL_DATA','OPEN_BOX']
  else:
    print 'Error in write_analysis function: thisSearch is not valid'
    sys.exit(1)
  PDdir = pipedown + '/'
  if opts.symlink_plots:
    symlink_plot_directory(opts.datadir + images_dir,\
                      opts.physdir + images_dir)
  else:
    copy_plot_directory(opts.datadir + images_dir,\
                      opts.physdir + images_dir)
  if opts.symlink_plots:
    symlink_plot_directory(opts.datadir + PDdir,\
        opts.physdir + PDdir)
  else:
    copy_plot_directory(opts.datadir + PDdir,\
        opts.physdir + PDdir,open_box = opts.full_data)



  #table and venn diagram  
  page = heading(page, "General information", "see details here")
  try:
    segs = get_segments_tag('SELECTED_SEGS')
    keys = ("segments","H1","H2","L1","G1","V1","T1") 
    msg = "The segment files above were created with no data quality flags set"
    page.p(msg)

    page.p("The times analyzed according to hipe are:")
    coincs = get_coincident_segments(thisSearch)  

    ## here is the table. The first column is another table with the durations, 
#    page.add("<table><tr><td>\n")
    page.table()
    page.tr();
    page.td('coincidence'); 
    page.td('duration(s)'); 
    page.tr.close()
    for key in coincs.keys():        
      page.tr()
      if thisSearch=='playground':
        file = webdir +thisSearch +'/segments/'+key+"_play_segs_analyzed.txt"
      elif thisSearch=='full_data' or thisSearch=='full_data_slide':
        file = webdir +'full_data/segments/'+key+"_segs_analyzed.txt"
      page.td(e.a(key,href=file,rel="external"))
  
      page.td(coincs.get(key))
      page.tr.close()
    page.table.close()
    page.add("</td>\n")
  except:
    page.div.close()
    logText(logfile, "Error in generating analysed segments.", "error")
    pass

  try:
#    page.add("<td>")
    #create figure for the venn diagram
#    data = ""
#    for coinc in ("H1","H1H2","H2","H2L1","L1","H1L1","H1H2L1"):
#        data = data +coincs.get(coinc) + " "
#   The Venn diagram is currently non-functional, please fix!

#    success = create_venn(data, thisSearch)
#    if success is True: 
#      alt = "Venn diagramm"
#    else: 
#      alt = "No picture generated, check your matlab installation. you must \
#          have vennX.m available (from matapps/src/inspiral/matSpiral package) "
#    # and incorporate into html
#
#    comment = "Venn diagram showing distribution of"
#    for coinc in  coincs.keys():
#      comment = comment + " "+ coinc
#    comment += " coincidences "
#   Venn diagram currently unfunctional
#    page = add_figure(page, fnames =[thisSearch+"/venn_"+thisSearch+".png"],\
#      caption=comment, size="full", alt=[alt])
#    page.add("</td></tr></table>")
    
    page.add("Here below is the detail of the data and ligo-data section \
        of ihope.ini ")
    try:
      page = add_config_section(page, "data")
      page = add_config_section(page, "ligo-data")
      page.div.close()
    except:
      page.div.close()  
      pass
  except:
    page.div.close()  
    logText(logfile, "Error in generating Venn diagram.", "error")
    pass
 
  # A keyword to identify the section
  # title will be the name of the section.
  # tags is a tag to search for the relevant cache file (will use the first one found)
  # imagetags. if not none, will only pu a subset of images corresponding to the list provided. 
  # captions is the text to be written in the caption
  # Comments: If not none a comment will be appended to the web page
  # configs : is not none, the part of the ini file corresponding to this config name will be written in the web page
  # images_dirs : What directory to look in.

  sections = []
  htmlOnly = []
  htmlAndPlot = []
  allIfosFlag = []
  titles = {}
  tags = {}
  htmltags = {}
  imagetags = {}
  captions = {}
  comments = {}
  configs = {}
  images_dirs = {}

  section = 'inspiralrange'
  sections.append(section)
  allIfosFlag.append(section)
  titles[section] = "Inspiral range plots"
  tags[section] = '*inspiralrange*'
  imagetags[section] = ['range_plot','range_hist','range_mass']
  captions[section] = """ Inspiral Horizon distance for a \
       (1.4,1.4) solar mass system with SNR=8 (first plot), and \
       histograms(second plot). The last plot shows the \
       expected horizon distance for any total mass, using an SNR=8."""
  comments[section] = None
  configs[section] = None
  images_dirs[section] = images_dir

  section = 'numtemplates'
  sections.append(section)
  allIfosFlag.append(section)
  titles[section] = "Variation in template bank and triggered template size"
  tags[section] = '*plotnumtemplates*'
  imagetags[section] = None
  captions[section] = "Variation in template bank and triggered template bank size"
  comments[section] = None
  configs[section] = 'tmpltbank'
  images_dirs[section] = images_dir

  section = 'inspiral1'
  sections.append(section)
  titles[section] = "First inspiral stage (using SIRE_FIRST*_ files)"
  tags[section] = '*plotinspiral_FIRST_*'
  imagetags[section] = ['snr_vs_time','snr_cum_hist','snr_histogram']
  captions[section] = "Trigger rate at first inspiral stage"
  comments[section] = None
  configs[section] = 'inspiral'
  images_dirs[section] = images_dir

  section = 'thinca1'
  sections.append(section)
  titles[section] = "First coincidence stage (using COIRE_FIRST*_ and COIRE_SLIDE_FIRST*_ files)"
  tags[section] = '*plotthinca_FIRST_*'
  imagetags[section] = ['dummy_on_purpose_so_that_no_images_is_passed_to_the_web_page']
  captions[section] = "Trigger rate at first coincidence stage"
  comments[section] = "<span style=\"color:red\">This section has no images (on purpose)</span>. \
   To obtain the first thinca stage images, clik on the links in the caption here below"
  configs[section] = 'thinca'
  images_dirs[section] = images_dir

  section = 'inspiral2'
  sections.append(section)
  titles[section] = 'Second inspiral stage using INSPIRAL_SECOND files (CAT2 and CAT3 vetoes applied and CBC hardware injections removed).'
  tags[section] = '*plotinspiral_SECOND_*CAT_4*'
  imagetags[section] = ['snr_vs_time','snr_vs_chisq_log','snr_histogram']
  captions[section] = "Trigger rate at second inspiral stage"
  comments[section] = None
  configs[section] = 'veto-inspiral'
  images_dirs[section] = images_dir

  section = 'ifar2'
  sections.append(section)
  titles[section] = "IFAR with CAT2 vetoes applied (including CBC hardware injections)."
  tags[section] = '*plotifar*'+ifartag[0]+'*CAT_2*'+ifartag[1]+'*'
  imagetags[section] = ['cumhist_uncombined_ifar','cumhist_combined_ifar']
  captions[section] = "Combined and uncombined IFAR vs cumulative number of triggers"
  comments[section] = None
  configs[section] = 'plotifar'
  images_dirs[section] = PDdir

  section = 'ifar3'
  sections.append(section)
  titles[section] = "IFAR with CAT2 vetoes applied and CBC hardware injections removed." 
  tags[section] = '*plotifar*'+ifartag[0]+'*CAT_3*'+ifartag[1]+'*'
  imagetags[section] = ['cumhist_uncombined_ifar','cumhist_combined_ifar']
  captions[section] = "Combined and uncombined IFAR vs cumulative number of triggers"
  comments[section] = None
  configs[section] = 'plotifar'
  images_dirs[section] = PDdir

  section = 'ifar4'
  sections.append(section)
  titles[section] = "IFAR with CAT2 and CAT3 vetoes applied and CBC hardware injections removed." 
  tags[section] = '*plotifar*'+ifartag[0]+'*CAT_4*'+ifartag[1]+'*'
  imagetags[section] = ['cumhist_uncombined_ifar','cumhist_combined_ifar']
  captions[section] = "Combined and uncombined IFAR vs cumulative number of triggers"
  comments[section] = None
  configs[section] = 'plotifar'
  images_dirs[section] = PDdir

  section = 'plotcumhist2'
  sections.append(section)
  titles[section] = "Effective SNR with CAT2 vetoes applied (including CBC hardware injections)."
  tags[section] = '*plotcumhist*'+ifartag[0]+'*CAT_2*'+ifartag[2]+'*'
  imagetags[section] = ['cumhist']
  captions[section] = 'Effective SNR vs number of triggers'
  comments[section] = None
  configs[section] = None
  images_dirs[section] = PDdir

  section = 'plotcumhist3'
  sections.append(section)
  titles[section] = "Effective SNR with CAT2 vetoes applied and CBC hardware injections removed."
  tags[section] = '*plotcumhist*'+ifartag[0]+'*CAT_3*'+ifartag[2]+'*'
  imagetags[section] = ['cumhist']
  captions[section] = 'Effective SNR vs number of triggers'
  comments[section] = None
  configs[section] = None
  images_dirs[section] = PDdir

  section = 'plotcumhist4'
  sections.append(section)
  titles[section] = "Effective SNR with CAT2 and CAT3 vetoes applied and CBC hardware injections removed."
  tags[section] = '*plotcumhist*'+ifartag[0]+'*CAT_4*'+ifartag[2]+'*'
  imagetags[section] = ['cumhist']
  captions[section] = 'Effective SNR vs number of triggers'
  comments[section] = None
  configs[section] = None
  images_dirs[section] = PDdir

  section = 'plotslides'
  sections.append(section)
  titles[section] = "Trigger rates at second coincidence stage (CAT2 and CAT3 vetoes applied and CBC hardware injections removed)"
  tags[section] = '*plotslides*'+ifartag[0]+'*CAT_4*'+ifartag[2]+'*'
  imagetags[section] = ['rates','durations_per_slide']
  captions[section] = 'Trigger rate at second coincidence stage'
  comments[section] = None
  configs[section] = None
  images_dirs[section] = PDdir

  if thisSearch=='full_data':
    section = 'zerofar2'
    sections.append(section)
    titles[section] = "FAR vs SNR for loud events with CAT2 vetoes applied" 
    tags[section] = '*lalapps_cbc_plotrates*CAT_2*'
    imagetags[section] = ['cumulative_F1']
    captions[section] = "SNR vs cumulative rate of triggers"
    comments[section] = None
    configs[section] = 'zerofar'
    images_dirs[section] = PDdir

    section = 'zerofar3'
    sections.append(section)
    titles[section] = "FAR vs SNR for loud events with CAT2 vetoes applied and CBC hardware injections removed." 
    tags[section] = '*lalapps_cbc_plotrates*CAT_3*'
    imagetags[section] = ['cumulative_F1']
    captions[section] = "SNR vs cumulative rate of triggers"
    comments[section] = None
    configs[section] = 'zerofar'
    images_dirs[section] = PDdir

    section = 'zerofar4'
    sections.append(section)
    titles[section] = "FAR vs SNR for loud events with CAT2 and CAT3 vetoes applied and CBC hardware injections removed." 
    tags[section] = '*lalapps_cbc_plotrates*CAT_4*'
    imagetags[section] = ['cumulative_F1']
    captions[section] = " vs cumulative number of triggers"
    comments[section] = None
    configs[section] = 'zerofar'
    images_dirs[section] = PDdir

  section = 'printlc2zl'
  sections.append(section)
  htmlOnly.append(section)
  titles[section] = "Loudest zero lag events with CAT2 vetoes applied (including CBC hardware injections)."
  htmltags[section] = '*'+ifartag[0]+'*CAT_2_VETO_LOUDEST*'+ifartag[1]+'*_EVENTS_BY_COMBINED_FAR*'
  captions[section] = 'Loudest events after Category 2 vetoes (including hardware injections).'
  comments[section] = None
  configs[section] = None
  images_dirs[section] = PDdir

  section = 'printlc3zl'
  sections.append(section)
  htmlOnly.append(section)
  titles[section] = "Loudest zero lag events with CAT2 vetoes applied and CBC hardware injections removed."
  htmltags[section] = '*'+ifartag[0]+'*CAT_3_VETO_LOUDEST*'+ifartag[1]+'*_EVENTS_BY_COMBINED_FAR*'
  captions[section] = 'Loudest events after Category 3 vetoes (including hardware injections).'
  comments[section] = None
  configs[section] = None
  images_dirs[section] = PDdir

  section = 'printlc4zl'
  sections.append(section)
  htmlOnly.append(section)
  titles[section] = "Loudest zero lag events with CAT2 and CAT3 vetoes applied and CBC hardware injections removed."
  htmltags[section] = '*'+ifartag[0]+'*CAT_4_VETO_LOUDEST*'+ifartag[1]+'*_EVENTS_BY_COMBINED_FAR*'
  captions[section] = 'Loudest events after Category 4 vetoes (hardware injections removed).'
  comments[section] = None
  configs[section] = None
  images_dirs[section] = PDdir

  section = 'printlcsl2'
  sections.append(section)
  htmlOnly.append(section)
  titles[section] = "Loudest time slide events with CAT2 vetoes applied (including CBC hardware injections)."
  htmltags[section] = '*'+ifartag[0]+'*CAT_2_VETO_LOUDEST_SLIDE_EVENTS_BY_COMBINED_FAR*'
  captions[section] = 'Loudest time slide events after Category 2 vetoes.'
  comments[section] = None
  configs[section] = None
  images_dirs[section] = PDdir

  section = 'printlcsl3'
  sections.append(section)
  htmlOnly.append(section)
  titles[section] = "Loudest time slide events with CAT2 vetoes applied and CBC hardware injections removed."
  htmltags[section] = '*'+ifartag[0]+'*CAT_3_VETO_LOUDEST_SLIDE_EVENTS_BY_COMBINED_FAR*'
  captions[section] = 'Loudest time slide events after Category 3 vetoes.'
  comments[section] = None
  configs[section] = None
  images_dirs[section] = PDdir


  section = 'printlcsl4'
  sections.append(section)
  htmlOnly.append(section)
  titles[section] = "Loudest time slide events with CAT2 and CAT3 vetoes applied and CBC hardware injections removed."
  htmltags[section] = '*'+ifartag[0]+'*CAT_4_VETO_LOUDEST_SLIDE_EVENTS_BY_COMBINED_FAR*'
  captions[section] = 'Loudest time slide events after Category 4 vetoes.'
  comments[section] = None
  configs[section] = None
  images_dirs[section] = PDdir

  page = create_sections(page,opts,sections,tags,captions,imagetags,\
                    images_dirs,titles,comments,configs,htmltags,htmlOnly,\
                    htmlAndPlot,allIfosFlag)

  return page


# ***************************************************************************
def write_summaryFiles(page, opts):
  """
  Creates a section to provide the relevant summary files such as COIRE_SUMMARY
  files

  @param page: the HTML document page.
  @type page: markup.py document
  @param opts: the main options variable
  @type opts: write_ihope_page options
  @return: an update of HTML document
  """
  webdir = opts.webdir
  ini = opts.config_file
  directories = []
  if opts.playground:
    directories.append("playground")
  if opts.full_data or opts.full_data_slide:  
    directories.append("full_data")
  if opts.injection:
    injections = hipecp.items("injections")
    injections.sort()
    for this in injections: 
      directories.append(this[0])

  page.add("This section provides a list of COIRE and SIRE summary files.<br/>")

  for thisSearch in directories:
    page.add('<h1>' + thisSearch  +' section</h1>')
    page.add("Summary files are available in the search directory")
    page = add_href(page, href= webdir +thisSearch, name="here")

    path = opts.physdir + "/" + thisSearch
    mkdir(path)
    # link all the SUMMARY SIRE text files
    page = heading(page, "Final SIRE Summary files ")
    prefix = opts.datadir + thisSearch + '/'
    ifos = get_ifos()
    for ifo in ifos:
      if 'inj' in thisSearch:
        fname = ifo + '-SIRE_INJECTIONS*SUMMARY_FIRST_' + thisSearch.upper() \
            + opts.txt_suffix     
      else:
        fname = ifo + '-SIRE_SUMMARY_FIRST_' + thisSearch.upper() \
            + opts.txt_suffix     
      page.add( "<p>SIRE Summary file for "+ifo+".</p>")
      flist = glob.glob(prefix + fname)
      if len(flist) == 1: 
        file = open(flist[0], "r")
        text = file.read()
        file.close()
        page.add("<pre>"+fname+"<br/>"+text+"</pre>")
      else:
        page.add('<pre>Could not find '+fname+'</pre>')

    # copy all SUMMARY SIRE file 
    try:
      dest = path + "/sire"
      mkdir(dest)
      src = opts.datadir+thisSearch+'/'+"*SIRE*SUMMARY*xml* " 
      if opts.no_copy is False:
        command = 'cp ' +src +' ' + dest
        make_external_call(command, opts.debug, opts.debug, True)   
      page.add("All SIRE files are available in sire directory")
      page = add_href(page, href= webdir +thisSearch+"/sire", name="here")
    except: 
      page.add("<pre>Problem while copying the sire summary files.</pre>")
    page.div.close()

    # link and copy all the SUMMARY COIRE text files
    
    if 'inj' in thisSearch:
      page = heading(page, "Final COIRE Summary files ")
      for ifostring in get_ifo_coinc():
        if len(ifostring)>2:
          fnames = []
          inj_fname = ifostring +'-COIRE_INJECTIONS*SUMMARY_SECOND_' + \
              ifostring +'_'+thisSearch.upper()+'_CAT_3_VETO' +opts.txt_suffix
          fnames.append(inj_fname)
          for fname in fnames:
            page.add( "<p>COIRE Summary file "+ ifostring +".</p>")
            flist = glob.glob(prefix + fname)
            if len(flist) == 1: 
              file = open(flist[0], "r")
              text = file.read()
              file.close()
              page.add("<pre>" + fname + "<br/>" + text + "</pre>")
            else: 
              page.add('<pre>Could not find '+fname+'</pre>')
      try:
        dest = path + "/coire"
        mkdir(dest)
        src = opts.datadir+thisSearch+'/'+"*COIRE*SUMMARY*xml* "
        if opts.no_copy is False:
          command = 'cp ' +src +' ' + dest
          make_external_call(command, opts.debug, opts.debug, True)
        page.add("All COIRE files are available in coire directory")
        page = add_href(page, href= webdir +thisSearch+"/coire", name="here")
      except:
        page.add("<pre>Problem while copying the coire summary files.</pre>")
      page.div.close()

    if thisSearch == "full_data" or thisSearch == "playground":
      if thisSearch == "playground" or opts.full_data:
        page = heading(page, "Final COIRE zero lag files")
        for ifostring in get_ifo_coinc():
          if len(ifostring)>2:
            fnames = []
            fname = ifostring +'-COIRE_SUMMARY_SECOND_' + ifostring \
                + '_' + thisSearch.upper() + '_CAT_3_VETO' + opts.txt_suffix
            fnames.append(fname)
            for fname in fnames:
              page.add( "<p>COIRE Zero Lag Summary file "+ ifostring +".</p>")
              flist = glob.glob(prefix + fname)
              if len(flist) == 1:
                file = open(flist[0], "r")
                text = file.read()
                file.close()
                page.add("<pre>" + fname + "<br/>" + text + "</pre>")
              else:
                page.add('<pre>Could not find '+fname+'</pre>')
        try:
          dest = path + "/coire_zero_lag"
          mkdir(dest)
          src = opts.datadir+thisSearch+'/'+"*COIRE_SUMMARY*xml* "
          if opts.no_copy is False:
            command = 'cp ' +src +' ' + dest
            make_external_call(command, opts.debug, opts.debug, True)
          page.add("All COIRE zero lag files are available in coire directory")
          page = add_href(page, href= webdir +thisSearch+"/coire_zero_lag",\
	         name="here")
        except:
          page.add("<pre>Problem while copying the coire summary files.</pre>")
        page.div.close()

      if thisSearch == "playground" or opts.full_data_slide:
        page = heading(page, "Final COIRE slide files")
        for ifostring in get_ifo_coinc():
          if len(ifostring)>2:
            fnames = []
            slidefname = ifostring +'-COIRE_SLIDE_SUMMARY_SECOND_' + ifostring \
                + '_' + thisSearch.upper() + '_CAT_3_VETO' + opts.txt_suffix
            fnames.append(slidefname)
            for fname in fnames:
              page.add( "<p>COIRE Slide Summary file "+ ifostring +".</p>")
              flist = glob.glob(prefix + fname)
              if len(flist) == 1:
                file = open(flist[0], "r")
                text = file.read()
                file.close()
                page.add("<pre>" + fname + "<br/>" + text + "</pre>")
              else:
                page.add('<pre>Could not find '+fname+'</pre>')
        try:
          dest = path + "/coire_slide"
          mkdir(dest)
          src = opts.datadir+thisSearch+'/'+"*COIRE_SLIDE_SUMMARY*xml* "
          if opts.no_copy is False:
            command = 'cp ' +src +' ' + dest
            make_external_call(command, opts.debug, opts.debug, True)
          page.add("All COIRE slide files are available in coire directory")
          page = add_href(page, href= webdir +thisSearch+"/coire_slide",\
	         name="here")
        except:
          page.add("<pre>Problem while copying the coire summary files.</pre>")
	page.div.close()

  return page
 

# ***************************************************************************
def write_injection(page, opts,injection,pipedown='pipedown'):
  """
  Creates the injection section

  @param page: the HTML document page.
  @type page: markup.py document
  @param opts: the main options variable
  @type opts: write_ihope_page options
  @return: an update of HTML document
  """
  webdir = opts.webdir
  ini = opts.config_file
  ifos = get_ifos()
  images_dir = injection + '_summary_plots/'
  if opts.symlink_plots:
    symlink_plot_directory(opts.datadir + images_dir,\
                      opts.physdir + images_dir)
  else:
    copy_plot_directory(opts.datadir + images_dir,\
                      opts.physdir + images_dir)

  PDdir = pipedown + '/'
  if opts.symlink_plots:
    symlink_plot_directory(opts.datadir + PDdir,\
        opts.physdir + PDdir)
  else:
    copy_plot_directory(opts.datadir + PDdir,\
        opts.physdir + PDdir,open_box = opts.full_data)

  page.add("This section summarizes the analysis of the injection runs.<br/>")


  
  #---------------------------- the injection section
  # title will be the name of the section.
  # tags is a tag to search for the relevant cache file (will use the first one found)
  # imagetags. if not none, will only pu a subset of images corresponding to the list provided. 
  # captions is the text to be written in the caption
  # configs : is not none, the part of the ini file corresponding to this config name will be written in the web page
 
  sections = []
  htmlOnly = []
  htmlAndPlot = []
  allIfosFlag = []
  titles = {}
  tags = {}
  htmltags = {}
  imagetags = {}
  captions = {}
  comments = {}
  configs = {}
  images_dirs = {}

  section = 'foundmissed1'
  sections.append(section)
  titles[section] = "Efficiency, Found and Missed plots (first stage)"
  tags[section] = '*plotinspmissed_FIRST*'
  imagetags[section] = ['dummy_on_purpose_so_that_no_images_is_passed_to_the_web_page']
  captions[section] = "Found and Missed injections links (first stage)"
  comments[section] = "<span style=\"color:red\">This section has no images \
        (on purpose)</span>. To obtain the efficiency and found/missed plots \
        at the first stage, click on the links in the caption here below"
  configs[section] = 'plotinspmissed'
  images_dirs[section] = images_dir

  section = 'foundmissed2'
  sections.append(section)
  titles[section] = "Closeby missed injections (Second coincidence stage CAT2 and CAT3 vetoes applied and CBC hardware injections removed)."
  tags[section] = '*plotinspmissed_SECOND*CAT_4*'
  imagetags[section] = ['map-8','map-9']
  captions[section] = 'Found and Missed injections. Effective distance versus total mass'
  comments[section] = None
  configs[section] = 'plotinspmissed'
  images_dirs[section] = images_dir

  section = 'efficiency2'
  sections.append(section)
  titles[section] = "Efficiencies (Second stage, CAT2 and CAT3 vetoes applied and CBC hardware injections removed)"
  tags[section] = '*plotinspmissed_SECOND*CAT_4*'
  imagetags[section] = ['efficiency-eff_dist']
  captions[section] = 'Efficiencies (Second stage)'
  comments[section] = None
  configs[section] = 'plotinspmissed'
  images_dirs[section] = images_dir

  section = 'inspinj1'
  sections.append(section)
  titles[section] = "Accuracies (First stage)"
  tags[section] = '*plotinspinj_FIRST*'
  imagetags[section] = ['mchirp_vs_snr_accuracy_scatter_mtotal']
  captions[section] = 'Accuracy of recovered chirp mass vs SNR (first stage)'
  comments[section] = '<span style=\"color:red\">The colours in these plots \
   represent the chirp mass of the injection. </span>'
  configs[section] = 'plotinspinj'
  images_dirs[section] = images_dir

  section = 'inspinj2'
  sections.append(section)
  titles[section] = "Accuracies (Second stage, CAT2 and CAT3 vetoes applied and CBC hardware injections removed)"
  tags[section] = '*plotinspinj_SECOND*CAT_4*'
  imagetags[section] = ['mchirp_vs_snr_accuracy_scatter_mtotal']
  captions[section] = 'Accuracy of recovered chirp mass vs SNR (second stage)'
  comments[section] = '<span style=\"color:red\">The colours in these plots \
   represent the chirp mass of the injection. </span>'
  configs[section] = 'plotinspinj'
  images_dirs[section] = images_dir

  section = 'ethinca'
  sections.append(section)
  titles[section] = "The Ethinca Parameter (CAT2 and CAT3 vetoes applied and CBC hardware injections removed)"
  tags[section] = '*plotethinca_SECOND*CAT_4*'
  imagetags[section] = ['ethinca_vs_combined']
  captions[section] = "Ethinca versus combined SNR in the different ifo combinations"
  comments[section] = None
  configs[section] = 'plotethinca'
  images_dirs[section] = images_dir

  section = 'chisq'
  sections.append(section)
  titles[section] = "The chisq parameter (CAT2 and CAT3 vetoes applied and CBC hardware injections removed)" 
  tags[section] = '*plotsnrchi*CAT_4*'
  imagetags[section] = ['chisq_inj_newsnr_lines','rsq']
  captions[section] = "chi-square versus SNR and rsq versus SNR"
  comments[section] = "Chi-squared vs SNR with contours of constant newSNR, \
    the second plot shows the rsq veto."
  configs[section] = 'plotsnrchi'
  images_dirs[section] = images_dir

  section = 'injsqf'
  sections.append(section)
  htmlOnly.append(section)
  titles[section] = "Found injections that are not louder than the loudest slide (CAT2 and CAT3 vetoes applied and CBC hardware injections removed)"
  htmltags[section] = '*'+injection.upper()+'_CAT_4_VETO_QUIETEST_FOUND_COMPARED_TO_PLAYGROUND*'
  captions[section] = 'This table lists found injections that are not louder than the loudest event in the (full data) time slides. Equivalently, found injections that have a FAR greater than zero. These tables are shown after category 4 vetoes are applied.'
  comments[section] = None
  configs[section] = None
  images_dirs[section] = PDdir

  section = 'injscm'
  sections.append(section)
  htmlOnly.append(section)
  titles[section] = "Nearby missed injections (CAT2 and CAT3 vetoes applied and CBC hardware injections removed)"
  htmltags[section] = '*'+injection.upper()+'_CAT_4_VETO_CLOSEST_MISSED_INJECTIONS_SUMMARY*'
  captions[section] = 'This table shows a list of the closest injections that were not marked as "found" by our pipeline. These tables are shown after category 4 vetoes are applied.'
  comments[section] = None
  configs[section] = None
  images_dirs[section] = PDdir

  section = 'pipedownacc1'
  sections.append(section)
  titles[section] = "Accuracies I (Pipedown, CAT2 and CAT3 vetoes applied and CBC hardware injections removed)"
  tags[section] = '*ligolw_cbc_plotfm_fm_dist_v_dt*'+injection.upper()+'_CAT_4_VETO*'
  imagetags[section] = ['ligolw_cbc_plotfm_fm_dist_v_dt']
  captions[section] = 'Injected decisive distance vs difference between recovered and injected end times of "found" injections. Stars are injections found with zero combined far; circles are injections found with non-zero combined far. The decisive distance is the second smallest effective distance.'
  comments[section] = None
  configs[section] = None
  images_dirs[section] = PDdir

  section = 'pipedownacc2'
  sections.append(section)
  titles[section] = "Accuracies II (Pipedown, CAT2 and CAT3 vetoes applied and CBC hardware injections removed)"
  tags[section] = '*ligolw_cbc_plotfm_fm_lin*'+injection.upper()+'_CAT_4_VETO*'
  imagetags[section] = ['ligolw_cbc_plotfm_fm_lin']
  captions[section] = 'Recovered parameter accuracies vs. difference between recovered and injected end times and injection tags. Stars are triggers with zero combined far, circles are triggers with non-zero far.'
  comments[section] = None
  configs[section] = None
  images_dirs[section] = PDdir

  section = 'pipedownfm1'
  sections.append(section)
  titles[section] = "Found and missed plots I (Pipedown, CAT2 and CAT3 vetoes applied and CBC hardware injections removed)"
  tags[section] = '*ligolw_cbc_plotfm_fm_dist_v_end_time*'+injection.upper()+'_CAT_4_VETO*'
  imagetags[section] = ['ligolw_cbc_plotfm_fm_dist_v_end_time']
  captions[section] = 'Found and missed injections: injected decisive distance vs. injected end time by mass bin. Stars are injections found with zero combined far; circles are injections found with non-zero combined far. Red crosses are missed injections (vetoed injections are excluded). The decisive distance is the second smallest effective distance.'
  comments[section] = None
  configs[section] = None
  images_dirs[section] = PDdir

  section = 'pipedownfm2'
  sections.append(section)
  titles[section] = "Found and missed plots II (Pipedown, CAT2 and CAT3 vetoes applied and CBC hardware injections removed)"
  tags[section] = '*ligolw_cbc_plotfm_fm_dist_v_param*'+injection.upper()+'_CAT_4_VETO*'
  imagetags[section] = ['ligolw_cbc_plotfm_fm_dist_v_param']
  captions[section] = 'Found and missed injections: injected decisive distance vs. injected mass. Stars are injections found with zero combined far; circles are injections found with non-zero combined far. Red crosses are missed injections (vetoed injections excluded). The decisive distance is the second smallest effective distance.'
  comments[section] = None
  configs[section] = None
  images_dirs[section] = PDdir

  page = create_sections(page,opts,sections,tags,captions,imagetags,\
                    images_dirs,titles,comments,configs,htmltags,htmlOnly,\
                    htmlAndPlot,allIfosFlag)

  return page

# ***************************************************************************

def write_hw_injection(page, opts):
  """
  Creates the injection section

  @param page: the HTML document page.
  @type page: markup.py document
  @param opts: the main options variable
  @type opts: write_ihope_page options
  @return: an update of HTML document
  """
  webdir = opts.webdir
  copy_plot_directory(opts.datadir + "hardware_injection_summary",\
      opts.physdir + "hardware_inj")  

  page.add("This section summarizes the analysis of the hardware injections.<br/>")
  
  sections = []
  htmlOnly = []
  htmlAndPlot = []
  allIfosFlag = []
  titles = {}
  tags = {}
  htmltags = {}
  imagetags = {}
  captions = {}
  comments = {}
  configs = {}
  images_dirs = {}

  section = 'hwinjpage1'
  sections.append(section)
  htmlAndPlot.append(section)
  titles[section] = "Hardware injection summary page at CAT 1"
  tags[section] = '*HWINJ_SUMMARY-*'
  htmltags[section] = '*HWINJ_SUMMARY-*'
  imagetags[section] = None
  captions[section] = ''
  comments[section] = None
  configs[section] = None
  images_dirs[section] = "hardware_inj"

  section = 'hwinjpage2'
  sections.append(section)
  htmlAndPlot.append(section)
  titles[section] = "Hardware injection summary page at CAT 2"
  tags[section] = '*HWINJ_SUMMARY_CAT_2*'
  htmltags[section] = '*HWINJ_SUMMARY_CAT_2*'
  imagetags[section] = None
  captions[section] = ''
  comments[section] = None
  configs[section] = None
  images_dirs[section] = "hardware_inj"

  page = create_sections(page,opts,sections,tags,captions,imagetags,\
                    images_dirs,titles,comments,configs,htmltags,htmlOnly,\
                    htmlAndPlot,allIfosFlag)

  return page

# ***************************************************************************
def write_about(page, opts):
  """
  Creates the section "About". 

  @param page: the HTML document page.
  @type page: markup.py document
  @param opts: the main options variable
  @type opts: write_ihope_page options
  @return: an update of HTML document
  """
  webdir = opts.webdir
  
  page.p("This page was automatically generated with write_ihope_page using\
      the following ini file")
  tmp  = open(opts.config)
  tmp2 = tmp.read()
  tmp.close() 
  # the < character in a <pre> HTML tag is read and interpreted, So, we need to 
  # change it to &lt
  page.pre(tmp2.replace('<', '&lt;'))
  # page.pre(tmp.read())

  page.p("and the following command line arguments:")
  text=""
  for arg in sys.argv[:]:
    text = text +  arg +" "
  page.pre( text )
  
  page.p(__Id__[4:len(__Id__)-5])
  return page


# ***************************************************************************
def add_config_section(page, section):
  """
  Copy and paste a section of the ihope.ini into the HTML page within verbatim
  tags

  @param page: the HTML document page.
  @type page: markup.py document
  @param section: the name of the section to copy and paste
  @type section: string
  @return: an update of HTML document
  """
  # section may be empty or set to None
  try:
    ini  = hipecp.items(section)
  except:
    return page

  page.add("<pre>")
  page.add("------------------------------------" +section)
  for i in  ini :
    page.add(i[0]+' = '+i[1])
  page.add("</pre>")
  return page


# ***************************************************************************
def heading(page, title="None", label="Switch details on/off", heading="h3"):
  """
  Add a hx HTML section within the document with the ability to toggle it on/off
   
  @param page: the HTML document page.
  @type page: markup.py document
  @param title: the name of the section
  @type title: string
  @param label: the label to put on the toggle button 
  @type label: string
  @param heading: the HTML heading (h3 by default)
  @type heading: string
  @return: an update of HTML document
  """
  #increment block number
  global count_block
  global h2_num,h3_num
  global fig_num

  input=str(count_block)
  count_block=count_block+1

  fig_num = 1

  if heading == 'h2':
    section_num = str(h2_num)+"."
  elif heading == 'h3':
    section_num = str(h2_num)+"."+str(h3_num)+"."
  else:
    raise ValueError, "heading must be either h2 or h3 in heading function"


  page.add("<"+heading+">"+ section_num+ title)
  h3_num = h3_num + 1
  text = label 

  page.add("</"+heading+">")
  page.div(id="div_"+input , style='display:none') 

  logText(logfile, '      Enters sub-section: '+title )

  return page 

# ***************************************************************************
def add_caption(page, caption):
  """
  Add a caption to the HTML document. Should be used with the figure only. 
  Take care of the figure numbering.

  @param page: the HTML document page.
  @type page: markup.py document
  @param caption: the name of the section to copy and paste
  @type caption: string
  @return: an update of HTML document
  """
  global fig_num
  page.p("<b>Figure "+str(fig_num) + ":</br> "+caption)
  fig_num = fig_num + 1

  return page


# ***************************************************************************
def add_figure(page,webDir,fnames="test", caption="add a caption", size=None, \
    alt="no figure found",source="not specified",html=None,html_file = None, \
    all_ifos = False):
  """
  Add a figure to the HTML document, taking care of the numbering.
 
  @param page: the HTML document page.
  @type page: markup.py document
  @param caption: the name of the section to copy and paste
  @type caption: string
  @param size: the size of the image ("full", "half", "third")
  @type size: string
  @param alt: an HTML alt 
  @type alt: a string
  @param source: the source of the figure
  @type source: string
  @return: an update of HTML document
  """
  global fig_num
  page.add("<!-- insert a figure -->\n<div class=\"figure\">")
  this_count = 0
 
  # set the size of the images
  if size==None:
    if len(fnames)==1: 
      size="full"
    elif len(fnames)==2:
      size="half"
    else : size="third"

  # Determine the ifo
  ifo = None
  if html:
    for combo in ['H1','H2','L1','V1','H1L1','H1H2','H2L1','H2V1','H1V1','L1V1','H1L1V1','H1H2L1','H2L1V1','H1H2V1','H1H2L1V1']:
      if html.startswith(combo):
        ifo = combo
  if ifo and not all_ifos:
    page.h5(ifo + ' were operating')

# Source is the thumnail, target is the full-resolution image. 

  for fnam in fnames:
    source=webDir+"/"+fnam.strip()
    title  = ' title="%s" ' % alt[this_count].strip()
    title  = title.replace("_thumb","")
    target = source.replace("_thumb","")
    if size == 'full':
      source = source.replace("_thumb","")
    page.add("\t<a href=\"" + target+"\"" +  title + " rel=\"external\">\n" )
    try:
      thisalt = alt[this_count]
      thisalt.replace("\n","")
      page.img(class_=size ,src=source, alt=thisalt )
    except:
      page.img(class_=size ,src=source )
      
    page.add("</a>")
    this_count = this_count +1

  section_num = str(h2_num)+"."+str(h3_num-1)+"."+str(fig_num)
  page.add("\t<p><b>Figure " + section_num + "</b>:  " + caption + "</p>")
  fig_num +=1

  # add a yellow box (depends on yourcss style though) that contains the 
  # filename of the pictures.
  if fnames is not None:
    page.add("<pre>Source: "+webDir+"\n" )
    for fnam in fnames:
      page.add(fnam.replace("\n",""))
    page.add("</pre>")

  if html_file:
    this = open(html_file, "r")
    for line in this:
      page.add(line)
    page.add("<pre>Source: "+ html_file + "\n </pre>" )
    this.close()

  page.div.close()

  #check that the file exsits. If not, returns an error message
  return page

# ***************************************************************************
def create_venn(data, tag):
  """
  Create a venn diagram for the 2 or 3 ifos case
  data has to be H1H2 H1 H2L1 H2 H1L1 L1 H1H2L1 array 

    >>> data = 
    >>> tag = 
    >>> create_venn(data, tag)

  @param data: an input dictionary
  @param data: numbers
  @param tag: 
  @type tag: string
  @return: a boolean. True is succesfule, False otherwise.
  """

  try:
    mscript = open("temp.m", "w")
    mscript.write("data = [")
    mscript.write( data +"];" ) 
    mscript.write(" vennX(data\'/3600/24,0.01);")
    mscript.write(" k=colormap(jet); k = [1 1 1; k];colormap(k);")
    mscript.write("saveas(gcf,\'venn_"+tag+".png\')")
    mscript.close()
    command=("matlab -nodisplay -nodesktop -nojvm -nosplash   < temp.m > /tmp/test ;  rm -f temp.m; mv venn_"+tag+".png "+opts.physdir+"/"+tag+"/")
    if not opts.debug:
      make_external_call(command, opts.debug, opts.debug, True)  
    return True 
  except:
    logText(logfile, """WARNING   The matlab command to create the venn diagram failed. 
                Check that matlab is properly set,and vennX.m is available
                (see matapps/src/searches/inspiral/matspiral/utilities")  
          """, "error")
    return False

# ***************************************************************************
def get_coincident_segments(tag):
  """
  @param tag: "playground" or "full_data" 
  @type tag: string
  return: a dictionary with the coincident time for each ifo combination
  """
  ifos = get_ifos()
  thisdata = {}
  thisdata['segments'] = ['duration(days)' ,'duration (s)']
  output={}
  ifo_coincs = get_ifo_coinc()  

  try:
    for coinc in ifo_coincs:
      if tag=="playground":
        command = "awk \'{sum=sum+$4} END {print sum}\' "\
          +opts.datadir+ tag +"/"+coinc+"_play_segs_analyzed.txt"
      elif tag=="full_data" or tag=="full_data_slide":
        command = "awk \'{sum=sum+$4} END {print sum}\' "\
          +opts.datadir+"full_data/"+coinc+"_segs_analyzed.txt"

      output[coinc], status = \
        make_external_call(command, False,opts.debug, True)

    logText(logfile, '...Get the analyzed segments duration...')
  except:
    logText(logfile , 'problem(s) while parsnig the coincident segments', \
      "error")

  return output

  


# ***************************************************************************
def get_segments_tag(tag):
  """
  reads segment files and return table of duration for each ifo
 
    >>> get_segments_tag("SELECTED_SEGS")

  @param tag: a tag to the segments
  @type tag: string
  return: 
  """
  # variables for file location
  datadir = opts.datadir
  ifos = get_ifos()
  thisdata = {}
  directory = datadir+'/segments/'
  thisdata['segments'] = ['filename', 'duration(days)' ,'duration (s)']

  # the name of the files (except for the ifo)
  this_tag  = '-' + tag + '-' + opts.gps_start_time + '-' + opts.duration + '.xml'

  for ifo in ifos:
    # read segments from science segments xml file
    scifile          = datadir + '/segments/' + ifo + '-SCIENCE_SEGMENTS-' + opts.gps_start_time + '-' + opts.duration + '.xml'
    scifile_xmldoc, digest  = utils.load_fileobj(open(scifile, 'r'), contenthandler=ContentHandler)
    scifile_segtable = table.get_table(scifile_xmldoc, lsctables.SegmentTable.tableName)
    scifile_seglist  = segments.segmentlist()
    for seg in scifile_segtable:
      endtime   = seg.end_time + seg.end_time_ns * 10**(-9)
      starttime = seg.start_time + seg.start_time_ns * 10**(-9)
      scifile_seglist.append(segments.segment(starttime, endtime))

    # read segments from tag xml file
    file                 = datadir + '/segments/' + ifo + this_tag
    file_xmldoc, digest  = utils.load_fileobj(open(file, 'r'), contenthandler=ContentHandler)
    file_segtable        = table.get_table(file_xmldoc, lsctables.SegmentTable.tableName)
    file_seglist         = segments.segmentlist()
    for seg in file_segtable:
      endtime   = seg.end_time + seg.end_time_ns * 10**(-9)
      starttime = seg.start_time + seg.start_time_ns * 10**(-9)
      file_seglist.append(segments.segment(starttime, endtime))

    # find the intersection of science segments and tag segments
    scifile_seglist &= file_seglist
    ifotime = abs(scifile_seglist)
    thisdata[ifo] = [os.path.basename(file),  str(float(ifotime)/86400.0), ifotime]

  return thisdata




# ***************************************************************************
def get_ifo_coinc():
  """
  return: list of valid coincidences with respect to the ini file

  """
  # get the ifo requested
  ifos = get_ifos()
  # get number of ifos to look at
  numifos = []  
  for option in ["one-ifo","two-ifo","three-ifo","four-ifo", "five-ifo"]:
    if hipecp.has_option("ifo-details",option): 
	tmp = option.split("-")
 	numifos.append(tmp[0])
  # now let us fill the different possible ifo combinations
  output = []
  for num in numifos:
    if num=="one":
      for ifo in ifos:
        output.append(ifo)
    elif num=="two":
      for ifo1 in ifos:
        for ifo2 in ifos:
          if ifo1 < ifo2:
            output.append(ifo1+ifo2)
    elif num=="three":
      for ifo1 in ifos:
        for ifo2 in ifos:
          for ifo3 in ifos:
            if ifo1 < ifo2 and ifo2 < ifo3:
              output.append(ifo1+ifo2+ifo3)
    elif num=="four":
      for ifo1 in ifos:
        for ifo2 in ifos:
          for ifo3 in ifos:
            for ifo4 in ifos:
              if ifo1 < ifo2 and ifo2 < ifo3 and ifo3 < ifo4:
                output.append(ifo1+ifo2+ifo3+ifo4)
  return output


  


# ***************************************************************************
def get_ifos():
  """
  read the ifos used in the ifo-details section of the ihope ini file
  return: list of ifos
  """
  ifos=[] 
  for option in ["g1-data","h1-data","h2-data","l1-data", "v1-data"]:
    if hipecp.has_option("ifo-details",option): ifos.append(option[0:2].upper() )

  return ifos


# ***************************************************************************
# ***************************************************************************
def get_version(executable): 
  """
  Search for the tag and version of an executable using the --version argument.

  @param executable: the name of an executable
  @type executable: string
  return: the tag name if any and the version of the executable
  """
  output=[]
  try:
    # this is surely useless now to make a distinction with trigbank 
    cmd = opts.ihope_directory + "/" + opts.gps_start_time + "-" + \
    opts.gps_end_time + "/executables/" + executable + " --version"
    output,status = make_external_call(cmd, opts.debug, opts.debug, True)
  except:  
    output= '(not found)' 

  output = output.split('\n')
  
  #searching for the version, which may be empty(undefined)
  version = '. Undefined version (does this executable have --version?).'
  for entry in output:
    if entry.startswith('$lalappsCommitID'):
      version = 'Version is '+ entry.split(' ')[1] +'.'

  return version

# ***************************************************************************
def copy_segments():
  """
  This function copies the segments list into the web page directory 
  """

  msg =   "Copying segments into the web page directory (in yoururl/segments)"
  logText(logfile, msg)
  # first we create this directory      
  mkdir(opts.physdir+'/segments')
  mkdir(opts.physdir+'/catlists')
  # parsing the ini file, find the cat file and thenread the ./segments directory
  try:
    location = opts.physdir + '/segments/'
    if opts.no_copy is False:
      command = 'cp '+opts.datadir +'/segments/* ' + location
      dummy,status = make_external_call(command, opts.debug, opts.debug, True)
      command = 'cp '+opts.datadir +'/segments/*cat*.txt ' + opts.physdir + '/catlists/' 
      dummy,status = make_external_call(command, opts.debug, opts.debug, True)
  except:
    logText(logfile, "Could not copy the segment files", "error")
    pass

 #  the selected segment files
  for thisSearch in ['playground', 'full_data']:
    mkdir(opts.physdir+thisSearch)
    try :
      stdout ='Copying the selected segments files into /'+thisSearch
      for this in get_ifo_coinc():
          if thisSearch=='playground':
            seg = '/'+thisSearch+'/'+this +'_play_segs_analyzed.txt'
          elif thisSearch=='full_data':
            seg = '/'+thisSearch+'/'+this +'_segs_analyzed.txt'

          if opts.no_copy is False:
            dest = opts.physdir + thisSearch + '/segments'
            mkdir(dest)
            command = 'cp '+opts.datadir + seg + ' ' + dest
            dummy,status = \
                make_external_call(command, opts.debug, opts.debug, True)

          if status>0:
            stdout += " WARNING: could not copy a selected segment file "
            stdout += "from " + thisSearch + " for " + this + "ifo combination"
    except:
      stdout +=" WARNING: problem while copying a selected segment"
      stdout += "(from " + thisSearch + "). passing..."    
      logText(logfile,  stdout, "warning")
      pass
    else: 
      logText(logfile,  stdout)

def copy_files(dataDir,webDir,name):
  '''
  A function to use python to copy files matching name from dataDir to webDir
  '''
  # Loop over matching files
  for file in glob.glob(os.path.join(dataDir,name)):
    # Copy the file
    shutil.copy(file, webDir)

def copy_plot_directory(dataDir,webDir,open_box = True):
  '''
  A function to copy the XXX_summary_plots directory to the html directory
  '''
  # The open_box flag is an old flag that was used to not copy open box plots
  # At the moment this flag is unused but is kept as it may be used again.
  if not os.path.isdir(webDir):
    mkdir(webDir)
    print >> sys.stdout, ' .../ copying all' + dataDir + ' files to ' + webDir
    # Copy the html and cache files
    copy_files(dataDir,webDir,'*html')
    copy_files(dataDir,webDir,'*cache')

    # Copy the images directory
    from_path = os.path.join(dataDir, 'Images')
    if os.path.isdir(from_path):
      to_path = os.path.join(webDir, 'Images')
      if os.path.exists(to_path):
        shutil.rmtree(to_path)
      shutil.copytree(from_path, to_path)
    else:
      # Why does the directory not exist?
      errMsg = "Directory %s does not exist." %(from_path,)
      logText(logfile, errMsg, "warning")

    # Copy the omega_scans directory
    from_path = os.path.join(dataDir, 'omega_scans')
    if os.path.isdir(from_path):
      to_path = os.path.join(webDir, 'omega_scans')
      if os.path.exists(to_path):
        shutil.rmtree(to_path)
      shutil.copytree(from_path, to_path)
    else:
      # Why does the directory not exist?
      errMsg = "Directory %s does not exist." %(from_path,)
      logText(logfile, errMsg, "warning")

    print >> sys.stdout, '... Done'

def symlink_plot_directory(dataDir,webDir):
  '''
  A function to symlink, rather than copy, the plot directory to the html
  directory
  '''
  if not os.path.isdir(webDir):
    if dataDir[-1] == '/':
      dataDir = dataDir[:-1]
    if webDir[-1] == '/':
      webDir = webDir[:-1]
    command = 'ln -s ' + dataDir + ' ' + webDir
    print command
    print >> sys.stdout, ' .../symlinking all' + dataDir + ' files to ' + webDir
    make_external_call(command, opts.debug, opts.debug, True)
    print >> sys.stdout, '... Done'

    
# ***************************************************************************
def fom(page, opts, cachefile_tag=None, caption="fix me",\
	 image_tag=None,directory="playground_summary_plots",all_ifos=False):
  """
  This function reads a cachefile, copy the files to the relevant directory, and  update the HTML document to add figures and pertinent sections. 

  
  @param page: the HTML document page.
  @type page: markup.py document
  @param opts: the main options variable
  @type opts: write_ihope_page options
  @param cachefile_tag: the pattern of the cachefile to look at
  @type cachefile_tag: string
  @param caption: a list of caption
  @type caption: list of string 
  @param image_tag: a list of tag (including *) to select files within a
  @param directory:  the directory to look at
  @type directory: string
  cachefile
  """

  dataDir = opts.datadir+directory+'/'
  webDir = opts.physdir + directory + '/'

  # create a list of cachefiles, each of which containing
  # the cachefile_tag argument in its name
  thisglob = webDir + cachefile_tag +'cache'
  cachefileList =  glob.glob(thisglob)
  
  if opts.verbose is True: 
    print "        Searching for files with this(ese) tag(s): " +str(image_tag)

  # for each cachefile we create a div figure section with all 
  # figures whose name is part of the image_tag list
  cachefileList.sort()
  for eachcache in cachefileList:
    # read this cachefile
    this = open(eachcache, "r")
    fnameList = []
    # Lots of lists but it ensures that the order is always the same
    fileNameList = []
    href = None
    #for each file contained in the cachefile
    if opts.debug is True :
      print >>sys.stdout, "        --> Copying files from " +eachcache

    for filename in this:
      filename = filename.replace("\n","")
      filename = filename.replace(" ","")
      fileNameList.append(filename)
    
    fileNameList.sort()
 
    for filename in fileNameList:
      # if the file is an image, we copy it in ./Images
      # if the file is an image, we copy the file. Then, if the pattern 
      # (in image_tag)  matches the filename, then we create am html section.
      if filename[-4:]=='.png':
        # we check if the pattern match the file
        if image_tag is not None:
          for eachPattern  in image_tag:
            if patternFoundInFilename(filename, eachPattern) is True:
              fnameList.append(filename.replace(".png", "_thumb.png"))
        else:
          fnameList.append(filename.replace(".png", "_thumb.png"))
      elif filename[-5:] == '.html' and ('followup' not in filename):
        href = filename.replace("/pictures", "")
        href = href.replace("\n", "")
        htmlFile = href
        href = directory + '/' + href

    if href:
      msg =" <a href=\"" + href + "\""
      msg+=" rel=\"external\" >"
      msg += """ <br/> --- <b>Click here (to open link in a new tab/window)</b> to get all pictures (full resolution) 
          as well as the pylal arguments used to generate the plots</a> """
      
      source =  eachcache.split('/') 
      source =  source[len(source)-1]
      page = add_figure(page, directory, fnames=fnameList, \
          caption=(caption+' '+msg), source=source,size=None, alt=fnameList,\
          html = htmlFile,all_ifos=all_ifos)
      # close the cachefile
    this.close()
  
#  page.div.close()

  return page
# ***************************************************************************

def html_and_plot(page, opts, cachefile_tag=None,\
          caption='fix_me',image_tag=None,\
          directory=None,html_tag=None,all_ifos=False):
  """
  This function does fom and html_insert together!
  """
  dataDir = opts.datadir+directory+'/'
  webDir = opts.physdir + directory + '/'
  # create a list of cachefiles, each of which containing
  # the cachefile_tag argument in its name
  thisglob = webDir + cachefile_tag +'cache'
  cachefileList =  glob.glob(thisglob)
  # And do the same for the html_tag glob.
  thisglob = webDir + html_tag + 'html'
  htmlfileList = glob.glob(thisglob)

  # Now we need to match them up
  fileListDict = []
  cacheFiles = {}
  htmlFiles = {}
  cachefileList.sort()
  for file in cachefileList: 
    ifo = None
    for combo in ['H1','H2','L1','V1','H1L1','H1H2','H2L1','H2V1','H1V1','L1V1','H1L1V1','H1H2L1','H2L1V1','H1H2V1','H1H2L1V1']:
      if os.path.basename(file).startswith(combo):
        ifo = combo
    if ifo:
      cacheFiles[ifo] = file
      fileListDict.append(ifo)
    else:
      print >> sys.stderr,"WARNING: file doesn't seem to match an ifo combo!"

  for file in htmlfileList:
    ifo = None
    for combo in ['H1','H2','L1','V1','H1L1','H1H2','H2L1','H2V1','H1V1','L1V1','H1L1V1','H1H2L1','H2L1V1','H1H2V1','H1H2L1V1']:
      if os.path.basename(file).startswith(combo):
        ifo = combo
    if ifo:
      htmlFiles[ifo] = file
    else:
      print >> sys.stderr,"WARNING: file doesn't seem to match an ifo combo!"
    if ifo not in cacheFiles.keys():
      print >> sys.stderr,"WARNING: file doesn't seem to correspond to cache file list"

  if opts.verbose is True:
    print "        Searching for files with this(ese) tag(s): " +str(image_tag)

  for ifo in fileListDict:
    if ifo in htmlFiles.keys():
      # read the image cachefile
      this = open(cacheFiles[ifo], "r")
      fnameList = []
      # Lots of lists but it ensures that the order is always the same
      fileNameList = []
      href = None
      #for each file contained in the cachefile

      if opts.debug is True:
        print >>sys.stdout, "        --> Copying files from " +cacheFiles[ifo]

      for filename in this:
        filename = filename.replace("\n","")
        filename = filename.replace(" ","")
        fileNameList.append(filename)

      fileNameList.sort()

      for filename in fileNameList:
        # if the file is an image, we copy it in ./Images
        # if the file is an image, we copy the file. Then, if the pattern 
        # (in image_tag)  matches the filename, then we create am html section.
        if filename[-4:]=='.png':
          # we check if the pattern match the file
          if image_tag is not None:
            for eachPattern  in image_tag:
              if patternFoundInFilename(filename, eachPattern) is True:
                fnameList.append(filename.replace(".png", "_thumb.png"))
          else:
            fnameList.append(filename.replace(".png", "_thumb.png"))
        elif filename[-5:] == '.html' and ('followup' not in filename):
          href = filename.replace("/pictures", "")
          href = href.replace("\n", "")
          htmlFile = href
          href = directory + '/' + href

      if href:
        msg =" <a href=\"" + href + "\""
        msg+=" rel=\"external\" >"
        msg += """ <br/> --- <b>Click here (to open link in a new tab/window)</b> to get all pictures (full resolution) 
            as well as the pylal arguments used to generate the plots</a> """
  
        source =  cacheFiles[ifo].split('/')
        source =  source[len(source)-1]
        page = add_figure(page, directory, fnames=fnameList, \
            caption=(caption+' '+msg), source=source,size=None, alt=fnameList,\
            html = htmlFile,html_file = htmlFiles[ifo],all_ifos=all_ifos)
        # close the cachefile
        this.close()

  return page
  
###############################################################################

def html_insert(page,opts, html_tag=None,\
          caption='fix_me',directory='playground_summary_plots',all_ifos=False):
  """
  This function copies a set of html files verbatim into the page. 

  
  @param page: the HTML document page.
  @type page: markup.py document
  @param opts: the main options variable
  @type opts: write_ihope_page options
  @param html_tag: the pattern of the html files to look at
  @type html_tag: string
  @param caption: a list of caption
  @type caption: list of string 
  @param directory:  the directory to look at
  @type directory: string
  """

  dataDir = opts.datadir+directory+'/'
  webDir = opts.physdir + directory + '/'

  # create a list of html files to be used
  thisglob = webDir + html_tag +'html'
  htmlfileList =  glob.glob(thisglob)

  if opts.verbose is True:
    print "        Copying this html file : " +str(html_tag)

  # for each cachefile we create a div figure section with all 
  # figures whose name is part of the image_tag list
  global fig_num
  htmlfileList.sort()
  for htmlfile in htmlfileList:
    # read this file
    this = open(htmlfile, "r")
    page.add("<!-- insert a figure -->\n<div class=\"figure\">")
    this_count = 0
   # Determine the ifo
    ifo = None
    for combo in ['H1','H2','L1','V1','H1L1','H1H2','H2L1','H2V1','H1V1','L1V1','H1L1V1','H1H2L1','H2L1V1','H1H2V1','H1H2L1V1']:
      if os.path.basename(htmlfile).startswith(combo):
        ifo = combo
    if ifo and not all_ifos:
      page.h5('When ' + ifo + ' interferometers operating')
    # Add the html file into the page
    for line in this:
      page.add(line)
    # add a yellow box (depends on yourcss style though) that contains the 
    # html filename.
    section_num = str(h2_num)+"."+str(h3_num-1)+"."+str(fig_num)
    page.add("\t<p><b>Figure " + section_num + "</b>:  " + caption + "</p>")
    fig_num +=1
    page.add("<pre>Source: "+ htmlfile + "\n </pre>" )
    page.div.close()
    # close the file
    this.close()

  return page

# ***************************************************************************
def set_style():
  """
  Function to copy the style file as read from the write_ihope_page.ini file
  """

  tmp = []
  tmp.append(configcp.get("main", "style"))
  try:
    style = configcp.get("main","style")
    command = 'cp ' + style + " " +opts.physdir
    make_external_call(command, opts.debug, opts.debug, True )
    tmp=style.split('/')
    return tmp[len(tmp)-1]
    
  except:
    print sys.stderr()<< 'could not copy the style file'
    pass
    return ""

# ----------------------------------------------------------------------------
def parse_arguments():
  """
  Function to parse the arguments and check their validity.
  """
  usage =  """ %prog [options]
  Program to write webpage from upperlimit.py
  """

  parser = OptionParser( usage = usage, version = "%prog CVS "+__Id__ )

  parser.add_option("-C","--config-file",action="store",type="string",\
      metavar=" INI File",\
      help="""ini file with information about run directories. The ini file should look like
-------                        
[main].............................................
gps-start-time  = 847555570........................
gps-end-time    = 849974770........................
title           = \"Low mass CBC analysis\"........
ihope-ini-file  = ihope.ini........................
ihope-directory = /archive/home/cokelaer/S5/Month1/full_analysis/....
home-directory  = /archive/home/...................
url             = ldas-jobs.ligo.caltech.edu.......
username        = cokelaer.........................
style           = /archive/home/cokelaer/style.css.
output          = index.shtml""")
  parser.add_option("-A","--open-the-box",action="store_true",\
      default=False,dest="full_data", metavar="DOANALYSIS",\
      help="" )
  parser.add_option("-T","--symlink-plots",action="store_true",\
      default=False, metavar="DOTUNING",\
      help="" )
  parser.add_option("-S","--skip-summary",action="store_false",\
      default=True,dest="summary", metavar="SOSUMMARY",\
      help="Skip generating brief summary pages")
  parser.add_option("-f","--skip-summary-files",action="store_false",\
      default=True,dest="summary_files", metavar="SOSUMMARY",\
      help="Skip generating the summary files section.")
  parser.add_option("-U","--skip-upperlimit",action="store_false",\
      default=True,dest="upperlimit", metavar="DOUPPERLIMIT",\
      help="" )
  parser.add_option("-I","--skip-injection",action="store_false",\
      default=True,dest="injection", metavar="DOUPPERLIMIT",\
      help="" )
  parser.add_option("-H","--skip-hardware-injection",action="store_false",\
      default=True,dest="hardware_injection", metavar="SKIPHWINJ",\
      help="" )
  parser.add_option("-P","--skip-playground",action="store_false",\
      default=True,dest="playground", metavar="DOPLAYGROUND",\
      help="" )
  parser.add_option("-B","--skip-full-data-slide",action="store_false",\
      default=True,dest="full_data_slide", metavar="DOPLAYGROUND",\
      help="" )
  parser.add_option("-F","--skip-followup",action="store_false",\
      default=True,dest="followup", metavar="DOFOLLOWUP",\
      help="" )              
  parser.add_option("-N","--no-copy",action="store_true",\
       default=False,dest="no_copy",metavar="NOCOPY",\
       help="do not copy any file (useful for debugging only) ")
  parser.add_option("-D","--debug",action="store_true",\
       default=False,dest="debug",metavar="DODEBUG",\
       help="More verbose than --verbose ")
  parser.add_option("-V","--verbose",action="store_true",\
      default=False, dest="verbose",metavar="VERBOSE",\
      help="talkative script! One can also add the debug option." )

  (opts,args) = parser.parse_args()

  if opts.config_file is None:
    raise ValueError

  return opts,args

# ----------------------------------------------------------------------------

opts,args = parse_arguments()
#############################################################################
#  MAIN PART                                                                #
#############################################################################
fig_num = 1
count_block = 0  # counter for the unique id of the toggle boxes
h2_num = 0       # counter for the section h2
h3_num = 1       # counter for the section h3
config   =  opts.config_file
opts.config = config # save the name of the ini file, why ?
configcp = glue.pipeline.DeepCopyableConfigParser()
configcp.read(config)
maxdiv = 200
print >>sys.stdout, "|------------------- Initialisation"
# First, we open an xml file, for the log file
logfile_name = __name__+".xml"
print >>sys.stdout,"Openning the log file (" +logfile_name+")."
logfile = open(logfile_name, "w")
logfile.write("""<?xml version="1.0" encoding="ISO-8859-1"?>
<?xml-stylesheet type="text/xsl" href="write_ihope_page.xsl"?>
<log>""")


#---------------------------------------
# then, we parse the write_ihope_page.ini file
print >>sys.stdout,"Parsing the ini file: " + opts.config
try:
  opts.config_file 	= configcp.get("main", "ihope-ini-file")
  opts.gps_start_time 	= configcp.get("main", "gps-start-time")
  opts.gps_end_time 	= configcp.get("main", "gps-end-time")
  opts.ihope_directory 	= configcp.get("main", "ihope-directory")
  opts.title	 	= configcp.get("main", "title").replace( '"', '' )
  opts.output	 	= configcp.get("main", "output")
  opts.html_directory 	= configcp.get("main", "html-directory")
  Lead                  = configcp.get("main", "lead").replace( '"', '' )
  Second                = configcp.get("main", "second").replace( '"', '' )
  Notes                 = configcp.get("main", "notes").replace( '"', '' )
except:
  print >> sys.stderr, "ERROR::The ini file does not have the proper field in the [main] section" 
  print >> sys.stderr, """       Consider adding one of those fields if missing: ihope-ini-file, \
	gps-start-time,gps-end-time, ihope-directory, title,url, username, output"""
  raise
  
#------------------------------------
#sub-products of the ini file parsing
opts.gpsdir =  '/'+str(opts.gps_start_time)+'-'+str(opts.gps_end_time)+'/'   # note the / at the end
opts.duration = str(int(opts.gps_end_time) - int(opts.gps_start_time))
opts.daydur = float(opts.duration) / 86400
opts.webdir = './'
opts.datadir = os.path.abspath(opts.ihope_directory) + opts.gpsdir 
opts.physdir = opts.html_directory + opts.gpsdir
opts.txt_suffix = '-'+opts.gps_start_time+'-'+opts.duration+'.txt'


#----------------------
# openning the html file
print >>sys.stdout,"Openning the HTML (" + opts.output+")"
try:
  html_file = file(opts.output,"w")
except:
  msg = "Cannot open %" % opts.output
  print >>sys.stderr, msg
  raise


#-----------------------------------------
# here is the directory we want to extract information from
msg = "Entering this directory (where we will get all the relevant information)" + opts.datadir
print >> sys.stdout, msg
if not  os.path.isdir(opts.datadir):
  raise  "%s is not a valid directory. Check your gps time." % opts.datadir
# which physical name is 
msg = "We will copy all images and HTML documents into this directory " +opts.physdir
logText(logfile, msg,"done")
# Make the physdir and all necessary directories above it
try:
  os.makedirs(opts.physdir)
except OSError:
  pass
mkdir(opts.physdir+'/Images')

# get the style file and copy it (must happen after the html directories have been created)
opts.style = set_style()

#-----------------------------------------
# now we can parse the ihope.ini file itself
msg =   "Parsing the ihope ini file (" + opts.config_file+")"
logText(logfile, msg)
hipe   = opts.datadir+'/'+opts.config_file
hipecp = glue.pipeline.DeepCopyableConfigParser()
hipecp.read(hipe)
make_external_call( 'cp '+opts.config_file + ' ' + opts.physdir, False, opts.debug, True)

#-----------------------------------------
# now we copy the segments to the web directory and other relevant information
copy_segments()


#-----------------------------------------
###### create the section labels  
html_sections={}

# the string "full_data" is hardcoded somewhere else, so it must remain as it is
html_order = ['general', 'datainfo', 'fullsummary','play_full_slide_summary',\
    'playsummary','playground', 'full_data', \
    'full_data_slide', 'upperlimit_full','upperlimit_play','injection',\
    'hardware_injections', 'followup', 'upperlimit', \
    'summary_files', 'logfile', 'about']

html_sections['general'] = "General Information"
html_sections['datainfo'] = "Data Information"
defaultPage = 'datainfo'
if opts.summary is True and opts.playground is True:
  html_sections['playsummary'] = "Summary of playground"
  defaultPage = 'playsummary'
if opts.summary is True and opts.full_data_slide is True:
  html_sections['play_full_slide_summary'] = \
      "Summary of playground with full data slides"
  defaultPage='play_full_slide_summary'
if opts.summary is True and opts.full_data is True:
  html_sections['fullsummary'] = "Summary of full data"
  defaultPage = 'fullsummary'
if opts.playground is True: html_sections['playground']	= "Playground"
if opts.full_data_slide is True: 
  html_sections['full_data_slide'] = "Full Data Slide And Playground Zero lag"
if opts.injection is True: 
  if hipecp.has_section("injections"):
    html_sections['injection'] = "Injection"
    for inj in hipecp.options("injections"):
      html_sections[inj] = inj.upper() 
    html_sections['allinj'] = 'All Injections Combined'
  else:
    opts.injection = False
    print >> sys.stderr, "No injections section found in ini file."
    print >> sys.stderr, "Now running with --skip-injections"
if opts.hardware_injection is True: 
  html_sections['hardware_injections'] = "Hardware Injection"
if opts.full_data is True: html_sections['full_data'] = "Full Data"
if opts.followup is True: html_sections['followup'] = "Followup"
if opts.upperlimit is True: 
  html_sections['upperlimit_play'] = "Search Sensitivity"
if opts.upperlimit and opts.full_data: 
  html_sections['upperlimit_full'] = "Upper Limits"
if opts.summary_files:
  html_sections['summary_files'] = "Summary files"
html_sections['logfile'] = "Log File"
html_sections['about'] = "About"

title = opts.title+ ": from "+str(opts.gps_start_time)\
  +" to "+str(opts.gps_end_time) 
script = {}
script['toggle.js'] = 'javascript'
script['//ajax.googleapis.com/ajax/libs/jquery/1.2.6/jquery.min.js']\
    = 'javascript'
create_toggle()
# Finally, we create the html document 
msg =   "Creating HTML document"
logText(logfile, msg)

mainpage = initialize_page(title,opts.style,script)

mainpage.div(id_="wrapper")
mainpage.div(id_="menubar")
mainpage.div(id_="menu")
# The main part of the code, which loops over all subsections
subPages = {}
for each_section in html_order:
  if each_section in html_sections:
    if each_section == 'injection':
      injs = hipecp.options("injections")
      injs.sort()
      for inj in injs:
        subPages[inj]=\
            write_sub_page(opts,inj,html_sections,opts.style,script)
      subPages[inj]=\
          write_sub_page(opts,'allinj',html_sections,opts.style,script)
    else:
      subPages[each_section]=\
          write_sub_page(opts,each_section,html_sections,opts.style,script)
mainpage.div.close()
mainpage.div.close()
mainpage.div(id_="ihope")
mainpage.add('<h2> i h o p e </h2>')
if opts.full_data:
  mainpage.add('<img src="https://www.lsc-group.phys.uwm.edu/ligovirgo/cbc/public/segments/S5/OpenBox.jpg">')
else:
  mainpage.add('<img src="https://www.lsc-group.phys.uwm.edu/ligovirgo/cbc/public/segments/S5/thomasLegacy.jpg">')
mainpage.div.close()
mainpage.div(id_='header')
mainpage.add('<h1>' + title  +' </h1>')
leads = 'Lead: ' + Lead + '&nbsp;&nbsp&nbsp Second: '+ Second
mainpage.add('<h3> ' + leads + ' </h3>')
mainpage.add('<h3>' + Notes + ' </h3>')
mainpage.div.close()
mainpage.div(id_='iframecontent')
mainpage.add('<p id="placeholder">Please select a report section on the left.</p>')
mainpage.add('<script type="text/javascript">')
mainpage.add("""loadFrame('""" + defaultPage + """.html')""")
mainpage.add('</script>')
mainpage.div.close()
mainpage.div.close()

# what is the False for ? I cannot remenber
html_file.write(mainpage(False))
html_file.close()
# just remove the <html lang="en"> line tat is not a strict HTML code
cmd = 'sed -e \'s/<html lang=\"en\">//\' '+ opts.output + '>/tmp/tmp.html'
output,status = make_external_call(cmd, True, True, False)
if status==0:
  cmd = 'mv /tmp/tmp.html '+opts.output
  make_external_call(cmd, True, True, True)
if opts.full_data:
  cmd = 'sed -i \'s/<body>/<body\ class=\"openbox\">/g\' ' + opts.output + '>/tmp/tmp.html'
  make_external_call(cmd, True, True, True)

# that's it for the html creation. let us copy it to the requested directory  

print '---------------------FINISHED ---------------------'
print '--- HTML file created. '
print '--- Copying html documents in ' +opts.physdir
make_external_call('mv  '+opts.output +' ' + opts.physdir, opts.debug, opts.debug, True)
make_external_call( 'mv toggle.js '+ opts.physdir, opts.debug, opts.debug,  True)

#let us close the log file
logfile.write("</log>")
logfile.close()
logfile = __name__+".xml"
output, status = make_external_call( 'grep WARNING '+ logfile +'| wc - | awk \'{print $1}\' - ', opts.debug, opts.debug, True)

if status==0:
  if int(output)==0:
    print 'No warnings'
  else:
    print '\n\n\nThere are warnings : '+str(int(output))+' . Check the log file '+logfile
  
  output, status = make_external_call('mv '+logfile + " "+opts.physdir, True,True,True) 
else:
  print 'Could not find the log file ' +logfile
  
 
#Finally create the xsl for the log xml file
logfile = open(__name__+".xsl", "w")
logfile.write("""<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

<xsl:template match="/">
  <html>
  <body>
  <center><h1>Log file summary </h1></center>
  <xsl:apply-templates/>
  </body>
  </html>
</xsl:template>

<xsl:template match="section">
<h2 color="blue">Section: <xsl:value-of select="."/></h2>
<br />
</xsl:template>

<xsl:template match="done">
<center>
  <div>
    <table bgcolor="green" border="2" width="80%">
      <tr>
        <td width="80%"><xsl:value-of select="."/></td>
        <td bgcolor="white"> passed</td>
      </tr>
    </table>
  </div>
</center>
</xsl:template>
<xsl:template match="warning">
<center>
  <div>
    <table bgcolor="orange" border="2" width="80%">
      <tr>
        <td width="80%"><xsl:value-of select="."/></td>
        <td bgcolor="white"> fix me</td>
      </tr>
    </table>
  </div>
</center>
</xsl:template>
<xsl:template match="error">
<center>
  <div>
    <table bgcolor="red" border="2" width="80%">
      <tr>
        <td width="80%"><xsl:value-of select="."/></td>
        <td bgcolor="white"> skipped</td>
      </tr>
    </table>
  </div>
</center>
</xsl:template>



</xsl:stylesheet>
""")

logfile.close()
output, status = make_external_call('mv '+__name__+".xsl" + " "+opts.physdir, True,True,True) 


