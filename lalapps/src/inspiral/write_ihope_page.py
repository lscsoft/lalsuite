#!/usr/bin/python

# $Id$
__author__ = "Thomas Cokelaer <thomas.cokelaer@astro.cf.ac.uk>"
__version__ = "$Revision$"
__date__ = "$Date$"
__Id__ = "$Id$"
__name__ = "write_ihope_page"

import sys
import copy
import os
from optparse import *
import glob
from glue  import markup
from glue.markup import oneliner as e
import ConfigParser



def functionId(nFramesUp):
    """ Create a string naming the function n frames up on the stack.
    """
    try:
      co = sys._getframe(nFramesUp+1).f_code
      msg = "%s (%s @ %d)" % (co.co_name, co.co_filename, co.co_firstlineno)
      if msg.startswith("?") is False:
        print "-->ERROR in function: " + msg
    except:
      msg=[]
    return msg

# ***************************************************************************
# ***************************************************************************
def logText(logfile, text, tag="done", debug=True):
  
  if tag=="warning":
    msg= "<"+tag+">"+text+"</"+tag+">\n"
    logfile.write(msg)
    print >>sys.stdout,text

  elif tag=="error":
    msg = "<"+tag+">"+text
    logfile.write(text)
    i =1 
    while len(msg)>0:
      msg = functionId(i)
      if msg.startswith("?") is False:
        logfile.write("\n"+msg)
      else :
        logfile.write("</"+tag+">\n")
        msg=[]
      i=i+1

  else :
    msg = "<"+tag+">"+text+"</"+tag+">\n"
    logfile.write( msg )
    print >>sys.stderr,text


def make_external_call(command, show_stdout=False, show_command=False, show_error=True):
  """
  Run a program on the shell and print informative messages on failure.
  """
  if show_command: print command

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
      print  this_output[1:]  #first character is \n

  stdin.close()
  out.close()
  err.close()  
  return this_output, status

# ***************************************************************************
# ***************************************************************************
def mkdir( newdir ):
  """
  Create a directory, don't complain if it exists

  newdir = name of directory to be created
  """
  if os.path.isdir(newdir): pass
  elif os.path.isfile(newdir):
    raise OSError("a file with the same name as the desired " \
                  "dir, '%s', already exists." % newdir)
  else: os.mkdir(newdir)


# ***************************************************************************
# ***************************************************************************
def create_toggle():
  fname = open("toggle.js", "w")
  fname.write("function toggleVisible(division) {\n"+\
    "  if (document.getElementById(\"div_\" + division).style.display == \"none\") {\n"+\
    "    document.getElementById(\"div_\" + division).style.display = \"block\";\n"+\
    "    document.getElementById(\"input_\" + division).checked = true;\n"+\
    "  } else {\n"+\
    "    document.getElementById(\"div_\" + division).style.display = \"none\";\n"+\
    "    document.getElementById(\"input_\" + division).checked = false;\n"+\
    "   }\n"+\
    " }\n"
   )
  fname.close()

def write_results(page, opts, section):
  """
  """
  print >>sys.stdout,"Creating section "+section +"..."
  if section=='general': 	page = write_general(page, opts)
  elif section=='toc': 		page = write_toc(page, opts)
  elif section=='summary': 	page = write_summary(page, opts)
  elif section=='playground': 	page = write_analysis(page, opts,thisSearch='playground')
  elif section=='tuning': 	page = write_tuning(page, opts)
  elif section=='injection': 	page = write_injection(page, opts)
  elif section=='analysis': 	page = write_analysis(page, opts,thisSearch='full_data')
  elif section=='upperlimit': 	page = write_upperlimit(page, opts)
  elif section=='logfile': 	page = write_logfile(page, opts)
  elif section=='about': 	page = write_about(page, opts)
  else:
    msg = "this section ("+ section +") doest not exist. "
    logText(logfile, msg, "error")

  return page
# ***************************************************************************
def write_logfile(page , opts):
  """
  This function creates an entry with the logfile information.
 
  @param page:
  @param opts:
  """
  # get the directory of the url
  dir = opts.webdir +'/'

  page.p(e.br())
  page.add("<!-- beginning of the general section -->")
  page.div(class_="encadre")
  page.div(class_="contenu")
  page = write_title(page, html_sections['logfile'], "zola")
  page.div(id="div_zola", style='display:none')
  page.a("logfile", href=dir +  __name__+'.xml')
  page.div.close()
  page.div.close()
  page.div.close()

  return page
# ***************************************************************************
def write_toc(page , opts):
  """ 
  This function creates the table of contents, using the html_sections 
  variable defined in the main part of this program.

  @param page:
  @param opts:
  """
  # get the directory of the url
  dir = opts.webdir +'/'

  # and the section names
  items = html_sections 
  order = html_order 
  # create a frame for the toc, and add contents 
  page.add("<!-- beginning of toc section-->")
  page.div(class_="encadre")
  page.div(class_="contenu")
  page.h2()
  page.add("<a name=\"toc\"> Table of Contents  </a>")
  page.h2.close()
  page.ol()
  for this_item in order:
    try:
      page.li(e.a(this_item, href=dir + opts.output + '#' + items[this_item]))
      msg = "added "+this_item
      logText(logfile, msg, "done")
    except:
      msg = 'skip '+this_item
      logText(logfile, msg, "warning")
#      pass

  page.ol.close()
  page.div.close()
  page.div.close()

  return page

# ***************************************************************************
# ***************************************************************************
def write_general(page,opts):
  webdir = opts.webdir
  ini = opts.config_file

  page.p(e.br())
  page.add("<!-- beginning of the general section -->")
  page.div(class_="encadre")
  page.div(class_="contenu")
  page = write_title(page, html_sections['general'], "rimbaud")
  page.div(id="div_rimbaud", style='display:none')
  text=  'This page summarizes the analysis of the data segment from GPS time %s up to %s' % (opts.gps_start_time,  opts.gps_end_time)
  page.add(text)

  # The search parameters 
  try:
    executables = ("inspiral", "tmpltbank", "sire", "thinca", "trigbank", "coire" , "inspinj")
    page.h3("The search used the following resources:")
    
    page.ol() #starts enumerate
    page.add("<li> LAL/LALApps "+get_version("inspiral")) #first section
    page.ul()
    for exe in executables:
      text = "<b>lalapps_"+exe+ "</b>    "+get_version(exe)
      page.li(text)   
    page.ul.close()
    page.li.close()

    page.add("<li>Segment information:") #section section
    page.ul()
    for this in get_ifos():
        seg = this +'-SELECTED_SEGS-'+opts.gps_start_time+'-'+opts.duration+'.txt'
        this_ref = webdir + '/segments/'+seg
        page.li(e.a(seg,href=this_ref) )
    page.ul.close()
    page.li.close()

    text=("The configuration is contained in the file <a href=\"" + \
        webdir + "/" + ini + "\">" + ini + "</a>")
    page.li(text) #third section
    text = "A list of category files stored in this directory ("  \
   	+ "<a href=\""+ webdir + "/catlists/\"> catlists</a>"  +\
  	"), and listed in the ini file."

    page.li(text) # fourth section
    page.ol.close()
  except:
    logText(logfile, """Problem with the executable: cannot find them ? """, "warning")
    pass
  else:
    logText(logfile,  "...Get the executables version...")


  # The ifo requested
  page = heading(page, "This search concerned the following combinaison of ifos")
  page = add_config_section(page, "ifo-details")
  page.div.close()


  page.div.close()
  page.div.close()
  page.div.close()
  
  return page


# ***************************************************************************
# ***************************************************************************
def write_title(page, text, tag):
  page.h2()
  page.add("<a name=\""+text.replace(" ", "_")+"\">")
  input = "input_"+tag
  page.input(id=input, type="checkbox", checked="checked", onclick="toggleVisible('"+tag+"')")
  page.add(text)
  page.add("</a>")
  page.h2.close() 
  page.div(e.a("return to top", href="#toc" ), class_="back")

  return page


# ***************************************************************************
# ***************************************************************************
def write_summary(page,opts):
  # first, get some information 
  webdir = opts.webdir
  datadir = opts.datadir
  ini = opts.config_file
  ifos = get_ifos()

  page.p(e.br()) # to put a space between different section
  page.add("<!-- beginning data summary section-->")
  page.div(class_="encadre") #(1)
  page.div(class_="contenu") #(2)
  # title with a div to toggle on/off the content 
  page = write_title(page, html_sections['summary'], "verlaine")
  
  page.div(id="div_verlaine", style='display:none') #(3)
  page = heading(page, "Selected segments", "Switch details on/off")  #(4)
  page.p("""The segments files above were created with no data quality flags 
set. The times analyzed according to hipe are (SELECTED_SEGS):""")
  page.p(e.br())
  page.table()
  segs = get_segments_tag("SELECTED_SEGS")
  keys = ("segments","H1","H2","L1","G1","V1","T1") 
  for key in keys:
    if segs.has_key(key):
      page.tr()
      page.td(key)
      page.td(segs[key])
      page.tr.close()
  page.table.close()
  page.div.close() #(3)

  i=0
  for cat in ("CATEGORY_1_VETO_SEGS","CATEGORY_2_VETO_SEGS","CATEGORY_3_VETO_SEGS","CATEGORY_4_VETO_SEGS"):
    i=i+1
    try:
      page = heading(page, "Category " + str(i) + " veto segments")  #(4)
      
      page.add("This category includes the following flags : ")
      for ifo in ifos:
        command = 'awk \'{print $1}\' ' +  datadir +'/segments/' + ifo + "cat"+str(i)+".txt"
        flags, status = make_external_call(command, opts.debug, opts.debug, True)
        page.pre(flags)

      page.add("The veto times (category "+str(i)+") according to hipe ( "+cat+")")
      page.table()
      segs = get_segments_tag(cat)
      keys = ("segments","H1","H2","L1","G1","V1","T1") 
      for key in keys:
         if segs.has_key(key):
           page.tr()
           page.td(key)
           page.td(segs[key])
           page.tr.close()
      page.table.close()
      page.div.close() # (3)

    except:
      logText(logfile, "Problems parsing category veto segment list "+str(i), "warning")
     

  for tag in ["SCIENCE_SEGMENTS", "RDS_C03_L2"] :
    page = heading(page, "Science segments ") #(4)
    page.add("The science segments (category 1) according to hipe are "+tag+": <br/>")
    page.table()
    segs = get_segments_tag(tag)
    keys = ("segments","H1","H2","L1","G1","V1","T1") 
    for key in keys:
      if segs.has_key(key):
        page.tr()
        page.td(key)
        page.td(segs[key])
        page.tr.close()
    page.table.close()
    page.div.close() #(3)


  page.div.close() #(2)
  page.div.close() #(1)
  page.div.close() #(0)
  
  return page




# ***************************************************************************
# ***************************************************************************
def write_injection(fname,dir):
  fname.write("<a name=\"eff\"><h2>Efficiency Plots</h2></a>\n\n")

  for ifos in opts.ifos:
    fname.write("<h3>" + ifos + "</h3>")

    if "L" in ifos:
      plot_name = "plotnumgalaxies/" + ifos + "/" + ifos + "_" + \
          opts.x_value + "_" + opts.y_value
      logText(logfile, plot_name, "warning")
      add_figure(fname,dir,fnames=[plot_name + "_efficiency.png",\
          plot_name + "_mc_errors.png"],\
          size="half",caption="efficiency plot, and Monte Carlo errors")

    else:
      plot_name = "plotnumgalaxies/" + ifos + "/" + ifos + "_" + opts.x_value
      add_figure(fname,dir,fnames=[plot_name + "_efficiency.png"],\
          size="full",caption="efficiency plot")

      add_figure(fname,dir,fnames=[plot_name + "_hanford_calibration.png",\
        plot_name + "_monte_carlo.png",plot_name + "_waveform.png"],\
          size="third",caption="efficiency plots with systematic errors")

# ***************************************************************************
# ***************************************************************************
def write_upperlimit(page, opts):
  
  webdir = opts.webdir
  ini = opts.config_file
  ifos = get_ifos()
  mkdir(opts.physdir+"/playground")
  page.p(e.br())
  page.div(class_="encadre")
  page.div(class_="contenu")
  page = write_title(page, html_sections['upperlimit'], "zola")
  page.div(id="div_zola", style='display:none')
  page.add("Upper Limit Results")

  # gaussian
  page = heading(page, "Gaussian Mass Distribution", "see details here", "h3")
  try:
    page.pre()
    page.add("Compute Posterior output")
    tmp_file = file(opts.datadir + "analysis/plotnumgalaxies/" + "ul-output-gaussian.log")
    page.add(tmp_file.read())
    fname.write("</pre>")

    page = add_figure(fname,dir,fnames=["plotnumgalaxies/Gaussian-posterior-pdf.png", \
      "plotnumgalaxies/Gaussian-posterior-cdf.png"], \
      size="half",\
      caption="posterior and cumulative posterior for the rate")
  except:pass
  page.div.close()

  # uniform
  page = add_input_h3(page, "Uniform Mass Distribution")
  
  page.add("Compute Posterior output:\n")
  #tmp_file = file(opts.actual_directory + "/plotnumgalaxies/" + \
  #    "ul-output-uniform.log")
  #fname.write(tmp_file.read())

  try:
    add_figure(fname,dir,fnames=[ "plotnumgalaxies/Uniform-rate-v-mass.png"], \
      size="full", caption="rate upper limit vs total mass")
    uniFiles = glob.glob(opts.actual_directory + \
      "/plotnumgalaxies/Uniform-*-mtotal*pdf.png")
    for thisFile in uniFiles:
      # keep only the file name
      thisFile = thisFile.split("/")[-1]
      # get the mass range
      mlo = thisFile.split("-")[1]
      mhi = thisFile.split("-")[2]
      add_figure(fname,dir,fnames=["plotnumgalaxies/" + thisFile,\
          "plotnumgalaxies/" + thisFile.strip("pdf.png") + "cdf.png"],\
          size="third", \
          caption="posterior and cumulative posterior in mass range " + \
          mlo + " to " + mhi)
  except:pass
  page.div.close()# close the h3
  
  #---------------------------- loudest events
  page = add_input_h3(page, "loudest events")
  page.div.close()# close the h3

  # the end
  page.div.close()
  page.div.close()
  page.div.close()
  return page

# ***************************************************************************
# ***************************************************************************
def write_tuning(page,opts):
  webdir = opts.webdir
  ini = opts.config_file
  ifos = get_ifos()
  mkdir(opts.physdir+"/playground")
  page.p(e.br())
  page.div(class_="encadre")
  page.div(class_="contenu")
  page = write_title(page, html_sections['tuning'], "ronsard")
  page.div(id="div_ronsard", style='display:none')
  page.add("This section summarizes the tuning with e-thinca, r-sqaure, chi-square and h1-h2 cuts.<br/>")
  
  #---------------------------- EThinca
  page = heading(page, "E-thinca for injections") #(4)
  for injtag in ['bns_inj', 'bhns_inj', 'bbh_inj']:
    page.h4(injtag)
    for combi in ["H1H2L1", "H1H2", "H1L1", "H2L1"]:
      #lists = glob.glob(opts.datadir+"/"+ injtag +"/"+ifo+"-SIRE_INJECTIONS_*FOUND-*.txt")
#      summary = open(lists[0])
      caption = combi +" E-thinca parameter at the first stage of the coincidence analysis ("+injtag +"tag)."      
#      caption = caption + "<pre> -----------------------"+lists[0]+"\n" +summary.read()+"</pre>"
      if combi=="H1H2L1":
        page = add_figure(page, fnames=  [\
            "injections/ethinca/"+combi+"_first_coinc_"+injtag+"_H1_H2_ethinca_vs_H1_snr.png", \
            "injections/ethinca/"+combi+"_second_coinc_"+injtag+"_H1_H2_ethinca_vs_H1_snr.png", \
            "injections/ethinca/"+combi+"_second_coinc_"+injtag+"_H1_H2_ethinca_vs_H1_snr_andtotMass.png", \
            "injections/ethinca/"+combi+"_first_coinc_"+injtag+"_H1_L1_ethinca_vs_H1_snr.png", \
            "injections/ethinca/"+combi+"_second_coinc_"+injtag+"_H1_L1_ethinca_vs_H1_snr.png", \
            "injections/ethinca/"+combi+"_second_coinc_"+injtag+"_H1_L1_ethinca_vs_H1_snr_andtotMass.png", \
            "injections/ethinca/"+combi+"_first_coinc_"+injtag+"_H2_L1_ethinca_vs_H2_snr.png", \
            "injections/ethinca/"+combi+"_second_coinc_"+injtag+"_H2_L1_ethinca_vs_H2_snr.png",\
            "injections/ethinca/"+combi+"_second_coinc_"+injtag+"_H2_L1_ethinca_vs_H2_snr_andtotMass.png"], \
            caption=caption, size="third",alt=["","","","","","","", "", ""])
      else:
        page = add_figure(page, fnames=  [\
            "injections/ethinca/"+combi+"_first_coinc_"+injtag+"_"+combi[0:2]+"_"+combi[2:4]+"_ethinca_vs_"+combi[0:2]+"_snr.png", \
            "injections/ethinca/"+combi+"_second_coinc_"+injtag+"_"+combi[0:2]+"_"+combi[2:4]+"_ethinca_vs_"+combi[0:2]+"_snr.png",\
            "injections/ethinca/"+combi+"_second_coinc_"+injtag+"_"+combi[0:2]+"_"+combi[2:4]+"_ethinca_vs_"+combi[0:2]+"_snr_andtotMass.png"],\
            caption=caption, size="third",alt=["", "", ""])
  page.div.close()


  page = heading(page, "Chi-square tuning") #(4)
  page.div.close()
  page = heading(page, "R-square tuning") #(4)
  page.div.close()
  page = heading(page, "H1H2 distance cut") #(4)
  page.div.close()


  # the end
  page.div.close()
  page.div.close()
  page.div.close()
  return page


# ***************************************************************************
# ***************************************************************************
def write_analysis(page,opts, thisSearch='playground'):
  
  webdir = opts.webdir
  ini = opts.config_file
  ifos = get_ifos()
  mkdir(opts.physdir+"/"+thisSearch)


  page.p(e.br())
  page.div(class_="encadre")
  page.div(class_="contenu")
  if thisSearch=='playground':
    page = write_title(page, html_sections[thisSearch], "hugo")
    images_dir='/playground_summary_plots/'
  else:
    page = write_title(page, html_sections[thisSearch], "hugo")
    images_dir='/analysis_summary_plots/'
  page.div(id="div_hugo", style='display:none')
  page.add("This section summarizes the analysis of the "+thisSearch+" data.<br/>")
  
  #table and venn diagramm
  
  try:
    page = heading(page, "General information", "see details here")
    segs = get_segments_tag('SELECTED_SEGS')
    print segs
    keys = ("segments","H1","H2","L1","G1","V1","T1") 
    page.p("The segment files above were created with no data quality flags set")
    page.p("The times analyzed accoring to hipe are:")
    coincs = get_coincident_segments(thisSearch)  
    ## here is the table. The first column is another table with the durations, 
    ## and the second one is the venn diagram  
    page.add("<table><tr><td>\n")
    page.table()
    page.tr();
    page.td('coincidence'); 
    page.td('duration(s)'); 
    page.tr.close()
    for key in coincs.keys():        
      page.tr()
      if thisSearch=='playground':
        file = webdir +'/'+thisSearch +'/'+key+"_play_segs_analyzed.txt"
      else:
        file = webdir +'/'+thisSearch +'/'+key+"_segs_analyzed.txt"
      page.td(e.a(key,href=file))
  
      page.td(coincs.get(key))
      page.tr.close()
    page.table.close()
    page.add("</td>\n")
  except:
    logText(logfile, "Error in generating analysed segments.", "error")
    pass

  try:
    page.add("<td>")
    #create figure for the venn diagram
    data = ""
    for coinc in ("H1","H1H2","H2","H2L1","L1","H1L1","H1H2L1"):
        data = data +coincs.get(coinc) + " "
    create_venn(data, thisSearch)
    # and incorporate into html
    comment = "Venn diagramm showing distribution of"
    for coinc in  coincs.keys():
      comment = comment + " "+ coinc
    page = add_figure(page, fnames =[thisSearch+"/venn_"+thisSearch+".png"], caption=comment, size="half")
    page.add("</td></tr></table>")
    page.add("here below is the detail of the data and ligo-data section of ihope.ini ")
    page = add_config_section(page, "data")
    page = add_config_section(page, "ligo-data")
    page.div.close()  
  except:
    logText(logfile, "Error in generating Venn diagram.", "error")
    pass
 
  
  #----------------------------  The horizon distance  
  page = heading(page, "Inspiral range plots" )
  caption = """ "Inspiral Horizon distance for a \
       (1.4,1.4) solar mass system with SNR=8 (first sub-figure), and \
       histograms(second sub-figure). The last sub-figure shows the \
       expected horizon distance for any total mass, using an SNR=8."""  
  tag = '*inspiralrange*'
  page = fom(page, opts, tag, caption, images_dir)
  page.div.close()

  #---------------------------- the number of templates
  page = heading(page, "Variation in template bank and triggered template \
	bank size")
  page = add_config_section(page, "tmpltbank")
  tag = '*plotnumtemplate*'
  caption = "Variation in template bank and triggered template bank size"
  page = fom(page, opts, tag, caption, images_dir)
  page.div.close()
  
  #---------------------------- the first inspiral stage (no veto)
  page = heading(page, "First inspiral stage")
  page = add_config_section(page, "inspiral")
  tag = '*plotinspiral_FIRST_*PLAYGROUND*'
  caption = "Trigger rate at first inspiral stage"
  page = fom(page, opts, tag, caption, images_dir)
  page.div.close()
  
  #---------------------------- the second  inspiral stage (no veto)
  page = heading(page, "Second inspiral stage")
  page = add_config_section(page, "inspiral")
  tag = '*plotinspiral_SECOND*PLAYGROUND*'
  caption = "Trigger rate at second inspiral stage"
  page = fom(page, opts, tag, caption, images_dir)
  page.div.close()
  #---------------------------- the first inspiral stage (no veto)
  page = heading(page, "Second thinca step (all ifo combinaison)")
  page = add_config_section(page, "thinca")
  tag = '*plotthinca_SECOND*PLAYGROUND*'
  caption = "Trigger rate at second inspiral stage"
  page = fom(page, opts, tag, caption, images_dir)
  page.div.close()

  

  
  #---------------------------- the inspiral trigger rate, category 1 and 2
#  page.p("The parameters for the inspiral jobs are ")
#  page = add_config_section(page, "inspiral")
#  page.p("We can look more closely at the first inspiral stage results, and in \
#	particular at the trigger distribution (clustered, category veto 1 and 2")
#  # add figures and summary files
#  for ifo in get_ifos():
#    try :
#      page.h4(ifo)
#      summary = open(opts.datadir+"/"+ thisSearch  +"/"+ifo+"-SIRE-"+opts.gps_start_time+"-"+opts.duration+".txt")
#      comment = ifo+" snr versus time (left-top), and  cumulative snr (left-bottom) clustered case no veto"
#      comment += "<center><pre>-------------------summary file (category 1) \n"  +summary.read()+"</pre></center>"
#      summary.close()
#    except:
#      logText(logfile, 'SIRE files not found at first inspiral stage', "warning")
#      pass
#  
#    try: 
#     # print opts.datadir+"/"+ thisSearch  +"/"+ifo+"-SIRE_CAT_2_VETO-"+opts.gps_start_time+"-"+opts.duration+".txt"
#      summary = open(opts.datadir+"/"+ thisSearch  +"/"+ifo+"-SIRE_CAT_2_VETO-"+opts.gps_start_time+"-"+opts.duration+".txt")
#     
#      comment += ifo+" snr versus time (right-top), and  cumulative snr (right-bottom) clustered case no veto"
#      comment += "<center><pre>-------------------summary file (category 2) \n"  +summary.read()+"</pre></center>"
#      summary.close()
#      logText(logfile, 'parsing SIRE files at first inspiral stage (cat veto 2)')
#    except:
#      logText(logfile, 'SIRE files not found at first inspiral stage (cat veto 2)', "warning")
#      pass
#
#    page.p("The " +ifo + " inspiral files after the first inspiral stage are clustered (no veto \
#	applied), and the results are shown here below in the next set of plots.")
#    page = add_config_section(page, ifo.lower()+"-inspiral")
#    page = add_figure(page, fnames = [\
#	thisSearch+"/plotinspiral_"+thisSearch+"_"+ifo+"_end_time_vs_snr.png", \
#	thisSearch+"/plotinspiral_"+thisSearch+"_cat_2_"+ifo+"_end_time_vs_snr.png", \
#	thisSearch+"/plotinspiral_"+thisSearch+"_"+ifo+"_snr_cum_hist.png" ,\
#	thisSearch+"/plotinspiral_"+thisSearch+"_cat_2_"+ifo+"_snr_cum_hist.png" \
#	], caption=comment, size="half")
#
#  page.div(id="todo")
#  page.p("todo:distribution in tau0, tau3 ? which tool")
#  page.div.close()# close the todo
#  page.div.close()# close the div openned in add_input_h3 function
#  
#  #---------------------------- the first coincidence stagee, no veto
#  page = heading(page, "First coincidence stage (no veto)")
#  page = add_config_section(page, "thinca")
#  page.add("the number of slide is "+get_numslide("./"+thisSearch)+" time slides")
#  page = add_config_section(page, "thinca-slide")
#  fnames=[]
#  for combi in ("H1H2L1", "H1H2", "H1L1", "H2L1"):
#    page.h4(combi)
#    fnames=[]
#    alt=[]
#    fname_prefix = thisSearch+"/plotthinca_"+thisSearch+"_first_coinc_"
#    fname_suffix ="-"+ str(opts.gps_start_time)  +"-"+str(opts.duration)+".png"
#    alt.append(combi+"hist_slide_trigs")
#    fnames.append(fname_prefix + combi + "_hist_slide_trigs"+fname_suffix) 
#    alt.append(combi+"plot_slide_trigs")
#    fnames.append(fname_prefix + combi + "_plot_slide_trigs"+fname_suffix) 
#    alt.append(combi+"cum_slide_trigs")
#    fnames.append(fname_prefix + combi + "_cum_hist_snr"+fname_suffix) 
#    alt.append(combi+"snr versus snr")
#    if combi=="H1H2":
#      fnames.append(fname_prefix +  "H1_vs_H2_snr"+fname_suffix) 
#    if combi=="H1L1":
#      fnames.append(fname_prefix +  "H1_vs_L1_snr"+fname_suffix) 
#    if combi=="H2L1":
#      fnames.append(fname_prefix +  "H2_vs_L1_snr"+fname_suffix) 
#    

#    comment = "First stage coincidence " +combi + " case. From top to \
#	botton, and left to right, we have a histogram of the number \
#	of triggers in each slide, the distribution of number of triggers\
#	 in each slide and zero lag, the cumulative histogram of number of\
#	 events versus SNR, and the SNR scatter plots for 2 ifos."    
#    try :
#      summary = open(opts.datadir+"/"+ thisSearch  +"/"+combi+"-COIRE-"+opts.gps_start_time+"-"+opts.duration+".txt")
#      comment += "<center><pre>-------------------coire summary file (category 1) \n"  +summary.read()+"</pre></center>"
#      summary.close()
#    except:
#      pass
#    page = add_figure(page, fnames = fnames, caption=comment, size="half", alt=alt)
#  page.div.close()# close the div openned in add_input_h3 function
  


  #---------------------------- the second stage inspiral results, first stage, no veto
#  page = heading(page, "Inspiral jobs (second stage) no veto")
#  page = add_config_section(page, "inspiral")
  # add figures and summary files
#  for ifo in get_ifos():
#    page.h4(ifo)
#    comment = ifo+" snr versus time (left), and  cumulative snr (right) clustered case no veto"
#    page.p("The " +ifo + "inspiral files after the first inspiral stage are clustered (no veto \
#	applied), and the results are shown here below in the next set of plots.")
#    page = add_config_section(page, ifo.lower()+"-inspiral")
#    page = add_figure(page, fnames = [\
#	"playground/plotinspiral_"+thisSearch+"_second_stage_"+ifo+"_end_time_vs_snr.png", \
#	"playground/plotinspiral_"+thisSearch+"_second_stage_"+ifo+"_snr_cum_hist.png", \
#	"playground/plotinspiral_"+thisSearch+"_second_stage_"+ifo+"_snr_over_chi_histogram.png"\
#	], caption=comment, size="third", alt=["end time versus SNR" , "snr cumulative histogram" ,"snr over chi-square histogram"])
#  page.div.close()# close the div openned in add_input_h3 function
  

  #---------------------------- the second coincidence stagee, no veto
 # page = heading(page, "Second coincidence stage (no veto)")
 # page = add_config_section(page, "thinca-2")
 # page.add("the number of slide is "+get_numslide("./"+thisSearch)+" time slides")
 # page = add_config_section(page, "thinca-slide")
 # fnames=[]
 # for combi in ("H1H2L1", "H1H2", "H1L1", "H2L1"):
 #   page.h4(combi)
 #   fnames=[]
 #   alt=[]
 #   fname_prefix = thisSearch +"/plotthinca_"+thisSearch+"_second_coinc_"
 #   fname_suffix ="-"+ str(opts.gps_start_time)  +"-"+str(opts.duration)+".png"
 #   alt.append(combi+"hist_slide_trigs")
 #   fnames.append(fname_prefix + combi + "_hist_slide_trigs"+fname_suffix) 
 #   alt.append(combi+"plot_slide_trigs")
 #   fnames.append(fname_prefix + combi + "_plot_slide_trigs"+fname_suffix) 
 #   alt.append(combi+"cum_slide_trigs")
 #   fnames.append(fname_prefix + combi + "_cum_hist_effective_snr"+fname_suffix) 
 #   alt.append(combi+"snr versus snr")
 #   alt.append(combi+"effective snr versus effective snr")
 #   if combi=="H1H2":
 #     fnames.append(fname_prefix +  "H1_vs_H2_snr"+fname_suffix) 
 #     fnames.append(fname_prefix +  "H1_vs_H2_effective_snr"+fname_suffix) 
 #   if combi=="H1L1":
 #     fnames.append(fname_prefix +  "H1_vs_L1_snr"+fname_suffix) 
 #     fnames.append(fname_prefix +  "H1_vs_L1_effective_snr"+fname_suffix) 
 #   if combi=="H2L1":
 #     fnames.append(fname_prefix +  "H2_vs_L1_snr"+fname_suffix) 
 #     fnames.append(fname_prefix +  "H2_vs_L1_effective_snr"+fname_suffix) 
 #   if combi=="H1H2L1":
 #     fnames.append(fname_prefix +  "H1H2L1_H1_vs_L1_snr"+fname_suffix) 
 #     fnames.append(fname_prefix +  "H1H2L1_H1_vs_L1_effective_snr"+fname_suffix) 
    

 #   comment = "Second stage coincidence " +combi + " case. From top to \
#	botton, and left to right, we have a histogram of the number \
#	of triggers in each slide, the distribution of number of triggers\
#	 in each slide and zero lag, the cumulative histogram of number of\
#	 events versus SNR, and the SNR scatter plots for 2 ifos."    
#    try :
#      logText(logfile, opts.datadir+"/"+ thisSearch  +"/"+combi+"-COIRE_"+combi+"-"+opts.gps_start_time+"-"+opts.duration+".txt")
#      summary =open(opts.datadir+"/"+ thisSearch  +"/"+combi+"-COIRE_"+combi+"-"+opts.gps_start_time+"-"+opts.duration+".txt")
#      comment += "<center><pre>-------------------coire second stage, summary file (category 1) \n"  +summary.read()+"</pre></center>"
#      summary.close()
#    except:
#      page.div(id="todo")
#      page.p("summary file from sire not found")
#      page.div.close()
#      pass
#    page = add_figure(page, fnames = fnames, caption=comment, size="half", alt=alt)
#  page.div(id="todo")
#  page.p("todo:add --html-output option in plotthinca, have the cache file forthe first coincidence stage available ")
#  page.div.close()# close the todo
#
#  page.h4("All trigges together")
#  comment = "Cumulative histogram versus statistic (effective SNR) in top-left panel, SNR versus time (top right). slide histograms (bottom)"
#  page.p("The " +ifo + "inspiral files after the first inspiral stage are clustered (no veto \
#	applied), and the results are shown here below in the next set of plots.")
#  page = add_figure(page, fnames = [\
#	thisSearch+"/plotthinca_"+thisSearch+"_second_coinc_cum_hist_effective_snr-847555570-2419200.png",\
#	thisSearch+"/plotthinca_"+thisSearch+"_second_coinc_effective_snr_vs_time-847555570-2419200.png",\
#	thisSearch+"/plotthinca_"+thisSearch+"_second_coinc_hist_slide_trigs-847555570-2419200.png",\
#	thisSearch+"/plotthinca_"+thisSearch+"_second_coinc_plot_slide_trigs-847555570-2419200.png"\
#	], caption=comment, size="half")
#
#  page.div.close()# close the div openned in add_input_h3 function
#
#
#  
#  #---------------------------- loudest events
#  page = heading(page, "Summary")
#  combis = ['H1H2L1', 'H1H2', 'H1L1', 'H2L1']
#  vetos = ['', '_CAT_2_VETO-', '_CAT_3_VETO-', '_CAT_4_VETO-']
# # try:
#  for combi in combis:  
#    page.table();
#    page.tr(e.td(combi, width=5))
#    page.tr()
#    page.td('');
#    page.td(vetos)
#    page.tr.close()
#    page.tr()
#    page.td('zero lag')    
#    for veto in vetos:
#      text = opts.datadir  + '/'+thisSearch+'/'+ combi+'-COIRE_'+combi+veto+'*txt'
#      file = glob.glob(text)
#      try:
#        command = 'grep clusters  '+file[0] + '|awk \'{print $NF}\''
#        data, status = make_external_call(command, opts.debug, opts.debug, True)
#        page.td(data)
#      except:
#        page.td('not found')
#        logText(logfile, 'WARNING' + text +' not found', "warning") 
#    page.tr.close()
#    page.tr()
#    page.td('time slides')    
#    for veto in vetos:
#      try:
#        file = glob.glob(opts.datadir  + '/'+thisSearch+'/'+ combi+'-COIRE_SLIDE_'+combi+veto+'*txt')
#        command = 'grep clusters  '+file[0] + '|awk \'{print $NF}\''
#        data, status = make_external_call(command, opts.debug, opts.debug, True)
#        page.td(data)
#      except:
#        page.td('file not found')
#        logText(logfile, 'WARNING' + text +' not found', "warning") 
#
#    page.tr.close()
#    page.table.close()
#  except:
#      page.div(id="todo")
#      page.p("Summary file from sire not found")
#      page.div.close()
#      pass
 # page.div.close()# close the h3


  #---------------------------- loudest events

  # the end
  page.div.close()
  page.div.close()
  page.div.close()
  return page

# ***************************************************************************
# ***************************************************************************
def write_injection(page, opts):
  webdir = opts.webdir
  ini = opts.config_file
  ifos = get_ifos()
  mkdir(opts.physdir+"/injections")
  page.p(e.br())

  page.div(class_="encadre")
  page.div(class_="contenu")
  page = write_title(page, html_sections['injection'], "baudelaire")
  page.div(id="div_baudelaire", style='display:none')
  page.add("This section summarizes the analysis of the injection runs.<br/>")

  #---------------------------- found and missed first inspiral stage
  page = heading(page, "Found versus missed plots (first inspiral stage)")
  for injtag in ['bbh_inj', 'bns_inj', 'bhns_inj']:
    page.h4(injtag)
    for ifo in get_ifos():  
      lists = glob.glob(opts.datadir+"/"+ injtag +"/"+ifo+"-SIRE_INJECTIONS_*FOUND-*.txt")
      try:
        summary = open(lists[0])
        this = summary.read()
      except:
        this = "No summary file found"
        lists = ["No summary file found"]
        pass
      caption = ifo +" Found and missed injections at the first stage of the inspiral analysis ("+injtag +"tag)."      
      caption = caption + "<pre> -----------------------"+lists[0]+"\n" +this+"</pre>"
      page = add_figure(page, fnames=  [\
        "injections/"+ifo+"_first_insp_"+injtag+"_dist-mchirp.png", \
        "injections/"+ifo+"_first_insp_"+injtag+"_time-dist.png"], \
        caption=caption, size="half",alt=["", ""])
  page.div.close()
  
  #---------------------------- found and missed first coinc stage
  page = heading(page, "Found versus missed plots (first coincidence stage)")
  for injtag in ['bbh_inj', 'bns_inj', 'bhns_inj']:
    page.h4(injtag)
    for combi in get_ifo_coinc():  
      if len(combi)>2:
#        blabla(logfile, opts.datadir+"/"+ injtag +"/"+combi+"-COIRE_INJECTIONS_*_FOUND_"+combi+"-*.txt")
        lists = glob.glob(opts.datadir+"/"+ injtag +"/"+combi+"-COIRE_INJECTIONS_*_FOUND_"+combi+"-*.txt")
        try:
          summary = open(lists[0])
          this = summary.read()
        except:
          this = "No summary file found"
          lists = ["No summary file found"]
          pass
        caption = combi +" Found and missed injections at the first coincidence stage ("+injtag +"tag)."      
        caption = caption + "<pre> -----------------------"+lists[0]+"\n" +this+"</pre>"
        page = add_figure(page, fnames=  [\
          "injections/"+combi+"_first_coinc_bbh_inj_missed_found_chirp_dist_vs_mchirp.png", \
          "injections/"+combi+"_first_coinc_bbh_inj_missed_found_vs_end_time.png", \
    	  "injections/"+combi+"_first_coinc_bbh_inj_missed_found_vs_mchirp.png"], \
          caption=caption, size="half",alt=["", "", ""])
  page.div.close()
  #---------------------------- close by missed
  page = heading(page,  "close by missed")
  page = add_figure(page, fnames = ["H1H2L1-found_missed.png"], caption="", size="half",alt=[""])
  page = add_figure(page, fnames = ["H1H2-found_missed.png"], caption="", size="half",alt=[""])
  page = add_figure(page, fnames = ["H1L1-found_missed.png"], caption="", size="half",alt=[""])
  page.div.close()
  
  #---------------------------- parameter efficiency
  page = heading(page, "Parameter efficiency, second coincidence stage, no veto")
  tags = ['eff_dist_frac', 'end_time','end_time_vs_snr','end_time_vs_eff_dist', 'eta','eta_vs_snr', 'mchirp_vs_eff_dist','mchirp']
  tags = ['eff_dist_frac', 'end_time_vs_snr','eta_vs_snr', 'mchirp_vs_snr']
  for tag in tags:
    page.h4(tag)
    for combi in ['H1H2L1','H1H2','H1L1','H2L1']:
    #for combi in ['H1H2L1']:
      page.h5(combi)
      fnames=[]
      alt=[]
      for injection_tag in ['bns_inj', 'bhns_inj', 'bbh_inj']:
        this = 'plotinspinj_'+injection_tag +'_'+combi+'_'+tag+'_accuracy.png'
        alt.append(this)
        fnames.append('injections/plotinspinj/'+this)
      page = add_figure(page, fnames =fnames, caption="test", size="third",alt=alt)
  page.div.close()

  page = heading(page, "Efficiency, second coincidence stage, veto and no veto")
  page.div.close()
  
  # the end
  page.div.close()
  page.div.close()
  page.div.close()

  return page


def write_about(page, opt):
  webdir = opts.webdir
  page.p(e.br())
  page.div(class_="encadre")
  page.div(class_="contenu")
  page = write_title(page, html_sections['about'], "balzac")
  page.div(id="div_balzac", style='display:none')
  page.p("This page was automatically generated with write_ihope_page using the following ini file")
  tmp  = open(opts.config)
  page.pre(tmp.read())
  tmp.close() 

  page.p("and the following command line arguments:")
  text=""
  for arg in sys.argv[:]:
    text = text +  arg +" "
  page.pre( text )
  
  page.p(__Id__[4:len(__Id__)-5])
  page.div.close()
  page.div.close()
  page.div.close()
  return page


# ***************************************************************************
# ***************************************************************************
def get_numslide(run):
  logText(logfile, opts.datadir + "/" + run + "/inspiral_hipe_" + run )
  return str(29)
  

# ***************************************************************************
# ***************************************************************************
def add_config_section(page, section):
  ini  = hipecp.items(section)
  page.add("<pre>")
  page.add("------------------------------------" +section)
  for i in  ini :
    page.add(i[0]+' = '+i[1])
  page.add("</pre>")
  return page

# ***************************************************************************
# ***************************************************************************
def add_input_h3(page, title):
  #increment block number
  global count_block 
  input=str(count_block)
  count_block=count_block+1

  page.add("<h3> "+title)
  page.input(id="input_"+input, type="checkbox", onclick="toggleVisible('"+input+"')")
  page.add("See details here</h3>")
  page.div(id="div_"+input , style='display:none')  
  return page

def heading(page, title="None", label="Switch details on/off", heading="h3"):
  #increment block number
  global count_block 
  input=str(count_block)
  count_block=count_block+1

  page.add("<"+heading+">"+ title)
#  page.input(id="input_"+input, type="checkbox", onclick="toggleVisible('"+input+"')", size="4")
  text = label 

  page.input(id="input_"+input, type="button", onclick="toggleVisible('"+input+"')", value=text ,class_="toggle")
  page.add("</"+heading+">")
  page.div(id="div_"+input , style='display:none') 


  logText(logfile, '      Enters sub-section: '+title )
  return page 

# ***************************************************************************
# ***************************************************************************
def add_caption(page, caption):
  global fig_num
  page.p("Figure "+str(fig_num) + ": "+caption)
  fig_num = fig_num + 1
  return page

def add_figure(page,fnames="test", caption="add a caption", size="full", alt="no figure found"):
  global fig_num
  dir = opts.webdir
  page.add("<!-- insert a figure -->\n<div class=\"figure\">")
  this_count = 0
  
  for fnam in fnames:
    source=dir+"/"+fnam
    
    title = " title=\""+alt[this_count]+"\""
    page.add("\t<a href=\"" + source.replace("_small", "")+"\"" +  title + ">\n" )
    try:
      page.img(class_=size ,src=fnam, alt=alt[this_count] )
    except:
      page.img(class_=size ,src=fnam )
      
    page.add("</a>")
    this_count = this_count +1
  page.add("\t<p class=\"figure\">Figure " + str(fig_num) + ":  " + caption + "</p>\n</div>\n\n")
  fig_num +=1

  #check that the file exsits. If not, returns an error message
  return page

# ***************************************************************************
# ***************************************************************************
####
def create_venn(data, tag):
  """
  Create a venn diagramm for the 2 or 3 ifos case
  data has to be H1H2 H1 H2L1 H2 H1L1 L1 H1H2L1 array 
  """
  try:
    mscript = open("temp.m", "w")
    mscript.write("data = [")
    mscript.write( data +"];" ) 
    mscript.write(" vennX(data\'/3600/24,0.01);")
    mscript.write(" k=colormap(jet); k = [1 1 1; k];colormap(k); saveas(gcf,\'venn_"+tag+".png\')")
    mscript.close()
    command=("matlab -nodisplay -nodesktop -nojvm -nosplash   < temp.m > /tmp/test ;  rm -f temp.m; mv venn_"+tag+".png "+opts.physdir+"/"+tag+"/")
    if not opts.debug:
      make_external_call(command, opts.debug, opts.debug, True)   
  except:
    logText(logfile, """WARNING   The matlab command to create the venn diagram failed. 
                Check that matlab is properly set,and vennX.m is available
                (see matapps/src/searches/inspiral/matspiral/utilities")  
          """, "error")

# ***************************************************************************
# ***************************************************************************
def sortDict(d):
    items = d.items()
    items.sort()
    return [value for key, value in items]

# ***************************************************************************
# ***************************************************************************
def get_coincident_segments(tag):
  """
  return duration of coincident segments
  """
  ifos = get_ifos()
  thisdata = {}
  thisdata['segments'] = ['duration(days)' ,'duration (s)']
  duration = str((int(opts.gps_end_time) - int(opts.gps_start_time)))
  output={}
  ifo_coincs = get_ifo_coinc()  

  try:
    for coinc in ifo_coincs:
      if tag=="playground":
        command = "awk \'{sum=sum+$4} END {print sum}\' "+opts.datadir+ tag +"/"+coinc+"_play_segs_analyzed.txt"
      elif tag=="analysis":
        command = "awk \'{sum=sum+$4} END {print sum}\' "+opts.datadir+ tag +"/"+coinc+"_segs_analyzed.txt"

      output[coinc], status = make_external_call(command, opts.debug,opts.debug, True)
    logText(logfile, '...Get the analyzed segments duration...')
  except:
    logText(logfile , 'problem(s) while parsnig the coincident segments', "error")

  return output

  


# ***************************************************************************
# ***************************************************************************
def get_segments_tag(tag):
  """
  reads segment files and return table of duration for each ifo
  """
  datadir = opts.datadir
  ifos = get_ifos()
  thisdata = {}
  thisdata['segments'] = ['duration(days)' ,'duration (s)']
  duration = str((int(opts.gps_end_time) - int(opts.gps_start_time)))
  duration.split('.')
  this_tag  = '-' + tag + '-' + opts.gps_start_time + '-' + str(duration) + '.txt'
  if tag=="RDS_C03_L2":
    this_tag  = '_' + tag + '-' + opts.gps_start_time + '-' + str(duration) + '.txt'
  command = 'awk \'{print NF}\' ' +  datadir +'/segments/' + ifos[1] + this_tag
  try:
    ncols, status = make_external_call(command, False, opts.debug, True)
    ncols = ncols[len(ncols)-2]
    if float(ncols)==4:  
      for ifo in ifos:
        command = 'awk \'{sum=sum+$4} END {print sum/3600/24}\' ' +  datadir +'/segments/' + ifo + this_tag
        output_days, status = make_external_call(command, opts.debug, opts.debug, True)
        command = 'awk \'{sum=sum+$4} END {print sum}\' ' +  datadir +'/segments/' + ifo + this_tag
        output_seconds, status = make_external_call(command, opts.debug, opts.debug, True)
        thisdata[ifo] = [ output_days, output_seconds]
    elif float(ncols)==2:
      for ifo in ifos:
        command = 'awk \'{sum=sum+$2-$1} END {print sum/3600/24}\' ' +  datadir +'/segments/' + ifo + this_tag
        output_days, status = make_external_call(command,opts.debug,opts.debug, True)
        command = 'awk \'{sum=sum+$2-$1} END {print sum}\' ' +  datadir +'/segments/' + ifo + this_tag
        output_seconds, status = make_external_call(command, opts.debug, opts.debug, True)
        thisdata[ifo] = [ output_days, output_seconds]
  except:
    logText(logfile,  'Error while parsing the segment duration files', "error")
 # else:
#    logText(logfile,  '...Get the segment duration for '+tag+'...')

  return thisdata


# ***************************************************************************
# ***************************************************************************
def get_ifo_coinc():
  """
  return list of valid coincidences with respect to the ini file
  """
  # get the ifo requested
  ifos = get_ifos()
  # get number of ifos to look at
  numifos = []  
  for option in ["one-ifo","two-ifo","three-ifo","four-ifo", "five-ifo"]:
    if hipecp.has_option("ifo-details",option): 
	tmp = option.split("-")
 	numifos.append(tmp[0])
  # now let us fill the different possible ifo combinaisons
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
# ***************************************************************************
def get_ifos():
  """
  read the ifos used in the ifo-details section of the ihope ini file
  """
  ifos=[] 
  for option in ["g1-data","h1-data","h2-data","l1-data", "v1-data"]:
    if hipecp.has_option("ifo-details",option): ifos.append(option[0:2].upper() )

  return ifos


# ***************************************************************************
# ***************************************************************************
def get_version(executable): 
  output=[]
  try:
    if executable=="trigbank":
      pathname = hipecp.get("condor", "trigbank") 
    else:
      pathname = hipecp.get("condor", executable) 
    s =  os.path.split(pathname)
    if len(s[0])==0:
      cmd = opts.ihope_directory + s[1] + '  --version'
    else:
      cmd = s[1] + ' --version '
    output,status = make_external_call(cmd, opts.debug, opts.debug, True)
  except:  
    output= '(not found)' 
    pass

  if output.find('CVS Tag')>0:
    output = output.split()
    output = ' has tag version '+output[len(output)-2]
  elif output.find('CVS Version')>0:
    output = output.split()
    output = ' no tag found but the version is '+output[len(output)-7] + output[len(output)-6] 
  else:
    output = ' not tag found, no version found...check me'
  return output

# ***************************************************************************
# ***************************************************************************
#### function to copy the segment files in the appropriate directory
def copy_segments():


  # parsing the ini file, find the cat file and thenread the ./segments directory
  valid="done"
  try:
    stdout= '...... 1 - the category files into /catlists'
    cats = hipecp.items("segments")
    location = opts.physdir+'/catlists/'
    mkdir(location)
    for file in cats:
      if file[0].find('analyze')<=0:
        command = 'cp '+opts.datadir +'/segments/'+  file[1] +' ' + location
        dummy,status = make_external_call(command, False, opts.debug, True)
        if status>0:
          stdout += " WARNING:: could not copy a category segment list"
          valid="warning"
          break
  except:
    stdout +=" ERROR:: problem while copying category segment list. passing..."    
    logText(logfile,  stdout, "warning") # we do not want to quit here. just a warning
    pass
  else: 
    logText(logfile,  stdout, valid)

 #  the segment files
  valid="done"
  try :
    stdout ='...... 2 - the segments files into /segments'
    for this in get_ifos():
        seg = this +'-SELECTED_SEGS-'+opts.gps_start_time+'-'+opts.duration+'.txt'
        try:
          mkdir(opts.physdir+'/segments')
        except:raise
  
        command = 'cp '+opts.datadir+'/segments/'+  seg + ' '+opts.physdir+'/segments/'
        dummy,status = make_external_call(command, True, opts.debug, True)
        if status>0:
          stdout += " WARNING:: could not copy a segment list"
          valid="warning"
          break
  except:
    stdout +=" WARNING:: problem while copying segment list. passing..."    
    logText(logfile,  stdout, "warning") # we do not want to quit here. just a warning
    pass
  else: 
    logText(logfile,  stdout, valid)

 #  the selected segment files
  for thisSearch in ['playground', 'analysis']:
    try :
      stdout ='...... 3 - the selected segments files into /'+thisSearch
      for this in get_ifo_coinc():
          if thisSearch=='playground':
            seg = '/'+thisSearch+'/'+this +'_play_segs_analyzed.txt'
          elif thisSearch=='full_data':
            seg = '/'+thisSearch+'/'+this +'_segs_analyzed.txt'
          mkdir(opts.physdir+'/'+thisSearch)
          command = 'cp '+opts.datadir + seg + ' '+opts.physdir+'/' +thisSearch +'/'
          dummy,status = make_external_call(command, opts.debug, opts.debug, True)
          if status>0:
            stdout += " WARNING:: could not copy a selected segment (playground)"
            break
    except:
      stdout +=" WARNING:: problem while copying a selected segment (playground). passing..."    
      logText(logfile,  stdout, "warning")
      pass
    else: 
      logText(logfile,  stdout)

def fom(page, opts, tag, caption, directory="playground_summary_plots"):
  """
  """

  dir = opts.datadir+'/'+directory+'/' 
  # create the div that will contain the figure 
  page.div(class_="figure")
  # search for the cache file that contain the tag in
  # the directory considered (given as fifth argument)
  thisglob = dir + tag +'cache'
  filenames =  glob.glob(thisglob)
  # and parse this cache files searching for images
  for eachcache in filenames:
    this = open(eachcache, "r")
    fnameList = []
    for filename in this:
       #For each images, there is also a thumb file
      if filename.find('png')>=0:
        fnameList.append(filename.replace(".png", "_thumb.png"))
        command = 'cp ' + dir + filename+' ' +opts.physdir +'/Images/'       
        make_external_call(command.replace("\n", " "), False, opts.debug, True)
        command = 'cp ' + dir + filename.replace(".png","_thumb.png")+' ' +opts.physdir  + '/Images/' 
        make_external_call(command.replace("\n", " "), False, opts.debug, True)
    this.close()
    
    this = open(eachcache, "r")
    for filename in this:
      if filename.find('html')>=0:
        command = 'cp ' + dir + filename+' ' +opts.physdir        
        make_external_call(command.replace("\n", " "), False, opts.debug, True)
        msg =" <a href=\""+filename.replace("/pictures", "")+\
            "\">click here to get all pictures and arguments used to generate the plots</a> "
    this.close()
   
    page = add_figure(page, fnames=fnameList, caption=(caption+msg), size="third", alt=fnameList)

    this.close()
  page.div.close()

  return page
# ***************************************************************************
# ***************************************************************************
def set_style():
  """
  """
  tmp = []
  tmp.append(configcp.get("main", "style"))
  make_external_call('cp ' + configcp.get("main", "style")+ " " +opts.physdir, opts.debug, opts.debug, True )

  # 10 style file should be allright ;-)
  for i in range(1,10,1):
    try:
      tmp.append (configcp.get("main", "style"+str(i)))    
      make_external_call('cp ' + configcp.get("main", "style"+str(i))+ " " +opts.physdir, opts.debug, opts.debug, True)
    except:
      pass
  style = tmp

  count = 0
  for this in style:
    tmp =  this.split('/')
    style[count] = tmp[len(tmp)-1]
    count = count +1

  return style

# ----------------------------------------------------------------------------
def parse_arguments():
  usage =  """ %prog [options]
  Program to write webpage from upperlimit.py
  """

  parser = OptionParser( usage = usage, version = "%prog CVS "+__Id__ )

  parser.add_option("-C","--config-file",action="store",type="string",\
      metavar=" INI File",\
      help="ini file with information about run directories" )
  parser.add_option("-A","--skip-analysis",action="store_false",\
      default=True,dest="analysis", metavar="DOANALYSIS",\
      help="" )
  parser.add_option("-T","--skip-tuning",action="store_false",\
      default=True,dest="tuning", metavar="DOTUNING",\
      help="" )
  parser.add_option("-U","--skip-upperlimit",action="store_false",\
      default=True,dest="upperlimit", metavar="DOUPPERLIMIT",\
      help="" )
  parser.add_option("-I","--skip-injection",action="store_false",\
      default=True,dest="injection", metavar="DOUPPERLIMIT",\
      help="" )
  parser.add_option("-S","--skip-summary",action="store_false",\
      default=True,metavar="DOSUMMARY",\
      help="" )
  parser.add_option("-P","--skip-playground",action="store_false",\
      default=True,dest="playground", metavar="DOPLAYGROUND",\
      help="" )
  parser.add_option("-D","--debug",action="store_true",\
       default=False,dest="debug",metavar="DODEBUG",\
       help="" )
  parser.add_option("-V","--verbose",action="store_true",\
      default=False, dest="verbose",metavar="VERBOSE",\
      help="" )

  (opts,args) = parser.parse_args()

  if opts.config_file is None:
    raise ValueError,\
"""
------------------------------------------------------------------------------
the arguments --config-file must be used and followed by an ini file, an
example of which is  :

[main]
gps-start-time  = 847555570
gps-end-time    = 849974770
title           = "Low mass CBC analysis"
ihope-ini-file  = ihope.ini
ihope-directory = /archive/home/cokelaer/S5/Month1/full_analysis/
url             = ldas-jobs.ligo.caltech.edu
username        = cokelaer
output          = index.shtml
style           = /archive/home/cokelaer/style.css
"""

  return opts,args

# ----------------------------------------------------------------------------

#try:
#  from glue import markup
#  from markup import oneliner as e
#except: 
#  raise ImportError("Require markup.py to generate the html page")
#
opts,args = parse_arguments()


#############################################################################
#  MAIN PART                                                                #
#############################################################################
fig_num = 1
count_block = 0
config   =  opts.config_file
opts.config = config # save the name of the ini file, why ?
configcp = ConfigParser.ConfigParser()
configcp.read(config)

#----------------------------------------
# First, we open an xml file, for the log
print >>sys.stdout,"Openning the log file"
try :
  logfile = open(__name__+".xml", "w")
  logfile.write("""<?xml version="1.0" encoding="ISO-8859-1"?>
<?xml-stylesheet type="text/xsl" href="write_ihope_page.xsl"?>
<log>
""")
except:
  print >>sys.stderr, "Could not create the log file...quitting"
  raise 

#---------------------------------------
print >>sys.stdout,"Parsing the ini file: " + opts.config
#parsing the ini file
try:
  opts.config_file 	= configcp.get("main", "ihope-ini-file")
  opts.gps_start_time 	= configcp.get("main", "gps-start-time")
  opts.gps_end_time 	= configcp.get("main", "gps-end-time")
  opts.ihope_directory 	= configcp.get("main", "ihope-directory")
  opts.title	 	= configcp.get("main", "title")
  opts.url	 	= configcp.get("main", "url")
  opts.username	 	= configcp.get("main", "username")
  opts.output	 	= configcp.get("main", "output")
except:
  print >>sys.stderr, "ERROR::The ini file does not have the proper field in the [main] section" 
  print >> sys.stderr,"       Consider adding those fields if missing: ihope-ini-file, gps-start-time,gps-end-time, ihope-directory, title,url, username, output"
  raise
  
#------------------------------------
#sub-products of the ini file parsing
print >>sys.stdout,"Creating sub-product variables"
opts.gpsdir =  '/'+str(opts.gps_start_time)+'-'+str(opts.gps_end_time)+'/'
opts.duration = str(int(opts.gps_end_time) - int(opts.gps_start_time))
opts.webdir = 'http://' + opts.url + '/~' + opts.username  + opts.gpsdir
opts.datadir = opts.ihope_directory + opts.gpsdir
opts.physdir = '/archive/home/'+opts.username+'/public_html/'+ opts.gpsdir
  
#------------------------------------------------------------------------------------
#read the different css style available and copy the style files in the web directory.
print >>sys.stdout,"Getting the CSS style"
try:
  opts.style = set_style()
  for this in opts.style:
    print >> sys.stdout, "   Found " + this
except:
  logText(logfile, "Problem with the style file. (either not copied or not found)", "warning")
  raise


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
# here is the directory we want to look at
msg = "Entering " + opts.datadir
print >> sys.stdout, msg
if not  os.path.isdir(opts.datadir):
  raise  "%s is not a valid directory. Check your gps time." % opts.datadir
# which physical name is 
if not os.path.isdir(opts.physdir):
  "%s is not a valid directory. Trying to create it." % opts.physdir
  try:
    os.mkdir(opts.physdir)
    logText(logfile, "##Warning:: directory "+opts.physdir+" created","warning" )
    os.mkdir(opts.physdir+'/Images')
    logText(logfile, "##Warning:: directory "+opts.physdir+"/Images created","warning" )
  except: pass
else: 
  print opts.physdir +'already exists'
if not os.path.isdir(opts.physdir+'Images'):
  "%s is not a valid directory. Trying to create it." % opts.physdir
  try:
    os.mkdir(opts.physdir+'/Images')
    logText(logfile, "##Warning:: directory "+opts.physdir+"/Images created","warning" )
  except: pass
else: 
  print opts.physdir+'/Images' +'already exists'

#-----------------------------------------
# now we can parse the ihope.ini file
msg =   "...Parsing the ihope ini file (" + opts.config_file+")"
logText(logfile, msg)
hipe   = opts.ihope_directory+'/'+opts.config_file
hipecp = ConfigParser.ConfigParser()
hipecp.read(hipe)
make_external_call( 'cp '+opts.config_file + ' ' + opts.physdir, False, opts.debug, True)

#-----------------------------------------
# now we copy the segments to the web directory
msg =   "Copying segments"
logText(logfile, msg)
copy_segments()


#-----------------------------------------
###### create the section labels  
html_sections={}
html_order = ['toc', 'general', 'summary', 'playground', 'injection', \
    'tuning', 'analysis', 'upperlimit', 'logfile', 'about']

html_sections['toc'] = "Table of Contents"
html_sections['general'] = "General Information"
html_sections['summary'] = "Data Information"
if opts.playground is True: html_sections['playground']	= "Playground"
if opts.injection is True: html_sections['injection'] = "Injection"
if opts.tuning is True: html_sections['tuning']	= "Tuning"
if opts.analysis is True: html_sections['analysis'] = "Full Data"
if opts.upperlimit is True: html_sections['upperlimit']	= "Upper Limit"
html_sections['logfile'] = "Log File"
html_sections['about'] = "About"


title = opts.title+ ": from "+str(opts.gps_start_time)+" to "+str(opts.gps_end_time) 
script = {}
script['toggle.js'] = 'javascript'
create_toggle()
# Finally, we create the html document 
msg =   "Creating HTML document"
logText(logfile, msg)

page = markup.page(mode="strict_html")
page._escape = False
doctype="""<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">"""
doctype+="""\n<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">"""

page.init(title=title, css=opts.style, script=script , doctype=doctype)
page.h1(opts.title +" (" + opts.gps_start_time +"-" +opts.gps_end_time+")")

##### we write the toc, and general section
for each_section in html_order:
  try:
    logText(logfile, html_sections[each_section], "section")
    page = write_results(page, opts, each_section)
  except:
    logText(logfile, each_section, "section")
    msg = "skipped  "+each_section + " section"
    logText(logfile, msg, "warning")
    print >>sys.stdout, msg
    pass

 

##### the endi


page.add("<table><tr><td>")
page.add("<a href=\"http://validator.w3.org/check?uri=referer\">")
page.add("<img src=\"http://validator.w3.org/images/valid_icons/valid-xhtml10\" alt=\"Valid XHTML 1.0!\" height=\"31\" width=\"88\"/></a> ")
page.add("<a href=\"http://jigsaw.w3.org/css-validator/\">")
page.add("<img src=\"http://www.w3.org/Icons/valid-css\" alt=\"Valid CSS!\" height=\"31\" width=\"88\" /></a> ")
page.add("</td></tr></table>")

html_file.write(page(False))
html_file.close()
logfile.write("</log>")
logfile.close()

print '---------------------FINISHED ---------------------'
print '--- HTML file created. '
print '--- Copying html documents in ' +opts.physdir
make_external_call('mv  '+opts.output +' ' + opts.physdir, opts.debug, opts.debug, True)
make_external_call( 'mv toggle.js '+ opts.physdir, opts.debug, opts.debug,  True)

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



