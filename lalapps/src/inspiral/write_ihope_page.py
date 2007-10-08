#!/usr/bin/python

# $Id$
__author__ = "Thomas Cokelaer <thomas.cokelaer@astro.cf.ac.uk>"
__version__ = "$Revision$"
__date__ = "$Date$"
__Id__ = "$Id$"


import sys
import copy
import os
from optparse import *
import glob
import markup
from markup import oneliner as e
import ConfigParser


# ***************************************************************************
# ***************************************************************************
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
      print this_output

  stdin.close()
  out.close()
  err.close()  
  return this_output

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


# ***************************************************************************
# ***************************************************************************
def write_toc(page , opts):
  page.add("<!-- beginning of toc section-->")
  dir = opts.webdir 
  items = html_sections 
  page.div(id="encadre")
  page.div(id="contenu")
  page.a(e.h2("Table of Contents"), name="toc")
  page.ol()
  for this_item in items:
    page.li(e.a(this_item, href=dir+"/"+opts.output+'#'+this_item))
  page.ol.close()

  
  page.div.close()
  page.div.close()
  return page

# ***************************************************************************
# ***************************************************************************
def write_general(page,opts):
  webdir = opts.webdir
  ini = opts.config_file

  page.br()
  page.add("<!-- beginning of the general section -->")
  page.div(id="encadre")
  page.div(id="contenu")
  page = write_title(page, html_sections[0], "rimbaud")
  page.div(id="div_rimbaud", style='display:block')
  text=  'This page summarizes the analysis of the data segment from GPS time %s up to %s' % (opts.gps_start_time,  opts.gps_end_time)
  page.add(text)

  page.h3("The search used the following resources:")
  page.br()
  page.ol()
  page.li("LAL/LALApps has "+get_version("inspiral"))

  executables = ("inspiral", "tmpltbank", "sire", "thinca", "trigbank", "coire" , "inspinj")
  page.ul()
  for exe in executables:
    text = "<b>lalapps_"+exe+ "</b>:    "+get_version(exe)
    page.li(text)   
  page.ul.close()

  page.li("Segment information:")
  page.ul()
  for this in get_ifos():
      seg = this +'-SELECTED_SEGS.txt'
      this_ref = webdir + '/segments/'+seg
      page.li(e.a(seg,href=this_ref) )
  page.ul.close()
  text=("The configuration is contained in the file <a href=\"" + \
      webdir + "/" + ini + "\">" + ini + "</a>")
  page.li(text)
  text = "A list of category files stored in this directory ("  \
	+ "<a href=\""+ webdir + "/catlists/\"> catlists</a>"  +\
	"), and listed in the ini file."
  page.li(text)
  page.ul.close()
  page.ol.close()

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
  page.add("<a name=\""+text+"\">")
  page.h2()
  input = "input_"+tag
  page.input(id=input, type="checkbox", checked="on", onclick="toggleVisible('"+tag+"')")
  page.add(text)
  page.h2.close() 
  page.add("</a>")
  page.div(e.a("return to top", href="#toc" ), id="back")
#  page.div(e.a("close all", href="#toc" ), id="back")

  return page


# ***************************************************************************
# ***************************************************************************
def write_data_summary(page,opts):
  webdir = opts.webdir
  datadir = opts.datadir
  ini = opts.config_file
  ifos = get_ifos()

  page.br()
  page.add("<!-- beginning data summary section-->")
  # starts with 2 divs that needs to be closed
  page.div(id="encadre") #(1)
  page.div(id="contenu") #(2)
  # title with a div to toggle on/off the content 
  page = write_title(page, html_sections[1], "verlaine")
  page.div(id="div_verlaine", style='display:block') #(3)
  # start an heading that opens a div (4)
  page.p("This section summarizes the data segments information.")
  page = heading(page, "Selected segments", "Switch details on/off", "h3") 
  page.add("The segments files above were created with no data quality flags set")
  page.br()
  page.add("The times analyzed according to hipe are (SELECTED_SEGS):")
  page.br()
  page.table()
  segs = get_segments()
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
        command = 'awk \'{print $1}\' ' +  datadir +'/../' + ifo + "cat"+str(i)+".txt"
        flags = make_external_call(command, 0, 0)
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
      print "Problems parsing category veto segment list "+str(i)
      page.add("unknown or not set. ")
      
  page = heading(page, "Science segments ") #(4)
  page.add("The science segments (category 1) according to hipe are (SCIENCE): <br/>")
  page.table()
  segs = get_segments_tag("SCIENCE_SEGMENTS")
  keys = ("segments","H1","H2","L1","G1","V1","T1") 
  for key in keys:
    if segs.has_key(key):
      page.tr()
      page.td(key)
      page.td(segs[key])
      page.tr.close()
  page.table.close()
  page.div.close() #(3)
  
  page  = heading(page, "RDS_C03_L2 segments") #(4)
  page.add("The RDS_C)3_L2 segments according to hipe are (SCIENCE): <br/>")
  page.table()
  segs = get_segments_tag("RDS_C03_L2")
  keys = ("segments","H1","H2","L1","G1","V1","T1") 
  for key in keys:
    if segs.has_key(key):
      page.tr()
      page.td(key)
      page.td(segs[key])
      page.tr.close()
  page.table.close()
  page.div.close() # (3)


  page.div.close() #(2)
  page.div.close() #(1)
  page.div.close() #(0)
  
  return page


# ***************************************************************************
# ***************************************************************************
def add_figure(fname,dir,fnames,size,caption):
  global fig_num
  fig_num +=1
  fname.write("<div class=\"figure\"> \n<p>")
  for fnam in fnames:
    fname.write("<a href=\"" + "dir" + "/" + "fnam" +"\">" + \
      "<img class=\"" + size + "\" src=\"" + "dir" + "/" + fnam + "\"></a>")
  fname.write("<p>Figure " + str(fig_num) + ":  " + caption + "\n</div>\n\n")


# ***************************************************************************
# ***************************************************************************
def write_injections_efficiencies(fname,dir):
  fname.write("<a name=\"eff\"><h2>Efficiency Plots</h2></a>\n\n")

  for ifos in opts.ifos:
    fname.write("<h3>" + ifos + "</h3>")

    if "L" in ifos:
      plot_name = "plotnumgalaxies/" + ifos + "/" + ifos + "_" + \
          opts.x_value + "_" + opts.y_value
      print plot_name
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
  page.br()
  page.div(id="encadre")
  page.div(id="contenu")
  page = write_title(page, html_sections[6], "zola")
  page.div(id="div_zola", style='display:block')
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
  page.div(id="todo")
  page.p( "to do : many things")
  page.div.close()# close the todo
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
  page.br()
  page.div(id="encadre")
  page.div(id="contenu")
  page = write_title(page, html_sections[4], "ronsard")
  page.div(id="div_ronsard", style='display:block')
  page.add("This section summarizes the tuning .<br/>")
  
  #---------------------------- loudest events
  page = add_input_h3(page, "injections and time slides")
  page.div(id="todo")
  page.p( "add list of loudest events. Not mix this with follow up")
  page.div.close()# close the todo
  page.div.close()# close the h3

  # the end
  page.div.close()
  page.div.close()
  page.div.close()
  return page


# ***************************************************************************
# ***************************************************************************
def write_fulldata(page,opts):
  print "....Processing  full data section"

  this_tag = "analysis"
  webdir = opts.webdir
  ini = opts.config_file
  ifos = get_ifos()
  mkdir(opts.physdir+"/" +this_tag)
  page.br()
  page.div(id="encadre")
  page.div(id="contenu")
  page = write_title(page, html_sections[5], "rabelais")
  page.div(id="div_rabelais", style='display:block')
  page.add("This section summarizes the analysis of the full data.<br/>")
  
  #table and venn diagramm
  page = add_input_h3(page, "General information")
  page.table()
  segs = get_segments()
  keys = ("segments","H1","H2","L1","G1","V1","T1") 
  page.add("The segment files above were created with no data quality flags set")
  page.add("The times analyzed accoring to hipe are:")
  coincs = get_coincident_segments(this_tag)  
  ## here is the table. The first column is another table with the durations, 
  ## and the second one is the venn diagram  
  page.add("<table><tr><td>\n")
  page.table()
  page.tr(); page.td('coincidence'); page.td('duration(s)'); page.tr.close()
  for key in coincs.keys():    
    page.tr()
    page.td(key) 
    page.td(coincs.get(key))
    page.tr.close()
  page.table.close()
  page.add("</td><td>\n")
  #create figure for the venn diagram
  data = ""
  for coinc in ("H1","H1H2","H2","H2L1","L1","H1L1","H1H2L1"):
      data = data +coincs.get(coinc) + " "
  create_venn(data, this_tag)
  # and incorporate into html
  comment = "Venn diagramm showing distribution of"
  for coinc in  coincs.keys():
    comment = comment + " "+ coinc
  page = add_figure(page, fnames =[this_tag+"/venn_"+this_tag+".png"], caption=comment, size="half")
  page.add("</td></tr></table>")
  page.add("here below is the detail of the data and ligo-data section of ihope.ini ")
  page = add_config_section(page, "data")
  page = add_config_section(page, "ligo-data")
  page.div.close()  #end of h3
 
  #----------------------------  The horizon distance   
  page = add_input_h3(page, "Inspiral Horizon distance")
  comment = "Inspiral Horizon Distance for a (1.4,1.4) system with a SNR=8"
  page = add_figure(page, fnames =[\
	this_tag+"/horizon_"+this_tag+"_range_hist.png", \
	this_tag+"/horizon_"+this_tag+"_range_plot.png", \
	this_tag+"/horizon_"+this_tag+"_range_mass.png"],
  caption=comment, size="half")
  page.div.close()#end of h3

  #---------------------------- the number of templates
  page = add_input_h3(page, "Variation in template bank and triggered template bank size")
  page = add_config_section(page, "tmpltbank")
  # add figures and summary files
  comment = "Variation in template bank and triggered template bank size"
  page = add_figure(page, fnames =[this_tag+"/plotnumtemplates.png"], caption=comment, size="half")
  page.div.close() #end of h3


  # the end
  page.div.close()
  page.div.close()
  page.div.close()
  print "....Processing  full data section over"
  return page
     
     
# ***************************************************************************
# ***************************************************************************
def write_playground_summary(page,opts):
  print "....Processing  playground section"
  webdir = opts.webdir
  ini = opts.config_file
  ifos = get_ifos()
  mkdir(opts.physdir+"/playground")


  page.br()
  page.div(id="encadre")
  page.div(id="contenu")
  page = write_title(page, html_sections[2], "hugo")
  page.div(id="div_hugo", style='display:block')
  page.add("This section summarizes the analysis of the playground data.<br/>")
  
  #table and venn diagramm
  page = heading(page, "General information", "see details here", "h3")
  page.table()
  segs = get_segments()
  keys = ("segments","H1","H2","L1","G1","V1","T1") 
  page.add("The segment files above were created with no data quality flags set")
  page.add("The times analyzed accoring to hipe are:")
  coincs = get_coincident_segments('playground')  
  ## here is the table. The first column is another table with the durations, 
  ## and the second one is the venn diagram  
  page.add("<table><tr><td>\n")
  page.table()
  page.tr(); page.td('coincidence'); page.td('duration(s)'); page.tr.close()
  for key in coincs.keys():    
    page.tr()
    page.td(key) 
    page.td(coincs.get(key))
    page.tr.close()
  page.table.close()
  page.add("</td><td>\n")
  #create figure for the venn diagram
  data = ""
  for coinc in ("H1","H1H2","H2","H2L1","L1","H1L1","H1H2L1"):
      data = data +coincs.get(coinc) + " "
  create_venn(data, "playground")
  # and incorporate into html
  comment = "Venn diagramm sihowing distribution of"
  for coinc in  coincs.keys():
    comment = comment + " "+ coinc
  page = add_figure(page, fnames =["/playground/venn_playground.png"], caption=comment, size="half")
  page.add("</td></tr></table>")
  page.add("here below is the detail of the data and ligo-data section of ihope.ini ")
  page = add_config_section(page, "data")
  page = add_config_section(page, "ligo-data")
  page.div.close()  
 
  #----------------------------  The horizon distance  
  page = heading(page, "Inspiral range plots" )
  page.add("")  
  page.div(class_="figure")
  page.add("<!--#include virtual=\"playground/plotinspiralrange_playground.html\"-->")
  page = add_caption(page, "Inspiral Horizon distance for a \
	(1.4,1.4) solar mass system with SNR=8 (first sub-figure), and \
	histograms(second sub-figure). The last sub-figure shows the \
	expected horizon distance for any total mass, using an SNR=8."   )
  page.div.close()
  page.div.close()

  #---------------------------- the number of templates
  page = heading(page, "Variation in template bank and triggered template \
	bank size")
  page = add_config_section(page, "tmpltbank")
  page.div(class_="figure")
  page.add("<!--#include virtual=\"playground/plotnumtemplates_playground.html\"-->")
  page = add_caption(page, "Variation in template bank and triggered \
	 template bank size")
  page.div.close()
  page.div.close()
  
  #---------------------------- the inspiral results, first stage, no veto
  page = heading(page, "Inspiral jobs (first stage)")
  page.p("This section summarizes the results of the first inspiral stage. \
	First, let us look at the rate of triggers (unclustered, no veto)")
  #---------------------------- the inspiral trigger rate, no veto
  page.div(class_="figure")
  page.add("<!--#include virtual=\"playground/plotnumtriggers_playground.html\"-->")
  page = add_caption(page, "Variation in inspiral file size")
  page.div.close()
  #---------------------------- the inspiral trigger rate, category 1 and 2
  page.p("The parameters for the inspiral jobs are ")
  page = add_config_section(page, "inspiral")
  page.p("We can look more closely at the first inspiral stage results, and in \
	particular at the trigger distribution (clustered, category veto 1 and 2")
  # add figures and summary files
  for ifo in get_ifos():
    page.h4(ifo)
    summary = open(opts.datadir+"/"+ "playground"  +"/"+ifo+"-SIRE-"+opts.gps_start_time+"-"+opts.duration+".txt")
    comment = ifo+" snr versus time (left), and  cumulative snr (right) clustered case no veto"
    comment += "<center><pre>-------------------summary file (category 1) \n"  +summary.read()+"</pre></center>"
    summary.close()
    
    summary = open(opts.datadir+"/"+ "playground"  +"/"+ifo+"-SIRE_CAT_2_VETO-"+opts.gps_start_time+"-"+opts.duration+".txt")
    comment = ifo+" snr versus time (left), and  cumulative snr (right) clustered case no veto"
    comment += "<center><pre>-------------------summary file (category 1) \n"  +summary.read()+"</pre></center>"
    summary.close()

    page.p("The " +ifo + " inspiral files after the first inspiral stage are clustered (no veto \
	applied), and the results are shown here below in the next set of plots.")
    page = add_config_section(page, ifo.lower()+"-inspiral")
    page = add_figure(page, fnames = [\
	"playground/plotinspiral_playground_"+ifo+"_end_time_vs_snr.png", \
	"playground/plotinspiral_playground_cat_2_"+ifo+"_end_time_vs_snr.png", \
	"playground/plotinspiral_playground_"+ifo+"_snr_cum_hist.png" ,\
	"playground/plotinspiral_playground_cat_2_"+ifo+"_snr_cum_hist.png" \
	], caption=comment, size="half")
  page.div(id="todo")
  page.p("distribution in tau0, tau3 ? which tool")
  page.div.close()# close the todo
  page.div.close()# close the div openned in add_input_h3 function
  
  #---------------------------- the first coincidence stagee, no veto
  page = heading(page, "First coincidence stage (no veto)")
  page = add_config_section(page, "thinca")
  page.add("the number of slide is "+get_numslide("./playground")+" time slides")
  page = add_config_section(page, "thinca-slide")
  fnames=[]
  for combi in ("H1H2L1", "H1H2", "H1L1", "H2L1"):
    fnames=[]
    fnames.append("/playground/"+combi+"_first_coinc_playground_cum_hist_effective_snr.png") 
    fnames.append("/playground/"+combi+"_first_coinc_playground_hist_slide_trigs.png") 
    fnames.append("/playground/"+combi+"_first_coinc_playground_plot_slide_trigs.png") 
    fnames.append("/playground/"+combi+"_first_coinc_playground_cum_hist_snr.png") 
    

    comment = "first stage coincidence" +combi + "case"    
    page = add_figure(page, fnames = fnames, caption=comment, size="half")
  page.div.close()# close the div openned in add_input_h3 function
  

  #---------------------------- the second stage inspiral results, first stage, no veto
  page = add_input_h3(page, "Inspiral jobs (second stage) no veto")
  page = add_config_section(page, "inspiral")
  # add figures and summary files
  for ifo in get_ifos():
    page.h4(ifo)
    comment = ifo+" snr versus time (left), and  cumulative snr (right) clustered case no veto"
    page.p("The " +ifo + "inspiral files after the first inspiral stage are clustered (no veto \
	applied), and the results are shown here below in the next set of plots.")
    page = add_config_section(page, ifo.lower()+"-inspiral")
    page = add_figure(page, fnames = [\
	"playground/plotinspiral_"+ifo+"_end_time_vs_snr.png", \
	"playground/plotinspiral_veto"+ifo+"_end_time_vs_snr.png", \
	"playground/plotinspiral_"+ifo+"_snr_vs_chisq.png", \
	"playground/plotinspiral_"+ifo+"_snr_cum_hist.png" \
	"playground/plotinspiral_veto"+ifo+"_snr_cum_hist.png" \
	"playground/plotinspiral_"+ifo+"_snr_vs_chisq.png" \
	], caption=comment, size="half")
  page.div(id="todo")
  page.p("add the sire summary files?")
  page.div.close()# close the todo
  page.div.close()# close the div openned in add_input_h3 function
  

  #---------------------------- the second coincidence stagee, no veto
  page = add_input_h3(page, "Second coincidence stage no veto")
  page = add_config_section(page, "thinca-2")
  page.add("the number of slide is "+get_numslide("./playground")+" time slides")
  page = add_config_section(page, "thinca-slide")
  fnames=[]
  for combi in ("H1H2L1", "H1H2", "H1L1", "H2L1"):
    fnames=[]
    if combi=="H1H2L1":
      fnames.append("second_coinc_H_vs_L1_snr.png") 
      fnames.append("second_coinc_"+combi+"_plot_slide_trigs.png") 
    elif combi=="H1H2":
      fnames.append("second_coinc_H1_vs_H2_snr.png") 
      fnames.append("second_coinc_"+combi+"_plot_slide_trigs.png") 
    elif combi=="H1L1":
      fnames.append("second_coinc_H1_vs_L1_snr.png") 
      fnames.append("second_coinc_"+combi+"_plot_slide_trigs.png") 
    elif combi=="H2L1":
      fnames.append("second_coinc_H2_vs_L1_snr.png") 
      fnames.append("second_coinc_"+combi+"_plot_slide_trigs.png") 
    comment = "second stage coincidence" +combi + "case"
    page = add_figure(page, fnames = fnames, caption=comment, size="half")
  page.div.close()# close the div openned in add_input_h3 function
  
  #---------------------------- loudest events
  page = add_input_h3(page, "loudest events")
  page.div(id="todo")
  page.p( "add list of loudest events. Not mix this with follow up")
  page.div.close()# close the todo
  page.div.close()# close the h3

  # the end
  page.div.close()
  page.div.close()
  page.div.close()
  print "....Processing  playground section done"
  return page

# ***************************************************************************
# ***************************************************************************
def write_injection(page, opts):
  webdir = opts.webdir
  ini = opts.config_file
  ifos = get_ifos()
  mkdir(opts.physdir+"/injections")
  page.br()

  page.div(id="encadre")
  page.div(id="contenu")
  page = write_title(page, html_sections[3], "baudelaire")
  page.div(id="div_baudelaire", style='display:block')
  page.add("This section summarizes the analysis of the injection runs.<br/>")

  #---------------------------- found and missedi first inspiral stage
  page = add_input_h3(page, "Found versus missed plots (first inspiral stage)")
  for injtag in ['bbh_inj', 'bns_inj', 'bhns_inj']:
    page.h4(injtag)
    for ifo in get_ifos():  
      lists = glob.glob(opts.datadir+"/"+ injtag +"/"+ifo+"-SIRE_INJECTIONS_*FOUND-*.txt")
      summary = open(lists[0])
      caption = ifo +" Found and missed injections at the first stage of the inspiral analysis ("+injtag +"tag)."      
      caption = caption + "<pre> -----------------------"+lists[0]+"\n" +summary.read()+"</pre>"
      page = add_figure(page, fnames=  [\
        "injections/"+ifo+"_first_insp_"+injtag+"_dist-mchirp.png", \
        "injections/"+ifo+"_first_insp_"+injtag+"_time-dist.png"], \
        caption=caption, size="half")
  page.div.close()
  
  #---------------------------- found and missed first coinc stage
  page = add_input_h3(page, "Found versus missed plots (first coincidence stage)")
  for injtag in ['bbh_inj', 'bns_inj', 'bhns_inj']:
    page.h4(injtag)
    for combi in get_ifo_coinc():  
      if len(combi)>2:
        print opts.datadir+"/"+ injtag +"/"+combi+"-COIRE_INJECTIONS_*_FOUND_"+combi+"-*.txt"
        lists = glob.glob(opts.datadir+"/"+ injtag +"/"+combi+"-COIRE_INJECTIONS_*_FOUND_"+combi+"-*.txt")
        print lists
        summary = open(lists[0])
        caption = combi +" Found and missed injections at the first coincidence stage ("+injtag +"tag)."      
        caption = caption + "<pre> -----------------------"+lists[0]+"\n" +summary.read()+"</pre>"
        page = add_figure(page, fnames=  [\
          "injections/"+combi+"_first_coinc_bbh_inj_missed_found_chirp_dist_vs_mchirp.png", \
          "injections/"+combi+"_first_coinc_bbh_inj_missed_found_vs_end_time.png", \
    	  "injections/"+combi+"_first_coinc_bbh_inj_missed_found_vs_mchirp.png"], \
          caption=caption, size="half")
  page.div.close()
 
  #---------------------------- close by missed
  page = add_input_h3(page, "close by missed")
  page = add_figure(page, fnames = "H1H2L1-found_missed.png", caption="", size="half")
  page = add_figure(page, fnames = "H1H2-found_missed.png", caption="", size="half")
  page = add_figure(page, fnames = "H1L1-found_missed.png", caption="", size="half")
  page.div.close()
  
  #---------------------------- parameter efficiency
  page = add_input_h3(page, "Parameter efficiency")
  page = add_figure(page, fnames = "H1H2L1-found_missed.png", caption="", size="half")
  page = add_figure(page, fnames = "H1H2-found_missed.png", caption="", size="half")
  page = add_figure(page, fnames = "H1L1-found_missed.png", caption="", size="half")
  page.div.close()

  
  # the end
  page.div.close()
  page.div.close()
  page.div.close()

  return page


def write_about(page, opt):
  webdir = opts.webdir
  page.br()
  page.div(id="encadre")
  page.div(id="contenu")
  page = write_title(page, html_sections[8], "balzac")
  page.div(id="div_balzac", style='display:block')
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
  print opts.datadir + "/" + run + "/inspiral_hipe_" + run
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
  return page 

# ***************************************************************************
# ***************************************************************************
def error_section(page, section):
  page.br()
  page.div(id="encadre")
  page.div(id="contenu")
  page.br()
  page.div(id="todo")
  page.p("Error(s) during the processing of " +section +" section. Skipping this section.")
  page.div.close()
  page.br()
  page.div.close()
  page.div.close()
  return page

# ***************************************************************************
# ***************************************************************************
def add_caption(page, caption):
  global fig_num
  page.p("Figure "+str(fig_num) + ": "+caption)
  fig_num = fig_num + 1
  return page

def add_figure(page,fnames, caption, size):
  global fig_num
  dir = opts.webdir
  page.add("<!-- insert a figure -->\n<div class=\"figure\">")
  for fnam in fnames:
    page.add("\t<a href=\"" + dir + "/" + fnam +"\">\n" + \
      "  <img class=\"" + size + "\" src=\"" + dir + "/" + fnam + "\"\>\n  </a>")
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
    command=("matlab -nodisplay -nodesktop -nojvm -nosplash   < temp.m > /tmp/test ; rm -f temp.m; mv venn_"+tag+".png "+opts.physdir+"/"+tag+"/")
    if not opts.debug:
      make_external_call(command, True, True, True)    
  except:
    print """   The matlab command to create the venn diagram failed. 
                Check that matlab is properly set,and vennX.m is available
                (see matapps/src/searches/inspiral/matspiral/utilities")  
          """

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
  print "....Get "+ tag +" -------------------------------------------"
  ifos = get_ifos()
  thisdata = {}
  thisdata['segments'] = ['duration(days)' ,'duration (s)']
  duration = str((int(opts.gps_end_time) - int(opts.gps_start_time)))
  output={}
  ifo_coincs = get_ifo_coinc()  
  for coinc in ifo_coincs:
    if tag=="playground":
      command = "awk \'{sum=sum+$4} END {print sum}\' "+opts.datadir+ tag +"/"+coinc+"_play_segs_analyzed.txt"
    elif tag=="analysis":
      command = "awk \'{sum=sum+$4} END {print sum}\' "+opts.datadir+ tag +"/"+coinc+"_segs_analyzed.txt"

    output[coinc] = make_external_call(command, 0,0)


  return output

  


# ***************************************************************************
# ***************************************************************************
def get_segments_tag(tag):
  """
  reads segment files and return table of duration for each ifo
  """
  datadir = opts.datadir
  print "....Get "+ tag +" -------------------------------------------"
  ifos = get_ifos()
  thisdata = {}
  thisdata['segments'] = ['duration(days)' ,'duration (s)']
  duration = str((int(opts.gps_end_time) - int(opts.gps_start_time)))
  duration.split('.')
  this_tag  = '-' + tag + '-' + opts.gps_start_time + '-' + str(duration) + '.txt'
  if tag=="RDS_C03_L2":
    this_tag  = '_' + tag + '-' + opts.gps_start_time + '-' + str(duration) + '.txt'
  
  command = 'awk \'{print NF}\' ' +  datadir +'/' + ifos[1] + this_tag
  ncols = make_external_call(command, 0, 0)
  ncols = ncols[len(ncols)-2]
  if float(ncols)==4:  
    for ifo in ifos:
      print '......Get segment durations for '+ifo
      command = 'awk \'{sum=sum+$4} END {print sum/3600/24}\' ' +  datadir +'/' + ifo + this_tag
      output_days = make_external_call(command, 0, 0, False)
      command = 'awk \'{sum=sum+$4} END {print sum}\' ' +  datadir +'/' + ifo + this_tag
      output_seconds = make_external_call(command, 0, 0, False)
      thisdata[ifo] = [ output_days, output_seconds]
  elif float(ncols)==2:
    for ifo in ifos:
      print '......Get segment durations for '+ifo
      command = 'awk \'{sum=sum+$2-$1} END {print sum/3600/24}\' ' +  datadir +'/' + ifo + this_tag
      output_days = make_external_call(command, 0, 0, False)
      command = 'awk \'{sum=sum+$2-$1} END {print sum}\' ' +  datadir +'/' + ifo + this_tag
      output_seconds = make_external_call(command, 0, 0, False)
      thisdata[ifo] = [ output_days, output_seconds]

  return thisdata


# ***************************************************************************
# ***************************************************************************
def get_segments():
  """
  reads the segment lists. Same as get_segments_tag but filename are slightly different.
  """
  print "....Get  selected segments ---------------------------------------------"
  datadir = opts.datadir
  ifos = get_ifos()
  thisdata = {}
  thisdata['segments'] = ['duration(days)' ,'duration (s)']
  for ifo in ifos:
    print '......Get segment durations for '+ifo
    command = 'awk \'{sum=sum+$4} END {print sum/3600/24}\' ' +  datadir +'/'+ifo+'-SELECTED_SEGS.txt '
    output_days = make_external_call(command, 0, 0)
    command = 'awk \'{sum=sum+$4} END {print sum}\' ' +  datadir +'/'+ifo+'-SELECTED_SEGS.txt '
    output_seconds = make_external_call(command, 0, 0)
    thisdata[ifo] = [ output_days, output_seconds]
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
  try:
    if executable=="trigbank":
      pathname = hipecp.get("condor", "trigtotmplt") 
    else:
      pathname = hipecp.get("condor", executable) 
    s =  os.path.split(pathname)
    if len(s[0])==0:
      cmd = opts.ihope_directory + s[1] + '  --version'
    else:
      cmd = s[1] + ' --version '
    output = make_external_call(cmd, False, False)
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
  # First, we copy the file (style, data) that are needed to create the web page
  try :
    print '...Copying category files '
    cats = hipecp.items("segments")
    location = opts.physdir+'/catlists/'
    mkdir(location)
    for file in cats:
      if file[0].startswith('analyze')<=0:
        command = 'cp '+ file[1] +' ' + location
        try: 
          make_external_call(command, 0,0)
        except:pass
  except:
    print "WARNING:: could not copy the category segment list. continuing..."
  # the segment files
  try :
    print '...Copying segment files '
    for this in get_ifos():
        seg = this +'-SELECTED_SEGS.txt'
        mkdir(opts.physdir+'/segments')
        make_external_call('cp '+opts.datadir + seg + ' '+opts.physdir+'/segments/', 0,0)
  except:
    print "WARNING:: could not copy the segments list. continuing..."
# ***************************************************************************
# ***************************************************************************
def set_style():
  """
  """
  tmp = []
  tmp.append(configcp.get("main", "style"))
  make_external_call('cp ' + configcp.get("main", "style")+ " " +opts.physdir)
  for i in range(1,10,1):
    try:
      tmp.append (configcp.get("main", "style"+str(i)))    
      make_external_call('cp ' + configcp.get("main", "style"+str(i))+ " " +opts.physdir)
    except:pass
  style = tmp

  count = 0
  for this in style:
    tmp =  this.split('/')
    style[count] = tmp[len(tmp)-1]
    count = count +1

  return style

# ***************************************************************************
# ***************************************************************************
#############################################################################
usage =  """ %prog [options]
Program to write webpage from upperlimit.py
"""

parser = OptionParser( usage = usage, version = "%prog CVS "+__Id__ )

parser.add_option("-C","--config-file",action="store",type="string",\
    default='write_ihope_page.ini', metavar=" INI File",\
    help="ini file with information about run directories" )
parser.add_option("-A","--skip-analysis",action="store_false",\
    default=True, dest="analysis",metavar="DOANALYSIS",\
    help="" )
parser.add_option("-T","--skip-tuning",action="store_false",\
    default=True, dest="tuning",metavar="DOTUNING",\
    help="" )
parser.add_option("-U","--skip-upperlimit",action="store_false",\
    default=True, dest="upperlimit",metavar="DOUPPERLIMIT",\
    help="" )
parser.add_option("-S","--skip-summary",action="store_false",\
    default=True, dest="summary",metavar="DOSUMMARY",\
    help="" )
parser.add_option("-P","--skip-playground",action="store_false",\
    default=True, dest="playground",metavar="DOPLAYGROUND",\
    help="" )
parser.add_option("-D","--debug-mode",action="store_true",\
    default=False, dest="debug",metavar="DODEBUG",\
    help="" )



(opts,args) = parser.parse_args()

#############################################################################
#  MAIN PART                                                                #
#############################################################################
fig_num = 0
count_block = 0
config   =  opts.config_file
opts.config = config # save the name of the ini file, why ?
configcp = ConfigParser.ConfigParser()
configcp.read(config)



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
  #read hte different css style available and copy the style files in the web directory.
except:
  print "ERROR::The ini file "+ opts.config_file + " does not have the proper field in the [main] section"
  print "ERROR::Consider adding ihope-ini-file, gps-start-time,gps-end-time, ihope-directory, title,url, username, output and style fields"
  raise
  
#sub-products of the ini file parsing
opts.gpsdir =  '/'+str(opts.gps_start_time)+'-'+str(opts.gps_end_time)+'/'
opts.duration = str(int(opts.gps_end_time) - int(opts.gps_start_time))
opts.webdir = 'http://' + opts.url + '/~' + opts.username  + opts.gpsdir
opts.datadir = opts.ihope_directory + opts.gpsdir
opts.physdir = '/archive/home/'+opts.username+'/public_html/'+ opts.gpsdir
  
#read hte different css style available and copy the style files in the web directory.
try:
  opts.style = set_style()
except:
  print "WARNING::problem with the style file. (either not copied or not found)"
  raise

html_file = file(opts.output,"w")


# here is the directory we want to look at
if not  os.path.isdir(opts.datadir):
  raise  "%s is not a valid directory. Check your gps time." % opts.datadir
# which physical name is 
if not os.path.isdir(opts.physdir):
  "%s is not a valid directory. Trying to create it." % opts.physdir
  try:
    os.mkdir(opts.physdir)
    print "##Warning:: directory "+opts.physdir+" created"
  except: pass


# now we can parse the ihope.ini file 
print "...Parsing the ihope ini file" 
hipe   = opts.ihope_directory+'/'+opts.config_file
hipecp = ConfigParser.ConfigParser()
hipecp.read(hipe)
make_external_call( 'cp '+opts.config_file + ' ' + opts.physdir)

copy_segments()


###### create the section labels  
print '...Creating html documents'
html_sections = ("General information",  \
		"Data summary", \
		"Playground", \
		"Injection", \
		"Tuning", \
		"Full data", \
		"Upper limit", \
		"Openning the box",
		"About this page")
title = opts.title+ ": from "+str(opts.gps_start_time)+" to "+str(opts.gps_end_time) 
script = {}
script['javascript'] = 'toggle.js'
create_toggle()
# Finallt, we create the html document 
page = markup.page()
page._escape = False
page.init(title=title, css=opts.style, script=script)
page.h1(opts.title +" (" + opts.gps_start_time +"-" +opts.gps_end_time+")")

##### we write the toc, and general section
try: 
  page = write_toc(page, opts)
except:
  page.div("problem in toc section. skippingg this section. "); 
  pass

try:
  page = write_general(page, opts)
except:
  page.add("problem in "+html_sections[0]+" section. skipping this section. "); 
  if opts.debug:
    raise
  else:
    pass

if opts.summary:
  try:
    page = write_data_summary(page, opts)
  except: 
    page.add("problem in "+html_sections[1]+" section. skipping this section. "); 
    pass

if opts.playground:
  try:
    page = write_playground_summary(page, opts)
  except:
    page = error_section(page, html_sections[2]); 
    if opts.debug:
      raise
    else:
      pass
  
#page = write_injection(page, opts)

if opts.tuning:
  try:
    page = write_tuning(page, opts)
  except: 
    page.add("problem in "+html_sections[4]+" section. skipping this section. "); 
    pass

if opts.analysis:
  try:
    page = write_fulldata(page, opts)
  except: 
    page.add("problem in "+html_sections[5]+" section. skipping this section. "); 
    pass

if opts.upperlimit:
  try:
    page = write_upperlimit(page, opts)
  except:  
    page.add("problem in "+html_sections[6]+" section. skipping this section. "); 
    pass

#try:
#  page = write_openbox(page, opts)
#except:
#  page = error_section(page, html_sections[7]);
#  pass

#try:
page = write_about(page, opts)
#except:   page = error_section(page, html_sections[8]); pass
 

##### the end
html_file.write(page(False))
html_file.close()

print 'html file created. Parsing to convert in nice html output'



print '...Copying html documents in ' +opts.physdir
make_external_call('mv  '+opts.output +' ' + opts.physdir, '','',True)
make_external_call( 'mv toggle.js '+ opts.physdir, '', 'True')



