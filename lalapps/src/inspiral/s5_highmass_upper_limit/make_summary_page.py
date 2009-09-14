#!/usr/bin/env /usr/bin/python

import shutil
import os
import sys
import glob
import ConfigParser
import subprocess
from optparse import OptionParser

class wiki(object):
  def __init__(self,open_box=False,fname="wiki.txt"):
    if open_box: fname = "open_box_" + fname
    self.fname = fname
    self.file = open(fname,"w")

  def image_link(self,path,webserver):
    thumb = "thumb_" + path
    command = 'convert ' + path + ' -resize 300x300 -antialias ' + thumb
    print command
    popen = subprocess.Popen(command.split())
    popen.communicate()
    status = popen.returncode
    s = '[[ImageLink('+webserver+'/'+thumb+','+webserver+'/'+path+',width=300][,alt=none])]]'
    self.file.write(s)

  def image_table(self,image_list, webserver):
    if not image_list: return
    for j, i in enumerate(image_list):
      if not (j) % 3: self.file.write("\n\n||")
      self.image_link(i, webserver)
      self.file.write("||")
    self.file.write("\n\n")
  
  def image_glob(self, pat):
    image_list = []
    for image in glob.glob(pat):
      if 'thumb' in image: continue
      else: image_list.append(image)
    image_list.sort()
    return image_list

  def section(self,title):
    s = "=== "+title.strip()+" ===\n"
    self.file.write(s)

  def write(self,val):
    self.file.write(val)

  def finish(self):
    self.file.close()  

def parse_command_line():
  parser = OptionParser(version = "%prog CVS $Id$", usage = "%prog [options] [file ...]", description = "%prog computes mass/mass upperlimit")
  parser.add_option("--webserver", help = "Set the webserver path.  Required.  Example https://ldas-jobs.ligo.caltech.edu/~channa/highmass_months_23-24_summary_page")
  parser.add_option("--open-box", action = "store_true", help = "Produce open box page")
  parser.add_option("--output-name-tag", default = "", metavar = "name", help = "Set the basename for image search")
  opts, filenames = parser.parse_args()

  if not opts.webserver:
    print >>sys.stderr, "must specify a webserver"
    sys.exit(1)
  return opts, filenames

###########################################################
# MAIN
###########################################################

opts, filenames = parse_command_line()

webserver = opts.webserver
open_box = opts.open_box
base_name = opts.output_name_tag

if open_box: print >>sys.stderr, "WARNING: OPENING THE BOX"

page = wiki(open_box, fname=base_name+"wiki.txt")

page.section("Injection Parameters")

image_list = page.image_glob(base_name+'_sim_dist*.png') 
print base_name+'_sim_dist*.png'
page.image_table(image_list,webserver)

page.section("Found / Missed")
image_list = page.image_glob(base_name+'_deff_vs_mchirp_*.png')
page.image_table(image_list,webserver)
image_list = page.image_glob(base_name+'_deff_vs_t_*.png')
page.image_table(image_list,webserver)

page.section("Parameter Accuracy")
image_list = page.image_glob(base_name+'_mchirp_acc_frac_*.png')
page.image_table(image_list,webserver)
image_list = page.image_glob(base_name+'_eta_acc_frac_*.png')
page.image_table(image_list,webserver)
image_list = page.image_glob(base_name+'_t_acc_*.png')
page.image_table(image_list,webserver)

page.section("Playground Chi-squared")
image_list = page.image_glob(base_name+'_playground_chi2_vs_rho_*.png')
page.image_table(image_list,webserver)

page.section("Playground Effective SNR scatter")
image_list = page.image_glob(base_name+'_playground_rho_*.png')
page.image_table(image_list,webserver)

page.section("Playground SNR")
image_list = page.image_glob(base_name+'_playground_count_vs_snr_*.png')
page.image_table(image_list,webserver)

page.section("Playground Ifar")
image_list = page.image_glob(base_name+'_playground_count_vs_ifar*.png')
page.image_table(image_list,webserver)

try:
  for l in open(base_name+"playground_summary_table.txt").readlines(): page.write(l)
except: print >>sys.stderr, "WARNING: couldn't find playground summary, continuing"

if open_box:
    print >>sys.stderr, "WARNING: OPENING THE BOX"

    page.section("Full Data Chi-squared")
    image_list = page.image_glob(base_name+'_chi2_vs_rho_*.png')
    page.image_table(image_list,webserver)

    page.section("Full Data Effective SNR scatter")
    image_list = page.image_glob(base_name+'_rho_*.png')
    page.image_table(image_list,webserver)

    page.section("Full Data SNR")
    image_list = page.image_glob(base_name+'_count_vs_snr_*.png')
    page.image_table(image_list,webserver)

    page.section("Full Data Ifar")
    image_list = page.image_glob(base_name+'_count_vs_ifar_*.png')
    page.image_table(image_list,webserver)

    try:
      for l in open(base_name+"summary_table.txt").readlines(): page.write(l)
    except: print >>sys.stderr, "WARNING: couldn't find summary, continuing"

    # UPPER LIMIT PLOTS
    #ifos_list = [f.replace('volume_time.png','') for f in page.image_glob('*volume_time.png')]
    #ifos_string = ",".join(ifos_list)
    page.section("Volume x time ")
    try:
      filenames = page.image_glob(base_name+'*range_summary.txt')
      print filenames
      files = [open(f).readlines() for f in filenames]
      for f in filenames:
        page.write("|| || %s ||" % (f.replace('range_summary.txt',''),) )
      page.write("\n")
      for i in range(len(files[0])):
        for f in files: 
          page.write(f[i].strip())
        page.write("\n")
    except: print >>sys.stderr, "WARNING: couldn't find Range summary " + f + ", continuing"

    page.write("\n")
    image_list = page.image_glob(base_name+'*volume_time.png') 
    page.image_table(image_list,webserver)

    page.section("error on Volume x time ")
    image_list = page.image_glob(base_name+'*fractional_error.png')
    page.image_table(image_list,webserver)

    page.section("lambda ")
    image_list = page.image_glob(base_name+'*lambda.png')
    page.image_table(image_list,webserver)

    page.section("Selected posteriors ")
    image_list = page.image_glob(base_name+'*posterior.png') 
    page.image_table(image_list,webserver)

    page.section("upper limit ")
    image_list = page.image_glob(base_name+'*upper_limit.png') 
    page.image_table(image_list,webserver)

    page.section("Combined upper limit")
    image_list = [base_name+'combinedupper_limit.png', base_name+'combinedposterior.png']
    page.image_table(image_list,webserver)

try: 
  for l in open(base_name+"plotsummary.txt").readlines(): page.write(l)
except: print >>sys.stderr, "WARNING couldn't find plotsummary.txt"

page.finish()
