#!/usr/bin/env /usr/bin/python

import shutil
import os
import sys
import glob
import ConfigParser
import subprocess

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

try: webserver = sys.argv[1]
except: 
  print >>sys.stderr, "YOU MUST SPECIFY A WEBSERVER AS THE FIRST ARGUMENT (e.g. https://ldas-jobs.ligo.caltech.edu/~channa/highmass_months_23-24_summary_page)"
  sys.exit(1)

open_box = False
try: # see if you want to open the box
  if sys.argv[2] == "open":
    print >>sys.stderr, "WARNING: OPENING THE BOX"
    open_box = True
except: pass

page = wiki(open_box)

page.section("Injection Parameters")

image_list = page.image_glob('cbc_plotsummary_5_sim_dist*.png') 
page.image_table(image_list,webserver)

page.section("Found / Missed")
image_list = page.image_glob('cbc_plotsummary_0_deff_vs_mchirp_*.png')
page.image_table(image_list,webserver)
image_list = page.image_glob('cbc_plotsummary_0_deff_vs_t_*.png')
page.image_table(image_list,webserver)

page.section("Parameter Accuracy")
image_list = page.image_glob('cbc_plotsummary_1_mchirp_acc_frac_*.png')
page.image_table(image_list,webserver)
image_list = page.image_glob('cbc_plotsummary_1_eta_acc_frac_*.png')
page.image_table(image_list,webserver)
image_list = page.image_glob('cbc_plotsummary_1_t_acc_*.png')
page.image_table(image_list,webserver)

page.section("Playground Chi-squared")
image_list = page.image_glob('cbc_plotsummary_2_playground_chi2_vs_rho_*.png')
page.image_table(image_list,webserver)

page.section("Playground Effective SNR scatter")
image_list = page.image_glob('cbc_plotsummary_3_playground_rho_*.png')
page.image_table(image_list,webserver)

page.section("Playground SNR")
image_list = page.image_glob('cbc_plotsummary_4_playground_count_vs_snr_*.png')
page.image_table(image_list,webserver)

page.section("Playground Ifar")
image_list = page.image_glob('cbc_plotsummary_4_playground_count_vs_ifar*.png')
page.image_table(image_list,webserver)

try:
  for l in open("playground_summary_table.txt").readlines(): page.write(l)
except: print >>sys.stderr, "WARNING: couldn't find playground summary, continuing"

if open_box:
    print >>sys.stderr, "WARNING: OPENING THE BOX"

    page.section("Full Data Chi-squared")
    image_list = page.image_glob('cbc_plotsummary_2_chi2_vs_rho_*.png')
    page.image_table(image_list,webserver)

    page.section("Full Data Effective SNR scatter")
    image_list = page.image_glob('cbc_plotsummary_3_rho_*.png')
    page.image_table(image_list,webserver)

    page.section("Full Data SNR")
    image_list = page.image_glob('cbc_plotsummary_4_count_vs_snr_*.png')
    page.image_table(image_list,webserver)

    page.section("Full Data Ifar")
    image_list = page.image_glob('cbc_plotsummary_4_count_vs_ifar_*.png')
    page.image_table(image_list,webserver)

    try:
      for l in open("summary_table.txt").readlines(): page.write(l)
    except: print >>sys.stderr, "WARNING: couldn't find summary, continuing"

    # UPPER LIMIT PLOTS
    ifos_list = [f.replace('volume_time.png','') for f in page.image_glob('*volume_time.png')]
    ifos_string = ",".join(ifos_list)
    page.section("Volume x time " + ifos_string)
    try:
      filenames = page.image_glob('*range_summary.txt"')
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
    image_list = page.image_glob('*volume_time.png') 
    page.image_table(image_list,webserver)

    page.section("error on Volume x time " + ifos_string)
    image_list = page.image_glob('*fractional_error.png')
    page.image_table(image_list,webserver)

    page.section("lambda " + ifos_string)
    image_list = page.image_glob('*lambda.png')
    page.image_table(image_list,webserver)

    page.section("Selected posteriors " + ifos_string)
    image_list = page.image_glob('*posterior.png') 
    page.image_table(image_list,webserver)

    page.section("upper limit " + ifos_string)
    image_list = page.image_glob('*upper_limit.png') 
    page.image_table(image_list,webserver)

    page.section("Combined upper limit")
    image_list = ['combinedupper_limit.png', 'combinedposterior.png']
    page.image_table(image_list,webserver)

try: 
  for l in open("plotsummary.txt").readlines(): page.write(l)
except: print >>sys.stderr, "WARNING couldn't find plotsummary.txt"

page.finish()
