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
    command = 'convert ' + path + ' -resize 300x300 ' + thumb
    print command
    popen = subprocess.Popen(command.split())
    popen.communicate()
    status = popen.returncode
    s = '[[ImageLink('+webserver+'/'+thumb+','+webserver+'/'+path+',width=300][,alt=none])]]'
    self.file.write(s)

  def image_table(self,image_list, webserver):
    if not image_list: return
    self.file.write("||")
    for i in image_list:
      self.image_link(i, webserver)
      self.file.write("||")
    self.file.write("\n")

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
image_list = ['cbc_plotsummary_5_sim_dist_m1_m2_EOBNRpseudoFourPN.png','cbc_plotsummary_5_sim_dist_m1_m2_IMRPhenomAtwoPN.png']
page.image_table(image_list,webserver)

page.section("Found / Missed")
image_list = ['cbc_plotsummary_0_deff_vs_mchirp_H1H2L1.png','cbc_plotsummary_0_deff_vs_mchirp_H1L1.png','cbc_plotsummary_0_deff_vs_mchirp_H2L1.png']
page.image_table(image_list,webserver)
image_list = ['cbc_plotsummary_0_deff_vs_t_H1H2L1.png','cbc_plotsummary_0_deff_vs_t_H1L1.png','cbc_plotsummary_0_deff_vs_t_H2L1.png']
page.image_table(image_list,webserver)

page.section("Parameter Accuracy")
image_list = ['cbc_plotsummary_1_mchirp_acc_frac_IMRPhenomAtwoPN_H1.png','cbc_plotsummary_1_mchirp_acc_frac_IMRPhenomAtwoPN_H2.png','cbc_plotsummary_1_mchirp_acc_frac_IMRPhenomAtwoPN_L1.png']
page.image_table(image_list,webserver)
image_list = ['cbc_plotsummary_1_eta_acc_frac_IMRPhenomAtwoPN_H1.png','cbc_plotsummary_1_eta_acc_frac_IMRPhenomAtwoPN_H2.png','cbc_plotsummary_1_eta_acc_frac_IMRPhenomAtwoPN_L1.png']
page.image_table(image_list,webserver)
image_list = ['cbc_plotsummary_1_t_acc_IMRPhenomAtwoPN_H1.png','cbc_plotsummary_1_t_acc_IMRPhenomAtwoPN_H2.png','cbc_plotsummary_1_t_acc_IMRPhenomAtwoPN_L1.png']
page.image_table(image_list,webserver)
image_list = ['cbc_plotsummary_1_mchirp_acc_frac_EOBNRpseudoFourPN_H1.png','cbc_plotsummary_1_mchirp_acc_frac_EOBNRpseudoFourPN_H2.png','cbc_plotsummary_1_mchirp_acc_frac_EOBNRpseudoFourPN_L1.png']
page.image_table(image_list,webserver)
image_list = ['cbc_plotsummary_1_eta_acc_frac_EOBNRpseudoFourPN_H1.png','cbc_plotsummary_1_eta_acc_frac_EOBNRpseudoFourPN_H2.png','cbc_plotsummary_1_eta_acc_frac_EOBNRpseudoFourPN_L1.png']
page.image_table(image_list,webserver)
image_list = ['cbc_plotsummary_1_t_acc_EOBNRpseudoFourPN_H1.png','cbc_plotsummary_1_t_acc_EOBNRpseudoFourPN_H2.png','cbc_plotsummary_1_t_acc_EOBNRpseudoFourPN_L1.png']
page.image_table(image_list,webserver)

page.section("Playground Chi-squared")
image_list = ['cbc_plotsummary_2_playground_chi2_vs_rho_H1.png', 'cbc_plotsummary_2_playground_chi2_vs_rho_H2.png', 'cbc_plotsummary_2_playground_chi2_vs_rho_L1.png']
page.image_table(image_list,webserver)

page.section("Playground Effective SNR scatter")
image_list = ['cbc_plotsummary_3_playground_rho_H1_vs_L1.png', 'cbc_plotsummary_3_playground_rho_H1_vs_H2.png', 'cbc_plotsummary_3_playground_rho_H2_vs_L1.png']
page.image_table(image_list,webserver)

page.section("Playground Ifar")
image_list = ['cbc_plotsummary_4_playground_count_vs_ifar_H1H2L1.png','cbc_plotsummary_4_playground_count_vs_ifar_H1L1.png','cbc_plotsummary_4_playground_count_vs_ifar_H2L1.png']
page.image_table(image_list,webserver)

page.section("Playground SNR")
image_list = ['cbc_plotsummary_4_playground_count_vs_snr_H1H2L1.png','cbc_plotsummary_4_playground_count_vs_snr_H1L1.png','cbc_plotsummary_4_playground_count_vs_snr_H2L1.png']
page.image_table(image_list,webserver)

try:
  for l in open("playground_summary_table.txt").readlines(): page.write(l)
except: print >>sys.stderr, "WARNING: couldn't find playground summary, continuing"

if open_box:
    print >>sys.stderr, "WARNING: OPENING THE BOX"

    page.section("Full Data Chi-squared")
    image_list = ['cbc_plotsummary_2_chi2_vs_rho_H1.png', 'cbc_plotsummary_2_chi2_vs_rho_H2.png', 'cbc_plotsummary_2_chi2_vs_rho_L1.png']
    page.image_table(image_list,webserver)

    page.section("Full Data Effective SNR scatter")
    image_list = ['cbc_plotsummary_3_rho_H1_vs_L1.png', 'cbc_plotsummary_3_rho_H1_vs_H2.png', 'cbc_plotsummary_3_rho_H2_vs_L1.png']
    page.image_table(image_list,webserver)

    page.section("Full Data Ifar")
    image_list = ['cbc_plotsummary_4_count_vs_ifar_H1H2L1.png','cbc_plotsummary_4_count_vs_ifar_H1L1.png','cbc_plotsummary_4_count_vs_ifar_H2L1.png']
    page.image_table(image_list,webserver)

    page.section("Full Data SNR")
    image_list = ['cbc_plotsummary_4_count_vs_snr_H1H2L1.png','cbc_plotsummary_4_count_vs_snr_H1L1.png','cbc_plotsummary_4_count_vs_snr_H2L1.png']
    page.image_table(image_list,webserver)
    try:
      for l in open("summary_table.txt").readlines(): page.write(l)
    except: print >>sys.stderr, "WARNING: couldn't find summary, continuing"

    # UPPER LIMIT PLOTS
    page.section("Volume x time H1H2L1, H1L1, H2L1")
    try:
      files = [open("H1H2L1range_summary.txt").readlines(), open("H1L1range_summary.txt").readlines(), open("H2L1range_summary.txt").readlines()]
      page.write("||<-2> H1H2L1||||<-2> H1L1||||<-2> H2L1||\n")
      for i in range(len(files[0])):
        for f in files: 
          page.write(f[i].strip())
        page.write("\n")
    except: print >>sys.stderr, "WARNING: couldn't find Range summary " + f + ", continuing"
    page.write("\n")
    image_list = ['H1H2L1volume_time.png', 'H1L1volume_time.png','H2L1volume_time.png']
    page.image_table(image_list,webserver)

    page.section("error on Volume x time H1H2L1, H1L1, H2L1")
    image_list = ['H1H2L1fractional_error.png', 'H1L1fractional_error.png','H2L1fractional_error.png']
    page.image_table(image_list,webserver)

    page.section("lambda H1H2L1, H1L1, H2L1")
    image_list = ['H1H2L1lambda.png', 'H1L1lambda.png','H2L1lambda.png']
    page.image_table(image_list,webserver)

    page.section("Selected posteriors H1H2L1, H1L1, H2L1")
    image_list = ['H1H2L1posterior.png', 'H1L1posterior.png','H2L1posterior.png']
    page.image_table(image_list,webserver)

    page.section("upper limit H1H2L1, H1L1, H2L1")
    image_list = ['H1H2L1upper_limit.png', 'H1L1upper_limit.png','H2L1upper_limit.png']
    page.image_table(image_list,webserver)

    page.section("Combined upper limit")
    image_list = ['combinedupper_limit.png', 'combinedposterior.png']
    page.image_table(image_list,webserver)

try: 
  for l in open("plotsummary.txt").readlines(): page.write(l)
except: print >>sys.stderr, "WARNING couldn't find plotsummary.txt"

page.finish()
