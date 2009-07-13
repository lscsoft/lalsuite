#!/usr/bin/env /usr/bin/python

import shutil
import os
import sys
import glob
import ConfigParser
import subprocess

class wiki(object):
  def __init__(self,fname="wiki.txt"):
    self.fname = fname
    self.file = open(fname,"w")

  def image_link(self,path,webserver):
    thumb = "thumb_" + path
    command = 'convert ' + path + ' -resize 400x400 ' + thumb
    print command
    popen = subprocess.Popen(command.split())
    popen.communicate()
    status = popen.returncode
    s = '[[ImageLink('+webserver+'/'+thumb+','+webserver+'/'+path+',width=400][,alt=none])]]'
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

  def finish(self):
    self.file.close()  

try: webserver = sys.argv[1]
except: 
  print >>sys.stderr, "YOU MUST SPECIFY A WEBSERVER AS THE FIRST ARGUMENT (e.g. https://ldas-jobs.ligo.caltech.edu/~channa/highmass_months_23-24_summary_page)"
  sys.exit(1)

page = wiki()

page.section("Found / Missed")
image_list = ['cbc_plotsummary_0_deff_vs_mchirp_H1H2L1.png','cbc_plotsummary_0_deff_vs_mchirp_H1L1.png','cbc_plotsummary_0_deff_vs_mchirp_H2L1.png']
page.image_table(image_list,webserver)
image_list = ['cbc_plotsummary_0_deff_vs_t_H1H2L1.png','cbc_plotsummary_0_deff_vs_t_H1L1.png','cbc_plotsummary_0_deff_vs_t_H2L1.png']
page.image_table(image_list,webserver)

page.section("Parameter Accuracy")
image_list = ['cbc_plotsummary_1_mchirp_acc_frac_H1.png','cbc_plotsummary_1_mchirp_acc_frac_H2.png','cbc_plotsummary_1_mchirp_acc_frac_L1.png']
page.image_table(image_list,webserver)
image_list = ['cbc_plotsummary_1_eta_acc_frac_H1.png','cbc_plotsummary_1_eta_acc_frac_H2.png','cbc_plotsummary_1_eta_acc_frac_L1.png']
page.image_table(image_list,webserver)
image_list = ['cbc_plotsummary_1_t_acc_H1.png','cbc_plotsummary_1_t_acc_H2.png','cbc_plotsummary_1_t_acc_L1.png']
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

try: # see if you want to open the box
  if sys.argv[2] == "open":

    print >>sys.stderr, "WARNING: OPENING THE BOX"

    page.section("Full Data Chi-squared")
    image_list = ['cbc_plotsummary_2_chi2_vs_rho_H1.png', 'cbc_plotsummary_2_chi2_vs_rho_H2.png', 'cbc_plotsummary_2_chi2_vs_rho_L1.png']
    page.image_table(image_list,webserver)

    page.section("Full Data Effective SNR scatter")
    image_list = ['cbc_plotsummary_3_rho_H1_vs_L1.png', 'cbc_plotsummary_3_rho_H1_vs_H2.png', 'cbc_plotsummary_3_rho_H2_vs_L1.png']
    page.image_table(image_list,webserver)

    page.section("Full Data Ifar")
    image_list = ['cbc_plotsummary_4_nonplayground_count_vs_ifar_H1H2L1.png','cbc_plotsummary_4_nonplayground_count_vs_ifar_H1L1.png','cbc_plotsummary_4_nonplayground_count_vs_ifar_H2L1.png']
    page.image_table(image_list,webserver)

    page.section("Full Data SNR")
    image_list = ['cbc_plotsummary_4_nonplayground_count_vs_snr_H1H2L1.png','cbc_plotsummary_4_nonplayground_count_vs_snr_H1L1.png','cbc_plotsummary_4_nonplayground_count_vs_snr_H2L1.png']
    page.image_table(image_list,webserver)
except: pass

page.finish()
