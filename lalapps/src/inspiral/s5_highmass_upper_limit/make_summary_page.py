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

  def image_link(self,path):
    thumb = "thumb_" + path
    command = 'convert ' + path + ' -resize 400x400 ' + thumb
    print command
    popen = subprocess.Popen(command.split())
    popen.communicate()
    status = popen.returncode
    s = '[[ImageLink('+path+','+thumb+',width=400][,alt=none])]]'
    self.file.write(s)

  def image_table(self,image_list):
    if not image_list: return
    self.file.write("||")
    for i in image_list:
      self.image_link(i)
      self.file.write("||")
    self.file.write("\n")

  def section(self,title):
    s = "=== "+title.strip()+" ===\n"
    self.file.write(s)

  def finish(self):
    self.file.close()  

page = wiki()

page.section("Found / Missed")
image_list = ['cbc_plotsummary_0_deff_vs_mchirp_H1H2L1.png','cbc_plotsummary_0_deff_vs_mchirp_H1L1.png','cbc_plotsummary_0_deff_vs_mchirp_H2L1.png']
page.image_table(image_list)
image_list = ['cbc_plotsummary_0_deff_vs_t_H1H2L1.png','cbc_plotsummary_0_deff_vs_t_H1L1.png','cbc_plotsummary_0_deff_vs_t_H2L1.png']
page.image_table(image_list)




page.section("Playground Chi-squared")
image_list = ['cbc_plotsummary_2_playground_chi2_vs_rho_H1.png', 'cbc_plotsummary_2_playground_chi2_vs_rho_H2.png', 'cbc_plotsummary_2_playground_chi2_vs_rho_L1.png']
page.image_table(image_list)

page.section("Playground Effective SNR scatter")
image_list = ['cbc_plotsummary_3_playground_rho_H1_vs_L1.png', 'cbc_plotsummary_3_playground_rho_H1_vs_H2.png', 'cbc_plotsummary_3_playground_rho_H2_vs_L1.png']
page.image_table(image_list)

page.section("Playground Ifar")
image_list = ['cbc_plotsummary_4_playground_count_vs_ifar_H1H2L1.png','cbc_plotsummary_4_playground_count_vs_ifar_H1L1.png','cbc_plotsummary_4_playground_count_vs_ifar_H2L1.png']
page.image_table(image_list)

page.section("Playground SNR")
image_list = ['cbc_plotsummary_4_playground_count_vs_snr_H1H2L1.png','cbc_plotsummary_4_playground_count_vs_snr_H1L1.png','cbc_plotsummary_4_playground_count_vs_snr_H2L1.png']
page.image_table(image_list)

if sys.argv[1] == "open":

  print >>sys.stderr, "WARNING: OPENING THE BOX"

  page.section("Full Data Chi-squared")
  image_list = ['cbc_plotsummary_2_chi2_vs_rho_H1.png', 'cbc_plotsummary_2_chi2_vs_rho_H2.png', 'cbc_plotsummary_2_chi2_vs_rho_L1.png']
  page.image_table(image_list)

  page.section("Full Data Effective SNR scatter")
  image_list = ['cbc_plotsummary_3_rho_H1_vs_L1.png', 'cbc_plotsummary_3_rho_H1_vs_H2.png', 'cbc_plotsummary_3_rho_H2_vs_L1.png']
  page.image_table(image_list)

  page.section("Full Data Ifar")
  image_list = ['cbc_plotsummary_4_nonplayground_count_vs_ifar_H1H2L1.png','cbc_plotsummary_4_nonplayground_count_vs_ifar_H1L1.png','cbc_plotsummary_4_nonplayground_count_vs_ifar_H2L1.png']
  page.image_table(image_list)

  page.section("Full Data SNR")
  image_list = ['cbc_plotsummary_4_nonplayground_count_vs_snr_H1H2L1.png','cbc_plotsummary_4_nonplayground_count_vs_snr_H1L1.png','cbc_plotsummary_4_nonplayground_count_vs_snr_H2L1.png']
  page.image_table(image_list)




page.finish()
