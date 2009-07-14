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
    command = 'convert ' + path + ' -resize 300x300 ' + thumb
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
page.section("Volume x time H1H2L1, H1L1, H2L1")
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
page.section("Combined upper limit")
image_list = ['combinedupper_limit.png', 'combinedposterior.png']

page.image_table(image_list,webserver)
page.finish()

