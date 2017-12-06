#!/usr/bin/env python

"""
Script to be setup and launch a cron job for the automated running of the known pulsar analysis.

The script requires an initial configuration file. Each time the script will be re-run by cron it will itself
automatically update the times in the configuration file.

(C) Matthew Pitkin 2015
"""

# make print statements python 3-proof
from __future__ import print_function, division

import sys
import os
import ast
import calendar
import numpy as np
import subprocess as sp
import socket
import smtplib

import argparse
from ConfigParser import ConfigParser

from pylal import git_version

# try importing astropy
try:
  from astropy.time import Time
except ImportError:
  print("Could not import astropy: make sure astropy is installed (e.g. 'pip install astropy') and in the PYTHONPATH", file=sys.stderr)
  sys.exit(1)


# try import python-crontab
try:
  from crontab import CronTab
except ImportError:
  print("Could not import python-crontab: make sure it is installed (e.g. 'pip install python-crontab') and in the PYTHONPATH", file=sys.stderr)
  sys.exit(1)


__author__ = "Matthew Pitkin <matthew.pitkin@ligo.org>"
__version__ = "git id %s"%git_version.id
__date__ = git_version.date


# number of seconds for each allowed time step value
tsdic = {'hourly': 1440., 'daily': 86400., 'weekly': 7.*86400., 'monthly': 86400.*calendar.monthrange(Time.now().datetime.year, Time.now().datetime.month)[1]}


# function to remove crontab job
def remove_cron(cronid):
  try:
    cron = CronTab(user=True)
    cron.remove_all(comment=cronid)
    cron.write()
  except:
    print("Warning... count not remove crontab job with comment ID '%s'" % cronid)

  return


# function to send an email
def send_email(FROM, TO, SUBJECT, MESSAGE, server, quitserver=True):
  # set a notication email template
  emailtemplate = "From: {0}\nTo: {1}\nSubject: {3}\n\n{4}"
  message = emailtemplate.format(FROM, TO, SUBJECT, MESSAGE)

  try:
    server.sendmail(FROM, TO, message)
  except:
    print("Warning... unable to send email.")

  if quitserver:
    server.quit()


# main function
if __name__=='__main__':
  description = """This script will setup and run a cron job for an automated known pulsar pipeline.
A configuration .ini file is required.
  """

  parser = argparse.ArgumentParser( description = description, version = __version__ )
  parser.add_argument("inifile", help="The configuration (.ini) file")

  # parse input options
  opts = parser.parse_args()

  inifile = opts.inifile

  startcron = False # variable to say whether to create the crontab (if this is the first time the script is run then this will be changed to True later)
  cronid = 'knopeJob' # an ID for the crontab job
  croncommand = '%s -i {0}' % sys.argv[0] # set the cron command (which will re-run this script)

  # open and parse config file
  cp = ConfigParser()
  try:
    cp.read(inifile)
  except:
    print("Error... cannot parse configuration file '%s'" % inifile, file=sys.stderr)
    sys.exit(1)

  # if configuration file has previous_endtimes option then the cronjob must have started
  if not cp.has_option('times', 'previous_endtimes'): # make sure to start the crontab job
    startcron = True

  # open and parse the run configuation file
  if cp.has_option('configuration', 'file'):
    runconfig = cp.get('configuration', 'file') # Get main configuration ini template for the run

    if not os.path.isfile(runconfig):
      print("Error... run configuration file '%s' does not exist!" % runconfig, file=sys.stderr)
      if not startcron: remove_cron(cronid) # remove cron job
      sys.exit(1)
  else:
    print("Error... must specify a run configuration '.ini' file", file=sys.stderr)
    if startcron: remove_cron(cronid)
    sys.exit(1)

  cprun = ConfigParser()
  try:
    cprun.read(runconfig)
  except:
    print("Error... could not read run configuration '.ini' file", file=sys.stderr)
    sys.exit(1)

  # try getting email address information (for notifications of e.g. job failures, re-running rescue DAGs, etc)
  email = None
  server = None
  if cp.has_option('notification', 'email'):
    try:
      email = cp.get('notification', 'email')
    except:
      print("Warning... could not get email address from configuration file. No notifications will be sent.")

    if '@' not in email:
      print("Warning... '%s' is not a valid email address. No notifications will be sent.")
      email = None

  # set email server
  if email != None:
    try:
      server = smtplib.SMTP('localhost')
    except:
      print("Warning... could not get SMTP server. No notication emails will be sent.")
      email = None

  # set email sender (if this fails maybe just hardcode sender to be 'matthew.pitkin@ligo.org')
  FROM = None
  if email != None:
    try:
      HOST = socket.getfqdn()
      USER = os.environ['USER']
      FROM = USER+'@'+HOST
    except:
      FROM = 'matthew.pitkin@ligo.org'

  # Get the current time
  now = Time.now()
  gpsnow = int(now.replicate('gps').value)

  # check for DAG completion in previous analyses
  prevdags = None
  try:
    rundir = cprun.get('analysis', 'run_dir') # get run directory where DAG will have been created
  except:
    rundir = os.getcwd() # if 'run_dir' was not set then the current working directory will have been used
  if not startcron:
    try:
      prevdags = ast.literal_eval(cp.get('configuration', 'previous_dags'))
    except:
      errmsg = "Error... no previous DAG file(s) have been set."
      print(errmsg, file=sys.stderr)
      remove_cron(cronid) # remove cron job
      if email != None:
        subject = sys.argv[0] + ': Error message'
        send_email(FROM, email, subject, errmsg, server)
      sys.exit(1)

    # check if the last dag has completed yet
    lockfile = prevdags[-1] + '.lock'
    if os.path.isfile(lockfile):
      # DAG has not yet complete, so wait before trying again (by updating the crontab entry)
      try:
        # reset to run again later
        if timestep in ['hourly', 'daily']: # if hourly or daily just wait until the next run
          print("Previous DAG not finished. Re-running later")
          sys.exit(0)
        else: # add a day to the crontab job and re-run then
          cron = CronTab(user=True)
          for job in cron.find_comment(cronid)
            thisjob = job # cron job

          # get a detlaT for a day
          t1 = Time('2010-01-01 00:00:00')
          t2 = Time('2010-01-02 00:00:00')
          dt = t2-t1
          newcrontime = now + dt
          if timestep == 'weekly':
            thisjob.dow.on(newcrontime.datetime.strftime('%a').upper()) # get new day of the week
          else:
            thisjob.day.on(newcrontime.datetime.day) # get new month of the year
          cron.write()
          sys.exit(0)
      except:
        errmsg = "Error... could not reset the crontab to wait for DAG completion."
        print(errmsg, file=sys.stderr)
        remove_cron(cronid)
        if email != None:
          subject = sys.argv[0] + ': Error message'
          send_email(FROM, email, subject, errmsg, server)
        sys.exit(1)

    # check if there is a rescue DAG, and if so run it, and wait
    rescuefile = prevdags[-1] + '.rescue'
    if rescuefile in os.listdir(os.path.dirname(prevdags[-1])):
      # if more than 2 rescue DAGs have already been run then just abort as there's probably some problem
      if rescuefile+'003' in os.listdir(os.path.dirname(prevdags[-1])):
        errmsg = "Error... rescue DAG has been run twice and there are still failures. Automation code is aborting. Fix the problem and then retry"
        print(errmsg, file=sys.stderr)
        remove_cron(cronid) # remove cron job
        if email != None:
          subject = sys.argv[0] + ': Error message'
          send_email(FROM, email, subject, errmsg, server)
        sys.exit(1)

      # run rescue DAG
      from subprocess import Popen
      x = Popen(['condor_submit_dag', prevdags[-1]])
      x.wait()
      if x.returncode != 0:
        errmsg = "Error... unable to submit rescue DAG for '%s'. Automation code is aborting." % prevdags[-1]
        remove_cron(cronid) # remove cron job
        if email != None:
          subject = sys.argv[0] + ': Error message'
          send_email(FROM, email, subject, errmsg, server)
        sys.exit(1)

      # wait until re-running
      try:
        # reset to run again later
        if timestep in ['hourly', 'daily']: # if hourly or daily just wait until the next run
          print("Previous DAG not finished. Re-running later")
          sys.exit(0)
        else: # add a day to the crontab job and re-run then
          cron = CronTab(user=True)
          for job in cron.find_comment(cronid)
            thisjob = job # cron job

          # get a detlaT for a day
          t1 = Time('2010-01-01 00:00:00')
          t2 = Time('2010-01-02 00:00:00')
          dt = t2-t1
          newcrontime = now + dt
          if timestep == 'weekly':
            thisjob.dow.on(newcrontime.datetime.strftime('%a').upper()) # get new day of the week
          else:
            thisjob.day.on(newcrontime.datetime.day) # get new month of the year
          cron.write()
          sys.exit(0)
      except:
        errmsg = "Error... could not reset the crontab to wait for rescue DAG completion."
        print(errmsg, file=sys.stderr)
        remove_cron(cronid)
        if email != None:
          subject = sys.argv[0] + ': Error message'
          send_email(FROM, email, subject, errmsg, server)
        sys.exit(1)

  # Get the start time of the automated analysis - if not present default to the current time
  if cp.has_option('times', 'starttime'):
    try:
      starttime = cp.getint('times', 'starttime')
    except:
      errmsg = "Error... could not parse 'starttime' in '[times]'. A start time is required."
      print(errmsg, file=sys.stderr)
      if email != None:
        subject = sys.argv[0] + ': Error message'
        send_email(FROM, email, subject, errmsg, server)
      if not startcron: remove_cron(cronid)
      sys.exit(1)

  # check start time is in the past
  if startime >= gpsnow:
    errmsg = "Error... start time (%f) must be in the past!" % starttime
    print(errmsg, file=sys.stderr)
    if email != None:
      subject = sys.argv[0] + ': Error message'
      send_email(FROM, email, subject, errmsg, server)
    if not startcron: remove_cron(cronid)
    sys.exit(1)

  # Get the end time of the automated analysis - if not present default to infinity (never stop!)
  if cp.has_option('times', 'endtime'):
    try:
      endtime = cp.getint('times', 'endtime')
    except:
      print("Warning... could not parse 'endtime' in '[times]'. Defaulting to Infinity.")
      endtime = np.inf
  else:
    # defaulting to infinity
    endtime = np.inf

  # check end time is after start time
  if endtime <= starttime:
    errmsg = "Error... start time is after end time!"
    print(errmsg, file=sys.stderr)
    if not startcron: remove_cron(cronid) # remove cron job
    if email != None:
      subject = sys.argv[0] + ': Error message'
      send_email(FROM, email, subject, errmsg, server)
    sys.exit(1)

  # get the time increment for re-running the pipeline on (can be "hourly", "daily", "weekly" or "monthly")
  if cp.has_option('times', 'steps'):
    timestep = cp.get('times', 'steps')

    if timestep not in ['hourly', 'daily', 'weekly', 'monthly']:
      errmsg = "Error... 'steps' value in '[times'] must be 'hourly', 'daily', 'weekly' or 'monthly'"
      print(errmsg, file=sys.stderr)
      if not startcron: remove_cron(cronid) # remove cron job
      if email != None:
        subject = sys.argv[0] + ': Error message'
        send_email(FROM, email, subject, errmsg, server)
      sys.exit(1)
  else:
    errmsg = "Error... must specify a time step 'steps' value in '[times]'"
    print(errmsg, file=sys.stderr)
    if not startcron: remove_cron(cronid) # remove cron job
    if email != None:
      subject = sys.argv[0] + ': Error message'
      send_email(FROM, email, subject, errmsg, server)
    sys.exit(1)

  # get a lag time (this will add a lag to gpsnow - if there is a lag between data creation and replication on the various sites a lag may be required, so that the data exists)
  if cp.has_option('times', 'lag'):
    try:
      timelag = gp.getint('times', 'lag')
    except:
      timelag = 0
  else:
    timelag = 0 # default to no lag

  # check if this is not the first run of the script
  if not startcron:
    try:
      prev_ends = ast.literal_eval(cp.get('times', 'previous_endtimes'))
    except:
      errmsg = "Error... cannot parse previous end times list"
      print(errmsg, file=sys.stderr)
      if startcron: remove_cron(cronid) # remove cron job
      if email != None:
        subject = sys.argv[0] + ': Error message'
        send_email(FROM, email, subject, errmsg, server)
      sys.exit(1)

    # update start time in the configuration file to the previous end time
    newstart = prev_ends[-1]

    # if new start time is after end time stop the cronjob and exit
    if newstart >= endtime:
      remove_cron(cronid)
      sys.exit(0) # end the script

    newend = newstart + int(tsdic[timestep])

    # check if new end time is past the overall end time
    if newend >= endtime: newend = endtime
    else:
      # check if the current time is later than the new end time
      if newend < gpsnow-timelag: newend = gpsnow-timelag # set end time to now

    prev_ends.append(newend)

    cp.set('times', 'previous_endtimes', '['+', '.join(str(z) for z in prev_ends)+']') # output as list
  else: # create previous end times
    newstart = starttime
    newend = newstart + int(tsdic[timestep])

    # check if the current time is later than the new end time
    if newend >= endtime: newend = endtime
    else:
      if newend < gpsnow-timelag: newend = gpsnow-timelag # set end time to now

    cp.set('times', 'previous_endtimes', '['+str(newend)+']')

  # Get the script for running the full pipeline
  if cp.has_option('configuration', 'exec'):
    runscript = cp.get('configuration', 'exec')

    if not (os.path.isfile(runscript) and os.access(runscript, os.X_OK)):
      errmsg = "Error... run script '%s' does not exist or is not executable" % runscript
      print(errmsg, file=sys.stderr)
      if email != None:
        subject = sys.argv[0] + ': Error message'
        send_email(FROM, email, subject, errmsg, server)
      sys.exit(1)
  else:
    errmsg = "Error... a run script executable 'exec' is required in the '[configuration]' section."
    print(errmsg, file=sys.stderr)
    if not startcron: remove_cron(cronid) # remove cron job
    if email != None:
      subject = sys.argv[0] + ': Error message'
      send_email(FROM, email, subject, errmsg, server)
    sys.exit(1)

  # edit start and end times for the main run configuration script
  if cprun.has_section('analysis'):
    cprun.set('analysis', 'starttime', str(newstart))  # set start time
    cprun.set('analysis', 'endtime', str(newend))      # set end time
    cprun.set('analysis', 'autonomous', 'True')        # set 'autonomous' to true
    cprun.set('analysis', 'submit_dag', 'True')        # set to make sure Condor DAG is submitted

    # create file name for DAG
    dagname = 'automated_run_%s-%s.dag' % (str(newstart), str(newend))
    cprun.set('analysis', 'dag_name', dagname) # add this dag file name to the automation code configuration script (to be used to check for DAG completion)

    if prevdags != None:
      # add on new DAG file to list
      prevdags.append(os.path.join(rundir, dagname))
      cp.set('configuation', 'previous_dags', '['+', '.join(z for z in prevdags)+']') # output as list
    else: # output DAG file to previous_dags list
      cp.set('configuation', 'previous_dags', '['+dagname+']')

    # add the initial start time
    cprun.set('analysis', 'autonomous_initial_start', str(starttime))

    # write updated parameters to the file
    fc = open(runconfig, 'w')
    cprun.write(fc)
    fc.close()
  else:
    errmsg = "Error... run configuration file '%s' has no '[analysis]' section!" % runconfig
    print(errmsg, file=sys.stderr)
    if not startcron: remove_cron(cronid)
    if email != None:
      subject = sys.argv[0] + ': Error message'
      send_email(FROM, email, subject, errmsg, server)
    sys.exit(1)

  # Write out updated configuration file
  fc = open(inifile, 'w')
  cp.write(fc)
  fc.close()

  ### RUN ANALYSIS SCRIPT ###
  p = sp.Popen('{0} -i {1}'.format(runscript, runconfig), shell=True)
  out, err = p.communicate()
  if p.returncode != 0:
    errmsg = "Error... problem running main script '%s'.: %s, %s" % (runscript, out, err)
    print(errmsg, file=sys.stderr)
    if not startcron: remove_cron(cronid)
    if email != None:
      subject = sys.argv[0] + ': Error message'
      send_email(FROM, email, subject, errmsg, server)
    sys.exit(1)
  ###########################

  # create crontab job
  if startcron:
    try:
      cron = CronTab(user=True)
      job = cron.new(command=croncommand.format(runconfig), comment=cronid)

      # set job time - this will start at the next time step (as we've just run the first step)
      day = now.datetime.day
      month = now.datetime.month
      year = now.datetime.year
      hour = now.datetime.hour
      minute = now.datetime.minute
      dow = now.datetime.strftime('%a').upper() # day of the week

      if timestep == 'hourly': # required for 'hourly'
        job.minute.on(minute)
      elif timestep == 'daily': # required for 'daily'
        job.minute.on(minute)
        job.hour.on(hour)
      elif timestep == 'weekly': # required for 'weekly'
        job.minute.on(minute)
        job.hour.on(hour)
        job.dow.on(dow)
      elif timestep == 'monthly': # required for 'monthly'
        job.minute.on(minute)
        job.hour.on(hour)
        job.day.on(day)
      else:
        print("Error... unrecognised 'timestep' option '%s'" % timestep, file=sys.stderr)
        sys.exit(1)

      cron.write()
    except:
      errmsg = "Error... could not create crontab job"
      print(errmsg, file=sys.stderr)
      if email != None:
        subject = sys.argv[0] + ': Error message'
        send_email(FROM, email, subject, errmsg, server)
      sys.exit(1)

  sys.exit(0)
