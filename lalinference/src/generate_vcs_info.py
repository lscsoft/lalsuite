# generate_vcs_info.py - determine git version info
#
# Based on generateGitID.sh by Reinhard Prix
#
# Copyright (C) 2009,2010, Adam Mercer <adam.mercer@ligo.org>
# Copyright (C) 2009,2010, Nickolas Fotopoulos <nvf@gravity.phys.uwm.edu>
# Copyright (C) 2008,2009, John T. Whelan <john.whelan@ligo.org>
# Copyright (C) 2008, Reinhard Prix <reinhard.ligo.org>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or (at
# your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

#
# preamble
#

# metadata
__author__ = 'Adam Mercer <adam.mercer@ligo.org>'

# import required system modules
import exceptions
import os
import sys
import time
import optparse
import filecmp
import shutil
try:
  import subprocess
except ImportError:
  sys.exit("Python-2.4, or higher is required")

#
# class definitions
#

# version info class
class git_info(object):
  def __init__(self):
    id = None
    date = None
    branch = None
    tag = None
    author = None
    committer = None
    status = None

# git invocation error exception handler
class GitInvocationError(exceptions.LookupError):
  pass

#
# helper methods
#

# command line option parsing
def parse_args():
  # usage information
  usage = '%prog [options] project src_dir build_dir'

  # define option parser
  parser = optparse.OptionParser(usage=usage)
  parser.add_option('--sed', action='store_true', default=False,
      help='output sed commands for version replacement')
  parser.add_option('--sed-file', action='store_true', default=False,
      help='output sed file for version replacement')

  # parse command line options
  opts, args = parser.parse_args()

  # check for positional arguments
  if (len(args) != 3):
    parser.error('incorrect number of command line options specified')

  # check that both --sed and --file are not specified
  if opts.sed and opts.sed_file:
    parser.error('cannot specify both --sed and --sed-file')

  return opts, args[0], args[1], args[2]

#
# process management methods
#

# return output from running given command
def call_out(command):
  """
  Run the given command (with shell=False) and return a tuple of
  (int returncode, str output). Strip the output of enclosing whitespace.
  """
  # start external command process
  p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

  # get outputs
  out, _ = p.communicate()

  return p.returncode, out.strip()

def check_call_out(command):
  """
  Run the given command (with shell=False) and return the output as a
  string. Strip the output of enclosing whitespace.
  If the return code is non-zero, throw GitInvocationError.
  """
  # start external command process
  p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

  # get outputs
  out, _ = p.communicate()

  # throw exception if process failed
  if p.returncode != 0:
    raise GitInvocationError, 'failed to run "%s"' % " ".join(command)

  return out.strip()

#
# git version method
#

def in_git_repository():
  """
  Return True if git is available and we are in a git repository; else
  return False.

  NB: Unfortunately there is a magic number without any documentation to back
  it up. It turns out that git status returns non-zero exit codes for all sorts
  of success conditions, but I cannot find any documentation of them. 128 was
  determined empirically. I sure hope that it's portable.
  """
  ret_code, git_path = call_out(('which', 'git'))
  return (ret_code != 1) and not git_path.startswith('no') \
      and (call_out((git_path, 'status'))[0] != 128)

# generate sed command file containing version info
def generate_git_version_info():
  # info object
  info = git_info()

  git_path = check_call_out(('which', 'git'))

  # determine basic info about the commit
  # %H -- full git hash id
  # %ct -- commit time
  # %an, %ae -- author name, email
  # %cn, %ce -- committer name, email
  git_id, git_udate, git_author_name, git_author_email, \
    git_committer_name, git_committer_email = \
    check_call_out((git_path, 'log', '-1',
    '--pretty=format:%H,%ct,%an,%ae,%cn,%ce')).split(",")

  git_date = time.strftime('%Y-%m-%d %H:%M:%S +0000',
    time.gmtime(float(git_udate)))
  git_author = '%s <%s>' % (git_author_name, git_author_email)
  git_committer = '%s <%s>' % (git_committer_name, git_committer_email)

  # determine branch
  branch_match = check_call_out((git_path, 'rev-parse',
    '--symbolic-full-name', 'HEAD'))
  if branch_match == "HEAD":
    git_branch = None
  else:
    git_branch = os.path.basename(branch_match)

  # determine tag
  status, git_tag = call_out((git_path, 'describe', '--exact-match',
    '--tags', git_id))
  if status != 0:
    git_tag = None

  # refresh index
  check_call_out((git_path, 'update-index', '-q', '--refresh'))

  # check working copy for changes
  status_output = subprocess.call((git_path, 'diff-files', '--quiet'))
  if status_output != 0:
    git_status = 'UNCLEAN: Modified working tree'
  else:
    # check index for changes
    status_output = subprocess.call((git_path, 'diff-index', '--cached',
      '--quiet', 'HEAD'))
    if status_output != 0:
      git_status = 'UNCLEAN: Modified index'
    else:
      git_status = 'CLEAN: All modifications committed'

  # determine version strings
  info.id = git_id
  info.date = git_date
  info.branch = git_branch
  info.tag = git_tag
  info.author = git_author
  info.committer = git_committer
  info.status = git_status

  return info

#
# main program entry point
#

if __name__ == "__main__":
  # parse command line options
  options, project, src_dir, build_dir = parse_args()

  # filenames
  basename = '%sVCSInfo.h' % project
  infile = '%s.in' % basename  # used after chdir to src_dir
  tmpfile= '%s/%s.tmp' % (build_dir, basename)
  srcfile = '%s/%s' % (src_dir, basename)
  dstfile = '%s/%s' % (build_dir, basename)

  # copy vcs header to build_dir, if appropriate
  if os.access(srcfile, os.F_OK) and not os.access(dstfile, os.F_OK):
    shutil.copy(srcfile, dstfile)

  # change to src_dir
  os.chdir(src_dir)

  # generate version information output, if appropriate
  # NB: First assume that git works and we are in a repository since
  # checking it is expensive and rarely necessary.
  try:
    info = generate_git_version_info()
  except GitInvocationError:
    # We expect a failure if either git is not available or we are not
    # in a repository. Any other failure is unexpected and should throw an
    # error.
    if not in_git_repository():
        sys.exit(0)
    else:
        sys.exit("Unexpected failure in discovering the git version")

  if options.sed_file:
    # output sed command file to stdout
    print 's/@ID@/%s/g' % info.id
    print 's/@DATE@/%s/g' % info.date
    print 's/@BRANCH@/%s/g' % info.branch
    print 's/@TAG@/%s/g' % info.tag
    print 's/@AUTHOR@/%s/g' % info.author
    print 's/@COMMITTER@/%s/g' % info.committer
    print 's/@STATUS@/%s/g' % info.status
  elif options.sed:
    # generate sed command line options
    sed_cmd = ('sed',
               '-e', 's/@ID@/%s/' % info.id,
               '-e', 's/@DATE@/%s/' % info.date,
               '-e', 's/@BRANCH@/%s/' % info.branch,
               '-e', 's/@TAG@/%s/' % info.tag,
               '-e', 's/@AUTHOR@/%s/' % info.author,
               '-e', 's/@COMMITTER@/%s/' % info.committer,
               '-e', 's/@STATUS@/%s/' % info.status,
               infile)

    # create tmp file
    # FIXME: subprocess.check_call becomes available in Python 2.5
    sed_retcode = subprocess.call(sed_cmd, stdout=open(tmpfile, "w"))
    if sed_retcode:
      raise GitInvocationError, "Failed call (modulo quoting): " \
          + " ".join(sed_cmd) + " > " + tmpfile

    # only update vcs header if appropriate
    if os.access(dstfile, os.F_OK) and filecmp.cmp(dstfile, tmpfile):
      os.remove(tmpfile)
    else:
      os.rename(tmpfile, dstfile)
  else:
    # output version info
    print 'Id: %s' % info.id
    print 'Date: %s' % info.date
    print 'Branch: %s' % info.branch
    print 'Tag: %s' % info.tag
    print 'Author: %s' % info.author
    print 'Committer: %s' % info.committer
    print 'Status: %s' % info.status

# vim: syntax=python tw=72 ts=2 et
