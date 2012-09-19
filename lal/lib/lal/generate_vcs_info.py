# generate_vcs_info.py - determine git version info
#
# Based on generateGitID.sh by Reinhard Prix
#
# Copyright (C) 2009,2010,2012, Adam Mercer <adam.mercer@ligo.org>
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
import os
import sys
import time
import optparse
import filecmp
import shutil
import subprocess

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
    builder = None
    build_date = None

# git invocation error exception handler
class GitInvocationError(LookupError):
  pass

#
# helper methods
#

# command line option parsing
def parse_args():
  # usage information
  usage = '%prog project src_dir build_dir'

  # define option parser
  parser = optparse.OptionParser(usage=usage)

  # parse command line options
  opts, args = parser.parse_args()

  # check for positional arguments
  if (len(args) != 3):
    parser.error('incorrect number of command line options specified')

  return args[0], args[1], args[2]

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
    raise GitInvocationError('failed to run "%s"' % " ".join(command))

  # convert byte objects to strings, if appropriate
  if not isinstance(out, str) and sys.version_info >= (3,):
    out = str(out, encoding='utf8')

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

  # get builder
  retcode, builder_name = call_out((git_path, 'config', 'user.name'))
  if retcode:
    builder_name = 'Unknown User'
  retcode, builder_email = call_out((git_path, 'config', 'user.email'))
  if retcode:
    builder_email = ''
  git_builder = '%s <%s>' % (builder_name, builder_email)

  # determine current time and treat it as the build time
  git_build_date = time.strftime('%Y-%m-%d %H:%M:%S +0000', time.gmtime())

  # determine version strings
  info.id = git_id
  info.date = git_date
  info.branch = git_branch
  info.tag = git_tag
  info.author = git_author
  info.committer = git_committer
  info.status = git_status
  info.builder = git_builder
  info.build_date = git_build_date

  return info

#
# main program entry point
#

if __name__ == "__main__":
  # parse command line options
  project, src_dir, build_dir = parse_args()

  # filenames, header
  header_basename = '%sVCSInfo.h' % project
  header_infile = '%s.in' % header_basename  # used after chdir to src_dir
  header_tmpfile= '%s/%s.tmp' % (build_dir, header_basename)
  header_srcfile = '%s/%s' % (src_dir, header_basename)
  header_dstfile = '%s/%s' % (build_dir, header_basename)

  # filename, module
  if project == 'LALApps':
    module_basename = 'git_version.py'
    module_infile = '%s.in' % module_basename
    module_tmpfile = '%s/%s.tmp' % (build_dir, module_infile)
    module_srcfile = '%s/%s' % (src_dir, module_basename)
    module_dstfile = '%s/%s' % (build_dir, module_basename)

  # copy vcs header to build_dir, if appropriate
  if os.access(header_srcfile, os.F_OK) and not os.access(header_dstfile, os.F_OK):
    shutil.copy(header_srcfile, header_dstfile)

  # copy vcs module to build_dir, if appropriate
  if project == 'LALApps':
    if os.access(module_srcfile, os.F_OK) and not os.access(module_dstfile, os.F_OK):
      shutil.copy(module_srcfile, module_dstfile)

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

  # construct sed command options
  sed_opts = ['sed',
             '-e', 's/@ID@/%s/' % info.id,
             '-e', 's/@DATE@/%s/' % info.date,
             '-e', 's/@BRANCH@/%s/' % info.branch,
             '-e', 's/@TAG@/%s/' % info.tag,
             '-e', 's/@AUTHOR@/%s/' % info.author,
             '-e', 's/@COMMITTER@/%s/' % info.committer,
             '-e', 's/@STATUS@/%s/' % info.status,
             '-e', 's/@BUILDER@/%s/' % info.builder,
             '-e', 's/@BUILD_DATE@/%s/' % info.build_date]

  # construct sed command for vcs header
  header_sed_cmd = list(sed_opts)
  header_sed_cmd.append(header_infile)

  # construct sed command for vcs module
  if project == 'LALApps':
    module_sed_cmd = list(sed_opts)
    module_sed_cmd.append(module_infile)

  # create tmp files
  try:
    subprocess.check_call(header_sed_cmd, stdout=open(header_tmpfile, 'w'))
    if project == 'LALApps':
      subprocess.check_call(module_sed_cmd, stdout=open(module_tmpfile, 'w'))
  except CalledProcessError:
    raise GitInvocationError("Failed call (modulo quoting): " \
        + " ".join(sed_cmd) + " > " + tmpfile)

  # update if appropriate
  if os.access(header_dstfile, os.F_OK) and filecmp.cmp(header_dstfile, header_tmpfile):
    os.remove(header_tmpfile)
    if project == 'LALApps':
      os.remove(module_tmpfile)
  else:
    print('  GEN    %s' % header_basename)
    os.rename(header_tmpfile, header_dstfile)
    if project == 'LALApps':
      print('  GEN    %s' % module_basename)
      os.rename(module_tmpfile, module_dstfile)

# vim: syntax=python tw=72 ts=2 et
