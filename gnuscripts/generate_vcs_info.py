# generate_vcs_info.py - determine git version info
#
# Based on generateGitID.sh by Reinhard Prix
#
# Copyright (C) 2012-2013 Karl Wette <karl.wette@ligo.org>
# Copyright (C) 2009-2014, Adam Mercer <adam.mercer@ligo.org>
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
import re

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

#
# process management methods
#

def _call_out(work_dir, command):

  # start external command process
  p = subprocess.Popen(command, cwd=work_dir, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

  # get outputs
  out, _ = p.communicate()

  # convert byte objects to strings, if appropriate
  if not isinstance(out, str) and sys.version_info >= (3,):
    out = str(out, encoding='utf8')

  return p.returncode, out.strip()

def _check_call_out(work_dir, command):

  # call _call_out
  returncode, out = _call_out(work_dir, command)

  # throw exception if process failed
  if returncode != 0:
    raise subprocess.CalledProcessError(returncode, command)

  return out.strip()

#
# git version methods
#

def get_vcs_info(repo_dir, git_path='git'):
  """
  Returns an instance of class git_info() containing version info
  derived from repository in 'repo_dir' by running 'git_path'.
  """

  # info object
  info = git_info()

  # determine basic info about the commit
  # %H -- full git hash id
  # %ct -- commit time
  # %an, %ae -- author name, email
  # %cn, %ce -- committer name, email
  git_id, git_udate, git_author_name, git_author_email, \
    git_committer_name, git_committer_email = \
    _check_call_out(repo_dir, (git_path, 'log', '-1',
      '--pretty=format:%H,%ct,%an,%ae,%cn,%ce')).split(",")

  git_date = time.strftime('%Y-%m-%d %H:%M:%S +0000',
    time.gmtime(float(git_udate)))
  git_author = '%s <%s>' % (git_author_name, git_author_email)
  git_committer = '%s <%s>' % (git_committer_name, git_committer_email)

  # determine branch
  branch_match = _check_call_out(repo_dir, (git_path, 'rev-parse',
    '--symbolic-full-name', 'HEAD'))
  if branch_match == "HEAD":
    git_branch = None
  else:
    git_branch = os.path.basename(branch_match)

  # determine tag
  status, git_tag = _call_out(repo_dir, (git_path, 'describe', '--exact-match',
    '--tags', git_id))
  if status != 0:
    git_tag = None

  # refresh index
  _check_call_out(repo_dir, (git_path, 'update-index', '-q', '--refresh'))

  # check working copy for changes
  status, output = _call_out(repo_dir, (git_path, 'diff-files', '--quiet'))
  if status != 0:
    git_status = 'UNCLEAN: Modified working tree'
  else:
    # check index for changes
    status, output = _call_out(repo_dir, (git_path, 'diff-index', '--cached',
      '--quiet', 'HEAD'))
    if status != 0:
      git_status = 'UNCLEAN: Modified index'
    else:
      git_status = 'CLEAN: All modifications committed'

  # get builder
  retcode, builder_name = _call_out(repo_dir, (git_path, 'config', 'user.name'))
  if retcode:
    builder_name = 'Unknown User'
  retcode, builder_email = _call_out(repo_dir, (git_path, 'config', 'user.email'))
  if retcode:
    builder_email = ''
  git_builder = '%s <%s>' % (builder_name, builder_email)

  # determine version strings
  info.id = git_id
  info.date = git_date
  info.branch = git_branch
  info.tag = git_tag
  info.author = git_author
  info.committer = git_committer
  info.status = git_status
  info.builder = git_builder

  return info

def generate_vcs_info(dst_file, src_file, git_path='git'):
  """
  Generates a version info file 'dst_file' from a template 'src_file',
  by running 'git_path' on the repository containing 'src_file'.
  Returns True if 'dst_file' was modified, False otherwise.
  """

  # get git repository directory
  repo_dir = os.path.dirname(src_file)
  if len(repo_dir) == 0:
    repo_dir = '.'

  # create temporary destination filename
  dst_tmpfile = dst_file + '.tmp'

  # generate version information output, if appropriate
  # NB: assume that git works and we are in a repository
  info = get_vcs_info(repo_dir, git_path)

  # construct dictionary of replacements
  replacements = {
    '@ID@' : info.id,
    '@DATE@' : info.date,
    '@BRANCH@' : info.branch,
    '@TAG@' : info.tag,
    '@AUTHOR@' : info.author,
    '@COMMITTER@' : info.committer,
    '@STATUS@' : info.status,
    '@BUILDER@' : info.builder,
    }

  # open source temporary destination file
  src = open(src_file, 'r')
  dst = open(dst_tmpfile, 'w')

  # read source file, perform replacements, write to temporary destination file
  for line in src:
    for repl in replacements:
      line = line.replace(repl, str(replacements[repl]))
    dst.write(line)

  # close files
  src.close()
  dst.close()

  # update if appropriate
  if os.access(dst_file, os.F_OK) and filecmp.cmp(dst_file, dst_tmpfile):
    os.remove(dst_tmpfile)
    modified = False
  else:
    os.rename(dst_tmpfile, dst_file)
    modified = True

  return modified

#
# main program entry point
#

if __name__ == "__main__":

  # usage information
  usage = '%prog [--git-path=/path/to/git] [--am-v-gen=$(AM_V_GEN)] output_file source_file'

  # define option parser
  parser = optparse.OptionParser(usage=usage)
  parser.add_option("--git-path", dest="git_path", default="git")
  parser.add_option("--am-v-gen", dest="am_v_gen", default="")

  # parse command line options
  opts, args = parser.parse_args()
  git_path = opts.git_path
  am_v_gen = opts.am_v_gen

  # check for positional arguments
  if (len(args) != 2):
    parser.error('incorrect number of command line options specified')
  dst_file = args[0]
  src_file = args[1]

  # generate version info
  modified = generate_vcs_info(dst_file, src_file, git_path)

  # print generation message
  if modified and len(am_v_gen) > 0:
    prefix = re.findall("[ ]*GEN[ ]*", am_v_gen)
    print('%s %s' % (prefix[0], dst_file))

# vim: syntax=python tw=72 ts=2 et
