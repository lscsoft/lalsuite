#!/usr/bin/env python
#
# generates and submits a Metronome build of lalsuite from a given git URL and id

import os
import sys
import commands

from optparse import OptionParser


# define common options (e.g., -v -q -h) 
parser = OptionParser(version="%prog $Id$")
parser.add_option("-v", "--verbose",
                  action="store_true", dest="verbose", default=False,
                  help="print verbose status messages to stdout")
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose",
                  help="don't print status messages to stdout")

# define custom options
parser.add_option("-r", "--git-repo", dest="git_repo",
                  default="git://ligo-vcs.phys.uwm.edu/lalsuite.git",
                  help="git repository", metavar="URL (git:// or file://)")
parser.add_option("-i", "--git-id", dest="git_id",
                  help="git id", metavar="ID")
parser.add_option("-b", "--git-branch", dest="git_branch",
                  help="git branch", metavar="NAME")
parser.add_option("--harness-git-repo", dest="harness_git_repo",
                  help="build harness git repository",
                  metavar="URL (git:// or file://)")
parser.add_option("--harness-git-id", dest="harness_git_id",
                  help="build harness git id", metavar="ID")
parser.add_option("--harness-git-branch", dest="harness_git_branch",
                  help="build harness git branch", metavar="NAME")
parser.add_option("--template-dir", dest="template_dir",
                  help="Metronome submit file template location", metavar="DIR")
(options, args) = parser.parse_args()

# copy options to local vars so we can modify them
git_repo = options.git_repo
git_id = options.git_id
git_branch = options.git_branch
harness_git_repo = options.harness_git_repo
harness_git_id = options.harness_git_id
harness_git_branch = options.harness_git_branch
template_dir = options.template_dir

# if unspecified, harness should use the same repo/branch/tag as code
if harness_git_repo is None:
    harness_git_repo = git_repo
if harness_git_id is None:
    harness_git_id = git_id
if harness_git_branch is None:
    harness_git_branch = git_branch

if template_dir is None:
    template_dir = "%s/share/nmi" % os.getenv("GLUE_PREFIX")

# sanity-check our options
if git_id is None:
    print >> sys.stderr, "ERROR: git_id (-i) required"
    sys.exit(2)
if git_branch is None:
    print >> sys.stderr, 'WARNING: no git_branch specified, using "unknown_branch" in metadata'
    git_branch = "unknown_branch"

# TODO: confirm that the given git id exists on the given git branch
# (currently the git branch is unused except as informational
# metadata, so user error can lead to wrong and confusing metadata).
# See <http://stackoverflow.com/questions/2444100> for how to do this.

# propagate relevant options to Metronome via special environment vars
# (this allows us to use a single static Metronome submit file
# template, located in glue/lib/nmi, for all lalsuite builds)
os.environ["_NMI_GIT_REPO"] = git_repo
os.environ["_NMI_GIT_ID"] = git_id
os.environ["_NMI_GIT_BRANCH"] = git_branch
os.environ["_NMI_HARNESS_GIT_REPO"] = harness_git_repo
os.environ["_NMI_HARNESS_GIT_ID"] = harness_git_id
os.environ["_NMI_HARNESS_GIT_BRANCH"] = harness_git_branch
os.environ["_NMI_HARNESS_INPUT_SPEC_FILE"] = "%s/lalsuite-build-scripts.location" % template_dir

if options.verbose: os.system("env|fgrep _NMI|sort")

# we need to tell Metronome how & where to find the build harness code
# depending on whether it was specified as a git:// or file:// URL
if 'git://' not in harness_git_repo and 'file://' not in harness_git_repo:
    print >> sys.stderr, "ERROR: invalid harness repository URL (\"%s\"): only git:// and file:// schemes are currently supported" % harness_git_repo
    sys.exit(2)

(status, output) = commands.getstatusoutput("nmi_submit %s/lalsuite-build.spec" % template_dir)

print output
sys.exit(status)
