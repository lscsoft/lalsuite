#!/bin/sh

# A == reference (presumed good) unit test
# B == candidate unit test

# exit immediately if any command exits with a non-zero status.
set -e
# treat unset variables as an error when performing parameter expansion.
set -u

# implicit parameters (should be set earlier in the workflow)
export _NMI_WORKFLOW_NODENAME
export _NMI_GIT_BRANCH_A
export _NMI_GIT_BRANCH_B
export _NMI_GIT_ID_A
export _NMI_GIT_ID_B

export _NMI_HARNESS_GIT_REPO=${_NMI_HARNESS_GIT_REPO:="git://ligo-vcs.phys.uwm.edu/lalsuite.git"}
export _NMI_HARNESS_GIT_BRANCH=${_NMI_HARNESS_GIT_BRANCH:="master"}

if [[ $# -lt 2 ]]; then
    echo usage: $0 runid_a runid_b
    exit 1
fi

export _NMI_RUNID_A=$1
export _NMI_RUNID_B=$2

# assumed to be defined already

echo
env | grep ^_NMI_ | sort
echo

nmi_submit --no-wait cmdfile
