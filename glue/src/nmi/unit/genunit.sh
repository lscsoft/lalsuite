#!/bin/sh

# TODO: this ought to look in a lib dir under the glue install
TEMPLATE_DIR=$LSCSOFT_SRCDIR/lalsuite/glue/src/nmi/unit

#source $GLUE_LOCATION/etc/glue-user-env.sh
#source /opt/lscsoft/glue/etc/glue-user-env.sh

TMPDIR=${TMPDIR:-/tmp}

# exit immediately if any command exits with a non-zero status.
set -e
# treat unset variables as an error when performing parameter expansion.
set -u

USAGE="usage: $0 xml_file_to_reproduce lalapps_build_gid"
if [[ $# -lt 2 ]]; then
    echo $USAGE
    exit 1
fi

if [[ ! -f $1 ]]; then
    echo "error: $1 not found"
    echo $USAGE
    exit 1
fi

# create a temp dir for the submit files
SUBMIT_DIR=$TMPDIR/$USER/$(basename $0).${_NMI_GIT_ID}.$$.$(date +%s)
mkdir -p $SUBMIT_DIR
cd $SUBMIT_DIR

# make a copy of Metronome run-specification (aka submit) and
# input-specification files
# TODO: these should come from the current lscsoft install
# (e.g., $GLUE_PREFIX/lib/nmi) instead
cp $TEMPLATE_DIR/{cmdfile,scripts.git,remote_*.sh} .

# make sure given build GID exists
#nmi_gid2runid $2 > /dev/null

export _NMI_XML_FILE=$1
export _NMI_XML_BASENAME=$(basename $_NMI_XML_FILE)

# implicit parameters (should be set earlier in the workflow)
echo _NMI_GIT_ID=$_NMI_GIT_ID > /dev/null
echo _NMI_GIT_BRANCH=$_NMI_GIT_BRANCH > /dev/null
export _NMI_GIT_ID _NMI_GIT_BRANCH

# extract process (executable) from the reference XML file
NAME=$(ligolw_print -t process -c program $_NMI_XML_FILE)
export _NMI_LAL_EXE=lalapps_$NAME

# extract shorter workflow node name from XML filename (Miron would disapprove)
# TODO: can/should this be extracted more reliably from the metadata inside?
export _NMI_WORKFLOW_NODENAME=$(echo $_NMI_XML_BASENAME | sed 's/[^-]*-//; s/-.*//')

# generate an NMI input specification file pointing to the given build gid
export _NMI_BUILD_SOURCE_FILE=build.$2.$$.nmi
cp -f $TEMPLATE_DIR/build.nmi.template $_NMI_BUILD_SOURCE_FILE
echo "input_runids = $2" >> $_NMI_BUILD_SOURCE_FILE

# print some handy debugging info
echo SUBMIT_DIR=$SUBMIT_DIR
echo LSCSOFT_SRCDIR=$LSCSOFT_SRCDIR
env | grep ^_NMI_ | sort

nmi_submit --no-wait cmdfile
