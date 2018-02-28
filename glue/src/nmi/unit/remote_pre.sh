#!/bin/sh

# this task has been separated from the remote_task (execute.sh) in
# order to insulate that task's lalapps code from the production glue
# environment we use for ligolw_print here

# 2013-04-17 pfcouvar: remove since we now use a standard OS install of lalsuite
#. /opt/lscsoft/glue/etc/glue-user-env.sh

# exit immediately if any command exits with a non-zero status.
set -e
# treat unset variables as an error when performing parameter expansion.
set -u
# print (unexpanded) shell input lines as they are read
set -v
# print (expanded) commands before they're executed
set -x

#debugging
find .

# since we specified "untar_results = symlink" in the nmi input file, we need to untar it here
#tar zxvf results.tar.gz

# extract params from the reference XML file
# NOTE: the sed command strips out single-quote chars surrounding
# filename values (which can occur in, e.g., coire's XML output)
ligolw_print -t process_params -c param -c value -d " " $NMI_ligo_reference_xml | tr \\n " " | sed s/"'"//g > params.tmp

# extract input files from the reference XML file
# NOTE: the sed command strips out single-quote chars surrounding
# filename values (which can occur in, e.g., coire's XML output)
ligolw_print -t search_summvars -c name -c string -d " " $NMI_ligo_reference_xml | grep ^input_file | cut -f 2 -d " " | sed s/"'"//g > input_files.tmp

# if --bank-file or --veto-file params are defined,
# add their values (XML filenames) to our symlink list
# NOTE: the sed command strips out single-quote chars surrounding
# the filename value (which can occur in, e.g., coire's XML output)
ligolw_print -t process_params -c param -c value -d " " $NMI_ligo_reference_xml | egrep '^--..-veto-file|^--bank-file|^--injection-file|^--glob' | cut -f 2 -d " " | sed s/"'"//g > symlinks.tmp

# extract process (executable) from the reference XML file
#echo -n lalapps_ > exe.tmp
#ligolw_print -t process -c program $NMI_ligo_reference_xml >> exe.tmp

