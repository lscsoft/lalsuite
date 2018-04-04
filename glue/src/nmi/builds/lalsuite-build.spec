# Metronome submit template for CBC builds

# metadata for the Metronome DB
description = build of lalsuite from repo=$(GIT_REPO), branch=$(GIT_BRANCH), id=$(GIT_ID)
project = cbc
project_version = $(GIT_BRANCH)
component = lalsuite
component_version = $(GIT_ID)
run_type = build

# build parameters
platforms = x86_64_sl_6.1
inputs = $(HARNESS_INPUT_SPEC_FILE)
remote_task = remote_task.sh
remote_post = remote_post.py
remote_post_args = --verbose

# send the submitting user an email if any task fails
notify = $(USER)@syr.edu
notify_fail_only = true

# stream "live" stdout/err from tasks back to the submit node
stream_user_io = true

# we set these attributes so that they get propagated to the
# environment of our input spec files and/or remote_* scripts above
git_repo = $(GIT_REPO)
git_id = $(GIT_ID)
git_branch = $(GIT_BRANCH)

# 2013-04-17 pfcouvar: if we submit without RequestMemory and run on a
# slot with <4GB of RAM, the build will run but quickly grow beyond
# that and get preempted (which is normal/good Condor behavior if we
# didn't specify how much we need up front and used more than the slot
# offers) -- but for some reason the job will then just sit idle
# indefinitely.  We need to just specify enough up front.
platform_job_requestmemory = 5000
