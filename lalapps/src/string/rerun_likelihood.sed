# usage:
#
# sed -f rerun_likelihood.sed <filename.dag >new.dag
#
# edits the .dag file so that all jobs are marked DONE except the
# meas_likelihood jobs and their children which are all marked as not DONE

# find JOB statements
/^[[:space:]]*JOB[[:space:]]/ {
# fine meas_likelihood and calc_likelihood jobs and
# remove DONE if present
/lalapps_string_meas_likelihood/ {
s?\(.*\)[[:space:]]*DONE[[:space:]]*$?\1?
b done
}
/lalapps_string_calc_likelihood/ {
s?\(.*\)[[:space:]]*DONE[[:space:]]*$?\1?
b done
}
# append DONE to the end if not present
/DONE[[:space:]]*$/! s?$? DONE?
}
:done
