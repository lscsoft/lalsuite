# ----------------------------------------------------------------------
# LALSuite GitLab-CI: set up build environment
# ----------------------------------------------------------------------

# combine ./configure flags into one variable
# - concatenate all variables with names LCI_CONFIGURE_FLAGS...
# - into one variable named LCI_CONFIGURE_FLAGS
export LCI_CONFIGURE_FLAGS=
for var in $(printf "%s\n" ${!LCI_CONFIGURE_FLAGS*} | sort); do
    if [ "X${LCI_CONFIGURE_FLAGS}" = X ]; then
        LCI_CONFIGURE_FLAGS="${!var}"
    else
        LCI_CONFIGURE_FLAGS="${LCI_CONFIGURE_FLAGS} ${!var}"
    fi
done

# set source data to Gitlab-CI pipeline creation date
# - how hard can it be to parse a date in Unix...?
SOURCE_DATE_EPOCH=$(LANG=C date -j -f "%Y-%m-%dT%H:%M:%SZ" "${CI_PIPELINE_CREATED_AT}" +"%s" 2>/dev/null || true)
if [ "X${SOURCE_DATE_EPOCH}" = X ]; then
    SOURCE_DATE_EPOCH=$(LANG=C date -u -D "%Y-%m-%dT%H:%M:%SZ" -d "${CI_PIPELINE_CREATED_AT}" '+%s' 2>/dev/null || true)
    if [ "X${SOURCE_DATE_EPOCH}" = X ]; then
        SOURCE_DATE_EPOCH=$(LANG=C date -u -d "${CI_PIPELINE_CREATED_AT}" +"%s" 2>/dev/null || true)
        if [ "X${SOURCE_DATE_EPOCH}" = X ]; then
            echo "ERROR: could not set SOURCE_DATE_EPOCH"
            exit 1
        fi
    fi
fi
export SOURCE_DATE_EPOCH

# print environment
echo "========================= build environment =========================="
export -p
echo "----------------------------------------------------------------------"
