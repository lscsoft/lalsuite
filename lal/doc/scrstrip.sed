# get rid of extraneous lines
/^#verbatim/d
/^#\/verbatim/d
/^#ignore/d
/^#\/ignore/d

# remove irrelevant markup commands
s/^##verse/##/
s/^##\/verse/##/

# get rid of extra #s at beginning of line
/^##/s/^#[#]*/#/
