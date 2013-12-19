# Environment script generation
# Author: Karl Wette, 2011

# print an error message to standard error, remove output file, and exit
function msg(str) {
    print "generate_user_env.awk: " str >"/dev/stderr"
    if (output != "") {
        system("rm -f " output)
    }
    exit 1
}

# script setup
BEGIN {
    # check that required variables were set on command line
    if (SED == "") {
        msg("no value for 'SED' given on command line")
    }
    if (output == "") {
        msg("no value for 'output' given on command line")
    }
    # remove output file
    system("rm -f " output)
}

# comments are passed through
/^#/ {
    print $0 >>output
}

# source another environment setup file
#   syntax: source FILE
$1 == "source" {
    for (i = 2; i <= NF; ++i) {
        sourcefiles[$i] = $i
    }
}

# set a environment variable to a given value or path
#   syntax: set NAME value
#           set NAME /dir1 /dir2
$1 == "set" {
    name = $2
    if (name == "") {
        msg("no name given to prepend")
    }
    setvars[name] = ""
    value = ""
    for (i = 3; i <= NF; ++i) {
        found = index(":" setvars[name] ":", ":" $i ":") > 0 ||
            index(":" value ":", ":" $i ":") > 0
        if (!found) {
            value = value ":" $i
        }
    }
    if (value != "") {
        setvars[name] = substr(value, 2)
    }
    delete pathvars[name]
}

# prepend a value to a path environment variable
#   syntax: prepend PATH /dir1 /dir2 /dir3 ...
$1 == "prepend" {
    name = $2
    if (name == "") {
        msg("no name given to prepend")
    }
    value = ""
    for (i = 3; i <= NF; ++i) {
        found = index(":" pathvars[name] ":", ":" $i ":") > 0 ||
            index(":" value ":", ":" $i ":") > 0
        if (!found) {
            value = value $i ":"
        }
    }
    if (value != "") {
        if (pathvars[name] == "") {
            pathvars[name] = sprintf("${%s}", name)
        }
        pathvars[name] = value pathvars[name]
    }
    delete setvars[name]
}

# append a value to a path environment variable
#   syntax: append PATH /dir1 /dir2 /dir3 ...
$1 == "append" {
    name = $2
    if (name == "") {
        msg("no name given to prepend")
    }
    value = ""
    for (i = 3; i <= NF; ++i) {
        found = index(":" pathvars[name] ":", ":" $i ":") > 0 ||
            index(":" value ":", ":" $i ":") > 0
        if (!found) {
            value = value ":" $i
        }
    }
    if (value != "") {
        if (pathvars[name] == "") {
            pathvars[name] = sprintf("${%s}", name)
        }
        pathvars[name] = pathvars[name] value
    }
    delete setvars[name]
}

# output environment variables
END {
    # if output file ends in 'csh', use C shell syntax, otherwise Bourne shell syntax
    csh = (output ~ /\.csh$/)
    # output source files
    for (sourcefile in sourcefiles) {
        if (csh) {
            print "source " sourcefile ".csh" >>output
        }
        else {
            print ". " sourcefile ".sh" >>output
        }
    }
    # output set variables
    for (name in setvars) {
        if (csh) {
            print "setenv " name " \"" setvars[name] "\"" >>output
        }
        else {
            print name "=\"" setvars[name] "\"" >>output
            print "export " name >>output
        }
    }
    # output prepend/append variables
    for (name in pathvars) {
        sed_script = "s|::|:!@%COLON%@!:|;"
        split(pathvars[name], pathvar, ":")
        for (i in pathvar) {
            if (substr(pathvar[i], 1, 1) != "$") {
                sed_script = sed_script "s|:" pathvar[i] ":|:|;"
            }
        }
        sed_script = sed_script "s|::*|:|g;s|^:||;s|:$||;s|!@%COLON%@!||g;"
        if (csh) {
            print "if ( ! ${?" name "} ) setenv " name >>output
            print "setenv " name " `echo \":${" name "}:\" | " SED " '" sed_script "'`" >>output
            print "setenv " name " \"" pathvars[name] "\"" >>output
        }
        else {
            print name "=`echo \":${" name "}:\" | " SED " '" sed_script "'`" >>output
            print name "=\"" pathvars[name] "\"" >>output
            print "export " name >>output
        }
    }
}
