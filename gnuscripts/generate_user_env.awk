# Environment script generation
# Author: Karl Wette, 2011

# before processing any file
BEGIN {
    prompt = "user_env.awk:"
    # check that required variables were set on command line
    if (package == "") {
        print prompt, "no value for 'package' given on command line" >"/dev/stderr"
        exit 1
    }
    if (output == "") {
        print prompt, "no value for 'output' given on command line" >"/dev/stderr"
        exit 1
    }
    sourcecount = 0
}

# source another environment setup file
#   syntax: source FILE
$1 == "source" {
    ++sourcecount
    sourcefiles[sourcecount] = $2
}

# set a environment variable to a given value or path
#   syntax: set NAME value
#           set NAME /dir1 /dir2
$1 == "set" {
    name = $2
    if (name == "") {
        print prompt, "no name given to prepend" > "/dev/stderr"
        exit 1
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
    setvars[name] = substr(value, 2)
    delete pathvars[name]
}

# prepend a value to a path environment variable
#   syntax: prepend PATH /dir1 /dir2 /dir3 ...
$1 == "prepend" {
    name = $2
    if (name == "") {
        print prompt, "no name given to prepend" > "/dev/stderr"
        exit 1
    }
    if (pathvars[name] == "") {
        pathvars[name] = sprintf("${%s}", name)
    }
    value = ""
    for (i = 3; i <= NF; ++i) {
        found = index(":" pathvars[name] ":", ":" $i ":") > 0 ||
            index(":" value ":", ":" $i ":") > 0
        if (!found) {
            value = value $i ":"
        }
    }
    pathvars[name] = value pathvars[name]
    delete setvars[name]
}

# append a value to a path environment variable
#   syntax: append PATH /dir1 /dir2 /dir3 ...
$1 == "append" {
    name = $2
    if (name == "") {
        print prompt, "no name given to prepend" > "/dev/stderr"
        exit 1
    }
    if (pathvars[name] == "") {
        pathvars[name] = sprintf("${%s}", name)
    }
    value = ""
    for (i = 3; i <= NF; ++i) {
        found = index(":" pathvars[name] ":", ":" $i ":") > 0 ||
            index(":" value ":", ":" $i ":") > 0
        if (!found) {
            value = value ":" $i
        }
    }
    pathvars[name] = pathvars[name] value
    delete setvars[name]
}

# output environment variables
END {
    # if output file ends in 'csh', use C shell syntax, otherwise Bourne shell syntax
    csh = (output ~ /\.csh$/)
    # print usage
    printf "# source this file to access '%s'\n", package >output
    printf "# usage:\n" >>output
    printf "#   source %s\n", output >>output
    # output source files
    for (i = 1; i <= sourcecount; ++i) {
        if (csh) {
            printf "if ( -f %s.csh ) then\n", sourcefiles[i] >>output
            printf "   source %s.csh $argv:q\n", sourcefiles[i] >>output
            printf "endif\n" >>output
        }
        else {
            printf "if [ -f %s.sh ]; then\n", sourcefiles[i] >>output
            printf "   source %s.sh \"$@\"\n", sourcefiles[i] >>output
            printf "fi\n" >>output
        }
    }
    # output set variables
    for (name in setvars) {
        if (csh) {
            printf "setenv %s \"%s\"\n", name, setvars[name] >>output
        }
        else {
            printf "%s=\"%s\"\n", name, setvars[name] >>output
            printf "export %s\n", name >>output
        }
    }
    # output prepend/append variables
    #   - use sed, if available, to first remove path
    #     elements from path variable if already present
    for (name in pathvars) {
        if (csh) {
            printf "if ( ! ${?%s} ) setenv %s\n", name, name >>output
            printf "which sed >&/dev/null\n" >>output
            printf "if ( ! $status ) then\n" >>output
        }
        else {
            printf "type -p sed >/dev/null\n" >>output
            printf "if [ $? -eq 0 ]; then\n" >>output
        }
        split(pathvars[name], pathvar, ":")
        for (i in pathvar) {
            if (substr(pathvar[i], 1, 1) != "$") {
                if (csh) {
                    printf "   setenv %s `echo \":${%s}:\" | sed 's|:%s:|:|;s|^:||;s|:$||'`\n", name, name, pathvar[i] >>output
                }
                else {
                    printf "   %s=`echo \":${%s}:\" | sed 's|:%s:|:|;s|^:||;s|:$||'`\n", name, name, pathvar[i] >>output
                }
            }
        }
        if (csh) {
            printf "   setenv %s `echo \"%s\" | sed 's|^:||;s|:$||'`\n", name, pathvars[name] >>output
            printf "else\n" >>output
            printf "   setenv %s \"%s\"\n", name, pathvars[name] >>output
            printf "endif\n" >>output
        }
        else {
            printf "   %s=`echo \"%s\" | sed 's|^:||;s|:$||'`\n", name, pathvars[name] >>output
            printf "else\n" >>output
            printf "   %s=\"%s\"\n", name, pathvars[name] >>output
            printf "fi\n" >>output
            printf "export %s\n", name >>output
        }
    }
}
