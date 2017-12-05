# Environment script generation
# Author: Karl Wette, 2011--2014

# src[name] holds environment scripts to source, srclen its length
# env[name] holds the value of environment variable 'name'
# set[name,path] is defined if 'path' has already been set in 'name'
# sed[name] holds sed scripts to remove duplicates from 'name'

# input records are separated by semi-colons
BEGIN {
  RS = ";"
  srclen = 0
}

# first filter out any whitespace-only lines
/^[ \n\t]*$/ {
  next
}

# source an environment script
#   syntax: source script
$1 == "source" {
  if (NF != 2) {
    print "generate_user_env.awk: syntax error in '" $0 "'" >"/dev/stderr"
    exit 1
  }
  src[srclen] = $2
  ++srclen
  next
}

# set an environment variable to a given value
#   syntax: set NAME value
$1 == "set" {
  if (NF != 3) {
    print "generate_user_env.awk: syntax error in '" $0 "'" >"/dev/stderr"
    exit 1
  }
  name = toupper($2)
  env[name] = $3
  next
}

# prepend a value to an environment variable
#   syntax: prepend PATH path1 path2 path3 ...
$1 == "prepend" {
  if (NF == 1) {
    print "generate_user_env.awk: syntax error in '" $0 "'" >"/dev/stderr"
    exit 1
  }
  if (NF >= 3) {
    name = toupper($2)
    for (i = 3; i <= NF; ++i) {
      if ( !( (name,$i) in set ) ) {
        set[name, $i] = 1
        env[name] = env[name] $i ":"
        sed[name] = sed[name] "s|" $i ":||g;"
      }
    }
    env[name] = env[name] "$" name
  }
  next
}

# append a value to an environment variable
#   syntax: append PATH path1 path2 path3 ...
$1 == "append" {
  if (NF == 1) {
    print "generate_user_env.awk: syntax error in '" $0 "'" >"/dev/stderr"
    exit 1
  }
  if (NF >= 3) {
    name = toupper($2)
    for (i = 3; i <= NF; ++i) {
      if ( !( (name,$i) in set ) ) {
        set[name, $i] = 1
        env[name] = ":" env[name] $i
        sed[name] = sed[name] "s|:" $i "||g;"
      }
    }
    env[name] = "$" name env[name]
  }
  next
}

# raise error for any other input
{
  print "generate_user_env.awk: unknown directive '" $0 "'" >"/dev/stderr"
  exit 1
}

# output environment variables in both C and Bourne shell syntax
END {
  envempty = 1
  for (i = 0; i < srclen; ++i) {
    envempty = 0
    print "csh:source " src[i] ".csh"
    print "sh:. " src[i] ".sh"
    print "fish:. " src[i] ".fish"
  }
  for (name in env) {
    envempty = 0
    if (name == "PATH") {
      fisharraytranslate = "| tr ':' '\\n'"
    } else {
      fisharraytranslate = ""
    }
    print "csh:if ( ! ${?" name "} ) setenv " name
    print "sh:export " name
    if (sed[name] != "") {
      print "csh:setenv " name " `echo \"$" name "\" | @SED@ -e '" sed[name] "'`"
      print "sh:" name "=`echo \"$" name "\" | @SED@ -e '" sed[name] "'`"
      print "fish:set " name " (echo \"$" name "\" | @SED@ -e 's| |:|g;" sed[name] "'" fisharraytranslate ")"
    }
    print "csh:setenv " name " \"" env[name] "\""
    print "sh:" name "=\"" env[name] "\""
    print "fish:set " name " (echo \"" env[name] "\" | @SED@ -e 's| |:|g;'" fisharraytranslate ")"
  }
  if (envempty) {
    print "generate_user_env.awk: no user environment script was generated" >"/dev/stderr"
    exit 1
  }
}
