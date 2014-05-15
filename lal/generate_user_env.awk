# Environment script generation
# Author: Karl Wette, 2011--2014

# env[name] holds the value of environment variable 'name'
# set[name,path] is defined if 'path' has already been set in 'name'
# sed[name] holds sed scripts to remove duplicates from 'name'

# set an environment variable to a given value
#   syntax: set NAME value
$1 == "set" {
  if (NF == 3) {
    name = toupper($2)
    env[name] = $3
  }
}

# prepend a value to an environment variable
#   syntax: prepend PATH path1 path2 path3 ...
$1 == "prepend" {
  if (NF >= 3) {
    name = toupper($2)
    for (i = 3; i <= NF; ++i) {
      if ( !( (name,$i) in set ) ) {
        set[name, $i] = 1
        env[name] = env[name] $i ":"
        sed[name] = sed[name] "s|" $i ":||g;"
      }
    }
    env[name] = env[name] "${" name "}"
  }
}

# append a value to an environment variable
#   syntax: append PATH path1 path2 path3 ...
$1 == "append" {
  if (NF >= 3) {
    name = toupper($2)
    for (i = 3; i <= NF; ++i) {
      if ( !( (name,$i) in set ) ) {
        set[name, $i] = 1
        env[name] = ":" env[name] $i
        sed[name] = sed[name] "s|:" $i "||g;"
      }
    }
    env[name] = "${" name "}" env[name]
  }
}

# output environment variables in both C and Bourne shell syntax
END {
  for (name in env) {
    print "csh:if ( ! ${?" name "} ) setenv " name
    print "sh:export " name
    if (sed[name] != "") {
      print "csh:setenv " name " `echo \"${" name "}\" | @SED@ -e '" sed[name] "'`"
      print "sh:" name "=`echo \"${" name "}\" | @SED@ -e '" sed[name] "'`"
    }
    print "csh:setenv " name " \"" env[name] "\""
    print "sh:" name "=\"" env[name] "\""
  }
}
