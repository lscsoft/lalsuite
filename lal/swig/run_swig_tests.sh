# Set up environment and run SWIG test scripts under libtool
# Author: Karl Wette, 2014

# record command line
printf "\nCommand:\n%s" "$0"
printf " '%s'" "$@"
printf "\n\n"

# parse command line
echo "Arguments:"
libtool="$1"
echo "libtool='${libtool}'"
shift
libs="$1"
echo "libs='${libs}'"
shift
pathname="$1"
echo "pathname='${pathname}'"
shift
pathval="$1"
echo "pathval='${pathval}'"
shift
cmd="$1"
echo "cmd='${cmd}'"
shift
args=`printf " '%s'" "$@"`
echo "args=\"${args}\""
echo

# set LAL debugging and message level
LAL_DEBUG_LEVEL="memdbg,msglvl1"
export LAL_DEBUG_LEVEL

# set language path
eval "${pathname}=${pathval}; export ${pathname}"

# list directories in language path
for dir in `echo ${pathval} | sed 's|:| |g'`; do
    echo "Contents of path directory $dir:"
    ls -l "$dir"
    echo
done

# recursively build unique list of all libtool libraries
oldlibs="${libs}"
while :; do

    # store the dependency libraries from 'oldlibs' in 'newlibs'
    newlibs=
    for oldlib in ${oldlibs}; do

        # extract 'dependency_libs' value from 'oldlib'
        oldlibdeps=`. ${oldlib}; echo ${dependency_libs}`

        # add each dependency library to 'newlibs', excluding system libraries
        for newlib in ${oldlibdeps}; do
            case ${newlib} in
                /lib/*|/usr/*|/opt/*)
                    ;;
                *.la)
                    if test "x${newlibs}" = x; then
                        newlibs="${newlib}"
                    else
                        newlibs="${newlibs} ${newlib}"
                    fi
                    ;;
            esac
        done

    done

    # if 'newlibs' is empty, no more libraries, so break
    if test "x${newlibs}" = x; then
        break
    fi

    # otherwise, add any library in 'newlibs' to 'libs', ignoring duplicates
    for newlib in ${newlibs}; do
        case " ${libs} " in
            *" ${newlib} "*)
                ;;
            *)
                libs="${libs} ${newlib}"
                ;;
        esac
    done

    # start over by looking for the dependency libraries of 'newlibs'
    oldlibs="${newlibs}"

done

# build libtool command line
ltcmd="${libtool} --mode=execute"
for lib in ${libs}; do
    ltcmd="${ltcmd} -dlopen ${lib}"
done

# print environment under libtool command
paths="/=\$/d;/^PATH=/p;/^${pathname}=/p;/LIBRARY_PATH=/p;/^LAL/p"
echo "Environment:"
eval "${ltcmd} ${SHELL} -c 'set'" | sed -n -e "${paths}"
echo

# print test script command
echo "Test script command:"
echo "${ltcmd} ${cmd} ${args}"
echo

# run test script
echo "Test script output:"
eval "${ltcmd} ${cmd} ${args}"
exitval="$?"

# print test script return status
echo
echo "Test script returned: ${exitval}"
exit ${exitval}
