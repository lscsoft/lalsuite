#!/bin/sh

# I don't know why which doesn't work on Solaris...
which()
{
  keep_IFS=${IFS}
  IFS=:
  for dir in ${PATH} ; do
    if test -x $dir/$1 ; then
      echo "$dir/$1"
      IFS=${keep_IFS}
      exit 0
    fi
  done
  IFS=${keep_IFS}
  exit 1
}

d=`pwd`
top_builddir=${top_builddir:-"$d"}
export top_builddir
file=${1:-"uninstall.sh"}
RM=`which rm`

$RM -f rm
cat > rm <<EOF
#!/bin/sh
echo "rm \$*" >> $d/$file
exit 0
EOF
chmod +x rm

cat > $file <<EOF
#!/bin/sh
EOF

orig_PATH=$PATH
PATH=$d:$PATH
export PATH

( cd $top_builddir && make uninstall )

PATH=$orig_PATH
export PATH
$RM -f rm

cat >> $file <<EOF
exit 0
EOF

chmod +x $file
exit 0
