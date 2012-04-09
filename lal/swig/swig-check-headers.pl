# check headers for correct SWIGLAL constructs
# Author: Karl Wette, 2011

require 5.008_008;
use strict;
use Getopt::Long;
use File::Basename;
use Cwd qw(getcwd abs_path);
use File::Spec::Functions qw(abs2rel catdir catfile);
use List::Util qw(sum);

sub parse_thru_cpp;

# exit status
my $status = 0;

# parse command line
my ($ifacefile, $headerdir);
GetOptions("interface=s" => \$ifacefile,
           "include=s"   => \$headerdir);
if (!defined($ifacefile) || !defined($headerdir)) {
    print "usage: $0 --interface <file.i> --include <dir>\n";
    exit 1;
}

# check environment
die "environment variable 'CPP' is not set" if !defined($ENV{CPP});

# check input
die "'$ifacefile' is not a file: $!" if (! -f $ifacefile);
die "'$headerdir' is not a directory: $!" if (! -d $headerdir);

# in top-level interface file, look for SWIGLAL_STRUCT_...
# macros, and headers being swiglal_included
my %structmacros;
my @headers;
open FILE, $ifacefile or die "could not open '$ifacefile': $!";
while (<FILE>) {
    chomp;
    if (/^(SWIGLAL_STRUCT_\w+)\((.*?)\);$/) {

        # save macro and struct name
        my $macroname = $1;
        my $structname = $2;

        # print error if there has already been a macro defined
        if (defined($structmacros{$structname})) {
            print "$ifacefile: multiple SWIGLAL_STRUCT_...($structname) macros\n";
            $status = 1;
        }

        # remember that macro has been defined
        $structmacros{$structname} = $macroname;

    }
    if (/^swiglal_include\((.*?)\)$/) {

        # add to list of headers
        push @headers, $1;

    }
}
close FILE;

foreach my $headerfile (@headers) {

    # check that the header exists
    my $file = abs2rel(abs_path(catfile($headerdir, $headerfile)), getcwd());
    die "'$file' is not a file: $!" if (! -f $file);

    # parse the header through the C preprocessor.
    # do not define SWIG so that all SWIG-related code should be removed.
    my $headernoswig = parse_thru_cpp $file;

    # print error for any SWIGLAL macros remaining in the header
    my %noswigmacros;
    while ($headernoswig =~ /(SWIGLAL_\w+)/sg) {
        $noswigmacros{$1} = 1;
    }
    foreach (keys %noswigmacros) {
        print "$file: '$_' macro remains when preprocessed without -DSWIG\n";
        $status = 1;
    }

    # parse the header through the C preprocessor.
    # define SWIG so that code that is excluded from the SWIG interface will be removed.
    my $headerswig = parse_thru_cpp $file, '-DSWIG';

    # find all structs in the header
    while ($headerswig =~ /struct\s+(\S*)\s*{(.*?)}\s*(\S*);/sg) {

        # get the struct (tag) name, and contents
        my $structtagname = $1;
        my $structcode = $2;
        my $structname = $3;

        # print error if struct has no tag name, or tag name does not
        # begin with 'tag', or name does not equal 'tag' + tag name
        if ($structtagname ne '') {
            if ($structtagname !~ /^tag/) {
                print "$file: tag name of struct '$structtagname' does not begin with 'tag'\n";
                $status = 1;
            }
            if ($structname ne '') {
                if ('tag'.$structname ne $structtagname) {
                    print "$file: tag name of struct '$structtagname' does not equal 'tag' + name '$structname'\n";
                    $status = 1;
                }
            }
            else {
                $structname = $structtagname;
                $structname =~ s/^tag//;
            }
        }
        else {
            print "$file: struct '$structname' has no tag name\n";
            $status = 1;
        }

        # print error if struct does not have a SWIGLAL_STRUCT macro,
        # or argument to SWIGLAL_STRUCT does not match struct name
        if ($structcode =~ /SWIGLAL_STRUCT\((.*?)\);/) {
            if ($1 ne $structname) {
                print "$file: struct $structname: SWIGLAL_STRUCT() argument '$1' must equal '$structname'\n";
                $status = 1;
            }
        }
        else {
            print "$file: struct $structname is missing macro SWIGLAL_STRUCT();\n";
            $status = 1;
        }

        # print error if struct does not have a SWIGLAL_STRUCT_... macro defined in top-level interface file
        if (!defined($structmacros{$structname})) {
            print "$ifacefile: missing macro SWIGLAL_STRUCT_...($structname);\n";
            $status = 1;
        }

        # process SWIGLAL_STRUCT_... macro
        delete $structmacros{$structname};

    }

}

# print error if SWIGLAL_STRUCT... macros remain in top-level interface file
foreach (keys %structmacros) {
    print "$ifacefile: unprocessed $structmacros{$_}($_) macro\n";
    $status = 1;
}

exit $status;


# parse $file through the C preprocessor to
# remove comments and preprocessor directives.
sub parse_thru_cpp {
    my ($file, @opts) = @_;

    my $cmd = "sed '/^ *# *include/d' '$file' | $ENV{CPP} @opts - 2>/dev/null";
    open FILE, "$cmd |" or die "could not execute '$cmd': $!";
    my $header = "";
    while (<FILE>) {
        $header .= $_;
    }
    close FILE;

    return $header;
}
