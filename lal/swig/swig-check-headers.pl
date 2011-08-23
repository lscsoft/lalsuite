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
my ($ifacefile, $headerdir, $verbose);
GetOptions("interface=s" => \$ifacefile,
           "include=s"   => \$headerdir,
           "verbose"     => \$verbose);
if (!defined($ifacefile) || !defined($headerdir)) {
    print "usage: $0 [--verbose] --interface <file.i> --include <dir>\n";
    exit 1;
}

# check environment
die "environment variable 'CPP' is not set" if !defined($ENV{CPP});

# check input
die "'$ifacefile' is not a file: $!" if (! -f $ifacefile);
die "'$headerdir' is not a directory: $!" if (! -d $headerdir);

# get list of headers being swiglal_included
my @headers;
open FILE, $ifacefile or die "could not open '$ifacefile': $!";
while (<FILE>) {
    chomp;
    push @headers, $1 if /^swiglal_include\((.*?)\)$/;
}
close FILE;

my ($numlines, $numswiglines, %numswigmacros);

foreach my $headerfile (@headers) {

    # check that the header exists
    my $file = abs2rel(abs_path(catfile($headerdir, $headerfile)), getcwd());
    die "'$file' is not a file: $!" if (! -f $file);

    # parse the header through the C preprocessor.
    # do not define SWIG so that all SWIG-related
    # code should be removed.
    my $headernoswig = parse_thru_cpp $file;

    # count lines of header code
    if ($verbose) {
        foreach (split $/, $headernoswig) {
            ++$numlines;
        }
    }

    # find any SWIGLAL macros remaining in the header
    my %noswigmacros;
    while ($headernoswig =~ /(SWIGLAL_\w+)/sg) {
        $noswigmacros{$1} = 1;
    }
    foreach (keys %noswigmacros) {
        print "$file: '$_' macro remains when preprocessed without -DSWIG\n";
        $status = 1;
    }

    # parse the header through the C preprocessor.
    # define SWIG so that code that is excluded from
    # the SWIG interface will be removed.
    my $headerswig = parse_thru_cpp $file, '-DSWIG';

    # count lines of header code and number of SWIGLAL macros
    if ($verbose) {
        foreach (split $/, $headerswig) {
            ++$numswiglines;
            next if /^\s*#/;
            ++$numswigmacros{$1} if /(SWIGLAL_\w+)\(/;
        }
    }

    # find all structs in the header
    while ($headerswig =~ /struct\s+(\S*)\s*{(.*?)}\s*(\S*);/sg) {

        # get the struct name (or tagname) and contents
        my $structname = $3 ? $3 : $1;
        my $structcode = $2;

        # print if struct does not have a SWIGLAL_STRUCT_LALALLOC macro
        if ($structcode !~ /SWIGLAL_STRUCT_(?:NO_)?LALALLOC\(\);/) {
            print "$file: missing SWIGLAL_STRUCT_LALALLOC in struct $structname\n";
            $status = 1;
        }

        # print if struct contains mismatched SWIGLAL_DYNAMIC_[12]ARRAY macros
        foreach my $n (qw(1 2)) {

            # name and arguments of last macro
            my ($lastmacro, $lastargs);

            while ($structcode =~ /SWIGLAL_DYNAMIC_${n}DARRAY_(BEGIN|END)\(([^\)]+?)\);/sg) {

                # name and arguments of current macro
                my $thismacro = $1;
                my $thisargs = $2;

                # arguments are equivalent up to whitespace
                $thisargs =~ s/\s//g;

                my $mismatched = 0;
                if ($thismacro eq 'BEGIN') {

                    # BEGIN macro is mismatched is the last macro was also a BEGIN
                    $mismatched = 1 if $lastmacro eq 'BEGIN';

                    # store as last macro
                    $lastmacro = $thismacro;
                    $lastargs = $thisargs;

                }
                else {

                    # END macro is mismatched if last macro is not a BEGIN or arguments do not match
                    $mismatched = 1 if $lastmacro ne 'BEGIN' || $lastargs ne $thisargs;

                    # erase last macro
                    $lastmacro = $lastargs = undef;

                }

                if ($mismatched) {
                    print "$file: mismatched SWIGLAL_DYNAMIC_${n}DARRAY macro in struct $structname\n";
                    $status = 1;
                    last;
                }

            }
        }

    }

}

# print header file statistics
if ($verbose) {
    my $numswigmacros = sum(values %numswigmacros);
    print  "$ifacefile:\n";
    printf "% 6i lines of header code\n", $numlines;
    printf "%+6i lines when processed with -DSWIG (%0.1f%%)\n", $numswiglines - $numlines, 100.0 * ($numswiglines - $numlines) / $numlines;
    printf "% 6i use(s) of 'SWIGLAL' macros (%0.1f%%)\n", $numswigmacros, 100.0 * $numswigmacros / $numlines;
    foreach (sort { $a cmp $b } keys %numswigmacros) {
        printf "% 6i use(s) of '$_' macro\n", $numswigmacros{$_};
    }
}

exit $status;


# parse $file through the C preprocessor to
# remove comments and preprocessor directives.
sub parse_thru_cpp {
    my ($file, @opts) = @_;

    my $cmd = "$ENV{CPP} @opts '$file' 2>/dev/null";
    open FILE, "$cmd |" or die "could not execute '$cmd': $!";
    my $header = "";
    my $incl = 0;
    while (<FILE>) {
        # since the C preprocessor might also include some
        # system headers, parse lines beginning with '#' for
        # the header file name; if found, turn on the $incl
        # flag until another '#' line is encountered
        $incl = /\"$file\"/ if /^#/;
        # include line only of $incl is true
        $header .= $_ if $incl;
    }
    close FILE;

    return $header;
}
