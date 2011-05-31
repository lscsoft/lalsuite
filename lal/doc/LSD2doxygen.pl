#!/usr/bin/perl

use strict;
use Switch;
use File::Temp qw/mktemp/;
use Pod::Usage;

# parse command line
my ($fname, $action, $diffcmd, $diffpipe, $noclean, $nolatex, $printfn);
unshift @ARGV, "--diffcmd=diffless";
while (my $arg = shift @ARGV) {
    if ($arg =~ /^--/p) {
        my $opt = ${^POSTMATCH};
        if ($opt =~ /^diffcmd=/p) {
            my $diffopt = ${^POSTMATCH};
            switch ($diffopt) {
                case "diff" {
                    $diffcmd = "diff";
                    $diffpipe = "";
                }
                case "diffless" {
                    $diffcmd = "diff -y -W 160";
                    $diffpipe = "| less";
                }
                case "xxdiff" {
                    $diffcmd = "xxdiff";
                    $diffpipe = "";
                    $action = "diff";
                }
                case "kompare" {
                    $diffcmd = "diff -u5";
                    $diffpipe = "| kompare -o -";
                    $action = "diff";
                }
            }
        }
        else {
            switch ($opt) {
                case "print" {
                    $printfn = 1;
                }
                case "noclean" {
                    $noclean = 1;
                }
                case "nolatex" {
                    $nolatex = 2;
                }
                else {
                    $action = $opt;
                }
            }
        }
    }
    else {
        $fname = $arg;
    }
}
pod2usage(                            -exitval => 0) if $action eq "help";
pod2usage(-message => "No filename!", -exitval => 2) if !defined($fname);
pod2usage(-message => "No action!",   -exitval => 2) if !defined($action);

# exit status
my $status = 0;

# temporary file names
my $tmp1 = mktemp("$fname.LSD2doxygen.XXXXX");
my $tmp2 = mktemp("$fname.LSD2doxygen.XXXXX");

# slurp the original file
print "===== $fname =====\n" if $printfn;
die "'$fname' is not a file!" if !(-f $fname);
my $origfile = "";
open FILE, "<$fname" or die "Could not open '$fname'!: $!";
while (<FILE>) {
    $origfile .= $_;
}
close FILE;

# start with the original file
my $file = $origfile;

# regular expression to match an #if* / #endif CPP block
# based on http://perldoc.perl.org/perlfaq6.html,
# "Can I use Perl regular expressions to match balanced text?"
my $ifdef = qr{(?<BLOCK>                          # name this match
                ^\#\s*?(?:if|ifdef|ifndef).*?\n   # match an #if* directive
                (?:                               # either:
                 (?:                              #    don't match an #if*/#endif directive
                  \n|                             #       blank line
                  [^\#].*?\n|                     #       any line not beginning with #
                  \#\s*?[^ie].*?\n|               #       any non #i*/#e* directive
                  \#\s*?el.*?\n|                  #       #else/#elif directives
                  \#\s*?in.*?\n                   #       #include directives
                 )++                              #    match without backtracking
                 |                                # or
                 (?-1)                            #    must have found a nested #if*/#endif, so
                )*                                #    recurse to beginning of capturing pattern
                \#\s*?endif.*?\n                  # closing #endif directive
               )
              }mx;

# regular expressions to match C / C++ comments, and things that aren't comments / #directives
# taken from http://perldoc.perl.org/perlfaq6.html,
# "How do I use a regular expression to strip C style comments from a file?"
my $c_comment      = qr{(?<BLOCK>/\*[^*]*\*+([^/*][^*]*\*+)*/)}s;
my $cpp_comment    = qr{(?<BLOCK>//([^\\]|[^\n][\n]?)*?\n)}s;
my $something_else = qr{(?<OTHER>
                         "(\\.|[^"\\])*"|         # double-quoted string
                         '(\\.|[^'\\])*'|         # single-quoted string
                         (`(?:[^`']++|(?-1))*')|  # M4-style-quoted string
                         .[^/"'`\\#]*)            # anything else
                       }sx;

# look for LSD documentation in #if 0 / #endif blocks
sub findIf0Blocks {
    my ($text) = @_;
    $text =~ s{$ifdef|$something_else}{
        my ($block, $other) = ($+{BLOCK}, $+{OTHER});
        my $retn;
        # if pattern BLOCK is defined ...
        if (defined($block)) {
            # get first line, middle, and last line of BLOCK
            my ($first, $middle, $last) = ($block =~ m!\A(.*?\n)((?:.*?\n)*)(.*?\n)\Z!sg);
            # if BLOCK begins with '#if 0', process it
            if ($first =~ m!^#if\s+0\s!) {
                $retn = cleanupLSD($block);
            }
            # otherwise look for LSD doc block in middle
            else {
                $retn = $first . findIf0Blocks($middle) . $last;
            }
        }
        # otherwise print anything else that was matched
        else {
            $retn = $other;
        }
        $retn
    }ge;
    return $text;
}
$file = findIf0Blocks($file);

# look for LSD documentation in comments
$file =~ s{$c_comment|$cpp_comment|$something_else}{
    my ($block, $other) = ($+{BLOCK}, $+{OTHER});
    # if pattern BLOCK is defined, process it
    # otherwise print anything else that was matched
    defined($block) ? cleanupLSD($block) : $other
}ge;

# act
switch ($action) {

    # replace the original file with the new doxygen one
    case "doit" {
        open FILE, ">$fname" or die "Could not open '$fname'!: $!";
        print FILE $file;
        close FILE;
    }

    # show what changes were made to the file
    case "diff" {
        if ($file eq $origfile) {
            print "$fname: no changes were made\n";
        }
        else {
            open DIFF, "| $diffcmd $fname - $diffpipe" or die "'$diffcmd' failed: $!";
            print DIFF $file;
            close DIFF;
        }
    }

    # check that file still generates the same C code
    # when parsed through the C preprocessor
    case ["check", "checkdiff"] {

        # parse the original file through the CPP
        my ($origcpp, $origcpperr) = parseThruCPP($fname);

        # replace it with the doxygenated file
        open FILE, ">$fname" or die "Could not open '$fname'!: $!";
        print FILE $file;
        close FILE;

        # parse the doxygenated file through the CPP
        my ($cpp, $cpperr) = parseThruCPP($fname);

        # restore the original file
        open FILE, ">$fname" or die "Could not open '$fname'!: $!";
        print FILE $origfile;
        close FILE;

        # any differences?
        if ($origcpp ne $cpp) {
            if ($action =~ /diff/) {

                # write processed files to temporary files
                open FILE, ">$tmp1" or die "Could not open '$tmp1'!: $!";
                print FILE $origcpp;
                close FILE;
                open FILE, ">$tmp2" or die "Could not open '$tmp2'!: $!";
                print FILE $cpp;
                close FILE;

                # print differences
                system "$diffcmd $tmp1 $tmp2 $diffpipe";

            }
            else {
                print "$fname: code has been modified!\n";
                $status = 1;
            }
        }

        # or errors?
        my @err;
        foreach (keys %$cpperr) {
            push @err, $_ if !defined($origcpperr->{$_});
        }
        if (@err) {
            print "$fname: errors in preprocessing:\n";
            map { print "   $_\n" } @err;
        }

    }
    else {
        die "Invalid action '$action'!";
        $status = 1;
    }
}

# remove temporary files
unlink $tmp1 if -f $tmp1;
unlink $tmp2 if -f $tmp2;

# exit
exit $status;

# parse file through the C preprocessor
sub parseThruCPP {
    my ($fname) = @_;
    my ($cpp, $err);

    # proprocess $fname and get output
    system "gcc -E $fname >$tmp1 2>$tmp2";
    open FILE, "<$tmp1" or die "Could not open '$tmp1'!: $!";
    while (<FILE>) {
        $cpp .= $_;
    }
    close FILE;
    open FILE, "<$tmp2" or die "Could not open '$tmp2'!: $!";
    while (<FILE>) {
        chomp;
        $err->{$_} = 1;
    }
    close FILE;

    return ($cpp, $err);
}

# clean up LSD documentation
sub cleanupLSD {
    my ($text) = @_;

    # regex for non-line-breaking whitespace
    my $n = "[^\\S\n]";

    # return no cleanup was asked for
    return $text if $noclean;

    # return if there are no LSD tags
    return $text if
        (($text =~ m!</?lal(?:LaTeX|Verbatim|ErrTable)[^>]*?>!) == 0);

    # get rid of LSD LaTeX and Verbatim tags
    $text =~ s!</?lal(?:LaTeX|Verbatim)[^>]*?>!!sg;

    # make embedded C comments safe
    while (($text =~ s!\A(.+)/\*!$1/\\*!sg) > 0) {}
    while (($text =~ s!\*/(.+)\Z!*\\/$1!sg) > 0) {}

    # replace first line #if / last line #endif directives with doxygen comments
    $text =~ s!\A#if[^\n]*!/**!;
    $text =~ s!#endif[^\n]*\Z!*/!;

    # replace long first / last line comments with doxygen comments
    $text =~ s!\A/($n|\*)+!/**!;
    $text =~ s!($n|\*)+/\Z!*/!;

    # get rid of any long string of divider characters
    $text =~ s!([-*%+=])\1{4,}!!sg;

    # use CVS Id tag as a hook to place a '\file' command
    $text =~ s!^(\s*\*?\s*)(Revision:\s*)?\$Id\$!\\file!mp;

    # get rid of other CVS tags
    $text =~ s!\$(?:Date|Revision)\$!!mg;

    # convert 'Author:' string to doxygen formatting
    $text =~ s!^(\s*\*?\s*)Author:!$1\\author!mp;

    # convert LSD error table to a doxygen group
    $text =~ s!<lalErrTable[^>]*?>(\s*)\*/!\\name Error Codes \*/ $1/*@\{*/!sp;
    $text =~ s!(/\*+)(\s*)</lalErrTable[^>]*?>!/*@\}*/$2$1!sp;

    # try to clean up embedded LaTeX, if asked for
    if (!$nolatex) {

        # regexes for balanced braces and brackets
        my $bbr  = qr!({(?:[^{}]++|(?-1))*})!;
        my $wbbr = qr!{((?:[^{}]*$bbr)*[^{}]*)}!;
        my $bbk  = qr!(\[(?:[^[\]]++|(?-1))*\])!;
        my $wbbk = qr!\[((?:[^[\]]*$bbk)*[^[\]]*)\]!;

        # regex substitution for illegal \ref characters
        my $illref = sub { $_[0] =~ s![:.()^\-]!_!g; $_[0] };

        # remove these LaTeX commands:
        # environments
        $text =~ s!\\(?:begin|end)$n*{(?:
                   center|document|obeylines
                   )}!!mgx;
        # two arguments
        $text =~ s!\\(?:
                   providecommand
                   )$bbr$bbr!!mgx;
        # two arguments, first optional
        #$text =~ s!\\(?:
        #idx
        #)$bbk?$bbr!!mgx;
        # one argument
        $text =~ s!\\(?:
                   vfill|vspace
                   )$bbr!!mgx;
        # no arguments
        $text =~ s!\\(?:
                   footnotesize|medskip|newpage|noindent
                  )$n*!!mgx;

        # remove these LaTeX commands but leave argument:
        # one argument
        $text =~ s!\\(?:
                   index
                   )$wbbr!$1!mgx;

        # flag these environments for manual intervention
        $text =~ s!\\(begin|end)$n*{(
                   figure|table
                  )}!(MANUAL INTERVENTION $1 $2)!mgpx;

        # convert formulae
        $text =~ s!\$\$(.+?)\$\$!\\f[$1\\f]!sg;
        $text =~ s!\$(.+?)\$!\\f\$$1\\f\$!sg;
        $text =~ s!\\begin$n*{displaymath}!\\f[!mg;
        $text =~ s!\\end$n*{displaymath}!\\f]!mg;
        $text =~ s{\\begin$n*{(equation\*?|eqnarray\*?)}(.*?)\\end$n*{\1}}{
            my $env = $1;
            my $eqn = $2;
            my $anch = '';
            $eqn =~ s{\\label$n*$wbbr}{
                $_ = $1;
                $anch .= '\anchor ' . &$illref($_) . ' ';
                '\label{' . $_ . '}'
            }sge;
            $anch . '\f{' . $env . '}{' . $eqn . '\f}'
        }sge;

        # convert descriptions
        sub desc {
            my ($text) = @_;
            $text =~ s{\\begin$n*{description}(?<LIST>.*?)\\end$n*{description}}{
                $_ = $+{LIST};
                while (/\\item\[/) {
                    s!\\item$wbbk(?<TEXT>.*?)(?<END>\n*\\item\[|\Z)!<dt>$1</dt><dd>$+{TEXT}</dd>$+{END}!sx;
                }
                '<dl>' . desc($_) . '</dl>'
            }sge;
            return $text;
        }
        $text = desc($text);

        # convert numbered and unnumbered lists
        sub list {
            my ($text) = @_;
            $text =~ s{\\begin$n*{(?<ENV>enumerate|itemize)}(?<LIST>.*?)\\end$n*{\k<ENV>}}{
                my $e = $+{ENV};
                $_ = $+{LIST};
                $e =~ s!enumerate!ol!;
                $e =~ s!itemize!ul!;
                while (/\\item/) {
                    s!\\item(?<TEXT>.*?)(?<END>\n*\\item|\Z)!<li>$+{TEXT}</li>$+{END}!sx;
                }
                "<$e>" . list($_) . "</$e>"
            }sge;
            return $text;
        }
        $text = list($text);

        # convert tables
        $text =~ s{\\begin$n*{tabular}$bbr?(?<TABLE>.*?)\\end$n*{tabular}}{
            $_ = $+{TABLE};
            $_ =~ s!\\hline!!sg;
            $_ =~ s!$n*\\\\$n*\n!</td></tr>\n<tr><td>!sg;
            $_ =~ s|$n*(?<!\\)&$n*|</td><td>|sg;
            '<table><tr><td>' . $_ . '</td></tr></table>'
        }sge;

        # convert verbatim
        $text =~ s!\\begin$n*{(?:verbatim|quote)}!\\code!mg;
        $text =~ s!\\end$n*{(?:verbatim|quote)}!\\endcode!mg;

        # can't convert pictures, but save the code
        $text =~ s!\\begin$n*{picture}!\\verbatim!mg;
        $text =~ s!\\end$n*{picture}!\\endverbatim!mg;

        # replace formatting commands
        $text =~ s!\{\\(tt|it|rm|sc|sl|bf|sf)$n+!\\text$1\{!sg;
        $text =~ s!\\verb(.)(.+?)\1!\\texttt{$2}!mg;
        $text =~ s!\\emph!\\textit!sg;
        $text =~ s!\\text(?:sc|sl|bf|sf)!\\texttt!sg;
        $text =~ s{\\text(tt|it)\s*$wbbr}{
            my $e = $1;
            $_ = $2;
            s/\\_/_/g;
            /^[\w_:]+$/ ?
                ($e eq 'tt' ? "\\c $_"      : "\\e $_"      ) :
                ($e eq 'tt' ? "<tt>$_</tt>" : "<em>$_</em>" )
        }sge;

        # special treatment of 'Synopsis/Prototypes/Description/Uses/Algorithm/Notes' sections: turn into 'heading'
        $text =~ s!\\(?:sub)*section\*?{(Synopsis|Description|Prototypes|Algorithm|Notes|Uses|Usage)}!\\heading{$1}!g;

        # rephrase (sub)section commands, turning labels (if present) into anchors
        $text =~ s{\\((?:sub)*section)\*?\s*$wbbr\n(?<LBL>\\label\s*$bbr)?}{
            my $level = $1;
            my $title = $2;
            my $lbl = $+{LBL};
            if (defined($lbl)) {
                $lbl =~ s!\\label\s*$wbbr!$1!;
            } else {
                $lbl = "TODOref";
            }
            $lbl = &$illref($lbl);
            $_ = '\\' . $level . ' ' . $lbl . ' '. $title . "\n";
            $_
        }sge;

        # replace paragraph command by 'heading'
        $text =~ s!\\paragraph\*?$wbbr!\\heading{$1}!mg;

        # preserve references
        $text =~ s{[~ ]*\(*\\(?:eq)?ref$wbbr\)*}{
            $_ = $1;
            '\TODOref{' . &$illref($_) . '}'
        }sge;

        ## intelligent guess about equation references
        $text =~ s!([Ee]qs?|[Ee]quations?)[.\\~]*\\TODOref!\1.\\eqref!mg;
        ## intelligent guess about figure references
        $text =~ s!([Ff]igs?|[Ff]igures?)[.\\~]*\\TODOref!\1.\\figref!mg;
        ## intelligent guess about table references
        $text =~ s!([Tt]ab?|[Tt]ables?)[.\\~]*\\TODOref!\1.\\tableref!mg;

        # replace probable filenames with references
        $text =~ s!<tt>(.*?\.[ch])</tt>!\\ref \1!mg;

        # replace citations by refs
        $text =~ s{\\cite\s*$wbbr}{
            my $ref = $1;
            $ref = &$illref($ref);
            '[\ref ' . $ref .']'
        }mge;

        # replace bibitems by anchors
        $text =~ s{\\bibitem\s*$wbbr}{
            my $ref = $1;
            $ref = &$illref($ref);
            '\anchor ' . $ref . ' <b>[' . $ref . "]</b> "
        }mge;
        # and get rid of 'bibliography'
        $text =~ s!\\begin{thebibliography}{.*}!(MANUAL INTERVENTION begin bibliography)!;
        $text =~ s!\\end{thebibliography}!(MANUAL INTERVENTION end bibliography)!;

        # replace LaTeX's "\_" by "_"
        $text =~ s!\\_!_!g;

        # replace miscellaneous LaTeX commands
        $text =~ s!\\lq!`!g;
        $text =~ s!``|''!"!g;

        # replace \href{} hyperlinks
        $text =~ s!\\href\s*{(.*)}{(.*)}!<a href="$1">$2</a>!g;

        # replace protected space "~" by "\ " which is understood by doxygen
        $text =~ s!([^\\])~!$1\\ !g;

        # remove any empty LaTeX comments
        $text =~ s!^$n*%$n*$!!mg;

        # mark input commands as 'TODO'
        $text =~ s!\\input$wbbr!\TODOinput{$1}!g
    }

    # get rid of empty comments
    $text =~ s!\A/\*/\Z!!;
    $text =~ s{/\*+(\s*)\*+/}{
        my $wsp = $1;
        $wsp =~ s![^\n]!!g;
        $wsp
    }sge;
    return $text if $text =~ m!\A\n*\Z!sg;

    # remove trailing whitespace
    $text =~ s!$n*$!!mg;

    return $text;
}

__END__

=pod

=head1 SYNOPSIS

LSD2doxygen {action} [options] "file"

where {action} (required) is one of:

=over

=item --help

Display this message.

=item --diff

Run the conversion on "file" and display a side-by-side diff of the changes made.

=item --check

Run the original "file" through the C preprocessor, do the conversion,
run the converted "file" through the C preprocessor, and check that
i) the preprocessed files are identical before/after conversion
ii) no change in error messages from the C preprocessor.

=item --checkdiff

Same as --check, but display a side-by-side diff of the preprocessed output before/after conversion

=item --doit

Do the conversion!

=back

and [options] are one of:

=over

=item --diffcmd={diff,diffless,xxdiff,kompare}

Select command to use for displaying diff (default is 'diffless')

=item --noclean

Disable LSD to doxygen conversion (i.e. no changes should be made).

=item --nolatex

Do not convert any LaTeX commands to doxygen.

=item --print

Print file name before starting.

=back

=cut
