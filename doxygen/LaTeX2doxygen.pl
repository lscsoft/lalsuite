#!/usr/bin/perl

use strict;
use Switch;
use Pod::Usage;

# display help
pod2usage(-exitval => 1) if @ARGV == 0;

# loop over LaTeX files given on command line
for my $texfname (@ARGV) {

    # check that LaTeX file exists
    die "'$texfname' is not a file!" if !(-f $texfname);
    die "'$texfname' is not a LaTeX file!" if !($texfname =~ /\.tex$/);

    # generate output Doxygen filename
    my $doxyfname = $texfname;
    $doxyfname =~ s/\.tex$/.dox/;

    # check that Doxygen file does not exist
    die "File '$doxyfname' already exists!" if (-f $doxyfname);

    # open files
    open TEX, "<$texfname" or die "Could not open '$texfname'!: $!";
    open DOXY, ">$doxyfname" or die "Could not open '$doxyfname'!: $!";

    # copy header of LaTeX file verbatim to doxygen file
    while (<TEX>) {
        last if $_ =~ /^ *\\begin{document}/;
        print DOXY "// $_";
    }

    # slurp the contents of the LaTeX document
    my $latex = "";
    while (<TEX>) {
        last if $_ =~ /^ *\\end{document}/;
        $latex .= $_;
    }
    close TEX;

    # convert to Doxygen
    my $doxygen = toDoxygen($latex);

    # write to output file
    print DOXY "/**\n\\page pagename\n\n$doxygen\n*/\n";
    close DOXY;

}

exit(0);
    

# convert LaTeX to doxygen documentation
sub toDoxygen {
    my ($text) = @_;

    # regex for non-line-breaking whitespace
    my $n = "[^\\S\n]";

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
    {

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
        $text =~ s!\\\[!\\f[!sg;
        $text =~ s!\\\]!\\f]!sg;
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
        $text =~ s!\\(?:sub)*section\*?{(Synopsis|Description|Prototypes|Algorithm|Notes|Uses|Usage)}!\n### $1 ###\n!g;

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
        $text =~ s!\\paragraph\*?$wbbr!\n### $1 ###\n!mg;

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

        # replace protected space "~" by " "
        $text =~ s!([^\\])~!$1 !g;

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

LaTeX2doxygen "file.tex" ...

Converts each LaTeX "file.tex" to Doxygen, which is written to "file.dox"

=cut
