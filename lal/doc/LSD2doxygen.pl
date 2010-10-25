#!/usr/bin/perl

use strict;
use Switch;
use File::Temp qw/mktemp/;
use Pod::Usage;

my $diffcmd = "diff -y -W 160";

# parse command line
my ($fname, $action, $noclean, $nolatex, $printfn);
while (my $arg = shift @ARGV) {
    if ($arg =~ /^--/p) {
	my $opt = ${^POSTMATCH};
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
    else {
	$fname = $arg;
    }
}
pod2usage(                            -exitval => 0) if $action eq "help";
pod2usage(-message => "No filename!", -exitval => 2) if !defined($fname);
pod2usage(-message => "No action!",   -exitval => 2) if !defined($action);

# temporary file names
my $tmp1 = mktemp("$fname.LSD2doxygen.XXXXX");
my $tmp2 = mktemp("$fname.LSD2doxygen.XXXXX");

# slurp the original file
print "$fname\n" if $printfn;
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
	    print "No changes were made to '$fname'\n";
	}
	else {
	    open DIFF, "| $diffcmd $fname -" or die "'$diffcmd' failed: $!";
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
		system "$diffcmd $tmp1 $tmp2";
		
	    }
	    else {
		print "Code has been modified in '$fname'!\n";
	    }
	}
	
	# or errors?
	my @err;
	foreach (keys %$cpperr) {
	    push @err, $_ if !defined($origcpperr->{$_});
	}
	if (@err) {
	    print "Errors in preprocessing '$fname'!\n";
	    map { print "   $_\n" } @err;
	}
	
    }
    else {
	die "Invalid action '$action'!";
    }
}

# remove temporary files
unlink $tmp1 if -f $tmp1;
unlink $tmp2 if -f $tmp2;

# parse file through the C preprocessor
sub parseThruCPP {
    my ($fname) = @_;
    my ($cpp, $err);

    # proprocess $fname and get output
    system "cpp -E $fname >$tmp1 2>$tmp2";
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
    
    # return no cleanup was asked for
    return $text if $noclean;

    # get rid of any LSD directives, return if there are none
    return $text if 
	(($text =~ s!</?lal(?:LaTeX|Verbatim|ErrTable)[^>]*?>!!sg) == 0);

    # make embedded C comments safe
    while (($text =~ s!\A(.+)/\*!$1/-*!sg) > 0) {}
    while (($text =~ s!\*/(.+)\Z!*-/$1!sg) > 0) {}

    # replace the first and last non-blank lines with doxygen comments
    $text =~ s!\A(\n*)[^\n]*?\n!$1/**\n!sg;
    $text =~ s!\n[^\n]*?(\n*)\Z!\n*/$1!sg;

    # get rid of CVS tags
    $text =~ s!\$(?:Id|Date|Revision)\$!!mg;

    # convert Author: comments to doxygen
    $text =~ s!^(\s*\*?\s*)Author:!$1\\author!mg;

    # try to clean up embedded LaTeX, if asked for
    if (!$nolatex) {
	
	# regexes for:
	# non-line-breaking whitespace
	my $n = "[^\\S\n]";
	# balanced braces
	my $bbr  = qr!({(?:[^{}]++|(?-1))*})!;
	my $wbbr = qr!{((?:[^{}]*$bbr)*[^{}]*)}!;
	# balanced brackets
	my $bbk  = qr!(\[(?:[^[\]]++|(?-1))*\])!;
	my $wbbk = qr!\[((?:[^[\]]*$bbk)*[^[\]]*)\]!;

	# remove these LaTeX commands:
	# environments
	foreach (qw(center document figure obeylines table wrapfigure)) {
	    $text =~ s!\\(?:begin|end)$n*{$_}!!mg;
	}
	# two arguments
	foreach (qw(providecommand)) {
	    $text =~ s!\\$_$bbr$bbr!!mg;
	}
	# two arguments, first optional
	foreach (qw(idx)) {
	    $text =~ s!\\$_$bbk?$bbr!!mg;
	}
	# one argument
	foreach (qw(index input label vfill vspace)) {
	    $text =~ s!\\$_$bbr!!mg;
	}
	# no arguments
	foreach (qw(footnotesize medskip newpage noindent)) {
	    $text =~ s!\\$_ *!!mg;
	}

	# convert formulae
	$text =~ s{(\$\$?)(.+?)\1}{
	    $_ = $2;
	    $_ =~ /\n/ ? '\f[' . $_ . '\f]' : '\f$' . $_ . '\f$'
	}sge;
	$text =~ s!\\begin$n*{(?:equation|displaymath)}!\\f[!mg;
	$text =~ s!\\end$n*{(?:equation|displaymath)}!\\f]!mg;
	$text =~ s!\\begin$n*{eqnarray\*?}!\\f{eqnarray*}{!mg;
	$text =~ s!\\end$n*{eqnarray\*?}!\\f}!mg;

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
	$text =~ s!\{\\(tt|it|rm|sc|sl|bf|sf)$n*!\\text$1\{!sg;
	$text =~ s!\\verb(.)(.+?)\1!\\texttt{$2}!mg;
	$text =~ s!\\emph!\\textit!sg;
	$text =~ s!\\text(?:sc|sl|bf|sf)!\\texttt!sg;
	$text =~ s{\\text(tt|it)$wbbr}{
	    my $e = $1;
	    $_ = $2;
	    s/\\_/_/g;
	    /^[\w_:]+$/ ?
		($e eq 'tt' ? "\\c $_"      : "\\e $_"      ) :
		($e eq 'tt' ? "<tt>$_</tt>" : "<em>$_</em>" )
	}sge;
	
	# replace subsection commands
	$text =~ s!\\(?:sub)*section\*?$wbbr!\\par $1!mg;
	$text =~ s!\\paragraph\*?$wbbr!<b>$1</b>!mg;

	# replace citations
	$text =~ s{\\(?:cite|ref)$wbbr}{
	    $_ = $1;
	    s/://g;
	    '\ref ' . $_
	}mge;

	# miscellaneous
	$text =~ s!\\lq!`!g;
	$text =~ s|(?<!\\)@|\\@|g;

    }

    # get rid of empty comments
    $text =~ s{/\*+(\s*)\*+/}{
	my $wsp = $1;
	$wsp =~ s![^\n]!!g;
	$wsp
    }sge;
    return $text if $text =~ m!\A\n*\Z!sg;

    return $text;
}

__END__

=pod

=head1 SYNOPSIS

LSD2doxygen {action} [debug-options] "file"

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

and [debug-options] are one of:

=over

=item --noclean

Disable LSD to doxygen conversion (i.e. no changes should be made).

=item --nolatex

Do not convert any LaTeX commands to doxygen.

=item --print

Print file name before starting.

=back

=cut
