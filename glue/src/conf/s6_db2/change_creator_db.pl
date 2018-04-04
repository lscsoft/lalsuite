#!/usr/bin/env perl

if($#ARGV!=0)
{
    print "Usage: change_creator_db.pl <site>\n";
    exit(1);
}

$this=$ARGV[0];

print "The site is $this\n";

$tableFileName="creator_db_table.txt";

open(TABLE,"<$tableFileName");

while(<TABLE>)
{
    ($site,$number) = split;
    $table{$site}=$number;
}

$number=$table{$this};
if(!$number)
{
    print "Unknown site.\n";
    print "Here is the list of known sites:\n";
    print join(" ",keys(%table)),"\n";
    exit(2);
}

print "Setting creator_db to $table{$this}\n";

@files=glob("*.sql.in *.sh.in");
map(chomp,@files);

foreach $file (@files)
{
    open(IN,"<$file");
    @lines=<IN>;
    close(IN);

    $_ = $file;
    s/\.in$//;
    $outfile = $_;

    open(OUT,">$outfile");
    foreach $line (@lines)
    {
	if($line =~ m/creator_db/ && $line =~ m/DEFAULT/)
	{
	    print OUT "      creator_db         INTEGER NOT NULL WITH DEFAULT $number,\n";
	}
        elsif($line =~ m/_cdb/ && $line =~ m/DEFAULT/)
        {
            $_ = $line;
            s/DEFAULT 1/DEFAULT $number/;
            print OUT $_;
        }
        elsif($line =~ m/_cdb/ && $line =~ m/ELSE/)
        {
            $_ = $line;
            s/ELSE 1/ELSE $number/;
            print OUT $_;
        }
        elsif(($line =~ m/connect/ || $line =~ m/update/) && $line =~ m/seg_/)
        {
            $_ = $line;
            s/seg_lho/seg_$this/;
            print OUT $_;
        }
	else
	{
	    print OUT $line;
	}
    }
    close(OUT);
}


