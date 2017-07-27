use strict;
use warnings;
use POSIX;
use Getopt::Std;

sub create
{
	my $out = $_[0];
	my $num = $_[1];

	my $mid = floor($num/2);
	my $b = 1.0;
	my $g = 0.0;
	my $r = 0.0;
	my $delta = 1/$mid;

	open(OUT,">$out");
	for(my $i=0; $i<$num; $i++)
	{
		if($i<$mid)
		{
			$g+=$delta;
			$r+=$delta;
		}
		else
		{
			$b-=$delta;
			$g-=$delta;
		}

		print OUT "$b	$g	$r\n";
	}
	close(OUT);
}

if(scalar@ARGV==0)
{
        print "\n\nScript creates 'bwr' colormap.\nCalling: perl create_bwr_colormap.pl \n\nOptions:\n-h this help message\n-o cmap.dat\n-n number_of_rows\n\n";
        exit;
}

my %args = ();
getopts('ho:n:',\%args);

if($args{h})
{
        print "\n\nScript creates 'bwr' colormap.\nCalling: perl create_bwr_colormap.pl \n\nOptions:\n-h this help message\n-o cmap.dat\n-n number_of_rows\n\n";
        exit;
}

my $out = "cmap.dat";
my $num = 64;

if($args{o})
{ $out = $args{o}; }
if($args{n})
{ $num = $args{n}; }

create($out,$num);


